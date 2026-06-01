# ============================================================
# Data Cleaning and Integration in Official Statistics
# Part 2 Project: Statistical Matching on Palmer Penguins
# ============================================================

# This script applies a micro non-parametric statistical matching approach
# to the Palmer Penguins dataset.
#
# The objective is to reconstruct missing joint information between:
# - body mass, observed only in file A;
# - bill length, observed only in file B.
#
# The original dataset is complete, which allows us to artificially hide
# one variable and then evaluate the quality of the imputation.

# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------

library(palmerpenguins)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(knitr)

# Set a seed to make the random split and random donor selection reproducible
set.seed(123)

# ------------------------------------------------------------
# 2. Load and prepare the complete dataset
# ------------------------------------------------------------

# We keep only the variables needed for the statistical matching exercise.
# Observations with missing values are removed to simplify the experiment.

penguins_clean <- penguins %>%
  select(
    species,
    island,
    sex,
    flipper_length_mm,
    body_mass_g,
    bill_length_mm
  ) %>%
  na.omit()

# ------------------------------------------------------------
# 3. Build two disjoint files A and B
# ------------------------------------------------------------

# The complete dataset is split horizontally:
# - file A contains one subset of penguins;
# - file B contains the remaining subset of penguins.
#
# The two files therefore contain different individuals.

n <- nrow(penguins_clean)

sample_A_id <- sample(
  1:n,
  size = floor(n / 2),
  replace = FALSE
)

A_full <- penguins_clean[sample_A_id, ]
B_full <- penguins_clean[-sample_A_id, ]

# ------------------------------------------------------------
# 4. Create the statistical matching setting
# ------------------------------------------------------------

# The dataset is then split vertically:
#
# X = common variables observed in both A and B:
#     species, island, sex, flipper length
#
# Y = variable observed only in A:
#     body mass
#
# Z = variable observed only in B:
#     bill length
#
# In the experiment, the true bill length in A is hidden during matching,
# but stored separately for evaluation.

A <- A_full %>%
  select(
    species,
    island,
    sex,
    flipper_length_mm,
    body_mass_g
  )

B <- B_full %>%
  select(
    species,
    island,
    sex,
    flipper_length_mm,
    bill_length_mm
  )

# True hidden values of bill length in A, used only for evaluation
A_true_bill <- A_full$bill_length_mm

# ------------------------------------------------------------
# 5. Define the hot deck matching function
# ------------------------------------------------------------

# The function imputes bill_length_mm in the recipient file A using
# observed bill_length_mm values from the donor file B.
#
# Two specifications are implemented:
#
# - "rich": donors are restricted to the same species, sex and island.
#           Among them, the donor with the closest flipper length is selected.
#
# - "poor": donors are restricted only to the same species.
#           One donor is selected randomly.
#
# The rich specification is the baseline method used in the project.
# The poor specification is used later for sensitivity analysis.

hot_deck_match <- function(recipient_file, donor_file, spec = "rich") {
  
  imputed_values <- numeric(nrow(recipient_file))
  
  for (i in 1:nrow(recipient_file)) {
    
    # Current recipient unit in file A
    recipient <- recipient_file[i, ]
    
    # Rich specification:
    # use a more restrictive donation class based on species, sex and island
    if (spec == "rich") {
      possible_donors <- donor_file %>%
        filter(
          species == recipient$species,
          sex == recipient$sex,
          island == recipient$island
        )
    }
    
    # Poor specification:
    # use only species as the matching variable
    if (spec == "poor") {
      possible_donors <- donor_file %>%
        filter(
          species == recipient$species
        )
    }
    
    # Fallback rule:
    # if no donor is found in the initial class, relax the matching rule
    # and keep only species.
    if (nrow(possible_donors) == 0) {
      possible_donors <- donor_file %>%
        filter(species == recipient$species)
    }
    
    # Final fallback:
    # if still no donor is available, use the whole donor file.
    if (nrow(possible_donors) == 0) {
      possible_donors <- donor_file
    }
    
    # In the rich specification, select the donor with the closest flipper length.
    if (spec == "rich") {
      best_donor <- possible_donors %>%
        mutate(
          distance = abs(flipper_length_mm - recipient$flipper_length_mm)
        ) %>%
        slice_min(distance, n = 1, with_ties = TRUE) %>%
        slice_sample(n = 1)
    }
    
    # In the poor specification, select one donor randomly.
    if (spec == "poor") {
      best_donor <- possible_donors %>%
        slice_sample(n = 1)
    }
    
    # Assign the donor's bill length to the recipient
    imputed_values[i] <- best_donor$bill_length_mm
  }
  
  # Return the recipient file completed with the imputed bill length
  recipient_file %>%
    mutate(bill_length_imputed = imputed_values)
}

# ------------------------------------------------------------
# 6. Define an evaluation function
# ------------------------------------------------------------

# Since the original data were complete, the true hidden values of bill length
# are known. This allows direct evaluation of the imputation.
#
# We compute:
# - RMSE: root mean squared error;
# - MAE: mean absolute error;
# - true correlation between body mass and true bill length;
# - synthetic correlation between body mass and imputed bill length;
# - absolute correlation error.

evaluate_matching <- function(A_file, B_file, A_true_values, spec = "rich") {
  
  matched_data <- hot_deck_match(A_file, B_file, spec = spec) %>%
    mutate(
      bill_length_true = A_true_values,
      error = bill_length_imputed - bill_length_true
    )
  
  true_cor <- cor(
    matched_data$body_mass_g,
    matched_data$bill_length_true,
    use = "complete.obs"
  )
  
  synthetic_cor <- cor(
    matched_data$body_mass_g,
    matched_data$bill_length_imputed,
    use = "complete.obs"
  )
  
  tibble(
    specification = spec,
    rmse = sqrt(mean(matched_data$error^2, na.rm = TRUE)),
    mae = mean(abs(matched_data$error), na.rm = TRUE),
    true_correlation = true_cor,
    synthetic_correlation = synthetic_cor,
    correlation_error = abs(true_cor - synthetic_cor)
  )
}

# ------------------------------------------------------------
# 7. Apply the baseline statistical matching procedure
# ------------------------------------------------------------

# The baseline method is the rich distance hot deck:
# same species, same sex, same island, then nearest flipper length.

A_matched <- hot_deck_match(A, B, spec = "rich") %>%
  mutate(
    bill_length_true = A_true_bill,
    error = bill_length_imputed - bill_length_true
  )

# Compute baseline evaluation indicators
baseline_results <- evaluate_matching(
  A_file = A,
  B_file = B,
  A_true_values = A_true_bill,
  spec = "rich"
)

print(baseline_results)

# ------------------------------------------------------------
# 8. Compare true and imputed distributions
# ------------------------------------------------------------

# A good statistical matching procedure should preserve the marginal
# distribution of the imputed variable.

A_matched %>%
  select(bill_length_true, bill_length_imputed) %>%
  pivot_longer(
    cols = everything(),
    names_to = "variable",
    values_to = "bill_length_mm"
  ) %>%
  ggplot(aes(x = bill_length_mm, linetype = variable)) +
  geom_density(linewidth = 1) +
  labs(
    title = "True vs imputed distribution of bill length",
    x = "Bill length (mm)",
    y = "Density",
    linetype = ""
  ) +
  theme_minimal()

# ------------------------------------------------------------
# 9. Compare true and synthetic relationships between Y and Z
# ------------------------------------------------------------

# The core objective of statistical matching is not only to impute
# plausible individual values, but also to reconstruct the relationship
# between variables that are not jointly observed.
#
# Here, we compare:
# - the true relationship between body mass and bill length;
# - the synthetic relationship between body mass and imputed bill length.

A_matched %>%
  select(body_mass_g, bill_length_true, bill_length_imputed) %>%
  pivot_longer(
    cols = c(bill_length_true, bill_length_imputed),
    names_to = "variable",
    values_to = "bill_length_mm"
  ) %>%
  ggplot(aes(x = body_mass_g, y = bill_length_mm)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ variable) +
  labs(
    title = "True vs synthetic relationship",
    x = "Body mass (g)",
    y = "Bill length (mm)"
  ) +
  theme_minimal()

# ------------------------------------------------------------
# 10. Sensitivity analysis
# ------------------------------------------------------------

# After applying the baseline procedure, we test whether matching quality
# depends on:
# - the variables used to define similarity;
# - the available sample size.
#
# The sample_size parameter refers to the total number of penguins drawn
# before splitting the data into A and B.
#
# Example:
# - sample_size = 100 means about 50 penguins in A and 50 in B.
# - sample_size = 300 means about 150 penguins in A and 150 in B.

run_one_simulation <- function(sample_size, specification) {
  
  # Draw a subsample from the complete dataset
  sampled_data <- penguins_clean %>%
    slice_sample(n = sample_size)
  
  # Split this subsample into two disjoint files
  sample_A_id <- sample(
    1:nrow(sampled_data),
    size = floor(sample_size / 2),
    replace = FALSE
  )
  
  A_full_sim <- sampled_data[sample_A_id, ]
  B_full_sim <- sampled_data[-sample_A_id, ]
  
  # Build file A: common variables + body mass
  A_sim <- A_full_sim %>%
    select(
      species,
      island,
      sex,
      flipper_length_mm,
      body_mass_g
    )
  
  # Build file B: common variables + bill length
  B_sim <- B_full_sim %>%
    select(
      species,
      island,
      sex,
      flipper_length_mm,
      bill_length_mm
    )
  
  # Apply matching and return evaluation indicators
  evaluate_matching(
    A_file = A_sim,
    B_file = B_sim,
    A_true_values = A_full_sim$bill_length_mm,
    spec = specification
  ) %>%
    select(rmse, mae, correlation_error)
}

# Run the simulation:
# - 5 sample sizes;
# - 20 repetitions for each setting;
# - 2 specifications: poor and rich.

simulation_results <- expand_grid(
  sample_size = c(50, 100, 150, 250, 300),
  repetition = 1:20,
  specification = c("poor", "rich")
) %>%
  mutate(
    result = pmap(
      list(sample_size, specification),
      ~ run_one_simulation(
        sample_size = ..1,
        specification = ..2
      )
    )
  ) %>%
  unnest(result)

# Summarise the results across repetitions
simulation_summary <- simulation_results %>%
  group_by(sample_size, specification) %>%
  summarise(
    mean_rmse = mean(rmse, na.rm = TRUE),
    mean_mae = mean(mae, na.rm = TRUE),
    mean_correlation_error = mean(correlation_error, na.rm = TRUE),
    .groups = "drop"
  )

print(simulation_summary)

# ------------------------------------------------------------
# 11. Plot sensitivity results
# ------------------------------------------------------------

# This graph shows how the average RMSE changes with sample size
# and matching specification.

ggplot(
  simulation_summary,
  aes(x = sample_size, y = mean_rmse, linetype = specification)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(
    title = "Matching quality by sample size",
    x = "Total sample size before splitting into A and B",
    y = "Average RMSE",
    linetype = "Specification"
  ) +
  theme_minimal()