
# ============================================================
# Data Cleaning and Integration in Official Statistics
# Project script
# ============================================================

# ------------------------------------------------------------
# 1. Load packages and data
# ------------------------------------------------------------

library(dplyr)
library(NHANES)
library(tidyr)
library(validate)
library(validatetools)
library(VIM)
library(deductive)
library(errorlocate)
library(ggplot2)

# Load the NHANES dataset included in the NHANES package
data(NHANES)

# ------------------------------------------------------------
# 2. Build the working dataset
# ------------------------------------------------------------

# The project focuses on adults from the 2011-2012 survey wave.
# Duplicated individuals are removed using the ID variable.
# Only a subset of variables useful for the cleaning exercise is kept.

df_project <- NHANES %>%
  filter(SurveyYr == "2011_12", Age >= 20) %>%
  distinct(ID, .keep_all = TRUE) %>%
  select(
    ID,
    Age,
    AgeDecade,
    Gender,
    Education,
    MaritalStatus,
    Race1,
    nPregnancies,
    nBabies,
    Weight,
    Height,
    BMI,
    Pulse,
    TotChol,
    Testosterone,
    PhysActive
  ) %>%
  mutate(across(where(is.factor), as.character))

# ------------------------------------------------------------
# 3. Define empirical thresholds for outlier detection
# ------------------------------------------------------------

# We use the interquartile range rule to identify unusually low
# or high values for total cholesterol and pulse.

# Cholesterol thresholds
res_chol <- summary(df_project$TotChol)

up_chol <- res_chol[["3rd Qu."]] +
  1.5 * IQR(df_project$TotChol, na.rm = TRUE)

lo_chol <- res_chol[["1st Qu."]] -
  1.5 * IQR(df_project$TotChol, na.rm = TRUE)

# Pulse thresholds
res_pulse <- summary(df_project$Pulse)

up_pulse <- res_pulse[["3rd Qu."]] +
  1.5 * IQR(df_project$Pulse, na.rm = TRUE)

lo_pulse <- res_pulse[["1st Qu."]] -
  1.5 * IQR(df_project$Pulse, na.rm = TRUE)

# ------------------------------------------------------------
# 4. Define validation rules
# ------------------------------------------------------------

# The validation rules combine:
# - range checks;
# - missing value checks;
# - logical consistency rules;
# - physiological plausibility checks;
# - consistency between BMI, weight and height.

rules <- validator(
  age = Age >= 20 & Age <= 80,
  
  age_decade = is.na(Age) | !is.na(AgeDecade),
  
  marital_status = !is.na(MaritalStatus),
  
  education = !is.na(Education),
  
  physical_activity = !is.na(PhysActive),
  
  male_not_pregnant =
    Gender != "male" | (is.na(nPregnancies) & is.na(nBabies)),
  
  female_pregnant =
    Gender != "female" | (nPregnancies >= nBabies),
  
  testosterone =
    Gender != "male" | Testosterone > 150,
  
  bmi =
    abs(BMI - (Weight / (Height / 100)^2)) < 0.1,
  
  pulse =
    Pulse > lo_pulse & Pulse < up_pulse,
  
  cholesterol =
    TotChol > lo_chol & TotChol < up_chol
)

# Remove redundant validation rules if any
rules <- validatetools::remove_redundancy(rules)

# Apply the validation rules to the original dataset
cf <- confront(df_project, rules)

# Summarise the number of records that pass or fail each rule
check <- summary(cf)
print(check[1:7])

# ------------------------------------------------------------
# 5. Impute missing demographic variables
# ------------------------------------------------------------

# AgeDecade can be reconstructed for individuals aged 80,
# because NHANES top-codes age at 80 and the corresponding
# age decade is necessarily "70+".

df_project <- df_project %>%
  mutate(
    AgeDecade_imp = ifelse(Age == 80 & is.na(AgeDecade), TRUE, FALSE),
    AgeDecade = ifelse(AgeDecade_imp == TRUE, "70+", AgeDecade)
  )

# Marital status and education are imputed using hot deck imputation.
# Donor classes are defined using age decade, gender and race.

df_project <- hotdeck(
  df_project,
  variable = c("MaritalStatus", "Education"),
  domain_var = c("AgeDecade", "Gender", "Race1"),
  impNA = TRUE
)

# Visualise the remaining missingness after demographic imputation
aggr(
  df_project[, c(
    "MaritalStatus",
    "Education",
    "MaritalStatus_imp",
    "Education_imp",
    "AgeDecade",
    "AgeDecade_imp"
  )],
  col = c("skyblue", "tomato", "green"),
  delimiter = "_imp",
  numbers = TRUE,
  prop = FALSE,
  cex.axis = 0.8,
  main = "After demographic hot deck imputation"
)

# ------------------------------------------------------------
# 6. Impute anthropometric variables and reconstruct BMI
# ------------------------------------------------------------

# Missing weight and height values are imputed using hot deck imputation.
# Donor classes are defined using age decade and gender.

df_project <- hotdeck(
  df_project,
  variable = c("Weight", "Height"),
  domain_var = c("AgeDecade", "Gender"),
  impNA = TRUE
)

# If BMI is missing, it is reconstructed from weight and height.

df_project <- df_project %>%
  mutate(
    BMI_imp = if_else(is.na(BMI), TRUE, FALSE),
    BMI = if_else(BMI_imp == TRUE, Weight / (Height / 100)^2, BMI)
  )

# ------------------------------------------------------------
# 7. Visualise pulse outliers
# ------------------------------------------------------------

# The dashed horizontal lines correspond to the lower and upper
# thresholds computed with the IQR rule.

df_project %>%
  ggplot(aes(x = Age, y = Pulse, color = Gender)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = lo_pulse, linetype = "dashed", color = "red") +
  geom_hline(yintercept = up_pulse, linetype = "dashed", color = "red") +
  labs(
    title = "Pulse by gender",
    x = "Age",
    y = "Pulse"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# ------------------------------------------------------------
# 8. Treat outliers and impute biological variables
# ------------------------------------------------------------

# Implausible testosterone values for men are set to missing.
# Pulse and cholesterol outliers are also set to missing.
# They are then imputed using k-nearest-neighbour imputation.

df_project <- df_project %>%
  mutate(
    Testosterone = ifelse(
      Gender == "male" & Testosterone <= 150,
      NA_real_,
      Testosterone
    ),
    Pulse = ifelse(
      !is.na(Pulse) & (Pulse <= lo_pulse | Pulse >= up_pulse),
      NA_real_,
      Pulse
    ),
    TotChol = ifelse(
      !is.na(TotChol) & (TotChol <= lo_chol | TotChol >= up_chol),
      NA_real_,
      TotChol
    )
  )

df_project <- kNN(
  df_project,
  variable = c("Pulse", "TotChol", "Testosterone"),
  dist_var = c("Age", "Gender", "Weight", "Height"),
  k = 5,
  imp_var = TRUE
)

# Re-apply the validation rules after this cleaning step
cf <- confront(df_project, rules)
check <- summary(cf)
print(check[1:7])

# ------------------------------------------------------------
# 9. First imputation attempt for pregnancy-related variables
# ------------------------------------------------------------

# Pregnancy-related variables are imputed using kNN.
# Predictors include gender, age, race, education and BMI.

df_project <- kNN(
  df_project,
  variable = c("nPregnancies", "nBabies"),
  dist_var = c("Gender", "Age", "Race1", "Education", "BMI"),
  k = 5,
  imp_var = TRUE
)

# Re-check validation rules after the first pregnancy imputation
cf <- confront(df_project, rules)
check <- summary(cf)
print(check[1:7])

# ------------------------------------------------------------
# 10. Correct logical inconsistencies in pregnancy variables
# ------------------------------------------------------------

# Men should not have pregnancy or baby counts.
# For women, the number of pregnancies should not be lower
# than the number of babies.
# Inconsistent values are set back to missing before re-imputation.

df_project <- df_project %>%
  mutate(
    nPregnancies = ifelse(
      Gender == "male",
      NA_integer_,
      nPregnancies
    ),
    nBabies = ifelse(
      Gender == "male",
      NA_integer_,
      nBabies
    ),
    nPregnancies = ifelse(
      Gender == "female" & nPregnancies < nBabies,
      NA_integer_,
      nPregnancies
    ),
    nBabies = ifelse(
      Gender == "female" & is.na(nPregnancies),
      NA_integer_,
      nBabies
    )
  )

# Re-impute pregnancy-related variables with hot deck imputation.
# Donor classes are defined using age decade, gender and race.

df_project <- hotdeck(
  df_project,
  variable = c("nBabies", "nPregnancies"),
  domain_var = c("AgeDecade", "Gender", "Race1"),
  impNA = TRUE
)

# ------------------------------------------------------------
# 11. Final validation check
# ------------------------------------------------------------

# The final confrontation checks whether the cleaned dataset
# satisfies the validation rules.

cf <- confront(df_project, rules)
check <- summary(cf)
print(check[1:7])

