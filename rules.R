# --- PHYSIQUE ---
Age >= 20 & Age <= 80
is.na(Age) | !is.na(AgeDecade)

# --- DÉMOGRAPHIE & CONDITIONNELLES ---
# On utilise la logique : !CONDITION | RESULTAT
!is.na(MaritalStatus)
!is.na(Education)
!is.na(PhysActive)

# --- GENRE ET REPRODUCTION ---
Gender != "male" | (is.na(nPregnancies) & is.na(nBabies))
Gender != "female" | (nPregnancies >= nBabies)

# --- JEUNES FEMMES ---
!(Gender == "female" & Age < 20) | is.na(nPregnancies)
!(Gender == "female" & Age < 20) | is.na(nBabies)

# --- BIOLOGIE & TENSION ---
Gender != "male" | Testosterone > 150
abs(BMI - (Weight / (Height/100)^2)) < 0.1
Pulse > 50 & Pulse < 110
TotChol > 1.5 & TotChol < 15.0




