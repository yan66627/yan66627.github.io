# =============================================================================
# Blood Pressure Control & Medication Use — Regression Analysis
# Author: Yan Jin
#
# Analyzes blood pressure control in hypertensive patients from a survey dataset.
# Fits logistic regression (BP control) and linear regression (SBP/DBP) with
# adjustment for demographics and comorbidities. Generates publication-ready
# summary tables using gtsummary.
#
# Dependencies: dplyr, gtsummary
# =============================================================================

library(dplyr)
library(gtsummary)

# -----------------------------------------------------------------------------
# 1. Load and subset data to hypertensive patients
# -----------------------------------------------------------------------------
data_raw <- read.csv("Preprocessed_data.csv")

analysis_data <- data_raw %>%
  filter(htn_accaha == "Yes") %>%
  select(
    svy_id, svy_weight_mec,
    bp_med_use,                          # Exposure: antihypertensive medication use
    bp_control_accaha, bp_sys_mean, bp_dia_mean,   # Outcomes
    demo_age_years, demo_gender, demo_race,        # Demographics
    cc_diabetes, cc_ckd, cc_cvd_any               # Comorbidities
  )

# -----------------------------------------------------------------------------
# 2. Recode variables
# -----------------------------------------------------------------------------
analysis_data <- analysis_data %>%
  mutate(
    bp_med_use  = ifelse(bp_med_use == "Yes", 1, 0),
    bp_control  = ifelse(bp_control_accaha == "Yes", 1, 0),
    diabetes    = ifelse(cc_diabetes == "Yes", 1, 0),
    ckd         = ifelse(cc_ckd == "Yes", 1, 0),
    cvd         = ifelse(cc_cvd_any == "Yes", 1, 0),
    gender      = factor(demo_gender),
    race        = factor(demo_race)
  )

# Complete cases only
analysis_data_complete <- analysis_data %>%
  filter(!is.na(bp_med_use), !is.na(bp_control),
         !is.na(bp_sys_mean), !is.na(bp_dia_mean),
         !is.na(demo_age_years))

cat("Analytic sample size:", nrow(analysis_data_complete), "\n")

# -----------------------------------------------------------------------------
# 3. Table 1 — Descriptive statistics by medication use
# -----------------------------------------------------------------------------
table1 <- analysis_data_complete %>%
  select(demo_age_years, gender, race, diabetes, ckd, cvd,
         bp_med_use, bp_control) %>%
  mutate(
    bp_med_use = factor(bp_med_use, levels = c(0, 1), labels = c("No", "Yes")),
    bp_control = factor(bp_control, levels = c(0, 1), labels = c("Uncontrolled", "Controlled")),
    diabetes   = factor(diabetes,   levels = c(0, 1), labels = c("No", "Yes")),
    ckd        = factor(ckd,        levels = c(0, 1), labels = c("No", "Yes")),
    cvd        = factor(cvd,        levels = c(0, 1), labels = c("No", "Yes"))
  ) %>%
  tbl_summary(
    by        = bp_med_use,
    statistic = list(demo_age_years ~ "{mean} \u00b1 {sd}",
                     all_categorical() ~ "{n} ({p}%)"),
    missing   = "no"
  ) %>%
  add_overall() %>%
  bold_labels()

table1

# -----------------------------------------------------------------------------
# 4. Logistic regression — BP control (binary outcome)
# -----------------------------------------------------------------------------
logit_model <- glm(
  bp_control ~ bp_med_use + demo_age_years + gender + race + diabetes + ckd + cvd,
  family = binomial(),
  data   = analysis_data_complete
)

tbl_regression(logit_model, exponentiate = TRUE) %>%
  bold_labels() %>%
  bold_p(t = 0.05)

# -----------------------------------------------------------------------------
# 5. Linear regression — Systolic blood pressure (continuous outcome)
# -----------------------------------------------------------------------------
sbp_model <- lm(
  bp_sys_mean ~ bp_med_use + demo_age_years + gender + race + diabetes + ckd + cvd,
  data = analysis_data_complete
)

tbl_regression(sbp_model) %>%
  bold_labels() %>%
  bold_p(t = 0.05)
