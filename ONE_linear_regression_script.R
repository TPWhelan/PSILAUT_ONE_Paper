# --- Packages ---
library(tidyverse)
library(broom)

# =========================
# CONFIG: edit just this bit
# =========================
DATA_CSV              <- "data.csv"
OUTCOME_VAR           <- "outcome_variable"     # outcome measured at multiple doses
BASELINE_PREDICTOR    <- "baseline_metric"         # baseline-only predictor (taken at Dose 0)
TIMEVARYING_COVAR     <- c("covariate", "covariate")    # covariate measured at both doses (e.g., motion), c("Mean_FD", "Visit")
BASELINE_COVARS_SCALE <- character(0)     # baseline covariates to include SCALED c("Age") if none use character(0) 
BASELINE_COVARS_RAW   <- character(0)    # baseline covariates to include UNSCALED (e.g., factor) c("Sex")  
DOSE_BASE             <- 0            # baseline dose
DOSE_FU               <- 5            # follow-up dose


# --- Load your long-format dataset ---
dat <- read.csv(DATA_CSV)


# Ensure types
dat <- dat %>%
  mutate(
    ID   = as.factor(ID),
    Dose = as.integer(Dose)
  ) %>%
  mutate(across(all_of(BASELINE_COVARS_RAW), as.factor))  # e.g., Sex

# --- Build wide table with Dose 0 & Dose 5 values and covariates from Dose 0 ---
wide_core <- dat %>%
  filter(Dose %in% c(DOSE_BASE, DOSE_FU)) %>%
  mutate(.dose_label = as.character(Dose)) %>%
  select(ID, .dose_label, all_of(c(OUTCOME_VAR, TIMEVARYING_COVAR))) %>%
  pivot_wider(
    names_from  = .dose_label,
    values_from = all_of(c(OUTCOME_VAR, TIMEVARYING_COVAR)),
    names_glue  = "{.value}_{.dose_label}"
  )

# normalize any stray trailing ".0" in names
names(wide_core) <- sub("(.*_\\d+)\\.0$", "\\1", names(wide_core))

# ensure expected columns exist; fail early if not
expected_cols <- c(
  paste0(OUTCOME_VAR, "_", DOSE_BASE),
  paste0(OUTCOME_VAR, "_", DOSE_FU),
  paste0(TIMEVARYING_COVAR, "_", DOSE_BASE),
  paste0(TIMEVARYING_COVAR, "_", DOSE_FU)
)
missing_now <- setdiff(expected_cols, names(wide_core))
if (length(missing_now)) {
  stop(sprintf(
    "After pivot, the following expected columns are missing: %s\nAvailable columns: %s",
    paste(missing_now, collapse = ", "),
    paste(names(wide_core), collapse = ", ")
  ))
}

# 2) Baseline covariates/predictor from Dose 0
pred_0_name <- paste0(BASELINE_PREDICTOR, "_0")
covars_0 <- dat %>%
  filter(Dose == DOSE_BASE) %>%
  transmute(
    ID,
    !!pred_0_name := .data[[BASELINE_PREDICTOR]],
    across(all_of(BASELINE_COVARS_SCALE)),
    across(all_of(BASELINE_COVARS_RAW))
  )

# --- Combine & compute difference ---
outcome_0 <- paste0(OUTCOME_VAR, "_", DOSE_BASE)
outcome_1 <- paste0(OUTCOME_VAR, "_", DOSE_FU)
fd_0      <- paste0(TIMEVARYING_COVAR, "_", DOSE_BASE)
fd_1      <- paste0(TIMEVARYING_COVAR, "_", DOSE_FU)
diff_name <- paste0(OUTCOME_VAR, "_diff")

fo_cl4_psilaut_block_wide <- wide_core %>%
  left_join(covars_0, by = "ID") %>%
  mutate(!!diff_name := .data[[outcome_1]] - .data[[outcome_0]])

# Optional: drop rows missing outcome or key predictors
needed_cols <- c(diff_name, pred_0_name,
                 BASELINE_COVARS_SCALE, BASELINE_COVARS_RAW,
                 fd_0, fd_1)
fo_cl4_psilaut_block_wide <- fo_cl4_psilaut_block_wide %>%
  filter(if_all(all_of(needed_cols), ~ !is.na(.x)))

# --- Linear model ---
lhs <- paste0("scale(", diff_name, ")")
rhs_terms <- c(
  paste0("scale(", pred_0_name, ")"),
  if (length(BASELINE_COVARS_SCALE)) paste0("scale(", BASELINE_COVARS_SCALE, ")"),
  if (length(BASELINE_COVARS_RAW))   BASELINE_COVARS_RAW,
  paste0("scale(", fd_1, ")"),
  paste0("scale(", fd_0, ")")
)
form <- as.formula(paste(lhs, "~", paste(unlist(rhs_terms), collapse = " + ")))

fit <- lm(data = fo_cl4_psilaut_block_wide, formula = form)

# --- Output ---
summary(fit)
tidy_fit <- broom::tidy(fit, conf.int = TRUE)
print(tidy_fit)

# Partial r for the baseline predictor
coef_name <- paste0("scale(", pred_0_name, ")")
t_val <- coef(summary(fit))[coef_name, "t value"]
df_res <- df.residual(fit)
r_partial <- t_val / sqrt(t_val^2 + df_res)
r_partial


# FDR correction
p_values <- c()
p_fdr <- p.adjust(p_values, method = "BH")
p_fdr