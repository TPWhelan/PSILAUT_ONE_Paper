# =============================================================================
# Linear mixed-effects models with two-sided permutation tests for a
# repeated-measures dose x group design (placebo, 2 mg and 5 mg psilocybin;
# autistic and non-autistic participants).
#
# Model:   outcome ~ Dose * Group + Visit + Mean_FD + (1 | ID)
#
# Permutation schemes:
#   - Dose effects and dose x group interactions: Freedman-Lane permutation of
#     residuals from a reduced model, permuted within participants.
#   - Group effects: group labels permuted between participants, keeping each
#     participant's sessions together.
#
# Input: a CSV file in long format, one row per session, with columns
#   ID       participant identifier
#   Group    0 = non-autistic, 1 = autistic
#   Dose     0, 2 or 5 (mg)
#   Visit    visit number (1, 2, 3)
#   Mean_FD  mean framewise displacement for the scan
#   <outcome> the connectivity measure to analyse (set OUTCOME_VAR below)
#
# A simulated dataset with this structure (dummy_data.csv) is provided so the
# script can be run end to end; it contains no study data.
#
# Requires: tidyverse, lme4, lmerTest, emmeans
# =============================================================================

#############################
## SETTINGS
#############################

DATA_PATH      <- "dummy_data.csv"
OUTCOME_VAR    <- "Network"                      # outcome column to analyse
ID_VAR         <- "ID"
GROUP_VAR      <- "Group"
DOSE_VAR       <- "Dose"
EXTRA_COVARS   <- c("Visit", "Mean_FD")          # covariates 

DOSE_LEVELS    <- c(0, 2, 5)
DOSE_LABELS    <- c("Placebo", "2mg", "5mg")
GROUP_LEVELS   <- c(0, 1)
GROUP_LABELS   <- c("Non_Autistic", "Autistic")

# Optional exclusion of individual sessions, as a dplyr filter expression,
# e.g. "!(ID == 'P001' & Dose == 0)". Leave "" to keep all sessions.
EXCLUDE_EXPR   <- ""

N_PERM           <- 10                        # number of permutations
SEED_BASE        <- 42                           # random seed
PAIRWISE_ADJ     <- "fdr"                        # adjustment for p_adj columns (within each table only)
PRINT_PARAMETRIC <- TRUE                         # also print parametric (lmerTest/emmeans) results
COMPLETERS_ONLY  <- FALSE                        # TRUE = sensitivity analysis restricted to
                                                 # participants with a session at every dose

#############################
## PACKAGES
#############################
suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(lmerTest)
  library(emmeans)
})

# Satterthwaite degrees of freedom, matching lmerTest's summary() output
emm_options(lmer.df = "satterthwaite")

#############################
## DATA
#############################
data <- read.csv(DATA_PATH, stringsAsFactors = FALSE)

if (nzchar(EXCLUDE_EXPR)) {
  data <- dplyr::filter(data, !! rlang::parse_expr(EXCLUDE_EXPR))
}

need_cols <- unique(c(OUTCOME_VAR, ID_VAR, GROUP_VAR, DOSE_VAR, EXTRA_COVARS))
data <- data %>% dplyr::select(all_of(need_cols)) %>% tidyr::drop_na()

data <- data %>%
  mutate(
    !!ID_VAR    := factor(.data[[ID_VAR]]),
    !!GROUP_VAR := factor(.data[[GROUP_VAR]], levels = GROUP_LEVELS, labels = GROUP_LABELS),
    !!DOSE_VAR  := factor(.data[[DOSE_VAR]],  levels = DOSE_LEVELS,  labels = DOSE_LABELS)
  )

if (COMPLETERS_ONLY) {
  n_doses    <- tapply(data[[DOSE_VAR]], data[[ID_VAR]], function(x) length(unique(x)))
  completers <- names(n_doses)[which(n_doses == length(DOSE_LEVELS))]
  data <- droplevels(data[data[[ID_VAR]] %in% completers, ])
  cat("\nCOMPLETERS ONLY:", length(completers), "participants\n")
  print(table(data[!duplicated(data[[ID_VAR]]), GROUP_VAR]))
}

# Every dose must have sessions in both groups
dose_by_group <- table(data[[DOSE_VAR]], data[[GROUP_VAR]])
cat("\nSessions by dose (rows) and group (columns):\n"); print(dose_by_group)
if (any(dose_by_group == 0)) {
  stop("At least one dose has no sessions in a group. Check the data, EXCLUDE_EXPR and OUTCOME_VAR.")
}

#############################
## MODEL FORMULAS
#############################
collapse_plus <- function(terms) {
  terms <- terms[terms != ""]
  if (length(terms) == 0) "1" else paste(terms, collapse = " + ")
}

rhs_terms <- function(include_interaction = TRUE, include_group = TRUE) {
  if (include_interaction) {
    base <- sprintf("%s*%s", DOSE_VAR, GROUP_VAR)
  } else {
    base <- c(DOSE_VAR, if (include_group) GROUP_VAR else NULL)
  }
  collapse_plus(c(base, EXTRA_COVARS))
}

form_full_int <- as.formula(sprintf("%s ~ %s + (1|%s)", OUTCOME_VAR, rhs_terms(TRUE, TRUE), ID_VAR))
form_red_int  <- as.formula(sprintf("%s ~ %s + (1|%s)", OUTCOME_VAR, rhs_terms(FALSE, TRUE), ID_VAR))
form_full_add <- form_red_int
form_no_dose  <- as.formula(sprintf("%s ~ %s + (1|%s)", OUTCOME_VAR,
                                    collapse_plus(c(GROUP_VAR, EXTRA_COVARS)), ID_VAR))

#############################
## HELPER FUNCTIONS
#############################

# t value for a named fixed-effect coefficient
get_t_for_coef <- function(model, coef_name){
  tab <- summary(model)$coefficients
  if (!coef_name %in% rownames(tab)) {
    stop(sprintf("Coefficient '%s' not found. Available: %s",
                 coef_name, paste(rownames(tab), collapse = ", ")))
  }
  as.numeric(tab[coef_name, "t value"])
}

# Names of the dose main-effect coefficients (e.g. "Dose2mg", "Dose5mg")
find_dose_main_coef_names <- function(model, dose_var = DOSE_VAR){
  rn <- rownames(summary(model)$coefficients)
  dose_main <- rn[grepl(paste0("^", dose_var), rn) & !grepl(":", rn)]
  names(dose_main) <- gsub(paste0("^", dose_var), "", dose_main)
  dose_main
}

# Dose coefficients: Freedman-Lane, residuals of a model without dose terms
# permuted within participants. With non-autistic participants as the reference
# group, these coefficients are the dose effects in the non-autistic group.
perm_dose_coef <- function(dat, coef_name, id_var = ID_VAR, B = N_PERM, seed = 1){
  set.seed(seed)
  m_full    <- lmer(form_full_int, data = dat, REML = TRUE)
  t_obs     <- abs(get_t_for_coef(m_full, coef_name))
  m_no_dose <- lmer(form_no_dose,  data = dat, REML = TRUE)

  resid_null  <- resid(m_no_dose)
  fitted_null <- fitted(m_no_dose)
  split_idx   <- split(seq_len(nrow(dat)), dat[[id_var]])

  t_perm <- numeric(B)
  for (b in seq_len(B)) {
    resid_perm <- resid_null
    for (ii in split_idx) resid_perm[ii] <- sample(resid_null[ii], replace = FALSE)
    dat$Y_perm <- fitted_null + resid_perm
    m_full_p  <- lmer(update(form_full_int, Y_perm ~ .), data = dat, REML = TRUE)
    t_perm[b] <- abs(get_t_for_coef(m_full_p, coef_name))
  }
  list(coef = coef_name, t_obs = t_obs, p_perm = mean(t_perm >= t_obs), t_perm = t_perm)
}

# Group effect: group labels permuted between participants
perm_group_coef <- function(dat, formula, B = N_PERM, seed = 1){
  set.seed(seed)
  m_full     <- lmer(formula, data = dat, REML = TRUE)
  coef_group <- paste0(GROUP_VAR, levels(dat[[GROUP_VAR]])[2])
  t_obs      <- abs(get_t_for_coef(m_full, coef_group))

  subj_map <- dat %>% distinct(!!sym(ID_VAR), !!sym(GROUP_VAR)) %>%
    rename(ID = all_of(ID_VAR), Group = all_of(GROUP_VAR))

  t_perm <- numeric(B)
  for (b in seq_len(B)) {
    map_b <- subj_map %>% mutate(Group_perm = sample(Group)) %>% select(ID, Group_perm)
    d2 <- dat %>% left_join(map_b, by = setNames("ID", ID_VAR)) %>%
      mutate(!!GROUP_VAR := .data$Group_perm) %>% select(-Group_perm)
    m_full_p  <- lmer(formula, data = d2, REML = TRUE)
    t_perm[b] <- abs(get_t_for_coef(m_full_p, coef_group))
  }
  list(coef = coef_group, t_obs = t_obs, p_perm = mean(t_perm >= t_obs), t_perm = t_perm)
}

# Joint dose x group interaction (2 df): Freedman-Lane permutation of the
# likelihood-ratio statistic (models fitted by maximum likelihood)
perm_FL_LRT <- function(form_full, form_red, dat, id_var = ID_VAR, B = N_PERM, seed = 1){
  set.seed(seed)
  m_red  <- lmer(form_red,  data = dat, REML = FALSE)
  m_full <- lmer(form_full, data = dat, REML = FALSE)
  obs    <- anova(m_red, m_full)$Chisq[2]

  resid_red  <- resid(m_red)
  fitted_red <- fitted(m_red)
  split_idx  <- split(seq_len(nrow(dat)), dat[[id_var]])

  perm <- numeric(B)
  for (b in seq_len(B)) {
    resid_perm <- resid_red
    for (ii in split_idx) resid_perm[ii] <- sample(resid_red[ii], replace = FALSE)
    dat$Y_perm <- fitted_red + resid_perm
    m_full_p <- lmer(update(form_full, Y_perm ~ .), data = dat, REML = FALSE)
    m_red_p  <- lmer(update(form_red,  Y_perm ~ .), data = dat, REML = FALSE)
    perm[b]  <- anova(m_red_p, m_full_p)$Chisq[2]
  }
  list(stat_obs = obs, p_perm = mean(perm >= obs), perm = perm)
}

# Dose x group interaction at each dose (change from placebo compared between
# groups): Freedman-Lane, residuals of the additive model permuted within participants
perm_interaction_dosewise <- function(dat, B = N_PERM, seed = 1){
  set.seed(seed)
  m_red  <- lmer(form_red_int,  data = dat, REML = TRUE)
  m_full <- lmer(form_full_int, data = dat, REML = TRUE)

  em_form <- as.formula(paste("~", GROUP_VAR, "|", DOSE_VAR))
  group_diff <- function(model) {
    gd <- as.data.frame(contrast(emmeans(model, em_form), method = "revpairwise"))
    setNames(gd$estimate, as.character(gd[[DOSE_VAR]]))
  }
  gd_obs <- group_diff(m_full)
  c2_obs <- gd_obs[["2mg"]] - gd_obs[["Placebo"]]
  c5_obs <- gd_obs[["5mg"]] - gd_obs[["Placebo"]]

  resid_red  <- resid(m_red)
  fitted_red <- fitted(m_red)
  split_idx  <- split(seq_len(nrow(dat)), dat[[ID_VAR]])

  c2_perm <- numeric(B); c5_perm <- numeric(B)
  for (b in seq_len(B)) {
    resid_perm <- resid_red
    for (ii in split_idx) resid_perm[ii] <- sample(resid_red[ii], replace = FALSE)
    dat$Y_perm <- fitted_red + resid_perm
    m_full_p <- lmer(update(form_full_int, Y_perm ~ .), data = dat, REML = TRUE)
    gd_p <- group_diff(m_full_p)
    c2_perm[b] <- gd_p[["2mg"]] - gd_p[["Placebo"]]
    c5_perm[b] <- gd_p[["5mg"]] - gd_p[["Placebo"]]
  }
  p2 <- mean(abs(c2_perm) >= abs(c2_obs))
  p5 <- mean(abs(c5_perm) >= abs(c5_obs))
  tibble(
    Contrast = c("Dose2mg:Group (vs Placebo)", "Dose5mg:Group (vs Placebo)"),
    Est_obs  = c(c2_obs, c5_obs),
    p_perm   = c(p2, p5),
    p_adj    = p.adjust(c(p2, p5), method = PAIRWISE_ADJ)
  )
}

# Pairwise dose contrasts within one group: model fitted in that group only,
# emmeans contrast t ratios, Freedman-Lane residual permutation within participants
perm_within_group_pairs <- function(dat, group_level, B = N_PERM, seed = 1){
  set.seed(seed)
  d <- dat %>% dplyr::filter(.data[[GROUP_VAR]] == group_level)

  dose_counts <- table(d[[DOSE_VAR]])
  cat(sprintf("Sample sizes for %s: %s\n", group_level,
              paste(names(dose_counts), dose_counts, sep = "=", collapse = ", ")))

  rhs_full  <- collapse_plus(c(DOSE_VAR, EXTRA_COVARS))
  rhs_null  <- collapse_plus(EXTRA_COVARS)
  form_full <- as.formula(sprintf("%s ~ %s + (1|%s)", OUTCOME_VAR, rhs_full, ID_VAR))
  form_null <- as.formula(sprintf("%s ~ %s + (1|%s)", OUTCOME_VAR, rhs_null, ID_VAR))

  m_full <- lmer(form_full, data = d, REML = TRUE)
  m_null <- lmer(form_null, data = d, REML = TRUE)

  pairs_obs     <- summary(pairs(emmeans(m_full, as.formula(paste("~", DOSE_VAR))), adjust = "none"))
  contrasts_obs <- as.character(pairs_obs$contrast)
  t_obs         <- pairs_obs$t.ratio

  resid_null  <- resid(m_null)
  fitted_null <- fitted(m_null)
  split_idx   <- split(seq_len(nrow(d)), d[[ID_VAR]])

  t_perm <- matrix(NA_real_, nrow = B, ncol = length(t_obs),
                   dimnames = list(NULL, contrasts_obs))
  failures <- 0L

  for (b in seq_len(B)) {
    resid_perm <- resid_null
    for (ii in split_idx) resid_perm[ii] <- sample(resid_null[ii], replace = FALSE)
    d$Y_perm <- fitted_null + resid_perm

    pairs_p <- try({
      m_full_p <- lmer(as.formula(sprintf("Y_perm ~ %s + (1|%s)", rhs_full, ID_VAR)),
                       data = d, REML = TRUE)
      summary(pairs(emmeans(m_full_p, as.formula(paste("~", DOSE_VAR))), adjust = "none"))
    }, silent = TRUE)
    if (inherits(pairs_p, "try-error")) { failures <- failures + 1L; next }

    t_perm[b, ] <- pairs_p$t.ratio[match(contrasts_obs, as.character(pairs_p$contrast))]
  }

  if (failures > 0) {
    cat(sprintf("Warning: %d/%d permutations failed for %s\n", failures, B, group_level))
  }

  p_perm <- colMeans(abs(t_perm) >= matrix(abs(t_obs), nrow = B, ncol = length(t_obs), byrow = TRUE),
                     na.rm = TRUE)
  tibble(
    Group    = as.character(group_level),
    Contrast = contrasts_obs,
    Est_obs  = pairs_obs$estimate,
    t_obs    = t_obs,
    p_perm   = as.numeric(p_perm),
    p_adj    = p.adjust(as.numeric(p_perm), method = PAIRWISE_ADJ)
  )
}

#############################
## PRIMARY MODEL
#############################
omni_model <- lmer(form_full_int, data = data, REML = TRUE)

cat("\nFixed effects (beta, SE, Satterthwaite df, t, parametric p):\n")
print(summary(omni_model)$coefficients)

if (PRINT_PARAMETRIC) {
  cat("\nParametric Type III ANOVA (Satterthwaite):\n")
  print(anova(omni_model, type = 3))

  cat("\nParametric within-group dose contrasts (emmeans; uncorrected):\n")
  print(summary(pairs(emmeans(omni_model, as.formula(paste("~", DOSE_VAR, "|", GROUP_VAR))),
                      adjust = "none")))
}

#############################
## PERMUTATION TESTS
#############################

# 1) Dose x group interaction, joint test (2 df)
res_int <- perm_FL_LRT(form_full_int, form_red_int, dat = data, seed = SEED_BASE + 1)
cat("\nPermutation p, dose x group interaction (joint 2-df test):", signif(res_int$p_perm, 4), "\n")

# 1b) Dose x group interaction at each dose
res_int_dosewise <- perm_interaction_dosewise(data, seed = SEED_BASE + 11)
cat("\nPermutation p, dose x group interaction at each dose:\n")
print(res_int_dosewise)

# 2) Dose coefficients (dose effects in the reference group, non-autistic)
coef_names_vec <- find_dose_main_coef_names(omni_model)
res_dose_list  <- lapply(seq_along(coef_names_vec), function(i)
  perm_dose_coef(data, coef_names_vec[i], seed = SEED_BASE + 20 + i))
dose_coef_tbl <- tibble(
  Coefficient = as.character(coef_names_vec),
  t_obs       = sapply(res_dose_list, `[[`, "t_obs"),
  p_perm      = sapply(res_dose_list, `[[`, "p_perm")
) %>% mutate(p_adj = p.adjust(p_perm, method = PAIRWISE_ADJ))
cat("\nPermutation p, dose coefficients (|t| shown):\n")
print(dose_coef_tbl)

# 3) Group effect, additive model
res_group <- perm_group_coef(data, form_full_add, seed = SEED_BASE + 3)
cat("\nPermutation p, group effect (additive model): |t| =", signif(res_group$t_obs, 4),
    "| p =", signif(res_group$p_perm, 4), "\n")

# 3b) Group effect at placebo, interaction model
res_baseline_group <- perm_group_coef(data, form_full_int, seed = SEED_BASE + 6)
cat("\nPermutation p, group effect at placebo (interaction model): |t| =",
    signif(res_baseline_group$t_obs, 4), "| p =", signif(res_baseline_group$p_perm, 4), "\n")

# 4) Pairwise dose contrasts within each group
cat("\nWithin-group pairwise dose contrasts:\n")
pw_nonaut <- perm_within_group_pairs(data, GROUP_LABELS[1], seed = SEED_BASE + 51)
pw_aut    <- perm_within_group_pairs(data, GROUP_LABELS[2], seed = SEED_BASE + 52)
print(pw_nonaut)
print(pw_aut)
