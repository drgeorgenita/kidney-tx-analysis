# Kidney transplant recipients aged <60 vs >=60 years:
# baseline characteristics and clinical outcomes with OR and RR

source("R/helpers.R")
install_missing(c("openxlsx", "dplyr", "tableone", "flextable", "officer",
                  "survival", "survminer"))
library(openxlsx); library(dplyr); library(tableone)
library(flextable); library(officer); library(survival); library(survminer)

in_file <- "data/cohort.xlsx"
out <- "output/age_outcomes"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
op <- function(f) file.path(out, f)

COLS <- c("<60" = "#2166ac", ">=60" = "#d6604d")

# ---- Data -------------------------------------------------------------------
raw <- read.xlsx(in_file, na.strings = c("", "NA", "N/A"))

cd <- function(x, g) as.integer(as.character(x) %in% g)

df <- raw %>%
  filter(!is.na(over_60)) %>%
  mutate(
    Group = factor(over_60, levels = c(0, 1), labels = c("<60", ">=60")),

    rec_sex         = factor(rec_sex),
    ethnicity       = factor(ethnicity),
    PRD_group       = factor(PRD_group),
    dialysis_status = factor(dialysis_status),
    imm_risk        = factor(imm_risk),
    UKKRRI_quartile = factor(UKKRRI_quartile, levels = c("R1", "R2", "R3", "R4")),
    HLA_MM_level    = factor(HLA_MM_level),
    UKKDR_quartlie  = factor(UKKDR_quartlie,
                             levels = c("D1", "D2", "D3", "D4", "live donor")),
    tx_type         = factor(tx_type),

    # baseline comorbidities reported as n (%)
    rec_HTN            = factor(rec_HTN),
    rec_CV_disease     = factor(rec_CV_disease),
    rec_diabetes       = factor(rec_diabetes),
    rec_smoking_status = factor(rec_smoking_status),
    pre_tx_DSA         = factor(pre_tx_DSA),

    SSI_binary      = as.integer(!is.na(CDC_SSI_30_day) &
                                   CDC_SSI_30_day != 0 & CDC_SSI_30_day != "0"),
    SSI_superficial = as.integer(CDC_SSI_30_day == "superficial incisional SSI"),
    SSI_organ       = as.integer(CDC_SSI_30_day == "organ/space SSI"),
    SSI_deep        = as.integer(CDC_SSI_30_day == "deep SSI/abscess"),

    CD_grade = factor(post_op_clavien_dindo_grade,
                      levels = c("0", "1", "2", "3a", "3b", "4a")),
    CD_0     = cd(post_op_clavien_dindo_grade, "0"),
    CD_1     = cd(post_op_clavien_dindo_grade, "1"),
    CD_2     = cd(post_op_clavien_dindo_grade, "2"),
    CD_3a    = cd(post_op_clavien_dindo_grade, "3a"),
    CD_3b    = cd(post_op_clavien_dindo_grade, "3b"),
    CD_3ab   = cd(post_op_clavien_dindo_grade, c("3a", "3b")),
    CD_4a    = cd(post_op_clavien_dindo_grade, "4a"),
    CD_4ab   = cd(post_op_clavien_dindo_grade, c("4a", "4b")),
    CD_3plus = cd(post_op_clavien_dindo_grade, c("3a", "3b", "4a", "4b")),

    across(c(age_at_tx, rec_BMI, charlson_comorbidity_index, CIT_min, donor_age,
             cr_1m, cr_3m, cr_6m, cr_12m, egfr_1m, egfr_3m, egfr_6m, egfr_12m,
             LOS_days, survival_days, life_of_tx_days), as.numeric),

    across(c(DGF, PNF, ITU_admission, graft_failure, recipient_deceased,
             CMV_viraemia, BK_viraemia, EBV_viraemia,
             UTI_30_day, Ischaemic_CV_even_30_day, DVT_PE_30_day,
             LRTI_HAP_30_day, readmission_30_day, intervention_for_SSI,
             reintervention_30_days, fall_12m, malignancy_12mo,
             de_novoDSA, BPAR, rejection), as.numeric)
  )

n_all <- nrow(df); n_lt <- sum(df$over_60 == 0); n_ge <- sum(df$over_60 == 1)
lab_all <- sprintf("Total (n=%d)", n_all)
lab_lt  <- sprintf("Age <60 (n=%d)", n_lt)
lab_ge  <- sprintf("Age ≥60 (n=%d)", n_ge)
cat(sprintf("n=%d | <60: %d | >=60: %d\n", n_all, n_lt, n_ge))

# ---- Table 1: baseline characteristics --------------------------------------
vars_demo <- c(
  "age_at_tx", "rec_sex", "ethnicity",
  "charlson_comorbidity_index", "rec_BMI",
  "rec_HTN", "rec_CV_disease", "rec_diabetes", "rec_smoking_status",
  "PRD_group", "dialysis_status", "imm_risk", "pre_tx_DSA",
  "UKKRRI_quartile", "HLA_MM_level",
  "donor_age", "tx_type", "CIT_min", "UKKDR_quartlie"
)
vars_demo <- vars_demo[vars_demo %in% names(df)]
cat_vars  <- vars_demo[sapply(df[vars_demo], function(x) is.factor(x) || is.character(x))]
cont_vars <- vars_demo[sapply(df[vars_demo], is.numeric)]

# Shapiro-Wilk: normal -> mean (SD) and t-test; otherwise median (IQR) and Mann-Whitney
shapiro_p <- sapply(cont_vars, function(v) {
  x <- df[[v]][!is.na(df[[v]])]
  if (length(x) < 3) return(NA_real_)
  tryCatch(shapiro.test(x)$p.value, error = function(e) NA_real_)
})
print(round(shapiro_p, 4))
nonnormal_vars <- names(shapiro_p)[!is.na(shapiro_p) & shapiro_p < 0.05]

# Fisher exact where any expected cell < 5, otherwise chi-squared
exact_vars <- cat_vars[sapply(cat_vars, function(v) {
  tryCatch(any(suppressWarnings(chisq.test(table(df[[v]], df$Group)))$expected < 5),
           error = function(e) TRUE)
})]

t1 <- CreateTableOne(vars = vars_demo, strata = "Group", data = df,
                     factorVars = cat_vars, addOverall = TRUE)
t1_printed <- print(t1, nonnormal = nonnormal_vars, exact = exact_vars,
                    showAllLevels = FALSE, quote = FALSE, noSpaces = TRUE,
                    printToggle = TRUE, smd = FALSE, test = TRUE)

get_test <- function(v) {
  if (v %in% cat_vars) return(if (v %in% exact_vars) "Fisher exact" else "Chi-squared")
  if (v %in% nonnormal_vars) "Mann-Whitney U" else "Student t-test"
}
test_col <- sapply(rownames(t1_printed), function(rn) {
  rn_clean <- trimws(gsub("\\s*\\(.*\\)", "", rn))
  matched <- vars_demo[sapply(vars_demo, function(v) startsWith(tolower(rn_clean), tolower(v)))]
  if (length(matched)) get_test(matched[1]) else ""
})

t1_df <- as.data.frame(t1_printed, stringsAsFactors = FALSE)
t1_df[["test"]] <- NULL
t1_df <- cbind(Variable = rownames(t1_df), t1_df, `Statistical test` = test_col)

# age defines the groups, so its p-value is not reported
age_row <- which(startsWith(tolower(trimws(gsub("\\s*\\(.*\\)", "", t1_df$Variable))),
                            "age_at_tx"))
if (length(age_row)) {
  t1_df$p[age_row] <- ""
  t1_df[["Statistical test"]][age_row] <- "Grouping variable"
}

p_num <- suppressWarnings(as.numeric(t1_df$p))
sig1  <- which(t1_df$p == "<0.001" | (!is.na(p_num) & p_num < 0.05))

ft1 <- flextable(t1_df) %>%
  set_header_labels(Variable = "Variable", Overall = lab_all, `<60` = lab_lt,
                    `>=60` = lab_ge, p = "p-value",
                    `Statistical test` = "Statistical test") %>%
  bold(part = "header") %>%
  fontsize(size = 10, part = "all") %>%
  font(fontname = "Calibri", part = "all") %>%
  border_inner_h(part = "body") %>%
  add_footer_lines(paste0(
    "Values: mean (SD) for normally distributed continuous variables, median (IQR) ",
    "otherwise; n (%) for categorical. Normality assessed by Shapiro-Wilk. ",
    "Continuous: Student t-test or Mann-Whitney U; categorical: chi-squared, or ",
    "Fisher exact where any expected cell <5. Highlighted rows: p<0.05. ",
    "Recipient age defines the study groups; its p-value is not reported.")) %>%
  fontsize(size = 9, part = "footer") %>%
  italic(part = "footer") %>%
  autofit()
if (length(sig1)) ft1 <- bg(ft1, i = sig1, bg = "#FFF2CC")

read_docx() %>%
  body_add_par("Table 1. Baseline demographics and transplant characteristics",
               style = "heading 1") %>%
  body_add_flextable(ft1) %>%
  print(target = op("table1_baseline.docx"))

# ---- Table 2: clinical outcomes with OR and RR ------------------------------
# OR and RR for >=60 relative to <60, Wald intervals on the log scale
calc_or <- function(e_ge, n_ge, e_lt, n_lt) {
  a <- e_ge; b <- n_ge - e_ge; c <- e_lt; d <- n_lt - e_lt
  if (a == 0 || b == 0 || c == 0 || d == 0) return(list(est = NA, lo = NA, hi = NA))
  est <- (a * d) / (b * c); se <- sqrt(1/a + 1/b + 1/c + 1/d)
  list(est = est, lo = exp(log(est) - 1.96 * se), hi = exp(log(est) + 1.96 * se))
}
calc_rr <- function(e_ge, n_ge, e_lt, n_lt) {
  if (e_lt == 0 || n_lt == 0 || n_ge == 0) return(list(est = NA, lo = NA, hi = NA))
  est <- (e_ge / n_ge) / (e_lt / n_lt)
  se  <- sqrt(1/max(e_ge, 0.5) - 1/n_ge + 1/max(e_lt, 0.5) - 1/n_lt)
  list(est = est, lo = exp(log(est) - 1.96 * se), hi = exp(log(est) + 1.96 * se))
}
fmt_est <- function(x) {
  if (is.na(x$est)) return("-")
  sprintf("%.2f (%.2f–%.2f)", x$est, x$lo, x$hi)
}
sig_stars <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**",
                         ifelse(p < 0.05, "*", "ns")))

fisher_row <- function(var, label, section) {
  x0 <- df[[var]][df$over_60 == 0]; x1 <- df[[var]][df$over_60 == 1]
  n0 <- sum(!is.na(x0)); n1 <- sum(!is.na(x1))
  e0 <- sum(x0 == 1, na.rm = TRUE); e1 <- sum(x1 == 1, na.rm = TRUE)
  p  <- fisher.test(matrix(c(e1, n1 - e1, e0, n0 - e0), nrow = 2))$p.value
  data.frame(Section = section, Variable = label,
             lt60 = sprintf("%d/%d (%.1f%%)", e0, n0, 100 * e0 / n0),
             ge60 = sprintf("%d/%d (%.1f%%)", e1, n1, 100 * e1 / n1),
             OR_CI = fmt_est(calc_or(e1, n1, e0, n0)),
             RR_CI = fmt_est(calc_rr(e1, n1, e0, n0)),
             p = fmt_p(p), Sig = sig_stars(p),
             pct_lt = 100 * e0 / n0, pct_ge = 100 * e1 / n1,
             stringsAsFactors = FALSE)
}
mw_row <- function(var, label, section) {
  x0 <- df[[var]][df$over_60 == 0]; x1 <- df[[var]][df$over_60 == 1]
  q0 <- quantile(x0, c(.25, .5, .75), na.rm = TRUE)
  q1 <- quantile(x1, c(.25, .5, .75), na.rm = TRUE)
  p  <- wilcox.test(x0, x1, exact = FALSE)$p.value
  data.frame(Section = section, Variable = label,
             lt60 = sprintf("%g (%g–%g)", q0[2], q0[1], q0[3]),
             ge60 = sprintf("%g (%g–%g)", q1[2], q1[1], q1[3]),
             OR_CI = "-", RR_CI = "-", p = fmt_p(p), Sig = sig_stars(p),
             pct_lt = NA, pct_ge = NA, stringsAsFactors = FALSE)
}

t2 <- rbind(
  fisher_row("SSI_binary",           "SSI, any type (30d)",     "SSI (30 days)"),
  fisher_row("SSI_superficial",      "  Superficial incisional", "SSI (30 days)"),
  fisher_row("SSI_organ",            "  Organ/space",            "SSI (30 days)"),
  fisher_row("SSI_deep",             "  Deep incisional/abscess","SSI (30 days)"),
  fisher_row("intervention_for_SSI", "  Intervention for SSI",   "SSI (30 days)"),

  fisher_row("UTI_30_day",               "Urinary tract infection (30d)",            "30-day complications"),
  fisher_row("Ischaemic_CV_even_30_day", "Ischaemic cardiovascular event (30d)",     "30-day complications"),
  fisher_row("DVT_PE_30_day",            "DVT / pulmonary embolism (30d)",           "30-day complications"),
  fisher_row("LRTI_HAP_30_day",          "LRTI / hospital-acquired pneumonia (30d)", "30-day complications"),
  fisher_row("reintervention_30_days",   "Surgical reintervention (30d)",            "30-day complications"),
  fisher_row("readmission_30_day",       "Unplanned readmission (30d)",              "30-day complications"),

  fisher_row("CD_0",     "Grade 0: no complication",              "Clavien-Dindo"),
  fisher_row("CD_1",     "Grade 1: minor, no treatment",          "Clavien-Dindo"),
  fisher_row("CD_2",     "Grade 2: pharmacological treatment",    "Clavien-Dindo"),
  fisher_row("CD_3a",    "Grade 3a: intervention without GA",     "Clavien-Dindo"),
  fisher_row("CD_3b",    "Grade 3b: intervention under GA",       "Clavien-Dindo"),
  fisher_row("CD_3ab",   "Grade 3 (3a+3b): any intervention",     "Clavien-Dindo"),
  fisher_row("CD_4a",    "Grade 4a: single organ failure / ICU",  "Clavien-Dindo"),
  fisher_row("CD_4ab",   "Grade 4 (4a+4b): organ failure / ICU",  "Clavien-Dindo"),
  fisher_row("CD_3plus", "Any grade ≥3",                     "Clavien-Dindo"),

  mw_row("cr_1m",    "Creatinine 1 month (μmol/L)",   "Renal function"),
  mw_row("cr_3m",    "Creatinine 3 months (μmol/L)",  "Renal function"),
  mw_row("cr_6m",    "Creatinine 6 months (μmol/L)",  "Renal function"),
  mw_row("cr_12m",   "Creatinine 12 months (μmol/L)", "Renal function"),
  mw_row("egfr_1m",  "eGFR 1 month (mL/min/1.73m²)",   "Renal function"),
  mw_row("egfr_3m",  "eGFR 3 months (mL/min/1.73m²)",  "Renal function"),
  mw_row("egfr_6m",  "eGFR 6 months (mL/min/1.73m²)",  "Renal function"),
  mw_row("egfr_12m", "eGFR 12 months (mL/min/1.73m²)", "Renal function"),
  mw_row("LOS_days", "Length of stay (days)",                "Renal function"),

  fisher_row("DGF",                "Delayed graft function",         "Primary outcomes"),
  fisher_row("PNF",                "Primary non-function",           "Primary outcomes"),
  fisher_row("ITU_admission",      "ICU admission",                  "Primary outcomes"),
  fisher_row("graft_failure",      "Graft failure",                  "Primary outcomes"),
  fisher_row("recipient_deceased", "Patient death",                  "Primary outcomes"),
  fisher_row("fall_12m",           "Falls (12 months)",              "Primary outcomes"),
  fisher_row("malignancy_12mo",    "De novo malignancy (12 months)", "Primary outcomes"),

  fisher_row("de_novoDSA",   "De novo DSA",                          "Immunology and virology"),
  fisher_row("rejection",    "Rejection on biopsy",                  "Immunology and virology"),
  fisher_row("BPAR",         "Biopsy-proven acute rejection (BPAR)", "Immunology and virology"),
  fisher_row("CMV_viraemia", "CMV viraemia",                         "Immunology and virology"),
  fisher_row("BK_viraemia",  "BK viraemia",                          "Immunology and virology"),
  fisher_row("EBV_viraemia", "EBV viraemia",                         "Immunology and virology")
)

print(t2[t2$Sig != "ns", c("Variable", "lt60", "ge60", "OR_CI", "RR_CI", "p")],
      row.names = FALSE)

sig2 <- which(t2$Sig != "ns")
ft2 <- flextable(t2[, c("Section", "Variable", "lt60", "ge60", "OR_CI", "RR_CI", "p", "Sig")]) %>%
  set_header_labels(Section = "Section", Variable = "Variable", lt60 = lab_lt,
                    ge60 = lab_ge, OR_CI = "OR (95% CI)", RR_CI = "RR (95% CI)",
                    p = "p-value", Sig = "Sig.") %>%
  merge_v(j = "Section") %>%
  valign(j = "Section", valign = "top") %>%
  bold(part = "header") %>%
  align(j = c("lt60", "ge60", "OR_CI", "RR_CI", "p", "Sig"), align = "center", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  font(fontname = "Calibri", part = "all") %>%
  border_inner_h(part = "body") %>%
  add_footer_lines(paste0(
    "Categorical variables: n/N (%); continuous variables: median (IQR). ",
    "OR, odds ratio; RR, relative risk (age ≥60 vs <60). ",
    "- = not estimable (zero events in one group) or continuous variable. ",
    "Fisher exact test for categorical and Mann-Whitney U for continuous variables. ",
    "Highlighted rows: p<0.05.")) %>%
  fontsize(size = 9, part = "footer") %>%
  italic(part = "footer") %>%
  autofit()
if (length(sig2)) {
  ft2 <- ft2 %>% bg(i = sig2, bg = "#FFF2CC") %>%
    bold(i = sig2, j = c("p", "Sig")) %>%
    color(i = sig2, j = c("p", "Sig"), color = "#C00000")
}

read_docx() %>%
  body_add_par("Table 2. Clinical outcomes, age <60 vs ≥60", style = "heading 1") %>%
  body_add_flextable(ft2) %>%
  print(target = op("table2_outcomes.docx"))

# ---- Survival ---------------------------------------------------------------
km_plot <- function(fit, ylab, title) {
  ggsurvplot(fit, data = df, palette = unname(COLS), risk.table = TRUE,
             conf.int = TRUE, pval = TRUE, pval.method = TRUE,
             legend.labs = c("Age <60", "Age >=60"), legend.title = "",
             xlab = "Time from transplant (days)", ylab = ylab, title = title,
             risk.table.height = 0.28, ggtheme = theme_bw(base_size = 13))
}
km_pt <- survfit(Surv(survival_days, recipient_deceased) ~ Group, data = df)
km_gr <- survfit(Surv(life_of_tx_days, graft_failure) ~ Group, data = df)
ggsave(op("km_patient_survival.png"),
       print(km_plot(km_pt, "Patient survival probability", "Patient survival")),
       width = 9, height = 6.5, dpi = 300)
ggsave(op("km_graft_survival.png"),
       print(km_plot(km_gr, "Graft survival probability", "Graft survival")),
       width = 9, height = 6.5, dpi = 300)

# ---- Figures ----------------------------------------------------------------
bar_by_group <- function(dat, x, title, angle = 0) {
  counts <- dat %>% filter(!is.na(.data[[x]])) %>%
    count(Group, level = .data[[x]]) %>%
    group_by(Group) %>% mutate(pct = n / sum(n) * 100)
  ggplot(counts, aes(x = level, y = pct, fill = Group)) +
    geom_col(position = "dodge", width = 0.65, colour = "white") +
    scale_fill_manual(values = COLS) +
    labs(title = title, x = "Clavien-Dindo grade", y = "Proportion (%)", fill = "") +
    theme_bw(base_size = 13) +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = angle, hjust = if (angle) 1 else 0.5))
}

ggsave(op("clavien_dindo.png"), bar_by_group(df, "CD_grade", "Clavien-Dindo grade"),
       width = 8, height = 5, dpi = 300)

df$CD_combined <- factor(case_when(
  df$post_op_clavien_dindo_grade == "0" ~ "Grade 0",
  df$post_op_clavien_dindo_grade == "1" ~ "Grade 1",
  df$post_op_clavien_dindo_grade == "2" ~ "Grade 2",
  df$post_op_clavien_dindo_grade %in% c("3a", "3b") ~ "Grade 3 (3a+3b)",
  df$post_op_clavien_dindo_grade %in% c("4a", "4b") ~ "Grade 4 (4a+4b)"
), levels = c("Grade 0", "Grade 1", "Grade 2", "Grade 3 (3a+3b)", "Grade 4 (4a+4b)"))
ggsave(op("clavien_dindo_combined.png"),
       bar_by_group(df, "CD_combined", "Clavien-Dindo grade (combined)", angle = 15),
       width = 8, height = 5, dpi = 300)

egfr_long <- data.frame(
  Group     = rep(df$Group, 4),
  Timepoint = factor(rep(c("1 month", "3 months", "6 months", "12 months"), each = n_all),
                     levels = c("1 month", "3 months", "6 months", "12 months")),
  eGFR      = c(df$egfr_1m, df$egfr_3m, df$egfr_6m, df$egfr_12m)
)
p_egfr <- ggplot(egfr_long, aes(x = Timepoint, y = eGFR, fill = Group)) +
  geom_boxplot(outlier.alpha = 0.5, width = 0.6) +
  scale_fill_manual(values = COLS) +
  labs(title = "eGFR after transplantation", x = "",
       y = "eGFR (mL/min/1.73m²)", fill = "") +
  theme_bw(base_size = 13) + theme(legend.position = "top")
ggsave(op("egfr_boxplot.png"), p_egfr, width = 8, height = 5, dpi = 300)

outcome_bars <- c(
  SSI = "SSI, any type (30d)", UTI = "Urinary tract infection (30d)",
  `Ischaemic CV` = "Ischaemic cardiovascular event (30d)",
  `DVT/PE` = "DVT / pulmonary embolism (30d)",
  `LRTI/HAP` = "LRTI / hospital-acquired pneumonia (30d)",
  Readmission = "Unplanned readmission (30d)", DGF = "Delayed graft function",
  PNF = "Primary non-function", ICU = "ICU admission", Death = "Patient death",
  Falls = "Falls (12 months)", Malignancy = "De novo malignancy (12 months)",
  `EBV viraemia` = "EBV viraemia"
)
ob <- t2[match(outcome_bars, t2$Variable), c("pct_lt", "pct_ge")]
outcomes_long <- data.frame(
  Variable = factor(rep(names(outcome_bars), 2), levels = rev(names(outcome_bars))),
  Group    = factor(rep(c("<60", ">=60"), each = length(outcome_bars)), levels = c("<60", ">=60")),
  pct      = c(ob$pct_lt, ob$pct_ge)
)
p_out <- ggplot(outcomes_long, aes(x = Variable, y = pct, fill = Group)) +
  geom_col(position = "dodge", width = 0.65, colour = "white") +
  coord_flip() +
  scale_fill_manual(values = COLS) +
  labs(title = "Incidence of clinical outcomes", x = "", y = "Incidence (%)", fill = "") +
  theme_bw(base_size = 13) + theme(legend.position = "top")
ggsave(op("outcomes_barplot.png"), p_out, width = 9, height = 6, dpi = 300)

cat("Outputs written to", out, "\n")
