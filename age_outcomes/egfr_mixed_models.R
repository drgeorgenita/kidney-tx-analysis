# eGFR after kidney transplantation: linear mixed-effects models on repeated
# measurements at 1, 3, 6 and 12 months.

# Run for the full cohort or for one recipient age group:
#   Rscript age_outcomes/egfr_mixed_models.R            # full cohort
#   Rscript age_outcomes/egfr_mixed_models.R 70-79      # recipients aged 70-79

source("R/helpers.R")
install_missing(c("readxl", "janitor", "lme4", "lmerTest", "broom.mixed", "gtsummary",
                  "officer", "rvg", "flextable", "dplyr", "tidyr", "ggdist",
                  "beeswarm", "performance", "car"))
library(readxl); library(lme4); library(lmerTest); library(broom.mixed)
library(gtsummary); library(officer); library(rvg); library(flextable)
library(dplyr); library(tidyr); library(ggdist); library(performance)

theme_gtsummary_compact()
theme_gtsummary_journal()
set_flextable_defaults(font.family = "Arial")

args <- commandArgs(trailingOnly = TRUE)
SUBGROUP <- if (length(args)) args[1] else "all"
is_sub <- SUBGROUP != "all"

in_file <- "data/cohort_long.xlsx"
out <- file.path("output", if (is_sub) paste0("egfr_age_", SUBGROUP) else "egfr")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

prs <- read_pptx()
slide_n <- 0

egfr_lab <- "eGFR (mL/min/1.73m²)"

# ---- Model summaries --------------------------------------------------------
tidy_fixed <- function(x, ...) {
  out <- broom.mixed::tidy(x, ...)
  out[out$effect == "fixed", ]
}

# Variance components, ICC and Nakagawa R2 (falls back to an intercept-only
# approximation if performance cannot compute them)
re_components <- function(model) {
  vc <- as.data.frame(VarCorr(model))
  is_id <- tolower(vc$grp) == "id"
  tau00 <- vc$vcov[is_id & vc$var1 == "(Intercept)" & is.na(vc$var2)]
  tau00 <- if (length(tau00)) tau00[1] else vc$vcov[1]
  sig2 <- sigma(model)^2
  slope_idx <- which(is_id & vc$var1 != "(Intercept)" & is.na(vc$var2))
  has_slope <- length(slope_idx) > 0
  corr_idx <- which(is_id & !is.na(vc$var2))

  icc <- tryCatch({
    v <- as.numeric(performance::icc(model)$ICC_adjusted)
    if (length(v) == 1 && is.finite(v)) v else tau00 / (tau00 + sig2)
  }, error = function(e) tau00 / (tau00 + sig2))

  r2 <- tryCatch(performance::r2_nakagawa(model), error = function(e) NULL)
  r2m <- suppressWarnings(as.numeric(r2$R2_marginal))[1]
  r2c <- suppressWarnings(as.numeric(r2$R2_conditional))[1]
  if (!is.finite(r2m) || !is.finite(r2c)) {
    fev <- var(as.vector(model.matrix(model) %*% fixef(model)))
    tot <- fev + tau00 + sig2
    r2m <- fev / tot
    r2c <- (fev + tau00) / tot
  }
  list(sig2 = sig2, tau00 = tau00,
       tau11 = if (has_slope) vc$vcov[slope_idx[1]] else NA_real_,
       slope_nm = if (has_slope) vc$var1[slope_idx[1]] else NA_character_,
       rho01 = if (length(corr_idx)) vc$sdcor[corr_idx[1]] else NA_real_,
       icc = icc, n_id = n_distinct(model@frame$id), n_obs = nobs(model),
       r2m = r2m, r2c = r2c, has_slope = has_slope,
       aic = AIC(model), bic = BIC(model))
}

# Fixed and random effects in one table (layout as sjPlot::tab_model)
model_ft <- function(model, outcome_label = "eGFR") {
  tryCatch({
    fe <- broom.mixed::tidy(model, conf.int = TRUE, effects = "fixed")
    fe_rows <- data.frame(
      Predictors = fe$term,
      Estimates  = ifelse(is.na(fe$estimate), "", sprintf("%.2f", fe$estimate)),
      CI = ifelse(is.na(fe$conf.low), "", sprintf("%.2f to %.2f", fe$conf.low, fe$conf.high)),
      p  = fmt_p(fe$p.value),
      is_sig = !is.na(fe$p.value) & fe$p.value < 0.05,
      stringsAsFactors = FALSE)

    cc <- re_components(model)
    re_lab <- c("σ²", "τ₀₀ id")
    re_val <- c(sprintf("%.2f", cc$sig2), sprintf("%.2f", cc$tau00))
    if (cc$has_slope) {
      re_lab <- c(re_lab, sprintf("τ₁₁ %s", cc$slope_nm), "ρ₀₁")
      re_val <- c(re_val, sprintf("%.2f", cc$tau11), sprintf("%.2f", cc$rho01))
    }
    re_lab <- c(re_lab, "ICC", "N id", "Observations", "Marginal R² / Conditional R²")
    re_val <- c(re_val, sprintf("%.2f", cc$icc), cc$n_id, cc$n_obs,
                sprintf("%.3f / %.3f", cc$r2m, cc$r2c))

    d_all <- rbind(
      fe_rows,
      data.frame(Predictors = "Random Effects", Estimates = "", CI = "", p = "", is_sig = FALSE),
      data.frame(Predictors = re_lab, Estimates = re_val, CI = "", p = "", is_sig = FALSE)
    )
    nfe <- nrow(fe_rows)

    flextable(d_all[, c("Predictors", "Estimates", "CI", "p")]) |>
      add_header_row(values = c("", outcome_label, "", ""), colwidths = c(1, 1, 1, 1)) |>
      align(i = 1, part = "header", align = "center") |>
      bold(part = "header") |>
      italic(i = 2, part = "header") |>
      bold(i = nfe + 1, j = 1, part = "body") |>
      bg(i = nfe + 1, bg = "#F0F0F0", part = "body") |>
      bold(i = which(d_all$is_sig), j = "p", part = "body") |>
      border_remove() |>
      hline_top(part = "header", border = fp_border(width = 1.5)) |>
      hline(i = 2, part = "header", border = fp_border(width = 0.5)) |>
      hline_bottom(part = "header", border = fp_border(width = 1)) |>
      hline(i = nfe, border = fp_border(width = 0.5, style = "dashed"), part = "body") |>
      hline_bottom(part = "body", border = fp_border(width = 1)) |>
      fontsize(size = 9, part = "body") |>
      fontsize(size = 10, part = "header") |>
      width(j = 1, width = 2.2) |> width(j = 2, width = 1.0) |>
      width(j = 3, width = 1.6) |> width(j = 4, width = 0.7)
  }, error = function(e) { cat("  model_ft failed:", conditionMessage(e), "\n"); NULL })
}

# Random intercept vs random intercept + slope, side by side
compare_ri_rs_ft <- function(m_ri, m_rs, name_ri = "Random intercept",
                             name_rs = "Random intercept + slope") {
  tryCatch({
    f_ri <- broom.mixed::tidy(m_ri, conf.int = TRUE, effects = "fixed")
    f_rs <- broom.mixed::tidy(m_rs, conf.int = TRUE, effects = "fixed")
    terms <- union(f_ri$term, f_rs$term)
    gp <- function(d, col) vapply(terms, function(t) {
      i <- match(t, d$term); if (is.na(i)) NA_real_ else d[[col]][i]
    }, numeric(1))
    fmt_e <- function(e) ifelse(is.na(e), "", sprintf("%.2f", e))
    fmt_c <- function(lo, hi) ifelse(is.na(lo), "", sprintf("%.2f to %.2f", lo, hi))
    p_ri <- gp(f_ri, "p.value"); p_rs <- gp(f_rs, "p.value")

    fe <- data.frame(
      Predictors = terms,
      Est1 = fmt_e(gp(f_ri, "estimate")), CI1 = fmt_c(gp(f_ri, "conf.low"), gp(f_ri, "conf.high")),
      p1 = fmt_p(p_ri),
      Est2 = fmt_e(gp(f_rs, "estimate")), CI2 = fmt_c(gp(f_rs, "conf.low"), gp(f_rs, "conf.high")),
      p2 = fmt_p(p_rs),
      sig1 = !is.na(p_ri) & p_ri < 0.05, sig2 = !is.na(p_rs) & p_rs < 0.05,
      stringsAsFactors = FALSE)

    reval <- function(cc) c(
      sprintf("%.2f", cc$sig2), sprintf("%.2f", cc$tau00),
      if (cc$has_slope) sprintf("%.2f", cc$tau11) else "-",
      if (cc$has_slope && is.finite(cc$rho01)) sprintf("%.2f", cc$rho01) else "-",
      sprintf("%.2f", cc$icc), cc$n_id, cc$n_obs,
      sprintf("%.3f", cc$r2m), sprintf("%.3f", cc$r2c),
      sprintf("%.1f", cc$aic), sprintf("%.1f", cc$bic))
    re_lab <- c("σ²", "τ₀₀ id", "τ₁₁ (slope)",
                "ρ₀₁", "ICC", "N id", "Observations",
                "Marginal R²", "Conditional R²", "AIC", "BIC")
    re <- data.frame(Predictors = re_lab,
                     Est1 = reval(re_components(m_ri)), CI1 = "", p1 = "",
                     Est2 = reval(re_components(m_rs)), CI2 = "", p2 = "",
                     sig1 = FALSE, sig2 = FALSE)
    hdr <- data.frame(Predictors = "Random Effects / Fit", Est1 = "", CI1 = "", p1 = "",
                      Est2 = "", CI2 = "", p2 = "", sig1 = FALSE, sig2 = FALSE)
    d_all <- rbind(fe, hdr, re)
    nfe <- nrow(fe)

    ft <- flextable(d_all[, c("Predictors", "Est1", "CI1", "p1", "Est2", "CI2", "p2")]) |>
      set_header_labels(Est1 = "Estimate", CI1 = "95% CI", p1 = "p",
                        Est2 = "Estimate", CI2 = "95% CI", p2 = "p") |>
      add_header_row(values = c("", name_ri, name_rs), colwidths = c(1, 3, 3)) |>
      align(i = 1, part = "header", align = "center") |>
      bold(part = "header") |>
      bold(i = nfe + 1, j = 1, part = "body") |>
      bg(i = nfe + 1, bg = "#F0F0F0", part = "body") |>
      border_remove() |>
      hline_top(part = "header", border = fp_border(width = 1.5)) |>
      hline_bottom(part = "header", border = fp_border(width = 1)) |>
      hline(i = nfe, border = fp_border(width = 0.5, style = "dashed"), part = "body") |>
      hline_bottom(part = "body", border = fp_border(width = 1)) |>
      vline(j = 4, border = fp_border(width = 0.5, color = "#CCCCCC"), part = "all") |>
      fontsize(size = 9, part = "body") |>
      fontsize(size = 10, part = "header")
    if (any(d_all$sig1)) ft <- bold(ft, i = which(d_all$sig1), j = "p1", part = "body")
    if (any(d_all$sig2)) ft <- bold(ft, i = which(d_all$sig2), j = "p2", part = "body")
    autofit(ft)
  }, error = function(e) { cat("  compare_ri_rs_ft failed:", conditionMessage(e), "\n"); NULL })
}

save_check <- function(model, name) {
  tryCatch({
    png(file.path(out, paste0("diagnostics_", name, ".png")),
        width = 1500, height = 1050, res = 120)
    print(check_model(model))
    dev.off()
  }, error = function(e) {
    try(dev.off(), silent = TRUE)
    cat(sprintf("  check_model skipped (%s): %s\n", name, conditionMessage(e)))
  })
}

# ---- Plots ------------------------------------------------------------------
box_moment <- function(data, xvar, xlabel, title, palette = NULL) {
  tryCatch({
    df <- data |> filter(!is.na(.data[[xvar]]), !is.na(Timepoint))
    if (nrow(df) == 0) return(NULL)
    lvls <- levels(factor(df[[xvar]]))
    if (is.null(palette)) {
      palette <- if (grepl("age.group", xvar)) pal_age
        else if (grepl("TxType", xvar)) pal_txtype
        else if (grepl("UKKRRI.quartile", xvar)) pal_ukkrri
        else if (grepl("UKKDR.quartlie", xvar)) pal_ukkdr
        else if (grepl("WHO", xvar)) pal_whops
        else if (length(lvls) == 2) pal_binary
        else scales::hue_pal()(length(lvls))
    }
    ggplot(df, aes(x = eGFR, y = .data[[xvar]], fill = .data[[xvar]])) +
      geom_boxplot(outlier.shape = NA, linewidth = 0.4, colour = "#555555",
                   alpha = 0.55, width = 0.6) +
      geom_swarm(colour = dot_col, alpha = dot_alpha, size = dot_size) +
      scale_fill_manual(values = setNames(rep_len(palette, length(lvls)), lvls)) +
      facet_grid(. ~ Timepoint) +
      labs(title = title, x = egfr_lab, y = xlabel) +
      theme_thesis()
  }, error = function(e) { cat(sprintf("  box_moment failed (%s)\n", title)); NULL })
}

scatter_egfr <- function(data, xvar, xlabel, title) {
  data |>
    filter(!is.na(.data[[xvar]])) |>
    ggplot(aes(x = .data[[xvar]], y = eGFR)) +
    geom_point(colour = col_scatter_pt, alpha = 0.55, size = 1.8) +
    geom_smooth(colour = col_scatter_line, method = "lm", se = FALSE, linewidth = 1.0) +
    labs(title = title, x = xlabel, y = egfr_lab) +
    theme_thesis()
}

clean_term <- function(x) {
  x <- sub("^age\\.group", "age group ", x)
  x <- sub("^age\\.at\\.tx", "recipient age", x)
  x <- sub("^donor\\.age", "donor age", x)
  x <- sub("^donor\\.HTN", "donor HTN ", x)
  x <- sub("^proteinurie", "donor proteinuria ", x)
  x <- sub("^dialysis\\.status", "dialysis status ", x)
  x <- sub("^CatTx", "transplant category ", x)
  x <- sub("^CIT\\.hours", "CIT hours", x)
  x <- sub("^UKKRRI\\.quartile", "UKKRRI quartile ", x)
  x <- sub("^UKKDR\\.quartlie", "UKKDRI quartile ", x)
  trimws(gsub("[_.]", " ", x))
}

forest_plot <- function(model, title) {
  tryCatch({
    fe <- broom.mixed::tidy(model, conf.int = TRUE, effects = "fixed") |>
      filter(term != "(Intercept)") |>
      mutate(sig = !is.na(p.value) & p.value < 0.05,
             label = clean_term(term),
             label = factor(label, levels = rev(label)))
    if (nrow(fe) == 0) return(NULL)
    ggplot(fe, aes(x = estimate, y = label, colour = sig)) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = col_ref, linewidth = 0.5) +
      geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.16, linewidth = 0.7) +
      geom_point(size = 3.1) +
      scale_colour_manual(values = c("FALSE" = col_ns, "TRUE" = col_sig),
                          labels = c("FALSE" = "p ≥ 0.05", "TRUE" = "p < 0.05"),
                          name = NULL) +
      labs(title = title, subtitle = "Estimate with 95% CI",
           x = "Coefficient (mL/min/1.73m²)", y = NULL) +
      theme_thesis() +
      theme(legend.position = "bottom", axis.text.y = element_text(hjust = 1))
  }, error = function(e) { cat("  forest plot failed:", conditionMessage(e), "\n"); NULL })
}

# Ranked patient-level random intercepts (BLUPs) with 95% intervals
re_plot <- function(model, title) {
  tryCatch({
    re_df <- broom.mixed::tidy(model, effects = "ran_vals", conf.int = TRUE)
    re_df <- re_df[re_df$term == "(Intercept)", , drop = FALSE]
    if (nrow(re_df) == 0) return(NULL)
    if (!"conf.low" %in% names(re_df)) {
      re_df$conf.low  <- re_df$estimate - 1.96 * re_df$std.error
      re_df$conf.high <- re_df$estimate + 1.96 * re_df$std.error
    }
    re_df <- re_df[order(re_df$estimate), ]
    re_df$rank  <- seq_len(nrow(re_df))
    re_df$above <- re_df$estimate >= 0
    ggplot(re_df, aes(x = rank, y = estimate, colour = above)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = col_ref, linewidth = 0.5) +
      geom_linerange(aes(ymin = conf.low, ymax = conf.high), linewidth = 0.3, alpha = 0.30) +
      geom_point(size = 0.9, alpha = 0.85) +
      scale_colour_manual(values = c("FALSE" = "#D55E00", "TRUE" = "#0072B2"), guide = "none") +
      labs(title = title, subtitle = "Patient-level random intercepts with 95% CI, ranked",
           x = "Patients (ranked)", y = "Random intercept (mL/min/1.73m²)") +
      theme_thesis() +
      theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank())
  }, error = function(e) { cat("  re_plot failed:", conditionMessage(e), "\n"); NULL })
}

# descriptive plot followed by the random intercepts of the univariable model
plot_with_re <- function(p, label, rhs, re_title, data = d) {
  add_plot(p, label)
  m <- tryCatch(lmer(as.formula(paste("eGFR ~", rhs, "+ (1 | id)")), data = data),
                error = function(e) NULL)
  if (!is.null(m)) add_plot(re_plot(m, re_title), paste(label, "random effects"))
}

# univariable random-intercept models for a set of predictors
uv_lmer <- function(vars, labels, caption, data = d) {
  vars <- vars[vars %in% names(data)]
  data |>
    dplyr::select(eGFR, id, all_of(vars)) |>
    tbl_uvregression(y = eGFR, method = lmer, formula = "{y} ~ {x} + (1 | id)",
                     tidy_fun = tidy_fixed, label = labels[vars]) |>
    modify_header(label = "**Predictor**", n_obs = "**N**",
                  estimate = "**Beta (95% CI)**", p.value = "**p-value**") |>
    bold_p(t = 0.05) |>
    modify_caption(caption)
}

# ---- Data -------------------------------------------------------------------
d <- read_excel(in_file, sheet = 1, na = c("", "NA", "N/A", "NULL")) |>
  janitor::clean_names() |>
  rename(
    age.at.tx = age_at_tx, rec.sex = rec_sex, over.60 = over_60, age.group = age_group,
    rec.BMI = rec_bmi, rec.smoking.status = rec_smoking_status,
    rec.diabetes = rec_diabetes, dialysis.status = dialysis_status,
    rec.liver.disease = rec_liver_disease,
    charlson.comorbidity.index = charlson_comorbidity_index,
    WHO.PS = who_ps, ASA = asa, TxType = tx_type, anatomie = anatomy,
    HLA.MM.level = hla_mm_level, donor.age = donor_age, donor.sex = donor_sex,
    donor.eGFR.admission = donor_e_gfr_admission,
    UKKRRI.quartile = ukkrri_quartile, UKKDR.quartlie = ukkdr_quartlie,
    ABDR.MM.sum = abdr_mm_sum, induction.IS = induction_is, imm.risk = imm_risk,
    LOS.days = los_days, CIT.min = cit_min, cr.1m = cr_1m,
    UKKRRI = ukkrri, UKKDRI = ukkdri, Moment = moment, eGFR = e_gfr,
    donor_HTN = donor_htn
  ) |>
  mutate(
    ComplexUro = as.numeric(complex_uro_path_reflux_hydron_aug_bladder),
    proteinurie_num = suppressWarnings(as.numeric(as.character(donor_proteinuria))),
    across(c(age.at.tx, LOS.days, CIT.min, donor.age, charlson.comorbidity.index,
             rec.BMI, UKKRRI, UKKDRI, donor.eGFR.admission, cr.1m), as.numeric),
    Moment = factor(as.character(Moment), levels = c("I", "II", "III", "IV"),
                    labels = c("egfr.1m", "egfr.3m", "egfr.6m", "egfr.12m")),
    Moment_num = c(1, 3, 6, 12)[as.integer(Moment)],
    Timepoint = factor(Moment, labels = c("1 month", "3 months", "6 months", "12 months")),
    eGFR = as.numeric(eGFR),
    id = factor(id),
    TxType = factor(case_when(TxType == "DCD" ~ "DCD",
                              TxType == "DBD" ~ "DBD",
                              grepl("^LD", TxType, ignore.case = TRUE) ~ "Live Donor"),
                    levels = c("DCD", "DBD", "Live Donor")),
    HLA.MM.level = factor(HLA.MM.level),
    WHO.PS = factor(WHO.PS),
    ASA = factor(ASA),
    rec.smoking.status = factor(rec.smoking.status),
    rec.diabetes = factor(rec.diabetes),
    rec.liver.disease = factor(rec.liver.disease),
    dialysis.status = factor(ifelse(dialysis.status == "pre", "Pre-emptive", "On dialysis (HD/PD)"),
                             levels = c("Pre-emptive", "On dialysis (HD/PD)")),
    proteinurie = factor(ifelse(proteinurie_num > 0, "Yes", "No"), levels = c("No", "Yes")),
    imm.risk = factor(imm.risk),
    UKKRRI.quartile = factor(UKKRRI.quartile, levels = c("R1", "R2", "R3", "R4")),
    UKKDR.quartlie = factor(case_when(
      UKKDR.quartlie %in% c("D1", "D2", "D3", "D4") ~ UKKDR.quartlie,
      grepl("live", UKKDR.quartlie, ignore.case = TRUE) ~ "live donor"),
      levels = c("live donor", "D1", "D2", "D3", "D4")),
    age.group = factor(age.group, levels = c("<50", "50-59", "60-69", "70-79")),
    anatomie = factor(case_when(
      grepl("^standard", anatomie, ignore.case = TRUE) ~ "Standard",
      grepl("not requiring", anatomie, ignore.case = TRUE) ~ "Non-standard, no back table",
      grepl("requiring back", anatomie, ignore.case = TRUE) ~ "Non-standard, back table required"),
      levels = c("Standard", "Non-standard, no back table", "Non-standard, back table required")),
    logCreat = log(cr.1m),
    CIT.10min = CIT.min / 10,
    CIT.hours = CIT.min / 60,
    etnic = factor(ifelse(grepl("white", tolower(ethnicity)) & !grepl("non", tolower(ethnicity)),
                          "White", "Non-white"), levels = c("Non-white", "White")),
    donor.HTN = factor(ifelse(donor_HTN == 1, "Yes", "No"), levels = c("No", "Yes")),
    induction.IS = factor(ifelse(grepl("Alemtuzumab", induction.IS, ignore.case = TRUE),
                                 "Alemtuzumab", as.character(induction.IS)),
                          levels = c("Basiliximab", "Alemtuzumab", "ATG"))
  ) |>
  filter(!is.na(eGFR))

if (is_sub) d <- d |> filter(age.group == SUBGROUP) |> droplevels()

n_pts <- n_distinct(d$id)
n_obs <- nrow(d)
cat(sprintf("%s: %d observations from %d patients\n",
            if (is_sub) paste("Age", SUBGROUP) else "Full cohort", n_obs, n_pts))
print(table(d$Timepoint))
print(table(d$TxType, useNA = "always"))

labels <- list(
  age.at.tx = "Recipient age (years)", age.group = "Age group", over.60 = "Age ≥60",
  rec.sex = "Sex", etnic = "Ethnicity (White vs Non-white)",
  rec.BMI = "BMI", rec.smoking.status = "Smoking status", rec.diabetes = "Diabetes",
  ComplexUro = "Complex urology", rec.liver.disease = "Liver disease",
  charlson.comorbidity.index = "Charlson CCI", ASA = "ASA grade",
  WHO.PS = "WHO performance status",
  dialysis.status = "Dialysis status (pre-emptive vs on dialysis)",
  LOS.days = "Recipient LOS (days)",
  TxType = "Transplant type", CIT.10min = "CIT (per 10 min)", CIT.hours = "CIT (hours)",
  HLA.MM.level = "HLA mismatch level", ABDR.MM.sum = "ABDR mismatch sum",
  induction.IS = "Induction immunosuppression", imm.risk = "Immunological risk",
  UKKRRI = "UKKRRI (continuous)", UKKRRI.quartile = "UKKRRI quartile",
  UKKDRI = "UKKDRI (continuous)", UKKDR.quartlie = "UKKDRI quartile",
  donor.age = "Donor age (years)", donor.sex = "Donor sex",
  donor.HTN = "Donor hypertension", proteinurie = "Donor proteinuria",
  anatomie = "Donor anatomy", CatTx = "Transplant category (live vs deceased)"
)

# ---- 1. Time ----------------------------------------------------------------
# Models compared by ML; refitted with REML for inference
fit_moment             <- lmer(eGFR ~ Moment + (1 | id), data = d, REML = FALSE)
fit_moment_cont        <- lmer(eGFR ~ Moment_num + (1 | id), data = d, REML = FALSE)
fit_moment_cont_slopes <- lmer(eGFR ~ Moment_num + (1 + Moment_num | id), data = d, REML = FALSE)
print(compare_performance(fit_moment, fit_moment_cont, fit_moment_cont_slopes))

time_delta_aic <- AIC(fit_moment_cont) - AIC(fit_moment_cont_slopes)

fit_moment             <- update(fit_moment, REML = TRUE)
fit_moment_cont        <- update(fit_moment_cont, REML = TRUE)
fit_moment_cont_slopes <- update(fit_moment_cont_slopes, REML = TRUE)

# random slopes retained only if non-singular and ML AIC improves by >= 2
time_rs_singular <- isSingular(fit_moment_cont_slopes, tol = 1e-4)
time_use_rs <- !time_rs_singular && is.finite(time_delta_aic) && time_delta_aic >= 2
fit_time <- if (time_use_rs) fit_moment_cont_slopes else fit_moment_cont
time_structure <- if (time_use_rs) {
  "random intercepts and slopes"
} else if (time_rs_singular) {
  "random intercepts (random-slopes model singular)"
} else {
  "random intercepts (random slopes did not improve AIC by ≥2)"
}
cat(sprintf("Time model: %s (delta AIC %.1f)\n", time_structure, time_delta_aic))

save_check(fit_moment, "time_factor")
save_check(fit_moment_cont, "time_continuous")
save_check(fit_moment_cont_slopes, "time_slopes")

p1 <- d |>
  filter(!is.na(Timepoint)) |>
  ggplot(aes(x = eGFR, y = Timepoint, fill = Timepoint)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.4, colour = "#555555", alpha = 0.55, width = 0.6) +
  geom_swarm(colour = dot_col, alpha = dot_alpha, size = dot_size) +
  scale_fill_manual(values = pal_age) +
  labs(title = "eGFR by timepoint", y = "Timepoint", x = egfr_lab) +
  theme_thesis()
add_plot(p1, "timepoint")
add_plot(re_plot(fit_time, "Random effects: time model"), "timepoint random effects")

# largest change between 1 and 12 months
dd_delta <- inner_join(
  d |> filter(Moment == "egfr.1m") |> dplyr::select(id, eGFR_1m = eGFR),
  d |> filter(Moment == "egfr.12m") |> dplyr::select(id, eGFR_12m = eGFR),
  by = "id") |>
  mutate(change = eGFR_12m - eGFR_1m) |>
  arrange(desc(abs(change)))
cat(sprintf("Patients with 1 and 12 month eGFR: %d\n", nrow(dd_delta)))
print(head(dd_delta, 10))

p1ext <- d |>
  filter(id %in% head(dd_delta$id, 8)) |>
  ggplot(aes(x = Moment_num, y = eGFR, colour = id)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ id) +
  scale_color_brewer(palette = "Dark2") +
  labs(title = "Patients with the largest eGFR change", x = "Time (months)", y = egfr_lab) +
  theme_thesis() +
  theme(strip.text = element_text(size = 9))
add_plot(p1ext, "extreme trajectories")

# ---- 2. Demographics --------------------------------------------------------
demo_vars <- c("age.at.tx", "age.group", "over.60", "rec.sex", "etnic")
# age group is constant within an age subgroup
if (is_sub) demo_vars <- setdiff(demo_vars, c("age.group", "over.60"))
tbl_s2 <- uv_lmer(demo_vars, labels, sprintf("Demographics | n=%d patients, %d observations",
                                             n_pts, n_obs))

plot_with_re(scatter_egfr(d, "age.at.tx", "Recipient age (years)", "eGFR vs recipient age"),
             "recipient age", "age.at.tx", "Random effects: recipient age")
if (!is_sub) {
  plot_with_re(box_moment(d, "age.group", "Age group", "eGFR by age group"),
               "age group", "age.group", "Random effects: age group")
}

# ---- 3. Comorbidities -------------------------------------------------------
comorb_vars <- c("rec.BMI", "rec.smoking.status", "rec.diabetes", "ComplexUro",
                 "rec.liver.disease", "charlson.comorbidity.index", "ASA", "WHO.PS",
                 "dialysis.status", "LOS.days")
tbl_s3 <- uv_lmer(comorb_vars, labels, sprintf("Comorbidities | n=%d patients", n_pts))

plot_with_re(box_moment(d, "WHO.PS", "WHO performance status", "eGFR by WHO performance status"),
             "WHO PS", "WHO.PS", "Random effects: WHO performance status")
plot_with_re(scatter_egfr(d, "rec.BMI", "BMI", "eGFR vs recipient BMI"),
             "BMI", "rec.BMI", "Random effects: BMI")
plot_with_re(scatter_egfr(d, "charlson.comorbidity.index", "Charlson CCI",
                          "eGFR vs Charlson Comorbidity Index"),
             "CCI", "charlson.comorbidity.index", "Random effects: Charlson CCI")
plot_with_re(box_moment(d, "dialysis.status", "Dialysis status", "eGFR by dialysis status"),
             "dialysis", "dialysis.status", "Random effects: dialysis status")
plot_with_re(scatter_egfr(d, "LOS.days", "Recipient LOS (days)", "eGFR vs recipient length of stay"),
             "recipient LOS", "LOS.days", "Random effects: recipient LOS")

# ---- 4. Transplant, ischaemia and immunology --------------------------------
tx_vars <- c("TxType", "CIT.10min", "CIT.hours", "HLA.MM.level", "ABDR.MM.sum",
             "induction.IS", "imm.risk")
tbl_s4 <- uv_lmer(tx_vars, labels, "Transplant, ischaemia and immunological factors")

plot_with_re(box_moment(d, "TxType", "Transplant type", "eGFR by transplant type"),
             "transplant type", "TxType", "Random effects: transplant type")

traj_sum <- d |>
  filter(!is.na(TxType)) |>
  group_by(TxType, Timepoint) |>
  summarise(mean_eGFR = mean(eGFR, na.rm = TRUE),
            se_eGFR = sd(eGFR, na.rm = TRUE) / sqrt(sum(!is.na(eGFR))),
            .groups = "drop")
traj_lab <- traj_sum |> group_by(TxType) |>
  slice_max(order_by = as.integer(Timepoint), n = 1) |> ungroup()
pal_tx3 <- setNames(pal_txtype, c("DCD", "DBD", "Live Donor"))

p4_traj <- ggplot(traj_sum, aes(x = Timepoint, y = mean_eGFR, group = TxType, colour = TxType)) +
  geom_ribbon(aes(ymin = mean_eGFR - se_eGFR, ymax = mean_eGFR + se_eGFR, fill = TxType),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1.4) +
  geom_point(aes(fill = TxType), shape = 21, size = 3.6, stroke = 1.3, colour = "white") +
  geom_text(data = traj_lab, aes(label = TxType), hjust = 0, nudge_x = 0.08,
            size = 3.3, fontface = "bold") +
  scale_colour_manual(values = pal_tx3) +
  scale_fill_manual(values = pal_tx3) +
  scale_x_discrete(expand = expansion(mult = c(0.04, 0.24))) +
  labs(title = "Mean eGFR by transplant type", subtitle = "Mean ± SE",
       x = "Timepoint", y = "Mean eGFR (mL/min/1.73m²)",
       caption = "DCD, donation after circulatory death; DBD, donation after brain death.") +
  theme_thesis()
add_plot(p4_traj, "transplant type trajectory")

plot_with_re(scatter_egfr(d, "CIT.min", "CIT (minutes)", "eGFR vs CIT (minutes)"),
             "CIT minutes", "CIT.min", "Random effects: CIT (minutes)")
plot_with_re(scatter_egfr(d, "CIT.hours", "CIT (hours)", "eGFR vs CIT (hours)"),
             "CIT hours", "CIT.hours", "Random effects: CIT (hours)")
plot_with_re(box_moment(d, "induction.IS", "Induction", "eGFR by induction immunosuppression"),
             "induction", "induction.IS", "Random effects: induction")

# ---- 5. UK risk indices -----------------------------------------------------
idx_vars <- c("UKKRRI", "UKKRRI.quartile", "UKKDRI", "UKKDR.quartlie")
tbl_s5 <- uv_lmer(idx_vars, labels, "UKKRRI and UKKDRI")

plot_with_re(scatter_egfr(d, "UKKRRI", "UKKRRI", "eGFR vs UKKRRI"),
             "UKKRRI", "UKKRRI", "Random effects: UKKRRI")
plot_with_re(scatter_egfr(d, "UKKDRI", "UKKDRI", "eGFR vs UKKDRI"),
             "UKKDRI", "UKKDRI", "Random effects: UKKDRI")
plot_with_re(box_moment(d, "UKKRRI.quartile", "UKKRRI quartile", "eGFR by UKKRRI quartile"),
             "UKKRRI quartile", "UKKRRI.quartile", "Random effects: UKKRRI quartile")
plot_with_re(box_moment(d, "UKKDR.quartlie", "UKKDRI quartile", "eGFR by UKKDRI quartile"),
             "UKKDRI quartile", "UKKDR.quartlie", "Random effects: UKKDRI quartile")

# ---- 6. Donor factors -------------------------------------------------------
donor_vars <- c("donor.age", "donor.sex", "donor.HTN", "proteinurie", "anatomie")
tbl_s6 <- uv_lmer(donor_vars, labels, "Donor factors")

plot_with_re(scatter_egfr(d, "donor.age", "Donor age (years)", "eGFR vs donor age"),
             "donor age", "donor.age", "Random effects: donor age")
plot_with_re(box_moment(d, "donor.HTN", "Donor hypertension", "eGFR by donor hypertension"),
             "donor HTN", "donor.HTN", "Random effects: donor hypertension")
plot_with_re(box_moment(d, "proteinurie", "Donor proteinuria", "eGFR by donor proteinuria"),
             "donor proteinuria", "proteinurie", "Random effects: donor proteinuria")
p_anat <- box_moment(d, "anatomie", "Anatomy", "eGFR by donor anatomy",
                     palette = c("#0072B2", "#56B4E9", "#E69F00"))
if (!is.null(p_anat)) p_anat <- p_anat + theme(axis.text.y = element_text(size = 7))
plot_with_re(p_anat, "anatomy", "anatomie", "Random effects: anatomy")

# ---- 7. Creatinine at 1 month -----------------------------------------------
d_cr <- d |> filter(!is.na(cr.1m))
ctrl <- lmerControl(optimizer = "bobyqa")
m_creat    <- lmer(eGFR ~ cr.1m + (1 | id), data = d_cr, control = ctrl)
m_logcreat <- lmer(eGFR ~ logCreat + (1 | id), data = d_cr, control = ctrl)

add_plot(scatter_egfr(d_cr, "cr.1m", "Creatinine (µmol/L)", "eGFR vs creatinine"),
         "creatinine")
add_plot(re_plot(m_creat, "Random effects: creatinine"), "creatinine random effects")
add_plot(scatter_egfr(d_cr, "logCreat", "log(creatinine)", "eGFR vs log(creatinine)"),
         "log creatinine")
add_plot(re_plot(m_logcreat, "Random effects: log(creatinine)"), "log creatinine random effects")

tbl_s7 <- tbl_merge(
  list(tbl_regression(m_creat, tidy_fun = tidy_fixed, intercept = TRUE,
                      label = list(cr.1m ~ "Creatinine at 1 month")),
       tbl_regression(m_logcreat, tidy_fun = tidy_fixed, intercept = TRUE,
                      label = list(logCreat ~ "log(creatinine at 1 month)"))),
  tab_spanner = c("**Creatinine**", "**log(creatinine)**")) |>
  modify_caption(sprintf("Creatinine models | n=%d patients", n_distinct(d_cr$id)))

# ---- Variance explained by each predictor -----------------------------------
marg_r2 <- function(xvar) {
  tryCatch({
    m <- lmer(as.formula(paste0("eGFR ~ ", xvar, " + (1 | id)")),
              data = filter(d, !is.na(.data[[xvar]])))
    100 * as.numeric(r2_nakagawa(m)$R2_marginal)
  }, error = function(e) NA_real_)
}

r2_df <- data.frame(
  var = c("UKKDR.quartlie", "CIT.hours", "TxType", "donor.age", "donor.HTN", "proteinurie",
          "UKKRRI.quartile", "age.group", "dialysis.status", "charlson.comorbidity.index",
          "WHO.PS", "rec.BMI", "LOS.days", "Moment_num"),
  label = c("UKKDRI quartile", "Cold ischaemia time", "Transplant type", "Donor age",
            "Donor hypertension", "Donor proteinuria", "UKKRRI quartile",
            "Recipient age group", "Dialysis status", "Charlson Comorbidity Index",
            "WHO performance status", "Body mass index", "Recipient LOS",
            "Measurement timepoint"),
  group = c(rep("Donor / logistic factors", 6), rep("Recipient factors", 7), "Time"),
  stringsAsFactors = FALSE)
r2_df$r2 <- vapply(r2_df$var, marg_r2, numeric(1))
r2_df <- r2_df[!is.na(r2_df$r2), ]
r2_df <- r2_df[order(r2_df$r2), ]
r2_df$label <- factor(r2_df$label, levels = r2_df$label)
r2_df$group <- factor(r2_df$group, levels = c("Donor / logistic factors", "Recipient factors", "Time"))

p_r2 <- ggplot(r2_df, aes(x = r2, y = label, fill = group)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = sprintf("%.1f%%", r2)), hjust = -0.15, size = 3.2, colour = "grey25") +
  scale_fill_manual(values = c("Donor / logistic factors" = "#1f3b63",
                               "Recipient factors" = "#d9912a", "Time" = "#9fb6cf"),
                    name = NULL) +
  scale_x_continuous("Variance in eGFR explained, univariable marginal R² (%)",
                     limits = c(0, max(20, ceiling(max(r2_df$r2) / 2.5) * 2.5 + 2.5)),
                     breaks = seq(0, 100, 2.5), expand = expansion(mult = c(0, 0.04))) +
  labs(y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.x = element_line(colour = "grey80"), legend.position = "bottom")
add_plot(p_r2, "variance explained")

# ---- Multivariable models ---------------------------------------------------
# fitted by ML for AIC comparison, refitted by REML for estimates
fit_pair <- function(f, data) {
  ml <- tryCatch(lmer(f, data = data, control = ctrl, REML = FALSE),
                 error = function(e) { cat("  fit failed:", conditionMessage(e), "\n"); NULL })
  list(ml = ml, reml = if (is.null(ml)) NULL else update(ml, REML = TRUE))
}

mv_table <- function(uni, model, lab, caption) {
  multi <- tbl_regression(model, tidy_fun = tidy_fixed, intercept = TRUE, label = lab) |>
    bold_p(t = 0.05)
  tbl_merge(list(uni, multi), tab_spanner = c("**Univariable**", "**Multivariable**")) |>
    modify_caption(caption)
}

mv_outputs <- function(m, name, title, uni, lab, caption_fmt, data) {
  if (is.null(m$reml)) return(NULL)
  save_check(m$reml, name)
  add_plot(forest_plot(m$reml, title), paste(name, "forest"))
  add_plot(re_plot(m$reml, paste("Random effects:", name)), paste(name, "random effects"))
  mv_table(uni, m$reml, lab,
           sprintf(caption_fmt, nobs(m$reml), n_distinct(data$id)))
}

# Model A: recipient and donor risk indices
d_ukk <- d |> filter(!is.na(UKKRRI.quartile), !is.na(UKKDR.quartlie))
lab_A <- list(UKKRRI.quartile ~ "UKKRRI quartile", UKKDR.quartlie ~ "UKKDRI quartile")
uni_A <- uv_lmer(c("UKKRRI.quartile", "UKKDR.quartlie"), labels, "", data = d_ukk)

mA1 <- fit_pair(eGFR ~ UKKRRI.quartile + UKKDR.quartlie + (1 | id), d_ukk)
mA2 <- fit_pair(eGFR ~ UKKRRI.quartile + UKKDR.quartlie + Moment_num + (1 + Moment_num | id), d_ukk)

tbl_mA1 <- mv_outputs(mA1, "Model A1", "Model A1: UKKRRI + UKKDRI quartile", uni_A, lab_A,
                      "Model A1 (random intercept) | n=%d observations, %d patients | Ref: R1, live donor",
                      d_ukk)
tbl_mA2 <- mv_outputs(mA2, "Model A2", "Model A2: UKKRRI + UKKDRI quartile + time", uni_A,
                      c(lab_A, list(Moment_num ~ "Time (months)")),
                      "Model A2 (random slopes) | n=%d observations, %d patients | Ref: R1, live donor",
                      d_ukk)

# Model B: clinical predictors
d_mB <- d |>
  mutate(CatTx = factor(ifelse(TxType == "Live Donor", "Live", "Deceased"),
                        levels = c("Live", "Deceased"))) |>
  filter(!is.na(donor.age), !is.na(donor.HTN), !is.na(CIT.hours),
         !is.na(proteinurie), !is.na(dialysis.status))
lab_B <- list(age.at.tx ~ "Recipient age (years)", donor.age ~ "Donor age (years)",
              donor.HTN ~ "Donor hypertension", proteinurie ~ "Donor proteinuria",
              dialysis.status ~ "Dialysis status", CatTx ~ "Transplant category",
              `CatTx:CIT.hours` ~ "Transplant category × CIT (hours)")
uni_B <- uv_lmer(c("age.at.tx", "donor.age", "donor.HTN", "proteinurie", "dialysis.status", "CatTx"),
                 labels, "", data = d_mB)

mB1 <- fit_pair(eGFR ~ age.at.tx + donor.age + donor.HTN + proteinurie +
                  dialysis.status + CatTx * CIT.hours + (1 | id), d_mB)
mB2 <- fit_pair(eGFR ~ age.at.tx + Moment_num + donor.age + donor.HTN + proteinurie +
                  dialysis.status + CatTx * CIT.hours + (1 + Moment_num | id), d_mB)

for (m in list(mB1$reml, mB2$reml)) {
  if (!is.null(m)) tryCatch(print(car::vif(m)), error = function(e) NULL)
}

tbl_mB1 <- mv_outputs(mB1, "Model B1", "Model B1: clinical predictors", uni_B, lab_B,
                      "Model B1 (random intercept) | n=%d observations, %d patients",
                      d_mB)
tbl_mB2 <- mv_outputs(mB2, "Model B2", "Model B2: clinical predictors + time", uni_B,
                      c(lab_B, list(Moment_num ~ "Time (months)")),
                      "Model B2 (random slopes) | n=%d observations, %d patients",
                      d_mB)

aic_table <- function(models, descr) {
  rows <- Filter(Negate(is.null), Map(function(m, nm, ds) {
    if (is.null(m)) return(NULL)
    data.frame(Model = nm, Predictors = ds, N_obs = nobs(m),
               N_patients = n_distinct(m@frame$id), AIC = round(AIC(m), 1))
  }, models, names(models), descr))
  if (!length(rows)) return(NULL)
  df <- do.call(rbind, rows)
  print(df)
  best <- which.min(df$AIC)
  df$Best <- ifelse(seq_len(nrow(df)) == best, "Lowest AIC", "")
  flextable(df) |>
    set_header_labels(N_obs = "N obs", N_patients = "N patients", Best = "") |>
    bold(part = "header") |>
    bold(i = best, part = "body") |>
    bg(i = best, bg = "#E8F5E9", part = "body") |>
    add_footer_lines("Fitted by maximum likelihood on the same observations. ΔAIC > 2 indicates a meaningfully better fit.") |>
    italic(part = "footer") |>
    fontsize(size = 8, part = "footer") |>
    autofit()
}
ft_aic_A <- aic_table(list(`A1 (RI)` = mA1$ml, `A2 (RS)` = mA2$ml),
                      c("UKKRRI + UKKDRI quartile", "UKKRRI + UKKDRI quartile + time"))
ft_aic_B <- aic_table(list(`B1 (RI)` = mB1$ml, `B2 (RS)` = mB2$ml),
                      c("Clinical predictors", "Clinical predictors + time"))

print(prs, target = file.path(out, "egfr_figures.pptx"))

# ---- Word document ----------------------------------------------------------
to_ft <- function(tbl) {
  tbl |> as_flex_table() |>
    fontsize(size = 9, part = "body") |>
    fontsize(size = 10, part = "header") |>
    autofit()
}

# multivariable tables: rows coloured by the multivariable p-value
to_ft_mv <- function(tbl) {
  ft <- to_ft(tbl)
  bd <- ft$body$dataset
  p_cols <- grep("p\\.value|p_value", names(bd), value = TRUE)
  if (length(p_cols)) {
    p_raw <- trimws(as.character(bd[[tail(p_cols, 1)]]))
    p_num <- suppressWarnings(as.numeric(p_raw))
    p_num[grepl("^<", p_raw)] <- 0
    sig  <- which(!is.na(p_num) & p_num < 0.05)
    nsig <- which(!is.na(p_num) & p_num >= 0.05)
    if (length(sig))  ft <- color(ft, i = sig,  color = "#C0392B", part = "body")
    if (length(nsig)) ft <- color(ft, i = nsig, color = "#2E86AB", part = "body")
  }
  ft |>
    add_footer_lines("Red, p < 0.05; blue, p ≥ 0.05 (multivariable p-value)") |>
    italic(part = "footer") |>
    fontsize(size = 8, part = "footer") |>
    color(color = "grey40", part = "footer")
}

re_table <- function(model) {
  vc <- as.data.frame(VarCorr(model))
  tau00 <- vc$vcov[tolower(vc$grp) == "id" & vc$var1 == "(Intercept)" & is.na(vc$var2)][1]
  sig2 <- sigma(model)^2
  fe_var <- var(as.vector(model.matrix(model) %*% fixef(model)))
  total <- fe_var + tau00 + sig2
  flextable(data.frame(
    Statistic = c("σ²", "τ₀₀ id", "ICC", "N patients", "Observations",
                  "Marginal R² / Conditional R²"),
    Value = c(round(sig2, 2), round(tau00, 2), round(tau00 / (tau00 + sig2), 2),
              n_distinct(model@frame$id), nobs(model),
              sprintf("%.3f / %.3f", fe_var / total, (fe_var + tau00) / total)))) |>
    set_header_labels(Statistic = "Random effects", Value = "") |>
    bold(part = "header") |>
    bg(part = "header", bg = "#D5E8F0") |>
    fontsize(size = 9, part = "body") |>
    width(j = 1, width = 2.5) |> width(j = 2, width = 1.2)
}

br <- function(doc) body_add_par(doc, " ", style = "Normal")

# full table for each univariable model with at least one term at p<0.05
add_sig_models <- function(doc, predictors) {
  for (pred in predictors) {
    d_sub <- d |> filter(!is.na(.data[[pred]]))
    if (nrow(d_sub) < 20) next
    m <- tryCatch(lmer(as.formula(paste0("eGFR ~ `", pred, "` + (1 | id)")),
                       data = d_sub, control = ctrl), error = function(e) NULL)
    if (is.null(m)) next
    co <- as.data.frame(summary(m)$coefficients)
    co <- co[rownames(co) != "(Intercept)", , drop = FALSE]
    if (nrow(co) == 0 || all(co[, ncol(co)] >= 0.05, na.rm = TRUE)) next
    ft <- model_ft(m)
    if (!is.null(ft)) doc <- doc |> body_add_flextable(ft) |> br()
  }
  doc
}

add_section <- function(doc, heading, tbl, predictors) {
  doc |>
    body_add_break() |>
    body_add_par(heading, style = "heading 2") |>
    body_add_flextable(to_ft(tbl)) |>
    br() |>
    add_sig_models(predictors)
}

add_model <- function(doc, heading, note, m, tbl) {
  doc <- doc |> body_add_break() |>
    body_add_par(heading, style = "heading 2") |>
    body_add_par(note, style = "Normal")
  if (!is.null(m))   doc <- doc |> body_add_flextable(re_table(m)) |> br()
  if (!is.null(tbl)) doc <- doc |> body_add_flextable(to_ft_mv(tbl))
  doc
}

title <- if (is_sub) {
  sprintf("eGFR analysis: recipients aged %s (n=%d)", SUBGROUP, n_pts)
} else {
  sprintf("eGFR analysis: full cohort (n=%d)", n_pts)
}

doc <- read_docx() |>
  body_add_par(title, style = "heading 1") |>
  body_add_par(sprintf(paste("%d observations from %d patients. Time model: %s.",
                             "Patient-level predictors: random-intercept models."),
                       n_obs, n_pts, time_structure), style = "Normal") |>
  br() |>
  body_add_par("eGFR by time since transplant", style = "heading 2") |>
  body_add_par(paste("Time in months as a continuous variable. Random slopes were retained",
                     "only if the model was non-singular and improved ML AIC by at least 2."),
               style = "Normal") |>
  body_add_flextable(model_ft(fit_time)) |>
  br() |>
  body_add_par("Random intercept vs random slopes (time model)", style = "heading 2") |>
  body_add_par(paste("Identical fixed effects; the models differ only in the random-effects",
                     "structure. Random slopes are considered for time only, as it is the",
                     "only predictor that varies within patients."),
               style = "Normal")
cmp_ft <- compare_ri_rs_ft(fit_moment_cont, fit_moment_cont_slopes)
if (!is.null(cmp_ft)) doc <- doc |> body_add_flextable(cmp_ft) |> br()

doc <- add_section(doc, "Demographics", tbl_s2, demo_vars)
doc <- add_section(doc, "Comorbidities", tbl_s3, comorb_vars)
doc <- add_section(doc, "Transplant, ischaemia and immunological factors", tbl_s4,
                   setdiff(tx_vars, "CIT.10min"))
doc <- add_section(doc, "UKKRRI and UKKDRI", tbl_s5, idx_vars)
doc <- add_section(doc, "Donor factors", tbl_s6, c(donor_vars, "donor.eGFR.admission"))

doc <- doc |> body_add_break() |>
  body_add_par("Creatinine at 1 month", style = "heading 2") |>
  body_add_flextable(to_ft_mv(tbl_s7)) |>
  br()

doc <- add_model(doc, "Model A1: UKKRRI + UKKDRI quartile",
                 "Reference: UKKRRI R1, UKKDRI live donor. Random intercept.",
                 mA1$reml, tbl_mA1)
doc <- add_model(doc, "Model A2: UKKRRI + UKKDRI quartile + time",
                 "Reference: UKKRRI R1, UKKDRI live donor. Random intercept and slope for time.",
                 mA2$reml, tbl_mA2)
if (!is.null(ft_aic_A)) {
  doc <- doc |> body_add_break() |>
    body_add_par("AIC: models A1 and A2", style = "heading 2") |>
    body_add_flextable(ft_aic_A)
}
doc <- add_model(doc, "Model B1: clinical predictors",
                 paste("Reference: live donor, no donor hypertension, no donor proteinuria,",
                       "pre-emptive. Recipient age continuous. Random intercept. Transplant",
                       "category and CIT are strongly correlated; interpret the interaction",
                       "with caution."),
                 mB1$reml, tbl_mB1)
doc <- add_model(doc, "Model B2: clinical predictors + time",
                 paste("Reference as model B1. Time as a fixed effect with a random slope",
                       "per patient."),
                 mB2$reml, tbl_mB2)
if (!is.null(ft_aic_B)) {
  doc <- doc |> body_add_break() |>
    body_add_par("AIC: models B1 and B2", style = "heading 2") |>
    body_add_flextable(ft_aic_B)
}

print(doc, target = file.path(out, "egfr_tables.docx"))
cat("Outputs written to", out, "\n")
