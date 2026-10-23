# Length of stay after kidney transplantation: negative binomial regression.
# Univariable models by predictor group and two multivariable models.

source("R/helpers.R")
install_missing(c("readxl", "MASS", "performance", "see", "gtsummary", "broom",
                  "officer", "rvg", "flextable", "dplyr", "ggdist", "beeswarm"))
suppressPackageStartupMessages({
  library(MASS); library(readxl); library(performance); library(see)
  library(gtsummary); library(broom); library(officer); library(rvg)
  library(flextable); library(dplyr); library(ggdist)
})

theme_gtsummary_compact()
set_flextable_defaults(font.family = "Arial", font.size = 9)

in_file <- "data/cohort.xlsx"
out <- "output/length_of_stay"
dir.create(out, showWarnings = FALSE, recursive = TRUE)

prs <- read_pptx()
slide_n <- 0

los_breaks <- c(5, 10, 15, 20, 25, 50, 75)

box_h <- function(data, xvar, ylab, title, colors) {
  data |> filter(!is.na(.data[[xvar]]), !is.na(LOS_days)) |>
    ggplot(aes(x = LOS_days, y = .data[[xvar]], fill = .data[[xvar]])) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.4, colour = "#555555",
                 alpha = 0.55, width = 0.6) +
    geom_swarm(colour = dot_col, alpha = dot_alpha, size = dot_size) +
    scale_fill_manual(values = colors) +
    scale_x_continuous(breaks = los_breaks, minor_breaks = NULL) +
    labs(title = title, x = "LOS (days)", y = ylab) +
    theme_thesis()
}

scatter_p <- function(data, xvar, xlabel, title) {
  data |> filter(!is.na(.data[[xvar]]), !is.na(LOS_days)) |>
    ggplot(aes(x = .data[[xvar]], y = LOS_days)) +
    geom_point(colour = col_scatter_pt, alpha = 0.55, size = 1.8) +
    geom_smooth(method = "lm", colour = col_scatter_line, se = FALSE, linewidth = 1) +
    scale_y_continuous(breaks = los_breaks, minor_breaks = NULL) +
    labs(title = title, x = xlabel, y = "LOS (days)") +
    theme_thesis()
}

irr_forest <- function(model, title) {
  tryCatch({
    co <- broom::tidy(model, conf.int = TRUE, exponentiate = TRUE) |>
      filter(term != "(Intercept)") |>
      mutate(sig = p.value < 0.05, term = factor(term, levels = rev(unique(term))))
    if (nrow(co) == 0) return(NULL)
    ggplot(co, aes(x = estimate, y = term, colour = sig)) +
      geom_vline(xintercept = 1, linetype = "dashed", colour = col_ref, linewidth = 0.5) +
      geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.16, linewidth = 0.7) +
      geom_point(size = 3.1) +
      scale_colour_manual(values = c("FALSE" = col_ns, "TRUE" = col_sig),
                          labels = c("FALSE" = "p ≥ 0.05", "TRUE" = "p < 0.05"),
                          name = NULL) +
      scale_x_log10() +
      labs(title = title, subtitle = "IRR with 95% CI (log scale)",
           x = "IRR (<1 shorter LOS, >1 longer LOS)", y = NULL) +
      theme_thesis() +
      theme(legend.position = "bottom", axis.text.y = element_text(hjust = 1))
  }, error = function(e) NULL)
}

nb_ft <- function(model, outcome_label = "LOS") {
  tryCatch({
    fe <- broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
    fe_rows <- data.frame(
      Predictors = fe$term,
      IRR = ifelse(is.na(fe$estimate), "", sprintf("%.3f", fe$estimate)),
      CI  = ifelse(is.na(fe$conf.low), "",
                   sprintf("%.3f to %.3f", fe$conf.low, fe$conf.high)),
      p   = fmt_p(fe$p.value),
      is_sig = !is.na(fe$p.value) & fe$p.value < 0.05,
      row_type = "fe", stringsAsFactors = FALSE)
    fit_rows <- data.frame(
      Predictors = c("Model fit", "AIC", "θ (dispersion)", "N observations"),
      IRR = c("", round(AIC(model), 1), round(model$theta, 3), nobs(model)),
      CI = "", p = "", is_sig = FALSE, row_type = c("hdr", "fit", "fit", "fit"),
      stringsAsFactors = FALSE)
    d_all <- rbind(fe_rows, fit_rows)
    nfe <- nrow(fe_rows)

    flextable(d_all[, c("Predictors", "IRR", "CI", "p")]) |>
      set_header_labels(IRR = "IRR", CI = "95% CI", p = "p-value") |>
      add_header_row(values = c("", outcome_label, "", ""), colwidths = c(1, 1, 1, 1)) |>
      align(i = 1, part = "header", align = "center") |>
      bold(part = "header") |>
      italic(i = 2, part = "header") |>
      bold(i = nfe + 1, j = 1, part = "body") |>
      bg(i = nfe + 1, bg = "#F0F0F0", part = "body") |>
      bold(i = which(d_all$is_sig), j = "p", part = "body") |>
      color(i = which(d_all$is_sig), j = c("IRR", "CI", "p"), color = col_sig, part = "body") |>
      color(i = which(!d_all$is_sig & d_all$row_type == "fe"), j = c("IRR", "CI", "p"),
            color = col_ns, part = "body") |>
      border_remove() |>
      hline_top(part = "header", border = fp_border(width = 1.5)) |>
      hline(i = 2, part = "header", border = fp_border(width = 0.5)) |>
      hline_bottom(part = "header", border = fp_border(width = 1)) |>
      hline(i = nfe, border = fp_border(width = 0.5, style = "dashed"), part = "body") |>
      hline_bottom(part = "body", border = fp_border(width = 1)) |>
      fontsize(size = 9, part = "body") |>
      fontsize(size = 10, part = "header") |>
      width(j = 1, width = 2.5) |> width(j = 2, width = 0.9) |>
      width(j = 3, width = 1.5) |> width(j = 4, width = 0.7)
  }, error = function(e) { cat("  nb_ft failed:", conditionMessage(e), "\n"); NULL })
}

# ---- Data -------------------------------------------------------------------
raw <- read_excel(in_file, na = c("", "NA", "N/A"))
names(raw)[grepl("cRF%", names(raw), fixed = TRUE)] <- "cRF_pct"
names(raw)[grepl("complex_uro", names(raw))]        <- "ComplexUro"

d <- raw |>
  filter(!is.na(over_60)) |>
  mutate(
    across(c(age_at_tx, LOS_days, CIT_min, donor_age, charlson_comorbidity_index,
             rec_BMI, donor_eGFR_admission, UKKDRI, UKKRRI), as.numeric),
    CIT_hours = CIT_min / 60,
    CIT_10min = CIT_min / 10,
    age_group = cut(age_at_tx, breaks = c(0, 49, 59, 69, 100),
                    labels = c("<50", "50-59", "60-69", "70-79")),
    over60f   = factor(over_60, levels = c(0, 1), labels = c("<60", "≥60")),
    eth_bin   = factor(ifelse(grepl("white", tolower(ethnicity)) &
                                !grepl("non", tolower(ethnicity)), "White", "Non-white"),
                       levels = c("Non-white", "White")),
    TxGroup   = factor(case_when(grepl("^LD", tx_type) ~ "LD",
                                 grepl("DBD", tx_type) ~ "DBD",
                                 grepl("DCD", tx_type) ~ "DCD"),
                       levels = c("DCD", "DBD", "LD")),
    CatTx     = factor(tx_categ, levels = c(0, 1), labels = c("Live", "Deceased")),
    anat_cat  = factor(case_when(
      grepl("standard anatomy/single",  anatomy, ignore.case = TRUE) ~ "Standard",
      grepl("not requiring back table", anatomy, ignore.case = TRUE) ~ "Non-standard, no back table",
      grepl("requiring back table",     anatomy, ignore.case = TRUE) ~ "Non-standard, back table required"),
      levels = c("Standard", "Non-standard, no back table", "Non-standard, back table required")),
    HLA_MM    = factor(HLA_MM_level),
    WHO_PS    = factor(WHO_PS),
    ASA       = factor(ASA),
    DonorHTNf = factor(suppressWarnings(as.numeric(as.character(donor_HTN))),
                       levels = c(0, 1), labels = c("No", "Yes")),
    donor_prot = suppressWarnings(as.numeric(as.character(donor_proteinuria))),
    DonorProteinuria = factor(ifelse(donor_prot > 0, "Yes", "No"), levels = c("No", "Yes")),
    donor_sex_f = factor(donor_sex),
    UKKDR_q   = factor(case_when(
      UKKDR_quartlie == "live donor" ~ "live donor",
      UKKDR_quartlie %in% c("D1", "D2", "D3", "D4") ~ UKKDR_quartlie),
      levels = c("live donor", "D1", "D2", "D3", "D4")),
    UKKRRI_q  = factor(ifelse(UKKRRI_quartile %in% c("R1", "R2", "R3", "R4"),
                              UKKRRI_quartile, NA),
                       levels = c("R1", "R2", "R3", "R4")),
    induction_IS = factor(case_when(
      grepl("Alemtuzumab", induction_IS, ignore.case = TRUE) ~ "Alemtuzumab",
      grepl("ATG",         induction_IS, ignore.case = TRUE) ~ "ATG",
      TRUE ~ "Basiliximab"),
      levels = c("Basiliximab", "Alemtuzumab", "ATG")),
    imm_risk     = factor(imm_risk, levels = c("standard", "intermediate", "high")),
    rec_diabetes = factor(rec_diabetes),
    rec_smoking  = factor(rec_smoking_status),
    dialysis_status = factor(ifelse(dialysis_status == "pre", "Pre-emptive",
                                    "On dialysis (HD/PD)"),
                             levels = c("Pre-emptive", "On dialysis (HD/PD)")),
    ComplexUro    = factor(suppressWarnings(as.numeric(as.character(ComplexUro))),
                           levels = c(0, 1), labels = c("No", "Yes")),
    ABDR_MM_sum   = suppressWarnings(as.numeric(as.character(ABDR_MM_sum))),
    donor_LOS_num = suppressWarnings(as.numeric(as.character(donor_LOS)))
  )

n <- nrow(d)
cat(sprintf("n=%d | LOS median %g (IQR %g-%g), range %g-%g\n", n,
            median(d$LOS_days, na.rm = TRUE),
            quantile(d$LOS_days, .25, na.rm = TRUE), quantile(d$LOS_days, .75, na.rm = TRUE),
            min(d$LOS_days, na.rm = TRUE), max(d$LOS_days, na.rm = TRUE)))

# ---- Distribution of LOS ----------------------------------------------------
los_med <- median(d$LOS_days, na.rm = TRUE)
los_q   <- quantile(d$LOS_days, c(0.25, 0.75), na.rm = TRUE)
bw_los  <- 3

hist_los <- d |> filter(!is.na(LOS_days)) |>
  ggplot(aes(x = LOS_days)) +
  geom_histogram(binwidth = bw_los, boundary = 0, fill = "#0072B2",
                 colour = "#FFFFFF", linewidth = 0.3, alpha = 0.9) +
  geom_vline(xintercept = los_med, linetype = "dashed", colour = "#444444", linewidth = 0.5) +
  annotate("text", x = los_med, y = Inf, vjust = 1.6, hjust = -0.1, size = 3.1,
           colour = "#444444",
           label = sprintf("median %g (IQR %g–%g)", los_med, los_q[1], los_q[2])) +
  scale_x_continuous(breaks = seq(0, 90, 10), minor_breaks = NULL) +
  labs(title = "Distribution of hospital length of stay",
       x = "LOS (days)", y = "Number of recipients") +
  theme_thesis()
add_plot(hist_los, "LOS distribution")

# marginal negative binomial fit overlaid on the observed counts
nb_marg  <- glm.nb(LOS_days ~ 1, data = filter(d, !is.na(LOS_days)))
nb_curve <- data.frame(x = 0:max(d$LOS_days, na.rm = TRUE))
nb_curve$y <- sum(!is.na(d$LOS_days)) * bw_los *
  dnbinom(nb_curve$x, size = nb_marg$theta, mu = exp(coef(nb_marg)[[1]]))

hist_los_nb <- d |> filter(!is.na(LOS_days)) |>
  ggplot(aes(x = LOS_days)) +
  geom_histogram(binwidth = bw_los, boundary = 0, fill = "#0072B2",
                 colour = "#FFFFFF", linewidth = 0.3, alpha = 0.85) +
  geom_line(data = nb_curve, aes(x = x, y = y), colour = "#D55E00", linewidth = 0.9) +
  annotate("text", x = 38, y = max(nb_curve$y) * 0.72, hjust = 0, size = 3.1,
           label = "fitted negative binomial", colour = "#D55E00") +
  scale_x_continuous(breaks = seq(0, 90, 10), minor_breaks = NULL) +
  labs(title = "Length of stay with fitted negative binomial",
       x = "LOS (days)", y = "Number of recipients") +
  theme_thesis()
add_plot(hist_los_nb, "LOS with NB fit")

# ---- Predictor plots --------------------------------------------------------
yn_cols <- c("No" = "#009E73", "Yes" = "#D55E00")
add_plot(box_h(d, "age_group", "Age group", "LOS by age group",
               setNames(pal_age, c("<50", "50-59", "60-69", "70-79"))), "age group")
add_plot(box_h(d, "over60f", "Age ≥60", "LOS by age ≥60",
               c("<60" = col_ns, "≥60" = col_sig)), "age >=60")
add_plot(scatter_p(d, "age_at_tx", "Recipient age (years)", "LOS vs recipient age"), "age")
add_plot(box_h(d, "TxGroup", "Transplant type", "LOS by transplant type",
               c("DCD" = "#D55E00", "DBD" = "#E69F00", "LD" = "#009E73")), "transplant type")
add_plot(box_h(d, "CatTx", "Transplant category", "LOS by donor type",
               c("Live" = "#009E73", "Deceased" = "#D55E00")), "transplant category")
add_plot(box_h(d, "DonorHTNf", "Donor hypertension", "LOS by donor hypertension", yn_cols),
         "donor HTN")
add_plot(box_h(d, "DonorProteinuria", "Donor proteinuria", "LOS by donor proteinuria", yn_cols),
         "donor proteinuria")
add_plot(box_h(d, "UKKRRI_q", "UKKRRI quartile", "LOS by UKKRRI quartile",
               setNames(pal_ukkrri, c("R1", "R2", "R3", "R4"))), "UKKRRI quartile")
add_plot(box_h(d, "UKKDR_q", "UKKDRI quartile", "LOS by UKKDRI quartile",
               setNames(pal_ukkdr, c("live donor", "D1", "D2", "D3", "D4"))), "UKKDRI quartile")
add_plot(scatter_p(d, "donor_age", "Donor age (years)", "LOS vs donor age"), "donor age")
add_plot(scatter_p(d, "CIT_hours", "CIT (hours)", "LOS vs cold ischaemia time"), "CIT")
add_plot(box_h(d, "dialysis_status", "Dialysis status", "LOS by dialysis status",
               c("Pre-emptive" = "#009E73", "On dialysis (HD/PD)" = "#D55E00")), "dialysis")
add_plot(scatter_p(d, "UKKRRI", "UKKRRI", "LOS vs UKKRRI"), "UKKRRI")
add_plot(scatter_p(d, "UKKDRI", "UKKDRI", "LOS vs UKKDRI"), "UKKDRI")

# ---- Univariable models -----------------------------------------------------
var_labels <- list(
  age_at_tx ~ "Recipient age (years)", age_group ~ "Age group",
  over60f ~ "Age ≥60", rec_sex ~ "Sex", eth_bin ~ "Ethnicity",
  rec_smoking ~ "Smoking status", rec_diabetes ~ "Diabetes",
  rec_BMI ~ "BMI", ComplexUro ~ "Complex urology",
  charlson_comorbidity_index ~ "Charlson CCI",
  ASA ~ "ASA grade", WHO_PS ~ "WHO performance status",
  TxGroup ~ "Transplant type (ref: DCD)",
  anat_cat ~ "Donor anatomy", HLA_MM ~ "HLA mismatch level",
  donor_age ~ "Donor age (years)", donor_sex_f ~ "Donor sex",
  DonorHTNf ~ "Donor hypertension",
  donor_eGFR_admission ~ "Donor eGFR at admission",
  UKKRRI ~ "UKKRRI (continuous)", UKKRRI_q ~ "UKKRRI quartile",
  UKKDR_q ~ "UKKDRI quartile", UKKDRI ~ "UKKDRI (continuous)",
  ABDR_MM_sum ~ "ABDR mismatch sum", DonorProteinuria ~ "Donor proteinuria",
  donor_LOS_num ~ "Donor LOS (days)",
  CIT_10min ~ "CIT (per 10 min)", CIT_hours ~ "CIT (hours)",
  induction_IS ~ "Induction immunosuppression", imm_risk ~ "Immunological risk",
  dialysis_status ~ "Dialysis status"
)

make_uv <- function(vars, caption, data = d) {
  lbl <- var_labels[sapply(var_labels, function(f) all.vars(f)[1] %in% vars)]
  data |>
    dplyr::select(LOS_days, all_of(vars)) |>
    tbl_uvregression(y = LOS_days, method = glm.nb, exponentiate = TRUE,
                     hide_n = FALSE, label = lbl) |>
    modify_header(label = "**Predictor**", n_obs = "**N**",
                  estimate = "**IRR (95% CI)**", p.value = "**p-value**") |>
    modify_footnote(estimate ~ "IRR, incidence rate ratio. IRR <1 shorter LOS, >1 longer LOS.") |>
    modify_caption(caption) |>
    bold_p(t = 0.05)
}

tbl_s1 <- make_uv(c("age_at_tx", "age_group", "over60f", "rec_sex", "eth_bin"),
                  sprintf("Demographics (n=%d)", n))
tbl_s2 <- make_uv(c("rec_smoking", "rec_diabetes", "rec_BMI", "ComplexUro",
                    "charlson_comorbidity_index", "ASA", "WHO_PS", "dialysis_status"),
                  sprintf("Comorbidities (n=%d)", n))
tbl_s3 <- make_uv(c("TxGroup", "anat_cat", "HLA_MM", "donor_age", "donor_sex_f",
                    "DonorHTNf", "donor_eGFR_admission"),
                  sprintf("Transplant, anatomy and donor characteristics (n=%d)", n))
tbl_s4 <- make_uv(c("UKKRRI", "UKKRRI_q", "UKKDR_q"),
                  sprintf("UKKRRI and UKKDRI (n=%d)", n))

# UKKDRI is defined for deceased donors only
d_dec <- d |> filter(!is.na(UKKDRI), TxGroup %in% c("DCD", "DBD"))
tbl_s4b <- make_uv("UKKDRI", sprintf("UKKDRI (continuous), deceased donors (n=%d)", nrow(d_dec)),
                   data = d_dec)

tbl_s5 <- make_uv(c("ABDR_MM_sum", "DonorProteinuria", "DonorHTNf", "donor_LOS_num",
                    "CIT_10min", "CIT_hours", "induction_IS", "imm_risk"),
                  sprintf("Ischaemia, donor and immunological factors (n=%d)", n))

# full model tables for predictors with p<0.05
sig_table <- function(varname) {
  tryCatch({
    m <- glm.nb(as.formula(paste0("LOS_days ~ `", varname, "`")),
                data = filter(d, !is.na(.data[[varname]])))
    co <- broom::tidy(m) |> filter(term != "(Intercept)")
    if (all(co$p.value >= 0.05, na.rm = TRUE)) NULL else nb_ft(m)
  }, error = function(e) NULL)
}
sig_vars <- list(
  s1 = "age_group",
  s2 = c("rec_BMI", "dialysis_status"),
  s3 = c("TxGroup", "CatTx", "donor_age", "DonorHTNf"),
  s4 = c("UKKRRI_q", "UKKDR_q"),
  s5 = c("ABDR_MM_sum", "donor_LOS_num", "CIT_hours", "DonorProteinuria")
)
sig_tables <- lapply(sig_vars, function(v) Filter(Negate(is.null), lapply(v, sig_table)))

# ---- Multivariable models ---------------------------------------------------
fit_nb <- function(f, data) {
  tryCatch(glm.nb(as.formula(f), data = data),
           error = function(e) { cat("  glm.nb failed:", conditionMessage(e), "\n"); NULL })
}

# Model A: recipient and donor risk indices
d_mA <- d |> filter(!is.na(UKKRRI_q), !is.na(UKKDR_q))
mA <- fit_nb("LOS_days ~ UKKRRI_q + UKKDR_q", d_mA)
lab_A <- list(UKKRRI_q ~ "UKKRRI quartile", UKKDR_q ~ "UKKDRI quartile")

tbl_mA <- ft_mA <- NULL
if (!is.null(mA)) {
  tbl_mA <- tbl_merge(
    list(d_mA |> dplyr::select(LOS_days, UKKRRI_q, UKKDR_q) |>
           tbl_uvregression(y = LOS_days, method = glm.nb, exponentiate = TRUE,
                            label = lab_A) |> bold_p(t = 0.05),
         tbl_regression(mA, exponentiate = TRUE, label = lab_A) |> bold_p(t = 0.05)),
    tab_spanner = c("**Univariable**", "**Multivariable**")) |>
    modify_caption(sprintf("Model A: UKKRRI + UKKDRI quartile | n=%d | AIC=%.1f | θ=%.3f",
                           nobs(mA), AIC(mA), mA$theta))
  add_plot(irr_forest(mA, "Model A: UKKRRI + UKKDRI quartile"), "model A forest")
  ft_mA <- nb_ft(mA)
}

# Model B: recipient age, donor factors and donor type by cold ischaemia time
d_mB <- d |> filter(!is.na(donor_age), !is.na(DonorHTNf), !is.na(DonorProteinuria))
mB <- fit_nb("LOS_days ~ age_at_tx + donor_age + DonorHTNf + DonorProteinuria + CatTx * CIT_hours",
             d_mB)
lab_B <- list(age_at_tx ~ "Recipient age (years)", donor_age ~ "Donor age (years)",
              DonorHTNf ~ "Donor hypertension", DonorProteinuria ~ "Donor proteinuria",
              CatTx ~ "Transplant category")

tbl_mB <- ft_mB <- NULL
if (!is.null(mB)) {
  tbl_mB <- tbl_merge(
    list(d_mB |> dplyr::select(LOS_days, age_at_tx, donor_age, DonorHTNf,
                               DonorProteinuria, CatTx) |>
           tbl_uvregression(y = LOS_days, method = glm.nb, exponentiate = TRUE,
                            label = lab_B) |> bold_p(t = 0.05),
         tbl_regression(mB, exponentiate = TRUE,
                        label = c(lab_B, list(`CatTx:CIT_hours` ~ "Transplant category × CIT (hours)"))) |>
           bold_p(t = 0.05)),
    tab_spanner = c("**Univariable**", "**Multivariable**")) |>
    modify_caption(sprintf("Model B | n=%d | AIC=%.1f | θ=%.3f",
                           nobs(mB), AIC(mB), mB$theta))
  add_plot(irr_forest(mB, "Model B: clinical predictors"), "model B forest")

  dx_B <- tryCatch({ set.seed(123); plot(check_model(mB)) },
                   error = function(e) NULL)
  if (!is.null(dx_B)) add_plot(dx_B, "model B diagnostics")
  ft_mB <- nb_ft(mB)
}

aic_df <- do.call(rbind, Filter(Negate(is.null), list(
  if (!is.null(mA)) data.frame(Model = "A", Predictors = "UKKRRI quartile + UKKDRI quartile",
                               n = nobs(mA), AIC = round(AIC(mA), 1), theta = round(mA$theta, 3)),
  if (!is.null(mB)) data.frame(Model = "B",
                               Predictors = "Recipient age + donor age + donor HTN + donor proteinuria + transplant category × CIT",
                               n = nobs(mB), AIC = round(AIC(mB), 1), theta = round(mB$theta, 3))
)))
ft_aic <- NULL
if (!is.null(aic_df)) {
  print(aic_df)
  ft_aic <- flextable(aic_df) |>
    set_header_labels(theta = "θ") |>
    bold(part = "header") |>
    add_footer_lines("AIC is comparable only between models fitted to the same observations.") |>
    autofit()
}

print(prs, target = file.path(out, "los_figures.pptx"))

# ---- Word document ----------------------------------------------------------
to_ft <- function(tbl) tbl |> as_flex_table() |> fontsize(size = 10, part = "header") |> autofit()

br  <- function(doc) body_add_par(doc, " ", style = "Normal")
add_section <- function(doc, heading, tbls, extra = list(), note = NULL) {
  doc <- body_add_par(doc, heading, style = "heading 2")
  if (!is.null(note)) doc <- body_add_par(doc, note, style = "Normal")
  for (t in tbls) doc <- doc |> body_add_flextable(to_ft(t)) |> br()
  for (ft in extra) doc <- doc |> body_add_flextable(ft) |> br()
  body_add_break(doc)
}

doc <- read_docx() |>
  body_add_par("Length of stay: negative binomial regression", style = "heading 1") |>
  body_add_par(sprintf("n=%d. IRR = exp(β); IRR <1 shorter LOS, IRR >1 longer LOS.", n),
               style = "Normal") |>
  br()
doc <- add_section(doc, "Demographics", list(tbl_s1), sig_tables$s1)
doc <- add_section(doc, "Comorbidities", list(tbl_s2), sig_tables$s2)
doc <- add_section(doc, "Transplant, anatomy and donor characteristics", list(tbl_s3),
                   sig_tables$s3,
                   note = "Reference: DCD (transplant type), standard anatomy. Donor age was not available for anonymised altruistic donors and imported organs.")
doc <- add_section(doc, "UKKRRI and UKKDRI", list(tbl_s4, tbl_s4b), sig_tables$s4,
                   note = "UKKDRI as a continuous variable is reported for deceased donors only.")
doc <- add_section(doc, "Ischaemia, donor and immunological factors", list(tbl_s5),
                   sig_tables$s5,
                   note = "Donor LOS is the donor's length of hospital stay, a component of the UKKDRI.")

doc <- body_add_par(doc, "Multivariable models", style = "heading 1")
if (!is.null(ft_aic)) doc <- doc |> body_add_flextable(ft_aic) |> br()

doc <- body_add_par(doc, "Model A: UKKRRI quartile + UKKDRI quartile", style = "heading 2") |>
  body_add_par("Reference: UKKRRI R1, UKKDRI live donor.", style = "Normal")
if (!is.null(ft_mA))  doc <- doc |> body_add_flextable(ft_mA) |> br()
if (!is.null(tbl_mA)) doc <- doc |> body_add_flextable(to_ft(tbl_mA)) |> br()

doc <- body_add_par(doc, "Model B: clinical predictors", style = "heading 2") |>
  body_add_par(paste("Reference: live donor, no donor hypertension, no donor proteinuria.",
                     "Transplant category and CIT are strongly correlated because CIT is short in live donation."),
               style = "Normal")
if (!is.null(ft_mB))  doc <- doc |> body_add_flextable(ft_mB) |> br()
if (!is.null(tbl_mB)) doc <- doc |> body_add_flextable(to_ft(tbl_mB))

print(doc, target = file.path(out, "los_tables.docx"))
cat("Outputs written to", out, "\n")
