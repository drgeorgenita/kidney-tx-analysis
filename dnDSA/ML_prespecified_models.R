# De novo DSA: validation of the ML predictors reported in previous study

source("R/helpers.R")
install_missing(c("readxl", "dplyr", "logistf", "pROC", "xgboost",
                  "officer", "rvg", "flextable", "scales"))
suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(logistf); library(pROC); library(xgboost)
  library(officer); library(rvg); library(flextable); library(scales)
})

in_file  <- "data/dsa_cohort.xlsx"
out      <- "output/dnDSA_prespecified"
OUTCOME  <- "de_novoDSA"
CIT_COL  <- "CIT_min"
SPLIT_SEED <- 1000; MODEL_SEED <- 1202
NROUNDS <- 20000; MAX_DEPTH <- 5; LR <- 0.01; N_PERM <- 200
dir.create(out, showWarnings = FALSE, recursive = TRUE)
op <- function(f) file.path(out, f)

theme_dsa <- function() theme_thesis() +
  theme(legend.position = "bottom", legend.title = element_blank())

# ---- Output collectors ------------------------------------------------------
doc  <- read_docx()
deck <- read_pptx()

add_tbl <- function(ft, heading, csv = NULL, df = NULL) {
  doc <<- doc |> body_add_par(heading, style = "heading 2") |>
    body_add_flextable(ft) |> body_add_par("")
  if (!is.null(csv)) write.csv(df, op(csv), row.names = FALSE)
}
add_fig <- function(gg, title, w = 9, h = 5.4) {
  deck <<- add_slide(deck, layout = "Title and Content", master = "Office Theme")
  deck <<- ph_with(deck, value = title, location = ph_location_type(type = "title"))
  deck <<- ph_with(deck, value = dml(ggobj = gg),
                   location = ph_location(left = 0.5, top = 1.3, width = w, height = h))
}
style_ft <- function(ft) ft |> theme_booktabs() |> fontsize(size = 9, part = "body") |>
  fontsize(size = 10, part = "header") |> bold(part = "header") |> autofit()
add_note <- function(ft, txt) ft |> add_footer_lines(txt) |>
  fontsize(size = 7, part = "footer") |> italic(part = "footer")
pf <- function(p) ifelse(is.na(p), "-", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
orstr <- function(o, l, h) ifelse(is.na(o), "-", sprintf("%.2f (%.2f–%.2f)", o, l, h))

# ---- Data -------------------------------------------------------------------
d <- read_excel(in_file)
d[[OUTCOME]] <- as.integer(as.character(d[[OUTCOME]]))
d <- d[!is.na(d[[OUTCOME]]), ]

if (CIT_COL %in% names(d)) {
  d$CIT_hours <- as.numeric(as.character(d[[CIT_COL]])) / 60
} else {
  message("Column ", CIT_COL, " not found; cold ischaemia time is omitted from the models")
}
d$depleting_induction <- factor(
  ifelse(grepl("Alem", d$induction_IS, ignore.case = TRUE), "Depleting",
  ifelse(grepl("Basil", d$induction_IS, ignore.case = TRUE), "Non-depleting", NA)),
  levels = c("Non-depleting", "Depleting"))
for (v in intersect(c("age_at_tx", "HLA_MM_level"), names(d)))
  d[[v]] <- as.numeric(as.character(d[[v]]))
n_tot <- nrow(d); n_ev <- sum(d[[OUTCOME]] == 1)
cat(sprintf("n=%d, dnDSA events=%d (%.1f%%)\n", n_tot, n_ev, 100 * n_ev / n_tot))

# ---- Univariable Firth models -----------------------------------------------
# One row per continuous predictor; a header row and one row per level for
# categorical predictors. AUC (DeLong CI) is for the univariable model.
firth_uni <- function(var, type, label, ref = NULL, levels = NULL) {
  if (!var %in% names(d)) { message("not in data: ", var); return(NULL) }
  sub <- d[, c(OUTCOME, var)]; names(sub) <- c("y", "x")
  sub <- sub[stats::complete.cases(sub), ]
  if (type == "cont") sub$x <- as.numeric(as.character(sub$x))
  else sub$x <- factor(as.character(sub$x), levels = if (!is.null(ref)) c(ref, levels) else NULL)
  if (nrow(sub) < 5 || length(unique(sub$y)) < 2) return(NULL)
  f <- logistf(y ~ x, data = sub)
  cf <- names(coef(f))
  auc_s <- tryCatch({
    ci <- as.numeric(pROC::ci.auc(pROC::roc(sub$y, as.numeric(predict(f, type = "response")),
                                            quiet = TRUE)))
    sprintf("%.3f (%.3f–%.3f)", ci[2], ci[1], ci[3])
  }, error = function(e) "-")
  R <- function(level, N, ev, OR, lo, hi, p, is_ref, is_hdr)
    data.frame(var = var, label = label, level = level, N = N, ev = ev, OR = OR, lo = lo,
               hi = hi, p = p, is_ref = is_ref, is_hdr = is_hdr, auc = auc_s,
               stringsAsFactors = FALSE)
  if (type == "cont") {
    return(R(NA, nrow(sub), sum(sub$y), exp(coef(f)[["x"]]), exp(f$ci.lower[["x"]]),
             exp(f$ci.upper[["x"]]), f$prob[["x"]], FALSE, FALSE))
  }
  lv <- levels(sub$x)
  rbind(R(NA, nrow(sub), sum(sub$y), NA, NA, NA, NA, FALSE, TRUE),
        do.call(rbind, lapply(seq_along(lv), function(i) {
          tn <- paste0("x", lv[i]); sel <- sub$x == lv[i]
          if (i == 1 || !tn %in% cf) R(lv[i], sum(sel), sum(sub$y[sel]), NA, NA, NA, NA, i == 1, FALSE)
          else R(lv[i], sum(sel), sum(sub$y[sel]), exp(coef(f)[[tn]]), exp(f$ci.lower[[tn]]),
                 exp(f$ci.upper[[tn]]), f$prob[[tn]], FALSE, FALSE)
        })))
}
run_spec <- function(s) do.call(firth_uni, stats::setNames(
  list(s[[1]], s[[2]], s[[3]],
       if (length(s) >= 4) s[[4]] else NULL,
       if (length(s) >= 5) s[[5]] else NULL),
  c("var", "type", "label", "ref", "levels")))

# variable, type, label, reference level, other levels
SCREEN <- list(
  list("age_at_tx", "cont", "Recipient age (per year)"),
  list("over_60", "cat", "Recipient age ≥60", "0", "1"),
  list("rec_sex", "cat", "Sex (female vs male)", "M", "F"),
  list("non_white", "cat", "Ethnicity (non-white vs white)", "0", "1"),
  list("dialysis_status", "cat", "Dialysis status", "HD", c("PD", "pre")),
  list("rec_BMI", "cont", "Body mass index (per unit)"),
  list("rec_smoking_status", "cat", "Smoking status", "0", c("1", "2")),
  list("rec_HTN", "cat", "Hypertension", "0", "1"),
  list("rec_diabetes", "cat", "Diabetes", "0", "1"),
  list("charlson_comorbidity_index", "cont", "Charlson Comorbidity Index (per point)"),
  list("ASA", "cont", "ASA grade (per grade)"),
  list("WHO_PS", "cont", "WHO performance status (per level)"),
  list("UKKRRI_quartile", "cat", "UKKRRI quartile", "R1", c("R2", "R3", "R4")),
  list("donor_age", "cont", "Donor age (per year)"),
  list("donor_sex", "cat", "Donor sex", "F", "M"),
  list("UKKDRI_quartile", "cat", "UKKDRI quartile", "D1", c("D2", "D3", "D4")),
  list("HLA_MM_level", "cont", "HLA mismatch level (per level)"),
  list("tx_categ", "cat", "Transplant category (deceased vs live)", "0", "1"),
  list("CIT_hours", "cont", "Cold ischaemia time (per hour)"),
  list("induction_IS", "cat", "Induction regimen", "Alemtuzumab1", c("Alemtuzumab2", "Basiliximab")),
  list("depleting_induction", "cat", "Depleting vs non-depleting induction", "Non-depleting", "Depleting"),
  list("recipient_CMV", "cat", "Recipient CMV serostatus", "0", "1"),
  list("donor_CMV", "cat", "Donor CMV serostatus", "0", "1"),
  list("recipient_EBV", "cat", "Recipient EBV serostatus", "0", "1"),
  list("donor_EBV", "cat", "Donor EBV serostatus", "0", "1"),
  list("flow_cytoT", "cat", "Flow-cytometry crossmatch (T cell)", "negative", "Non-negative"),
  list("flow_cytoB", "cat", "Flow-cytometry crossmatch (B cell)", "negative", "positive"),
  list("pregnancy_history", "cat", "Pregnancy history", "0", "1"),
  list("prev_transfusions", "cat", "Previous transfusions", "0", "1")
)
# recipient and donor risk indices are reported together or not at all
idx <- c("UKKRRI_quartile", "UKKDRI_quartile")
if (!all(idx %in% names(d))) {
  SCREEN <- Filter(function(x) !x[[1]] %in% idx, SCREEN)
  message("Risk indices omitted: ", paste(setdiff(idx, names(d)), collapse = ", "), " not in data")
}

uni <- do.call(rbind, lapply(SCREEN, run_spec))
ok <- !is.na(uni$p)
uni$p_BH <- NA_real_
uni$p_BH[ok] <- p.adjust(uni$p[ok], method = "BH")

uni_disp <- data.frame(
  Predictor = ifelse(is.na(uni$level), uni$label, paste0(" ", uni$level)),
  N = as.character(uni$N),
  Events = as.character(uni$ev),
  `OR (95% CI)` = ifelse(uni$is_hdr, "", ifelse(uni$is_ref, "ref", orstr(uni$OR, uni$lo, uni$hi))),
  p = ifelse(uni$is_hdr | uni$is_ref, "", pf(uni$p)),
  `p (BH)` = ifelse(uni$is_hdr | uni$is_ref, "", pf(uni$p_BH)),
  `AUC (95% CI)` = ifelse(uni$is_hdr | is.na(uni$level), uni$auc, ""),
  check.names = FALSE, stringsAsFactors = FALSE)
sig_uni <- which(!uni$is_hdr & !uni$is_ref & !is.na(uni$p) & uni$p < 0.05)
ft_uni <- style_ft(flextable(uni_disp))
if (length(sig_uni)) ft_uni <- bold(ft_uni, i = sig_uni, j = "p")
ft_uni <- add_note(ft_uni, sprintf(paste(
  "Firth penalised logistic regression (n=%d, %d events). p (BH): Benjamini-Hochberg",
  "adjusted across all reported comparisons. AUC with DeLong 95%% CI."), n_tot, n_ev))
add_tbl(ft_uni, "Univariable associations (Firth)", "univariable_firth_BH.csv", uni_disp)

# ---- Predictor selection diagram --------------------------------------------
fc_box <- function(xmin, xmax, ymin, ymax, title, body, fill = "#EAF1F8", tcol = "#1f3b63") {
  list(
    annotate("rect", xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
             fill = fill, colour = "#9FB6CF", linewidth = 0.5),
    annotate("text", x = (xmin + xmax) / 2, y = ymax - 2.5, label = title,
             fontface = "bold", size = 3.4, colour = tcol, vjust = 1, family = "Arial"),
    annotate("text", x = (xmin + xmax) / 2, y = ymax - 7, label = body,
             size = 2.9, colour = "#333333", vjust = 1, lineheight = 0.98, family = "Arial"))
}
arrow_seg <- function(x, xend, y, yend) annotate("segment", x = x, xend = xend, y = y, yend = yend,
  arrow = arrow(length = unit(0.15, "cm"), type = "closed"), colour = "#5A6B80", linewidth = 0.5)

g_flow <- ggplot() +
  fc_box(4, 58, 76, 99, "Candidate predictors screened univariably (n = 15)",
         paste("Demographics: age, sex, ethnicity, dialysis status",
               "Comorbidity / indices: BMI, smoking, diabetes, CCI, ASA, WHO PS, UKKRRI",
               "Donor / transplant: HLA mismatch, donor age & sex, transplant category, CIT",
               "Viral / induction: recipient CMV, recipient EBV, induction regimen",
               "Crossmatch / sensitisation: crossmatch, pregnancy, transfusions", sep = "\n")) +
  arrow_seg(31, 31, 76, 51) +
  arrow_seg(31, 62, 64, 64) +
  fc_box(62, 100, 50, 78, "Not carried forward to modelling",
         paste("Not tested in the previous study",
               "(comorbidity indices, viral serostatus, crossmatch)",
               "Unavailable in this cohort: amino-acid mismatch",
               "Constant by design: prior sensitisation",
               "(cohort unsensitised, cRF 0%)", sep = "\n"),
         fill = "#F5F0E8", tcol = "#8a5a00") +
  fc_box(4, 58, 26, 50, "Pre-specified predictor set for modelling (n = 6)",
         paste("Recipient age  ·  HLA mismatch level  ·  Sex",
               "Transplant category  ·  Induction regimen  ·  CIT", "",
               "Pre-specified from the previous study,",
               "not selected by univariable significance", sep = "\n"),
         fill = "#E3EFE7", tcol = "#0A6B4A") +
  scale_x_continuous(limits = c(0, 100)) +
  scale_y_continuous(limits = c(24, 101)) +
  labs(title = "Predictor selection for the pre-specified models") +
  theme_void(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 13, colour = "#111111",
                                  hjust = 0, family = "Arial", margin = margin(b = 6)),
        plot.margin = margin(10, 12, 10, 12))
add_fig(g_flow, "Predictor selection", w = 9, h = 6)

# ---- Concordance with the previous study ------------------------------------
# previous-study p-values (BH-adjusted): whole cohort, unsensitised subgroup
PRESPEC <- list(
  list("HLA_MM_level", "cont", "HLA mismatch level (per level)", "<0.001 (BH 0.01)", "0.75 (BH 0.82)"),
  list("depleting_induction", "cat", "Depleting vs non-depleting induction",
       "<0.001 (BH 0.01)", "0.27 (BH 0.51)", "Non-depleting", "Depleting"),
  list("CIT_hours", "cont", "Cold ischaemia time (per hour)", "<0.001 (BH 0.02)", "0.39 (BH 0.60)"),
  list("tx_categ", "cat", "Type of transplant (deceased vs live)",
       "0.22 (BH 0.27)", "0.49 (BH 0.60)", "0", "1"),
  list("age_at_tx", "cont", "Recipient age (per year)", "0.03 (BH 0.06)", "<0.05 (BH 0.30)"),
  list("rec_sex", "cat", "Sex (female vs male)", "0.01 (BH 0.01)", "0.48 (BH 0.60)", "M", "F"),
  list("recipient_CMV", "cat", "Recipient CMV serostatus",
       "not tested", "not tested", "0", "1")
)
con <- do.call(rbind, lapply(PRESPEC, function(s) {
  r <- run_spec(s[c(1:3, if (length(s) >= 6) 6:7)])
  if (is.null(r)) return(NULL)
  r <- r[!r$is_ref & !r$is_hdr, , drop = FALSE][1, ]
  r$prev_whole <- s[[4]]; r$prev_unsens <- s[[5]]
  r
}))
con$p_BH <- p.adjust(con$p, method = "BH")

con_disp <- data.frame(
  Predictor = con$label,
  `OR (95% CI)` = orstr(con$OR, con$lo, con$hi),
  `AUC (95% CI)` = con$auc,
  `This cohort p` = pf(con$p),
  `p (BH)` = pf(con$p_BH),
  `BH p<0.05` = ifelse(con$p_BH < 0.05, "yes", "no"),
  `Previous study, whole cohort` = con$prev_whole,
  `Previous study, unsensitised` = con$prev_unsens,
  check.names = FALSE, stringsAsFactors = FALSE)
add_tbl(add_note(style_ft(flextable(con_disp)), sprintf(paste(
  "Predictors from the previous study tested in this unsensitised cohort (Firth). p (BH):",
  "Benjamini-Hochberg adjustment across these %d predictors. Previous study: whole cohort",
  "n=460 (115 events), unsensitised subgroup n=348 (35 events); BH-adjusted p-values."), nrow(con))),
  "Concordance with the previous study", "concordance_prevstudy.csv", con_disp)

# ---- Pre-specified multivariable model --------------------------------------
CONT <- c("age_at_tx", "HLA_MM_level", "CIT_hours")
REFS <- list(rec_sex = "M", tx_categ = "0", induction_IS = "Alemtuzumab1",
             recipient_CMV = "0", over_60 = "0")

mv_fit <- function(vars) {
  miss <- setdiff(vars, names(d))
  if (length(miss)) message("not in data, dropped: ", paste(miss, collapse = ", "))
  vars <- intersect(vars, names(d))
  dd <- d[, c(OUTCOME, vars)]
  for (v in setdiff(vars, CONT)) {
    dd[[v]] <- factor(dd[[v]])
    if (v %in% names(REFS) && REFS[[v]] %in% levels(dd[[v]]))
      dd[[v]] <- relevel(dd[[v]], ref = REFS[[v]])
  }
  dd <- dd[stats::complete.cases(dd), ]
  list(fit = logistf(as.formula(paste(OUTCOME, "~", paste(vars, collapse = " + "))),
                     data = dd, control = logistf.control(maxit = 1000)),
       dat = dd, n = nrow(dd))
}

# coefficient table for selected terms, with BH adjustment across those terms
mv_coefs <- function(m, rows) {
  cf <- coef(m$fit)
  out <- do.call(rbind, lapply(rows, function(r) {
    i <- match(r[[1]], names(cf))
    if (is.na(i)) return(data.frame(term = r[[1]], label = r[[2]], OR = NA, lo = NA, hi = NA, p = NA))
    data.frame(term = r[[1]], label = r[[2]], OR = exp(cf[i]), lo = exp(m$fit$ci.lower[i]),
               hi = exp(m$fit$ci.upper[i]), p = m$fit$prob[i])
  }))
  out$p_BH <- NA_real_
  out$p_BH[!is.na(out$p)] <- p.adjust(out$p[!is.na(out$p)], method = "BH")
  rownames(out) <- NULL
  out
}

ML_VARS <- c("tx_categ", "induction_IS", "age_at_tx", "rec_sex", "HLA_MM_level", "CIT_hours")
m_pre <- mv_fit(ML_VARS)
# ordered as in the XGBoost importance ranking
rows_pre <- list(
  list("age_at_tx", "Recipient age (per year)"),
  list("HLA_MM_level", "HLA mismatch level (per level)"),
  list("rec_sexF", "Female vs male"),
  list("tx_categ1", "Deceased vs live donor"),
  list("induction_ISBasiliximab", "Basiliximab vs Alemtuzumab x1"),
  list("induction_ISAlemtuzumab2", "Alemtuzumab x2 vs x1"),
  list("CIT_hours", "Cold ischaemia time (per hour)"))
co_pre <- mv_coefs(m_pre, rows_pre)
auc_pre <- as.numeric(pROC::auc(pROC::roc(m_pre$dat[[OUTCOME]],
                                          as.numeric(predict(m_pre$fit, type = "response")),
                                          quiet = TRUE)))

mv_disp <- rbind(
  data.frame(Predictor = co_pre$label, `OR (95% CI)` = orstr(co_pre$OR, co_pre$lo, co_pre$hi),
             p = pf(co_pre$p), `p (BH)` = pf(co_pre$p_BH), check.names = FALSE),
  data.frame(Predictor = "Model AUC-ROC", `OR (95% CI)` = sprintf("%.3f", auc_pre),
             p = "", `p (BH)` = "", check.names = FALSE))
add_tbl(add_note(style_ft(flextable(mv_disp)), sprintf(paste(
  "Firth penalised logistic regression, n=%d. p (BH): Benjamini-Hochberg adjusted across",
  "the %d model terms. AUC is apparent (fitted and evaluated on the same data)."),
  m_pre$n, sum(!is.na(co_pre$p)))),
  "Pre-specified multivariable model (Firth)", "mv_prespecified.csv", mv_disp)

# ---- Forest plots -----------------------------------------------------------
XLIM_OR <- c(0.1, 20)
forest_plot <- function(df, title, subtitle, sort = TRUE) {
  lab_sig <- c("BH-adjusted p<0.05", "BH-adjusted p≥0.05")
  df$sig <- ifelse(!is.na(df$p) & df$p < 0.05, lab_sig[1], lab_sig[2])
  # intervals near 0 or infinity indicate a zero-event category
  df$sep <- is.na(df$lo) | is.na(df$hi) | df$lo < 0.02 | df$hi > 100
  dropped <- unique(df$label[df$sep])
  df <- df[!df$sep, , drop = FALSE]
  if (sort) df <- df[order(df$OR), , drop = FALSE]
  df$label <- factor(df$label, levels = if (sort) df$label else rev(df$label))
  cap <- "Intervals beyond the axis are truncated."
  if (length(dropped)) cap <- paste0("Not shown (zero-event category): ",
                                     paste(dropped, collapse = "; "), ". ", cap)
  ggplot(df, aes(OR, label, colour = sig)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = col_ref) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.22, linewidth = 0.7) +
    geom_point(size = 2.6) +
    scale_x_log10() + coord_cartesian(xlim = XLIM_OR) +
    scale_colour_manual(values = setNames(c(col_sig, col_ns), lab_sig)) +
    labs(title = title, subtitle = subtitle, x = "Odds ratio (log scale)", y = NULL, caption = cap) +
    theme_dsa()
}

co_pre_bh <- transform(co_pre, p = p_BH)
add_fig(forest_plot(co_pre_bh[!is.na(co_pre_bh$OR), ], "Pre-specified multivariable model",
                    "Adjusted Firth odds ratios (95% CI)", sort = FALSE),
        "Multivariable forest plot")

fp <- uni[!uni$is_ref & !is.na(uni$OR), ]
fp$label <- ifelse(is.na(fp$level), fp$label, paste0(fp$label, ": ", fp$level))
fp$p <- fp$p_BH
add_fig(forest_plot(fp[, c("label", "OR", "lo", "hi", "p")], "Univariable associations with dnDSA",
                    "Firth odds ratios (95% CI), Benjamini-Hochberg adjusted"),
        "Univariable forest plot", w = 9, h = 7)
add_fig(forest_plot(fp[fp$var %in% ML_VARS, c("label", "OR", "lo", "hi", "p")],
                    "Univariable associations: pre-specified predictors",
                    "Firth odds ratios (95% CI)"),
        "Univariable forest plot, pre-specified predictors", w = 8.5, h = 4.4)

# ---- Descriptive plots ------------------------------------------------------
grp_cols <- c("dnDSA−" = col_grey, "dnDSA+" = col_sig)
d$.grp <- factor(ifelse(d[[OUTCOME]] == 1, "dnDSA+", "dnDSA−"), levels = names(grp_cols))

bar_pred <- function(var, lab) {
  if (!var %in% names(d)) return(invisible())
  df <- d[!is.na(d[[var]]), ]
  tab <- as.data.frame(prop.table(table(df[[var]], df$.grp), 2) * 100)
  names(tab) <- c("level", "grp", "pct")
  add_fig(ggplot(tab, aes(level, pct, fill = grp)) +
            geom_col(position = "dodge", width = 0.7) +
            scale_fill_manual(values = grp_cols) +
            labs(title = lab, x = NULL, y = "% within group") + theme_dsa(),
          lab)
}
box_pred <- function(var, lab) {
  if (!var %in% names(d)) return(invisible())
  df <- d[!is.na(d[[var]]), ]
  df$v <- as.numeric(as.character(df[[var]]))
  add_fig(ggplot(df, aes(.grp, v, fill = .grp)) +
            geom_boxplot(width = 0.5, outlier.alpha = 0.4) +
            scale_fill_manual(values = grp_cols, guide = "none") +
            labs(title = lab, x = NULL, y = lab) + theme_dsa(),
          lab)
}
box_pred("age_at_tx", "Recipient age (years)")
box_pred("CIT_hours", "Cold ischaemia time (hours)")
box_pred("HLA_MM_level", "HLA mismatch level")
bar_pred("depleting_induction", "Induction (depleting vs non-depleting)")
bar_pred("tx_categ", "Type of transplant")
bar_pred("recipient_CMV", "Recipient CMV serostatus")
bar_pred("rec_sex", "Sex")

# ---- XGBoost and Firth on the same split ------------------------------------
ml_vars <- intersect(ML_VARS, names(d))
dat <- d[, c(OUTCOME, ml_vars)]
dat <- dat[stats::complete.cases(dat), ]
# same reference levels as the regression, so importance maps onto the same terms
for (v in setdiff(ml_vars, CONT)) {
  dat[[v]] <- factor(dat[[v]])
  if (v %in% names(REFS) && REFS[[v]] %in% levels(dat[[v]]))
    dat[[v]] <- relevel(dat[[v]], ref = REFS[[v]])
}
y_all <- as.integer(as.character(dat[[OUTCOME]]))
X_all <- model.matrix(~ ., data = as.data.frame(dat[, ml_vars]))

set.seed(SPLIT_SEED)
train <- sample(c(TRUE, FALSE), nrow(X_all), replace = TRUE)
Xtr <- X_all[train, ]; Xte <- X_all[!train, ]; yte <- y_all[!train]

fit_xgb <- function(y) {
  set.seed(MODEL_SEED)
  xgboost(Xtr, factor(y), nrounds = NROUNDS, max_depth = MAX_DEPTH,
          learning_rate = LR, verbosity = 0)
}
safe_auc <- function(l, s) if (length(unique(l)) < 2) NA_real_ else
  as.numeric(pROC::auc(pROC::roc(l, s, quiet = TRUE)))

mod <- fit_xgb(y_all[train])
pred_tr <- predict(mod, Xtr, type = "response")
pred_te <- predict(mod, Xte, type = "response")
auc_te <- safe_auc(yte, pred_te)

dat2 <- as.data.frame(dat); dat2$.y <- y_all
rhs <- paste(ml_vars, collapse = " + ")
firth_h <- logistf(as.formula(paste(".y ~", rhs)), data = dat2[train, ])
mm_te <- model.matrix(as.formula(paste("~", rhs)), data = dat2[!train, ])
cf <- coef(firth_h)
common <- intersect(colnames(mm_te), names(cf))
firth_pte <- as.numeric(plogis(mm_te[, common] %*% cf[common]))
firth_auc <- safe_auc(yte, firth_pte)

# label permutation: refit on shuffled outcomes with predictors unchanged
perm <- do.call(rbind, lapply(seq_len(N_PERM), function(i) {
  set.seed(2000 + i)
  yp <- sample(y_all)
  if (length(unique(yp[train])) < 2) return(NULL)
  mp <- tryCatch(fit_xgb(yp[train]), error = function(e) NULL)
  if (is.null(mp)) return(NULL)
  data.frame(train_auc = safe_auc(yp[train], predict(mp, Xtr, type = "response")),
             test_auc  = safe_auc(yp[!train], predict(mp, Xte, type = "response")))
}))
p_perm <- (sum(perm$test_auc >= auc_te, na.rm = TRUE) + 1) / (sum(!is.na(perm$test_auc)) + 1)
write.csv(perm, op("permutation_results.csv"), row.names = FALSE)

htab <- data.frame(Model = c("Firth penalised logistic regression", "XGBoost"),
                   `Held-out test AUC` = sprintf("%.3f", c(firth_auc, auc_te)),
                   check.names = FALSE)
add_tbl(add_note(style_ft(flextable(htab)),
                 sprintf("Same predictors and train/test split; %d test recipients, %d events.",
                         length(yte), sum(yte))),
        "Regression vs XGBoost", "ml_headtohead.csv", htab)

imp <- as.data.frame(xgb.importance(model = mod))
gain <- setNames(imp$Gain, imp$Feature)
cmp <- data.frame(Predictor = co_pre$label,
                  gain_num = ifelse(co_pre$term %in% names(gain), gain[co_pre$term], 0),
                  `Multivariable OR (95% CI)` = orstr(co_pre$OR, co_pre$lo, co_pre$hi),
                  p = pf(co_pre$p), `p (BH)` = pf(co_pre$p_BH),
                  check.names = FALSE, stringsAsFactors = FALSE)
cmp <- cmp[order(-cmp$gain_num), ]
cmp$Importance <- sprintf("%.3f", cmp$gain_num)
cmp <- cmp[, c("Predictor", "Importance", "Multivariable OR (95% CI)", "p", "p (BH)")]
add_tbl(add_note(style_ft(flextable(cmp)), paste(
  "XGBoost importance (gain, the share of the improvement in training fit attributable to",
  "each predictor, summing to one) alongside the pre-specified multivariable Firth odds",
  "ratios. Continuous predictors are per unit.")),
  "XGBoost importance and multivariable odds ratios", "ml_importance_vs_OR.csv", cmp)

cut_metrics <- function(thr, name) {
  yh <- as.integer(pred_te > thr)
  TP <- sum(yh == 1 & yte == 1); FP <- sum(yh == 1 & yte == 0)
  FN <- sum(yh == 0 & yte == 1); TN <- sum(yh == 0 & yte == 0)
  pct <- function(a, b) if (b > 0) sprintf("%.0f%%", 100 * a / b) else "-"
  data.frame(`Cut-off` = name, Threshold = if (is.finite(thr)) sprintf("%.3f", thr) else "-",
             Accuracy = pct(TP + TN, length(yte)), Sensitivity = pct(TP, TP + FN),
             Specificity = pct(TN, TN + FP), check.names = FALSE)
}
cuts <- c(Q25 = unname(quantile(pred_te, .25)), Median = median(pred_te), Mean = mean(pred_te),
          Q75 = unname(quantile(pred_te, .75)), Q90 = unname(quantile(pred_te, .90)))
cut_rows <- rbind(do.call(rbind, Map(cut_metrics, cuts, names(cuts))),
                  cut_metrics(Inf, "All negative (no-information rate)"))
add_tbl(style_ft(flextable(cut_rows)), "XGBoost classification by probability cut-off",
        "ml_cutoffs.csv", cut_rows)

# ---- Model figures ----------------------------------------------------------
roc_plot <- function(labels, scores, title, col = col_ns, fillcol = "#E8F1F8") {
  r   <- pROC::roc(labels, scores, quiet = TRUE)
  auc <- as.numeric(pROC::auc(r))
  ci  <- tryCatch(as.numeric(pROC::ci.auc(r))[c(1, 3)], error = function(e) c(NA, NA))
  df  <- data.frame(fpr = 1 - r$specificities, tpr = r$sensitivities)
  df  <- df[order(df$fpr, df$tpr), ]
  n   <- nrow(df)
  fill_df <- data.frame(fpr = c(as.vector(rbind(df$fpr, c(df$fpr[-1], df$fpr[n]))), 1, 0),
                        tpr = c(as.vector(rbind(df$tpr, df$tpr)), 0, 0))
  k   <- which.max(r$sensitivities + r$specificities - 1)
  opt <- data.frame(fpr = 1 - r$specificities[k], tpr = r$sensitivities[k],
                    lab = sprintf("Youden optimum\nSens %.2f  Spec %.2f",
                                  r$sensitivities[k], r$specificities[k]))
  auclab <- if (all(is.na(ci))) sprintf("AUC-ROC\n%.3f", auc) else
    sprintf("AUC-ROC\n%.3f (%.3f–%.3f)", auc, ci[1], ci[2])
  ggplot(df, aes(fpr, tpr)) +
    geom_polygon(data = fill_df, fill = fillcol, colour = NA) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = col_ref) +
    geom_step(colour = col, linewidth = 1.1, direction = "hv") +
    geom_point(data = opt, colour = col_sig, size = 3, shape = 1, stroke = 1.4) +
    geom_text(data = opt, aes(label = lab), hjust = -0.12, vjust = 0.4, size = 2.9,
              colour = col_sig, lineheight = 0.9) +
    annotate("label", x = 0.97, y = 0.06, hjust = 1, vjust = 0, size = 3.1, colour = col,
             fontface = "bold", label.size = 0.4, lineheight = 0.9, label = auclab) +
    scale_x_continuous("1 − specificity", labels = percent_format(accuracy = 1),
                       limits = c(0, 1), expand = c(0.004, 0)) +
    scale_y_continuous("Sensitivity", labels = percent_format(accuracy = 1),
                       limits = c(0, 1), expand = c(0.004, 0)) +
    coord_equal() + labs(title = title) + theme_thesis()
}

tryCatch(add_fig(roc_plot(m_pre$dat[[OUTCOME]], as.numeric(predict(m_pre$fit, type = "response")),
                          "Pre-specified multivariable model (apparent)", col = "#1f3b63"),
                 "ROC: multivariable model (apparent)", w = 6, h = 6),
         error = function(e) message("ROC skipped: ", conditionMessage(e)))

term_lab <- setNames(vapply(rows_pre, `[[`, "", 2), vapply(rows_pre, `[[`, "", 1))
imp_df <- data.frame(Lab = unname(term_lab),
                     Gain = ifelse(names(term_lab) %in% names(gain), gain[names(term_lab)], 0),
                     used = names(term_lab) %in% names(gain))
imp_df <- imp_df[order(imp_df$Gain), ]
imp_df$Lab <- factor(imp_df$Lab, levels = imp_df$Lab)
add_fig(ggplot(imp_df, aes(Gain, Lab, fill = used)) +
          geom_col(width = 0.68) +
          scale_fill_manual(values = c(`TRUE` = col_ns, `FALSE` = "#C9D6E5"), guide = "none") +
          scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
          labs(title = "XGBoost predictor importance",
               subtitle = "Pale bars: predictors not used by the model",
               x = "Importance (gain)", y = NULL) + theme_dsa(),
        "XGBoost importance")

# SHAP values (TreeSHAP) on the training set
tryCatch({
  relab <- function(v) ifelse(v %in% names(term_lab), term_lab[v], v)
  sc <- predict(mod, Xtr, predcontrib = TRUE)
  if (is.null(dim(sc))) sc <- matrix(sc, nrow = nrow(Xtr))
  if (is.null(colnames(sc))) {
    colnames(sc) <- if (ncol(sc) == ncol(Xtr) + 1) c(colnames(Xtr), "BIAS") else colnames(Xtr)
  }
  sc <- sc[, !colnames(sc) %in% c("BIAS", "(Intercept)"), drop = FALSE]
  fv <- as.data.frame(Xtr[, colnames(sc), drop = FALSE])
  n_train_ev <- sum(y_all[train])

  msh <- colMeans(abs(sc))
  sdf <- data.frame(Lab = relab(names(msh)), v = as.numeric(msh))
  sdf <- sdf[order(sdf$v), ]
  sdf$Lab <- factor(sdf$Lab, levels = sdf$Lab)
  add_fig(ggplot(sdf, aes(v, Lab)) + geom_col(fill = col_ns, width = 0.68) +
            scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
            labs(title = "SHAP importance (mean |SHAP|)",
                 subtitle = sprintf("Training set, %d events", n_train_ev),
                 x = "Mean |SHAP value|", y = NULL) + theme_dsa(),
          "SHAP importance", w = 8, h = 4.6)

  long <- do.call(rbind, lapply(colnames(sc), function(f) {
    x <- fv[[f]]; rg <- range(x, na.rm = TRUE)
    data.frame(feature = relab(f), shap = sc[, f],
               fval = if (diff(rg) > 0) (x - rg[1]) / diff(rg) else 0.5)
  }))
  long$feature <- factor(long$feature, levels = relab(names(sort(msh))))
  add_fig(ggplot(long, aes(shap, feature, colour = fval)) +
            geom_vline(xintercept = 0, colour = col_ref, linetype = "dashed") +
            geom_jitter(height = 0.18, width = 0, alpha = 0.65, size = 1.7) +
            scale_colour_gradient(low = "#2C7FB8", high = "#D55E00", breaks = c(0, 1),
                                  labels = c("low", "high"), name = "Feature\nvalue") +
            labs(title = "SHAP values by predictor",
                 subtitle = sprintf("Each point is one training recipient; %d events", n_train_ev),
                 x = "SHAP value", y = NULL) +
            theme_thesis() + theme(legend.position = "right"),
          "SHAP summary", w = 9, h = 5)
}, error = function(e) message("SHAP skipped: ", conditionMessage(e)))

add_fig(roc_plot(yte, pred_te, "Held-out test ROC: XGBoost", col = col_sig),
        "Test ROC: XGBoost", w = 6, h = 6)
add_fig(roc_plot(yte, firth_pte, "Held-out test ROC: penalised regression"),
        "Test ROC: Firth regression", w = 6, h = 6)

auc_ci <- function(l, s) tryCatch({
  ci <- as.numeric(pROC::ci.auc(pROC::roc(l, s, quiet = TRUE)))
  c(est = ci[2], lo = ci[1], hi = ci[3])
}, error = function(e) c(est = NA, lo = NA, hi = NA))
auc_df <- data.frame(model = factor(c("Firth penalised regression", "XGBoost"),
                                    levels = c("XGBoost", "Firth penalised regression")),
                     rbind(auc_ci(yte, firth_pte), auc_ci(yte, pred_te)))
add_fig(ggplot(auc_df, aes(est, model, colour = model)) +
          geom_vline(xintercept = 0.5, linetype = "dotted", colour = col_ref) +
          geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.16, linewidth = 0.9) +
          geom_point(size = 3.4) +
          geom_text(aes(label = sprintf("%.3f", est)), vjust = -1.1, size = 3.4, colour = "#222222") +
          scale_colour_manual(values = c("Firth penalised regression" = col_ns, "XGBoost" = col_sig),
                              guide = "none") +
          scale_x_continuous(limits = c(0.3, 1), breaks = seq(0.3, 1, 0.1)) +
          labs(title = "Held-out test AUC with 95% CI",
               subtitle = sprintf("DeLong CI; %d test recipients, %d events", length(yte), sum(yte)),
               x = "AUC-ROC", y = NULL) + theme_dsa(),
        "Test AUC with confidence intervals", w = 8, h = 3.2)

perm_long <- rbind(data.frame(auc = perm$train_auc, set = "Training"),
                   data.frame(auc = perm$test_auc, set = "Held-out test"))
add_fig(ggplot(perm_long, aes(auc, set)) +
          geom_vline(xintercept = 0.5, linetype = "dotted", colour = col_ref) +
          geom_jitter(height = 0.14, width = 0, size = 2.2, colour = col_ns, alpha = 0.55) +
          geom_vline(xintercept = auc_te, colour = col_sig, linewidth = 1) +
          annotate("text", x = auc_te, y = 2.4, hjust = -0.05, size = 3.3, colour = col_sig,
                   label = sprintf("observed test AUC %.3f", auc_te)) +
          scale_x_continuous(limits = c(0, 1)) +
          labs(title = "XGBoost performance with permuted outcome labels", x = "AUC", y = NULL,
               caption = sprintf("%d permutations; permutation p = %.3f", nrow(perm), p_perm)) +
          theme_dsa(),
        "Label-permutation test")

prob_df <- rbind(data.frame(p = pred_tr, set = "Training set"),
                 data.frame(p = pred_te, set = "Held-out test set"))
prob_df$set <- factor(prob_df$set, levels = c("Training set", "Held-out test set"))
add_fig(ggplot(prob_df, aes(p)) +
          geom_histogram(bins = 20, fill = col_ns, colour = "white", linewidth = 0.3) +
          geom_vline(xintercept = mean(y_all), linetype = "dashed", colour = col_ref) +
          facet_wrap(~ set, ncol = 2) +
          labs(title = "Predicted probability of de novo DSA",
               subtitle = "Dashed line: cohort event rate", x = "P(dnDSA+)", y = "Recipients") +
          theme_dsa(),
        "Predicted probabilities", w = 9, h = 4.4)

# ---- Write outputs ----------------------------------------------------------
print(doc,  target = op("dsa_prespecified_tables.docx"))
print(deck, target = op("dsa_prespecified_figures.pptx"))

cat(sprintf("\nn=%d, events=%d\n", n_tot, n_ev))
cat(sprintf("Test AUC: Firth %.3f, XGBoost %.3f; permutation p = %.3f\n", firth_auc, auc_te, p_perm))
cat(sprintf("Univariable: raw p<0.05 = %d, BH p<0.05 = %d\n",
            sum(uni$p < 0.05, na.rm = TRUE), sum(uni$p_BH < 0.05, na.rm = TRUE)))
cat(sprintf("Pre-specified predictors with BH p<0.05: %d of %d\n",
            sum(con$p_BH < 0.05, na.rm = TRUE), nrow(con)))
