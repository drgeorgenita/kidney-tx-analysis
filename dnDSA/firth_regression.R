# De novo DSA: univariable and multivariable Firth-penalised logistic regression

source("R/helpers.R")
install_missing(c("readxl", "logistf", "patchwork", "dplyr", "officer", "rvg",
                  "pROC", "scales", "flextable"))
suppressPackageStartupMessages({
  library(readxl); library(logistf); library(patchwork); library(dplyr)
  library(officer); library(rvg); library(pROC); library(scales); library(flextable)
})

in_file <- "data/dsa_cohort.xlsx"
out <- "output/dnDSA"
dir.create(out, showWarnings = FALSE, recursive = TRUE)

OUT <- "de_novoDSA"
d <- read_excel(in_file)
d[[OUT]] <- as.numeric(as.character(d[[OUT]]))
d <- d[!is.na(d[[OUT]]), , drop = FALSE]
n_tot <- nrow(d); n_ev <- sum(d[[OUT]] == 1)
message("n = ", n_tot, " (", n_ev, " dnDSA+)")

# ---- Predictors -------------------------------------------------------------
# "cont": per-unit effect (binary 0/1 variables coded as numbers are also per-unit)
# "cat":  factor with an explicit reference level
spec <- list(
  "Demographics" = list(
    list(var = "age_at_tx",       type = "cont", label = "Age at transplant (per year)"),
    list(var = "age_group",       type = "cat",  ref = "<50", levels = c("50-59", "60-69", "70-79"), label = "Age group"),
    list(var = "over_60",         type = "cat",  ref = "0",   levels = "1", label = "Age >=60 years"),
    list(var = "non_white",       type = "cat",  ref = "0",   levels = "1", label = "Non-white ethnicity"),
    list(var = "rec_sex",         type = "cat",  ref = "F",   levels = "M", label = "Recipient sex"),
    list(var = "dialysis_status", type = "cat",  ref = "HD",  levels = c("PD", "pre"), label = "Dialysis status")
  ),
  "Comorbidities and risk indices" = list(
    list(var = "rec_BMI",                    type = "cont", label = "Recipient BMI (per kg/m2)"),
    list(var = "rec_smoking_status",         type = "cont", label = "Smoking status (per category)"),
    list(var = "rec_HTN",                    type = "cont", label = "Hypertension (yes vs no)"),
    list(var = "rec_diabetes",               type = "cont", label = "Diabetes (yes vs no)"),
    list(var = "age_CCI",                    type = "cont", label = "Age component of CCI (per point)"),
    list(var = "charlson_comorbidity_index", type = "cont", label = "Charlson Comorbidity Index (per point)"),
    list(var = "ASA",                        type = "cont", label = "ASA grade (per grade)"),
    list(var = "WHO_PS",                     type = "cont", label = "WHO performance status (per point)"),
    list(var = "UKKRRI",                     type = "cont", label = "UKKRRI (per unit)"),
    list(var = "UKKRRI_quartile",            type = "cat",  ref = "R1", levels = c("R2", "R3", "R4"), label = "UKKRRI quartile")
  ),
  "Donor and transplant characteristics" = list(
    list(var = "HLA_MM_level", type = "cont", label = "HLA mismatch level (per level)"),
    list(var = "donor_age",    type = "cont", label = "Donor age (per year)"),
    list(var = "donor_sex",    type = "cat",  ref = "F",   levels = "M",   label = "Donor sex"),
    list(var = "tx_categ",     type = "cat",  ref = "0",   levels = "1",   label = "Transplant category (deceased vs live)"),
    list(var = "donor_rhesus", type = "cat",  ref = "neg", levels = "pos", label = "Donor rhesus"),
    list(var = "rec_rhesus",   type = "cat",  ref = "neg", levels = "pos", label = "Recipient rhesus")
  ),
  "Viral status and induction immunosuppression" = list(
    list(var = "donor_CMV",     type = "cat",  ref = "0", levels = "1", label = "Donor CMV seropositive"),
    list(var = "recipient_CMV", type = "cat",  ref = "0", levels = "1", label = "Recipient CMV seropositive"),
    list(var = "donor_EBV",     type = "cat",  ref = "0", levels = "1", label = "Donor EBV seropositive"),
    list(var = "recipient_EBV", type = "cat",  ref = "0", levels = "1", label = "Recipient EBV seropositive"),
    list(var = "induction_IS",  type = "cat",  ref = "Alemtuzumab1", levels = c("Alemtuzumab2", "Basiliximab"), label = "Induction regimen"),
    list(var = "CMV_log",       type = "cont", label = "CMV viral load (log, per unit)"),
    list(var = "EBV_log",       type = "cont", label = "EBV viral load (log, per unit)")
  ),
  "Crossmatch, post-transplant viraemia, prior sensitisation" = list(
    list(var = "flow_cytoT",        type = "cat", ref = "negative", levels = "Non-negative", label = "Flow cytometry T-cell XM"),
    list(var = "flow_cytoB",        type = "cat", ref = "negative", levels = "positive",     label = "Flow cytometry B-cell XM"),
    list(var = "EBV_viraemia",      type = "cat", ref = "0", levels = "1", label = "EBV viraemia"),
    list(var = "CMV_viraemia",      type = "cat", ref = "0", levels = "1", label = "CMV viraemia"),
    list(var = "Adeno_viraemia",    type = "cat", ref = "0", levels = "1", label = "Adenovirus viraemia"),
    list(var = "pregnancy_history", type = "cat", ref = "0", levels = "1", label = "Pregnancy history"),
    list(var = "prev_transfusions", type = "cat", ref = "0", levels = "1", label = "Previous transfusions")
  )
)

# complete-case data for one predictor, coded as specified
prep_sub <- function(s) {
  cc  <- !is.na(d[[s$var]]) & !is.na(d[[OUT]])
  sub <- data.frame(y = d[[OUT]][cc], x = d[[s$var]][cc])
  if (s$type == "cont") {
    sub$x <- as.numeric(as.character(sub$x))
  } else {
    sub$x <- factor(as.character(sub$x), levels = c(s$ref, s$levels))
    if (any(is.na(sub$x)))
      warning(sprintf("%s: %d rows with an undeclared level dropped", s$var, sum(is.na(sub$x))))
    sub <- sub[!is.na(sub$x), ]
  }
  sub
}

fit_one <- function(s, group) {
  sub <- prep_sub(s)
  fit <- logistf(y ~ x, data = sub)
  row <- function(level, level_label, is_ref, tn, sel) {
    data.frame(group, var = s$var, label = s$label, level = level, level_label = level_label,
               is_ref = is_ref,
               OR = if (is_ref) NA else exp(fit$coefficients[[tn]]),
               lo = if (is_ref) NA else exp(fit$ci.lower[[tn]]),
               hi = if (is_ref) NA else exp(fit$ci.upper[[tn]]),
               p  = if (is_ref) NA else fit$prob[[tn]],
               ev = sum(sub$y[sel]), N = sum(sel), stringsAsFactors = FALSE)
  }
  if (s$type == "cont") return(row(NA, s$label, FALSE, "x", rep(TRUE, nrow(sub))))
  rbind(row(s$ref, paste0(s$ref, " (ref)"), TRUE, NA, sub$x == s$ref),
        do.call(rbind, lapply(s$levels, function(lv)
          row(lv, lv, FALSE, paste0("x", lv), sub$x == lv))))
}

# AUC with DeLong 95% CI for the univariable model of each predictor
get_auc <- function(s) {
  na_row <- data.frame(var = s$var, AUC = NA_real_, AUC_lo = NA_real_, AUC_hi = NA_real_)
  tryCatch({
    sub <- prep_sub(s)
    fit <- logistf(y ~ x, data = sub)
    roc_obj <- pROC::roc(sub$y, fit$predict, quiet = TRUE)
    ci <- as.numeric(pROC::ci.auc(roc_obj, method = "delong"))
    data.frame(var = s$var, AUC = ci[2], AUC_lo = ci[1], AUC_hi = ci[3])
  }, error = function(e) na_row)
}

res <- do.call(rbind, lapply(names(spec), function(g)
  do.call(rbind, lapply(spec[[g]], fit_one, group = g))))
auc_df <- do.call(rbind, lapply(unlist(spec, recursive = FALSE), get_auc))
res <- cbind(res, auc_df[match(res$var, auc_df$var), c("AUC", "AUC_lo", "AUC_hi")])
res$AUC_str <- ifelse(is.na(res$AUC), "",
                      sprintf("%.3f (%.3f–%.3f)", res$AUC, res$AUC_lo, res$AUC_hi))
rownames(res) <- NULL
write.csv(res, file.path(out, "firth_univariable.csv"), row.names = FALSE)

# ---- Multivariable models ---------------------------------------------------
# full: all predictors with univariable p <= 0.10
# reduced: after backward selection
CAT_LEVELS <- list(
  tx_categ       = c("0", "1"),
  recipient_CMV  = c("0", "1"),
  recipient_EBV  = c("0", "1"),
  EBV_viraemia   = c("0", "1"),
  Adeno_viraemia = c("0", "1"),
  induction_IS   = c("Alemtuzumab1", "Alemtuzumab2", "Basiliximab"))

fit_mv <- function(rhs, maxit = NULL) {
  vars <- all.vars(as.formula(paste("~", rhs)))
  dat  <- d[, c(OUT, vars), drop = FALSE]
  for (v in intersect(vars, names(CAT_LEVELS)))
    dat[[v]] <- factor(as.character(dat[[v]]), levels = CAT_LEVELS[[v]])
  dat <- dat[stats::complete.cases(dat), , drop = FALSE]
  args <- list(formula = as.formula(paste(OUT, "~", rhs)), data = dat)
  if (!is.null(maxit)) {
    args$control   <- logistf.control(maxit = maxit)
    args$plcontrol <- logistpl.control(maxit = maxit)
  }
  list(fit = do.call(logistf, args), dat = dat)
}

mv_row <- function(m, var, lv, disp, group) {
  f <- m$fit; dat <- m$dat
  if (is.na(lv)) {
    tn <- var; ev <- sum(dat[[OUT]]); N <- nrow(dat)
  } else {
    tn <- paste0(var, lv); sel <- dat[[var]] == lv
    ev <- sum(dat[[OUT]][sel]); N <- sum(sel)
  }
  data.frame(group = group, var = paste0(var, "_", ifelse(is.na(lv), "c", lv)),
             label = disp, level = NA, level_label = disp, is_ref = FALSE,
             OR = exp(f$coefficients[[tn]]), lo = exp(f$ci.lower[[tn]]),
             hi = exp(f$ci.upper[[tn]]), p = f$prob[[tn]],
             ev = ev, N = N, stringsAsFactors = FALSE)
}

m_red <- fit_mv("tx_categ + recipient_CMV + induction_IS")
grp_red <- sprintf("Final multivariable model (Firth, n=%d)", nrow(m_red$dat))
mv_red <- rbind(
  mv_row(m_red, "tx_categ",      "1",            "Deceased vs live donor",        grp_red),
  mv_row(m_red, "recipient_CMV", "1",            "Recipient CMV+ vs CMV−",   grp_red),
  mv_row(m_red, "induction_IS",  "Alemtuzumab2", "Alemtuzumab x2 vs x1",          grp_red),
  mv_row(m_red, "induction_IS",  "Basiliximab",  "Basiliximab vs Alemtuzumab x1", grp_red))

m_full <- fit_mv(paste("tx_categ + recipient_CMV + recipient_EBV + induction_IS +",
                       "EBV_viraemia + Adeno_viraemia + ASA"), maxit = 1000)
grp_full <- sprintf("Full multivariable model (univariable p≤0.10, Firth, n=%d)", nrow(m_full$dat))
mv_full <- rbind(
  mv_row(m_full, "tx_categ",       "1",            "Deceased vs live donor",        grp_full),
  mv_row(m_full, "recipient_CMV",  "1",            "Recipient CMV+ vs CMV−",   grp_full),
  mv_row(m_full, "recipient_EBV",  "1",            "Recipient EBV+ vs EBV−",   grp_full),
  mv_row(m_full, "induction_IS",   "Alemtuzumab2", "Alemtuzumab x2 vs x1",          grp_full),
  mv_row(m_full, "induction_IS",   "Basiliximab",  "Basiliximab vs Alemtuzumab x1", grp_full),
  mv_row(m_full, "EBV_viraemia",   "1",            "EBV viraemia",                  grp_full),
  mv_row(m_full, "Adeno_viraemia", "1",            "Adenovirus viraemia",           grp_full),
  mv_row(m_full, "ASA",            NA,             "ASA grade (per grade)",         grp_full))

write.csv(rbind(mv_full, mv_red), file.path(out, "firth_multivariable.csv"), row.names = FALSE)

mv_auc <- function(m) {
  tryCatch({
    roc_obj <- pROC::roc(as.numeric(m$dat[[OUT]]), m$fit$predict, quiet = TRUE, direction = "<")
    ci <- as.numeric(pROC::ci.auc(roc_obj, method = "delong"))
    sprintf("%.3f (%.3f–%.3f)", ci[2], ci[1], ci[3])
  }, error = function(e) NA_character_)
}
auc_red  <- mv_auc(m_red)
auc_full <- mv_auc(m_full)
message("AUC reduced model: ", auc_red, " | full model: ", auc_full)

# ---- Formatting -------------------------------------------------------------
fmt_bound <- function(x)
  ifelse(is.na(x), "Inf",
  ifelse(abs(x) >= 100, sprintf("%.0f", x),
  ifelse(abs(x) > 0 & abs(x) < 0.01, sprintf("%.3f", x), sprintf("%.2f", x))))
or_string <- function(OR, lo, hi)
  ifelse(is.na(OR), "", sprintf("%.2f (%s to %s)", OR, fmt_bound(lo), fmt_bound(hi)))
p_string <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001",
                         ifelse(p < 0.10, sprintf("%.3f", p), sprintf("%.2f", p))))

# ---- Forest plots -----------------------------------------------------------
XLO <- 0.02; XHI <- 80
BREAKS <- c(0.05, 0.1, 0.25, 0.5, 1, 2, 4, 10, 40)
row_cols <- c("disp", "is_header", "is_ref", "OR", "lo", "hi", "p", "ev", "N")

# one row per continuous predictor; a header row plus one row per level for categorical
build_rows <- function(df) {
  out_rows <- lapply(unique(df$var), function(v) {
    vr <- df[df$var == v, ]
    if (nrow(vr) == 1 && is.na(vr$level)) {
      vr$disp <- vr$label; vr$is_header <- FALSE
      return(vr[, row_cols])
    }
    hdr <- data.frame(disp = vr$label[1], is_header = TRUE, is_ref = FALSE,
                      OR = NA, lo = NA, hi = NA, p = NA, ev = NA, N = NA)
    vr$disp <- paste0(" ", vr$level_label); vr$is_header <- FALSE
    rbind(hdr, vr[, row_cols])
  })
  r <- do.call(rbind, out_rows)
  r$ord <- seq_len(nrow(r))
  r$yf  <- factor(r$ord, levels = rev(r$ord))
  r
}

forest_plot <- function(df) {
  r <- build_rows(df)
  sig <- !is.na(r$p) & r$p < 0.05
  r$col <- ifelse(sig, col_sig, col_ns)
  pts <- r[!r$is_header & !is.na(r$OR), ]
  pts$lo_c <- pmax(pts$lo, XLO * 1.05)
  pts$hi_c <- pmin(pts$hi, XHI * 0.95)
  ylabs <- setNames(r$disp, as.character(r$yf))

  pf <- ggplot(pts) +
    geom_vline(xintercept = 1, linetype = "dashed", color = col_ref, linewidth = 0.5) +
    geom_errorbarh(aes(y = yf, xmin = lo_c, xmax = hi_c, color = col), height = 0, linewidth = 0.7) +
    geom_point(aes(y = yf, x = OR, color = col, fill = col), shape = 21, size = 3, stroke = 1.1)

  # arrows where the interval extends beyond the axis
  hi_off <- pts[!is.na(pts$hi) & pts$hi > XHI * 0.95, ]
  lo_off <- pts[!is.na(pts$lo) & pts$lo < XLO * 1.05, ]
  if (nrow(hi_off)) pf <- pf + geom_segment(data = hi_off,
    aes(y = yf, yend = yf, x = hi_c * 0.9, xend = hi_c, color = col),
    arrow = arrow(length = unit(0.07, "cm")), linewidth = 0.7)
  if (nrow(lo_off)) pf <- pf + geom_segment(data = lo_off,
    aes(y = yf, yend = yf, x = lo_c * 1.1, xend = lo_c, color = col),
    arrow = arrow(length = unit(0.07, "cm")), linewidth = 0.7)

  pf <- pf +
    scale_color_identity() + scale_fill_identity() +
    scale_x_log10(limits = c(XLO, XHI), breaks = BREAKS, labels = as.character(BREAKS)) +
    scale_y_discrete(labels = ylabs, drop = FALSE) +
    labs(x = "Odds ratio for de novo DSA (log scale, Firth-penalised)", y = NULL,
         caption = "Red, p<0.05; blue, p≥0.05") +
    theme_thesis() +
    theme(panel.grid.major.y = element_blank(),
          axis.text.y = element_text(hjust = 0), axis.ticks.y = element_blank(),
          plot.caption = element_text(hjust = 0, color = col_grey, size = 8))

  r$or_txt  <- ifelse(r$is_header, "", or_string(r$OR, r$lo, r$hi))
  r$evn_txt <- ifelse(r$is_header | is.na(r$N), "", paste0(r$ev, "/", r$N))
  r$txtcol  <- ifelse(is.na(r$OR), "#222222", r$col)
  r$fface   <- ifelse(sig, "bold", "plain")

  pt <- ggplot(r) +
    geom_text(aes(x = 0, y = yf, label = or_txt, color = txtcol, fontface = fface),
              hjust = 0, size = 3) +
    geom_text(aes(x = 1.05, y = yf, label = evn_txt), hjust = 0, size = 3, color = col_grey) +
    scale_color_identity() +
    scale_x_continuous(limits = c(0, 1.8)) +
    scale_y_discrete(drop = FALSE) +
    labs(title = "OR (95% CI)   ev/N") +
    theme_void(base_size = 11, base_family = "Arial") +
    theme(plot.title = element_text(size = 9, face = "bold", hjust = 0))

  pf + pt + plot_layout(widths = c(3, 1.6))
}

# ---- Figures ----------------------------------------------------------------
figs <- list()
addfig <- function(title, plot, w, h)
  figs[[length(figs) + 1]] <<- list(title = title, plot = plot, w = w, h = h)
fig_h <- function(df) 0.42 * nrow(build_rows(df)) + 2.2

# split a predictor group across slides at variable boundaries
ROWS_PER_SLIDE <- 13
pack_vars <- function(df, budget) {
  vars <- unique(df$var)
  rc <- vapply(vars, function(v) nrow(build_rows(df[df$var == v, , drop = FALSE])), integer(1))
  k <- max(1, ceiling(sum(rc) / budget))
  target <- sum(rc) / k
  chunks <- list(); cur <- character(0); cur_n <- 0
  for (i in seq_along(vars)) {
    if (cur_n > 0 && length(chunks) < k - 1 && (cur_n + rc[i] - target) > (target - cur_n)) {
      chunks[[length(chunks) + 1]] <- cur; cur <- character(0); cur_n <- 0
    }
    cur <- c(cur, vars[i]); cur_n <- cur_n + rc[i]
  }
  if (length(cur)) chunks[[length(chunks) + 1]] <- cur
  chunks
}
add_domain <- function(df, base_title) {
  ch <- pack_vars(df, ROWS_PER_SLIDE)
  for (i in seq_along(ch)) {
    sub <- df[df$var %in% ch[[i]], , drop = FALSE]
    ttl <- if (length(ch) > 1) sprintf("%s (%d/%d)", base_title, i, length(ch)) else base_title
    addfig(ttl, forest_plot(sub), 11, fig_h(sub))
  }
}

for (g in names(spec)) add_domain(res[res$group == g, , drop = FALSE], g)

key_spec <- list(
  c("recipient_CMV", "1",            "Recipient CMV seropositive"),
  c("recipient_EBV", "1",            "Recipient EBV seropositive"),
  c("induction_IS",  "Alemtuzumab2", "Alemtuzumab x2 vs x1"),
  c("tx_categ",      "1",            "Deceased vs live donor"),
  c("CMV_viraemia",  "1",            "CMV viraemia"),
  c("EBV_viraemia",  "1",            "EBV viraemia"),
  c("over_60",       "1",            "Recipient age >=60"))
key <- do.call(rbind, lapply(key_spec, function(k) {
  row <- res[res$var == k[1] & !is.na(res$level) & res$level == k[2], ]
  transform(row, label = k[3], level = NA, level_label = k[3])
}))
addfig("Strongest univariable associations", forest_plot(key), 11, fig_h(key))

roc_plot <- function(m, title, auc_str) {
  tryCatch({
    outcome <- as.numeric(m$dat[[OUT]])
    roc_obj <- pROC::roc(outcome, m$fit$predict, quiet = TRUE, direction = "<")
    roc_df  <- data.frame(spec = rev(roc_obj$specificities), sens = rev(roc_obj$sensitivities))
    best <- which.max(roc_obj$sensitivities + roc_obj$specificities - 1)
    best_pt <- data.frame(spec = roc_obj$specificities[best], sens = roc_obj$sensitivities[best])

    ggplot(roc_df, aes(x = 1 - spec, y = sens)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = col_ref, linewidth = 0.5) +
      geom_ribbon(aes(ymin = 0, ymax = sens), fill = col_ns, alpha = 0.08) +
      geom_step(direction = "hv", colour = col_ns, linewidth = 1.2) +
      geom_point(data = best_pt, colour = col_sig, size = 3.5, shape = 21,
                 fill = "white", stroke = 1.6) +
      annotate("text", x = 1 - best_pt$spec + 0.03, y = best_pt$sens - 0.04,
               label = sprintf("Youden optimum\nSens %.2f  Spec %.2f", best_pt$sens, best_pt$spec),
               hjust = 0, size = 3.0, colour = col_sig) +
      annotate("label", x = 0.62, y = 0.18, hjust = 0, size = 3.4, fontface = "bold",
               label = sprintf("AUC-ROC\n%s", ifelse(is.na(auc_str), "not available", auc_str)),
               colour = col_ns, fill = "#F0F6FF", label.size = 0.3,
               label.padding = unit(0.4, "lines")) +
      scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                         labels = percent_format(accuracy = 1), expand = c(0.01, 0.01)) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                         labels = percent_format(accuracy = 1), expand = c(0.01, 0.01)) +
      labs(title = title,
           subtitle = sprintf("Firth logistic regression, n = %d, %d events",
                              length(outcome), sum(outcome)),
           x = "1 − specificity", y = "Sensitivity") +
      theme_thesis() +
      theme(aspect.ratio = 1)
  }, error = function(e) NULL)
}

addfig(sprintf("%s: AUC %s", grp_full, auc_full), forest_plot(mv_full), 11, fig_h(mv_full))
p_roc_full <- roc_plot(m_full, "ROC: full multivariable model", auc_full)
if (!is.null(p_roc_full)) addfig("ROC: full multivariable model", p_roc_full, 7, 6.5)

addfig(sprintf("%s: AUC %s", grp_red, auc_red), forest_plot(mv_red), 11, fig_h(mv_red))
p_roc_red <- roc_plot(m_red, "ROC: final multivariable model", auc_red)
if (!is.null(p_roc_red)) addfig("ROC: final multivariable model", p_roc_red, 7, 6.5)

# dnDSA incidence by recipient CMV, donor type and induction
panel_spec <- list(
  list(var = "recipient_CMV", map = c("0" = "CMV−", "1" = "CMV+"), title = "Recipient CMV serostatus"),
  list(var = "tx_categ",      map = c("0" = "Live", "1" = "Deceased"),  title = "Transplant category"),
  list(var = "induction_IS",  map = c("Alemtuzumab1" = "Alem x1", "Alemtuzumab2" = "Alem x2",
                                      "Basiliximab" = "Basilix."),    title = "Induction regimen"))
bar_cols <- c("#0072B2", "#009E73", "#E69F00", "#CC79A7")
panels <- lapply(panel_spec, function(s) {
  cc <- !is.na(d[[s$var]]) & !is.na(d[[OUT]])
  x <- as.character(d[[s$var]][cc]); y <- d[[OUT]][cc]
  lev <- names(s$map)
  pd <- data.frame(cat = factor(s$map[lev], levels = s$map[lev]),
                   rate = sapply(lev, function(l) 100 * mean(y[x == l])),
                   ev = sapply(lev, function(l) sum(y[x == l])),
                   n = sapply(lev, function(l) sum(x == l)))
  ggplot(pd, aes(cat, rate, fill = cat)) +
    geom_col(width = 0.62, color = "white") +
    geom_text(aes(label = sprintf("%d/%d\n%.1f%%", ev, n, rate)), vjust = -0.2, size = 3) +
    scale_fill_manual(values = bar_cols[seq_along(lev)], guide = "none") +
    labs(x = NULL, y = "dnDSA incidence (%)", title = s$title) +
    expand_limits(y = max(pd$rate) * 1.35 + 3) +
    theme_thesis() +
    theme(panel.grid.major.x = element_blank())
})
addfig("De novo DSA incidence by recipient CMV, donor type and induction",
       wrap_plots(panels, nrow = 1), 13, 4.6)

# ---- dnDSA trajectories -----------------------------------------------------
# Peak MFI per visit is the highest bead across that visit's class I and class II
# results. Antibody labels list the specificities without MFI values.
clean_txt <- function(x) {
  x <- as.character(x)
  ifelse(is.na(x) | x %in% c("", "NA", "NULL", "na"), NA, x)
}
parse_mfi <- function(x) {
  x <- clean_txt(x)
  if (is.na(x)) return(numeric(0))
  as.numeric(unlist(regmatches(x, gregexpr("[0-9]{3,6}", x))))
}
antibody_label <- function(s) {
  s <- clean_txt(s)
  if (is.na(s)) return(NA_character_)
  s <- gsub("\\(\\s*cleared\\s*\\)", "", s, ignore.case = TRUE)
  s <- gsub("\\(\\s*[0-9]{3,6}[^)]*\\)", "|", s)
  u <- trimws(strsplit(s, "\\|")[[1]]); u <- u[u != ""]
  u <- gsub("DRB4\\*01:01 \\(DR53\\)", "DR53 (DRB4*01:01)", u)
  u <- gsub("^DQA1\\*05:01$", "DQ (DQA1*05:01)", u)
  paste(u, collapse = ", ")
}

dsa_pos <- d[as.character(d[[OUT]]) %in% c("1", "1.0"), ]
visits <- data.frame(
  label = c("2 weeks", "1 month", "3 months", "6 months", "12 months"), x = 1:5,
  c1 = paste0("DSA", 1:5, "_class_1(MFI)"),
  c2 = paste0("DSA", 1:5, "_class_2(MFI)"), stringsAsFactors = FALSE)
class_map <- c("1" = "I", "2" = "II", "3" = "I+II")
min_mfi <- 0

traj <- data.frame(); meta <- data.frame()
for (ri in seq_len(nrow(dsa_pos))) {
  rr <- dsa_pos[ri, ]
  pts <- do.call(rbind, lapply(seq_len(nrow(visits)), function(t) {
    m <- c(parse_mfi(rr[[visits$c1[t]]]), parse_mfi(rr[[visits$c2[t]]]))
    if (length(m)) data.frame(x = visits$x[t], mfi = max(m)) else NULL
  }))
  if (is.null(pts) || max(pts$mfi) < min_mfi) next
  pts$id <- as.character(rr$ID)
  traj <- rbind(traj, pts)
  meta <- rbind(meta, data.frame(id = as.character(rr$ID), label = antibody_label(rr$dnDSA_MFI),
                                 grp = class_map[[as.character(rr$dnDSA_class)]],
                                 first = min(pts$x), peak = max(pts$mfi)))
}

meta <- meta[order(match(meta$grp, c("I", "II", "I+II")), meta$first, -meta$peak), ]
meta$id <- factor(meta$id, levels = rev(meta$id))
traj$id <- factor(traj$id, levels = levels(meta$id))
traj <- merge(traj, meta[, c("id", "grp")], by = "id")
traj$grp <- factor(traj$grp, levels = c("I", "II", "I+II"))

p_traj <- ggplot(traj, aes(x, id)) +
  geom_line(aes(group = id, colour = grp), linewidth = 0.9, alpha = 0.5) +
  geom_point(aes(size = mfi, colour = grp), alpha = 0.95) +
  geom_text(aes(label = mfi), vjust = -1.6, size = 2.5, colour = "grey25") +
  scale_x_continuous(breaks = 1:5, labels = visits$label, limits = c(0.6, 5.4)) +
  scale_y_discrete(labels = setNames(meta$label, as.character(meta$id))) +
  scale_size_continuous(range = c(2, 13), name = "Peak MFI", breaks = c(2000, 6000, 12000)) +
  scale_colour_manual(values = c("I" = "#2c7fb8", "II" = "#1b9e77", "I+II" = "#d95f02"),
                      name = "HLA class", labels = c("Class I", "Class II", "Class I + II")) +
  labs(x = "Time since transplant", y = NULL) +
  theme_minimal(base_size = 11, base_family = "Arial") +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "right")

figs <- c(list(list(title = "De novo DSA trajectories by recipient", plot = p_traj,
                    w = 10.5, h = 7.2)), figs)

# ---- Write slides -----------------------------------------------------------
# each figure is scaled to fit the slide at its own aspect ratio
deck <- read_pptx()
sz <- slide_size(deck)
MARGIN <- 0.3; TITLE_BAND <- 1.05
for (f in figs) {
  deck <- add_slide(deck, layout = "Title Only", master = "Office Theme")
  deck <- ph_with(deck, value = f$title, location = ph_location_type(type = "title"))
  avail_w <- sz$width - 2 * MARGIN
  avail_h <- sz$height - TITLE_BAND - MARGIN
  a <- f$w / f$h
  if (avail_w / avail_h > a) { h <- avail_h; w <- h * a } else { w <- avail_w; h <- w / a }
  loc <- ph_location(left = (sz$width - w) / 2, top = TITLE_BAND + (avail_h - h) / 2,
                     width = w, height = h)
  deck <- tryCatch(ph_with(deck, value = dml(ggobj = f$plot), location = loc),
                   error = function(e) {
                     tmp <- tempfile(fileext = ".png")
                     ggsave(tmp, f$plot, width = w, height = h, units = "in", dpi = 200, bg = "white")
                     ph_with(deck, value = external_img(tmp, width = w, height = h), location = loc)
                   })
}
print(deck, target = file.path(out, "firth_figures.pptx"))

# ---- Word tables ------------------------------------------------------------
build_table <- function(df) {
  rows <- lapply(unique(df$var), function(v) {
    vr <- df[df$var == v, , drop = FALSE]
    auc_val <- ifelse(is.na(vr$AUC_str[1]), "", vr$AUC_str[1])
    mk <- function(pred, N, ev, or, p, auc, sig)
      data.frame(Predictor = pred, N = N, `Event N` = ev, `OR (95% CI)` = or,
                 `p-value` = p, `AUC (95% CI)` = auc, sig = sig,
                 check.names = FALSE, stringsAsFactors = FALSE)
    if (nrow(vr) == 1 && is.na(vr$level[1])) {
      return(mk(v, as.character(vr$N), as.character(vr$ev), or_string(vr$OR, vr$lo, vr$hi),
                p_string(vr$p), auc_val, !is.na(vr$p) & vr$p < 0.05))
    }
    rbind(mk(v, "", "", "", "", auc_val, FALSE),
          do.call(rbind, lapply(seq_len(nrow(vr)), function(i) {
            rw <- vr[i, ]
            mk(paste0(" ", rw$level), as.character(rw$N), as.character(rw$ev),
               if (isTRUE(rw$is_ref)) "-" else or_string(rw$OR, rw$lo, rw$hi),
               if (isTRUE(rw$is_ref)) "" else p_string(rw$p), "",
               !is.na(rw$p) & rw$p < 0.05)
          })))
  })
  do.call(rbind, rows)
}

doc <- read_docx()
for (g in names(spec)) {
  tb <- build_table(res[res$group == g, , drop = FALSE])
  sig <- which(tb$sig); tb$sig <- NULL
  ft <- flextable(tb) |> bold(part = "header") |>
    set_table_properties(layout = "autofit", width = 1)
  if (length(sig)) ft <- bold(ft, i = sig, part = "body")
  doc <- doc |> body_add_par(g, style = "heading 1") |> body_add_flextable(ft) |> body_add_par("")
}

add_mv_table <- function(doc, mvt, heading, auc_str) {
  tb <- data.frame(Predictor = c(mvt$label, "Model AUC-ROC (DeLong 95% CI)"),
                   `OR (95% CI)` = c(or_string(mvt$OR, mvt$lo, mvt$hi),
                                     ifelse(is.na(auc_str), "not computed", auc_str)),
                   `p-value` = c(p_string(mvt$p), ""),
                   check.names = FALSE, stringsAsFactors = FALSE)
  n <- nrow(tb)
  sig <- which(!is.na(mvt$p) & mvt$p < 0.05)
  ft <- flextable(tb) |> bold(part = "header") |>
    bold(i = n, part = "body") |>
    bg(i = n, bg = "#EBF3FB", part = "body") |>
    hline(i = n - 1, border = fp_border(width = 1, color = "#333333"), part = "body") |>
    set_table_properties(layout = "autofit", width = 1)
  if (length(sig)) ft <- bold(ft, i = sig, part = "body")
  doc |> body_add_par(heading, style = "heading 1") |> body_add_flextable(ft) |> body_add_par("")
}
doc <- add_mv_table(doc, mv_full, grp_full, auc_full)
doc <- add_mv_table(doc, mv_red, grp_red, auc_red)

# dnDSA cases
relabel_mfi_str <- function(s) {
  if (is.na(s)) return("-")
  s <- gsub("DRB4\\*01:01",  "DR53 (DRB4*01:01)", s)
  s <- gsub("DRB3\\*(\\S+)", "DR52 (DRB3*\\1)", s, perl = TRUE)
  s <- gsub("DRB5\\*(\\S+)", "DR51 (DRB5*\\1)", s, perl = TRUE)
  s <- gsub("DQA1\\*(\\S+)", "DQ (DQA1*\\1)", s, perl = TRUE)
  s <- gsub("DQB1\\*(\\S+)", "DQ (DQB1*\\1)", s, perl = TRUE)
  s <- gsub("DPA1\\*(\\S+)", "DP (DPA1*\\1)", s, perl = TRUE)
  s <- gsub("DPB1\\*(\\S+)", "DP (DPB1*\\1)", s, perl = TRUE)
  s
}
c4d_any <- function(row) any(c(row[["bx1_C4d"]], row[["bx2_C4d"]], row[["bx3_C4d"]]) == 1, na.rm = TRUE)
class_lbl <- c("1" = "Class I", "2" = "Class II", "3" = "Class I + II")

ab_tbl <- data.frame(
  Case               = seq_len(nrow(dsa_pos)),
  `HLA class`        = unname(class_lbl[as.character(dsa_pos$dnDSA_class)]),
  `Antibodies (MFI)` = sapply(dsa_pos$dnDSA_MFI, relabel_mfi_str),
  Rejection          = ifelse(is.na(dsa_pos$rejection_type) | dsa_pos$rejection_type == 0,
                              "-", as.character(dsa_pos$rejection_type)),
  `C4d positive`     = ifelse(apply(dsa_pos, 1, c4d_any), "Yes", "No"),
  check.names = FALSE, stringsAsFactors = FALSE)
rej_rows <- which(ab_tbl$Rejection != "-")
ft_ab <- flextable(ab_tbl) |> bold(part = "header") |>
  set_table_properties(layout = "autofit", width = 1)
if (length(rej_rows)) ft_ab <- ft_ab |> bold(i = rej_rows) |> color(i = rej_rows, color = col_sig)
doc <- doc |>
  body_add_par(sprintf("dnDSA-positive recipients (n = %d)", nrow(dsa_pos)), style = "heading 1") |>
  body_add_par("Red: biopsy-proven rejection. Class I: HLA-A/B/C; class II: HLA-DR/DQ/DP.",
               style = "Normal") |>
  body_add_flextable(ft_ab) |>
  body_add_par("")

# biopsy-proven rejection
rej_all <- d[!is.na(d$BPAR) & d$BPAR == 1, ]
clean_bx <- function(x) ifelse(is.na(x) | x == 0, "-", as.character(x))
c4d_lbl  <- function(x) ifelse(is.na(x), "-", ifelse(x == 1, "C4d+", "C4d−"))
rej_tbl <- data.frame(
  Case = seq_len(nrow(rej_all)),
  `dnDSA+` = ifelse(rej_all[[OUT]] == 1, "Yes", "No"),
  `Biopsy 1` = clean_bx(rej_all$bx1_findings1), `C4d (Bx1)` = c4d_lbl(rej_all$bx1_C4d),
  `Biopsy 2` = clean_bx(rej_all$bx2_findings1), `C4d (Bx2)` = c4d_lbl(rej_all$bx2_C4d),
  `Biopsy 3` = clean_bx(rej_all$bx3_findings1), `C4d (Bx3)` = c4d_lbl(rej_all$bx3_C4d),
  check.names = FALSE, stringsAsFactors = FALSE)
dsa_rows <- which(rej_tbl$`dnDSA+` == "Yes")
ft_rej <- flextable(rej_tbl) |> bold(part = "header") |>
  set_table_properties(layout = "autofit", width = 1)
if (length(dsa_rows)) ft_rej <- ft_rej |> bold(i = dsa_rows) |> color(i = dsa_rows, color = col_sig)
doc <- doc |>
  body_add_par(sprintf("Biopsy-proven rejection (n = %d)", nrow(rej_tbl)), style = "heading 1") |>
  body_add_par("Red: dnDSA-positive. - = not performed.", style = "Normal") |>
  body_add_flextable(ft_rej)

print(doc, target = file.path(out, "firth_tables.docx"))

cat(sprintf("\nStrongest univariable associations (%d events / %d)\n", n_ev, n_tot))
hl <- res[!is.na(res$p), ]
print(head(hl[order(hl$p), c("label", "level_label", "OR", "lo", "hi", "p", "ev", "N")], 8),
      row.names = FALSE)
