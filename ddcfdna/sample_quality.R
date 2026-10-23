# dd-cfDNA pre-analytical quality: purity and concentration by collection tube
# (EDTA vs Cell3 Preserver) and timepoint.

in_file <- "data/ddcfdna.csv"
out <- "output/ddcfdna"
dir.create(out, showWarnings = FALSE, recursive = TRUE)

raw <- read.csv(in_file, skip = 1, header = TRUE, stringsAsFactors = FALSE,
                na.strings = c("", "NA", "N/A", "NaN%", "MISSING", "#N/A"))
col <- function(pat) raw[[grep(pat, names(raw), ignore.case = TRUE, value = TRUE)[1]]]

as_pct <- function(x) {
  n <- suppressWarnings(as.numeric(trimws(gsub("%", "", as.character(x)))))
  ifelse(!is.na(n) & n <= 1, n * 100, n)
}

df_raw <- data.frame(
  cf_id   = suppressWarnings(as.integer(col("cf_DNA_ID"))),
  cohort  = col("^cohort$"),
  tube1   = col("timepoint1_cf_DNA_collection_blood_tube"),
  conc1   = suppressWarnings(as.numeric(col("cf_DNA_sample1_conc"))),
  purity1 = as_pct(col("timepoint1_cf_DNA_purity")),
  tube2   = col("timepoint2_cf_DNA_collection_blood_tube"),
  conc2   = suppressWarnings(as.numeric(col("timepoint2_cf_DNA_conc"))),
  purity2 = as_pct(col("timepoint2_cf_DNA_purity")),
  stringsAsFactors = FALSE
)
df_raw <- df_raw[!is.na(df_raw$cf_id), ]

is_cell3 <- function(x) grepl("cell3|preserver", x, ignore.case = TRUE)
is_edta  <- function(x) grepl("edta", x, ignore.case = TRUE)

# ---- Tube comparison, both timepoints pooled --------------------------------
pool <- function(v1, v2, f) {
  x <- c(df_raw[[v1]][f(df_raw$tube1)], df_raw[[v2]][f(df_raw$tube2)])
  x[!is.na(x)]
}
edta_pur   <- pool("purity1", "purity2", is_edta)
cell3_pur  <- pool("purity1", "purity2", is_cell3)
edta_conc  <- pool("conc1", "conc2", is_edta)
cell3_conc <- pool("conc1", "conc2", is_cell3)

summ <- function(x, lab, digits) {
  q <- quantile(x, c(0.25, 0.5, 0.75))
  cat(sprintf(paste0("%-6s median %.", digits, "f (IQR %.", digits, "f-%.", digits, "f), n=%d\n"),
              lab, q[2], q[1], q[3], length(x)))
}
cat("Purity (%)\n");  summ(edta_pur, "EDTA", 0);  summ(cell3_pur, "Cell3", 0)
cat("Concentration (ng/uL)\n"); summ(edta_conc, "EDTA", 3); summ(cell3_conc, "Cell3", 3)
print(wilcox.test(edta_pur, cell3_pur, exact = FALSE))
print(wilcox.test(edta_conc, cell3_conc, exact = FALSE))

# ---- Figure by tube and timepoint -------------------------------------------
# samples excluded at quality control are shown but not summarised
qc <- read.csv("data/ddcfdna_qc_exclusions.csv", stringsAsFactors = FALSE)
qc <- qc[qc$analysis == "quality", ]
excl_tp1 <- qc$cf_id[qc$timepoint == 1]
excl_tp2 <- qc$cf_id[qc$timepoint == 2]
in_study <- df_raw$cohort %in% c("original", "extended")

tp1 <- df_raw[in_study & df_raw$cf_id != 28, ]
tp1 <- data.frame(method = ifelse(is_cell3(tp1$tube1), "Cell3 TP1",
                                  ifelse(is_edta(tp1$tube1), "EDTA TP1", NA)),
                  purity = tp1$purity1, conc = tp1$conc1,
                  excluded = tp1$cf_id %in% excl_tp1)

tp2 <- df_raw[in_study & !is.na(df_raw$tube2) & !is.na(df_raw$conc2), ]
tp2 <- data.frame(method = ifelse(is_cell3(tp2$tube2), "Cell3 TP2",
                                  ifelse(is_edta(tp2$tube2), "EDTA TP2", NA)),
                  purity = tp2$purity2, conc = tp2$conc2,
                  excluded = tp2$cf_id %in% excl_tp2)

df <- rbind(tp1, tp2)
df <- df[!is.na(df$method), ]

grp_levels <- c("EDTA TP1", "EDTA TP2", "Cell3 TP1", "Cell3 TP2")
grp_labels <- c("EDTA\nTP1", "EDTA\nTP2", "Cell3\nTP1", "Cell3\nTP2")
pal   <- c("#2166AC", "#6BAED6", "#1A7834", "#74C476")
x_pos <- c(1, 2, 3.5, 4.5)
jit   <- function(n, w = 0.12) runif(n, -w, w)

draw_panel <- function(metric, ylab, log_scale, threshold, show_legend) {
  if (log_scale) {
    ylim    <- c(-2.3, 1.45)
    yticks  <- c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5)
    ylabels <- c("0.01", "0.03", "0.1", "0.32", "1", "3.2", "10", "32")
  } else {
    ylim    <- c(0, 116)
    yticks  <- seq(0, 100, 20)
    ylabels <- paste0(yticks, "%")
  }

  plot(NA, xlim = c(0.3, 5.2), ylim = ylim, xaxt = "n", yaxt = "n",
       xlab = "", ylab = "", bty = "l", las = 1)
  abline(h = yticks, col = "grey92", lwd = 0.6)
  if (!is.null(threshold)) {
    abline(h = if (log_scale) log10(threshold) else threshold, lty = 2, col = "#CC0000", lwd = 1.2)
  }
  axis(2, at = yticks, labels = ylabels, las = 1, cex.axis = 0.9, tcl = -0.3)
  axis(1, at = x_pos, labels = grp_labels, tick = FALSE, cex.axis = 0.95, padj = 0.5)
  title(ylab = ylab, line = 3.2, cex.lab = 1.05, font.lab = 2)

  for (i in seq_along(grp_levels)) {
    sub <- df[df$method == grp_levels[i] & !is.na(df[[metric]]), ]
    iv <- sub[[metric]][!sub$excluded]
    ev <- sub[[metric]][sub$excluded]
    if (log_scale) { iv <- log10(iv[iv > 0]); ev <- log10(ev[ev > 0]) }
    xp <- x_pos[i]

    if (length(iv) > 1) {
      q   <- quantile(iv, c(0.25, 0.5, 0.75))
      iqr <- q[3] - q[1]
      wlo <- max(min(iv), q[1] - 1.5 * iqr)
      whi <- min(max(iv), q[3] + 1.5 * iqr)
      rect(xp - 0.19, q[1], xp + 0.19, q[3],
           col = adjustcolor(pal[i], alpha.f = 0.28), border = pal[i], lwd = 1.2)
      segments(xp - 0.19, q[2], xp + 0.19, q[2], col = "black", lwd = 2.2)
      segments(xp, q[1], xp, wlo, col = pal[i], lwd = 1.2)
      segments(xp, q[3], xp, whi, col = pal[i], lwd = 1.2)
      segments(xp - 0.09, c(wlo, whi), xp + 0.09, c(wlo, whi), col = pal[i], lwd = 1.2)
    }
    if (length(iv)) {
      set.seed(i * 100)
      points(xp + jit(length(iv)), iv, pch = 16, col = adjustcolor(pal[i], alpha.f = 0.88), cex = 0.75)
    }
    if (length(ev)) {
      set.seed(i * 200)
      points(xp + jit(length(ev), w = 0.07), ev, pch = 2, col = "#CC0000", cex = 1.0, lwd = 1.8)
    }
    lbl <- if (length(ev) == 0) paste0("n=", length(iv)) else
      paste0("n=", length(iv), "\n(", length(ev), " excl.)")
    text(xp, ylim[2] * 0.975, lbl, cex = 0.68, col = "#444444", adj = c(0.5, 1))
  }

  if (show_legend) {
    y1 <- if (log_scale) -1.87 else 10
    y2 <- if (log_scale) -2.05 else 5
    points(x_pos[1] - 0.55, y1, pch = 16, col = "#555555", cex = 0.85)
    text(x_pos[1] - 0.38, y1, "Included", cex = 0.78, col = "#555555", adj = 0)
    points(x_pos[1] - 0.55, y2, pch = 2, col = "#CC0000", cex = 0.95, lwd = 1.8)
    text(x_pos[1] - 0.38, y2, "Excluded", cex = 0.78, col = "#CC0000", adj = 0)
  }
  mtext(if (show_legend) "B" else "A", side = 3, line = 0.5, adj = 0, cex = 1.2, font = 2, at = 0.3)
}

draw_figure <- function() {
  par(mfrow = c(1, 2), mar = c(4, 5, 2.5, 1.5), mgp = c(2.8, 0.6, 0), tcl = -0.3, family = "sans")
  draw_panel("purity", "dd-cfDNA purity (%)", FALSE, 80, FALSE)
  draw_panel("conc", "dd-cfDNA concentration (ng/uL, log10 scale)", TRUE, NULL, TRUE)
}

png(file.path(out, "cfdna_quality.png"), width = 3900, height = 1800, res = 300, bg = "white")
draw_figure()
invisible(dev.off())

pdf(file.path(out, "cfdna_quality.pdf"), width = 13, height = 6, bg = "white")
draw_figure()
invisible(dev.off())
