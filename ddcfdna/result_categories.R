# dd-cfDNA results grouped by clinical category
# (quiescent <=0.5%, intermediate 0.5-1.0%, elevated >=1.0%).

in_file <- "data/ddcfdna.csv"
out <- "output/ddcfdna"
dir.create(out, showWarnings = FALSE, recursive = TRUE)

raw <- read.csv(in_file, skip = 1, header = TRUE, stringsAsFactors = FALSE,
                na.strings = c("", "NA", "N/A", "NaN%", "MISSING", "#N/A"))
col <- function(pat) raw[[grep(pat, names(raw), ignore.case = TRUE, value = TRUE)[1]]]

as_pct <- function(x) {
  x <- suppressWarnings(as.numeric(gsub("%", "", x)))
  if (!all(is.na(x)) && max(x, na.rm = TRUE) <= 1) x * 100 else x
}

cf_id <- suppressWarnings(as.numeric(col("cf_DNA_ID")))
res1  <- as_pct(col("timepoint1_cf_DNA_sample1_result"))
res2  <- as_pct(col("timepoint2_cf_DNA_sample1_result"))

# samples excluded at quality control (cf_id, timepoint, analysis)
qc <- read.csv("data/ddcfdna_qc_exclusions.csv", stringsAsFactors = FALSE)
qc <- qc[qc$analysis == "results", ]
excl_tp1 <- qc$cf_id[qc$timepoint == 1]
excl_tp2 <- qc$cf_id[qc$timepoint == 2]

keep1 <- !is.na(cf_id) & !cf_id %in% excl_tp1 & !is.na(res1)
keep2 <- !is.na(cf_id) & !cf_id %in% excl_tp2 & !is.na(res2)
results <- c(res1[keep1], res2[keep2])

grps <- list(
  Quiescent    = results[results <= 0.5],
  Intermediate = results[results > 0.5 & results < 1.0],
  Elevated     = results[results >= 1.0]
)
for (g in names(grps)) {
  x <- grps[[g]]
  cat(sprintf("%-13s n=%d  median %.3f%% (IQR %.3f-%.3f%%)\n", g, length(x),
              median(x), quantile(x, 0.25), quantile(x, 0.75)))
}

cols  <- c("#2ECC71", "#F39C12", "#E74C3C")
light <- c("#A9DFBF", "#FAD7A0", "#F1948A")
labs  <- c("Quiescent\n≤0.5%", "Intermediate\n0.5–1.0%", "Elevated\n≥1.0%")

draw_plot <- function() {
  par(mar = c(5, 5.5, 2, 1.5), family = "sans")
  plot(NA, xlim = c(0.5, 3.8), ylim = c(-0.35, 5.0),
       xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "l")
  axis(2, at = seq(0, 5, 0.5), labels = paste0(seq(0, 5, 0.5), "%"),
       las = 1, cex.axis = 0.95, tcl = -0.3)
  title(ylab = "dd-cfDNA (%)", line = 4, cex.lab = 1.15, font.lab = 2)
  axis(1, at = 1:3, labels = labs, tick = FALSE, cex.axis = 1.0, padj = 0.5)
  abline(h = seq(0, 5, 0.5), col = "grey92", lwd = 0.5)
  abline(h = c(0.5, 1.0), lty = 2, col = "grey60", lwd = 0.9)

  for (i in seq_along(grps)) {
    x   <- grps[[i]]
    q   <- quantile(x, c(0.25, 0.5, 0.75))
    iqr <- q[3] - q[1]
    wlo <- max(min(x), q[1] - 1.5 * iqr)
    whi <- min(max(x), q[3] + 1.5 * iqr)
    inl <- x[x >= wlo & x <= whi]
    outl <- x[x < wlo | x > whi]

    rect(i - 0.22, q[1], i + 0.22, q[3], col = light[i], border = cols[i], lwd = 1.8)
    segments(i - 0.22, q[2], i + 0.22, q[2], col = cols[i], lwd = 2.8)
    segments(i, q[1], i, wlo, col = cols[i], lwd = 1.5)
    segments(i, q[3], i, whi, col = cols[i], lwd = 1.5)
    segments(i - 0.1, c(wlo, whi), i + 0.1, c(wlo, whi), col = cols[i], lwd = 1.5)

    set.seed(i * 42)
    points(i + runif(length(inl), -0.12, 0.12), inl,
           pch = 16, col = adjustcolor(cols[i], 0.7), cex = 0.75)
    if (length(outl)) {
      set.seed(i * 99)
      points(i + runif(length(outl), -0.08, 0.08), outl, pch = 18, col = cols[i], cex = 1.1)
    }
    text(i + 0.27, q[2], paste0(round(q[2], 2), "%"), cex = 0.78, col = cols[i], font = 2, adj = 0)
    text(i, -0.22, paste0("n=", length(x)), cex = 0.85, col = "#444444")
  }

  legend("topleft", legend = gsub("\n", " ", labs), fill = light, border = cols,
         bty = "n", cex = 0.88)
}

pdf(file.path(out, "ddcfdna_categories.pdf"), width = 8, height = 7, bg = "white")
draw_plot()
invisible(dev.off())

png(file.path(out, "ddcfdna_categories.png"), width = 8, height = 7,
    units = "in", res = 300, bg = "white")
draw_plot()
invisible(dev.off())
