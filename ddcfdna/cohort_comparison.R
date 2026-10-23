# dd-cfDNA study: baseline comparison of the original and extended cohorts.

in_file <- "data/ddcfdna.csv"

raw <- read.csv(in_file, skip = 1, header = TRUE, stringsAsFactors = FALSE,
                na.strings = c("", "NA", "N/A", "NaN%", "MISSING", "#N/A"))

col <- function(pat) raw[[grep(pat, names(raw), ignore.case = TRUE, value = TRUE)[1]]]
num <- function(pat) suppressWarnings(as.numeric(col(pat)))

cohort <- col("^cohort$")
c1 <- cohort == "original"
c2 <- cohort == "extended"
cat("Cohort 1:", sum(c1, na.rm = TRUE), " Cohort 2:", sum(c2, na.rm = TRUE), "\n\n")

# CIT < 200 min corresponds to living donation and is excluded
cit <- num("^CIT")
cit[!is.na(cit) & cit <= 200] <- NA

fmtp <- function(p) {
  if (is.na(p)) return("-")
  if (p < 0.001) return("p < 0.001")
  if (p < 0.01) return(sprintf("p = %.3f", p))
  sprintf("p = %.2f", p)
}

mwu <- function(x, label) {
  a <- x[c1 & !is.na(x)]
  b <- x[c2 & !is.na(x)]
  p <- wilcox.test(a, b, exact = FALSE)$p.value
  cat(sprintf("%-45s C1: median %.1f (IQR %.1f-%.1f) | C2: median %.1f (IQR %.1f-%.1f) | %s [Mann-Whitney U]\n",
              label, median(a), quantile(a, 0.25), quantile(a, 0.75),
              median(b), quantile(b, 0.25), quantile(b, 0.75), fmtp(p)))
}

fsh <- function(x_yes, label) {
  a_y <- sum(x_yes[c1], na.rm = TRUE); a_n <- sum(c1, na.rm = TRUE) - a_y
  b_y <- sum(x_yes[c2], na.rm = TRUE); b_n <- sum(c2, na.rm = TRUE) - b_y
  p <- fisher.test(matrix(c(a_y, a_n, b_y, b_n), nrow = 2))$p.value
  cat(sprintf("%-45s C1: %d/%d (%.0f%%) | C2: %d/%d (%.0f%%) | %s [Fisher exact]\n",
              label, a_y, a_y + a_n, 100 * a_y / (a_y + a_n),
              b_y, b_y + b_n, 100 * b_y / (b_y + b_n), fmtp(p)))
}

mwu(num("age_at_tx"), "Recipient age (years)")
mwu(cit, "Cold ischaemia time (min, deceased donors)")
cat("\n")

no_tx  <- num("no_tx")
hla_mm <- num("HLA_MM_level")
crf    <- num("cRF")
dgf    <- num("DGF")

fsh(col("rec_sex") == "M", "Male sex")
fsh(grepl("white", col("ethnicity"), ignore.case = TRUE), "White ethnicity")
fsh(!is.na(no_tx) & no_tx >= 2, "Retransplant")
fsh(!is.na(hla_mm) & hla_mm >= 4, "Level 4 HLA mismatch")
fsh(!is.na(crf) & crf > 20, "cRF >20%")
fsh(grepl("alemtuzumab", col("induction_IS"), ignore.case = TRUE), "Alemtuzumab induction")
fsh(!is.na(dgf) & dgf == 1, "Delayed graft function")
fsh(col("dialysis_status") == "HD", "Haemodialysis")
fsh(col("tx_type") == "DCD", "DCD donor")
