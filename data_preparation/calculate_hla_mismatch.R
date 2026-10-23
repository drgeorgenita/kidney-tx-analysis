# HLA mismatch level (1-4) from the A:B:DR mismatch string, e.g. "0:1:0".
# Run on its own to add the level to a spreadsheet, or source() to use add_mm_level().

library(readxl)
library(transplantr)

# Adds `MM level` derived from the `MM` column, reporting malformed strings
# and any disagreement with a level already recorded in the data
add_mm_level <- function(df) {
  df$MM <- toupper(trimws(as.character(df$MM)))
  bad <- !is.na(df$MM) & !grepl("^[0-2]:[0-2]:[0-2]$", df$MM)
  if (any(bad)) message("Rows with invalid MM format: ", paste(which(bad), collapse = ", "))

  recorded <- if ("MM level" %in% names(df)) suppressWarnings(as.numeric(df$`MM level`)) else NA
  df$`MM level` <- hla_mm_level_str(mm = df$MM, sep = TRUE)
  n_diff <- sum(recorded != df$`MM level`, na.rm = TRUE)
  if (n_diff > 0) message(n_diff, " rows where the recorded MM level differs from the derived level")
  df
}

if (sys.nframe() == 0) {
  in_file  <- "data/eligibility.xlsx"
  out_file <- "output/data_preparation/hla_mismatch_level.csv"
  dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)

  df <- add_mm_level(read_excel(in_file))
  print(table(df$`MM level`, useNA = "ifany"))
  write.csv(df, out_file, row.names = FALSE)
}
