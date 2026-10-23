# Shared plotting theme, palettes and export helpers.
# Sourced by the analysis scripts; run all scripts from the repository root.

library(ggplot2)

theme_thesis <- function(base_size = 11) {
  theme_minimal(base_size = base_size, base_family = "Arial") %+replace%
    theme(
      panel.background  = element_rect(fill = "#FFFFFF", colour = NA),
      plot.background   = element_rect(fill = "#FFFFFF", colour = NA),
      panel.grid.major  = element_line(colour = "#ECECEC", linewidth = 0.4),
      panel.grid.minor  = element_blank(),
      panel.border      = element_blank(),
      axis.line         = element_line(colour = "#7A7A7A", linewidth = 0.4),
      axis.ticks        = element_line(colour = "#7A7A7A", linewidth = 0.4),
      axis.ticks.length = unit(3, "pt"),
      axis.text         = element_text(colour = "#444444", size = base_size - 1),
      axis.title        = element_text(colour = "#222222", size = base_size),
      axis.title.x      = element_text(margin = margin(t = 7)),
      axis.title.y      = element_text(margin = margin(r = 7), angle = 90),
      plot.title        = element_text(colour = "#111111", size = base_size + 2,
                                       face = "bold", hjust = 0,
                                       margin = margin(b = 2)),
      plot.subtitle     = element_text(colour = "#6E6E6E", size = base_size - 1,
                                       hjust = 0, margin = margin(b = 10)),
      plot.caption      = element_text(colour = "#9A9A9A", size = base_size - 3,
                                       hjust = 0, margin = margin(t = 10)),
      plot.title.position   = "plot",
      plot.caption.position = "plot",
      strip.background  = element_blank(),
      strip.text        = element_text(colour = "#222222", size = base_size - 1,
                                       face = "bold", margin = margin(b = 4)),
      legend.background = element_blank(),
      legend.key        = element_blank(),
      legend.text       = element_text(colour = "#444444", size = base_size - 1),
      legend.title      = element_text(colour = "#222222", size = base_size - 1,
                                       face = "bold"),
      legend.position   = "none",
      plot.margin       = margin(14, 16, 12, 14)
    )
}

# Okabe-Ito colours for nominal groups, sequential blues for ordinal quartiles
col_sig    <- "#D55E00"
col_ns     <- "#0072B2"
col_ref    <- "#B4B4B4"
col_grey   <- "#777777"

pal_binary  <- c("#009E73", "#D55E00")
pal_compare <- c("#0072B2", "#D55E00")
pal_age     <- c("#0072B2", "#009E73", "#E69F00", "#CC79A7")
pal_txtype  <- c("#D55E00", "#E69F00", "#009E73")
pal_ukkrri  <- c("#9ECAE1", "#4292C6", "#2171B5", "#084594")
pal_ukkdr   <- c("#C6DBEF", "#9ECAE1", "#4292C6", "#2171B5", "#084594")
pal_whops   <- c("#009E73", "#E69F00", "#D55E00", "#CC79A7")

col_scatter_pt   <- "#0072B2"
col_scatter_line <- "#444444"

dot_col   <- "#9AA0A6"
dot_alpha <- 0.30
dot_size  <- 1.0

# Adds a ggplot to the global slide deck `prs` (created with officer::read_pptx).
# Plots are written as editable vector graphics, with a PNG fallback.
add_plot <- function(p, label = "", w = 12.7, h = 6.9) {
  if (is.null(p) || !inherits(p, "ggplot")) {
    cat(sprintf("  skipped: %s\n", label))
    return(invisible(NULL))
  }
  loc <- officer::ph_location(left = 0.3, top = 0.3, width = w, height = h)
  prs <<- officer::add_slide(prs, layout = "Blank", master = "Office Theme")
  status <- tryCatch({
    prs <<- officer::ph_with(prs, rvg::dml(ggobj = p, bg = "white"), location = loc)
    "vector"
  }, error = function(e1) {
    tryCatch({
      tmp <- tempfile(fileext = ".png")
      ggsave(tmp, p, width = w, height = h, units = "in", dpi = 200, bg = "white")
      prs <<- officer::ph_with(prs, officer::external_img(tmp, width = w, height = h),
                               location = loc)
      "png"
    }, error = function(e2) sprintf("failed: %s", conditionMessage(e2)))
  })
  slide_n <<- slide_n + 1
  cat(sprintf("  slide %d: %s [%s]\n", slide_n, label, status))
}

fmt_p <- function(p) {
  ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

install_missing <- function(pkgs) {
  new <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(new)) install.packages(new, repos = "https://cloud.r-project.org")
}
