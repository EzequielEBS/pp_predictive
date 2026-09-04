library(ggplot2)

# Shared plotting helpers for the power-prior comparison figures.
#
# This file exists to stop the same ~15-line theme() block and the same
# color palettes from being copy-pasted into every panel of
# code/linear_regression/ppc_lm_data.R, code/normal/ppc_normal_data.R,
# code/logistic_regression/actg_plot_post.r and
# code/poisson_regression/e16_plot_post.r. Source it once per script:
#
#   source("code/common/plot_theme.R")

# theme_pp(): the shared look for the conjugate-model prior/posterior/
# posterior-predictive panels (CRPS vs Hyvarinen vs full/no-borrowing vs
# NPP). Every one of those panels used the same base_size = 12 sizing and
# only varied legend.position, so that's the one thing left as an argument.
#
# Note this also fixes a latent ggplot2 deprecation: the legend background
# border used the old `size=` argument to element_rect(), which newer
# ggplot2 versions warn about (and may eventually drop) in favor of
# `linewidth=`.
theme_pp <- function(legend.position = "none", base_size = 12) {
  theme(
    text = element_text(size = base_size),
    axis.title = element_text(size = base_size + 2),
    axis.text = element_text(size = base_size),
    legend.title = element_text(size = base_size + 2),
    legend.text = element_text(size = base_size),
    strip.text = element_text(size = base_size - 1),
    legend.position = legend.position,
    legend.background = element_rect(
      fill = "white",
      color = "black",
      linewidth = 0.3,
      linetype = "solid"
    ),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )
}

# Color palette shared by the conjugate-model comparison plots (lm and
# normal). Each script still adds its own extra keys (e.g. "beta_curr" /
# "beta_hist" or "y0" / "y") on top of this via
# c(pal_pp_conjugate, "beta_curr" = "black", ...).
pal_pp_conjugate <- c(
  "CRPS"    = "#66A8D0",
  "Hyva"    = "#D06673",
  "eta = 0" = "#7f7f7f",
  "eta = 1" = "#7f7f7f",
  "NPP"     = "#D0C366"
)

# Color palette + linetypes shared by the GLM eta-comparison plots
# (actg_plot_post.r and e16_plot_post.r).
pal_pp_glm <- c(
  "big"       = "#D62828",  # red
  "inf_match" = "#F77F00",  # orange
  "inter"     = "#00B4D8",  # sky blue
  "rate"      = "#7B2D8B",  # purple
  "npp"       = "#023E8A",  # dark blue
  "small"     = "#2D6A4F",  # forest green
  "vnpp"      = "#E9C46A"   # gold
)

lty_pp_glm <- c(
  "big"       = "solid",
  "inf_match" = "solid",
  "inter"     = "solid",
  "rate"      = "solid",
  "npp"       = "dashed",
  "small"     = "solid",
  "vnpp"      = "dashed"
)
