## code/linear_regression/pp_npp_results_table.R
##
## Builds THREE gt-styled panels -- Average MSE, Average WIS, Average Bias --
## comparing PP vs. NPP against the KNOWN true beta/sigma used to generate
## the simulated data, from pp_npp_posterior_by_replicate.RData's
## quantile_df (produced by pp_npp_posterior_by_replicate.R), and combines
## them into ONE side-by-side HTML page (style modeled on a reference image
## the user shared: bold grey row-group headers, a "Method" spanner, and
## each panel colored with its own hue -- see the coloring section below).
## Row groups are "Regime" (congruence level: High/Small/No), with one row
## per beta coefficient / sigma inside each group -- NOT the (theta0,
## Discrepancy) row layout from the reference image, which comes from a
## different part of the project with a different scenario design; this
## data only has congruence + parameter to group/label rows by.
##
## GROUND TRUTH: true_params.RData (in the sibling onpp/simulated-data-
## regression/data/ folder) stores `beta` (a list, one entry per p in
## c(3, 10, 50, 100) -- see generate_data.r) and `sigma` (a single scalar,
## shared by every scenario). For our p = 3 scenarios (formula =
## y ~ X1 + X2), the true parameter vector is beta[[1]] = c(intercept, X1,
## X2). Congruence level (High/Small/No) only changes how far the
## HISTORICAL data's beta is scaled away from this true beta (see
## aux_fun.r::dkl_beta and generate_data.r's run_beta_block) -- it does NOT
## change the current/true data-generating parameters. So the SAME truth
## vector applies to every row of quantile_df regardless of congruence.
##
## ETA IS EXCLUDED from this table: it has no "true value" to score against
## -- it's an information-matching / discounting artifact, not a
## data-generating parameter. (Under PP its saved "posterior" is a
## degenerate point mass at eta_hat anyway -- see
## pp_npp_posterior_by_replicate.R's header.) Bias/MSE/WIS here only cover
## the beta coefficients and sigma, which do have known truth.
##
## POINT ESTIMATE for Bias/MSE: the posterior MEDIAN (quantile_level ==
## 0.5) -- the natural point forecast paired with WIS (Bracher et al.
## 2021), and the only point summary quantile_df guarantees is present
## (mean/sd were dropped when the checkpoint switched to full quantiles).
##
## WIS: computed per (congruence, replicate, method, parameter) directly
## from that row's saved quantile curve, via the standard interval-score
## decomposition (Bracher, Held, Lauer & Reich, 2021 -- "Evaluating
## epidemic forecasts in an interval format", eq. 4), generalized to
## whatever QUANTILE_LEVELS the posterior run used. Requires that grid to
## be symmetric around 0.5 and include 0.5 -- true for the default grid in
## pp_npp_posterior_by_replicate.R, and checked below.
##
## OUTPUT: one self-contained HTML page with three gt tables side by side
## (MSE, WIS, Bias), each with one row per (congruence, parameter) and one
## column per method (PP, NPP). Saved as figures/lm_results_table.html
## (always works) and, if webshot2 is installed, a combined
## figures/lm_results_table.png screenshot of that page.
##
## COLORING: each panel gets its OWN palette (MSE: blue, WIS: warm
## purple-orange, Bias: diverging blue-white-red, since Bias -- unlike
## MSE/WIS -- can be negative and the sign matters). Within a panel, each
## ROW is colored relative to its own PP/NPP values (not a single scale
## across the whole table) -- same convention as the row-wise data_color()
## loop in the user's own build_summary_gt_table(). For Bias specifically,
## the per-row domain is symmetric around 0 (+-max(|PP|,|NPP|)) rather than
## [min, max], so a cell's color reflects both its sign and its magnitude,
## not just its rank against the other method in that row.
##
## NOTE ON TESTING: the WIS formula itself (wis_score() below) was
## independently verified -- against a closed-form degenerate-distribution
## case and against an independent pinball-loss-average formula that is
## mathematically equivalent to the interval-score formula -- see the
## conversation this script was drafted in. The gt table rendering was NOT
## tested (gt/gtExtras aren't installed in the sandbox this was drafted
## in) -- double check the table's visual output once run for real.

library(dplyr)
library(tidyr)
library(gt)
# gtExtras/ggsci are NOT needed here -- every cell is colored explicitly via
# data_color() per row (see build_pp_npp_gt_table() below), so there's no
# separate whole-column base layer to also apply.

## ---------------------------------------------------------------------------
## Configuration
## ---------------------------------------------------------------------------

data_dir         <- Sys.getenv("ONPP_DATA_DIR", "../../onpp/simulated-data-regression/data")
RESULTS_PATH     <- "samples/linear_regression/pp_npp_posterior_by_replicate.RData"
TRUE_PARAMS_PATH <- file.path(data_dir, "true_params.RData")
FIG_DIR          <- "figures"

CONGRUENCE_ORDER <- c("High congruence", "Small congruence", "No congruence")
PARAMETER_ORDER  <- c("(Intercept)", "X1", "X2", "sigma")

if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

## ---------------------------------------------------------------------------
## Load results + ground truth
## ---------------------------------------------------------------------------

stopifnot(
  "Run pp_npp_posterior_by_replicate.R first -- checkpoint not found" =
    file.exists(RESULTS_PATH),
  "true_params.RData not found -- check ONPP_DATA_DIR" =
    file.exists(TRUE_PARAMS_PATH)
)

results         <- readRDS(RESULTS_PATH)
quantile_df     <- results$quantile_df
QUANTILE_LEVELS <- results$quantile_levels

## ---------------------------------------------------------------------------
## Defensive dedup: a checkpoint can end up with duplicate rows for a given
## (congruence, replicate, method, parameter, quantile_level) -- e.g. from
## re-running a code chunk interactively in an IDE (RStudio/Positron) before
## letting the full loop run, which appends to the in-memory quantile_df each
## time it's re-executed, and then gets persisted by the next save_checkpoint().
## wis_score() requires exactly one row per quantile level, so duplicates here
## would otherwise fail loudly and unhelpfully deep inside the WIS computation.
## Exact duplicates (same value every time) are safe to collapse automatically.
## TRUE conflicts (duplicate rows with DIFFERENT values) are not safe to guess
## at, so those stop the script with the offending combos printed instead.
## ---------------------------------------------------------------------------

dup_check <- quantile_df %>%
  group_by(congruence, replicate, method, parameter, quantile_level) %>%
  summarise(n = n(), n_distinct_values = n_distinct(value), .groups = "drop") %>%
  filter(n > 1)

if (nrow(dup_check) > 0) {
  conflicting <- dup_check %>% filter(n_distinct_values > 1)
  if (nrow(conflicting) > 0) {
    print(conflicting)
    stop(
      nrow(conflicting), " (congruence, replicate, method, parameter, quantile_level) ",
      "combination(s) above have duplicate rows with DIFFERENT values -- this can't be ",
      "safely deduplicated automatically. Likely cause: the checkpoint accumulated results ",
      "from more than one attempt (e.g. the loop was interrupted and restarted without a ",
      "clean resume). Inspect and fix samples/linear_regression/",
      "pp_npp_posterior_by_replicate.RData directly (clear the affected combo's rows and ",
      "re-run pp_npp_posterior_by_replicate.R for it) before re-running this script."
    )
  }
  message(
    nrow(dup_check), " (congruence, replicate, method, parameter, quantile_level) ",
    "combination(s) had duplicate but IDENTICAL rows (likely a re-run code chunk before ",
    "the main loop started) -- deduplicating automatically."
  )
  quantile_df <- distinct(quantile_df)
}

stopifnot(
  "quantile_df has no rows -- nothing to summarize" = nrow(quantile_df) > 0,
  "quantile_level = 0.5 (median) is required for the Bias/MSE point estimate" =
    0.5 %in% QUANTILE_LEVELS,
  "QUANTILE_LEVELS must be symmetric around 0.5 for the WIS formula below" =
    isTRUE(all.equal(sort(QUANTILE_LEVELS), sort(1 - QUANTILE_LEVELS)))
)

load(TRUE_PARAMS_PATH)  # loads `beta` (list) and `sigma` (scalar)
# p = 3 (intercept + X1 + X2) is beta[[1]] in generate_data.r's ps <- c(3, 10, 50, 100).
stopifnot(
  "Expected beta[[1]] to have length 3 (intercept, X1, X2) -- check generate_data.r's `ps`" =
    length(beta[[1]]) == 3
)
TRUTH <- c(
  "(Intercept)" = beta[[1]][1], "X1" = beta[[1]][2], "X2" = beta[[1]][3],
  "sigma" = sigma
)

## ---------------------------------------------------------------------------
## WIS (Weighted Interval Score) -- Bracher, Held, Lauer & Reich (2021)
## ---------------------------------------------------------------------------

#' WIS for one predictive quantile curve against one true value.
#'
#' @param levels quantile levels (any order; must include 0.5 and be
#'   symmetric around it)
#' @param values quantile VALUES at `levels`, same order/length
#' @param y the true value being scored against
wis_score <- function(levels, values, y) {
  values <- values[order(levels)]
  levels <- sort(levels)
  m <- values[levels == 0.5]
  stopifnot("wis_score(): levels must include exactly one 0.5 entry" = length(m) == 1)

  lower_levels <- levels[levels < 0.5]
  K <- length(lower_levels)
  if (K == 0) return(abs(y - m))  # only the median was supplied -- degenerate to absolute error

  interval_score <- vapply(lower_levels, function(tau) {
    alpha <- 2 * tau
    l <- values[levels == tau]
    u <- values[abs(levels - (1 - tau)) < 1e-8]
    stopifnot("wis_score(): quantile levels must be symmetric around 0.5 (missing upper pair)" =
                length(u) == 1)
    is_alpha <- (u - l) +
      (2 / alpha) * (l - y) * (y < l) +
      (2 / alpha) * (y - u) * (y > u)
    (alpha / 2) * is_alpha
  }, numeric(1))

  (0.5 * abs(y - m) + sum(interval_score)) / (K + 0.5)
}

## ---------------------------------------------------------------------------
## Per-replicate scoring, then aggregate to (congruence, parameter, method)
## ---------------------------------------------------------------------------

per_replicate <- quantile_df %>%
  filter(parameter %in% names(TRUTH)) %>%  # drop eta -- no ground truth to score against
  group_by(congruence, replicate, method, parameter) %>%
  summarise(
    true_value = TRUTH[[unique(parameter)]],
    median_est = value[quantile_level == 0.5],
    wis        = wis_score(quantile_level, value, TRUTH[[unique(parameter)]]),
    .groups    = "drop"
  ) %>%
  mutate(
    error    = median_est - true_value,
    sq_error = error^2
  )

summary_df <- per_replicate %>%
  group_by(congruence, parameter, method) %>%
  summarise(
    n    = n(),
    Bias = abs(mean(error)),
    MSE  = mean(sq_error),
    WIS  = mean(wis),
    .groups = "drop"
  )

## ---------------------------------------------------------------------------
## Reshape each metric into its own (Regime x Parameter) x (PP, NPP) table.
## ---------------------------------------------------------------------------

# Display labels for each parameter, as the beta/sigma symbols used in the
# thesis write-up rather than the raw R coefficient/column names. Plain
# Unicode (Greek letter + subscript digit), written as \u escapes rather
# than literal characters in the source file so it's not at the mercy of
# the file's on-disk encoding matching R's session encoding (a common
# mojibake source on Windows) -- renders correctly in any browser with no
# extra markup, JS, or MathJax dependency, and stays correct if the table
# is ever copied into plain text.
PARAMETER_LABELS <- c(
  "(Intercept)" = "\u03b2\u2080",  # beta_0
  "X1"          = "\u03b2\u2081",  # beta_1
  "X2"          = "\u03b2\u2082",  # beta_2
  "sigma"       = "\u03c3"          # sigma
)

# Display labels for congruence level, framed as discrepancy (how far the
# historical data's beta is from the true/current beta -- see the GROUND
# TRUTH note above) rather than the internal "congruence" naming.
CONGRUENCE_LABELS <- c(
  "High congruence"  = "No discrepancy",
  "Small congruence" = "Small discrepancy",
  "No congruence"    = "Large discrepancy"
)

make_metric_table <- function(summary_df, metric) {
  summary_df %>%
    select(congruence, parameter, method, value = all_of(metric)) %>%
    pivot_wider(names_from = method, values_from = value) %>%
    mutate(
      congruence = factor(congruence, levels = CONGRUENCE_ORDER),
      parameter  = factor(parameter, levels = PARAMETER_ORDER)
    ) %>%
    arrange(congruence, parameter) %>%
    mutate(
      congruence = CONGRUENCE_LABELS[as.character(congruence)],  # relabel AFTER sort order is locked in
      parameter  = PARAMETER_LABELS[as.character(parameter)]
    ) %>%
    rename(Regime = congruence, Parameter = parameter)
}

mse_df  <- make_metric_table(summary_df, "MSE")
wis_df  <- make_metric_table(summary_df, "WIS")
bias_df <- make_metric_table(summary_df, "Bias")

## ---------------------------------------------------------------------------
## gt table builder -- same styling as build_summary_gt_table() (bold grey
## row-group headers, "Method" spanner), adapted to our 2 methods (PP, NPP)
## and a Parameter column instead of Eps/Discrepancy. Each panel supplies
## its own row_palette (color hue) and row_domain (how a row's color scale
## is anchored) so MSE/WIS/Bias can each get a distinct look.
## ---------------------------------------------------------------------------

# Sequential row_domain: color anchored at [min, max] of that row's PP/NPP
# values -- appropriate for MSE and WIS, which are always >= 0.
sequential_domain <- function(row_vals) c(min(row_vals), max(row_vals))

# Diverging row_domain: color anchored at +-max(|PP|, |NPP|), i.e. symmetric
# around 0 -- appropriate for Bias, which can be positive or negative and
# whose SIGN matters, not just its rank against the other method.
diverging_domain <- function(row_vals) {
  m <- max(abs(row_vals))
  if (m == 0) c(-1, 1) else c(-m, m)  # avoid a zero-width domain when both are exactly 0
}

blue_palette <- function(n, direction = 1) {
  cols <- grDevices::colorRampPalette(c("#eaf1fb", "#1f3a7a"))(n)  # light -> deep blue
  if (direction == -1) cols <- rev(cols)
  cols
}
warm_palette <- function(n, direction = 1) {
  cols <- grDevices::colorRampPalette(c("#f7b267", "#c9184a", "#6a1b6a"))(n)  # orange -> magenta -> purple
  if (direction == -1) cols <- rev(cols)
  cols
}
diverging_palette <- function(n, direction = 1) {
  cols <- grDevices::colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(n)  # blue -> white -> red
  if (direction == -1) cols <- rev(cols)
  cols
}

# show_row_labels = FALSE is for panels placed to the right of the FIRST
# (fully-labeled) panel in the side-by-side layout: the Parameter column and
# the Regime/congruence row-group TEXT would just repeat what the first
# panel already shows, so this drops the Parameter column entirely and
# blanks the row-group label text. The row-GROUPING itself (and therefore
# the grey band + row heights) is deliberately kept even when the text is
# blanked, so every panel still has the exact same number/position of rows
# and bands -- that's what keeps a row in this panel lined up with the same
# row in the first panel once they're side by side.
build_pp_npp_gt_table <- function(df, title, row_palette, row_domain = sequential_domain,
                                   show_row_labels = TRUE) {
  method_cols <- c("PP", "NPP")

  tbl <- df |>
    gt(groupname_col = "Regime") |>
    cols_label(Parameter = "Parameter", PP = "PP", NPP = "NPP") |>
    fmt_number(
      columns  = all_of(method_cols),
      decimals = 4
    ) |>
    tab_header(title = title) |>
    tab_spanner(
      label   = "Method",
      columns = all_of(method_cols)
    ) |>
    tab_style(
      style     = cell_text(weight = "bold"),
      locations = cells_row_groups()
    ) |>
    tab_options(
      row_group.font.weight       = "bold",
      row_group.background.color  = "#f0f0f0",
      table.font.size             = 11,
      # gt tables default to an auto-centering left/right margin, which
      # looks like uneven/wrong spacing once several tables are placed side
      # by side in a flex container (the flex `gap` below ends up fighting
      # with each table's own margin). Zeroing it here makes the flex gap
      # the ONLY spacing between panels.
      table.margin.left           = 0,
      table.margin.right          = 0
    )

  if (!show_row_labels) {
    tbl <- tbl |>
      cols_hide(columns = "Parameter") |>
      text_transform(
        locations = cells_row_groups(),
        # A truly EMPTY string collapses the row-group band's height in
        # some browsers (no text -> no line-height to size the row by),
        # which is what broke row alignment between panels -- a single
        # non-breaking space keeps the cell "occupied" by one line of text
        # (same line-height as a real label) while remaining invisible.
        fn        = function(x) rep("\u00a0", length(x))
      )
  }

  # row-wise coloring, anchored per row_domain() -- same convention as
  # build_summary_gt_table()'s per-row loop, generalized to allow a
  # diverging (sign-aware) domain for Bias.
  for (i in seq_len(nrow(df))) {
    row_vals <- as.numeric(df[i, method_cols])
    tbl <- tbl |>
      data_color(
        columns = all_of(method_cols),
        rows    = i,
        palette = row_palette(100, direction = -1)[30:70],
        domain  = row_domain(row_vals)
      )
  }
  tbl
}

# Bias is the FIRST panel in the side-by-side layout (see combine_..._html()
# call below), so it's the one that keeps its Parameter column and Regime
# row-group text -- MSE and WIS, now to its right, have theirs hidden
# (show_row_labels = FALSE) since they'd just repeat what Bias already shows.
gt_bias <- build_pp_npp_gt_table(bias_df, "|Bias|", row_palette = diverging_palette,
                                  row_domain = diverging_domain)
gt_mse  <- build_pp_npp_gt_table(mse_df,  "MSE",  row_palette = blue_palette,
                                  show_row_labels = FALSE)
gt_wis  <- build_pp_npp_gt_table(wis_df,  "WIS",  row_palette = warm_palette,
                                  show_row_labels = FALSE)

## ---------------------------------------------------------------------------
## Combine the three panels into one side-by-side HTML page. gt tables don't
## natively lay out side by side within a single gt object, so each is
## rendered to self-contained HTML via gt::as_raw_html() and wrapped in a
## flex container -- the standard approach for combining independent gt
## tables into one page.
## ---------------------------------------------------------------------------

combine_side_by_side_html <- function(tables, path, gap_px = 20) {
  # Each raw-html gt table still carries a <style> block with its own
  # (now-zeroed) margins -- wrap each in a div with NO extra margin/padding
  # of its own, so `gap` on the flex container is the only source of space
  # between panels.
  #
  # `flex: 0 0 auto` (NOT `flex: 1 ...`) is the actual fix for the huge gaps
  # seen between panels: `flex-grow: 1` makes each item STRETCH to fill an
  # equal share of the container's width, and since a narrow gt table
  # doesn't itself stretch, the leftover stretched space inside each flex
  # item shows up as a big empty gap before the next one -- that's what
  # `flex: 1 1 0` was doing. `flex: 0 0 auto` sizes each div to its table's
  # natural width instead, so the actual `gap` value is the only space
  # between panels.
  panels <- vapply(tables, function(t) as.character(gt::as_raw_html(t)), character(1))
  divs <- paste0('<div style="flex: 0 0 auto; margin: 0;">', panels, "</div>")
  page <- paste0(
    '<!doctype html><html><head><meta charset="utf-8">',
    "<title>PP vs. NPP results</title></head>",
    '<body style="font-family: -apple-system, Helvetica, Arial, sans-serif; padding: 16px; margin: 0;">',
    sprintf(
      # `width: fit-content` shrinks the flex container down to the combined
      # natural width of its `flex: 0 0 auto` children instead of stretching
      # to fill its parent -- this is what fixes the layout in the HTML page
      # itself when opened directly in a browser. It does NOT, by itself, fix
      # the PNG export below: webshot2 screenshots the full browser viewport
      # (fixed at vwidth = 1700), so a narrower page still leaves blank space
      # in the image unless the screenshot is cropped to this div specifically
      # (see id="panels-wrapper" + the `selector` argument in the PNG export).
      '<div id="panels-wrapper" style="display:flex; flex-wrap:wrap; gap:%dpx; align-items:flex-start; justify-content:flex-start; width:fit-content;">',
      gap_px
    ),
    paste(divs, collapse = "\n"),
    "</div></body></html>"
  )
  writeLines(page, path)
}

## ---------------------------------------------------------------------------
## Save
## ---------------------------------------------------------------------------

HTML_PATH <- file.path(FIG_DIR, "lm_results_table.html")
PNG_PATH  <- file.path(FIG_DIR, "lm_results_table.png")

combine_side_by_side_html(list(gt_bias, gt_mse, gt_wis), HTML_PATH)

# PNG export needs webshot2 (screenshots the combined HTML page) -- wrapped
# so a missing/broken headless browser doesn't stop the (always-working)
# HTML export above.
#
# `selector = "#panels-wrapper"` crops the screenshot to that div's own
# bounding box instead of the full vwidth x vheight browser viewport. This
# is needed IN ADDITION to the `width: fit-content` CSS fix above: the CSS
# fix shrinks the div itself when the page is opened in a normal browser,
# but webshot2 always renders into a fixed-size 1700x1400 virtual window
# first -- without `selector`, the exported PNG would just be a screenshot
# of that whole window (i.e. the div plus all the blank space around it),
# regardless of how narrow the div's own CSS makes it.
tryCatch({
  if (!requireNamespace("webshot2", quietly = TRUE)) {
    stop("webshot2 not installed")
  }
  webshot2::webshot(
    url      = HTML_PATH,
    file     = PNG_PATH,
    selector = "#panels-wrapper",
    vwidth   = 1700,
    vheight  = 1400
  )
}, error = function(e) {
  message("PNG export skipped/failed (HTML page was still saved): ", conditionMessage(e))
})

cat("Saved combined MSE/WIS/Bias panels to ", HTML_PATH, "\n", sep = "")
print(gt_bias)
print(gt_mse)
print(gt_wis)
