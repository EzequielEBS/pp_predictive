## code/linear_regression/pp_npp_posterior_by_replicate.R
##
## For each (congruence, replicate) already estimated in
## samples/linear_regression/eta_by_replicate.RData, fits the posterior for
## the regression coefficients (beta), residual dispersion (sigma), and the
## power-prior discount (eta) two ways, BOTH via hdbayes (Stan), so the
## comparison uses the same inference engine on both sides:
##   (1) Power prior (PP) via hdbayes::glm.pp() with a0.vals FIXED at the
##       information-matching estimate eta_hat. Because eta is fixed (not
##       sampled), its "posterior" under PP is degenerate -- a point mass at
##       eta_hat -- so every quantile level for eta under PP gets the same
##       value. This is intentional: it's exactly what makes the PP vs. NPP
##       comparison meaningful (NPP treats eta's uncertainty honestly; PP
##       plugs in a point estimate and pretends there's no uncertainty about
##       it at all).
##   (2) Normalized power prior (NPP) via hdbayes::glm.npp() -- eta is a
##       full posterior, not fixed. Matches actg_sample.r's methodology.
##
## QUANTILES (not just mean/sd/95% CI) are saved for every parameter (each
## beta coefficient, sigma, eta), at the levels in QUANTILE_LEVELS below, so
## downstream code can compute the Weighted Interval Score (WIS) against the
## known simulation truth. QUANTILE_LEVELS is a parameter -- edit it to
## whatever grid your WIS computation expects. It defaults to the standard
## 11-interval / 23-level grid used by the COVID Forecast Hub convention
## (Bracher et al. 2021): 0.01, 0.025, 0.05, 0.10, ..., 0.90, 0.95, 0.975,
## 0.99, plus the median (0.5).
##
## "sigma" is hdbayes's "dispersion" draw column for family = gaussian().
## Double-check against the hdbayes docs for your installed version whether
## that's the residual SD or variance before feeding it into WIS -- this
## wasn't verified in the sandbox this script was drafted in.
##
## OUTPUT FORMAT: quantile_df is LONG -- one row per
## (congruence, replicate, method, parameter, quantile_level) -- which is
## the natural shape for a WIS join against ground truth (join on
## congruence/replicate/parameter, and quantile_level against the interval
## bounds WIS needs). "parameter" is the beta coefficient name (e.g.
## "(Intercept)", "X1", "X2"), "sigma", or "eta".
##
## COST WARNING -- this is substantial for BOTH methods, at full scale
## (NUM_SIM = 200 replicates x 3 congruence levels = 600 combos):
##   - PP:  1 Stan fit per combo                      ->    600 Stan fits
##   - NPP: (N_A0_GRID + 1) Stan fits per combo        -> 13,200 Stan fits
##   - Total: ~13,800 Stan fits, each 4 chains x (iter_warmup + iter_sampling)
## draws, ON TOP of the 600 already run for eta estimation. This can run for
## a very long time (plausibly days) depending on your machine. The script
## checkpoints after EVERY combo specifically because of this.
##
## CHECKPOINTING: results saved after every combo; a rerun skips (congruence,
## replicate) pairs already present in eta_df -- including ones where PP or
## NPP failed (pp_error / npp_error columns record which, if either, failed;
## neither is auto-retried -- see eta_estimation_by_replicate.R's header for
## the general pattern of clearing specific rows to force a retry).
## NOTE: this checkpoint's schema (quantile_df, long format) is NOT
## compatible with an older checkpoint saved before quantiles were added
## (which had a wide beta_df of post_mean/post_sd/q025/q975 only). Loading
## an old-format file will stop with an explicit message rather than
## silently losing information -- move/delete it to start fresh.
## Similarly, resuming with a DIFFERENT QUANTILE_LEVELS than the checkpoint
## was started with is refused, since that would silently mix two different
## quantile grids in the same quantile_df.

library(dplyr)
library(hdbayes)
library(posterior)
library(parallel)
library(readr)

set.seed(20260819)

## ---------------------------------------------------------------------------
## Configuration
## ---------------------------------------------------------------------------

data_dir <- Sys.getenv("ONPP_DATA_DIR", "../../onpp/simulated-data-regression/data")

CONGRUENCE_FILES <- c(
  "High congruence"  = "sim_data_high_cong_varying_beta_p3_n0_ge_n.csv",
  "Small congruence" = "sim_data_small_cong_varying_beta_p3_n0_ge_n.csv",
  "No congruence"    = "sim_data_no_cong_varying_beta_p3_n0_ge_n.csv"
)

formula <- y ~ X1 + X2
family  <- gaussian()

NUM_SIM   <- 200  # full scale, both methods -- see COST WARNING above
N_A0_GRID <- 21   # a0 grid for glm.npp.lognc (matches actg_sample.r)

# Quantile levels saved for every parameter (beta coefficients, sigma, eta).
# Edit this to match whatever grid your WIS computation expects -- default
# is the standard COVID-Forecast-Hub-style 11-interval / 23-level grid.
QUANTILE_LEVELS <- c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975)

ncores        <- detectCores() - 1
chains        <- 4
iter_warmup   <- 5000
iter_sampling <- 10000

BASE_SEED <- 23082026

ETA_RESULTS_PATH <- "samples/linear_regression/eta_by_replicate.RData"
RESULTS_PATH     <- "samples/linear_regression/pp_npp_posterior_by_replicate.RData"
if (!dir.exists(dirname(RESULTS_PATH))) dir.create(dirname(RESULTS_PATH), recursive = TRUE)

## ---------------------------------------------------------------------------
## Load eta estimates and raw data
## ---------------------------------------------------------------------------

stopifnot(
  "Run eta_estimation_by_replicate.R first -- eta_by_replicate.RData not found" =
    file.exists(ETA_RESULTS_PATH)
)
eta_results <- readRDS(ETA_RESULTS_PATH) %>%
  filter(is.na(error), replicate <= NUM_SIM)

congruence_data <- lapply(CONGRUENCE_FILES, function(fname) {
  read_csv(file.path(data_dir, fname), show_col_types = FALSE)
})
names(congruence_data) <- names(CONGRUENCE_FILES)

## ---------------------------------------------------------------------------
## Quantile-row helpers
## ---------------------------------------------------------------------------

# One draws vector -> one row per quantile level.
quantile_rows <- function(cong, rep, method, parameter, draws_vec, levels) {
  data.frame(
    congruence = cong, replicate = rep, method = method,
    parameter = parameter, quantile_level = levels,
    value = as.numeric(stats::quantile(draws_vec, probs = levels, names = FALSE, type = 7)),
    stringsAsFactors = FALSE
  )
}

# A fixed (non-sampled) point value -> the same value repeated at every
# quantile level (a degenerate distribution). Used for eta under PP, where
# eta is plugged in rather than sampled.
degenerate_quantile_rows <- function(cong, rep, method, parameter, point_value, levels) {
  data.frame(
    congruence = cong, replicate = rep, method = method,
    parameter = parameter, quantile_level = levels,
    value = point_value, stringsAsFactors = FALSE
  )
}

# A named draws matrix (one column per parameter) -> quantile_rows() for
# every column, stacked.
matrix_to_quantile_df <- function(cong, rep, method, draws_mat, levels) {
  do.call(rbind, lapply(colnames(draws_mat), function(pn) {
    quantile_rows(cong, rep, method, pn, draws_mat[, pn], levels)
  }))
}

## ---------------------------------------------------------------------------
## Checkpoint load / init
## ---------------------------------------------------------------------------

if (file.exists(RESULTS_PATH)) {
  results <- readRDS(RESULTS_PATH)
  if (is.null(results$quantile_df)) {
    stop(
      "Existing checkpoint at ", RESULTS_PATH, " uses the old summary-stats ",
      "schema (a wide beta_df of post_mean/post_sd/q025/q975) from before ",
      "full quantiles were added. It can't be resumed by this version of ",
      "the script -- those summary stats aren't enough to reconstruct full ",
      "quantiles for WIS. Move or delete the old file to start a fresh run."
    )
  }
  if (!isTRUE(all.equal(results$quantile_levels, QUANTILE_LEVELS))) {
    stop(
      "QUANTILE_LEVELS has changed since this checkpoint was started ",
      "(checkpoint: ", paste(results$quantile_levels, collapse = ", "),
      "; current: ", paste(QUANTILE_LEVELS, collapse = ", "), "). ",
      "Resuming would silently mix two different quantile grids in the ",
      "same quantile_df. Restore the original QUANTILE_LEVELS, or move/",
      "delete the checkpoint to start fresh with the new grid."
    )
  }
  quantile_df <- results$quantile_df
  eta_df      <- results$eta_df
  cat(sprintf("Resuming: %d quantile rows, %d combos already done\n",
              nrow(quantile_df), nrow(eta_df)))
} else {
  quantile_df <- data.frame(
    congruence = character(), replicate = integer(), method = character(),
    parameter = character(), quantile_level = double(), value = double(),
    stringsAsFactors = FALSE
  )
  eta_df <- data.frame(
    congruence = character(), replicate = integer(), eta_hat = double(),
    pp_elapsed_sec = double(), npp_elapsed_sec = double(),
    pp_error = character(), npp_error = character(),
    stringsAsFactors = FALSE
  )
}

combo_done <- function(cong, rep) {
  any(eta_df$congruence == cong & eta_df$replicate == rep)
}

# Prints a progress line immediately (flushed), instead of buffering until
# R decides to show it -- important here because a single combo (PP fit +
# the NPP a0 grid + the NPP fit) can run for many minutes with otherwise no
# visible output (refresh = 0 silences Stan's own per-chain progress).
progress <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), sprintf(...), "\n", sep = "")
  flush(stdout())
}

save_checkpoint <- function() {
  saveRDS(
    list(quantile_df = quantile_df, eta_df = eta_df, quantile_levels = QUANTILE_LEVELS),
    RESULTS_PATH
  )
}

## ---------------------------------------------------------------------------
## Main loop: one combo at a time, PP then NPP, both via hdbayes/Stan.
## Sequential outer loop -- glm.pp()/glm.npp() already parallelize their own
## chains (parallel_chains = ncores) and the NPP lognc grid uses its own
## cluster, so parallelizing this outer loop too would oversubscribe.
## ---------------------------------------------------------------------------

n_jobs <- nrow(eta_results)
cat(sprintf("Total combos: %d (%d congruence levels x %d replicates)\n",
            n_jobs, length(CONGRUENCE_FILES), NUM_SIM))

for (j in seq_len(n_jobs)) {
  cong    <- eta_results$congruence[j]
  rep     <- eta_results$replicate[j]
  eta_hat <- eta_results$eta_hat[j]

  if (combo_done(cong, rep)) next

  progress("[%d/%d] %s | replicate %d | starting (eta_hat = %.4f)",
            j, n_jobs, cong, rep, eta_hat)

  data_c <- congruence_data[[cong]]
  hist_data <- data_c %>% filter(data_id == "hist3", replicate == rep)
  curr_data <- data_c %>% filter(data_id == "current", replicate == rep)
  data_list <- list(curr_data, hist_data)
  seed_here <- BASE_SEED + which(names(CONGRUENCE_FILES) == cong) * 1e6 + rep

  ## ---- PP via hdbayes::glm.pp(a0.vals = eta_hat) ----
  progress("  PP: fitting (1 Stan fit, %d chains)...", chains)
  t0 <- proc.time()[["elapsed"]]
  pp_outcome <- tryCatch({
    set.seed(seed_here)  # glm.pp wraps cmdstanr sampling; matches the
                          # set.seed()-before-fit convention already used
                          # elsewhere in this repo (see actg_sample.r's note
                          # on cmdstanr not always honoring R's global RNG)
    fit_pp <- glm.pp(
      formula = formula, family = family, data.list = data_list,
      a0.vals = eta_hat,
      iter_warmup = iter_warmup, iter_sampling = iter_sampling,
      chains = chains, parallel_chains = ncores, refresh = 0
    )
    draws_mat <- fit_pp %>% select(-lp__) %>% as_draws_matrix()
    colnames(draws_mat)[colnames(draws_mat) == "dispersion"] <- "sigma"
    qrows <- rbind(
      matrix_to_quantile_df(cong, rep, "PP", draws_mat, QUANTILE_LEVELS),
      degenerate_quantile_rows(cong, rep, "PP", "eta", eta_hat, QUANTILE_LEVELS)
    )
    list(qrows = qrows, elapsed = proc.time()[["elapsed"]] - t0, error = NA_character_)
  }, error = function(e) {
    list(qrows = NULL, elapsed = proc.time()[["elapsed"]] - t0, error = conditionMessage(e))
  })
  if (!is.null(pp_outcome$qrows)) quantile_df <- rbind(quantile_df, pp_outcome$qrows)
  if (is.na(pp_outcome$error)) {
    progress("  PP: done (%.1fs)", pp_outcome$elapsed)
  } else {
    progress("  PP: FAILED (%.1fs): %s", pp_outcome$elapsed, pp_outcome$error)
  }

  ## ---- NPP via hdbayes::glm.npp() ----
  t0 <- proc.time()[["elapsed"]]
  npp_outcome <- tryCatch({
    a0_grid <- seq(0, 1, length.out = N_A0_GRID)
    logncfun <- function(a0, ...) {
      hdbayes::glm.npp.lognc(formula = formula, family = family, histdata = hist_data,
                              a0 = a0, ...)
    }
    progress("  NPP: fitting a0 grid (%d points, %d workers)...", N_A0_GRID, ncores)
    cl <- makeCluster(ncores)
    on.exit(stopCluster(cl), add = TRUE)
    clusterSetRNGStream(cl, seed_here)
    clusterExport(cl, varlist = c("formula", "family", "hist_data"), envir = environment())
    # Dispatched in batches of `ncores` (rather than one parLapply() call
    # over the whole grid) purely so progress can be reported as each batch
    # completes -- parLapply already assigns work statically across the
    # cluster, so batching this way costs negligible extra overhead relative
    # to each batch's Stan fit time, but turns "silence for the whole grid"
    # into a running count.
    a0_batches <- split(seq_along(a0_grid), ceiling(seq_along(a0_grid) / ncores))
    a0.lognc_list <- vector("list", length(a0_grid))
    for (batch_idx in a0_batches) {
      batch_res <- parLapply(cl = cl, X = a0_grid[batch_idx], fun = logncfun,
                              iter_warmup = iter_warmup, iter_sampling = 5000,
                              chains = chains, refresh = 0)
      a0.lognc_list[batch_idx] <- batch_res
      progress("  NPP: a0 grid %d/%d points done", max(batch_idx), length(a0_grid))
    }
    stopCluster(cl)
    a0.lognc <- data.frame(do.call(rbind, a0.lognc_list))

    progress("  NPP: a0 grid done, fitting final NPP model (%d chains)...", chains)
    fit_npp <- glm.npp(
      formula = formula, family = family, data.list = data_list,
      a0.lognc = a0.lognc$a0, lognc = matrix(a0.lognc$lognc, ncol = 1),
      iter_warmup = iter_warmup, iter_sampling = iter_sampling,
      chains = chains, parallel_chains = ncores, refresh = 0
    )

    fit_df <- fit_npp %>% as_draws_df()
    param_names <- setdiff(
      names(fit_df),
      c(".chain", ".iteration", ".draw", "lp__", "logit_a0s[1]", "a0_hist_1")
    )
    draws_mat <- as.matrix(fit_df[, param_names])
    colnames(draws_mat)[colnames(draws_mat) == "dispersion"] <- "sigma"
    eta_draws <- fit_df$a0_hist_1

    qrows <- rbind(
      matrix_to_quantile_df(cong, rep, "NPP", draws_mat, QUANTILE_LEVELS),
      quantile_rows(cong, rep, "NPP", "eta", eta_draws, QUANTILE_LEVELS)
    )
    list(
      qrows = qrows, npp_eta_median = median(eta_draws),
      elapsed = proc.time()[["elapsed"]] - t0, error = NA_character_
    )
  }, error = function(e) {
    list(qrows = NULL, npp_eta_median = NA_real_,
         elapsed = proc.time()[["elapsed"]] - t0, error = conditionMessage(e))
  })
  if (!is.null(npp_outcome$qrows)) quantile_df <- rbind(quantile_df, npp_outcome$qrows)
  if (is.na(npp_outcome$error)) {
    progress("  NPP: done (%.1fs), eta_median=%.4f", npp_outcome$elapsed, npp_outcome$npp_eta_median)
  } else {
    progress("  NPP: FAILED (%.1fs): %s", npp_outcome$elapsed, npp_outcome$error)
  }

  eta_df <- rbind(eta_df, data.frame(
    congruence = cong, replicate = rep, eta_hat = eta_hat,
    pp_elapsed_sec = pp_outcome$elapsed, npp_elapsed_sec = npp_outcome$elapsed,
    pp_error = pp_outcome$error, npp_error = npp_outcome$error,
    stringsAsFactors = FALSE
  ))
  save_checkpoint()

  pp_status <- if (is.na(pp_outcome$error)) {
    sprintf("ok (%.1fs)", pp_outcome$elapsed)
  } else {
    sprintf("FAILED: %s", pp_outcome$error)
  }
  npp_status <- if (is.na(npp_outcome$error)) {
    sprintf("ok (%.1fs), eta_median=%.4f", npp_outcome$elapsed, npp_outcome$npp_eta_median)
  } else {
    sprintf("FAILED: %s", npp_outcome$error)
  }
  progress("[%d/%d] %s | replicate %d | COMBO DONE | PP: %s | NPP: %s",
           j, n_jobs, cong, rep, pp_status, npp_status)
}

## ---------------------------------------------------------------------------
## Summary
## ---------------------------------------------------------------------------

cat("\n--- Summary ---\n")
cat(sprintf("Combos processed: %d / %d\n", nrow(eta_df), n_jobs))
cat(sprintf("PP failures: %d | NPP failures: %d\n",
            sum(!is.na(eta_df$pp_error)), sum(!is.na(eta_df$npp_error))))

if (0.5 %in% QUANTILE_LEVELS) {
  eta_medians <- quantile_df %>% filter(parameter == "eta", quantile_level == 0.5)
  if (nrow(eta_medians) > 0) {
    cat("\nMedian eta (quantile_level = 0.5) -- PP (= eta_hat, degenerate) vs NPP posterior median, by congruence:\n")
    print(
      eta_medians %>%
        group_by(congruence, method) %>%
        summarise(n = n(), median_eta = median(value), .groups = "drop")
    )
  }
} else {
  cat("\n(0.5 is not in QUANTILE_LEVELS -- skipping the median-eta summary print; quantile_df still has everything needed for WIS.)\n")
}
