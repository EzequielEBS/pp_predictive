## code/linear_regression/eta_estimation_by_replicate.R
##
## Fits the eta-estimation algorithm (sequential_K_for_eta, i.e. estimate_eta_glm
## with an adaptively grown K) for EVERY replicate, at EVERY congruence level,
## for the p=3 linear regression simulated data. This extends sample.r, which
## only ran replicate == 1 per congruence level.
##
## COST WARNING: this is num_sim x 3 congruence levels = up to 600 combos,
## and each combo requires its own Bayesian glm.pp fit (4 chains x
## (iter_warmup + iter_sampling) draws) BEFORE sequential_K_for_eta can even
## start. At your existing settings (5000 warmup + 10000 sampling, 4 chains)
## this can run for hours. epsilon (below) is the main cost/precision lever
## for the eta-estimation half of each combo -- smaller epsilon means
## sequential_K_for_eta will grow K (and rerun estimate_eta_glm) more before
## it's satisfied. Tune NUM_SIM and epsilon down for a first test run.
##
## CHECKPOINTING: results are saved to RESULTS_PATH after every single combo,
## and on a rerun any (congruence, replicate) pair already present there is
## skipped. You can safely stop this script (Ctrl+C, or let it crash) and
## rerun it later to pick up where it left off.
##
## NOTE: a FAILED combo is also recorded in the checkpoint (with its error
## message in the `error` column), so a rerun will NOT automatically retry
## it -- a failure surfaces in the end-of-run summary for you to look at,
## rather than being silently retried and possibly burning more compute on
## a systematic problem. To retry specific failures, remove their rows from
## the saved results_df (readRDS(RESULTS_PATH), filter out the ones you want
## to redo, saveRDS back) before rerunning.

library(dplyr)
library(tidyverse)
library(hdbayes)
library(posterior)
library(subsampling)
library(readr)
library(parallel)

source("code/aux_fun_inf_match_glm.r")

set.seed(20260819)

## ---------------------------------------------------------------------------
## Configuration
## ---------------------------------------------------------------------------

data_dir <- Sys.getenv("ONPP_DATA_DIR", "../../onpp/simulated-data-regression/data")

NUM_SIM <- 200   # replicates per congruence level (max available: 500). Set
                 # this lower (e.g. 20) for a first test run before committing
                 # to the full 200.

CONGRUENCE_FILES <- c(
  "High congruence"  = "sim_data_high_cong_varying_beta_p3_n0_ge_n.csv",
  "Small congruence" = "sim_data_small_cong_varying_beta_p3_n0_ge_n.csv",
  "No congruence"    = "sim_data_no_cong_varying_beta_p3_n0_ge_n.csv"
)

formula <- y ~ X1 + X2
family  <- gaussian()

ncores        <- detectCores() - 1
chains        <- 4
iter_warmup   <- 5000
iter_sampling <- 10000

# sequential_K_for_eta settings -- epsilon is the main cost/precision knob
# for the eta-estimation half of each combo (see COST WARNING above).
K0        <- 50
epsilon   <- 0.05
gamma     <- 0.95
BASE_SEED <- 23082026

RESULTS_PATH <- "samples/linear_regression/eta_by_replicate.RData"
if (!dir.exists(dirname(RESULTS_PATH))) dir.create(dirname(RESULTS_PATH), recursive = TRUE)

## ---------------------------------------------------------------------------
## Load data
## ---------------------------------------------------------------------------

congruence_data <- lapply(CONGRUENCE_FILES, function(fname) {
  read_csv(file.path(data_dir, fname), show_col_types = FALSE) %>%
    filter(replicate <= NUM_SIM)
})
names(congruence_data) <- names(CONGRUENCE_FILES)

## ---------------------------------------------------------------------------
## Checkpoint load / init
## ---------------------------------------------------------------------------

if (file.exists(RESULTS_PATH)) {
  results_df <- readRDS(RESULTS_PATH)
  cat(sprintf("Resuming: %d combos already completed in %s\n",
              nrow(results_df), RESULTS_PATH))
} else {
  results_df <- data.frame(
    congruence = character(), replicate = integer(),
    n_hist = integer(), n_curr = integer(),
    eta_hat = double(), lower = double(), upper = double(),
    K_final = integer(), elapsed_sec = double(),
    error = character(), stringsAsFactors = FALSE
  )
}

is_done <- function(cong, rep) {
  any(results_df$congruence == cong & results_df$replicate == rep)
}

save_checkpoint <- function() {
  saveRDS(results_df, RESULTS_PATH)
}

## ---------------------------------------------------------------------------
## Main loop: sequential over (congruence, replicate) -- deliberately NOT
## parallelized at this level. glm.pp() already parallelizes its chains
## (parallel_chains = ncores) and sequential_K_for_eta()/estimate_eta_glm()
## already spins up its own cluster internally (ncores workers) for the K
## subsample replicates. Parallelizing this outer loop too would nest
## clusters inside clusters and oversubscribe the machine.
## ---------------------------------------------------------------------------

jobs <- expand.grid(
  congruence = names(CONGRUENCE_FILES), replicate = seq_len(NUM_SIM),
  stringsAsFactors = FALSE
)

n_jobs <- nrow(jobs)
cat(sprintf("Total combos: %d (%d congruence levels x %d replicates)\n",
            n_jobs, length(CONGRUENCE_FILES), NUM_SIM))

for (j in seq_len(n_jobs)) {
  cong <- jobs$congruence[j]
  rep  <- jobs$replicate[j]

  if (is_done(cong, rep)) next

  t0 <- proc.time()[["elapsed"]]

  row <- tryCatch({
    data_c <- congruence_data[[cong]]
    hist_data <- data_c %>% filter(data_id == "hist3", replicate == rep)
    curr_data <- data_c %>% filter(data_id == "current", replicate == rep)

    data_list <- list(curr_data, hist_data)

    # a0.vals = 0: posterior for beta under curr_data alone (no borrowing),
    # matching sample.r -- these draws feed estimate_eta_glm's score-ratio
    # computation.
    fit0 <- glm.pp(
      formula = formula, family = family, data.list = data_list,
      a0.vals = 0,
      iter_warmup = iter_warmup, iter_sampling = iter_sampling,
      chains = chains, parallel_chains = ncores,
      refresh = 0
    )
    beta_draws <- fit0 %>% select(-lp__, -dispersion) %>% as_draws_matrix()

    # a distinct seed per (congruence, replicate) combo so runs aren't
    # silently correlated with each other
    seed_here <- BASE_SEED + which(names(CONGRUENCE_FILES) == cong) * 1e6 + rep

    est <- sequential_K_for_eta(
      formula = formula, curr_data = curr_data, hist_data = hist_data,
      beta_draws = beta_draws, family = family,
      K0 = K0, epsilon = epsilon, gamma = gamma,
      SEED0 = seed_here, verbose = FALSE
    )

    data.frame(
      congruence = cong, replicate = rep,
      n_hist = nrow(hist_data), n_curr = nrow(curr_data),
      eta_hat = est$hat_median, lower = est$lower, upper = est$upper,
      K_final = est$K_final, elapsed_sec = proc.time()[["elapsed"]] - t0,
      error = NA_character_, stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      congruence = cong, replicate = rep,
      n_hist = NA_integer_, n_curr = NA_integer_,
      eta_hat = NA_real_, lower = NA_real_, upper = NA_real_,
      K_final = NA_integer_, elapsed_sec = proc.time()[["elapsed"]] - t0,
      error = conditionMessage(e), stringsAsFactors = FALSE
    )
  })

  results_df <- rbind(results_df, row)
  save_checkpoint()

  status <- if (is.na(row$error)) {
    sprintf("eta_hat=%.4f (K=%d)", row$eta_hat, row$K_final)
  } else {
    sprintf("FAILED: %s", row$error)
  }
  cat(sprintf("[%d/%d] %s | replicate %d | %.1fs | %s\n",
              j, n_jobs, cong, rep, row$elapsed_sec, status))
}

## ---------------------------------------------------------------------------
## Summary
## ---------------------------------------------------------------------------

cat("\n--- Summary ---\n")
ok <- results_df %>% filter(is.na(error))
failed <- results_df %>% filter(!is.na(error))
cat(sprintf("Completed: %d / %d combos (%d failed)\n",
            nrow(ok), nrow(results_df), nrow(failed)))
if (nrow(ok) > 0) {
  print(
    ok %>%
      group_by(congruence) %>%
      summarise(
        n = n(),
        median_eta_hat = median(eta_hat),
        mean_K_final = mean(K_final),
        .groups = "drop"
      )
  )
}
if (nrow(failed) > 0) {
  cat("\nFailed combos (see results_df$error for details):\n")
  print(failed %>% select(congruence, replicate, error))
}
