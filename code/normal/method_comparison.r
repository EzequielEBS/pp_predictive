# code/normal/method_comparison.r
#
# End-to-end comparison of eta-estimation methods, and of their downstream
# effect on the power-prior (PP) posterior AND the normalized-power-prior
# (NPP) posterior, for the normal-mean / fixed-variance toy example
# (Section 3.2.1: theta0 = theta1 = theta, y0i ~ N(theta0, sigma^2),
# yi ~ N(theta1, sigma^2), prior theta ~ N(m0, v0)).
#
# Reformulated into two clearly separate phases:
#
#   PHASE 1 (generate_all_replicates()): simulate every (scenario,
#   replicate) dataset UP FRONT -- nothing here estimates anything, it just
#   produces a flat list of "jobs" (one historical + one current sample
#   each). Kept separate from phase 2 so the random data-generating step
#   never competes with, or is reordered by, the parallel loop below.
#
#   PHASE 2 (run_replicates(), via process_one_replicate()): for EACH job,
#   in the SAME loop --
#     (i)   compute eta four ways:
#             - eta_exact   : the "oracle" eta -- the same influence-
#                             matching ratio used by estimate_eta(), but
#                             evaluated at the TRUE generating means
#                             (mu0, mu1) instead of their MLEs;
#             - eta_hat     : estimate_eta() -- the MLE-based influence-
#                             matching estimator (aux_fun_normal.R);
#             - eta_hat_glm : estimate_eta_glm() -- the score-matching /
#                             subsampling estimator (aux_fun_inf_match_glm.r);
#             - eta_a0_star : a0_star_normal() -- the closed-form a0*, built
#                             from the average historical-vs-current
#                             Kullback-Leibler divergence at the TRUE
#                             generating means (mu0, mu1), like eta_exact;
#     (ii)  compute the EXACT PP posterior of theta, in closed form, for
#           each of the four eta values;
#     (iii) compute the NPP posterior of theta -- eta given a Beta(a_tilde,
#           b_tilde) prior instead of plugged in -- as joint (eta, theta)
#           draws from code/normal/normal_npp_eta_fixedvar.stan via Stan.
#   That loop can run in parallel across jobs (a single outer
#   parallel::makeCluster(), cross-platform).
#
# NPP with a FIXED/plugged-in eta is NOT separately implemented -- for a
# fixed eta its theta-posterior is mathematically identical to PP's (the
# NPP normalizing constant c(eta) does not depend on theta, so it cancels
# once the posterior is renormalized). NPP only genuinely differs from PP
# once eta itself is treated as unknown/random, which is exactly what (iii)
# above computes.

suppressMessages({
  library(dplyr)
  library(purrr)
  library(parallel)
})
if (requireNamespace("pbapply", quietly = TRUE)) library(pbapply)

source("code/normal/aux_fun_normal.R")
source("code/aux_fun_inf_match_glm.r")

SEED <- 20260819

#-------------------------------------------------------------------------------
# eta estimator (exact/oracle) and the exact PP / NPP posteriors
#-------------------------------------------------------------------------------

# "Oracle" eta: reuses edelta_normal() (the same Delta(x) integral from the
# screenshot / aux_fun_normal.R) so it is provably the same formula as
# estimate_eta(), just evaluated at the true means instead of their MLEs.
#   eta = sqrt( Delta(mu1) / Delta(mu0) )
exact_eta_normal <- function(post_par, v, n0, mu0, mu1, nu = .5) {
  delta_curr <- edelta_normal(n0, v, post_par, mu1, nu = nu)
  delta_hist <- edelta_normal(n0, v, post_par, mu0, nu = nu)
  exp(0.5 * (log(delta_curr) - log(delta_hist)))
}

# a0*: a closed-form eta estimate that plugs the AVERAGE Kullback-Leibler
# divergence between the historical-data model and the current-data model
# (evaluated at every historical unit) directly into
#
#   a0* = d / (2 * N0 * Dbar0 + d),
#   Dbar0 = (1/N0) sum_{i=1}^{N0} D(p_theta0(.|x0i) || p_theta(.|x0i)),
#   D(f||g) = int f log(f/g)   (Kullback-Leibler divergence),
#
# where d is the number of parameters and N0 the historical sample size.
# Uses the TRUE generating means (mu0, mu1) for theta0/theta -- the same
# "oracle" values exact_eta_normal() above uses -- rather than their MLEs,
# so this is a0*'s oracle counterpart (like eta_exact vs. eta_hat), not an
# estimator computed from the data itself.
#
# For THIS model (normal mean, KNOWN/fixed variance v, d = 1 parameter),
# p_theta(x) = N(x; theta, v), so D(N(mu0, v) || N(mu1, v)) has the usual
# closed form (mu0 - mu1)^2 / (2v) -- and since there are no per-unit
# covariates in this toy model, that value doesn't actually vary with i, so
# Dbar0 collapses to that single number regardless of N0 (the sum-then-
# average in the formula above is a no-op here, but is kept explicit in
# case this is ever reused with a model where D DOES vary by unit).
a0_star_normal <- function(n0, v, mu0, mu1, d = 1) {
  Dbar0 <- mean(rep((mu0 - mu1)^2 / (2 * v), n0))
  d / (2 * n0 * Dbar0 + d)
}

# Power-prior posterior of theta, exact conjugate-normal formula: raise the
# historical likelihood to the power eta (equivalently, scale its
# per-observation variance from v to v/eta), fold it into the N(m0, v0)
# prior to get the PP prior, then update that prior with the current-data
# likelihood. eta = 0 reduces to "ignore the historical data" (PP prior =
# the original N(m0, v0) prior); eta = 1 reduces to full pooling.
pp_posterior_normal <- function(eta, m0, v0, v, y0, y) {
  eta <- max(eta, 0)
  pp_prior <- if (eta == 0) {
    list(m_star = m0, v_star = v0)
  } else {
    post_par_fixed_var(m0, v0, v / eta, y0)
  }
  post <- post_par_fixed_var(pp_prior$m_star, pp_prior$v_star, v, y)
  list(m = post$m_star, v = post$v_star,
       prior_m = pp_prior$m_star, prior_v = pp_prior$v_star)
}

# NPP: joint (eta, theta) posterior draws, eta ~ Beta(a_tilde, b_tilde),
# fitted in code/normal/normal_npp_eta_fixedvar.stan:
#
#   theta | eta, D, D0 ~ N(mu_tilde_eta, tau2_tilde_eta),
#     tau2_tilde_eta^-1 = tau2_eta^-1 + n/sigma2,
#     mu_tilde_eta = tau2_tilde_eta * (mu_eta/tau2_eta + n*ybar/sigma2),
#
#   pi_NPP(eta | D, D0) propto pi_A(eta) * Z(mu_tilde_eta, tau2_tilde_eta)
#                                          / Z(mu_eta, tau2_eta),
#   Z(mu, tau2) = (2*pi*tau2)^(1/2) * exp(mu^2 / (2*tau2)),
#
# where (mu_eta, tau2_eta) is the same update rule applied to the original
# prior (m0, v0) using only the eta-power-discounted historical data D0.
# Only eta is sampled by Stan (NUTS/HMC) -- theta is analytically
# marginalized out of the target -- and theta | eta, D, D0 is then drawn
# EXACTLY (closed form) in `generated quantities`, once per posterior draw
# of eta. So (eta, theta) together, as Stan returns them, already ARE joint
# draws from the NPP posterior p(theta, eta | D, D0).
#
# `model` can be passed in already-compiled (see run_replicates() below,
# which compiles it ONCE and reuses it across every job) -- compiling from
# scratch on every call is fine for a single one-off fit, but far too slow
# once this is called once per replicate.
npp_joint_stan_normal <- function(a_tilde = 1, b_tilde = 1, m0 = 0, v0 = 1, v = 1,
                                   y0, y, iter_warmup = 2000, iter_sampling = 2000,
                                   chains = 4, seed = SEED,
                                   stan_file = "code/normal/normal_npp_eta_fixedvar.stan",
                                   model = NULL) {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop("npp_joint_stan_normal() needs cmdstanr (and a compiled cmdstan toolchain) -- ",
         "see code/normal/samples_normal.R for how this repo already uses it.")
  }
  if (is.null(model)) model <- cmdstanr::cmdstan_model(stan_file)
  stan_data <- list(n0 = length(y0), n = length(y), y0 = y0, y = y,
                     m0 = m0, v0 = v0, v = v, tilde_a = a_tilde, tilde_b = b_tilde)
  fit <- model$sample(data = stan_data, iter_warmup = iter_warmup,
                       iter_sampling = iter_sampling, chains = chains,
                       parallel_chains = chains, seed = seed, refresh = 0)
  draws <- fit$draws(c("eta", "theta")) %>% posterior::as_draws_df()
  list(draws = draws, fit = fit,
       eta_mean = mean(draws$eta), eta_sd = sd(draws$eta),
       theta_mean = mean(draws$theta), theta_sd = sd(draws$theta))
}

has_cmdstanr <- requireNamespace("cmdstanr", quietly = TRUE)
if (!has_cmdstanr) {
  message("cmdstanr not installed -- run_replicates() below will compute PP only ",
          "(no NPP draws), and the NPP rows/curves in the medians table and density ",
          "plots further down will be skipped.")
}

#-------------------------------------------------------------------------------
# scenarios
#-------------------------------------------------------------------------------

# scenarios: named list of mu0 values; mu1 is held fixed at 1 throughout,
# matching the convention used elsewhere in code/normal/ (sim_helpers.R,
# inf_match_comp.r, inf_match_conv.r) where "discrepancy" means |mu0 - mu1|.
default_scenarios <- list(
  "No discrepancy"    = 1,
  "Small discrepancy" = 0.5,
  "Large discrepancy" = 0
)

# One fresh seed per (scenario, replicate), drawn up front so results are
# reproducible regardless of apply/parallel ordering -- shared by
# generate_all_replicates() below (and, by extension, by every replicate's
# eta/PP/NPP computation in run_replicates(), since a job's seed is used to
# re-seed before each of its own random steps too).
make_replicate_seeds <- function(scenarios, n_rep, SEED) {
  set.seed(SEED)
  matrix(
    sample.int(.Machine$integer.max, length(scenarios) * n_rep),
    nrow = n_rep, ncol = length(scenarios),
    dimnames = list(NULL, names(scenarios))
  )
}

#-------------------------------------------------------------------------------
# PHASE 1: generate every (scenario, replicate) dataset up front
#-------------------------------------------------------------------------------

# Returns a flat list of "jobs", one per (scenario, replicate): the
# scenario label, mu0/mu1, the replicate index, the seed used to generate
# it, and the simulated data (kept as the full hist/curr data frame, since
# estimate_eta() needs that shape -- not just y0/y vectors). Nothing here
# estimates anything; that's entirely phase 2 (process_one_replicate()).
generate_all_replicates <- function(scenarios = default_scenarios, mu1 = 1,
                                     n_rep = 100, n = 100, r = 1.5, v = 1,
                                     m0 = 0, v0 = 1, SEED = 20260819) {
  n0 <- r * n
  seeds <- make_replicate_seeds(scenarios, n_rep, SEED)

  jobs <- lapply(names(scenarios), function(label) {
    mu0 <- scenarios[[label]]
    lapply(seq_len(n_rep), function(j) {
      seed <- seeds[j, label]
      set.seed(seed)
      data <- generate_normal_data(n0 = n0, n = n, mu0 = mu0, mu = mu1,
                                    sigma0 = v, sigma = v)
      list(scenario = label, mu0 = mu0, mu1 = mu1, replicate = j,
           seed = seed, data = data)
    })
  })
  unlist(jobs, recursive = FALSE)
}

#-------------------------------------------------------------------------------
# PHASE 2: eta estimates + PP posteriors + NPP posterior, per job
#-------------------------------------------------------------------------------

# One (scenario, replicate) job -> both PP and (optionally) NPP posteriors,
# on the SAME simulated data -- this is the "same loop" the eta estimates,
# the PP posteriors, and the NPP posterior are all computed in.
#
# NPP is computed only when `model` (a compiled normal_npp_eta_fixedvar.stan
# CmdStanModel) is supplied; pass model = NULL to get PP-only results (e.g.
# when cmdstanr isn't installed).
process_one_replicate <- function(job, v, m0, v0, nu, alpha,
                                   M_beta, K_glm, ncores_glm,
                                   a_tilde = 1, b_tilde = 1,
                                   iter_warmup = 500, iter_sampling = 500,
                                   chains = 2,
                                   stan_file = "code/normal/normal_npp_eta_fixedvar.stan",
                                   model = NULL) {
  # Re-seed here (not just once at data-generation time in phase 1) so
  # every downstream random step -- beta_draws, estimate_eta_glm()'s own
  # cluster RNG stream, Stan's `seed` argument -- is reproducible
  # regardless of which worker, or in what order, this job runs on.
  set.seed(job$seed)

  # NOTE: `data` is both the data frame and (inside dplyr::filter() calls
  # below) the name of its "hist"/"curr" indicator column -- this shadowing
  # is inherited as-is from aux_fun_normal.R::estimate_eta() and every
  # other script in code/normal/ that uses generate_normal_data().
  data <- job$data
  hist_data <- data %>% filter(data == "hist")
  curr_data <- data %>% filter(data == "curr")
  y0 <- hist_data %>% pull(y)
  y  <- curr_data %>% pull(y)
  n0 <- length(y0)

  # Posterior of theta from the current data alone (N(m0, v0) prior) --
  # this is the p(theta | D, sigma^2) ~ N(m, v) of Section 3.2.1, which
  # estimate_eta()/exact_eta_normal() score the historical data against,
  # and which beta_draws (below) are drawn from for estimate_eta_glm().
  post_par <- post_par_fixed_var(m0, v0, v, y)

  # (i) eta, four ways
  eta_exact <- exact_eta_normal(post_par, v, n0, job$mu0, job$mu1, nu = nu)
  eta_hat   <- estimate_eta(data, post_par, v, alpha = alpha, mle = TRUE, nu = nu)

  beta_draws <- matrix(
    rnorm(M_beta, mean = post_par$m_star, sd = sqrt(post_par$v_star)),
    ncol = 1
  )
  eta_hat_glm <- estimate_eta_glm(
    formula = y ~ 1, curr_data = curr_data, hist_data = hist_data,
    beta_draws = beta_draws, family = gaussian(), nu = nu,
    K = K_glm, SEED = job$seed, ncores = ncores_glm
  )$hat_median

  eta_a0_star <- a0_star_normal(n0, v, job$mu0, job$mu1, d = 1)

  etas <- c(exact = eta_exact, estimate_eta = eta_hat, estimate_eta_glm = eta_hat_glm,
            Reznik = eta_a0_star)

  # (ii) PP posterior of theta, for each of the four eta values
  pp_rows <- do.call(rbind, lapply(names(etas), function(nm) {
    eta <- unname(etas[nm])
    pp  <- pp_posterior_normal(eta, m0, v0, v, y0, y)
    data.frame(scenario = job$scenario, replicate = job$replicate,
               method = nm, eta = eta, pp_m = pp$m, pp_v = pp$v)
  }))

  # (iii) NPP: joint (eta, theta) draws via Stan, on the SAME (y0, y) --
  # only computed when a compiled model was supplied (run_replicates()
  # below compiles it once and passes it to every job).
  npp_draws <- if (!is.null(model)) {
    fit <- npp_joint_stan_normal(a_tilde = a_tilde, b_tilde = b_tilde,
                                  m0 = m0, v0 = v0, v = v, y0 = y0, y = y,
                                  iter_warmup = iter_warmup, iter_sampling = iter_sampling,
                                  chains = chains, seed = job$seed, stan_file = stan_file,
                                  model = model)
    draws <- as.data.frame(fit$draws)  # drop the posterior::draws_df class
    draws$scenario  <- job$scenario    # before rbind-ing across replicates
    draws$replicate <- job$replicate
    draws
  } else {
    NULL
  }

  list(pp = pp_rows, npp = npp_draws)
}

#-------------------------------------------------------------------------------
# driver: phase 1 (generate) then phase 2 (estimate), optionally in parallel
#-------------------------------------------------------------------------------

# ncores_glm defaults to 1, NOT max(1, detectCores() - 1): when ncores > 1,
# every job already runs on its own outer worker, and estimate_eta_glm()
# spins up its OWN cluster internally for its K_glm subsamples -- nesting a
# second layer of parallel workers inside every outer worker would
# over-subscribe the machine (ncores x ncores_glm processes) and pay
# cluster start/teardown overhead on every single replicate. Raise
# ncores_glm only if you deliberately want that nested parallelism (e.g.
# ncores = 1, so there is only one layer).
#
# Stan iterations/chains here are cut down from npp_joint_stan_normal()'s
# own single-fit defaults (4 chains x 2000+2000) since this now runs once
# PER REPLICATE (n_rep x length(scenarios) fits): see the note further down
# by the one-dataset comparison, which keeps the full-strength defaults for
# a single higher-quality reference fit per scenario instead.
run_replicates <- function(scenarios = default_scenarios, mu1 = 1,
                            n_rep = 100, n = 100, r = 1.5, v = 1,
                            m0 = 0, v0 = 1, nu = .5, alpha = 1/4,
                            M_beta = 10000, K_glm = 100, ncores_glm = 1,
                            a_tilde = 1, b_tilde = 1,
                            iter_warmup = 500, iter_sampling = 500, chains = 2,
                            stan_file = "code/normal/normal_npp_eta_fixedvar.stan",
                            SEED = 20260819,
                            ncores = max(1, detectCores() - 1)) {
  # PHASE 1 -- every dataset, generated up front.
  jobs <- generate_all_replicates(scenarios = scenarios, mu1 = mu1, n_rep = n_rep,
                                   n = n, r = r, v = v, m0 = m0, v0 = v0, SEED = SEED)

  compute_npp <- has_cmdstanr
  model <- if (compute_npp) cmdstanr::cmdstan_model(stan_file) else NULL

  # PHASE 2 -- eta + PP + NPP, all in the same call, for every job.
  run_job <- function(job) {
    process_one_replicate(job, v = v, m0 = m0, v0 = v0, nu = nu, alpha = alpha,
                           M_beta = M_beta, K_glm = K_glm, ncores_glm = ncores_glm,
                           a_tilde = a_tilde, b_tilde = b_tilde,
                           iter_warmup = iter_warmup, iter_sampling = iter_sampling,
                           chains = chains, stan_file = stan_file, model = model)
  }

  out <- if (ncores > 1) {
    cl <- makeCluster(ncores)
    on.exit(stopCluster(cl), add = TRUE)
    clusterEvalQ(cl, {
      suppressMessages({ library(dplyr); library(purrr); library(parallel) })
      source("code/normal/aux_fun_normal.R")
      source("code/aux_fun_inf_match_glm.r")
      requireNamespace("cmdstanr", quietly = TRUE)
      requireNamespace("posterior", quietly = TRUE)
      NULL
    })
    # process_one_replicate()/pp_posterior_normal()/exact_eta_normal()/
    # a0_star_normal()/npp_joint_stan_normal() are defined at top level in
    # this script, not in a sourced file -- unlike run_job()'s OWN captured
    # arguments (v, m0, model, ...), which travel automatically with its
    # closure, these need an explicit clusterExport() (matches how
    # estimate_eta_glm() itself already does this for its own
    # K_glm-subsample cluster).
    clusterExport(cl, varlist = c("process_one_replicate", "pp_posterior_normal",
                                   "exact_eta_normal", "a0_star_normal",
                                   "npp_joint_stan_normal"))
    lapply_fn <- if (requireNamespace("pbapply", quietly = TRUE)) {
      function(X, FUN) pbapply::pblapply(X, FUN, cl = cl)
    } else {
      function(X, FUN) parLapply(cl, X, FUN)
    }
    lapply_fn(jobs, run_job)
  } else {
    lapply_fn <- if (requireNamespace("pbapply", quietly = TRUE)) pbapply::pblapply else lapply
    lapply_fn(jobs, run_job)
  }

  results <- do.call(rbind, lapply(out, `[[`, "pp"))
  results$scenario <- factor(results$scenario, levels = names(scenarios))
  rownames(results) <- NULL

  npp_results <- if (compute_npp) {
    npp_results <- do.call(rbind, lapply(out, `[[`, "npp"))
    npp_results$scenario <- factor(npp_results$scenario, levels = names(scenarios))
    rownames(npp_results) <- NULL
    npp_results
  } else {
    NULL
  }

  list(results = results, npp_results = npp_results)
}

#-------------------------------------------------------------------------------
# run + verify + save
#-------------------------------------------------------------------------------
# n_rep/K_glm here are kept modest for a default smoke run -- estimate_eta_glm()
# is the dominant per-replicate cost (it fits 2 GLMs and spins up/tears down a
# 1-node cluster of K_glm score-matching subsamples every call -- see
# ncores_glm's note above run_replicates()), and NPP now adds a Stan fit on
# top of that for every replicate too. Increase n_rep/K_glm/iter_warmup/
# iter_sampling for a more precise comparison once you've confirmed the
# pipeline runs end-to-end; ncores defaults to using every core but one.

comparison <- run_replicates(
  scenarios = default_scenarios, mu1 = 1,
  n_rep = 100, n = 100, r = 1.5, v = 1, m0 = 0, v0 = 1,
  nu = .999999999, alpha = 1/4, M_beta = 1000, K_glm = 100, ncores_glm = 1,
  a_tilde = 1, b_tilde = 1, iter_warmup = 1000, iter_sampling = 1000, chains = 2,
  SEED = SEED, ncores = max(1, detectCores() - 1)
)
results     <- comparison$results
npp_results <- comparison$npp_results

# Summary: median eta and median PP posterior mean/variance by scenario x
# method (medians rather than means, matching the med_np()/med_app()
# convention used elsewhere in this codebase for eta summaries).
summary_tab <- results %>%
  group_by(scenario, method) %>%
  summarise(
    eta_median  = median(eta),
    eta_lower   = quantile(eta, .025),
    eta_upper   = quantile(eta, .975),
    pp_m_median = median(pp_m),
    pp_v_median = median(pp_v),
    .groups = "drop"
  )
print(summary_tab, n = Inf)

save(results, summary_tab, file = "samples/method_comparison_normal.RData")
if (requireNamespace("xtable", quietly = TRUE)) {
  print(xtable::xtable(summary_tab), include.rownames = FALSE)
}

if (has_cmdstanr) {
  # Per-replicate posterior median eta/theta, by scenario -- same shape as
  # summary_tab above, so the two can be read side by side: PP's posterior
  # (one point eta estimate per replicate) vs. NPP's posterior (a full
  # (eta, theta) draw set per replicate, eta itself estimated with
  # uncertainty via its Beta(a_tilde, b_tilde) prior).
  npp_summary_tab <- npp_results %>%
    group_by(scenario, replicate) %>%
    summarise(
      eta_mean  = mean(eta),
      theta_mean = mean(theta)
    ) %>%
    summarise(
      eta_median   = median(eta_mean),
      eta_lower    = quantile(eta_mean, .025),
      eta_upper    = quantile(eta_mean, .975),
      theta_median = median(theta_mean),
      theta_lower  = quantile(theta_mean, .025),
      theta_upper  = quantile(theta_mean, .975),
      .groups = "drop"
    )
  print(npp_summary_tab, n = Inf)

  save(results, npp_results, summary_tab, npp_summary_tab,
       file = "samples/method_comparison_normal.RData")
}

#===============================================================================
# Shared setup for the analysis sections below: theta's MSE/WIS-vs-truth
# table, and the medians table + density plots for both theta and eta.
#===============================================================================

MU1_TRUE <- 1  # matches default_scenarios' fixed mu1 (see its own comment above)

PP_METHOD_LABELS <- c(exact = "IMPP (exact)", estimate_eta = "IMPP (est. η)",
                       estimate_eta_glm = "IMPP (est. η GLM)", Reznik = "RPP (a0*)")

FIG_DIR <- "figures"
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

# Shared by every gt-table export below (theta's MSE/WIS panel and the
# medians table): gt tables don't natively lay out side by side within a
# single gt object, so each is rendered to self-contained HTML via
# gt::as_raw_html() and wrapped in a flex container -- same helper as
# code/linear_regression/pp_npp_results_table.R::combine_side_by_side_html().
# Defining these two here (rather than inside a requireNamespace("gt")
# block) is harmless even when gt isn't installed -- a function body isn't
# evaluated until it's called, and every call site below is itself guarded.
combine_side_by_side_html <- function(tables, path, page_title = "Table", gap_px = 20) {
  panels <- vapply(tables, function(t) as.character(gt::as_raw_html(t)), character(1))
  divs <- paste0('<div style="flex: 0 0 auto; margin: 0;">', panels, "</div>")
  page <- paste0(
    '<!doctype html><html><head><meta charset="utf-8">',
    sprintf("<title>%s</title></head>", page_title),
    '<body style="font-family: -apple-system, Helvetica, Arial, sans-serif; padding: 16px; margin: 0;">',
    sprintf(
      '<div id="panels-wrapper" style="display:flex; flex-wrap:wrap; gap:%dpx; align-items:flex-start; justify-content:flex-start; width:fit-content;">',
      gap_px
    ),
    paste(divs, collapse = "\n"),
    "</div></body></html>"
  )
  writeLines(page, path)
}

save_gt_png <- function(html_path, png_path) {
  tryCatch({
    if (!requireNamespace("webshot2", quietly = TRUE)) stop("webshot2 not installed")
    webshot2::webshot(url = html_path, file = png_path, selector = "#panels-wrapper",
                       vwidth = 1200, vheight = 900)
  }, error = function(e) {
    message("PNG export skipped/failed (HTML page was still saved): ", conditionMessage(e))
  })
}

#===============================================================================
# MSE / WIS analysis: theta vs. its known truth (mu1)
#===============================================================================
#
# theta's truth is mu1 (=1, the TRUE current-data-generating mean, fixed
# across every scenario/replicate -- see default_scenarios' own comment
# above), so an MSE/WIS-vs-truth table makes sense for it. (eta doesn't get
# one: unlike theta, eta has no fixed truth to score against replicate to
# replicate -- exact_eta_normal()'s own "oracle" value changes with every
# replicate's data, since it depends on that replicate's post_par$m_star --
# see the "Medians" section below instead.)
#
# Both PP (results) and NPP (npp_results) give a POINT estimate (PP: pp_m;
# NPP: the posterior median from its Stan draws) for Bias/MSE, and a
# predictive INTERVAL for WIS -- PP's from its closed-form Normal posterior
# N(pp_m, pp_v); NPP's from the empirical quantiles of its Stan draws.
#
# wis_score() below is copied verbatim from
# code/linear_regression/pp_npp_results_table.R (independently verified
# there against a closed-form degenerate case and an equivalent pinball-
# loss formula -- see that file's header for details). The gt-table-
# building code mirrors that same file's build_pp_npp_gt_table()/
# combine_side_by_side_html() pattern (generalized here to our variable
# number of method columns and single scenario grouping instead of
# congruence x parameter).

if (requireNamespace("gt", quietly = TRUE) && requireNamespace("tidyr", quietly = TRUE)) {

QUANTILE_LEVELS <- c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975)  # symmetric around .5, as wis_score() requires

#' WIS for one predictive quantile curve against one true value.
#' (verbatim copy of code/linear_regression/pp_npp_results_table.R's
#' wis_score() -- see that file for the derivation/verification notes.)
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

## ---- per-replicate score, PP (closed form) + NPP (Stan draws) -----------

theta_pp_scores <- results %>%
  rowwise() %>%
  mutate(
    wis = wis_score(QUANTILE_LEVELS,
                     qnorm(QUANTILE_LEVELS, mean = pp_m, sd = sqrt(pp_v)),
                     MU1_TRUE),
    error = pp_m - MU1_TRUE
  ) %>%
  ungroup() %>%
  transmute(scenario, replicate, method = unname(PP_METHOD_LABELS[method]), error, wis)

theta_npp_scores <- if (has_cmdstanr) {
  npp_results %>%
    group_by(scenario, replicate) %>%
    summarise(
      wis   = wis_score(QUANTILE_LEVELS, quantile(theta, probs = QUANTILE_LEVELS), MU1_TRUE),
      error = median(theta) - MU1_TRUE,
      .groups = "drop"
    ) %>%
    transmute(scenario, replicate, method = "NPP", error, wis)
} else {
  NULL
}

theta_summary <- bind_rows(theta_pp_scores, theta_npp_scores) %>%
  group_by(scenario, method) %>%
  summarise(Bias = abs(mean(error)), MSE = mean(error^2), WIS = mean(wis), .groups = "drop")

## ---- gt table builder, one row per scenario, one column per method -----
#
# Same convention as pp_npp_results_table.R::build_pp_npp_gt_table(): each
# ROW is colored relative to its own methods' values (not one scale across
# the whole table), via gt::data_color().
make_scenario_gt_table <- function(df, metric, title, method_cols, palette_fn) {
  wide <- df %>%
    select(scenario, method, value = all_of(metric)) %>%
    tidyr::pivot_wider(names_from = method, values_from = value) %>%
    mutate(scenario = factor(scenario, levels = names(default_scenarios))) %>%
    arrange(scenario) %>%
    rename(Scenario = scenario)

  tbl <- wide %>%
    gt::gt() %>%
    gt::fmt_number(columns = tidyselect::all_of(method_cols), decimals = 4) %>%
    gt::tab_header(title = title) %>%
    gt::tab_spanner(label = "Method", columns = tidyselect::all_of(method_cols)) %>%
    gt::tab_options(table.font.size = 11, table.margin.left = 0, table.margin.right = 0)

  for (i in seq_len(nrow(wide))) {
    row_vals <- as.numeric(wide[i, method_cols])
    tbl <- tbl %>%
      gt::data_color(
        columns = tidyselect::all_of(method_cols), rows = i,
        palette = palette_fn(100)[30:70],
        domain  = c(min(row_vals), max(row_vals))
      )
  }
  tbl
}

blue_palette <- function(n) grDevices::colorRampPalette(c("#eaf1fb", "#1f3a7a"))(n)
warm_palette <- function(n) grDevices::colorRampPalette(c("#f7b267", "#c9184a", "#6a1b6a"))(n)

theta_method_cols <- unname(c(PP_METHOD_LABELS, "NPP"))

gt_mse_theta <- make_scenario_gt_table(theta_summary, "MSE", "Average MSE (θ)",
                                        theta_method_cols, blue_palette)
gt_wis_theta <- make_scenario_gt_table(theta_summary, "WIS", "Average WIS (θ)",
                                        theta_method_cols, warm_palette)

## ---- combine MSE + WIS side by side, save as one PNG ---------------------
# (combine_side_by_side_html()/save_gt_png() are defined once, in the
# "Shared setup" section above, and reused here and by the medians table
# export below.)

THETA_HTML <- file.path(FIG_DIR, "method_comparison_theta_table.html")
THETA_PNG  <- file.path(FIG_DIR, "method_comparison_theta_table.png")

combine_side_by_side_html(list(gt_mse_theta, gt_wis_theta), THETA_HTML, page_title = "θ MSE / WIS")
save_gt_png(THETA_HTML, THETA_PNG)

cat("Saved θ MSE/WIS panels to ", THETA_HTML, "\n", sep = "")
print(gt_mse_theta)
print(gt_wis_theta)

} else {
  message("gt and/or tidyr not installed -- skipping the θ MSE/WIS table. ",
          "theta_summary-equivalent numbers are still in `results`/`npp_results` ",
          "directly if you want to compute them another way.")
}

#===============================================================================
# Medians (LaTeX table) and one-replicate density plots, for theta and eta
#===============================================================================
#
# An MSE/WIS-vs-truth framing doesn't fit eta well here (see the section
# above for why it works for theta but not eta), so this section gives, for
# eta specifically:
#
#   (a) a plain table of each method's MEDIAN eta by scenario -- no "truth"
#       involved, just a summary -- combining summary_tab (PP) and
#       npp_summary_tab (NPP), both already computed above, WIDE-format:
#       one column per method rather than one row per (scenario, method),
#       since with only eta (no theta alongside it any more) a method-per-
#       column layout reads more like a comparison at a glance. As a LaTeX
#       table AND as an HTML/PNG image (unlike the theta MSE/WIS panel
#       above, gt::data_color() is NOT used here -- there's no "truth" for
#       a color scale to highlight against, so the table stays plain);
#   (b) for ONE representative replicate per scenario, what each method's
#       actual (eta, theta) distribution LOOKS like: PP's closed-form
#       Normal curve (theta) / point estimate (eta) against NPP's Stan-
#       sampled density -- side by side, one panel per scenario.
#
# Plotting follows this repo's own established convention (see
# code/normal/ppc_normal_data.R and code/common/plot_theme.R): closed-form
# curves via geom_function(), Stan draws via geom_density(), point
# estimates via geom_vline(), combined with patchwork.

## ---- (a) table of median eta, one column per method -----------------------

if (requireNamespace("tidyr", quietly = TRUE)) {

medians_tab <- bind_rows(
  summary_tab %>%
    transmute(scenario, method = unname(PP_METHOD_LABELS[method]), eta_median),
  if (has_cmdstanr) {
    npp_summary_tab %>% transmute(scenario, method = "NPP", eta_median)
  } else {
    NULL
  }
)

eta_method_order <- c(unname(PP_METHOD_LABELS), "NPP")

medians_tab_display <- medians_tab %>%
  tidyr::pivot_wider(names_from = method, values_from = eta_median) %>%
  mutate(scenario = factor(scenario, levels = names(default_scenarios))) %>%
  arrange(scenario) %>%
  select(scenario, any_of(eta_method_order)) %>%
  rename(Scenario = scenario)

print(medians_tab_display, n = Inf)
if (requireNamespace("xtable", quietly = TRUE)) {
  print(xtable::xtable(medians_tab_display, digits = 4), include.rownames = FALSE)
}

# Same table, as an HTML/PNG image (no gt::data_color() -- this table has
# no "truth" to score against, so unlike the theta MSE/WIS panel above
# there's nothing meaningful for a color scale to highlight).
if (requireNamespace("gt", quietly = TRUE)) {
  eta_method_cols_present <- intersect(eta_method_order, names(medians_tab_display))

  gt_medians <- medians_tab_display %>%
    gt::gt() %>%
    gt::fmt_number(columns = tidyselect::all_of(eta_method_cols_present), decimals = 4) %>%
    gt::tab_header(title = "Median η by scenario and method") %>%
    gt::tab_options(table.font.size = 11, table.margin.left = 0, table.margin.right = 0)

  MEDIANS_HTML <- file.path(FIG_DIR, "method_comparison_medians_table.html")
  MEDIANS_PNG  <- file.path(FIG_DIR, "method_comparison_medians_table.png")

  combine_side_by_side_html(list(gt_medians), MEDIANS_HTML, page_title = "Medians")
  save_gt_png(MEDIANS_HTML, MEDIANS_PNG)

  cat("Saved medians panel to ", MEDIANS_HTML, "\n", sep = "")
  print(gt_medians)
} else {
  message("gt not installed -- skipping the medians table PNG/HTML export.")
}

} else {
  message("tidyr not installed -- skipping the median-eta table.")
}

## ---- (b) density plots for ONE replicate, θ and η, per scenario ---------

if (requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("patchwork", quietly = TRUE)) {
  suppressMessages({ library(ggplot2); library(patchwork) })
  if (file.exists("code/common/plot_theme.R")) source("code/common/plot_theme.R")

  REPLICATE_FOR_DENSITY <- 1  # any single replicate -- just needs to be the
                               # SAME one for PP and NPP within a scenario

  # NOTE: every value pulled out of a loop variable (mu_i/sd_i/eta_i/lab_i
  # below) is injected into aes()/geom_function() via !! or a fresh local()
  # environment, NOT referenced directly (e.g. aes(color = nm)) -- ggplot2
  # layers capture unevaluated expressions, so referencing a loop variable
  # directly would have every layer silently show the LOOP'S FINAL value
  # once the plot is actually built, not its own iteration's value. Both
  # fixes were verified independently (each produces a distinct curve/line
  # per method here, not all collapsed onto the last one) before use.
  make_theta_density_plot <- function(label) {
    pp_rows <- results %>% filter(scenario == label, replicate == REPLICATE_FOR_DENSITY)
    p <- ggplot()
    for (nm in names(PP_METHOD_LABELS)) {
      row   <- pp_rows %>% filter(method == nm)
      mu_i  <- row$pp_m
      sd_i  <- sqrt(row$pp_v)
      lab_i <- unname(PP_METHOD_LABELS[[nm]])
      p <- p + geom_function(
        fun = local({ mu <- mu_i; sdv <- sd_i; function(x) dnorm(x, mu, sdv) }),
        aes(color = !!lab_i), linewidth = 1
      )
    }
    if (has_cmdstanr) {
      npp_draws <- npp_results %>% filter(scenario == label, replicate == REPLICATE_FOR_DENSITY)
      p <- p + geom_density(data = npp_draws, aes(x = theta, color = "NPP"), linewidth = 1)
    }
    p +
      geom_vline(aes(xintercept = !!MU1_TRUE, color = "true"), linetype = "dashed", linewidth = 0.7) +
      labs(x = expression(theta), y = "", title = label) +
      scale_color_manual(
        name = NULL,
        values = c(setNames(c("#66A8D0", "#D06673", "#7f7f7f", "#5DA271"), unname(PP_METHOD_LABELS)),
                   "NPP" = "#D0C366", "true" = "black")
      ) +
      guides(color = guide_legend(nrow = 2)) +
      theme_bw() + theme_pp(legend.position = "bottom")
  }

  make_eta_density_plot <- function(label) {
    pp_rows <- results %>% filter(scenario == label, replicate == REPLICATE_FOR_DENSITY)
    p <- ggplot()
    if (has_cmdstanr) {
      npp_draws <- npp_results %>% filter(scenario == label, replicate == REPLICATE_FOR_DENSITY)
      p <- p + geom_density(data = npp_draws, aes(x = eta, color = "NPP"), linewidth = 1)
    }
    for (nm in names(PP_METHOD_LABELS)) {
      row   <- pp_rows %>% filter(method == nm)
      eta_i <- row$eta
      lab_i <- unname(PP_METHOD_LABELS[[nm]])
      p <- p + geom_vline(aes(xintercept = !!eta_i, color = !!lab_i), linewidth = 1)
    }
    p +
      xlim(c(0, 1)) +
      labs(x = expression(eta), y = "", title = label) +
      scale_color_manual(
        name = NULL,
        values = c(setNames(c("#66A8D0", "#D06673", "#7f7f7f", "#5DA271"), unname(PP_METHOD_LABELS)),
                   "NPP" = "#D0C366")
      ) +
      guides(color = guide_legend(nrow = 2)) +
      theme_bw() + theme_pp(legend.position = "bottom")
  }

  theta_plots <- lapply(names(default_scenarios), make_theta_density_plot)
  eta_plots   <- lapply(names(default_scenarios), make_eta_density_plot)

  # 2 rows (theta, eta) x length(default_scenarios) columns.
  density_grid <- patchwork::wrap_plots(c(theta_plots, eta_plots), ncol = length(default_scenarios))
  print(density_grid)

  FIG_DIR <- "figures"
  if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)
  ggsave(plot = density_grid, filename = file.path(FIG_DIR, "method_comparison_density.png"),
         width = 5 * length(default_scenarios), height = 8, dpi = 320)
} else {
  message("ggplot2 and/or patchwork not installed -- skipping the θ/η density plots.")
}

