expected_sq_norm_score <- function(X, beta, beta_hat,  family = gaussian(), off = NULL) {
  if (is.function(family)) family <- family()
  if (is.null(off)) off <- rep(0, nrow(X))

  # Quantities at beta
  eta      <- X %*% beta + off
  mu_b     <- family$linkinv(eta)
  dmu_deta <- family$mu.eta(eta)
  V_b      <- family$variance(mu_b)
  M_diag   <- dmu_deta / V_b       # diagonal of M(beta)

  # Quantities at beta_hat
  eta_hat  <- X %*% beta_hat + off
  mu_bhat  <- family$linkinv(eta_hat)
  V_bhat   <- family$variance(mu_bhat)

  # Term 1: (1/phi^2) || X^T M(beta) (mu(beta_hat) - mu(beta)) ||^2
  XtM_diff <- t(X) %*% (M_diag * (mu_bhat - mu_b))
  term1    <- sum(XtM_diff^2)

  # Term 2: (1/phi) Tr( X^T M(beta) V(beta_hat) M(beta)^T X )
  XM    <- t(X) %*% diag(as.vector(M_diag))
  term2 <- sum(diag(XM %*% diag(as.vector(V_bhat)) %*% t(XM)))

  term1 + term2
}

expected_sq_norm_score_vec <- function(X, beta_mat, beta_hat, family = gaussian(), off = NULL) {
  if (is.function(family)) family <- family()
  if (is.null(off)) off <- rep(0, nrow(X))

  n <- nrow(X)
  S <- nrow(beta_mat)   # number of samples

  # Quantities at beta_hat (computed once) — n x 1
  eta_hat <- X %*% beta_hat + off
  mu_bhat <- family$linkinv(eta_hat)
  V_bhat  <- as.vector(family$variance(mu_bhat))

  # Quantities at all betas simultaneously
  # eta_mat: n x S  (each column is eta for one beta sample)
  eta_mat      <- X %*% t(beta_mat) + off        # n x S
  mu_mat       <- family$linkinv(eta_mat)            # n x S
  dmu_vec      <- family$mu.eta(eta_mat)             # n x S
  V_vec        <- family$variance(mu_mat)            # n x S
  
  dmu_deta_mat <- matrix(dmu_vec,     nrow = n)      # n x S
  V_mat        <- matrix(V_vec,       nrow = n)      # n x S
  M_mat        <- dmu_deta_mat / V_mat   

  # Term 1: (1/phi^2) || X^T M(beta) (mu(beta_hat) - mu(beta)) ||^2
  # diff_mat: n x S — each column is mu_bhat - mu(beta_s)
  diff_mat  <- as.vector(mu_bhat) - mu_mat          # n x S
  # X^T (M * diff): p x S
  XtM_diff  <- t(X) %*% (M_mat * diff_mat)          # p x S
  term1     <- colSums(XtM_diff^2)           # S x 1

  # Term 2: (1/phi) Tr( X^T M V(beta_hat) M^T X )
  # = (1/phi) sum_i V_bhat_i * || X^T e_i ||^2 * M_i^2
  # (XM)^2 %*% V_bhat: p x S then sum over p
  MV_mat <- sweep(M_mat^2, 1, V_bhat, "*")
  XM_sq   <- (t(X)^2) %*% MV_mat       # p x S  (broadcast V_bhat over S)
  term2   <- colSums(XM_sq)                    # S x 1

  term1 + term2
}

med_app <- function(x, gamma = 0.95){
  ## Aproximação normal usando o método Delta.
  dens <- density(x)
  app.pdf <- approxfun(dens)
  med.hat <- median(x)
  n <- length(x)
  sd.approx <- 1/(4 * n * (app.pdf(med.hat))^2)
  pars <- c(med.hat, sqrt(sd.approx))
  approx.ci <- qnorm(p = c(1 - gamma, 1 + gamma)/2, mean = pars[1], sd = pars[2])
  return(
  list(
    hat_median = med.hat, 
    lower = approx.ci[1],
    upper = approx.ci[2]
  ))
}

estimate_eta_glm <- function(
  formula,
  curr_data,
  hist_data,
  beta_draws, 
  family = gaussian(), 
  nu = .5,
  offset.list = NULL,
  ssp = F,
  K = 100,
  SEED = 1234,
  g = 0.95,
  ncores = detectCores() - 1
) {
  # Estimate beta_hat and beta0_hat using IRLS
  formula <- as.formula(formula)
  outcome_var <- all.vars(formula)[1]
  y <- curr_data[[outcome_var]]
  X <- model.matrix(formula, data = curr_data)
  y0 <- hist_data[[outcome_var]]
  X0 <- model.matrix(formula, data = hist_data)
  if (!is.null(offset.list)) {
    beta_hat <- 
      glm(formula, data = curr_data, family = family, offset = offset.list[[1]])$coefficients %>%
      as.matrix(ncol = 1)
    beta0_hat <-
      glm(formula, data = hist_data, family = family, offset = offset.list[[2]])$coefficients %>%
      as.matrix(ncol = 1)
  } else {
    beta_hat <- glm(formula, data = curr_data, family = family)$coefficients %>%
      as.matrix(ncol = 1)
    beta0_hat <- glm(formula, data = hist_data, family = family)$coefficients %>%
      as.matrix(ncol = 1)
  }

  k0 <- floor(nrow(X0)^(nu))
  cl <- makeCluster(ncores)
  on.exit(stopCluster(cl), add = TRUE)   # moved up: guarantees cleanup even if
                                          # clusterExport/parLapply below errors
  clusterSetRNGStream(cl, SEED)
  if (ssp) clusterEvalQ(cl, library(subsampling))  # ssp.glm() needs this loaded
                                                     # on each worker, not just
                                                     # the master session

  clusterExport(cl, varlist =
    c("formula",
      "curr_data",
      "hist_data",
      "beta_draws",
      "family",
      "expected_sq_norm_score_vec",
      "X0",
      "beta_hat",
      "beta0_hat",
      "offset.list",
      "ssp",
      "k0"
    ),
    envir = environment()
  )

  etas <- parLapply(cl, 1:K, function(i) {
    if (ssp) {
      if (family$family == "binomial") {
        ssp.results <- ssp.glm(
          formula  = formula,
          data     = hist_data,
          n.plt    = nrow(hist_data),
          n.ssp    = k0,
          family   = "quasibinomial",
          sampling.method = "withReplacement"
        )
      } else{
        ssp.results <- ssp.glm(
          formula  = formula,
          data     = hist_data,
          n.plt    = nrow(hist_data),
          n.ssp    = k0,
          family   = family$family,
          sampling.method = "withReplacement"
        )
      }
      idx0 <- ssp.results$index
    } else {
      idx0 <- sample(nrow(X0), size = k0, replace = FALSE)
    }
    X0k0 <- X0[idx0, , drop = FALSE]
    if (!is.null(offset.list)) {
      off <- offset.list[[2]][idx0]
    } else {
      off <- NULL
    }
    num <- expected_sq_norm_score_vec(X0k0, beta_draws, beta_hat, family, off = off)
    den <- expected_sq_norm_score_vec(X0k0, beta_draws, beta0_hat, family, off = off)

    hat_eta <- exp(0.5 * (log(mean(num)) - log(mean(den))))
    return(hat_eta)
  }) %>%
    unlist()

  # c() (not list()) so hat_median/lower/upper end up alongside etas at the
  # TOP level of the returned list -- list(etas=etas, med_app(...)) would
  # nest med_app's list as an unnamed 2nd element instead, so res$hat_median
  # etc. would silently be NULL.
  return(c(list(etas = etas), med_app(etas, gamma = g)))
}


#' Given estimate_eta_glm()'s returned c(median, lower, upper) and the K it
#' was run with, recover the implied required K for a target relative
#' error epsilon at confidence gamma. Same reverse-engineering trick as
#' before: med_app()'s internal sd.approx = sigma^2 / K, recovered exactly
#' from the CI via qnorm's inverse (no need to see the raw etas).
n_required_from_eta_result <- function(res, K, epsilon, gamma) {
  m_hat <- res$hat_median
  upper <- res$upper
  z <- qnorm((1 + gamma) / 2)
  sd_hat <- (upper - m_hat) / z
  sigma2_hat <- sd_hat^2 * K
  n_required <- z^2 * sigma2_hat / (epsilon^2 * m_hat^2)
  list(m_hat = m_hat, sd_hat = sd_hat, sigma2_hat = sigma2_hat,
       n_required = n_required)
}
 
## ---- sequential driver: re-call estimate_eta_glm() with a bigger K --------

#' Grow K by directly re-running estimate_eta_glm() until the K it was run
#' with is consistent with what its own output says is needed.
#'
#' Each iteration is a FULL rerun of estimate_eta_glm() (refits both GLMs,
#' regenerates all K replicates from scratch) -- this is the cost of
#' treating it as a black box rather than exposing draw_one_eta()
#' separately. SEED is bumped each call so retries aren't silently
#' correlated with the previous run's random draws.
#'
#' @param formula,curr_data,hist_data,beta_draws,family,nu,offset.list,ssp,ncores
#'   passed straight through to estimate_eta_glm() -- same names, same
#'   defaults, so a call here mirrors a call to estimate_eta_glm() directly.
#' @param K0 starting K
#' @param epsilon target relative error
#' @param gamma target confidence level (passed through as g)
#' @param confirm_frac look-ahead confirmation batch size as a fraction of
#'   K_current (K is bumped by at least this much for one confirming rerun
#'   once the sufficiency check first passes)
#' @param max_iter safety cap on iterations
sequential_K_for_eta <- function(formula, curr_data, hist_data, beta_draws,
                                  family = gaussian(), nu = .5,
                                  offset.list = NULL, ssp = FALSE,
                                  ncores = detectCores() - 1,
                                  K0 = 100, epsilon, gamma = 0.95,
                                  SEED0 = 1234, confirm_frac = 0.1,
                                  max_iter = 15, verbose = TRUE) {
  K_current <- K0
  seed_offset <- 0
  history <- data.frame(iter = integer(), K = integer(), m_hat = double(),
                         n_required = double(), action = character(),
                         stringsAsFactors = FALSE)

  run_once <- function(K, seed) {
    estimate_eta_glm(
      formula = formula, curr_data = curr_data, hist_data = hist_data,
      beta_draws = beta_draws, family = family, nu = nu,
      offset.list = offset.list, ssp = ssp, K = K, SEED = seed, g = gamma,
      ncores = ncores
    )
  }

  for (iter in seq_len(max_iter)) {
    seed_offset <- seed_offset + 1
    res <- run_once(K_current, SEED0 + seed_offset)
    est <- n_required_from_eta_result(res, K_current, epsilon, gamma)
    sufficient <- est$n_required <= K_current
 
    if (!sufficient) {
      if (verbose) cat(sprintf(
        "iter %2d: K=%5d  m_hat=%.4f  n_required=%7.1f  [growing]\n",
        iter, K_current, est$m_hat, est$n_required))
      history <- rbind(history, data.frame(
        iter = iter, K = K_current, m_hat = est$m_hat,
        n_required = est$n_required, action = "grow",
        stringsAsFactors = FALSE))
      K_current <- ceiling(est$n_required)
      next
    }
 
    # Looks sufficient -- confirm with one more rerun at a slightly larger K
    K_confirm <- K_current + max(10, ceiling(confirm_frac * K_current))
    seed_offset <- seed_offset + 1
    res_check <- run_once(K_confirm, SEED0 + seed_offset)
    est_check <- n_required_from_eta_result(res_check, K_confirm, epsilon, gamma)
    confirmed <- est_check$n_required <= K_confirm
 
    if (verbose) cat(sprintf(
      "iter %2d: K=%5d  m_hat=%.4f  n_required=%7.1f  [checking sufficiency at K=%d -> n_required=%.1f, %s]\n",
      iter, K_current, est$m_hat, est$n_required, K_confirm,
      est_check$n_required, if (confirmed) "confirmed" else "not confirmed"))
 
    history <- rbind(history, data.frame(
      iter = iter, K = K_current, m_hat = est$m_hat,
      n_required = est$n_required,
      action = if (confirmed) "confirmed" else "look-ahead-failed",
      stringsAsFactors = FALSE))
 
    if (confirmed) {
      if (verbose) cat(sprintf("\nConverged after %d iterations. K_final = %d\n",
                                iter, K_confirm))
      # c() so etas/hat_median/lower/upper from res_check land at the TOP
      # level alongside K_final/history -- e.g. K_est$etas works directly,
      # matching how actg_sample.r uses the result.
      return(c(list(K_final = K_confirm, history = history), res_check))
    }
    K_current <- K_confirm
  }

  warning("max_iter reached without confirmed convergence; returning last result")
  c(list(K_final = K_current, history = history), res)
}
