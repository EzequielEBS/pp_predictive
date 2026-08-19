irls <- function(y, X, family = gaussian(), off = NULL, tol = 1e-8, maxit = 100) {
  if (is.function(family)) family <- family()
  if (is.null(off)) off <- rep(0, nrow(X))
  
  # Initialize
  mu  <- (y + mean(y)) / 2
  eta <- family$linkfun(mu)
  
  beta <- rep(0, ncol(X))
  
  for (i in seq_len(maxit)) {
    dmu_deta <- family$mu.eta(eta)
    V_mu     <- family$variance(mu)
    
    # Working response and weights
    z <- eta - off + (y - mu) / dmu_deta
    W <- as.vector(dmu_deta^2 / V_mu)
    
    # Weighted least squares
    beta_new <- lm.wfit(X, z, W)$coefficients
    
    # Update
    eta <- X %*% beta_new + off
    mu  <- family$linkinv(eta)
    
    if (max(abs(beta_new - beta)) < tol) break
    beta <- beta_new
  }
  
  beta_new %>% as.matrix(ncol = 1)
}

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

estimate_eta_glm <- function(
  formula,
  curr_data,
  hist_data,
  beta_draws, 
  family = gaussian(), 
  nu = .5,
  offset.list = NULL,
  ssp = F
) {
  # Estimate beta_hat and beta0_hat using IRLS
  formula <- as.formula(formula)
  outcome_var <- all.vars(formula)[1]
  y <- curr_data[[outcome_var]]
  X <- model.matrix(formula, data = curr_data)
  y0 <- hist_data[[outcome_var]]
  X0 <- model.matrix(formula, data = hist_data)
  # beta_hat <- irls(y, X, family, off = offset.list[[1]])
  # beta0_hat <- irls(y0, X0, family, off = offset.list[[2]])
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

  # print(sum((beta_hat - beta0_hat)^2))
  # print(beta_hat)
  # print(beta0_hat)

  k0 <- floor(nrow(X0)^(nu))
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
  X0k0 <- X0[idx0,]
  if (!is.null(offset.list)) {
    off <- offset.list[[2]][idx0]
  } else {
    off <- NULL
  }
  num <- expected_sq_norm_score_vec(X0k0, beta_draws, beta_hat, family, off = off)
  den <- expected_sq_norm_score_vec(X0k0, beta_draws, beta0_hat, family, off = off)

  hat_eta <- exp(0.5 * (log(mean(num)) - log(mean(den))))
  return(hat_eta)
}

compute_n <- function(
  l, 
  gamma,
  formula,
  curr_data,
  hist_data,
  beta_draws,
  family,
  max_iter = 1000,
  starting_M = 100,
  non_par = F,
  ncores = detectCores() - 1
) {
  M <- starting_M
  aux <- T
  while(aux) {
    print(paste0("M = ", M))
    cl <- makeCluster(ncores)

    clusterExport(cl, varlist = 
      c("estimate_eta_glm", 
        "formula", 
        "curr_data", 
        "hist_data", 
        "beta_draws", 
        "family",
        "expected_sq_norm_score_vec"
      ),
      envir = environment()
    )
    clusterEvalQ(cl, {
      library(dplyr)
      library(tidyr)
    })
    etas <- parLapply(cl, 1:M, function(i) {
      estimate_eta_glm(
        formula = formula,
        curr_data = curr_data,
        hist_data = hist_data,
        beta_draws = beta_draws, 
        family = family
      )
    }) %>% unlist()

    on.exit(stopCluster(cl))
    med_hat <- median(etas)
    dens <- density(etas)
    app_pdf <- approxfun(dens)
    sd_approx <- 1/(4 * M * (app_pdf(med_hat))^2)
    pars <- c(med_hat, sqrt(sd_approx))
    if (non_par) {
      bounds <- sort(etas)[qbinom(p = c(1 - gamma, 1 + gamma)/2, size = length(etas), prob = 0.5)]
    } else {
      bounds <- qnorm(p = c(1 - gamma, 1 + gamma)/2, mean = pars[1], sd = pars[2])
    }
    len <- bounds[2] - bounds[1]
    if (l < len & M < max_iter) {
      M <- 2 * M
      if (M > max_iter) {
        aux <- F
      }
    } else {
      aux <- F
    }
  }
  if (M > max_iter) {
    warning(paste0("Maximum number of iterations reached. M = ", M))
  }
  return(list(M = M, med_hat = med_hat, bounds = bounds))
}
