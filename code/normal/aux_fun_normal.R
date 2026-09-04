# The normal-mean conjugate update (unknown mean, unknown variance, Normal-
# Inverse-Gamma prior) implemented below is mathematically the p = 1,
# intercept-only special case of the linear-regression conjugate update in
# code/linear_regression/aux_fun_lm.R: take X = a column of 1s, mu = [m],
# S = [v], and V = w * I_n. Rather than maintaining two independent (and
# previously slightly-drifting) implementations of the same algebra, the
# four functions below are thin wrappers that delegate to the lm versions
# and unwrap the resulting 1x1 matrices back into scalars, so there's a
# single implementation to keep correct.
#
# Trade-off: pred_par_conj_lm() builds a full m x m predictive covariance
# matrix, so pred_par_conj_normal() below is O(length(y_tilde)^2) rather
# than the previous O(1) -- fine for the sample sizes used in this repo
# (tens to low hundreds of test points), but worth knowing if this is ever
# reused with a much larger test set.
source("code/linear_regression/aux_fun_lm.R")

post_par_conj_normal <- function(m, v, a, b, w, y) {
  n <- length(y)
  res <- post_par_conj_lm(
    mu = matrix(m, 1, 1), S = matrix(v, 1, 1), a = a, b = b,
    V = w * diag(n), X = matrix(1, n, 1), y = y
  )
  list(
    m_star = as.numeric(res$mu_star),
    v_star = as.numeric(res$S_star),
    a_star = res$a_star,
    b_star = res$b_star
  )
}

pred_par_conj_normal <- function(m_star, v_star, a_star, b_star, w_tilde, y_tilde) {
  n <- length(y_tilde)
  res <- pred_par_conj_lm(
    mu_star = matrix(m_star, 1, 1), S_star = matrix(v_star, 1, 1),
    a_star = a_star, b_star = b_star,
    V_tilde = w_tilde * diag(n), X_tilde = matrix(1, n, 1), y_tilde = y_tilde
  )
  list(
    nu_pred = res$nu_pred,
    m_pred  = as.numeric(res$mu_pred[1]),
    v_pred  = as.numeric(res$S_pred[1, 1])
  )
}

pp_hyper_conj_normal <- function(eta, m, v, a, b, y0) {
  n0 <- length(y0)
  res <- pp_hyper_conj_lm(
    eta, mu = matrix(m, 1, 1), S = matrix(v, 1, 1), a = a, b = b,
    X0 = matrix(1, n0, 1), y0 = y0
  )
  list(
    m_star = as.numeric(res$mu_star),
    v_star = as.numeric(res$S_star),
    a_star = res$a_star,
    b_star = res$b_star
  )
}

pp_post_par_conj_normal <- function(eta, m, v, a, b, y0, w, y) {
  pp_prior_par <- pp_hyper_conj_normal(eta, m, v, a, b, y0)
  m <- pp_prior_par$m_star
  v <- pp_prior_par$v_star
  a <- pp_prior_par$a_star
  b <- pp_prior_par$b_star
  post_par <- post_par_conj_normal(m, v, a, b, w, y)
  return(post_par)
}

generate_normal_data <- function(n0 = 100, n = 100, mu0 = 1, mu = 1, sigma0 = 1, sigma = 1) {
  hist_data <- data.frame(y = rnorm(n0, mu0, sigma0), data = "hist")
  curr_data <- data.frame(y = rnorm(n, mu, sigma), data = "curr")
  rbind(hist_data, curr_data)
}

post_par_fixed_var <- function(m0, v0, v, y) {
  n <- length(y)
  v_star <- 1/(1/v0 + n/v)
  m_star <- v_star * (m0/v0 + sum(y)/v)
  return(list(m_star = m_star,
              v_star = v_star))
}

edelta_normal <- function(n0, v, post_par, hat_theta, nu = .5) {
  m_star <- post_par$m_star
  v_star <- post_par$v_star
  q2 <- (m_star - hat_theta)^2
  k0 <- n0^(nu)
  delta <- k0/v + k0^2/v^2 * (v_star + q2)
  return(delta)
}

estimate_eta <- function(data, post_par, v, alpha = 1/4, mle = F, nu = .5) {
  y0 <- data %>% filter(data == "hist") %>% pull(y)
  y <- data %>% filter(data == "curr") %>% pull(y)
  n0 <- length(y0)
  n <- length(y)

  theta0_mle <- mean(y0)
  theta_mle <- mean(y)
  diff <- abs(theta0_mle - theta_mle)
  if ((diff < n^(-1/2 + alpha)) && !mle) {
    hat_theta0 <- theta_mle
  } else {
    hat_theta0 <- theta0_mle
  }
  # Compute the delta for the current data
  delta_curr <- edelta_normal(n0, v, post_par, theta_mle, nu = nu)

  # Compute the delta for the historical data
  delta_hist <- edelta_normal(n0, v, post_par, hat_theta0, nu = nu)

  # Estimate eta using the ratio of deltas
  eta_estimate <- exp(0.5 * (log(delta_curr) - log(delta_hist)))

  return(eta_estimate)
}

estimate_theta <- function(data, alpha = 1/4) {
  y0 <- data %>% filter(data == "hist") %>% pull(y)
  y <- data %>% filter(data == "curr") %>% pull(y)
  n0 <- length(y0)
  n <- length(y)

  theta0_mle <- mean(y0)
  theta_mle <- mean(y)
  diff <- abs(theta0_mle - theta_mle)
  if (diff < n^(-1/2 + alpha)) {
    hat_theta0 <- theta_mle
  } else {
    hat_theta0 <- theta0_mle
  }
  return(list(
    hat_theta0 = hat_theta0,
    theta_mle = theta_mle,
    theta0_mle = theta0_mle
  ))
}

# compute 95% confidence intervals for each scenario and each n
med_np <- function(x, gamma = 0.95){
  ## método não-paramétrico, baseado na binomial
  med.hat <- median(x)
  return(
    c(med.hat, sort(x)[qbinom(p = c(1 - gamma, 1 + gamma)/2, size = length(x), prob = 0.5)])
  )
}
