post_par_conj_lm <- function(mu, S, a, b, V, X, y) {
  # V: nxn data covariance matrix
  # mu: px1 prior mean
  # S: pxp prior covariance matrix
  # a: scalar, shape parameter of the inverse-gamma distribution
  # b: scalar, scale parameter of the inverse-gamma distribution
  # X: nxp design matrix
  # y: nx1 response vector
  
  n <- nrow(V)
  p <- length(mu)
  S_inv <- solve(S)
  V_inv <- solve(V)
  S_star <- solve(S_inv + t(X) %*% V_inv %*% X)
  mu_star <- S_star %*% (S_inv %*% mu + t(X) %*% V_inv %*% y)
  a_star <- a + n/2
  b_star <- b + 0.5 * (t(y) %*% V_inv %*% y + t(mu) %*% S_inv %*% mu - 
                         t(mu_star) %*% solve(S_star) %*% mu_star)
  return(list(mu_star = mu_star,
              S_star = S_star,
              a_star = a_star,
              b_star = as.numeric(b_star)))
}

pred_par_conj_lm <- function(mu_star, S_star, a_star, b_star, V_tilde, X_tilde, y_tilde) {
  # V_tilde: mxm covariance matrix for the new data
  # mu_star: px1 posterior mean from the previous function
  # S_star: pxp posterior covariance from the previous function
  # a_star: scalar, updated shape parameter of the inverse-gamma distribution
  # b_star: scalar, updated scale parameter of the inverse-gamma distribution
  # X_tilde: mxp design matrix for the new data
  # y_tilde: mx1 response vector for the new data
  
  m <- nrow(V_tilde)
  S_pred <- b_star / a_star * (solve(V_tilde) + X_tilde %*% S_star %*% t(X_tilde))
  mu_pred <- X_tilde %*% mu_star
  nu_pred <- 2 * a_star
  
  return(list(nu_pred = nu_pred,
              mu_pred = mu_pred,
              S_pred = S_pred))
}

pp_hyper_conj_lm <- function(eta, mu, S, a, b, X0, y0) {
  # eta: scalar, power parameter
  # mu: px1 prior mean
  # S: pxp prior covariance matrix
  # a: scalar, shape parameter of the inverse-gamma distribution
  # b: scalar, scale parameter of the inverse-gamma distribution
  # X0: nxp design matrix for the historical data
  # y0: nx1 response vector for the historical data
  
  n0 <- nrow(X0)
  if (eta == 0) {
    return(list(mu_star = mu,
                S_star = S,
                a_star = a,
                b_star = b))
  } else {
    V0 <- 1/eta * diag(n0)
    a0 <- a + 0.5*n0*(eta-1)
    post_pars <- post_par_conj_lm(mu, S, a0, b, V0, X0, y0)
    return(post_pars)
  }
}

pp_post_par_conj_lm <- function(eta, mu, S, a, b, X0, y0, V, X, y) {
  # eta: scalar, power parameter
  # mu: px1 prior mean
  # S: pxp prior covariance matrix
  # a: scalar, shape parameter of the inverse-gamma distribution
  # b: scalar, scale parameter of the inverse-gamma distribution
  # X0: nxp design matrix for the historical data
  # y0: nx1 response vector for the historical data
  # V: nxn covariance matrix for the current data
  # X: nxp design matrix for the current data
  # y: nx1 response vector for the current data
  
  pp_prior_par <- pp_hyper_conj_lm(eta, mu, S, a, b, X0, y0)
  mu <- pp_prior_par$mu_star
  S <- pp_prior_par$S_star
  a <- pp_prior_par$a_star
  b <- pp_prior_par$b_star
  n <- nrow(X)
  V <- diag(n)
  post_par <- post_par_conj_lm(mu, S, a, b, V, X, y)
  return(post_par)
}
