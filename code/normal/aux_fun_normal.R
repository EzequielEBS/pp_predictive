post_par_conj_normal <- function(m, v, a, b, w, y) {
  n <- length(y)
  v_star <- 1/(1/w*n + v)
  m_star <- v_star * (m/v + 1/w* sum(y))
  a_star <- a + n/2
  b_star <- b + 0.5 * (t(y) %*% y + 1/v * m^2 - 
                         1/v_star * m_star^2)
  return(list(m_star = m_star,
              v_star = v_star,
              a_star = a_star,
              b_star = as.numeric(b_star)))
}

pred_par_conj_normal <- function(m_star, v_star, a_star, b_star, w_tilde, y_tilde) {
  n <- length(y_tilde)
  v_pred <- b_star / a_star * (1/w_tilde + v_star)
  m_pred <- m_star
  nu_pred <- 2 * a_star
  
  return(list(nu_pred = nu_pred,
              m_pred = m_pred,
              v_pred = v_pred))
}

pp_hyper_conj_normal <- function(eta, m, v, a, b, y0) {
  n0 <- nrow(X0)
  if (eta == 0) {
    return(list(m_star = m,
                v_star = v,
                a_star = a,
                b_star = b))
  } else {
    w0 <- 1/eta
    a0 <- a + 0.5*n0*(eta-1)
    post_pars <- post_par_conj_normal(m, v, a0, b, w0, y0)
    return(post_pars)
  }
}

pp_post_par_conj_normal <- function(eta, m, v, a, b, y0, w, y) {
  pp_prior_par <- pp_hyper_conj_normal(eta, m, v, a, b, y0)
  m <- pp_prior_par$m_star
  v <- pp_prior_par$v_star
  a <- pp_prior_par$a_star
  b <- pp_prior_par$b_star
  w <- 1
  post_par <- post_par_conj_normal(m, v, a, b, w, y)
  return(post_par)
}