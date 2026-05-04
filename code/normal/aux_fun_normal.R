post_par_conj_normal <- function(m, v, a, b, w, y) {
  n <- length(y)
  v_star <- 1/(1/w*n + 1/v)
  m_star <- v_star * (m/v + 1/w* sum(y))
  a_star <- a + n/2
  b_star <- b + 0.5 * (1/w*t(y) %*% y + 1/v * m^2 - 
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
  n0 <- length(y0)
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

edelta_normal <- function(n0, v, post_par, hat_theta) {
  m_star <- post_par$m_star
  v_star <- post_par$v_star
  q2 <- (m_star - hat_theta)^2
  delta <- n0/v + n0^2/v^2 * (v_star + q2)
  return(delta)
}

estimate_eta <- function(data, post_par, v, alpha = 1/4, mle = F) {
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
  delta_curr <- edelta_normal(n0, v, post_par, theta_mle)
  
  # Compute the delta for the historical data
  delta_hist <- edelta_normal(n0, v, post_par, hat_theta0)
  
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
