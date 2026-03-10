library(mvtnorm)
library(parallel)
library(hdbayes)
library(cmdstanr)
library(dplyr)
library(posterior)

source("code/normal/aux_fun_normal.R")
load("data/sim_normal_data.RData")
y0 <- hist_data$y
n0 <- length(y0)
y_train <- curr_train$y
y_test <- curr_test$y
n_train <- length(y_train)
n_test <- length(y_test)

m <- 0
v <-  1
a <- 2 
b <- 1
a_tilde <- 1
b_tilde <- 1

#-------------------------------------------------------------------------------
# Sample NPP
#-------------------------------------------------------------------------------

# sample prior
sample_mu_pp <- function(eta) {
  pp_prior_par <- pp_hyper_conj_normal(eta, m, v, a, b, y0)
  m_star <- pp_prior_par$m_star
  v_star <- pp_prior_par$v_star
  a_star <- pp_prior_par$a_star
  b_star <- pp_prior_par$b_star
  tildemu <- LaplacesDemon::rst(1, mu = m_star, 
                                sigma = b_star / a_star * v_star, nu = 2*a_star)
  return(tildemu)
}

# sample post predictive
sample_postpred_pp <- function(eta) {
  pp_prior_par <- pp_hyper_conj_normal(eta, m, v, a, b, y0)
  m0 <- pp_prior_par$m_star
  v0 <- pp_prior_par$v_star
  a0 <- pp_prior_par$a_star
  b0 <- pp_prior_par$b_star
  post_par <- post_par_conj_normal(m0, v0, a0, b0, 1, y_train)
  m_star <- post_par$m_star
  v_star <- post_par$v_star
  a_star <- post_par$a_star
  b_star <- post_par$b_star
  w_tilde <- 1
  y_tilde <- y_test
  pred_par <- pred_par_conj_normal(m_star, v_star, a_star, b_star, w_tilde, y_tilde)
  
  tildey <- LaplacesDemon::rst(length(y_tilde), mu = pred_par$m_pred, 
                                sigma = pred_par$v_pred, nu = pred_par$nu_pred)
  return(tildey)
}

M <- 2000
draws_eta_prior <- rbeta(M, a_tilde, b_tilde)

draws_mu_npp_prior <- lapply(
  draws_eta_prior, sample_mu_pp
)
draws_mu_npp_prior <- do.call(rbind, draws_mu_npp_prior)

# sample posterior 
normal_model <- cmdstan_model("code/normal/normal_npp_eta.stan")

stan_data <- list(
  n = n_train,
  y = y_train,
  n0 = n0,
  y0 = y0,
  a = a,
  b = b,
  m = m,
  v = v,
  tilde_a = a_tilde,
  tilde_b = b_tilde
)

fit <- normal_model$sample(
  data = stan_data,
  iter_warmup = 2000,
  iter_sampling = 2000,
  chains = 4,
  parallel_chains = 4
)

draws_eta_post <- (fit$draws(c("eta")) %>% as_draws_df())$eta

sample_mu_pp_post <- function(eta) {
  pp_prior_par <- pp_hyper_conj_normal(eta, m, v, a, b, y0)
  m <- pp_prior_par$m_star
  v <- pp_prior_par$v_star
  a <- pp_prior_par$a_star
  b <- pp_prior_par$b_star
  w <- 1
  post_par <- post_par_conj_normal(m, v, a, b, w, y_train)
  m_star <- post_par$m_star
  v_star <- post_par$v_star
  a_star <- post_par$a_star
  b_star <- post_par$b_star
  
  tildemu <- LaplacesDemon::rst(1, mu = m_star, 
                                sigma = b_star / a_star * v_star, nu = 2*a_star)
  return(tildemu)
}

draws_mu_npp_post <- lapply(
  draws_eta_post, sample_mu_pp_post
)
draws_mu_npp_post <- do.call(rbind, draws_mu_npp_post)

draws_postpred_npp <- lapply(
  draws_eta_post, sample_postpred_pp
)
draws_postpred_npp <- do.call(rbind, draws_postpred_npp)

save(draws_eta_prior, file = 'samples/draws_eta_prior_normal.RData')
save(draws_mu_npp_prior, file = 'samples/draws_mu_npp_prior.RData')
save(draws_postpred_npp, file = 'samples/draws_postpred_npp_normal.RData')
save(draws_eta_post, file = 'samples/draws_eta_post_normal.RData')
save(draws_mu_npp_post, file = 'samples/draws_mu_npp_post.RData')

