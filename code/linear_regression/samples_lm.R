library(mvtnorm)
library(parallel)
library(hdbayes)
library(cmdstanr)
library(dplyr)
library(posterior)

source("code/linear_regression/aux_fun_lm.R")
load("data/sim_lm_data.RData")
formula <- y ~ X1
family <- gaussian()

res_hist <- hdbayes:::stack.data(formula = formula, data.list = list(hist_data))
res_curr <- hdbayes:::stack.data(formula = formula, data.list = list(curr_data))
y0 <- res_hist$y
X0 <- res_hist$X
n0 <- length(y0)
y <- res_curr$y
X <- res_curr$X
n <- length(y)
p <- ncol(X)


mu <- rep(0, p)
S <-  diag(p)
a <- 2 
b <- 1
a_tilde <- 1
b_tilde <- 1

#-------------------------------------------------------------------------------
# Sample NPP
#-------------------------------------------------------------------------------

# sample prior
sample_beta_pp <- function(eta) {
  pp_prior_par <- pp_hyper_conj_lm(eta, mu, S, a, b, X0, y0)
  mu_star <- pp_prior_par$mu_star
  S_star <- pp_prior_par$S_star
  a_star <- pp_prior_par$a_star
  b_star <- pp_prior_par$b_star
  V_tilde <- diag(n)
  
  tildebeta <- mvtnorm::rmvt(1, delta = mu_star, 
                             sigma = b_star / a_star * S_star,
                             df = 2*a_star)
  return(tildebeta)
}

# sample prior predictive
sample_priorpred_pp <- function(eta) {
  pp_prior_par <- pp_hyper_conj_lm(eta, mu, S, a, b, X0, y0)
  mu_star <- pp_prior_par$mu_star
  S_star <- pp_prior_par$S_star
  a_star <- pp_prior_par$a_star
  b_star <- pp_prior_par$b_star
  V_tilde <- diag(n)
  X_tilde <- X
  y_tilde <- y
  pred_par <- pred_par_conj_lm(mu_star, S_star, a_star, b_star, V_tilde, X_tilde, y_tilde)
  
  tildey <- mvtnorm::rmvt(1, delta = pred_par$mu_pred, sigma = pred_par$S_pred, df = pred_par$nu_pred)
  return(tildey)
}

m <- 2000
draws_eta_prior <- rbeta(m, a_tilde, b_tilde)

draws_beta_npp_prior <- lapply(
  draws_eta_prior, sample_beta_pp
)
draws_priorpred_npp <- lapply(
  draws_eta_prior, sample_priorpred_pp
)
draws_beta_npp_prior <- do.call(rbind, draws_beta_npp_prior)
draws_priorpred_npp <- do.call(rbind, draws_priorpred_npp)


# sample posterior 
lm_model <- cmdstan_model("code/linear_regression/lm_npp_eta.stan")

stan_data <- list(
  n = n,
  p = p,
  X = X,
  y = y,
  n0 = n0,
  X0 = X0,
  y0 = y0,
  a = a,
  b = b,
  mu = mu,
  V = S,
  tilde_a = a_tilde,
  tilde_b = b_tilde
)

fit <- lm_model$sample(
  data = stan_data,
  iter_warmup = 2000,
  iter_sampling = 2000,
  chains = 4,
  parallel_chains = 4
)

draws_eta_post <- (fit$draws(c("eta")) %>% as_draws_df())$eta

sample_beta_pp_post <- function(eta) {
  pp_prior_par <- pp_hyper_conj_lm(eta, mu, S, a, b, X0, y0)
  mu <- pp_prior_par$mu_star
  S <- pp_prior_par$S_star
  a <- pp_prior_par$a_star
  b <- pp_prior_par$b_star
  V <- diag(n)
  post_par <- post_par_conj_lm(mu, S, a, b, V, X, y)
  mu_star <- post_par$mu_star
  S_star <- post_par$S_star
  a_star <- post_par$a_star
  b_star <- post_par$b_star
  
  tildebeta <- mvtnorm::rmvt(1, delta = mu_star, 
                             sigma = b_star / a_star * S_star,
                             df = 2*a_star)
  return(tildebeta)
}

draws_beta_npp_post <- lapply(
  draws_eta_post, sample_beta_pp_post
)
draws_beta_npp_post <- do.call(rbind, draws_beta_npp_post)

save(draws_eta_prior, file = 'samples/draws_eta_prior.RData')
save(draws_beta_npp_prior, file = 'samples/draws_beta_npp_prior.RData')
save(draws_priorpred_npp, file = 'samples/draws_priorpred_npp.RData')
save(draws_eta_post, file = 'samples/draws_eta_post.RData')
save(draws_beta_npp_post, file = 'samples/draws_beta_npp_post.RData')

