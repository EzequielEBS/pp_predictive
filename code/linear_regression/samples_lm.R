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
res_curr_train <- hdbayes:::stack.data(formula = formula, data.list = list(curr_train))
res_curr_test <- hdbayes:::stack.data(formula = formula, data.list = list(curr_test))
y0 <- res_hist$y
X0 <- res_hist$X
n0 <- length(y0)
y_train <- res_curr_train$y
X_train <- res_curr_train$X
y_test <- res_curr_test$y
X_test <- res_curr_test$X
n_train <- length(y_train)
n_test <- length(y_test)
p <- ncol(X0)


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

# sample post predictive
sample_postpred_pp <- function(eta) {
  pp_prior_par <- pp_hyper_conj_lm(eta, mu, S, a, b, X0, y0)
  mu0 <- pp_prior_par$mu_star
  S0  <- pp_prior_par$S_star
  a0  <- pp_prior_par$a_star
  b0  <- pp_prior_par$b_star
  post_par <- post_par_conj_lm(mu0, S0, a0, b0, diag(n_train), X_train, y_train)
  mu_star <- post_par$mu_star
  S_star  <- post_par$S_star
  a_star  <- post_par$a_star
  b_star  <- post_par$b_star
  V_tilde <- diag(n_test)
  X_tilde <- X_test
  y_tilde <- y_test
  pred_par <- pred_par_conj_lm(mu_star, S_star, a_star, b_star, V_tilde, X_tilde, y_tilde)
  
  tildey <- mvtnorm::rmvt(1, delta = pred_par$mu_pred, sigma = pred_par$S_pred, df = pred_par$nu_pred)
  return(tildey)
}

m <- 2000
draws_eta_prior <- rbeta(m, a_tilde, b_tilde)

draws_beta_npp_prior <- lapply(
  draws_eta_prior, sample_beta_pp
)
draws_beta_npp_prior <- do.call(rbind, draws_beta_npp_prior)


# sample posterior 
lm_model <- cmdstan_model("code/linear_regression/lm_npp_eta.stan")

stan_data <- list(
  n = n_train,
  p = p,
  X = X_train,
  y = y_train,
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
  V <- diag(n_train)
  post_par <- post_par_conj_lm(mu, S, a, b, V, X_train, y_train)
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

draws_postpred_npp <- lapply(
  draws_eta_post, sample_postpred_pp
)
draws_postpred_npp <- do.call(rbind, draws_postpred_npp)

save(draws_eta_prior, file = 'samples/draws_eta_prior.RData')
save(draws_beta_npp_prior, file = 'samples/draws_beta_npp_prior.RData')
save(draws_postpred_npp, file = 'samples/draws_postpred_npp.RData')
save(draws_eta_post, file = 'samples/draws_eta_post.RData')
save(draws_beta_npp_post, file = 'samples/draws_beta_npp_post.RData')

