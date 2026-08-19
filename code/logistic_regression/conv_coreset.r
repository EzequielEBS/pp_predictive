library(tidyverse)
library(hdbayes)
library(posterior)
library(subsampling)

set.seed(06032026)

source("code/aux_fun_inf_match_glm.r")

generate_log_reg_data <- function(n, r, p, c, beta) {
  n0 <- floor(n * r)
  X0 <- cbind(1, matrix(rnorm(n0 * p), nrow = n0))
  y0 <- rbinom(n0, size = 1, prob = plogis(c*X0 %*% beta))
  
  X <- cbind(1, matrix(rnorm(n * p), nrow = n))
  y <- rbinom(n, size = 1, prob = plogis(X %*% beta))

  df_hist <- data.frame(X0, y = y0, data_id = "hist")
  df_curr <- data.frame(X, y = y, data_id = "current")
  return(rbind(df_hist, df_curr))
}

#------------------------------------------
# Small discrepancy setting
#------------------------------------------

n <- 100
p <- 3
r <- 1.5
c <- 2
beta <- matrix(c(0.5, 1, -1, 0.8), ncol = 1)

n_list <- c(10, 50, 100, 500, 1000, 5000, 10000, 50000, 100000)

formula <- y ~ X1 + X2 + X3 + X4
family <- binomial()

etas_no_disc <- lapply(n_list, function(n) {
  data <- generate_log_reg_data(n = n, r = r, p = p, c = c, beta = beta)
  hist_data <- data %>% filter(data_id == "hist")
  curr_data <- data %>% filter(data_id == "current")
  fit0 <- glm.pp(
    formula = formula,
    family = family,
    data.list = list(curr_data, hist_data),
    a0.vals = 0,
    iter_warmup = 5000,
    iter_sampling = 10000,
    chains = 4,
    parallel_chains = 4,
    refresh = 0
  )
  beta_draws <- fit0[, -1] %>% as_draws_matrix()
  estimate_eta_glm(
    formula = formula,
    curr_data = curr_data,
    hist_data = hist_data,
    beta_draws = beta_draws,
    family = family,
    ssp = TRUE
  )
})
