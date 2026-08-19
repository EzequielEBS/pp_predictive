library(hdbayes)
library(posterior)
library(dplyr)
library(parallel)
library(purrr)
library(arrow)
library(subsampling)

source("code/aux_fun_inf_match_glm.r")


## obtain number of cores
ncores        = 1
chains        = 4      ## number of Markov chains to run
iter_warmup   = 5000   ## warmup per chain for MCMC sampling
iter_sampling = 10000   ## number of samples post warmup per chain
formula = outcome ~ age + treatment + race + cd4
family  = binomial('logit')

hist_data <- actg019
curr_data <- actg036

# normalize age and cd4
curr_data$age <- (curr_data$age - mean(curr_data$age)) /
  (sd(curr_data$age))
curr_data$cd4 <- (curr_data$cd4 - mean(curr_data$cd4)) /
  (sd(curr_data$cd4))
hist_data$age <- (hist_data$age - mean(hist_data$age)) /
  (sd(hist_data$age))
hist_data$cd4 <- (hist_data$cd4 - mean(hist_data$cd4)) /
  (sd(hist_data$cd4))

data.list = list(curr_data, hist_data)

eta_inter <- 1/2
eta_rate <- 1/2 * nrow(curr_data) / nrow(hist_data)
eta_small <- 0.1
eta_big <- 0.9


fit0 = glm.pp(
  formula = formula, family = family, data.list = data.list,
  a0.vals = 0,
  iter_warmup = iter_warmup, iter_sampling = iter_sampling, 
  chains = chains, parallel_chains = ncores,
  refresh = 0
)
beta_draws <- fit0[, -1] %>% as_draws_matrix()

l <- 0.01
gamma <- 0.95
best_sample_size <- compute_n(
  l, 
  gamma,
  formula,
  curr_data,
  hist_data,
  beta_draws,
  family,
  max_iter = 2000,
  starting_M = 100,
  non_par = T
)
best_sample_size

ncores <- detectCores() - 1
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
etas <- parLapply(cl, 1:best_sample_size$M, function(i) {
  estimate_eta_glm(
    formula = formula,
    curr_data = curr_data,
    hist_data = hist_data,
    beta_draws = beta_draws, 
    family = family
  )
}) %>% unlist()
stopCluster(cl)

eta_inf_match <- median(etas)
eta_inf_match

etas_df <- data.frame(eta = etas)


fit_eta_inter <- glm.pp(
  formula = formula, family = family, data.list = data.list,
  a0.vals = eta_inter,
  iter_warmup = iter_warmup, iter_sampling = iter_sampling, 
  chains = chains, parallel_chains = ncores,
  refresh = 0
)
fit_eta_rate <- glm.pp(
  formula = formula, family = family, data.list = data.list,
  a0.vals = eta_rate,
  iter_warmup = iter_warmup, iter_sampling = iter_sampling, 
  chains = chains, parallel_chains = ncores,
  refresh = 0
)
fit_eta_small <- glm.pp(
  formula = formula, family = family, data.list = data.list,
  a0.vals = eta_small,
  iter_warmup = iter_warmup, iter_sampling = iter_sampling, 
  chains = chains, parallel_chains = ncores,
  refresh = 0
)
fit_eta_big <- glm.pp(
  formula = formula, family = family, data.list = data.list,
  a0.vals = eta_big,
  iter_warmup = iter_warmup, iter_sampling = iter_sampling, 
  chains = chains, parallel_chains = ncores,
  refresh = 0
)
fit_eta_inf_match <- glm.pp(
  formula = formula, family = family, data.list = data.list,
  a0.vals = eta_inf_match,
  iter_warmup = iter_warmup, iter_sampling = iter_sampling, 
  chains = chains, parallel_chains = ncores,
  refresh = 0
)


a0     = seq(0, 1, length.out = 21)

## wrapper to obtain log normalizing constant in parallel package
logncfun = function(a0, ...){
  hdbayes::glm.npp.lognc(
    formula = formula, family = family, histdata = hist_data, a0 = a0, ...
  )
}

cl = makeCluster(15)
clusterSetRNGStream(cl, 123)
clusterExport(cl, varlist = c('formula', 'family', 'hist_data'))

## call created function
a0.lognc = parLapply(
  cl = cl, X = a0, fun = logncfun, iter_warmup = iter_warmup, iter_sampling = 5000, 
  chains = chains, refresh = 0
)
stopCluster(cl)

a0.lognc = data.frame( do.call(rbind, a0.lognc) )

fit_npp = glm.npp(
  formula = formula, family = family, data.list = data.list,
  a0.lognc = a0.lognc$a0,
  lognc = matrix(a0.lognc$lognc, ncol = 1),
  iter_warmup = iter_warmup, iter_sampling = iter_sampling, 
  chains = chains, parallel_chains = ncores,
  refresh = 0
)

fit_vnpp <- read_parquet("samples/logistic_regression/samples_variational.parquet")

fit_eta_inter <- fit_eta_inter %>% select(-.chain, -.iteration, -.draw, -lp__)
fit_eta_inter$eta <- eta_inter
fit_eta_inter$method <- "inter"
fit_eta_rate <- fit_eta_rate %>% select(-.chain, -.iteration, -.draw, -lp__)
fit_eta_rate$eta <- eta_rate
fit_eta_rate$method <- "rate"
fit_eta_small <- fit_eta_small %>% select(-.chain, -.iteration, -.draw, -lp__)
fit_eta_small$eta <- eta_small
fit_eta_small$method <- "small"
fit_eta_big <- fit_eta_big %>% select(-.chain, -.iteration, -.draw, -lp__)
fit_eta_big$eta <- eta_big
fit_eta_big$method <- "big"
fit_eta_inf_match <- fit_eta_inf_match %>% select(-.chain, -.iteration, -.draw, -lp__)
fit_eta_inf_match$eta <- eta_inf_match
fit_eta_inf_match$method <- "inf_match"
fit_npp <- fit_npp %>% select(-.chain, -.iteration, -.draw, -lp__, -`logit_a0s[1]`)
fit_npp$method <- "npp"
fit_npp <- fit_npp %>% rename(eta = a0_hist_1)
fit_vnpp$method <- "vnpp"

all_fits <- bind_rows(
  fit_eta_inter, fit_eta_rate, fit_eta_small, fit_eta_big, fit_eta_inf_match,
  fit_npp, fit_vnpp
)
write_parquet(all_fits, "samples/logistic_regression/samples_all.parquet")
write_parquet(etas_df, "samples/logistic_regression/etas_inf_match.parquet")
