library(cmdstanr)
library(posterior)

SEED <- 20260819
set.seed(SEED)

load("data/sim_normal_data.RData")

model <- cmdstan_model("code/normal/normal.stan")

y0 <- hist_data$y
n0 <- length(y0)
n1 <- 10000
m0 <- 0
v0 <-  1
a0 <- 2
b0 <- 1

# best_a0_crps is not computed in this script -- it comes from
# code/normal/ppc_normal_data.R (the CRPS-optimal eta). Run that first in
# the same session, or set eta explicitly below.
stan_data <- list(n0 = n0,
    y0 = y0,
    n1 = n1,
    m0 = m0,
    v0 = v0,
    a0 = a0,
    b0 = b0,
    eta = best_a0_crps
  )

fit <- model$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  seed = SEED
)

gen_q <- model$generate_quantities(
  data = stan_data,
  fitted_params = fit
)



mean((gen_q %>% 
  as_draws_df())$`y1[1]`)
