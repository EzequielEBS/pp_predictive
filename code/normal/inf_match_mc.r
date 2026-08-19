library(purrr)
library(dplyr)

source("code/aux_fun_inf_match_glm.r")
source("code/normal/aux_fun_normal.R")

#-------------------------------------------------------------------------------
# No discrepancy
#-------------------------------------------------------------------------------

n <- 100
r <- 1.5
n0 <- r*n
mu0 <- 0
mu1 <- 1
v <- 1
data <- generate_normal_data(n0 = n0, n = n, mu0 = mu0, mu = mu1, sigma0 = v, sigma = v)
hist_data <- data %>% filter(data == "hist")
curr_data <- data %>% filter(data == "curr")
formula <- as.formula("y ~ 1")
y0 <- hist_data %>% pull(y)
y <- curr_data %>% pull(y)
X <- matrix(1, nrow = length(y), ncol = 1)
X0 <- matrix(1, nrow = length(y0), ncol = 1)
family <- gaussian()

m0 <- 0
v0 <- 1
post_par <- post_par_fixed_var(m0, v0, v, y)
beta_draws <- matrix(
  rnorm(10000, mean = post_par$m_star, sd = sqrt(post_par$v_star)),
  ncol = 1
)

hat_eta_glm <- estimate_eta_glm(
  formula = formula,
  curr_data = curr_data,
  hist_data = hist_data,
  beta_draws = beta_draws, 
  phi = v, 
  family = family
)
hat_eta_normal <- estimate_eta(data, post_par, v, mle = T)
print(hat_eta_glm)
print(hat_eta_normal)
