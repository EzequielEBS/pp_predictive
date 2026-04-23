library(tidyverse)
library(dplyr)
library(ggplot2)
library(scoringRules)
library(LaplacesDemon)
library(pbapply)
library(parallel)

#-------------------------------------------------------------------------------
# auxiliary functions
#-------------------------------------------------------------------------------

source("code/normal/generate_normal_data.R")
source("code/normal/aux_fun_normal.R")


# est_eta_sce <- function(data) {
#   y <- data %>% filter(data == "curr") %>% pull(y)
#   post_par <- post_par_conj_normal(m, v, a, b, w, y)
#   # d <- post_draws(n_draws = n_draws, post_par = post_par)
#   estimate_eta(data)
# }

#-------------------------------------------------------------------------------
# analysis
#-------------------------------------------------------------------------------

n_list <- c(10, 50, 100, 200, 500, 1000, 10000)
r <- 1.5
n0_list <- r*n_list
n_rep <- 1000

results_cong_sce <- pblapply(1:n_rep, function(j) {
  print(paste("Replicate", j))
  res <- lapply(seq_along(n_list), function(i) {
    n <- n_list[i]
    n0 <- n0_list[i]
    m0 <- 0
    v0 <- 1
    mu0 <- 1
    mu1 <- 1
    v <- 1
    data <- generate_normal_data(n0 = n0, n = n, mu0 = mu0, mu = mu1, sigma0 = v, sigma = v)
    y0 <- data %>% filter(data == "hist") %>% pull(y)
    y <- data %>% filter(data == "curr") %>% pull(y)
    post_par <- post_par_fixed_var(m0, v0, v, y)
    eta_estimate <- estimate_eta(data, post_par, v)
    return(eta_estimate)
  })
  return(unlist(res))
}) %>% 
  do.call(rbind, .) %>% as.data.frame()
colnames(results_cong_sce) <- paste0("n_", n_list)

results_scong_sce <- lapply(1:n_rep, function(j) {
  print(paste("Replicate", j))
  res <- lapply(seq_along(n_list), function(i) {
    n <- n_list[i]
    n0 <- n0_list[i]
    m0 <- 0
    v0 <- 1
    mu0 <- 0
    mu1 <- 1
    v <- 1
    data <- generate_normal_data(n0 = n0, n = n, mu0 = mu0, mu = mu1, sigma0 = v, sigma = v)
    y0 <- data %>% filter(data == "hist") %>% pull(y)
    y <- data %>% filter(data == "curr") %>% pull(y)
    post_par <- post_par_fixed_var(m0, v0, v, y)
    eta_estimate <- estimate_eta(data, post_par, v)
    return(eta_estimate)
  })
  return(unlist(res))
}) %>% 
  do.call(rbind, .) %>% as.data.frame()
colnames(results_scong_sce) <- paste0("n_", n_list)

results_incong_sce <- lapply(1:n_rep, function(j) {
    print(paste("Replicate", j))
  res <- lapply(seq_along(n_list), function(i) {
    n <- n_list[i]
    n0 <- n0_list[i]
    m0 <- 0
    v0 <- 1
    mu0 <- -1
    mu1 <- 1
    v <- 1
    data <- generate_normal_data(n0 = n0, n = n, mu0 = mu0, mu = mu1, sigma0 = v, sigma = v)
    y0 <- data %>% filter(data == "hist") %>% pull(y)
    y <- data %>% filter(data == "curr") %>% pull(y)
    post_par <- post_par_fixed_var(m0, v0, v, y)
    eta_estimate <- estimate_eta(data, post_par, v)
    return(eta_estimate)
  })
  return(unlist(res))
}) %>% 
  do.call(rbind, .) %>% as.data.frame()
colnames(results_incong_sce) <- paste0("n_", n_list)

# compute median for each scenario and each n
results_cong_sce_median <- apply(results_cong_sce, 2, median)
results_scong_sce_median <- apply(results_scong_sce, 2, median)
results_incong_sce_median <- apply(results_incong_sce, 2, median)

# compute 95% confidence intervals for each scenario and each n
med_np <- function(x, gamma = 0.95){
  ## método não-paramétrico, baseado na binomial
  med.hat <- median(x)
  return(
    c(med.hat, sort(x)[qbinom(p = c(1 - gamma, 1 + gamma)/2, size = length(x), prob = 0.5)])
  )
}
results_cong_sce_ci <- apply(results_cong_sce, 2, med_np)
results_scong_sce_ci <- apply(results_scong_sce, 2, med_np)
results_incong_sce_ci <- apply(results_incong_sce, 2, med_np)

# Plotting confidence intervals
plot_data <- data.frame(
  n = rep(n_list, 3),
  scenario = rep(c("No discrepancy", "Small discrepancy", "Large discrepancy"), each = length(n_list)),
  median = c(results_cong_sce_median, results_scong_sce_median, results_incong_sce_median),
  lower_ci = c(results_cong_sce_ci[2, ], results_scong_sce_ci[2, ], results_incong_sce_ci[2, ]),
  upper_ci = c(results_cong_sce_ci[3, ], results_scong_sce_ci[3, ], results_incong_sce_ci[3, ])
)
plot_data$scenario <- factor(
  plot_data$scenario,
  levels = c("No discrepancy", "Small discrepancy", "Large discrepancy")
)

ggplot(plot_data, aes(x = as.factor(n), y = median)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2) +
  facet_wrap(~ scenario, scales = "free_y", ncol = 1) +
  labs(x = "Sample Size (n)", y = expression(hat(eta))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.text = element_text(size = 12))

ggsave("figures/inf_match_conv.png", width = 10, height = 4)
