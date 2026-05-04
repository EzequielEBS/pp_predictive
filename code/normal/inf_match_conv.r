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

source("code/normal/aux_fun_normal.R")

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


theta_cong_sce <- pblapply(1:n_rep, function(j) {
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
    theta_estimates <- estimate_theta(data, alpha = 1/4)
    df <- data.frame(
      n = n, 
      rep = j, 
      hat_theta0 = theta_estimates$hat_theta0,
      theta0_mle = theta_estimates$theta0_mle,
      theta_mle = theta_estimates$theta_mle
    )
    return(df)
  })
  res <- do.call(rbind, res)
  return(res)
}) %>% 
  do.call(rbind, .) %>% as.data.frame()

theta_scong_sce <- pblapply(1:n_rep, function(j) {
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
    theta_estimates <- estimate_theta(data, alpha = 1/4)
    df <- data.frame(
      n = n, 
      rep = j, 
      hat_theta0 = theta_estimates$hat_theta0,
      theta0_mle = theta_estimates$theta0_mle,
      theta_mle = theta_estimates$theta_mle
    )
    return(df)
  })
  res <- do.call(rbind, res)
  return(res)
}) %>% 
  do.call(rbind, .) %>% as.data.frame()

theta_incong_sce <- pblapply(1:n_rep, function(j) {
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
    theta_estimates <- estimate_theta(data, alpha = 1/4)
    df <- data.frame(
      n = n, 
      rep = j, 
      hat_theta0 = theta_estimates$hat_theta0,
      theta0_mle = theta_estimates$theta0_mle,
      theta_mle = theta_estimates$theta_mle
    )
    return(df)
  }) %>% 
    do.call(rbind, .) %>% as.data.frame()
  return(res)
}) %>% 
  do.call(rbind, .) %>% as.data.frame()

plot_data_theta <- data.frame(
  scenario = 
    rep(c("No discrepancy", "Small discrepancy", "Large discrepancy"), each = nrow(theta_cong_sce))
)
plot_data_theta <- rbind(
  cbind(theta_cong_sce, scenario = "No discrepancy"),
  cbind(theta_scong_sce, scenario = "Small discrepancy"),
  cbind(theta_incong_sce, scenario = "Large discrepancy")
)

plot_data_theta <- plot_data_theta %>%
  group_by(scenario, n) %>%
  summarise(
    diff_hat_theta0 = median(abs(hat_theta0 - theta_mle)),
    diff_theta0_mle = median(abs(theta0_mle - theta_mle))
  ) %>% ungroup()

# plot_data_theta$scenario <- factor(
#   plot_data$scenario,
#   levels = c("No discrepancy", "Small discrepancy", "Large discrepancy")
# )
plot_data_theta <- plot_data_theta %>%
  pivot_longer(
    cols = c("diff_hat_theta0", "diff_theta0_mle"),
    names_to = "estimate_type",
    values_to = "estimate_value"
  )

plot_data_theta$scenario <- factor(
  plot_data_theta$scenario,
  levels = c("No discrepancy", "Small discrepancy", "Large discrepancy")
)
n_df <- data.frame(
  scenario = rep(c("No discrepancy", "Small discrepancy", "Large discrepancy"), each = length(n_list)),
  n = rep(n_list, 3),
  n_cutoff = n_list^(-1/4)
)
n_df$scenario <- factor(
  n_df$scenario,
  levels = c("No discrepancy", "Small discrepancy", "Large discrepancy")
)

ggplot(plot_data_theta,
  aes(x = as.factor(n), y = estimate_value, color = estimate_type)) +
  geom_point() +
  # geom_point(data = n_df, aes(x = as.factor(n), y = n_cutoff), color = "black", shape = 4, size = 3) +
  scale_color_manual(values = c("blue", "red"),
                      labels = c(bquote("||"~tilde(theta)[0] - hat(theta)~"||"), bquote("||"~hat(theta)[0] - hat(theta)~"||")),
                      name = "Median") +
  facet_wrap(~ scenario, scales = "free_y", ncol = 1) +
  labs(x = "Sample Size (n)", y = "") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.text = element_text(size = 12))
