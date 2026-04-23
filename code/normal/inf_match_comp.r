library(tidyverse)
library(dplyr)
library(ggplot2)
library(scoringRules)
library(LaplacesDemon)
library(patchwork)

#-------------------------------------------------------------------------------
# auxiliary functions
#-------------------------------------------------------------------------------

source("code/normal/generate_normal_data.R")
source("code/normal/aux_fun_normal.R")

comp_crps <- function(eta){
  single <- function(e){
    pp_prior_par <- post_par_fixed_var(m0, v0, v/e, y0)
    m0 <- pp_prior_par$m_star
    v0 <- pp_prior_par$v_star
    post_par <- post_par_fixed_var(m0, v0, v, y_train)
    m_star <- post_par$m_star
    v_star <- post_par$v_star
    y_tilde <- y_test
    pred_par <- list(m_pred = m_star,
                     v_pred = v_star + v)
    mean(crps_norm(y_tilde, pred_par$m_pred, sqrt(pred_par$v_pred)))
  }
  if (length(eta) > 1) sapply(eta, single) else single(eta)
}


#-------------------------------------------------------------------------------
# comparison
#-------------------------------------------------------------------------------

n <- 100
r <- 1
n0 <- r * n
m0 <- 0
v0 <- 1
v <- 1
mu0 <- 1
mu1 <- 1
data_cong <- generate_normal_data(n0 = n0, n = n, mu0 = mu0, mu = mu1, sigma0 = v, sigma = v)
y0_cong <- data_cong %>% filter(data == "hist") %>% pull(y)
y_cong <- data_cong %>% filter(data == "curr") %>% pull(y)
idx <- sample(seq_len(length(y_cong)), size = 0.3 * length(y_cong))
y_train <- y_cong[idx]
y_test <- y_cong[-idx]


y0 <- y0_cong
crps_optim <- optimize(comp_crps, interval = c(0,1))
best_a0_crps <- crps_optim$minimum
# hyva_optim <- optimize(comp_hyvarinen, interval = c(0,1))
# best_a0_hyvarinen <- hyva_optim$minimum

post_par_pool_cong <- post_par_fixed_var(m0, v0, v, data_cong$y)
post_par_cong <- post_par_fixed_var(m0, v0, v, y_cong)
eta_inf_cong <- estimate_eta(data_cong, post_par_cong, v)
pp_par_crps_cong <- if (best_a0_crps == 0) {
  list(m_star = m0, v_star = v0)
} else {
  post_par_fixed_var(m0, v0, v/best_a0_crps, y0_cong)
}
pp_par_inf_cong <- post_par_fixed_var(m0, v0, v/eta_inf_cong, y0_cong)
post_par_crps_cong <- post_par_fixed_var(
  pp_par_crps_cong$m_star, pp_par_crps_cong$v_star, v, y_cong
)
post_par_inf_cong <- post_par_fixed_var(
  pp_par_inf_cong$m_star, pp_par_inf_cong$v_star, v, y_cong
)
eta_true_cong <- sqrt(post_par_cong$v_star/v + 1/v * (post_par_cong$m_star - mu1)^2 + 1/n0) / 
  sqrt(post_par_cong$v_star/v + 1/v * (post_par_cong$m_star - mu0)^2 + 1/n0)

print(paste("True eta:", round(eta_true_cong, 3)))
print(paste("Estimated eta:", round(eta_inf_cong, 3)))
print(paste("CRPS eta:", round(best_a0_crps, 3)))

post_cong <- ggplot() +
  stat_function(fun = function(x) dnorm(x, mean = post_par_pool_cong$m_star, 
                                      sd = sqrt(post_par_pool_cong$v_star)),
                aes(color = "pool")) +
  stat_function(fun = function(x) dnorm(x, mean = post_par_cong$m_star, 
                                      sd = sqrt(post_par_cong$v_star)),
                aes(color = "curr")) +
  # stat_function(fun = function(x) dnorm(x, mean = post_par_crps_cong$m_star, 
  #                                     sd = sqrt(post_par_crps_cong$v_star)),
  #               aes(color = "crps")) +
  stat_function(fun = function(x) dnorm(x, mean = post_par_inf_cong$m_star, 
                                      sd = sqrt(post_par_inf_cong$v_star)),
                aes(color = "inf")) +
  geom_vline(aes(xintercept = mean(y0_cong), color = "y0", linetype = "y0")) +
  geom_vline(aes(xintercept = mean(y_cong),  color = "y",  linetype = "y")) +
  scale_color_manual(values = c(
                              "pool" = "#0079fbff",  # blue
                              "curr" = "#E45756",  # red
                              "crps" = "#72B7B2",  # teal
                              "inf"  = "#B279A2",  # muted purple
                              "y0"   = "#2E2E2E",  # dark gray (better than pure black)
                              "y"    = "#F58518"   # orange
                              ),
                      labels = c("pool" = "Pooled", 
                                  "curr" = "Curr", 
                                  "crps" = "CRPS", 
                                  "inf" = "Inf-Match",
                                  "y0" = "Hist Mean",
                                  "y" = "Curr Mean"
                              )
                            ) +
  scale_linetype_manual(values = c("y0" = "dashed", "y" = "dashed")) +
  labs(x = "", y = "", color = "") +
  xlim(1 - 1, 1 + 1) +
  theme_minimal() +
  guides(linetype = "none") +
  ggtitle("No discrepancy")
print(post_cong)

mu0 <- 0
mu1 <- 1
data_scong <- generate_normal_data(n0 = n0, n = n, mu0 = mu0, mu = mu1, sigma0 = v, sigma = v)
y0_scong <- data_scong %>% filter(data == "hist") %>% pull(y)
y_scong <- data_scong %>% filter(data == "curr") %>% pull(y)
idx <- sample(seq_len(length(y_scong)), size = 0.3 * length(y_scong))
y_train <- y_scong[idx]
y_test <- y_scong[-idx]

y0 <- y0_scong
crps_optim_scong <- optimize(comp_crps, interval = c(0,1))
best_a0_crps_scong <- crps_optim_scong$minimum

post_par_pool_scong <- post_par_fixed_var(m0, v0, v, data_scong$y)
post_par_scong <- post_par_fixed_var(m0, v0, v, y_scong)
eta_inf_scong <- estimate_eta(data_scong, post_par_scong, v)
pp_par_crps_scong <- if (best_a0_crps_scong == 0) {
  list(m_star = m0, v_star = v0)
} else {
  post_par_fixed_var(m0, v0, v/best_a0_crps_scong, y0_scong)
}
post_par_crps_scong <- post_par_fixed_var(
  pp_par_crps_scong$m_star, pp_par_crps_scong$v_star, v, y_scong
)
pp_par_inf_scong <- post_par_fixed_var(m0, v0, v/eta_inf_scong, y0_scong)
post_par_inf_scong <- post_par_fixed_var(
  pp_par_inf_scong$m_star, pp_par_inf_scong$v_star, v, y_scong
)
eta_true_scong <- sqrt(post_par_scong$v_star/v + 1/v * (post_par_scong$m_star - mu1)^2 + 1/n0) / 
  sqrt(post_par_scong$v_star/v + 1/v * (post_par_scong$m_star - mu0)^2 + 1/n0)

print(paste("True eta:", round(eta_true_scong, 3)))
print(paste("Estimated eta:", round(eta_inf_scong, 3)))
print(paste("CRPS eta:", round(best_a0_crps_scong, 3)))

post_scong <- ggplot() +
  stat_function(fun = function(x) dnorm(x, mean = post_par_pool_scong$m_star, 
                                      sd = sqrt(post_par_pool_scong$v_star)),
                aes(color = "pool")) +
  stat_function(fun = function(x) dnorm(x, mean = post_par_scong$m_star, 
                                      sd = sqrt(post_par_scong$v_star)),
                aes(color = "curr")) +
  # stat_function(fun = function(x) dnorm(x, mean = post_par_crps_scong$m_star, 
  #                                     sd = sqrt(post_par_crps_scong$v_star)),
  #               aes(color = "crps")) +
  stat_function(fun = function(x) dnorm(x, mean = post_par_inf_scong$m_star, 
                                      sd = sqrt(post_par_inf_scong$v_star)),
                aes(color = "inf")) +
  geom_vline(aes(xintercept = mean(y0_scong), color = "y0", linetype = "y0")) +
  geom_vline(aes(xintercept = mean(y_scong),  color = "y",  linetype = "y")) +
  scale_color_manual(values = c(
                              "pool" = "#0079fbff",  # blue
                              "curr" = "#E45756",  # red
                              "crps" = "#72B7B2",  # teal
                              "inf"  = "#B279A2",  # muted purple
                              "y0"   = "#2E2E2E",  # dark gray (better than pure black)
                              "y"    = "#F58518"   # orange
                              ),
                      labels = c("pool" = "Pooled", 
                                  "curr" = "Curr", 
                                  "crps" = "CRPS", 
                                  "inf" = "Inf-Match",
                                  "y0" = "Hist Mean",
                                  "y" = "Curr Mean"
                              )
                            ) +
  scale_linetype_manual(values = c("y0" = "dashed", "y" = "dashed")) +
  labs(x = "", y = "", color = "") +
  ggtitle("Small discrepancy") +
  xlim(0.5 - 1, 0.5 + 1) +
  theme_minimal() +
  guides(linetype = "none")

print(post_scong)

mu0 <- -1
mu1 <- 1
data_incong <- generate_normal_data(n0 = n0, n = n, mu0 = mu0, mu = mu1, sigma0 = v, sigma = v)
y0_incong <- data_incong %>% filter(data == "hist") %>% pull(y)
y_incong <- data_incong %>% filter(data == "curr") %>% pull(y)
idx <- sample(seq_len(length(y_incong)), size = 0.3 * length(y_incong))
y_train <- y_incong[idx]
y_test <- y_incong[-idx]

y0 <- y0_incong
crps_optim_incong <- optimize(comp_crps, interval = c(0,1))
best_a0_crps_incong <- crps_optim_incong$minimum

post_par_pool_incong <- post_par_fixed_var(m0, v0, v, data_incong$y)
post_par_incong <- post_par_fixed_var(m0, v0, v, y_incong)
eta_inf_incong <- estimate_eta(data_incong, post_par_incong, v)
pp_par_crps_incong <- if (best_a0_crps_incong == 0) {
  list(m_star = m0, v_star = v0)
} else {
  post_par_fixed_var(m0, v0, v/best_a0_crps_incong, y0_incong)
}
post_par_crps_incong <- post_par_fixed_var(
  pp_par_crps_incong$m_star, pp_par_crps_incong$v_star, v, y_incong
)
pp_par_inf_incong <- post_par_fixed_var(m0, v0, v/eta_inf_incong, y0_incong)
post_par_inf_incong <- post_par_fixed_var(
  pp_par_inf_incong$m_star, pp_par_inf_incong$v_star, v, y_incong
)
eta_true_incong <- sqrt(post_par_incong$v_star/v + 1/v * (post_par_incong$m_star - mu1)^2 + 1/n0) / 
  sqrt(post_par_incong$v_star/v + 1/v * (post_par_incong$m_star - mu0)^2 + 1/n0)

print(paste("True eta:", round(eta_true_incong, 3)))
print(paste("Estimated eta:", round(eta_inf_incong, 3)))
print(paste("CRPS eta:", round(best_a0_crps_incong, 3)))

post_incong <- ggplot() +
  stat_function(fun = function(x) dnorm(x, mean = post_par_pool_incong$m_star, 
                                      sd = sqrt(post_par_pool_incong$v_star)),
                aes(color = "pool")) +
  stat_function(fun = function(x) dnorm(x, mean = post_par_incong$m_star, 
                                      sd = sqrt(post_par_incong$v_star)),
                aes(color = "curr")) +
  # stat_function(fun = function(x) dnorm(x, mean = post_par_crps_incong$m_star, 
  #                                     sd = sqrt(post_par_crps_incong$v_star)),
  #               aes(color = "crps")) +
  stat_function(fun = function(x) dnorm(x, mean = post_par_inf_incong$m_star, 
                                      sd = sqrt(post_par_inf_incong$v_star)),
                aes(color = "inf")) +
  geom_vline(aes(xintercept = mean(y0_incong), color = "y0", linetype = "y0")) +
  geom_vline(aes(xintercept = mean(y_incong),  color = "y",  linetype = "y")) +
  scale_color_manual(values = c(
                              "pool" = "#0079fbff",  # blue
                              "curr" = "#E45756",  # red
                              "crps" = "#72B7B2",  # teal
                              "inf"  = "#B279A2",  # muted purple
                              "y0"   = "#2E2E2E",  # dark gray (better than pure black)
                              "y"    = "#F58518"   # orange
                              ),
                      labels = c("pool" = "Pooled", 
                                  "curr" = "Curr", 
                                  "crps" = "CRPS", 
                                  "inf" = "Inf-Match",
                                  "y0" = "Hist Mean",
                                  "y" = "Curr Mean"
                              )
                            ) +
  scale_linetype_manual(values = c("y0" = "dashed", "y" = "dashed")) +
  labs(x = "", y = "", color = "") +
  ggtitle("Large discrepancy") +
  xlim(-1 - 1, 1 + 1) +
  theme_minimal() +
  guides(linetype = "none")
print(post_incong)

post_plot <- post_cong + post_scong + post_incong + plot_layout(ncol = 3) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")
print(post_plot)

ggsave("figures/inf_match_comp.png", post_plot, width = 9, height = 4)

results <- data.frame(
  scenario = c("No discrepancy", "Small discrepancy", "Large discrepancy"),
  eta_inf = round(c(eta_inf_cong, eta_inf_scong, eta_inf_incong), 3),
  eta_true_inf = round(c(eta_true_cong, eta_true_scong, eta_true_incong), 3),
  eta_crps = round(c(best_a0_crps, best_a0_crps_scong, best_a0_crps_incong), 3)
)
print(xtable::xtable(results), include.rownames = FALSE)








ggplot() +
  stat_function(fun = function(x) dnorm(x, mean = -1, 
                                      sd = 1),
                aes(color = "pool")) +
  stat_function(fun = function(x) dnorm(x, mean = 1, 
                                      sd = 1),
                aes(color = "curr")) +
  scale_color_manual(values = c(
                              "pool" = "#0079fbff",  # blue
                              "curr" = "#E45756",  # red
                              "crps" = "#72B7B2",  # teal
                              "inf"  = "#B279A2",  # muted purple
                              "y0"   = "#2E2E2E",  # dark gray (better than pure black)
                              "y"    = "#F58518"   # orange
                              ),
                      labels = c("pool" = "Pooled", 
                                  "curr" = "Curr", 
                                  "crps" = "CRPS", 
                                  "inf" = "Inf-Match",
                                  "y0" = "Hist Mean",
                                  "y" = "Curr Mean"
                              )
                            ) +
  scale_linetype_manual(values = c("y0" = "dashed", "y" = "dashed")) +
  labs(x = "", y = "", color = "") +
  # ggtitle("Large discrepancy") +
  xlim(-1 - 3, 1 + 3) +
  theme_minimal() +
  guides(linetype = "none", color = "none")
                              
