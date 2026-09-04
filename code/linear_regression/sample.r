library(dplyr)
library(tidyverse)
library(hdbayes)
library(posterior)
library(subsampling)

set.seed(20260819)

# This simulated-data CSV lives outside this repo (in a sibling "onpp"
# checkout). The path below was previously hardcoded inline, which only
# works from one specific machine/working directory. Point data_dir at
# wherever that data actually lives for you (e.g. an absolute path, or set
# the ONPP_DATA_DIR environment variable) before running this script.
data_dir <- Sys.getenv("ONPP_DATA_DIR", "../../onpp/simulated-data-regression/data")

num_sim <- 200
data_high_cong_varying_beta_n0_ge_n <-
  read_csv(file.path(data_dir, "sim_data_high_cong_varying_beta_p3_n0_ge_n.csv")) %>%
  filter(replicate <= num_sim)
data_small_cong_varying_beta_n0_ge_n <-
  read_csv(file.path(data_dir, "sim_data_small_cong_varying_beta_p3_n0_ge_n.csv")) %>%
  filter(replicate <= num_sim)
data_no_cong_varying_beta_n0_ge_n <-
  read_csv(file.path(data_dir, "sim_data_no_cong_varying_beta_p3_n0_ge_n.csv")) %>%
  filter(replicate <= num_sim)

source("code/aux_fun_inf_match_glm.r")

ncores        = 4
chains        = 4      ## number of Markov chains to run
iter_warmup   = 5000   ## warmup per chain for MCMC sampling
iter_sampling = 10000   ## number of samples post warmup per chain
formula = y ~ X1 + X2
family  = gaussian()

#------------------------------------------
# No discrepancy setting
#------------------------------------------

hist_no_disc <- data_high_cong_varying_beta_n0_ge_n %>%
  filter(data_id == "hist3") %>%
  filter(replicate == 1)
curr_no_disc <- data_high_cong_varying_beta_n0_ge_n %>%
  filter(data_id == "current") %>%
  filter(replicate == 1)

data_list_no_disc <- list(curr_no_disc, hist_no_disc)

fit0_no_disc = glm.pp(
  formula = formula, family = family, data.list = data_list_no_disc,
  a0.vals = 0,
  iter_warmup = iter_warmup, iter_sampling = iter_sampling, 
  chains = chains, parallel_chains = ncores,
  refresh = 0
)
beta_draws_no_disc <- fit0_no_disc %>% select(-lp__, -dispersion) %>% as_draws_matrix()

etas_no_disc <- lapply(1:100, function(i) {
  estimate_eta_glm(
    formula = formula,
    curr_data = data_list_no_disc[[1]],
    hist_data = data_list_no_disc[[2]],
    beta_draws = beta_draws_no_disc,
    family = family
  )
})

#------------------------------------------
# Small discrepancy setting
#------------------------------------------

hist_small_disc <- data_small_cong_varying_beta_n0_ge_n %>%
  filter(data_id == "hist3") %>%
  filter(replicate == 1)
curr_small_disc <- data_small_cong_varying_beta_n0_ge_n %>%
  filter(data_id == "current") %>%
  filter(replicate == 1)

data_list_small_disc <- list(curr_small_disc, hist_small_disc)
fit0_small_disc = glm.pp(
  formula = formula, family = family, data.list = data_list_small_disc,
  a0.vals = 0,
  iter_warmup = iter_warmup, iter_sampling = iter_sampling, 
  chains = chains, parallel_chains = ncores,
  refresh = 0
)
beta_draws_small_disc <- fit0_small_disc %>% select(-lp__, -dispersion) %>% as_draws_matrix()
etas_small_disc <- lapply(1:100, function(i) {
  estimate_eta_glm(
    formula = formula,
    curr_data = data_list_small_disc[[1]],
    hist_data = data_list_small_disc[[2]],
    beta_draws = beta_draws_small_disc,
    family = family
  )
})

#------------------------------------------
# Large discrepancy setting
#------------------------------------------

hist_large_disc <- data_no_cong_varying_beta_n0_ge_n %>%
  filter(data_id == "hist3") %>%
  filter(replicate == 1)
curr_large_disc <- data_no_cong_varying_beta_n0_ge_n %>%
  filter(data_id == "current") %>%
  filter(replicate == 1)
data_list_large_disc <- list(curr_large_disc, hist_large_disc)
fit0_large_disc = glm.pp(
  formula = formula, family = family, data.list = data_list_large_disc,
  a0.vals = 0,
  iter_warmup = iter_warmup, iter_sampling = iter_sampling, 
  chains = chains, parallel_chains = ncores,
  refresh = 0
)
beta_draws_large_disc <- fit0_large_disc %>% select(-lp__, -dispersion) %>% as_draws_matrix()
etas_large_disc <- lapply(1:100, function(i) {
  estimate_eta_glm(
    formula = formula,
    curr_data = data_list_large_disc[[1]],
    hist_data = data_list_large_disc[[2]],
    beta_draws = beta_draws_large_disc,
    family = family
  )
})

#------------------------------------------
# Plotting
#------------------------------------------

etas_df <- data.frame(
  eta = c(unlist(etas_no_disc), unlist(etas_small_disc), unlist(etas_large_disc)),
  setting = factor(c(
    rep("No discrepancy", length(etas_no_disc) * length(etas_no_disc[[1]])),
    rep("Small discrepancy", length(etas_small_disc) * length(etas_small_disc[[1]])),
    rep("Large discrepancy", length(etas_large_disc) * length(etas_large_disc[[1]]))
  ), levels = c("No discrepancy", "Small discrepancy", "Large discrepancy"))
)

plot_etas <- ggplot(etas_df, aes(x = eta, fill = setting)) +
  geom_density(alpha = 0.5) +
  labs(x = expression(hat(eta)), y = "", fill = "Setting") +
  geom_vline(xintercept = median(unlist(etas_no_disc)), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = median(unlist(etas_small_disc)), linetype = "dashed", color = "orange") +
  geom_vline(xintercept = median(unlist(etas_large_disc)), linetype = "dashed", color = "skyblue") +
  scale_fill_manual(values = c("blue", "orange", "skyblue")) +
  theme_bw()
plot_etas

ggsave("figures/eta_estimates_lm.png", plot_etas, width = 6, height = 4) 
