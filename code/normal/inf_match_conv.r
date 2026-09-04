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

source("code/normal/sim_helpers.R")

set.seed(20260819)

#-------------------------------------------------------------------------------
# analysis
#-------------------------------------------------------------------------------

n_list <- c(10, 50, 100, 200, 500, 1000, 10000)
r <- 1.5
n_rep <- 1000

# Same three discrepancy scenarios as before (mu0 = 1 / 0.5 / 0, mu1 = 1
# throughout), now expressed as data instead of three copy-pasted loops --
# see code/normal/sim_helpers.R for the shared simulation logic.
plot_data <- run_convergence_sim(
  n_list, r = r, n_rep = n_rep,
  scenarios = list(
    "No discrepancy"    = c(mu0 = 1,   mu1 = 1),
    "Small discrepancy" = c(mu0 = 0.5, mu1 = 1),
    "Large discrepancy" = c(mu0 = 0,   mu1 = 1)
  )
)

p <- plot_convergence(plot_data)
p
ggsave("figures/inf_match_conv.png", plot = p, width = 10, height = 4)
