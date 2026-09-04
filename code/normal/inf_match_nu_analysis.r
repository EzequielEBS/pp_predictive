library(tidyverse)
library(dplyr)
library(ggplot2)
library(scoringRules)
library(LaplacesDemon)
library(pbapply)
library(parallel)
library(patchwork)

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
n0_list <- r*n_list
n_rep <- 1000
nu <- c(.99, .5, .01)

# Same three discrepancy scenarios as before (mu0 = 1 / 0 / -1, mu1 = 1
# throughout -- note these differ from the scenarios in inf_match_conv.r,
# which is preserved here rather than unified, since that's how the
# original scripts were parameterized).
scenarios <- list(
  "No discrepancy"    = c(mu0 = 1,  mu1 = 1),
  "Small discrepancy" = c(mu0 = 0,  mu1 = 1),
  "Large discrepancy" = c(mu0 = -1, mu1 = 1)
)

plots_nu <- pblapply(nu, function(nu_val) {
  plot_data <- run_convergence_sim(n_list, r = r, n_rep = n_rep, nu = nu_val, scenarios = scenarios)
  plot_convergence(plot_data, title = bquote(nu == .(nu_val)))
})

all_plots <- wrap_plots(plots_nu, ncol = 1)
all_plots
ggsave("figures/inf_match_nu_analysis.pdf", width = 10, height = 12)

# -------------------------------------------------------------------------------
# behavior near 1
#
# NOTE: preserved as in the original script -- plots_nu_near_one is computed
# but never saved or printed here. If this sweep was meant to produce a
# figure, add a ggsave()/print() call; left as a no-op rather than guessing
# an intended output file.
# -------------------------------------------------------------------------------

nu_near_one <- seq(0.8, 0.9, by = 0.01)

ncores <- max(1, parallel::detectCores() - 1)
cl <- makeCluster(ncores)
clusterSetRNGStream(cl, 20260819)

clusterEvalQ(cl, {
  library(dplyr)
  library(ggplot2)
  source("code/normal/sim_helpers.R")
})
clusterExport(cl,
  varlist = c(
    "n_list", "n0_list", "r", "n_rep", "nu_near_one", "scenarios"
  ),
  envir = .GlobalEnv
)

plots_nu_near_one <- pblapply(nu_near_one, function(nu_val) {
  plot_data <- run_convergence_sim(n_list, r = r, n_rep = n_rep, nu = nu_val, scenarios = scenarios)
  plot_convergence(plot_data, title = bquote(nu == .(nu_val)))
}, cl = cl)

stopCluster(cl)
