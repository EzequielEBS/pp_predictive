library(dplyr)

source("code/normal/aux_fun_normal.R")

# Runs the influence-matching eta estimator for one (mu0, mu1) discrepancy
# scenario across a grid of current-sample sizes n_list (with historical
# sizes n0_list, typically r * n_list), replicated n_rep times, and returns
# the median + 95% (non-parametric, binomial-based) CI of the estimated eta
# at each n.
#
# This is the one place the triple-nested "generate data / estimate eta /
# collect into a matrix" loop lives -- it used to be copy-pasted three times
# per file (once per scenario) in both inf_match_conv.r and
# inf_match_nu_analysis.r (which further repeated the whole triple for every
# nu value it swept over).
run_convergence_scenario <- function(n_list, n0_list, mu0, mu1, v = 1,
                                      n_rep = 1000, nu = .5, m0 = 0, v0 = 1) {
  results <- lapply(1:n_rep, function(j) {
    vapply(seq_along(n_list), function(i) {
      n <- n_list[i]
      n0 <- n0_list[i]
      data <- generate_normal_data(n0 = n0, n = n, mu0 = mu0, mu = mu1, sigma0 = v, sigma = v)
      y <- data %>% filter(data == "curr") %>% pull(y)
      post_par <- post_par_fixed_var(m0, v0, v, y)
      estimate_eta(data, post_par, v, mle = TRUE, nu = nu)
    }, numeric(1))
  })
  results <- do.call(rbind, results)
  colnames(results) <- paste0("n_", n_list)

  medians <- apply(results, 2, median)
  cis <- apply(results, 2, med_np)

  data.frame(
    n = n_list,
    median = medians,
    lower_ci = cis[2, ],
    upper_ci = cis[3, ]
  )
}

# Runs run_convergence_scenario() across a named list of scenarios (each a
# c(mu0=, mu1=) pair) and returns one combined, ready-to-plot data frame
# with a `scenario` factor column (levels in the order given).
run_convergence_sim <- function(n_list, r = 1.5, n_rep = 1000, nu = .5,
                                 scenarios) {
  n0_list <- r * n_list
  plot_data <- lapply(names(scenarios), function(label) {
    mu0 <- unname(scenarios[[label]]["mu0"])
    mu1 <- unname(scenarios[[label]]["mu1"])
    res <- run_convergence_scenario(n_list, n0_list, mu0, mu1, n_rep = n_rep, nu = nu)
    res$scenario <- label
    res
  })
  plot_data <- do.call(rbind, plot_data)
  plot_data$scenario <- factor(plot_data$scenario, levels = names(scenarios))
  plot_data
}

# Renders the standard faceted convergence plot (median +/- 95% CI vs n, one
# facet per scenario) from a data frame produced by run_convergence_sim().
#
# group = scenario on geom_line() matters: with x mapped to a factor and no
# other discrete aesthetic, ggplot2 groups by x itself, so each x-level
# becomes its own one-point "group" and the connecting line silently fails
# to draw (ggplot2 >= ~3.3 warns "each group consists of only one
# observation"). Confirmed via ggplot_build() that omitting this drops the
# line entirely -- this affected the original (pre-refactor) inline
# ggplot() calls in both inf_match_conv.r and inf_match_nu_analysis.r too.
plot_convergence <- function(plot_data, title = NULL) {
  p <- ggplot(plot_data, aes(x = as.factor(n), y = median, group = scenario)) +
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2) +
    facet_wrap(~ scenario, scales = "free_y", ncol = 1) +
    labs(x = "Sample Size (n)", y = expression(hat(eta))) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(strip.text = element_text(size = 12))
  if (!is.null(title)) p <- p + ggtitle(title)
  p
}
