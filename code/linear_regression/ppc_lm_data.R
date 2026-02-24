library(scoringRules)
library(mvtnorm)
library(ggplot2)
library(hdbayes)
library(dplyr)
library(MASS)
library(gridExtra)
library(LaplacesDemon)
library(patchwork)

# source("code/linear_regression/generate_lm_data.R")
load("samples/draws_beta_npp_post.RData")
load("samples/draws_beta_npp_prior.RData")
load("samples/draws_eta_post.RData")
load("samples/draws_eta_prior.RData")
load("samples/draws_priorpred_npp.RData")
load("samples/ess.RData")
load("data/sim_lm_data.RData")
source("code/linear_regression/aux_fun_lm.R")

formula <- y ~ X1
family      <- gaussian()
true_beta_hist <- c(-0.4, 0.5)
# true_beta_hist <- c(-0.2, 0.1)
true_sigma_hist <- 1
true_beta_curr <- c(-0.4, 0.5)
true_sigma_curr <- 1

res_hist <- hdbayes:::stack.data(formula = formula, data.list = list(hist_data))
res_curr <- hdbayes:::stack.data(formula = formula, data.list = list(curr_data))
y0 <- res_hist$y
X0 <- res_hist$X
n0 <- length(y0)
y <- res_curr$y
X <- res_curr$X
n <- length(y)
p <- ncol(X)


mu <- rep(0, p)
S <-  diag(p)
a <- 2 
b <- 1
a_tilde <- 1
b_tilde <- 1

hyva_t <- function(y, m, v, nu){
  out <- lapply(1:length(y), function(i){
    x_i <- y[i]
    m_i <- m[i]
    v_i <- v[i]
    tau_i <- sqrt(v_i)
    return(
      2 * (nu+1) * ((x_i-m_i)^2 - nu*1/v_i) / ((x_i-m_i)^2 + nu*1/v_i^2)^2 +
      (nu + 1)^2 * (x_i-m_i)^2 / ((x_i-m_i)^2 + nu*1/v_i^2)^2 
  )
  })
  return(mean(unlist(out)))
}

comp_crps <- function(eta) {
  single_crps <- function(e) {
    pp_prior_par <- pp_hyper_conj_lm(e, mu, S, a, b, X0, y0)
    mu_star <- pp_prior_par$mu_star
    S_star  <- pp_prior_par$S_star
    a_star  <- pp_prior_par$a_star
    b_star  <- pp_prior_par$b_star
    V_tilde <- diag(n)
    X_tilde <- X
    y_tilde <- y
    pred_par <- pred_par_conj_lm(mu_star, S_star, a_star, b_star, V_tilde, X_tilde, y_tilde)

    mean(crps_t(y_tilde, pred_par$nu_pred, pred_par$mu_pred, diag(pred_par$S_pred)))
  }

  vapply(eta, single_crps, numeric(1))
}

crps_optim <- optimize(comp_crps, interval = c(0,1))
best_a0_crps <- crps_optim$minimum
png("figures/crps_curve_lm.png", width = 6, height = 4, units = "in", res = 320)
curve(comp_crps(x), from = 0, to = 1, xlab = expression(eta), ylab = 'CRPS')
dev.off()

comp_hyvarinen <- function(eta) {
  single_hyva <- function(e) {
    pp_prior_par <- pp_hyper_conj_lm(e, mu, S, a, b, X0, y0)
    mu_star <- pp_prior_par$mu_star
    S_star  <- pp_prior_par$S_star
    a_star  <- pp_prior_par$a_star
    b_star  <- pp_prior_par$b_star
    V_tilde <- diag(n)
    X_tilde <- X
    y_tilde <- y

    pred_par <- pred_par_conj_lm(mu_star, S_star, a_star, b_star, V_tilde, X_tilde, y_tilde)

    S_pred  <- pred_par$S_pred
    mu_pred <- pred_par$mu_pred
    df      <- pred_par$nu_pred

    hyva_t(y_tilde, mu_pred, diag(S_pred), df)
  }

  vapply(eta, single_hyva, numeric(1))
}

hyva_optim <- optimize(comp_hyvarinen, interval = c(0,1))
best_a0_hyvarinen <- hyva_optim$minimum

png("figures/hyvarinen_curve_lm.png", width = 6, height = 4, units = "in", res = 320)
curve(comp_hyvarinen(x), from = 1e-4, to = 1, xlab = expression(eta), ylab = 'Hyvarinen score')
dev.off()

# Create legend names
name_crps <- "CRPS"
name_hyv  <- "Hyva"

# Wrap vectors as data frames
obs <- 1
# Observed value for vertical line
obs_val <- curr_data$y[obs]

name_a01 <- "eta = 1"
name_a00 <- "eta = 0"

# mytable <- cbind(
#   Dist = c("CRPS", "Hyvarinen", "NPP", "Full borrowing", "No borrowing"),  # remove "expression(...)" wrapper
#   pESS = c(round(ess_eta[[which.min(crps_list)]], 2),
#           round(ess_eta[[which.min(hyvarinen_list)]], 2),
#           round(ess_npp, 2),
#           round(ess_eta[[40]], 2),
#           round(ess_eta[[1]], 2))
# )

# tg <- tableGrob(
#   mytable,
#   rows = NULL,
#   theme = ttheme_default(
#     core = list(
#       fg_params = list(parse = TRUE, fontsize = 12)  # parse plotmath
#     ),
#     colhead = list(
#       fg_params = list(parse = TRUE, fontsize = 12)
#     )
#   )
# )

pp_par_crps <- pp_hyper_conj_lm(best_a0_crps, mu, S, a, b, X0, y0)
pp_par_hyva <- pp_hyper_conj_lm(best_a0_hyvarinen, mu, S, a, b, X0, y0)
pp_par_eta1 <- pp_hyper_conj_lm(1, mu, S, a, b, X0, y0)
pp_par_eta0 <- pp_hyper_conj_lm(0, mu, S, a, b, X0, y0)

pred_par_crps <- pred_par_conj_lm(pp_par_crps$mu_star, 
                                  pp_par_crps$S_star, 
                                  pp_par_crps$a_star, 
                                  pp_par_crps$b_star, 
                                  diag(n), 
                                  X, 
                                  y)
pred_par_hyva <- pred_par_conj_lm(pp_par_hyva$mu_star, 
                                  pp_par_hyva$S_star, 
                                  pp_par_hyva$a_star, 
                                  pp_par_hyva$b_star, 
                                  diag(n), 
                                  X, 
                                  y)
pred_par_eta1 <- pred_par_conj_lm(pp_par_eta1$mu_star,
                                  pp_par_eta1$S_star, 
                                  pp_par_eta1$a_star, 
                                  pp_par_eta1$b_star, 
                                  diag(n), 
                                  X, 
                                  y)
pred_par_eta0 <- pred_par_conj_lm(pp_par_eta0$mu_star,
                                  pp_par_eta0$S_star,
                                  pp_par_eta0$a_star,
                                  pp_par_eta0$b_star,
                                  diag(n),
                                  X,
                                  y)

# prior distributions
plot_beta1 <- ggplot() +
  geom_function(fun = function(x) {
    dst(x,
        mu=pp_par_crps$mu_star[1],
        sigma=pp_par_crps$S_star[1,1] * pp_par_crps$b_star / pp_par_crps$b_star,
        nu=2*pp_par_crps$a_star)
  },
  aes(color = "CRPS"),
  linewidth = 2.5,
  linetype = 'dotted'
  ) +
  geom_function(fun = function(x) {
    dst(x,
        mu=pp_par_hyva$mu_star[1],
        sigma=pp_par_hyva$S_star[1,1] * pp_par_hyva$b_star / pp_par_hyva$a_star,
        nu=pp_par_hyva$a_star * 2)
  },
  aes(color = "Hyva"),
  linewidth = 1
  ) +
  geom_function(fun = function(x) {
    dst(x,
        mu=pp_par_eta1$mu_star[1],
        sigma=pp_par_eta1$S_star[1,1] * pp_par_eta1$b_star / pp_par_eta1$a_star,
        nu=pp_par_eta1$a_star * 2)
  },
  aes(color = "eta = 1"),
  linetype = 'dashed',
  linewidth = 1.6
  ) +
  geom_function(fun = function(x) {
    dst(x,
        mu=pp_par_eta0$mu_star[1],
        sigma=pp_par_eta0$S_star[1,1] * pp_par_eta0$b_star / pp_par_eta0$a_star,
        nu=pp_par_eta0$a_star * 2)
  },
  aes(color = "eta = 0"),
  linetype = 'dotted',
  linewidth = 1.6
  )+
  geom_vline(aes(xintercept = true_beta_curr[1], color = 'beta_curr'), 
             linewidth = 0.5, linetype = "solid") +
  geom_vline(aes(xintercept = true_beta_hist[1], color = 'beta_hist'), 
             linewidth = 1.2, linetype = "dotted") +
  labs(x = expression(beta[0]), y = '') +
  theme_bw() +
  geom_density(data = data.frame(y = draws_beta_npp_prior[,1]),
               aes(x = y, color = "NPP"),
               linewidth = 1) +
  scale_color_manual(
    name = NULL,
    values = c(
      "CRPS"    = "#66A8D0",
      "Hyva"    = "#D06673",
      "eta = 0" = "#7f7f7f",
      "eta = 1" = "#7f7f7f",
      "NPP" = "#D0C366",
      "beta_curr" = "black",
      "beta_hist" = "brown"
    ),
    labels = c(
      "CRPS"    = "CRPS",
      "Hyva"    = "Hyvarinen",
      "eta = 0" = "No \nborrowing",
      "eta = 1" = "Full \nborrowing",
      "beta_curr" = expression(beta['curr']^~'true'),
      "beta_hist" = expression(beta['hist']^~'true')
    ),
    breaks = c("CRPS",
               "Hyva",
               "NPP",
               "eta = 1",
               "eta = 0",
               "beta_hist",
               "beta_curr")
  ) +
  xlim(c(true_beta_curr[1]-0.5,true_beta_curr[1]+0.5)) +
  theme(text = element_text(size = 12),        # Base text size
        axis.title = element_text(size = 14),  # Axis titles
        axis.text = element_text(size = 12),   # Axis tick labels
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 11),
        legend.position = c(0.8, 0.74),
        legend.background = element_rect(
          fill = "white",     # background color of the legend
          color = "black",    # border color
          size = 0.3,         # border thickness
          linetype = "solid"  # border type
        ),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA))

plot_beta2 <- ggplot() +
  geom_function(fun = function(x) {
    dst(x,
        mu=pp_par_crps$mu_star[2],
        sigma=pp_par_crps$S_star[2,2] * pp_par_crps$b_star / pp_par_crps$a_star,
        nu=pp_par_crps$a_star * 2)
  },
  aes(color = "CRPS"),
  linewidth = 2.5,
  linetype = 'dotted'
  ) +
  geom_function(fun = function(x) {
    dst(x,
        mu=pp_par_hyva$mu_star[2],
        sigma=pp_par_hyva$S_star[2,2] * pp_par_hyva$b_star / pp_par_hyva$a_star,
        nu=pp_par_hyva$a_star * 2)
  },
  aes(color = "Hyva"),
  linewidth = 1
  ) +
  geom_function(fun = function(x) {
    dst(x,
        mu=pp_par_eta1$mu_star[2],
        sigma=pp_par_eta1$S_star[2,2] * pp_par_eta1$b_star / pp_par_eta1$a_star,
        nu=pp_par_eta1$a_star * 2)
  },
  aes(color = "eta = 1"),
  linetype = 'dashed',
  linewidth = 1.6
  ) +
  geom_function(fun = function(x) {
    dst(x,
        mu=pp_par_eta0$mu_star[2],
        sigma=pp_par_eta0$S_star[2,2] * pp_par_eta0$b_star / pp_par_eta0$a_star,
        nu=pp_par_eta0$a_star * 2)
  },
  aes(color = "eta = 0"),
  linetype = 'dotted',
  linewidth = 1.6
  )+
  geom_vline(aes(xintercept = true_beta_curr[2], color = 'beta_curr'), 
             linewidth = 0.5, linetype = "solid") +
  geom_vline(aes(xintercept = true_beta_hist[2], color = 'beta_hist'), 
             linewidth = 1.2, linetype = "dotted") +
  labs(x = expression(beta[1]), y = '') +
  theme_bw() +
  geom_density(data = data.frame(y = draws_beta_npp_prior[,2]),
               aes(x = y, color = "NPP"),
               linewidth = 1) +
  scale_color_manual(
    name = NULL,
    values = c(
      "CRPS"    = "#66A8D0",
      "Hyva"    = "#D06673",
      "eta = 0" = "#7f7f7f",
      "eta = 1" = "#7f7f7f",
      "NPP" = "#D0C366",
      "beta_curr" = "black",
      "beta_hist" = "brown"
    ),
    labels = c(
      "CRPS"    = "CRPS",
      "Hyva"    = "Hyvarinen",
      "eta = 0" = "No \nborrowing",
      "eta = 1" = "Full \nborrowing",
      "beta_curr" = expression(beta['curr']^~'true'),
      "beta_hist" = expression(beta['hist']^~'true')
    ),
    breaks = c("CRPS",
               "Hyva",
               "NPP",
               "eta = 1",
               "eta = 0",
               "beta_hist",
               "beta_curr")
  )  +
  xlim(c(true_beta_curr[2]-0.5,true_beta_curr[2]+0.5)) +
  theme(text = element_text(size = 12),        # Base text size
        axis.title = element_text(size = 14),  # Axis titles
        axis.text = element_text(size = 12),   # Axis tick labels
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 11),
        legend.position = "none",
        legend.background = element_rect(
          fill = "white",     # background color of the legend
          color = "black",    # border color
          size = 0.3,         # border thickness
          linetype = "solid"  # border type
        ),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA))

plot_eta <- ggplot() +
  geom_function(fun = function(x) {
    dbeta(x, a_tilde, b_tilde)
  }, aes(color = "NPP"), linewidth = 1) +
  geom_vline(aes(xintercept = best_a0_crps, color = 'best_crps'), 
                  linewidth = 2.5,
                  linetype = 'dotted') +
  geom_vline(aes(xintercept = best_a0_hyvarinen, color = 'best_hyva'), 
                  linewidth = 0.5, linetype = 'solid') +
  labs(x = expression(eta), y = '') +
  theme_bw() +
  scale_color_manual(
    name = NULL,
    values = c(
      "NPP" = "#D0C366",
      "best_crps" = "#66A8D0",
      "best_hyva" = "#D06673"
    ),
    labels = c(
      "best_crps" = "CRPS",
      "best_hyva" = "Hyvarinen"
    ),
    breaks = c("best_crps",
               "best_hyva")
  ) +
  xlim(c(0,1)) +
  theme(text = element_text(size = 12),        # Base text size
        axis.title = element_text(size = 14),  # Axis titles
        axis.text = element_text(size = 12),   # Axis tick labels
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 11),
        legend.position = c(0.2, 0.1),
        legend.background = element_rect(
          fill = "white",     # background color of the legend
          color = "black",    # border color
          size = 0.3,         # border thickness
          linetype = "solid"  # border type
        ),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA))

plot_par <- plot_beta1 + plot_beta2 + plot_eta + plot_layout(ncol = 3)
plot_par
ggsave(plot = plot_par,"figures/prior_lm.png", width = 15, height = 5, dpi = 320)

# prior predictive 
ggplot() +
  geom_function(fun = function(x) {
    dst(x, 
        mu=pred_par_crps$mu_pred[obs],
        sigma=pred_par_crps$S_pred[obs,obs], 
        nu=pred_par_crps$nu_pred)
    },
  aes(color = "CRPS"),
  linewidth = 2.5,
  linetype = 'dotted'
  ) +
  geom_function(fun = function(x) {
    dst(x, 
        mu=pred_par_hyva$mu_pred[obs],
        sigma=pred_par_hyva$S_pred[obs,obs], 
        nu=pred_par_hyva$nu_pred)
  },
  aes(color = "Hyva"),
  linewidth = 1
  ) +
  geom_function(fun = function(x) {
    dst(x, 
        mu=pred_par_eta1$mu_pred[obs],
        sigma=pred_par_eta1$S_pred[obs,obs], 
        nu=pred_par_eta1$nu_pred)
  },
  aes(color = "eta = 1"),
  linetype = 'dashed',
  linewidth = 1.6
  ) +
  geom_function(fun = function(x) {
    dst(x, 
        mu=pred_par_eta0$mu_pred[obs],
        sigma=pred_par_eta0$S_pred[obs,obs], 
        nu=pred_par_eta0$nu_pred)
  },
  aes(color = "eta = 0"),
  linetype = 'dotted',
  linewidth = 1.6
  )+
  geom_function(fun = function(x) {
    dnorm(x, mean = sum(c(1,curr_data[obs,2])*true_beta_curr), sd = true_sigma_curr)
  },
  aes(color = "true"),
  linetype = 'dashed',
  linewidth = 1.6
  ) + 
  geom_vline(aes(xintercept = obs_val, color = 'Obs value'), linewidth = 0.5) +
  labs(x = expression(tilde(y)), y = '') +
  theme_bw() +
  geom_density(data = data.frame(y = draws_priorpred_npp[,obs]), 
               aes(x = y, color = "NPP"),
               linewidth = 1) +
  scale_color_manual(
    name = NULL,
    values = c(
      "CRPS"    = "#66A8D0",
      "Hyva"    = "#D06673",
      "eta = 0" = "#7f7f7f",
      "eta = 1" = "#7f7f7f",
      "NPP" = "#D0C366",
      "true" = "black",
      "Obs value" = "black"
    ),
    labels = c(
      "CRPS"    = "CRPS",
      "Hyva"    = "Hyvarinen",
      "eta = 0" = "No \nborrowing",
      "eta = 1" = "Full \nborrowing",
      "true" = "True \ndist",
      "Obs value" = expression(tilde(y)[obs])
    ),
    breaks = c("CRPS",
               "Hyva",
               "NPP",
               "eta = 1",
               "eta = 0",
               "true",
               "Obs value")
  ) +
  # geom_label(
  #   data = ess_labels,
  #   aes(x, y, label = label),
  #   inherit.aes = FALSE,
  #   fill = "white",   # box color
  #   color = "black",    # text color
  #   hjust = 0,
  #   size = 4,
  #   parse = TRUE
  # ) +
  # annotation_custom(tg, xmin=-5.4, xmax = -3.3, ymin=0.25, ymax=0.65)+
  xlim(c(obs_val-4,obs_val+4)) +
  theme(text = element_text(size = 12),        # Base text size
        axis.title = element_text(size = 14),  # Axis titles
        axis.text = element_text(size = 12),   # Axis tick labels
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 11),
        legend.position = c(0.9, 0.74),
        legend.background = element_rect(
          fill = "white",     # background color of the legend
          color = "black",    # border color
          size = 0.3,         # border thickness
          linetype = "solid"  # border type
        ),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA))

ggsave("figures/ppc_npp_lm.png", width = 8, height = 5, dpi = 320)
# ggsave("figures/ppc_npp_lm.pdf", width = 8, height = 5, dpi = 320, device = cairo_pdf)


# post distributions
pp_post_par_crps <- pp_post_par_conj_lm(best_a0_crps, mu, S, a, b, X0, y0, diag(n), X, y)
pp_post_par_hyva <- pp_post_par_conj_lm(best_a0_hyvarinen, mu, S, a, b, X0, y0, diag(n), X, y)
pp_post_par_eta1 <- pp_post_par_conj_lm(1, mu, S, a, b, X0, y0, diag(n), X, y)
pp_post_par_eta0 <- pp_post_par_conj_lm(0, mu, S, a, b, X0, y0, diag(n), X, y)

plot_post_beta1 <- ggplot() +
  geom_function(fun = function(x) {
    dst(x,
        mu=pp_post_par_crps$mu_star[1],
        sigma=pp_post_par_crps$S_star[1,1] * pp_post_par_crps$b_star / pp_post_par_crps$a_star,
        nu=pp_post_par_crps$a_star * 2)
  },
  aes(color = "CRPS"),
  linewidth = 2.5,
  linetype = 'dotted'
  ) +
  geom_function(fun = function(x) {
    dst(x,
        mu=pp_post_par_hyva$mu_star[1],
        sigma=pp_post_par_hyva$S_star[1,1] * pp_post_par_hyva$b_star / pp_post_par_hyva$a_star,
        nu=pp_post_par_hyva$a_star * 2)
  },
  aes(color = "Hyva"),
  linewidth = 1
  ) +
  geom_function(fun = function(x) {
    dst(x,
        mu=pp_post_par_eta1$mu_star[1],
        sigma=pp_post_par_eta1$S_star[1,1] * pp_post_par_eta1$b_star / pp_post_par_eta1$a_star,
        nu=pp_post_par_eta1$a_star * 2)
  },
  aes(color = "eta = 1"),
  linetype = 'dashed',
  linewidth = 1.6
  ) +
  geom_function(fun = function(x) {
    dst(x,
        mu=pp_post_par_eta0$mu_star[1],
        sigma=pp_post_par_eta0$S_star[1,1] * pp_post_par_eta0$b_star / pp_post_par_eta0$a_star,
        nu=pp_post_par_eta0$a_star * 2)
  },
  aes(color = "eta = 0"),
  linetype = 'dotted',
  linewidth = 1.6
  )+
  geom_vline(aes(xintercept = true_beta_curr[1], color = 'beta_curr'), 
             linewidth = 0.5, linetype = "solid") +
  geom_vline(aes(xintercept = true_beta_hist[1], color = 'beta_hist'), 
             linewidth = 1.2, linetype = "dotted") +
  labs(x = expression(beta[0]), y = '') +
  theme_bw() +
  geom_density(data = data.frame(y = draws_beta_npp_post[,1]),
               aes(x = y, color = "NPP"),
               linewidth = 1) +
  scale_color_manual(
    name = NULL,
    values = c(
      "CRPS"    = "#66A8D0",
      "Hyva"    = "#D06673",
      "eta = 0" = "#7f7f7f",
      "eta = 1" = "#7f7f7f",
      "NPP" = "#D0C366",
      "beta_curr" = "black",
      "beta_hist" = "brown"
    ),
    labels = c(
      "CRPS"    = "CRPS",
      "Hyva"    = "Hyvarinen",
      "eta = 0" = "No \nborrowing",
      "eta = 1" = "Full \nborrowing",
      "beta_curr" = expression(beta['curr']^~'true'),
      "beta_hist" = expression(beta['hist']^~'true')
    ),
    breaks = c("CRPS",
               "Hyva",
               "NPP",
               "eta = 1",
               "eta = 0",
               "beta_hist",
               "beta_curr")
  ) +
  xlim(c(true_beta_curr[1]-0.5,true_beta_curr[1]+0.5)) +
  theme(text = element_text(size = 12),        # Base text size
        axis.title = element_text(size = 14),  # Axis titles
        axis.text = element_text(size = 12),   # Axis tick labels
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 11),
        legend.position = c(0.8, 0.74),
        legend.background = element_rect(
          fill = "white",     # background color of the legend
          color = "black",    # border color
          size = 0.3,         # border thickness
          linetype = "solid"  # border type
        ),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA))

plot_post_beta2 <- ggplot() +
  geom_function(fun = function(x) {
    dst(x,
        mu=pp_post_par_crps$mu_star[2],
        sigma=pp_post_par_crps$S_star[2,2] * pp_post_par_crps$b_star / pp_post_par_crps$a_star,
        nu=pp_post_par_crps$a_star * 2)
  },
  aes(color = "CRPS"),
  linewidth = 2.5,
  linetype = 'dotted'
  ) +
  geom_function(fun = function(x) {
    dst(x,
        mu=pp_post_par_hyva$mu_star[2],
        sigma=pp_post_par_hyva$S_star[2,2] * pp_post_par_hyva$b_star / pp_post_par_hyva$a_star,
        nu=pp_post_par_hyva$a_star * 2)
  },
  aes(color = "Hyva"),
  linewidth = 1
  ) +
  geom_function(fun = function(x) {
    dst(x,
        mu=pp_post_par_eta1$mu_star[2],
        sigma=pp_post_par_eta1$S_star[2,2] * pp_post_par_eta1$b_star / pp_post_par_eta1$a_star,
        nu=pp_post_par_eta1$a_star * 2)
  },
  aes(color = "eta = 1"),
  linetype = 'dashed',
  linewidth = 1.6
  ) +
  geom_function(fun = function(x) {
    dst(x,
        mu=pp_post_par_eta0$mu_star[2],
        sigma=pp_post_par_eta0$S_star[2,2] * pp_post_par_eta0$b_star / pp_post_par_eta0$a_star,
        nu=pp_post_par_eta0$a_star * 2)
  },
  aes(color = "eta = 0"),
  linetype = 'dotted',
  linewidth = 1.6
  )+
  geom_vline(aes(xintercept = true_beta_curr[2], color = 'beta_curr'), 
             linewidth = 0.5, linetype = "solid") +
  geom_vline(aes(xintercept = true_beta_hist[2], color = 'beta_hist'), 
             linewidth = 1.2, linetype = "dotted") +
  labs(x = expression(beta[1]),
       y = '') +
  theme_bw() +
  geom_density(data = data.frame(y = draws_beta_npp_post[,2]),
               aes(x = y, color = "NPP"),
               linewidth = 1) +
  scale_color_manual(
    name = NULL,
    values = c(
      "CRPS"    = "#66A8D0",
      "Hyva"    = "#D06673",
      "eta = 0" = "#7f7f7f",
      "eta = 1" = "#7f7f7f",
      "NPP" = "#D0C366",
      "beta_curr" = "black",
      "beta_hist" = "brown"
    ),
    labels = c(
      "CRPS"    = "CRPS",
      "Hyva"    = "Hyvarinen",
      "eta = 0" = "No \nborrowing",
      "eta = 1" = "Full \nborrowing",
      "beta_curr" = expression(beta['curr']^~'true'),
      "beta_hist" = expression(beta['hist']^~'true')
    ),
    breaks = c("CRPS",
               "Hyva",
               "NPP",
               "eta = 1",
               "eta = 0",
               "beta_hist",
               "beta_curr")
  ) +
  xlim(c(true_beta_curr[2]-0.5,true_beta_curr[2]+0.5)) +
  theme(text = element_text(size = 12),        # Base text size
        axis.title = element_text(size = 14),  # Axis titles
        axis.text = element_text(size = 12),   # Axis tick labels
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 11),
        legend.position = "none",
        legend.background = element_rect(
          fill = "white",     # background color of the legend
          color = "black",    # border color
          size = 0.3,         # border thickness
          linetype = "solid"  # border type
        ),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA))

plot_post_eta <- ggplot() +
  geom_density(data = data.frame(y = draws_eta_post), aes(x = y, color = "NPP"), linewidth = 1) +
  geom_vline(aes(xintercept = best_a0_crps, color = 'best_crps'), 
                  linewidth = 2.5, linetype = 'dotted') +
  geom_vline(aes(xintercept = best_a0_hyvarinen, color = 'best_hyva'), linewidth = 0.5, linetype = 'dashed') +
  labs(x = expression(eta), y = '') +
  theme_bw() +
  scale_color_manual(
    name = NULL,
    values = c(
      "NPP" = "#D0C366",
      "best_crps" = "#66A8D0",
      "best_hyva" = "#D06673"
    ),
    labels = c(
      "best_crps" = "CRPS",
      "best_hyva" = "Hyvarinen"
    ),
    breaks = c("best_crps",
               "best_hyva")
  ) +
  xlim(c(0,1)) +
  theme(text = element_text(size = 12),        # Base text size
        axis.title = element_text(size = 14),  # Axis titles
        axis.text = element_text(size = 12),   # Axis tick labels
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 11),
        legend.position = c(0.2, 0.1),
        legend.background = element_rect(
          fill = "white",     # background color of the legend
          color = "black",    # border color
          size = 0.3,         # border thickness
          linetype = "solid"  # border type
        ),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA))

plot_par_post <- plot_post_beta1 + plot_post_beta2 + plot_post_eta + plot_layout(ncol = 3)
plot_par_post
ggsave(plot = plot_par_post,"figures/post_lm.png", width = 15, height = 5, dpi = 320)

