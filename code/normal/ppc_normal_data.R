library(scoringRules)
library(mvtnorm)
library(ggplot2)
library(hdbayes)
library(dplyr)
library(MASS)
library(gridExtra)
library(LaplacesDemon)
library(patchwork)

source("code/normal/generate_normal_data.R")
load("samples/draws_mu_npp_post.RData")
load("samples/draws_mu_npp_prior.RData")
load("samples/draws_eta_post_normal.RData")
load("samples/draws_eta_prior_normal.RData")
load("samples/draws_priorpred_npp_normal.RData")
load("data/sim_normal_data.RData")
source("code/normal/aux_fun_normal.R")

true_mu_hist <- 1
true_sigma_hist <- 1
true_mu_curr <- 1
true_sigma_curr <- 1

y0 <- hist_data$y
n0 <- length(y0)
y <- curr_data$y
n <- length(y)


m <- 0
v <-  1
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

comp_crps <- function(eta){
  single <- function(e){
    pp_prior_par <- pp_hyper_conj_normal(e, m, v, a, b, y0)
    m_star <- pp_prior_par$m_star
    v_star <- pp_prior_par$v_star
    a_star <- pp_prior_par$a_star
    b_star <- pp_prior_par$b_star
    w_tilde <- 1
    y_tilde <- y
    pred_par <- pred_par_conj_normal(m_star, v_star, a_star, b_star, w_tilde, y_tilde)
    mean(crps_t(y_tilde, pred_par$nu_pred, pred_par$m_pred, pred_par$v_pred))
  }
  if (length(eta) > 1) sapply(eta, single) else single(eta)
}

crps_optim <- optimize(comp_crps, interval = c(0,1))
best_a0_crps <- crps_optim$minimum

png("figures/crps_curve_normal.png", width = 8, height = 5, units = "in", res = 320)
curve(comp_crps(x), from = 1e-4, to = 1, xlab = expression(eta), ylab = 'CRPS')
dev.off()

comp_hyvarinen <- function(eta) {
  single <- function(e) {
    pp_prior_par <- pp_hyper_conj_normal(e, m, v, a, b, y0)
    m_star <- pp_prior_par$m_star
    v_star <- pp_prior_par$v_star
    a_star <- pp_prior_par$a_star
    b_star <- pp_prior_par$b_star
    w_tilde <- 1
    y_tilde <- y

    pred_par <- pred_par_conj_normal(m_star, v_star, a_star, b_star, w_tilde, y_tilde)

    n <- length(y_tilde)
    v_pred <- pred_par$v_pred
    m_pred <- pred_par$m_pred
    df <- pred_par$nu_pred

    hyva_t(y_tilde, rep(m_pred, n), rep(v_pred, n), df)
  }

  if (length(eta) > 1) sapply(eta, single) else single(eta)
}

hyva_optim <- optimize(comp_hyvarinen, interval = c(0,1))
best_a0_hyvarinen <- hyva_optim$minimum

png("figures/hyvarinen_curve_normal.png", width = 8, height = 5, units = "in", res = 320)
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

pp_par_crps <- pp_hyper_conj_normal(best_a0_crps, m, v, a, b, y0)
pp_par_hyva <- pp_hyper_conj_normal(best_a0_hyvarinen, m, v, a, b, y0)
pp_par_eta1 <- pp_hyper_conj_normal(1, m, v, a, b, y0)
pp_par_eta0 <- pp_hyper_conj_normal(0, m, v, a, b, y0)

pred_par_crps <- pred_par_conj_normal(pp_par_crps$m_star, 
                                  pp_par_crps$v_star, 
                                  pp_par_crps$a_star, 
                                  pp_par_crps$b_star, 
                                  1, 
                                  y)
pred_par_hyva <- pred_par_conj_normal(pp_par_hyva$m_star, 
                                  pp_par_hyva$v_star, 
                                  pp_par_hyva$a_star, 
                                  pp_par_hyva$b_star, 
                                  1, 
                                  y)
pred_par_eta1 <- pred_par_conj_normal(pp_par_eta1$m_star,
                                  pp_par_eta1$v_star, 
                                  pp_par_eta1$a_star, 
                                  pp_par_eta1$b_star, 
                                  1, 
                                  y)
pred_par_eta0 <- pred_par_conj_normal(pp_par_eta0$m_star,
                                  pp_par_eta0$v_star,
                                  pp_par_eta0$a_star,
                                  pp_par_eta0$b_star,
                                  1,
                                  y)

# prior distributions
plot_mu <- ggplot() +
  geom_function(fun = function(x) {
    dst(x,
        mu=pp_par_crps$m_star,
        sigma=pp_par_crps$v_star * pp_par_crps$b_star / pp_par_crps$b_star,
        nu=2*pp_par_crps$a_star)
  },
  aes(color = "CRPS"),
  linewidth = 2.5,
  linetype = 'dotted'
  ) +
  geom_function(fun = function(x) {
    dst(x,
        mu=pp_par_hyva$m_star,
        sigma=pp_par_hyva$v_star * pp_par_hyva$b_star / pp_par_hyva$a_star,
        nu=pp_par_hyva$a_star * 2)
  },
  aes(color = "Hyva"),
  linewidth = 1
  ) +
  geom_function(fun = function(x) {
    dst(x,
        mu=pp_par_eta1$m_star,
        sigma=pp_par_eta1$v_star * pp_par_eta1$b_star / pp_par_eta1$a_star,
        nu=pp_par_eta1$a_star * 2)
  },
  aes(color = "eta = 1"),
  linetype = 'dashed',
  linewidth = 1.6
  ) +
  geom_function(fun = function(x) {
    dst(x,
        mu=pp_par_eta0$m_star,
        sigma=pp_par_eta0$v_star * pp_par_eta0$b_star / pp_par_eta0$a_star,
        nu=pp_par_eta0$a_star * 2)
  },
  aes(color = "eta = 0"),
  linetype = 'dotted',
  linewidth = 1.6
  )+
  geom_vline(aes(xintercept = true_mu_curr[1], color = 'beta_curr'), 
             linewidth = 0.5, linetype = "solid") +
  geom_vline(aes(xintercept = true_mu_hist[1], color = 'beta_hist'), 
             linewidth = 1.2, linetype = "dotted") +
  labs(x = expression(mu), y = '') +
  theme_bw() +
  geom_density(data = data.frame(y = draws_mu_npp_prior[,1]),
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
      "beta_curr" = expression(mu['curr']^~'true'),
      "beta_hist" = expression(mu['hist']^~'true')
    ),
    breaks = c("CRPS",
               "Hyva",
               "NPP",
               "eta = 1",
               "eta = 0",
               "beta_hist",
               "beta_curr")
  ) +
  xlim(c(true_mu_curr[1]-1,true_mu_curr[1]+1)) +
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
plot_mu

plot_eta <- ggplot() +
  geom_function(fun = function(x) {
    dbeta(x, a_tilde, b_tilde)
  }, aes(color = "NPP"), linewidth = 1) +
  geom_vline(aes(xintercept = best_a0_crps, color = 'best_crps'), 
                  linewidth = 2.5, linetype = 'dotted') +
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

plot_par <- plot_mu + plot_eta + plot_layout(ncol = 2)
plot_par
ggsave(plot = plot_par,"figures/prior_normal.png", width = 15, height = 5, dpi = 320)


# prior predictive 
ggplot() +
  geom_function(fun = function(x) {
    dst(x, 
        mu=pred_par_crps$m_pred,
        sigma=pred_par_crps$v_pred, 
        nu=pred_par_crps$nu_pred)
  },
  aes(color = "CRPS"),
  linewidth = 2.5,
  linetype = 'dotted'
  ) +
  geom_function(fun = function(x) {
    dst(x, 
        mu=pred_par_hyva$m_pred,
        sigma=pred_par_hyva$v_pred, 
        nu=pred_par_hyva$nu_pred)
  },
  aes(color = "Hyva"),
  linewidth = 1
  ) +
  geom_function(fun = function(x) {
    dst(x, 
        mu=pred_par_eta1$m_pred,
        sigma=pred_par_eta1$v_pred, 
        nu=pred_par_eta1$nu_pred)
  },
  aes(color = "eta = 1"),
  linetype = 'dashed',
  linewidth = 1.6
  ) +
  geom_function(fun = function(x) {
    dst(x, 
        mu=pred_par_eta0$m_pred,
        sigma=pred_par_eta0$v_pred, 
        nu=pred_par_eta0$nu_pred)
  },
  aes(color = "eta = 0"),
  linetype = 'dotted',
  linewidth = 1.6
  )+
  geom_function(fun = function(x) {
    dnorm(x, mean = true_mu_curr, sd = true_sigma_curr)
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


ggsave("figures/ppc_npp_normal.png", width = 8, height = 5, dpi = 320)


# post distributions
pp_post_par_crps <- pp_post_par_conj_normal(best_a0_crps, m, v, a, b, y0, 1, y)
pp_post_par_hyva <- pp_post_par_conj_normal(best_a0_hyvarinen, m, v, a, b, y0, 1, y)
pp_post_par_eta1 <- pp_post_par_conj_normal(1, m, v, a, b, y0, 1, y)
pp_post_par_eta0 <- pp_post_par_conj_normal(0, m, v, a, b, y0, 1, y)

plot_post_mu <- ggplot() +
  geom_function(fun = function(x) {
    dst(x,
        mu=pp_post_par_crps$m_star,
        sigma=pp_post_par_crps$v_star * pp_post_par_crps$b_star / pp_post_par_crps$a_star,
        nu=pp_post_par_crps$a_star * 2)
  },
  aes(color = "CRPS"),
  linewidth = 2.5,
  linetype = 'dotted'
  ) +
  geom_function(fun = function(x) {
    dst(x,
        mu=pp_post_par_hyva$m_star,
        sigma=pp_post_par_hyva$v_star * pp_post_par_hyva$b_star / pp_post_par_hyva$a_star,
        nu=pp_post_par_hyva$a_star * 2)
  },
  aes(color = "Hyva"),
  linewidth = 1
  ) +
  geom_function(fun = function(x) {
    dst(x,
        mu=pp_post_par_eta1$m_star,
        sigma=pp_post_par_eta1$v_star * pp_post_par_eta1$b_star / pp_post_par_eta1$a_star,
        nu=pp_post_par_eta1$a_star * 2)
  },
  aes(color = "eta = 1"),
  linetype = 'dashed',
  linewidth = 1.6
  ) +
  geom_function(fun = function(x) {
    dst(x,
        mu=pp_post_par_eta0$m_star,
        sigma=pp_post_par_eta0$v_star * pp_post_par_eta0$b_star / pp_post_par_eta0$a_star,
        nu=pp_post_par_eta0$a_star * 2)
  },
  aes(color = "eta = 0"),
  linetype = 'dotted',
  linewidth = 1.6
  )+
  geom_vline(aes(xintercept = true_mu_curr[1], color = 'beta_curr'), 
             linewidth = 0.5, linetype = "solid") +
  geom_vline(aes(xintercept = true_mu_hist[1], color = 'beta_hist'), 
             linewidth = 1.2, linetype = "dotted") +
  labs(x = expression(mu), y = '') +
  theme_bw() +
  geom_density(data = data.frame(y = draws_mu_npp_post[,1]),
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
      "beta_curr" = expression(mu['curr']^~'true'),
      "beta_hist" = expression(mu['hist']^~'true')
    ),
    breaks = c("CRPS",
               "Hyva",
               "NPP",
               "eta = 1",
               "eta = 0",
               "beta_hist",
               "beta_curr")
  ) +
  # xlim(c(true_mu_curr[1]-0.5,true_mu_curr[1]+0.5)) +
  theme(text = element_text(size = 12),        # Base text size
        axis.title = element_text(size = 14),  # Axis titles
        axis.text = element_text(size = 12),   # Axis tick labels
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 11),
        legend.position = c(0.2, 0.74),
        legend.background = element_rect(
          fill = "white",     # background color of the legend
          color = "black",    # border color
          size = 0.3,         # border thickness
          linetype = "solid"  # border type
        ),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA))
plot_post_mu

plot_post_eta <- ggplot() +
  geom_density(data = data.frame(y = draws_eta_post), aes(x = y, color = "NPP"), linewidth = 1) +
  geom_vline(aes(xintercept = best_a0_crps, color = 'best_crps'), linewidth = 2.5, linetype = 'dotted') +
  geom_vline(aes(xintercept = best_a0_hyvarinen, color = 'best_hyva'), linewidth = 0.5, linetype = 'solid') +
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

plot_par_post <- plot_post_mu + plot_post_eta + plot_layout(ncol = 2)
plot_par_post
ggsave(plot = plot_par_post,"figures/post_normal.png", width = 15, height = 5, dpi = 320)
