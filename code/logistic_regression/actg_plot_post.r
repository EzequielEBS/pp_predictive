library(arrow)
library(ggplot2)
library(dplyr)
library(tidyr)
library(hdbayes)

hist_data <- actg019
curr_data <- actg036
n <- nrow(curr_data)
n0 <- nrow(hist_data)

samples <- read_parquet("samples/logistic_regression/samples_all.parquet")
samples_etas <- read_parquet("samples/logistic_regression/etas_inf_match.parquet")
samples_long <- samples %>%
  pivot_longer(cols = c("(Intercept)", "age", "race", "treatment", "cd4"), names_to = "parameter", values_to = "value")
samples_long <- samples_long %>%
  group_by(parameter) %>%
  filter(value <= quantile(value, 0.995)) %>%  # drop top 0.5% per panel
  ungroup()
samples_long <- samples_long |> filter(method != "small")

pal <- c(
  "big"       = "#D62828",  # red
  "inf_match" = "#F77F00",  # orange
  "inter"     = "#00B4D8",  # sky blue
  "rate"      = "#7B2D8B",  # purple
  "npp"       = "#023E8A",  # dark blue
  "small"     = "#2D6A4F",  # forest green
  "vnpp"      = "#E9C46A"   # gold
)

lty <- c(
  "big"       = "solid",
  "inf_match" = "solid",
  "inter"     = "solid",
  "rate"      = "solid",
  "npp"       = "dashed",
  "small"     = "solid",
  "vnpp"      = "dashed"
)

plot_etas <- ggplot(samples_etas, aes(x = eta)) +
  geom_density(fill = "lightblue") +
  geom_vline(xintercept = median(samples_etas$eta), 
    linetype = "dotted", color = "grey40", linewidth = 0.75) +
  theme_bw() +
  labs(x = expression(hat(eta)), y = "")

plot_beta <- ggplot(samples_long, 
  aes(x = value, color = method, linetype = method)) +
  geom_density(linewidth = 0.85, adjust = 1.2) +
  scale_color_manual(
    values = pal,
    name = NULL,
    labels = c(
      "inter" = expression(eta == 0.5),
      "rate" = bquote(eta == 0.5 * n / n[0] ~ "=" ~ .(round(0.5 * n / n0, 2))),
      # "small" = expression(eta == 0.1),
      "big" = expression(eta == 0.9),
      "inf_match" = bquote(eta[inf_match] == .( round(median(samples_etas$eta), 2) )),
      "npp" = "NPP",
      "vnpp" = "VNPP"
    ),
    breaks = c(
      "npp",
      "vnpp",
      "big",
      "inf_match",
      "inter",
      "rate"
    )
  ) +
  scale_linetype_manual(values = lty, guide = "none") +
  facet_wrap(~ parameter, scales = "free") +
  labs(title = "", x = "", y = "") +
  theme_bw() +
  theme(
    legend.position        = "inside",
    legend.position.inside = c(0.83, 0.25),
    legend.background      = element_rect(fill = "white", color = "grey80"),
    legend.key.width       = unit(1.5, "cm"),
    legend.title           = element_blank(),
    legend.text            = element_text(size = 13),  # <-- increase this
    strip.text             = element_text(size = 11)
  )
plot_beta

samples_or_long <- samples_long %>%
  mutate(value = ifelse(parameter != "(Intercept)", exp(value), value))
samples_or_long <- samples_or_long %>%
  group_by(parameter) %>%
  filter(value <= quantile(value, 0.995)) %>%  # drop top 0.5% per panel
  ungroup()

param_labels <- c(
  "(Intercept)" = "(Intercept)",
  "age"         = "OR(age)",
  "cd4"         = "OR(cd4)",
  "race"        = "OR(race)",
  "treatment"   = "OR(treatment)"
)

vline_df <- data.frame(parameter = unique(samples_or_long$parameter)) %>%
  filter(parameter != "(Intercept)")

plot_or <- ggplot(samples_or_long,  
  aes(x = value, color = method, linetype = method)) +
  geom_density(linewidth = 0.85, adjust = 1.2) +
  scale_color_manual(
    values = pal,
    name = NULL,
    labels = c(
      "inter" = expression(eta == 0.5),
      "rate" = bquote(eta == 0.5 * n / n[0] ~ "=" ~ .(round(0.5 * n / n0, 2))),
      # "small" = expression(eta == 0.1),
      "big" = expression(eta == 0.9),
      "inf_match" = bquote(eta[inf_match] == .( round(median(samples_etas$eta), 2) )),
      "npp" = "NPP",
      "vnpp" = "VNPP"
    ),
    breaks = c(
      "npp",
      "vnpp",
      "big",
      "inf_match",
      "inter",
      "rate"
    )
  ) +
  scale_linetype_manual(values = lty, guide = "none") +
  geom_vline(
    data        = vline_df,
    aes(xintercept = 1),
    linetype    = "dotted",
    color       = "grey40",
    inherit.aes = FALSE,
    linewidth   = 0.75
  ) +
  facet_wrap(~ parameter, scales = "free",
              labeller = labeller(parameter = param_labels)
  ) +
  labs(title = "", x = "", y = "") +
  theme_bw() +
  theme(
    legend.position        = "inside",
    legend.position.inside = c(0.83, 0.25),
    legend.background      = element_rect(fill = "white", color = "grey80"),
    legend.key.width       = unit(1.5, "cm"),
    legend.title           = element_blank(),
    legend.text            = element_text(size = 13),  # <-- increase this
    strip.text             = element_text(size = 11)
  )
plot_or

samples_eta_npp <- samples %>%
  filter(method %in% c("npp", "vnpp")) %>%
  select(method, starts_with("eta")) %>%
  pivot_longer(cols = starts_with("eta"), names_to = "parameter", values_to = "value")

plot_eta <- ggplot(samples_eta_npp, aes(x = value, color = method, linetype = method)) +
  geom_density(linewidth = 0.85, adjust = 1.2, bounds = c(0, 1)) +
  geom_vline(xintercept = median(samples_etas$eta),
    linetype = "dotted", color = "grey40", linewidth = 0.75) +
  scale_color_manual(
    values = pal,
    name = NULL,
    labels = c(
      "npp" = "NPP",
      "vnpp" = "VNPP"
    ),
    breaks = c("npp", "vnpp")
  ) +
  scale_linetype_manual(values = lty, guide = "none") +
  facet_wrap(~ parameter, scales = "free") +
  labs(title = "", x = expression(eta), y = "") +
  theme_bw() +
  theme(
    legend.position        = "inside",
    legend.position.inside = c(0.25, 0.25),
    legend.background      = element_rect(fill = "white", color = "grey80"),
    legend.key.width       = unit(1.5, "cm"),
    legend.title           = element_blank(),
    legend.text            = element_text(size = 11),  # <-- increase this
    strip.text             = element_text(size = 11)
  )
plot_eta

# save figures
ggsave("figures/actg_beta_comp.png", plot_beta, width = 10, height = 6)
ggsave("figures/actg_or_comp.png", plot_or, width = 10, height = 6)
ggsave("figures/actg_eta_comp.png", plot_eta, width = 4, height = 4)
ggsave("figures/actg_eta_density.png", plot_etas, width = 4, height = 4)


