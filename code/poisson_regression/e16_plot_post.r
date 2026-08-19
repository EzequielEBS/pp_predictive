library(arrow)
library(ggplot2)
library(dplyr)
library(tidyr)
library(hdbayes)

hist_data <- E1684
curr_data <- E1690
n <- nrow(curr_data)
n0 <- nrow(hist_data)

samples <- read_parquet("samples/poisson_regression/samples_all.parquet")
samples_etas <- read_parquet("samples/poisson_regression/etas_inf_match.parquet")
samples_long <- samples %>%
  pivot_longer(cols = c(
    "(Intercept)", 
    "interval(0.26, 0.54]",
    "interval(0.54, 0.89]",
    "interval(0.89, 1.63]",
    "interval(1.63, inf]",  
    "treatment",         
    "sex",                  
    "cage",                 
    "node_bin"            
    ), names_to = "parameter", values_to = "value")
samples_long <- samples_long %>%
  group_by(parameter) %>%
  filter(value <= quantile(value, 0.995)) %>%  # drop top 0.5% per panel
  ungroup()

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
      "small" = expression(eta == 0.1),
      "big" = expression(eta == 0.9),
      "inf_match" = bquote(eta[inf_match] == .( round(median(samples_etas$eta), 2) )),
      "npp" = "NPP",
      "vnpp" = "VNPP"
    ),
    breaks = c(
      "npp",
      "vnpp",
      "inf_match",
      "big",
      "rate",
      "inter",
      "small"
    )
  ) +
  scale_linetype_manual(values = lty, guide = "none") +
  facet_wrap(~ parameter, scales = "free") +
  labs(title = "", x = "", y = "") +
  theme_bw() +
  theme(
    legend.position        = "bottom",
    legend.position.inside = c(0.83, 0.25),
    legend.background      = element_rect(fill = "white", color = "grey80"),
    legend.key.width       = unit(1.5, "cm"),
    legend.title           = element_blank(),
    legend.text            = element_text(size = 13),  # <-- increase this
    strip.text             = element_text(size = 11)
  )
plot_beta


samples_eta_npp <- samples %>%
  filter(method %in% c("npp", "vnpp")) %>%
  select(method, starts_with("eta")) %>%
  pivot_longer(cols = starts_with("eta"), names_to = "parameter", values_to = "value")

plot_eta_npp <- ggplot(samples_eta_npp, aes(x = value, color = method, linetype = method)) +
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
ggsave("figures/e16_beta_comp.png", plot_beta, width = 10, height = 6)
ggsave("figures/e16_eta_comp.png", plot_eta_npp, width = 4, height = 4)
ggsave("figures/e16_eta_density.png", plot_etas, width = 4, height = 4)
