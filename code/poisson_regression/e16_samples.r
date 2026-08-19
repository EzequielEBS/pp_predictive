library(tidyverse)
library(survival)
library(hdbayes)
library(posterior)
library(arrow)


source("code/aux_fun_inf_match_glm.r")

hist <- E1684
curr <- E1690

## replace 0 failure times with 0.50 days
hist <- hist %>% mutate(failtime = if_else(failtime == 0, 0.50/365.25, failtime)) # 1 subject w/ failtime = 0
curr <- curr %>% mutate(failtime = if_else(failtime == 0, 0.50/365.25, failtime)) # 10 subjects w/ failtime = 0

## Center and scale age
hist$cage <- as.numeric( scale(hist$age, center = T, scale = T) )
curr$cage <- as.numeric( scale(curr$age, center = T, scale = T) )

estimate_eta_e16 <- function(
  formula,
  hist_data,
  curr_data,
  hist_raw,
  beta_draws, 
  family = poisson(), 
  nu = .5,
  offset.list = NULL
) {
  # Estimate beta_hat and beta0_hat using IRLS
  formula <- as.formula(formula)
  outcome_var <- all.vars(formula)[1]
  y <- curr_data[[outcome_var]]
  X <- model.matrix(formula, data = curr_data)
  y0 <- hist_data[[outcome_var]]
  X0 <- model.matrix(formula, data = hist_data)
  # beta_hat <- irls(y, X, family, off = offset.list[[1]])
  # beta0_hat <- irls(y0, X0, family, off = offset.list[[2]])
  if (!is.null(offset.list)) {
    beta_hat <- 
      glm(formula, data = curr_data, family = family, offset = offset.list[[1]])$coefficients %>%
      as.matrix(ncol = 1)
    beta0_hat <-
      glm(formula, data = hist_data, family = family, offset = offset.list[[2]])$coefficients %>%
      as.matrix(ncol = 1)
  } else {
    beta_hat <- glm(formula, data = curr_data, family = family)$coefficients %>%
      as.matrix(ncol = 1)
    beta0_hat <- glm(formula, data = hist_data, family = family)$coefficients %>%
      as.matrix(ncol = 1)
  }

  # print(sum((beta_hat - beta0_hat)^2))
  # print(beta_hat)
  # print(beta0_hat)

  k0 <- floor(nrow(hist_raw)^(nu))
  idx0 <- sample(nrow(hist_raw), size = k0, replace = FALSE)
  hist_data_k0 <- hist_raw[idx0, ]

  nbreaks <- 5
  probs   <- 1:nbreaks / nbreaks
  breaks  <- curr |> 
    filter(failcens == 1) |> 
    reframe(quant = quantile(failtime, probs = probs)) |> 
    unlist()

  breaks <- as.numeric(breaks[-nbreaks])

  split_fmla  <- Surv(failtime, failcens) ~ treatment + sex + cage + node_bin
  hist_pseudo_k0 <- survival::survSplit(formula = split_fmla, data = hist_data_k0, cut = breaks,
                                    episode = "interval", start = "start")

  breaks_labs <- paste("(", c(0, round(breaks, 2)), ", ", c(round(breaks, 2), "inf"), "]", sep="")
  levels_interval <- factor(seq_len(nbreaks), labels = breaks_labs)

  hist_pseudo_k0 <- hist_pseudo_k0 |>  
    mutate(exposure    = failtime - start,
          log_exposure = log(exposure),
          interval    = factor(interval, levels = seq_along(levels_interval), labels = breaks_labs))
  X0k0 <- model.matrix(formula, data = hist_pseudo_k0)

  if (!is.null(offset.list)) {
    off <- hist_pseudo_k0$log_exposure
  } else {
    off <- NULL
  }
  num <- expected_sq_norm_score_vec(X0k0, beta_draws, beta_hat, family, off = off)
  den <- expected_sq_norm_score_vec(X0k0, beta_draws, beta0_hat, family, off = off)

  hat_eta <- exp(0.5 * (log(mean(num)) - log(mean(den))))
  return(hat_eta)
}



nbreaks <- 5
probs   <- 1:nbreaks / nbreaks
breaks  <- curr |> 
  filter(failcens == 1) |> 
  reframe(quant = quantile(failtime, probs = probs)) |> 
  unlist()

breaks <- as.numeric(breaks[-nbreaks])

split_fmla  <- Surv(failtime, failcens) ~ treatment + sex + cage + node_bin
hist_pseudo <- survival::survSplit(formula = split_fmla, data = hist, cut = breaks,
                                   episode = "interval", start = "start")
curr_pseudo <- survival::survSplit(formula = split_fmla, data = curr, cut = breaks,
                                   episode = "interval", start = "start")

breaks_labs <- paste("(", c(0, round(breaks, 2)), ", ", c(round(breaks, 2), "inf"), "]", sep="")
levels_interval <- factor(seq_len(nbreaks), labels = breaks_labs)
ninterval   <- length(unique(curr_pseudo$interval))

hist_pseudo <- hist_pseudo |>  
  mutate(exposure    = failtime - start,
         log_exposure = log(exposure),
         interval    = factor(interval, levels = seq_along(levels_interval), labels = breaks_labs))

curr_pseudo <- curr_pseudo |>  
  mutate(exposure    = failtime - start,
         log_exposure = log(exposure),
         interval    = factor(interval, levels = seq_along(levels_interval), labels = breaks_labs))

fit.curr <- glm(failcens ~ interval + treatment + sex + cage + node_bin + offset(log_exposure),
                data = curr_pseudo, family = poisson(link = "log"))
fit.hist <- glm(failcens ~ interval + treatment + sex + cage + node_bin + offset(log_exposure),
                data = hist_pseudo, family = poisson(link = "log"))

data.list   <- list(curr_pseudo, hist_pseudo)
offset.list <- list(curr_pseudo$log_exposure, hist_pseudo$log_exposure)
fmla   <- failcens ~ interval + treatment + sex + cage + node_bin
family <- poisson(link = "log")
pars   <- names(coefficients(fit.curr))

a0 <- seq(0, 1, length.out = 20)

logncfun <- function(a0, ...){
  hdbayes::glm.npp.lognc(
    formula = fmla, family = family, histdata = hist_pseudo, a0 = a0, ...
  )
}

a0.lognc = lapply(
  X = a0, 
  FUN = logncfun,
  iter_warmup = 1000,
  iter_sampling = 2000,
  chains = 4,
  parallel_chains = 4,
  refresh = 0,
  show_messages = FALSE,
  show_exceptions = FALSE
)

a0.lognc <- data.frame( do.call(rbind, a0.lognc) )

fit_npp <- glm.npp(
  formula = fmla,
  family = family,
  data.list = data.list,
  offset.list = offset.list,
  a0.lognc = a0.lognc$a0,
  lognc = matrix(a0.lognc$lognc, ncol = 1),
  iter_warmup = 5000,
  iter_sampling = 10000,
  chains = 4,
  parallel_chains = 4,
  refresh = 0,
)

eta_inter <- 1/2
eta_rate <- 1/2 * nrow(data.list[[1]]) / nrow(data.list[[2]])
eta_small <- 0.1
eta_big <- 0.9
 

fit0 = glm.pp(
  formula = fmla,
  family = family,
  data.list = data.list,
  offset.list = offset.list,
  a0.vals = 0,
  iter_warmup = 5000,
  iter_sampling = 10000,
  chains = 4,
  parallel_chains = 1,
  refresh = 0,
)
beta_draws <- fit0[, -1] %>% as_draws_matrix()

etas <- lapply(1:100, function(i) {
  estimate_eta_e16(
    formula = fmla,
    curr_data = data.list[[1]],
    hist_data = data.list[[2]],
    beta_draws = beta_draws,
    family = family,
    offset.list = offset.list,
    hist_raw = hist,
  )
}) %>%
  unlist()
etas_df <- data.frame(eta = etas)

eta_inf_match <- median(etas)
eta_inf_match

fit_eta_inter <- glm.pp(
  formula = fmla,
  family = family,
  data.list = data.list,
  offset.list = offset.list,
  a0.vals = eta_inter,
  iter_warmup = 5000,
  iter_sampling = 10000,
  chains = 4,
  parallel_chains = 4,
  refresh = 0,
)
fit_eta_rate <- glm.pp(
  formula = fmla,
  family = family,
  data.list = data.list,
  offset.list = offset.list,
  a0.vals = eta_rate,
  iter_warmup = 5000,
  iter_sampling = 10000,
  chains = 4,
  parallel_chains = 4,
  refresh = 0,
)
fit_eta_small <- glm.pp(
  formula = fmla,
  family = family,
  data.list = data.list,
  offset.list = offset.list,
  a0.vals = eta_small,
  iter_warmup = 5000,
  iter_sampling = 10000,
  chains = 4,
  parallel_chains = 4,
  refresh = 0,
)
fit_eta_big <- glm.pp(
  formula = fmla,
  family = family,
  data.list = data.list,
  offset.list = offset.list,
  a0.vals = eta_big,
  iter_warmup = 5000,
  iter_sampling = 10000,
  chains = 4,
  parallel_chains = 4,
  refresh = 0,
)
fit_eta_inf_match <- glm.pp(
  formula = fmla,
  family = family,
  data.list = data.list,
  offset.list = offset.list,
  a0.vals = ifelse(eta_inf_match > 1, 1, eta_inf_match),
  iter_warmup = 5000,
  iter_sampling = 10000,
  chains = 4,
  parallel_chains = 4,
  refresh = 0,
)

fit_vnpp <- read_parquet("samples/poisson_regression/samples_variational.parquet")

fit_eta_inter <- fit_eta_inter %>% select(-.chain, -.iteration, -.draw, -lp__)
fit_eta_inter$eta <- eta_inter
fit_eta_inter$method <- "inter"
fit_eta_rate <- fit_eta_rate %>% select(-.chain, -.iteration, -.draw, -lp__)
fit_eta_rate$eta <- eta_rate
fit_eta_rate$method <- "rate"
fit_eta_small <- fit_eta_small %>% select(-.chain, -.iteration, -.draw, -lp__)
fit_eta_small$eta <- eta_small
fit_eta_small$method <- "small"
fit_eta_big <- fit_eta_big %>% select(-.chain, -.iteration, -.draw, -lp__)
fit_eta_big$eta <- eta_big
fit_eta_big$method <- "big"
fit_eta_inf_match <- fit_eta_inf_match %>% select(-.chain, -.iteration, -.draw, -lp__)
fit_eta_inf_match$eta <- eta_inf_match
fit_eta_inf_match$method <- "inf_match"
fit_npp <- fit_npp %>% select(-.chain, -.iteration, -.draw, -lp__, -`logit_a0s[1]`)
fit_npp$method <- "npp"
fit_npp <- fit_npp %>% rename(eta = a0_hist_1)
fit_vnpp$method <- "vnpp"

all_fits <- bind_rows(
  fit_eta_inter, fit_eta_rate, fit_eta_small, fit_eta_big, fit_eta_inf_match,
  fit_npp, fit_vnpp
)
write_parquet(all_fits, "samples/poisson_regression/samples_all.parquet")
write_parquet(etas_df, "samples/poisson_regression/etas_inf_match.parquet")
