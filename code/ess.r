library(mvtnorm)
library(parallel)
library(hdbayes)
library(cmdstanr)
library(dplyr)

load("data/sim_lm_data.RData")
formula <- y ~ X1
family <- gaussian()

# Finds the pseudo effective-sample-size (pESS) implied by a sequence D.m of
# "prior mass minus running-average mass" differences, one entry per
# candidate sample size m = 0, 1, ..., M. pESS is the (linearly interpolated)
# crossing point where D.m changes sign.
#
# This root-finding logic used to be copy-pasted (with a subtly different,
# tie-unsafe version) inside both ess_lm_normgamma_pp() and
# ess_lm_normgamma_npp(); it's factored out here so there is exactly one
# implementation to get right.
find_pess <- function(D.m) {
  candidates <- which(abs(D.m) == min(abs(D.m)))
  if (length(candidates) > 1) {
    warning("Multiple candidate roots found for pESS; using the first one.")
  }
  D.min.n <- candidates[1]
  D.min.v <- D.m[D.min.n]

  if (D.min.v < 0) {
    D.min.v.nxt <- D.m[D.min.n + 1]
    pESS <- D.min.n - 1 + (-D.min.v / (-D.min.v + D.min.v.nxt))
  } else if (D.min.v > 0) {
    if (D.min.n - 1 == 0) {
      pESS <- D.min.n - 1
    } else {
      D.min.v.prv <- D.m[D.min.n - 1]
      pESS <- D.min.n - 1 - (D.min.v / (D.min.v - D.min.v.prv))
    }
  } else {
    pESS <- D.min.n - 1
  }

  max(pESS, 0)
}

ess_lm_normgamma_pp <- function(M, N, X, beta, tau, S_eta, c, a, n,
                                 mc_cores = max(1, detectCores() - 1)) {
  p <- length(beta)
  D.m <- mclapply(0:M, function(m) {
    Dpplus <- sum(tau*diag(S_eta)) + p/(2*tau^2) + (a-1)/(tau^2)
    Dq0 <- tau/(c^2)*diag(S_eta) + p/(2*tau^2) + (a/c - 1)/(tau^2)
    Dqj <- lapply(1:N, function(j) {
      id <- sample(n, m, replace = TRUE)
      Xmt <- X[id,]
      tau*diag(t(Xmt) %*% Xmt)
    })
    Dq <- rowMeans(do.call(cbind, Dqj)) + m/(2*tau^2) + Dq0

    Dqplus <- sum(Dq)
    Dpplus - Dqplus
  },
  mc.cores = mc_cores)
  D.m <- unlist(D.m)

  find_pess(D.m)
}

ess_lm_normgamma_npp <- function(M,
                                 N,
                                 X0,
                                 y0,
                                 X,
                                 beta,
                                 tau,
                                 eta,
                                 mu_beta,
                                 S_beta,
                                 a,
                                 b,
                                 a1,
                                 b1,
                                 c,
                                 n,
                                 mc_cores = max(1, detectCores() - 1)) {
  p <- length(beta)
  n0 <- length(y0)
  invS_beta <- solve(S_beta)
  S_eta <- invS_beta + eta * t(X0) %*% X0
  invS_eta <- solve(S_eta)
  mu_eta <- invS_eta %*% (invS_beta %*% mu_beta + eta * t(X0) %*% y0)

  draws_peta <- lapply(1:10000, function(i){
    draw_tau <- rgamma(1, a, b)
    draw_beta <- rmvnorm(1, mean = mu_eta, sigma = 1/tau*invS_eta)
    return(c(draw_beta, draw_tau))
  })

  l0 <- lapply(draws_peta, function(draw) {
    draw_beta <- draw[1:p]
    draw_tau <- draw[p+1]
    dmvnorm(y0, mean = X0 %*% draw_beta, sigma = 1/tau*diag(1,n0), log = TRUE)
  })
  l0 <- unlist(l0)
  mean2_l0 <- mean(l0)^2

  l02 <- lapply(draws_peta, function(draw) {
    draw_beta <- draw[1:p]
    draw_tau <- draw[p+1]
    dmvnorm(y0, mean = X0 %*% draw_beta, sigma = 1/tau*diag(1,n0), log = TRUE)^2
  })
  l02 <- unlist(l02)
  mean_l02 <- mean(l02)

  D.m <- mclapply(0:M, function(m) {
    Dpplus <- sum(tau*diag(S_eta)) + p/(2*tau^2) + (a-1)/(tau^2) -
      mean2_l0 + mean_l02 +
      - (a1-1)/eta + (b1/c-1)/(1-eta)
    Dq0 <- tau/(c^2)*diag(S_eta) + p/(2*tau^2) + (a/c - 1)/(tau^2) -
      mean2_l0 + mean_l02 +
      - (a1/c-1)/eta + (b1/c-1)/(1-eta)
    Dqj <- lapply(1:N, function(j) {
      id <- sample(n, m, replace = TRUE)
      Xmt <- X[id,]
      tau*diag(t(Xmt) %*% Xmt)
    })
    Dq <- rowMeans(do.call(cbind, Dqj)) + m/(2*tau^2) + Dq0

    Dqplus <- sum(Dq)
    Dpplus - Dqplus
  },
  mc.cores = mc_cores)
  D.m <- unlist(D.m)

  find_pess(D.m)
}

set.seed(20260819)

res_hist          = hdbayes:::stack.data(formula = formula, data.list = list(hist_data))
res_curr          = hdbayes:::stack.data(formula = formula, data.list = list(curr_data))
y0            = res_hist$y
X0            = res_hist$X
y = res_curr$y
X = res_curr$X
n <- length(y)
p <- ncol(X)
a <- 2
b <- 1
mu_beta <- rep(0, p)
S_beta <- diag(p)
invS_beta <- solve(S_beta)

M <- n
N <- 100000
tau <- a/b
c <- 100

a0_list      <- seq(0, 1, length.out = 40)
ess_eta <- lapply(a0_list, function(eta) {
  S_eta <- invS_beta + eta * t(X0) %*% X0
  invS_eta <- solve(S_eta)
  mu_eta <- invS_eta %*% (invS_beta %*% mu_beta + eta * t(X0) %*% y0)
  beta <- mu_eta
  esss <- ess_lm_normgamma_pp(M, N, X, beta, tau, S_eta, c, a, n = n)
  return(esss)
})

plot(a0_list, ess_eta, type = 'b', xlab = expression(a[0]), ylab = 'ESS')

a1 <- 1
b1 <- 1
eta <- a1/(a1 + b1)
ess_npp <- ess_lm_normgamma_npp(M,
                     N,
                     X0,
                     y0,
                     X,
                     beta,
                     tau,
                     eta,
                     mu_beta,
                     S_beta,
                     a,
                     b,
                     a1,
                     b1,
                     c,
                     n = n)

save(ess_eta, ess_npp, file = "samples_ppc/ess.RData")
