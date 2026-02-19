library(mvtnorm)
library(parallel)
library(hdbayes)
library(cmdstanr)
library(dplyr)

# load data
# data(airquality)
# hist_data <- airquality %>%
#   filter(Month <= 6)
# curr_data <- airquality %>%
#   filter(Month > 6)
# 
# # remove rows with missing values
# hist_data <- na.omit(hist_data)
# curr_data <- na.omit(curr_data)
# 
# # create response variable
# hist_data$logWind <- log(hist_data$Wind)
# curr_data$logWind <- log(curr_data$Wind)
# 
# # normalize predictor variables
# hist_data$Ozone <- (hist_data$Ozone - mean(hist_data$Ozone)) / sd(hist_data$Ozone)
# curr_data$Ozone <- (curr_data$Ozone - mean(curr_data$Ozone)) / sd(curr_data$Ozone)
# hist_data$Solar.R <- (hist_data$Solar.R - mean(hist_data$Solar.R)) / sd(hist_data$Solar.R)
# curr_data$Solar.R <- (curr_data$Solar.R - mean(curr_data$Solar.R)) / sd(curr_data$Solar.R)
# hist_data$Temp <- (hist_data$Temp - mean(hist_data$Temp)) / sd(hist_data$Temp)
# curr_data$Temp <- (curr_data$Temp - mean(curr_data$Temp)) / sd(curr_data$Temp)
# 
# # set parameters for Stan data
# formula     <- logWind ~ Ozone + Solar.R + Temp + Day
# family      <- gaussian()


load("data/sim_lm_data.RData")
formula <- y ~ X1
family <- gaussian()

# ess_hist <- function(M, N, X0, X, beta, tau, c) {
#   n0 <- nrow(X0)
#   n <- nrow(X)
#   D.m <- mclapply(0:M, function(m) {
#     Dpplus <- sum(tau*diag(t(X0) %*% X0)) + n0/(2*tau^2)
#     Dq0 <- tau / c^2 * diag(t(X0) %*% X0) + n0/(2*tau^2)
#     
#     Dqj <- lapply(1:N, function(j) {
#       # Xmt <- cbind(rep(1, m),matrix(rnorm(m*(p-1)), nrow = m, ncol = p))
#       id <- sample(n, m, replace = TRUE)
#       Xmt <- X[id,]
#       # Ymt <- rmvnorm(1, Xm %*% beta, 1/tau*diag(1, m))
#       # betat <- rmvnorm(1, mu_eta, 1/tau * inv_S_eta)
#       tau*diag(t(Xmt) %*% Xmt)
#     })
#     Dq <- rowMeans(do.call(cbind, Dqj)) + m/(2*tau^2) + Dq0
#     
#     Dqplus <- sum(Dq)
#     deltam <- abs(Dpplus - Dqplus)
#   },
#   mc.cores = 14)
#   D.m <- unlist(D.m)
#   
#   D.min.n <- which(abs(D.m) == min(abs(D.m)))
#   D.min.v <- D.m[which(abs(D.m) == min(abs(D.m)))]
#   {
#     if (D.min.v < 0)       {
#       D.min.v.nxt <- D.m[D.min.n+1]
#       pESS <- D.min.n - 1 + (-D.min.v / (-D.min.v + D.min.v.nxt))
#     }
#     else if (D.min.v > 0)  {
#       D.min.v.prv <- D.m[D.min.n-1]
#       pESS <- D.min.n - 1 - (D.min.v / (D.min.v - D.min.v.prv))
#     }
#     else if (D.min.v == 0) {
#       pESS <- D.min.n -1
#     }
#   }
#   
#   return(pESS)
# }
# 
# ess_pi0 <- function(M, N, X0, X, beta, tau, mu_beta, S_beta, a, c) {
#   p <- ncol(X0)
#   n <- nrow(X)
#   invS_beta <- solve(S_beta)
#   D.m <- mclapply(0:M, function(m) {
#     Dpplus <- sum(tau*diag(invS_beta)) + p/(2*tau^2) + (a-1)/(tau^2)
#     Dq0 <- tau / c^2 * diag(invS_beta) + p/(2*tau^2) + (a/c-1)/(tau^2)
#     
#     Dqj <- lapply(1:N, function(j) {
#       # Xmt <- cbind(rep(1, m),matrix(rnorm(m*(p-1)), nrow = m, ncol = p))
#       id <- sample(n, m, replace = TRUE)
#       Xmt <- X[id,]
#       # Ymt <- rmvnorm(1, Xm %*% beta, 1/tau*diag(1, m))
#       # betat <- rmvnorm(1, mu_eta, 1/tau * inv_S_eta)
#       tau*diag(t(Xmt) %*% Xmt)
#     })
#     Dq <- rowMeans(do.call(cbind, Dqj)) + m/(2*tau^2) + Dq0
#     
#     Dqplus <- sum(Dq)
#     deltam <- abs(Dpplus - Dqplus)
#   },
#   mc.cores = 14)
#   D.m <- unlist(D.m)
#   
#   D.min.n <- which(abs(D.m) == min(abs(D.m)))
#   D.min.v <- D.m[which(abs(D.m) == min(abs(D.m)))]
#   {
#     if (D.min.v < 0)       {
#       D.min.v.nxt <- D.m[D.min.n+1]
#       pESS <- D.min.n - 1 + (-D.min.v / (-D.min.v + D.min.v.nxt))
#     }
#     else if (D.min.v > 0)  {
#       if (D.min.n - 1 == 0) {
#         pESS <- D.min.n - 1
#       } else {
#         D.min.v.prv <- D.m[D.min.n-1]
#         pESS <- D.min.n - 1 - (D.min.v / (D.min.v - D.min.v.prv))
#       }
#     }
#     else if (D.min.v == 0) {
#       pESS <- D.min.n -1
#     }
#   }
#   return(pESS)
# }


ess_lm_normgamma_pp <- function(M, N, X, beta, tau, S_eta, c, a) {
  p <- length(beta)
  D.m <- mclapply(0:M, function(m) {
    Dpplus <- sum(tau*diag(S_eta)) + p/(2*tau^2) + (a-1)/(tau^2)
    Dq0 <- tau/(c^2)*diag(S_eta) + p/(2*tau^2) + (a/c - 1)/(tau^2)
    Dqj <- lapply(1:N, function(j) {
      # Xmt <- cbind(rep(1, m),matrix(rnorm(m*(p-1)), nrow = m, ncol = p))
      id <- sample(n, m, replace = TRUE)
      Xmt <- X[id,]
      # Ymt <- rmvnorm(1, Xm %*% beta, 1/tau*diag(1, m))
      # betat <- rmvnorm(1, mu_eta, 1/tau * inv_S_eta)
      tau*diag(t(Xmt) %*% Xmt)
    })
    Dq <- rowMeans(do.call(cbind, Dqj)) + m/(2*tau^2) + Dq0

    Dqplus <- sum(Dq)
    deltam <- Dpplus - Dqplus
  },
  mc.cores = 14)
  D.m <- unlist(D.m)

  D.min.n <- which(abs(D.m) == min(abs(D.m)))
  D.min.v <- D.m[which(abs(D.m) == min(abs(D.m)))]
  {
    if (D.min.v < 0)       {
      D.min.v.nxt <- D.m[D.min.n+1]
      pESS <- D.min.n - 1 + (-D.min.v / (-D.min.v + D.min.v.nxt))
    }
    else if (D.min.v > 0)  {
      if (D.min.n - 1 == 0) {
        pESS <- D.min.n - 1
      } else {
        D.min.v.prv <- D.m[D.min.n-1]
        pESS <- D.min.n - 1 - (D.min.v / (D.min.v - D.min.v.prv))
      }
    }
    else if (D.min.v == 0) {
      pESS <- D.min.n -1
    }
  }
  if(pESS < 0) {
    pESS <- 0
  }
  return(pESS)
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
                                 c) {
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
      # Xmt <- cbind(rep(1, m),matrix(rnorm(m*(p-1)), nrow = m, ncol = p))
      id <- sample(n, m, replace = TRUE)
      Xmt <- X[id,]
      # Ymt <- rmvnorm(1, Xm %*% beta, 1/tau*diag(1, m))
      # betat <- rmvnorm(1, mu_eta, 1/tau * inv_S_eta)
      tau*diag(t(Xmt) %*% Xmt)
    })
    Dq <- rowMeans(do.call(cbind, Dqj)) + m/(2*tau^2) + Dq0
    
    Dqplus <- sum(Dq)
    deltam <- Dpplus - Dqplus
  },
  mc.cores = 14)
  D.m <- unlist(D.m)
  
  D.min.n <- which(abs(D.m) == min(abs(D.m)))
  D.min.v <- D.m[which(abs(D.m) == min(abs(D.m)))]
  {
    if (D.min.v < 0)       {
      D.min.v.nxt <- D.m[D.min.n+1]
      pESS <- D.min.n - 1 + (-D.min.v / (-D.min.v + D.min.v.nxt))
    }
    else if (D.min.v > 0)  {
      if (D.min.n - 1 == 0) {
        pESS <- D.min.n - 1
      } else {
        D.min.v.prv <- D.m[D.min.n-1]
        pESS <- D.min.n - 1 - (D.min.v / (D.min.v - D.min.v.prv))
      }
    }
    else if (D.min.v == 0) {
      pESS <- D.min.n -1
    }
  }
  if(pESS < 0) {
    pESS <- 0
  }
  return(pESS)
}

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
  esss <- ess_lm_normgamma_pp(M, N, X, beta, tau, S_eta, c, a)
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
                     c)

save(ess_eta, ess_npp, file = "samples_ppc/ess.RData")
