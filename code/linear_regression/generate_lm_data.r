n0 <- 50
n <- 50
beta0 <- c(-0.4, 0.5)
# beta0 <- c(-0.2, 0.1)
beta <- c(-0.4, 0.5)
sigma0 <- 1
sigma <- 1
X0 <- cbind(rep(1,n0), rnorm(n0))
X <- cbind(rep(1,n), rnorm(n))
hist_data <- data.frame(y = X0 %*% beta0 + rnorm(n0, 0, sigma0),
                        X1 = X0[,2])
curr_data <- data.frame(y = X %*% beta + rnorm(n, 0, sigma),
                        X1 = X[,2])
save(hist_data, curr_data, file="data/sim_lm_data.RData")
