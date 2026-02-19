n0 <- 50
n <- 50
mu0 <- 0.5
mu <- 0.5
sigma0 <- 1
sigma <- 1
hist_data <- data.frame(y = rnorm(n0, mu0, sigma0))
curr_data <- data.frame(y = rnorm(n, mu, sigma))
save(hist_data, curr_data, file="data/sim_normal_data.RData")
