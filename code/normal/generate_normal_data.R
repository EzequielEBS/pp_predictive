set.seed(20260310)
n0 <- 100
n <- 100
mu0 <- 1
mu <- 1
sigma0 <- 1
sigma <- 1
hist_data <- data.frame(y = rnorm(n0, mu0, sigma0))
curr_data <- data.frame(y = rnorm(n, mu, sigma))
n <- nrow(curr_data)
train_idx <- sample(seq_len(n), size = 0.7 * n)
curr_train <- curr_data[train_idx, , drop = FALSE]
curr_test <- curr_data[-train_idx, , drop = FALSE]

save(hist_data, curr_train, curr_test, file="data/sim_normal_data.RData")
