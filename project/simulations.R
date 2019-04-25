library(glmnet)
library(MASS)
library(ggplot2)
library(tibble)
library(dplyr)
library(latex2exp)

# Returns a random support with sparsity index `k' for a true parameter vector
# of size `p'.
random_support <- function(p, k) {
  S <- c(rep(1, k), rep(0, p - k))
  S <- sample(S)
  S
}

# Returns a vector of length `p' to be used as the true parameter vector. This
# vector is generated using the given support `S'.
true_beta <- function(p, S) {
  beta_min <- 0.5
  
  beta_star <- numeric(p)
  support_indices <- which(S == 1)
  for (i in support_indices) {
    beta_star[[i]] <- sample(c(-beta_min, beta_min), size = 1)
  }
  
  beta_star
}

# Returns the signed support of the vector `beta'.
signed_support <- function(beta) {
  sign(beta)
}

# Returns the linear sparsity index corresponding to size `p' of the parameter
# vector and parameter `alpha'.
linear_sparsity <- function(p, alpha) {
  ceiling(alpha * p)
}

# Returns the sublinear sparsity index corresponding to size `p' of the
# parameter vector and parameter `alpha'.
sublinear_sparsity <- function(p, alpha) {
  ceiling(alpha * p / log(alpha * p))
}

# Returns the fractional power sparsity index corresponding to size `p' of the
# parameter vector and parameters `alpha' and `delta'.
fractional_power_sparsity <- function(p, alpha, delta) {
  ceiling(alpha * p^delta)
}

# Returns a random `n' x `p' design matrix whose columns are generated from a
# zero-mean multivariate normal distribution with covariance matrix
# `cov_matrix'.
random_design <- function(n, p, cov_matrix) {
  X <- mvrnorm(n = n, mu = rep(0, p), Sigma = cov_matrix)
  X <- rbind(X)  # Force X to be a row vector when n = 1
  X
}

# Returns the sample size corresponding to control parameter `theta', `p', and
# sparsity index return `k'.
sample_size <- function(theta, p, k) {
  ceiling(2 * theta * k * log(p - k))
}

# Returns the regularization parameter corresponding to `n', `p', and `k', with
# noise level `sigma'.
lambda_n <- function(n, p, k, sigma) {
  sqrt((2 * sigma^2 * log(k) * log(p - k)) / n)
}

p_values <- c(128, 256, 512)
theta_seq <- seq(0.01, 2.4, length.out = 15)
n_trials <- 20
alpha <- 0.4
delta <- 0.75
sigma <- 0.5
mu <- 0.1

# Center plot titles.
theme_update(plot.title = element_text(hjust = 0.5))

success_prob_plot <- function(p, uniform, mixture = 1, sparsity_fn, ...) {
}

plot_data_list <- list()

for (p in p_values) {
  k <- fractional_power_sparsity(p, alpha, delta)
  #k <- sublinear_sparsity(p, alpha)
  #k <- linear_sparsity(p, alpha)
  
  S <- random_support(p, k)
  beta_star <- true_beta(p, S)
  S_signed <- signed_support(beta_star)
  
  #cov_matrix <- diag(p)
  cov_matrix <- toeplitz(mu^(0:(p-1)))
  
  probs <- numeric(length(theta_seq))
  for (i in seq_along(theta_seq)) {
    theta <- theta_seq[[i]]
    n <- sample_size(theta, p, k)
    
    # TODO: Should we try to avoid this from ever occurring?
    if (n == 1) {
      cat("n = 1 -- skipping")
      next
    }
    
    count <- 0
    for (j in 1:n_trials) {
      X <- random_design(n, p, cov_matrix)
      y <- X %*% beta_star + rnorm(n, sd = sigma)
      
      lambda <- lambda_n(n, p, k, sigma)
      # TODO: In order to speed up coefficient prediction, may want to
      # precompute the sequence of regularization parameters.
      fit <- glmnet(X, y, intercept = FALSE)
      beta <- as.numeric(predict(fit, s = lambda, type = "coefficients",
                                 exact = TRUE, x = X, y = y))[-1]
      count <- count + all(signed_support(beta) == S_signed)
    }
    
    probs[[i]] <- count / n_trials
  }
  
  plot_data_list[[p]] <- tibble(theta_seq, probs, p)
}

plot_data <- bind_rows(plot_data_list)
plot_data$p <- factor(plot_data$p, levels = p_values)
ggplot(plot_data, aes(theta_seq, probs, color = p)) +
  coord_cartesian(xlim = c(0, 2.4), ylim = c(0, 1)) +
  geom_point() +
  geom_line() +
  labs(color = TeX("$p$")) +
  xlab(TeX("Control parameter $\\theta(n, p, k)$")) +
  ylab("Probability of success") +
  theme(aspect.ratio = 1)
