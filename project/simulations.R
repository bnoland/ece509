library(glmnet)
library(MASS)

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
  # TODO
}

# Returns the fractional power sparsity index corresponding to size `p' of the
# parameter vector and parameters `alpha' and `delta'.
fractional_power_sparsity <- function(p, alpha, delta) {
  # TODO
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

#p_values <- c(128, 256, 512)
p_values <- 128
theta_seq <- seq(0.01, 2.4, length.out = 100)
n_trials <- 200
alpha <- 0.4
sigma <- 0.5

sparsity <- linear_sparsity

for (p in p_values) {
  k <- sparsity(p, alpha)
  S <- random_support(p, k)
  beta_star <- true_beta(p, S)
  S_signed <- signed_support(beta_star)
  
  probs <- numeric(length(theta_seq))
  for (i in seq_along(theta_seq)) {
    theta <- theta_seq[[i]]
    n <- sample_size(theta, p, k)
    
    # TODO
    if (n == 1) {
      cat("n = 1 -- skipping")
      next
    }
    
    count <- 0
    for (j in 1:n_trials) {
      X <- random_design(n, p, diag(p))
      y <- X %*% beta_star + rnorm(n)
      
      lambda <- lambda_n(n, p, k, sigma)
      fit <- glmnet(X, y, intercept = FALSE)#, standardize = FALSE)
      beta <- as.numeric(predict(fit, s = lambda, type = "coefficients",
                                 exact = TRUE, x = X, y = y))[-1]
      count <- count + all(signed_support(beta) == S_signed)
    }
    
    probs[[i]] <- count / n_trials
  }
  
  plot(theta_seq, probs)
}
