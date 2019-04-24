library(glmnet)

p_values <- c(128, 256, 512)

# Returns a random support for the true parameter vector for a given size `p' of
# the parameter vector.
random_support <- function(p) {
}

# Returns a vector of length `p' to be used as the true parameter vector. This
# vector is generated using the given support `S'.
true_beta <- function(p, S) {
}

# Returns the signed support of the vector `beta'.
signed_support <- function(beta) {
}

# Returns the linear sparsity index corresponding to size `p' of the parameter
# vector and parameter `alpha'.
linear_sparsity <- function(p, alpha) {
}

# Returns the sublinear sparsity index corresponding to size `p' of the
# parameter vector and parameter `alpha'.
sublinear_sparsity <- function(p, alpha) {
}

# Returns the fractional power sparsity index corresponding to size `p' of the
# parameter vector and parameters `alpha' and `delta'.
fractional_power_sparsity <- function(p, alpha, delta) {
}

# Returns a random `n' x `p' design matrix whose columns are generated from a
# zero-mean multivariate normal distribution with covariance matrix
# `cov_matrix'.
random_design <- function(n, p, cov_matrix) {
}

# Returns the sample size corresponding to control parameter `theta', `p', and
# `sparsity index return `k'.
sample_size <- function(theta, p, k) {
}

# Returns the regularization parameter corresponding to `n', `p', and `k', with
# noise level `sigma'.
lambda_n <- function(n, p, k, sigma) {
}

theta_seq <- seq(0, 2.4, length.out = 10)
n_trials <- 200
alpha <- 0.4
sigma <- 0.5

sparsity <- linear_sparsity

for (p in p_values) {
  S <- random_support(p)
  beta_star <- true_beta(p, S)
  S_signed <- signed_support(beta_star)
  
  probs <- numeric(length(theta_seq))
  for (i in seq_along(theta_seq)) {
    theta <- theta_seq[[i]]
    k <- sparsity(p, alpha)
    n <- sample_size(theta, p, k)
    
    count <- 0
    for (j in 1:n_trials) {
      X <- random_design(n, p, diag(p))
      #y <- X %*% beta_star + rnorm(n)
      
      lambda <- lambda_n(n, p, k, sigma)
      #fit <- glmnet(X, y, intercept = FALSE)
      # beta <- as.numeric(predict(fit, s = lambda, type = "coefficients",
      #                            exact = TRUE, x = X, y = y))
      count <- count + (signed_support(beta) == S_signed)
    }
    
    probs[[i]] <- count / length(probs)
  }
  
  plot(theta_seq, probs)
}
