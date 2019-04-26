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
# vector and parameter `gamma'.
linear_sparsity <- function(p, gamma) {
  ceiling(gamma * p)
}

# Returns the sublinear sparsity index corresponding to size `p' of the
# parameter vector and parameter `gamma'.
sublinear_sparsity <- function(p, gamma) {
  ceiling(gamma * p / log(gamma * p))
}

# Returns the fractional power sparsity index corresponding to size `p' of the
# parameter vector and parameters `gamma' and `delta'.
fractional_power_sparsity <- function(p, gamma, delta) {
  ceiling(gamma * p^delta)
}

# Returns a random `n' x `p' design matrix whose columns are generated from a
# zero-mean multivariate normal distribution with covariance matrix
# `cov_matrix'.
random_design <- function(n, p, cov_matrix) {
  X <- mvrnorm(n = n, mu = rep(0, p), Sigma = cov_matrix)
  X <- rbind(X)  # Force X to be a row vector when n = 1
  X
}

# Returns the sample size corresponding to `p', sparsity index `k', and control
# parameter `theta'.
sample_size <- function(theta, p, k) {
  ceiling(2 * theta * k * log(p - k))
}

# Returns the regularization parameter corresponding to `n', `p', and `k', with
# noise level `sigma'.
lambda_n <- function(n, p, k, sigma) {
  sqrt((2 * sigma^2 * log(k) * log(p - k)) / n)
}

# Returns a plot of probability of successful signed support recovery vs.
# control parameter for the values of p given in `p_values'. `title' is the plot
# title. `uniform' is true if the random designs are to be generated using a
# uniform Gaussian ensemble, and false for non-uniform. `alpha' is the elastic
# net mixing parameter. `sparsity_fn' is the function defining the sparsity
# index, and any extra parameters `...' are passed to this function.
success_prob_plot <- function(p_values, title, uniform, alpha,
                              sparsity_fn, ...) {
  theta_seq <- seq(0.01, 2.4, length.out = 15)  # Control parameters
  sigma <- 0.5  # Noise level
  mu <- 0.1     # Defines the Toeplitz covariance matrix in non-uniform case
  n_trials <- 200
  
  plot_data_list <- list()
  
  # Generate the plot data for each value of p.
  for (p in p_values) {
    k <- sparsity_fn(p, ...)
    S <- random_support(p, k)
    beta_star <- true_beta(p, S)
    S_signed <- signed_support(beta_star)
    
    # Select the appropriate covariance matrix for generating random designs.
    if (uniform)
      cov_matrix <- diag(p)
    else
      cov_matrix <- toeplitz(mu^(0:(p-1)))
    
    # Estimate the probability of successful signed support recovery for each
    # value of the control parameter.
    probs <- numeric(length(theta_seq))
    for (i in seq_along(theta_seq)) {
      theta <- theta_seq[[i]]
      n <- sample_size(theta, p, k)
      
      # TODO: Due to standardization, glmnet doesn't like this case. Should we
      # try to avoid this from ever occurring?
      if (n == 1) {
        warning("Got n = 1 with theta = ", signif(theta), "; skipping.",
                call. = FALSE, immediate. = TRUE)
        next
      }
      
      # Estimate the probability of successful signed support recovery by
      # repeatedly simulating a model with random design.
      count <- 0
      for (j in 1:n_trials) {
        X <- random_design(n, p, cov_matrix)
        y <- X %*% beta_star + rnorm(n, sd = sigma)
        
        lambda <- lambda_n(n, p, k, sigma)
        # TODO: In order to speed up coefficient prediction, may want to
        # precompute the sequence of regularization parameters.
        fit <- glmnet(X, y, intercept = FALSE, alpha = alpha)
        beta <- as.numeric(predict(fit, s = lambda, type = "coefficients",
                                   exact = TRUE, x = X, y = y,
                                   alpha = alpha))[-1]
        
        count <- count + all(signed_support(beta) == S_signed)
      }
      
      probs[[i]] <- count / n_trials
    }
    
    plot_data_list[[p]] <- tibble(theta_seq, probs, p)
  }
  
  # Merge all of the plot data into one big data frame.
  plot_data <- bind_rows(plot_data_list)
  plot_data$p <- factor(plot_data$p, levels = p_values)
  
  # Generate a plot of probability of success vs. the control parameter. Display
  # a line for each value of p.
  plot <- ggplot(plot_data, aes(theta_seq, probs, color = p)) +
    coord_cartesian(xlim = c(0, 2.4), ylim = c(0, 1)) +
    geom_point() +
    geom_line() +
    geom_vline(xintercept = 1, linetype = "dotted") +
    labs(title = title, color = TeX("$p$")) +
    xlab(TeX("Control parameter $\\theta(n, p, k)$")) +
    ylab("Probability of success") +
    theme(aspect.ratio = 1)
  
  invisible(plot)
}

# Writes the plot `plot' to disk with appropriate options. `filename' is the
# filename. Returns `plot' for use in pipe operator chains.
save_plot <- function(plot, filename) {
  ggsave(filename, plot, device = "png", width = 4, height = 4)
  plot
}

# Set global theme for plots.
theme_set(theme_bw())
theme_update(plot.title = element_text(hjust = 0.5))

# Values of p for which to run the simulations.
p_values <- c(128, 256, 512)

# Sparsity parameters.
gamma <- 0.4
delta <- 0.75

# Simulations from the paper ----------------------------------------------

# Uniform Gaussian ensemble.

success_prob_plot(
  p_values,
  title = "Linear sparsity",
  uniform = TRUE,
  alpha = 1,
  sparsity_fn = linear_sparsity,
  gamma
) %>%
  save_plot("uniform_linear_sparsity_alpha_1.png") %>%
  print

success_prob_plot(
  p_values,
  title = "Sublinear sparsity",
  uniform = TRUE,
  alpha = 1,
  sparsity_fn = sublinear_sparsity,
  gamma
) %>%
  save_plot("uniform_sublinear_sparsity_alpha_1.png") %>%
  print

success_prob_plot(
  p_values,
  title = "Fractional power sparsity",
  uniform = TRUE,
  alpha = 1,
  sparsity_fn = fractional_power_sparsity,
  gamma,
  delta
) %>%
  save_plot("uniform_fractional_power_sparsity_alpha_1.png") %>%
  print

# Non-uniform Gaussian ensemble.

success_prob_plot(
  p_values,
  title = "Linear sparsity",
  uniform = FALSE,
  alpha = 1,
  sparsity_fn = linear_sparsity,
  gamma
) %>%
  save_plot("nonuniform_linear_sparsity_alpha_1.png") %>%
  print

success_prob_plot(
  p_values,
  title = "Sublinear sparsity",
  uniform = FALSE,
  alpha = 1,
  sparsity_fn = sublinear_sparsity,
  gamma
) %>%
  save_plot("nonuniform_sublinear_sparsity_alpha_1.png") %>%
  print

success_prob_plot(
  p_values,
  title = "Fractional power sparsity",
  uniform = FALSE,
  alpha = 1,
  sparsity_fn = fractional_power_sparsity,
  gamma,
  delta
) %>%
  save_plot("nonuniform_fractional_power_sparsity_alpha_1.png") %>%
  print

# Custom simulations ------------------------------------------------------

# alpha = 0.75

success_prob_plot(
  p_values,
  title = "Linear sparsity",
  uniform = TRUE,
  alpha = 0.75,
  sparsity_fn = linear_sparsity,
  gamma
) %>%
  save_plot("uniform_linear_sparsity_alpha_075.png") %>%
  print

success_prob_plot(
  p_values,
  title = "Sublinear sparsity",
  uniform = TRUE,
  alpha = 0.75,
  sparsity_fn = sublinear_sparsity,
  gamma
) %>%
  save_plot("uniform_sublinear_sparsity_alpha_075.png") %>%
  print

success_prob_plot(
  p_values,
  title = "Fractional power sparsity",
  uniform = TRUE,
  alpha = 0.75,
  sparsity_fn = fractional_power_sparsity,
  gamma,
  delta
) %>%
  save_plot("uniform_fractional_power_sparsity_alpha_075.png") %>%
  print

# alpha = 0.25

success_prob_plot(
  p_values,
  title = "Linear sparsity",
  uniform = TRUE,
  alpha = 0.25,
  sparsity_fn = linear_sparsity,
  gamma
) %>%
  save_plot("uniform_linear_sparsity_alpha_025.png") %>%
  print

success_prob_plot(
  p_values,
  title = "Sublinear sparsity",
  uniform = TRUE,
  alpha = 0.25,
  sparsity_fn = sublinear_sparsity,
  gamma
) %>%
  save_plot("uniform_sublinear_sparsity_alpha_025.png") %>%
  print

success_prob_plot(
  p_values,
  title = "Fractional power sparsity",
  uniform = TRUE,
  alpha = 0.25,
  sparsity_fn = fractional_power_sparsity,
  gamma,
  delta
) %>%
  save_plot("uniform_fractional_power_sparsity_alpha_025.png") %>%
  print
