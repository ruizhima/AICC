dgp <- function(sample_size, variance_covariance_matrix, ar_coefficients){
  x <- runif(sample_size)
  epsilon <- rnorm(sample_size, sd = error_sd)
  y <- intercept + slope * x + epsilon
  return(y)
}
