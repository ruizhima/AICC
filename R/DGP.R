#' Data generating process described in Hurvich and Tsai (1993)
#'
#' @param sample_size Desired sample size
#' @param variance_covariance_matrix covariance matrix for the error term
#' @param ar_coefficients AR coefficients, no constant term. If AR(p) with m variables, it should be a vertical stack of all coefficient matrices (thus the dimension is mp*m)
#'
#' @return Returns a single series with sample size = sample_size
#' @export
#'
#' @examples
dgp <- function(sample_size, variance_covariance_matrix, ar_coefficients) {
  epsilon <- mvrnorm(sample_size,matrix(0,2,1),variance_covariance_matrix)
  y <- matrix(NaN,sample_size,2)
  lags <- nrow(ar_coefficients)/ncol(ar_coefficients)
  lag_1 <- ar_coefficients[1:2,1:2]
  if (lags == 1) {
    for (i in 1:sample_size) {
      if (i == 1) {
        y[i,1:2] = epsilon[i,1:2]
      }
      else {
        y[i,1:2] = lag_1 %*% y[i-1,1:2] + epsilon[i,1:2]
      }
    }
  }
  if (lags == 2) {
    lag_2 <- ar_coefficients[3:4,1:2]
    for (i in 1:sample_size) {
      if (i == 1) {
        y[i,1:2] = epsilon[i,1:2]
      }
      else if (i == 2) {
        y[i,1:2] = lag_1 %*% y[i-1,1:2] + epsilon[i,1:2]
      }
      else {
        y[i,1:2] = lag_1 %*% y[i-1,1:2] + lag_2 %*% y[i-2,1:2] + epsilon[i,1:2]
      }
    }
  }
  return(y)
}
