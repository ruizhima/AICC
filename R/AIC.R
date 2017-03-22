# This function calculate AIC
#' AIC
#'
#' @param covmat  covariance matrix of innovations
#' @param p  order of autoregression
#' @param m  dimension of the variable
#' @param n  number of observations
#'
#' @return
#' @export
#'
#' @examples
AIC <- function(covmat, p, m, n) {
  out <-n*(log(det(covmat))+m)+2*(p*m^2+m*(m+1)/2)
  return(out)
}
