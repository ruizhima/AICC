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
SIC <- function(covmat, p, m, n) {
  out <-n*(log(det(covmat))+m)+p*m^2*log(n)
  return(out)
}
