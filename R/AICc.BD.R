# This function calculates AICc.BD based on Brockwell and Davis (1991)
#' AICc (Brockwell and Davis (1991))
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
AICc.BD <-function(covmat, p, m, n){
  out <- n*log(det(covmat)+m)+2*n*(1+p*m^2)/(n-p*m^2-2)
  return(out)
}
