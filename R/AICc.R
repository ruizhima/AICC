# This function calculates AICc based on Hurvich and Tsai (1991)
#' AICc (Hurvich and Tsai (1991))
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
AICc <-function(covmat, p, m, n){
  out <- n*log(det(covmat))+(n*(n*m+p*m^2))/(n-p*m-m-1)
  return(out)
}
