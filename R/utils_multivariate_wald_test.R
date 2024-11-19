#' Multivariate Wald
#'
#' @param fit fit 
#' @param vars_to_test variables to test if they are zero
#'
#' @return p_value regarding chisquare test with df = difference in coeffficients
#' @export
#'
#' @examples
multivariate_wald <- function(fitted_model, coefs_to_test){
  
  #check if coefficients found
  if(!all( coefs_to_test %in% names(coef(fitted_model)))){stop("coefficients not found")}
  
  #design matrix
  X = model.matrix(fitted_model)
  
  #inverse variance
  V = diag(nrow = length(fitted_model$fitted.values),
           ncol = length(fitted_model$fitted.values),
           x = (fitted_model$fitted.values*(1- fitted_model$fitted.values)) )
  
  #betas
  betas = matrix(coefficients(fitted_model)[coefs_to_test])
  
  #inverse variance
  inverse_variance <- (t(X) %*% V %*% X)[coefs_to_test, coefs_to_test]
  
  variance <- stats::vcov(fitted_model)[coefs_to_test, coefs_to_test]
  inverse_variance <-  solve(variance)
  
  #check if variance well calculated
  # if(!all(solve(inverse_variance) == stats::vcov(fitted_model)[coefs_to_test, coefs_to_test])){"cov-var matrix not estimated properly"}
 
  #wald statistic
  W = t(betas) %*%  inverse_variance %*% betas
  
  p_value = pchisq(q = W,  df = length(coefs_to_test), lower.tail = FALSE )
  
  return(p_value)
  
}
