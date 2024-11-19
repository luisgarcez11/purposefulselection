#' Accuracy ROC
#'
#' @param fit estimate accuracy of a model
#'
#' @return
#' @export
#'
#' @examples
accuracy <- function(fit){
  
  roc_curve <- suppressMessages({pROC::roc(fit$y , fit$fitted.values ) })
  
  best_cut_point <- suppressMessages(roc_curve$thresholds[which.min(sqrt( (1-roc_curve$sensitivities)^2+(1-roc_curve$specificities)^2 ))])
  
  y_pred <- case_when(fit$fitted.values <= best_cut_point ~ 0,
                      fit$fitted.values > best_cut_point ~ 1)
  
  accuracy <- sum(y_pred == fit$y)/ length(fit$y)
  
  return(accuracy)
  
  
}
