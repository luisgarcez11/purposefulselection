#' Check linearity assumptions
#'
#' @param fit
#' @param vars
#'
#' @return
#' @export
#'
#' @examples
check_linearity <- function(fit, vars = c("age", "height")) {

  #linear predictions
  linear_pred <- fit$model_fit$linear.predictors

  #iterate variables
  for (var in vars) {
    p <- ggplot2::ggplot(
      data = fit$data,
      mapping = aes(x = fit$data[[var]], y = linear_pred)
    ) +
      ggplot2::geom_smooth() +
      ggplot2::geom_jitter() +
      ggplot2::labs(x = var, y = "Linear Predictor")


    print(p)
    invisible(readline(prompt="Press [enter] to continue:"))


  }

  message(paste("STEP 5: check linearity for",  paste(vars, collapse = ", "), "\n"))

  return(fit)


}
