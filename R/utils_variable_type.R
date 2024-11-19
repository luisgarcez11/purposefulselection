#' Get variable type
#'
#' @param x vector variable
#'
#' @return
#' @export
#'
#' @examples
variable_type = function(x){
  
  len_var = length(unique(x))
  categorical = suppressWarnings(all(is.na(as.numeric(x))))
  
  type = case_when(len_var <= 2 ~ "categorical",
                   categorical ~ "categorical",
                   len_var > 2 ~  "numerical")
  
  return(type)
  
}


#' Get numeric variables 
#'
#' @param .data vector variable
#'
#' @return
#' @export
#'
#' @examples
variable_types = function(.data){
  
  num_vars = names(.data)[sapply(.data, FUN = function(x) variable_type(x)) == "numerical"]
  cat_vars = names(.data)[sapply(.data, FUN = function(x) variable_type(x)) == "categorical"]
  
  return(list("num_vars" = num_vars, "cat_vars" = cat_vars))
  
}
