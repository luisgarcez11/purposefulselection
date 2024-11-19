
#' Add categorical share
#'
#' @param univariate_results univariate results
#' @param .data data used
#'
#' @return
#' @export
#'
#' @examples
add_categorical_share = function(univariate_results, .data){

  univariate_results <- univariate_results %>%
    mutate(share = purrr::map2_dbl(.x = term,
                                   .y = type,
                                   .f = ~case_when(.y == "numeric" ~ NA,
                                                   .y == "categorical" ~ prop.table(table(.data[[.x]])))))
    mutate(share = case_when(type == "numeric" ~ NA_integer_,
                             type == "categorical" ~ prop.table(table(.data))))

}
