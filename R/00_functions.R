
#' Likelihood ratio
#'
#' @param null_model
#' @param model
#'
#' @return
#' @export
#'
#' @examples
likelihood_ratio <- function(null_model, model){

  test <- lmtest::lrtest(null_model, model)

  return(list("lr_stat" = test$Chisq[2],
              "df" = test$Df[2],
              "p_value" = test$`Pr(>Chisq)`[2]))
}

#' Formula Combination
#'
#' @param vars
#'
#' @return
#' @export
#'
#' @examples
formula_combination <- function(vars, term_min = 1,
                                term_limit = length(vars),
                                operator = " + ",
                                term_based,
                                lookup_vars_data){

  #group terms into variables
  if(!term_based){
    vars <- lookup_vars_data %>%
      mutate(whole = purrr::map(.x = twin_terms,
                                .f = ~vars[which(vars %in% .x )])) %>%
      pull(whole) %>%
      unique() %>%
      lapply(., FUN = function(x) paste0(x, collapse = operator)) %>%
      unlist() %>%
      setdiff("")

    if(term_limit > length(vars) ){term_limit = length(vars) }
  }

  all_combinations_list <- lapply(term_min:term_limit,
                                  function(k) combn(vars, k, simplify = FALSE))
  all_combinations_formulas <- c()
  for(i in all_combinations_list){
    for( j in i){
      all_combinations_formulas <- c(all_combinations_formulas,
                                     paste(j, collapse = operator))
    }
  }
  return(all_combinations_formulas)
}



recateogrize <- function(fit, var, before, after, ref_value = NULL){

  testthat::expect_equal(length(before), length(after), info = "length should be the same before and after")

  rename_vector <- after
  names(rename_vector) <- before

  data <- fit$data %>%
    mutate_at(var, ~rename_vector[.])

  if(!is.null(ref_value) ){
    if(!ref_value %in% after){stop("ref_value must be within after variables")}
    data[[var]] <- relevel(as.factor(data[[var]]), ref = ref_value)}

  fit$data <- data
  fit$model_fit <- glm(formula = fit$model_fit$formula, family = "binomial", data = fit$data )
  fit$model_results <- broom::tidy(fit$model_fit)

  message(paste("RECAT: For variable", var, "replaced these categories", paste(before, collapse = ", ") ,  "with these ones:", paste(after, collapse = ", "), "\n" ))

  return(fit)


}


#' Lookup
#'
#' @param model
#' @param interactions render interactions in lookuptable
#'
#' @return
#' @export
#'
#' @examples
lookup_vars <- function(model, interactions = FALSE){

  lookup_vars_data <- tibble(term = character(), variable = character())

  for ( var1 in  names(model$model)){

    options_var1 <- names(table(model$model[[var1]]))

    if(length(options_var1) < 20){

    lookup_vars_data <- lookup_vars_data %>%
      bind_rows( tibble( term = paste(var1, options_var1, sep = ""), variable =  var1) )
    }

    lookup_vars_data <- lookup_vars_data %>%
      bind_rows( tibble( term = var1, variable =  var1) )

    if(interactions == TRUE){

    for ( var2 in  names(model$model)){

      options_var2 <- names(table(model$model[[var2]]))


      interaction_variable  = paste(var1, var2, sep = ":")

      if(length(options_var2) < 20){

      interaction_options <- apply(expand.grid(paste0(var1, options_var1, sep = ""),
                                               paste0(var2, options_var2, sep = ""), KEEP.OUT.ATTRS = FALSE),
                                   1, paste, collapse=":")

      lookup_vars_data <- lookup_vars_data %>%
        bind_rows(tibble(term = interaction_options,
                         variable =  interaction_variable))

      }


      lookup_vars_data <- lookup_vars_data %>%
        bind_rows(tibble(term = interaction_variable,
                         variable =  interaction_variable))




    }}
  }


  return(lookup_vars_data)

}


#' Update formula
#'
#' @param formula
#' @param term
#'
#' @return
#' @export
#'
#' @examples
update_formula <- function(formula, ...){

  form = stringr::str_remove_all(paste(deparse(formula), collapse = ""), "\"" ) %>%
    stringr::str_squish()

  form = paste(form, ...)

  return(formula(form))


}


#' Title
#'
#' @param m_results
#' @param data
#' @param lookup_table
#'
#' @return
#' @export
#'
#' @examples
add_lr_p_values = function(m_results, ref_model, data1, lookup_table){

  model_results1 <- m_results %>%
    filter(term != "(Intercept)" ) %>%
    left_join(lookup_table, by = "term") %>%
    mutate(formula = purrr::map(.x = twin_terms,
                                ~ sprintf("%s - %s",
                                          stringr::str_squish(paste(deparse(as.formula(ref_model$formula )),
                                                                    collapse = "")),
                                          paste0(.x, collapse = " - ") ))) %>%
    mutate(fit = purrr::map(.x = formula, ~glm(formula = .x, family = "binomial", data = data1 ) )) %>%
    mutate(results = purrr::map(.x = fit, ~broom::tidy(.x) %>%
                                  filter(term != "(Intercept)") %>%
                                  left_join(lookup_table, by = "term") )) %>%
    mutate(p_value_lr = purrr::map_dbl(.x = fit,
                                             ~ likelihood_ratio(ref_model, .x )$p_value ))

  return(model_results1)
}

