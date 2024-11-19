


#' Check significance for variables excluded in step1
#'
#' @param nonconfounded_fit object returned from check_confounding_for_excluded_in_step2
#'
#' @return
#' @export
#'
#' @examples
check_significance_for_excluded_in_step1 <- function(nonconfounded_fit, sig_level = 0.05){

  #model
  model_fit <- nonconfounded_fit$model_fit
  model_fit_form <- stringr::str_squish(paste(deparse(model_fit$formula), collapse = "")) %>%
    stringr::str_remove_all(pattern = "\"")

  #variable lookup
  lookup <- nonconfounded_fit$lookup

  #store all fits and combinations from adding excluded variables in step1
  check_vars_excluded_step1 <- tibble::tibble(
    added_variable = nonconfounded_fit$excluded_vars_step1) %>%
    left_join(nonconfounded_fit$lookup_vars_data, by = c("added_variable" = "term"))

  if(nonconfounded_fit$term_based){

    check_vars_excluded_step1 <- check_vars_excluded_step1 %>%
      mutate(formula = purrr::map(.x = added_variable,
                                  ~ sprintf("%s + %s",
                                            stringr::str_squish(paste(deparse(as.formula(nonconfounded_fit$model_fit$formula )),
                                                                      collapse = "")),
                                            .x ))) %>%
      mutate(fit = purrr::map(.x = formula, ~glm(formula = .x, family = "binomial", data = nonconfounded_fit$data ) )) %>%
      mutate(results = purrr::map(.x = fit, ~broom::tidy(.x) %>%
                                    filter(term != "(Intercept)") %>%
                                    left_join(nonconfounded_fit$lookup_vars_data, by = "term") ))%>%
      mutate(added_p_value_wald = purrr::map2_dbl(.x = results,
                                                  .y = added_variable,
                                                  ~ .x$p.value[which(.x$term == .y)] )) %>%
      mutate(added_p_value = added_p_value_wald)

  }else{

    check_vars_excluded_step1 <- check_vars_excluded_step1 %>%
      mutate(formula = purrr::map(.x = twin_terms,
                                  ~ sprintf("%s + %s",
                                            stringr::str_squish(paste(deparse(as.formula(nonconfounded_fit$model_fit$formula )),
                                                                      collapse = "")),
                                            paste0(.x, collapse = " + ") ))) %>%
      mutate(fit = purrr::map(.x = formula, ~glm(formula = .x, family = "binomial", data = nonconfounded_fit$data ) )) %>%
      mutate(results = purrr::map(.x = fit, ~broom::tidy(.x) %>%
                                    filter(term != "(Intercept)") %>%
                                    left_join(nonconfounded_fit$lookup_vars_data, by = "term") )) %>%
      mutate(added_p_value_lr = purrr::map_dbl(.x = fit,
                                               ~ likelihood_ratio(model_fit, .x )$p_value )) %>%
      mutate(added_p_value = added_p_value_lr)

  }

  added_variables <- c()

  #adding if any is significant from
  if(any(check_vars_excluded_step1$added_p_value < sig_level)){

    #formula with all significant added backs
    i = which(check_vars_excluded_step1$added_p_value < sig_level)
    formula1 = paste(model_fit_form,
                     paste(check_vars_excluded_step1$added_variable[i],
                           collapse = " + "),
                     sep = " + ")


    #fit multivariable
    model_fit <- glm(formula = formula1,
                     family = "binomial",
                     data = nonconfounded_fit$data )

    model_results <- broom::tidy(model_fit)

    #add lr p_values
    if(!nonconfounded_fit$term_based){
      model_results <- add_lr_p_values(m_results = model_results,
                                       ref_model = model_fit,
                                       data1 = nonconfounded_fit$data,
                                       lookup_table = nonconfounded_fit$lookup_vars_data) %>%
        mutate(p.value = p_value_lr)
    }

    #if any included model terms not significant they will be removed
    while(any(model_results$p.value[which(model_results$term %in% nonconfounded_fit$excluded_vars_step1)] > sig_level) ){

      #set terms to remove
      term_to_remove <- model_results %>%
        arrange(desc(p.value)) %>%
        left_join(nonconfounded_fit$lookup_vars_data) %>%
        filter(term %in% nonconfounded_fit$excluded_vars_step1 ) %>%
        filter(p.value > sig_level)

      if(nrow(term_to_remove) == 0 ){break}

      if(nonconfounded_fit$term_based){
        variable_to_remove <- term_to_remove$term
      }else{
        variable_to_remove <- term_to_remove$twin_terms[[1]]
      }


      #store model
      model_fit <- glm(formula = update_formula(formula1,  paste("-", paste0(variable_to_remove, collapse = " - ") )),
                       family = "binomial",
                       data = nonconfounded_fit$data )

      #add LR p-values
      if(nonconfounded_fit$sig_lr){
        model_results <- broom::tidy(model_fit) %>%
          add_lr_p_values(ref_model = model_fit,
                          data1 = nonconfounded_fit$data,
                          lookup_table = nonconfounded_fit$lookup_vars_data) %>%
          mutate(p.value = p_value_lr)
      }else{
        #add wald p-values
        model_results <- broom::tidy(model_fit) %>%
          filter(term != "(Intercept)" ) %>%
          left_join(nonconfounded_fit$lookup_vars_data, by = "term") %>%
          relocate("variable", .after = "term")
        }

    }

    model_results = model_fit %>%
      broom::tidy() %>%
      filter(term != "(Intercept)") %>%
      left_join(lookup_vars(model_fit), by = "term")%>%
      relocate("variable", .after = "term")

    added_variables <- setdiff(model_results$term, nonconfounded_fit$model_results$term)

    message(paste("added",
                  paste(added_variables, collapse = ", "),
                  "to the model, that was excluded in step1 (univariate analysis)"))

  }else{

    message("STEP4: PRELIMINARY MAIN EFFECTS MODEL: none of the variables excluded in step1 (univariate analysis) were added to the model \n")

    model_fit <- nonconfounded_fit$model_fit
    model_results <- nonconfounded_fit$model_results


  }

  #store info on this step
  check_vars_excluded_step1 <- check_vars_excluded_step1 %>%
    mutate(added = case_when(added_variable %in% added_variables ~ "Yes",
                             TRUE ~ "No"))

  #final model fit
  model_fit_form = stringr::str_remove_all(paste(deparse(model_fit$formula), collapse = ""), "\"" )

  #calculate lr and wald p-values
  model_results = model_fit %>%
    broom::tidy()

  #define wald and lr p-values
  if(nonconfounded_fit$sig_lr){
    model_results = model_results %>%
      add_lr_p_values(data1 = nonconfounded_fit$data,
                      ref_model = model_fit,
                      lookup_table = nonconfounded_fit$lookup_vars_data) %>%
      dplyr::rename("p.value_wald" = "p.value") %>%
      mutate(p.value = p_value_lr) }else{

      model_results = model_results %>% mutate(p.value_wald = p.value)
    }

  #model performance
  model_fit_auc <- suppressMessages(pROC::auc(model_fit$y, model_fit$fitted.values))
  model_fit_acc <- suppressMessages(accuracy(model_fit))


  #numeric_excluded
  num_terms_excluded = setdiff(nonconfounded_fit$num_terms, model_results$term)
  num_terms_included = setdiff(nonconfounded_fit$num_terms, num_terms_excluded)

  #categorical_excluded
  cat_terms_excluded = setdiff(nonconfounded_fit$cat_terms, model_results$term)
  cat_terms_included = setdiff(nonconfounded_fit$cat_terms, cat_terms_excluded)

  return(list("data" = nonconfounded_fit$data,
              "outcome_var" = nonconfounded_fit$outcome_var,
              "lookup" = nonconfounded_fit$lookup,
              "lookup_vars_data" = nonconfounded_fit$lookup_vars_data,
              "univariate_results" = nonconfounded_fit$univariate_results,
              "term_based" = nonconfounded_fit$term_based,
              "force_entry" = nonconfounded_fit$force_entry,
              "num_terms" = nonconfounded_fit$num_terms,
              "cat_terms" = nonconfounded_fit$cat_terms,
              "num_terms_excluded" = num_terms_excluded,
              "num_terms_included" = num_terms_included,
              "cat_terms_excluded" = cat_terms_excluded,
              "cat_terms_included" = cat_terms_included,


              "reduced_model_results" = nonconfounded_fit$reduced_model_results,
              "reduced_model_auc" = nonconfounded_fit$reduced_model_auc,
              "reduced_model_acc" = nonconfounded_fit$reduced_model_acc,

              "larger_model_results" = nonconfounded_fit$larger_model_results,
              "larger_model_auc" = nonconfounded_fit$larger_model_auc,
              "larger_model_acc" = nonconfounded_fit$larger_model_acc,

              "check_coef_variation" = nonconfounded_fit$check_coef_variation,
              "check_coef_variation_results" = nonconfounded_fit$model_results,
              "check_coef_variation_auc" = nonconfounded_fit$model_fit_auc,
              "check_coef_variation_acc" = nonconfounded_fit$model_fit_acc,
              "n_coef_added_back" = nonconfounded_fit$n_coef_added_back,

              "check_vars_excluded_step1" = check_vars_excluded_step1,
              "excluded_vars_step1" = nonconfounded_fit$excluded_vars_step1,
              "excluded_vars_step2" = nonconfounded_fit$excluded_vars_step2,
              "model_fit" = model_fit,
              "model_results" = model_results,
              "model_fit_auc" = model_fit_auc,
              "model_fit_acc" = model_fit_acc,

              #multivariate tests
              "complete_larger_lr_p_value" = nonconfounded_fit$complete_larger_lr_p_value,
              "complete_larger_wald_p_value" = nonconfounded_fit$complete_larger_wald_p_value,
              "larger_reduced_lr_p_value" = nonconfounded_fit$larger_reduced_lr_p_value,
              "larger_reduced_wald_p_value" = nonconfounded_fit$larger_reduced_wald_p_value,
              "reduced_conf_wald_p_value" = nonconfounded_fit$reduced_conf_wald_p_value,
              "reduced_conf_lr_p_value" = nonconfounded_fit$reduced_conf_lr_p_value
  ))

}
