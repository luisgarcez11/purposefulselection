
#' History fit
#'
#' @param fit returned object from check_significance_for_excluded_in_step1
#'
#' @return
#' @export
#'
#' @examples
history <- function(fit,
                    return_model = FALSE){

  message("compiling history...")

  variables_included <- unique(fit$univariate_results$term)
  if(any(is.na(variables_included))){variables_included <- variables_included[-which(is.na(variables_included))]}



  #set variable type
  num_terms = variable_types(fit$data)$num_vars
  cat_terms = variable_types(fit$data)$cat_vars

  #description table for numeric
  if(length(variables_included)>= 50){

    describe_numeric_data = tibble()
    limit = ceiling(length(variables_included)/50 )-1
    for ( i in 1:limit){
      m = 49
      lag = 50*i
      indexes = i:m + lag
      indexes = indexes[which(indexes< length(variables_included))]

      data <- as.data.frame(psych::describe(fit$data[,variables_included[indexes]], check = TRUE, omit = TRUE) ) %>%
        tibble::rownames_to_column("variable")
      describe_numeric_data = bind_rows(describe_numeric_data,data)
      }
  }else{
    describe_numeric_data <- as.data.frame(psych::describe(fit$data[,variables_included], check = TRUE, omit = TRUE) ) %>%
      tibble::rownames_to_column("variable")
  }


  #step0 - data description
  history_fit <- fit$univariate_results[, c("term", "variable", "type")] %>%
    mutate(force_entry = case_when(term %in% fit$force_entry ~"yes",
                                   variable %in% fit$force_entry ~"yes",
                                   TRUE ~ "no")) %>%
    left_join(describe_numeric_data %>%
                mutate_at("variable", ~stringr::str_remove_all(., "\\*")) %>%
                select(variable, skew, kurtosis, se, sd),
              by = c("term" = "variable"))  %>%
    mutate(n_occurr = purrr::map2_int(.x = term, .y = type, .f = ~case_when(.y == "categorical" ~ sum(fit$data[[.x]])))) %>%
    mutate(share = purrr::map2_dbl(.x = term,
                                   .y = type,
                                   .f = ~case_when(.y == "categorical" ~ sum(fit$data[[.x]])/nrow(fit$data)  ) ))
  # %>%
  #   mutate(largest_share = purrr::map_dbl(.x = term,
  #                                         .f = ~ round(max(prop.table(table(unlist(fit$data[,.x])))),2) )) %>%
  #   mutate(largest_share = case_when(type == "categorical" ~ largest_share,
  #                                    TRUE ~ NA_integer_))


  #step1 - univariate
  history_fit <- history_fit %>%
    left_join(fit$univariate_results[,c("term", "p.value")] %>%
                distinct(), by = "term") %>%
    dplyr::rename("univariate" = "p.value")

  if(fit$term_based == TRUE){
    history_fit <- history_fit %>%
      mutate(univariate_desc = case_when(term %in% fit$excluded_vars_step1 ~ "out",
                                         TRUE ~ "in"))
  }else{
    history_fit <- history_fit %>%
      mutate(univariate_desc = case_when(term %in% fit$excluded_vars_step1 ~ "out",
                                         TRUE ~ "in"))
    }

  #step2 - larger model
  history_fit <- history_fit %>%
    left_join(fit$larger_model_results[,c("term", "p.value")] %>%
                distinct(), by = "term") %>%
    dplyr::rename("larger_model" = "p.value") %>%
    mutate(larger_model_desc = univariate_desc)

  #step2 - reduced model
  history_fit <- history_fit %>%
    left_join(fit$reduced_model_results[,c("term", "p.value")] %>%
                distinct(), by = "term") %>%
    dplyr::rename("reduced_model" = "p.value")


    history_fit <- history_fit %>%
      mutate(reduced_model_desc = case_when(term %in% fit$excluded_vars_step2 ~ "out",
                                            univariate_desc == "out" ~ "out",
                                            TRUE ~ "in"))

  #step3 - check confouding
  if(any(!is.na(fit$check_coef_variation))){
  added_terms <- fit$check_coef_variation$added_term[which(fit$check_coef_variation$selected == "Yes")] %>%
    stringr::str_split(pattern = "\\+") %>%
    unlist() %>%
    stringr::str_squish()

  history_fit <- history_fit %>%
    left_join(fit$check_coef_variation_results[,c("term", "p.value")] %>%
                distinct(), by = "term") %>%
    dplyr::rename("check_conf_model" = "p.value")

  if(fit$term_based){
    history_fit <- history_fit %>% mutate(check_confounding_desc = case_when(term %in% added_terms ~ "in",
                                              reduced_model_desc == "out" ~ "out",
                                              TRUE ~ "in"))
  }else{
    history_fit <- history_fit %>% mutate(check_confounding_desc = case_when(variable %in% added_terms ~ "in",
                                                             reduced_model_desc == "out" ~ "out",
                                                             TRUE ~ "in"))
    }

  }else{
    history_fit <- history_fit %>%
      left_join(fit$check_coef_variation_results[,c("term", "p.value")] %>%
                  distinct(), by = "term") %>%
      dplyr::rename("check_conf_model" = "p.value") %>%
      mutate(check_confounding_desc = reduced_model_desc)
  }

  #step4 - check variables excluded on step1
  added_terms <- fit$check_vars_excluded_step1

  if(nrow(added_terms)> 0){


      history_fit <- history_fit %>%
        left_join(fit$check_vars_excluded_step1[,c("added_variable", "added_p_value", "added")],
                  by = c("term" = "added_variable")) %>%
        dplyr::rename("check_uni_exclusions" = "added_p_value" ) %>%
        mutate(check_uni_exclusions_desc = case_when( check_confounding_desc == "in" ~ "in",
                                                      check_confounding_desc == "out" & added == "Yes" ~ "in",
                                                      TRUE ~ "out")) %>%
        select(-added)


  }else{
    history_fit <- history_fit %>%
      mutate(check_uni_exclusions_desc = check_confounding_desc)
  }

  #step6- interations
  history_fit <- history_fit %>%
    full_join(fit$model_results[,c("term", "p.value")] %>%
                filter(term != "(Intercept)") %>%
                distinct(), by = "term") %>%
    dplyr::rename("final" = "p.value" ) %>%
    mutate(final_desc = case_when( check_uni_exclusions_desc == "in" ~ "in",
                                   check_uni_exclusions_desc == "out" ~ "out",
                                          TRUE ~ "in")) %>%
    mutate_if(is.double, ~round(., 4))


  #add AUC
  history_fit <- history_fit %>%
    add_row(variable = "AUC",
            larger_model_desc = as.character(round(fit$larger_model_auc,4)),
            reduced_model_desc = as.character(round(fit$reduced_model_auc,4)),
            check_confounding_desc = as.character(round(fit$check_coef_variation_auc , 4)),
            check_uni_exclusions_desc = as.character(round(fit$model_fit_auc , 4)),
            final_desc = as.character(round(fit$model_fit_auc , 4))
            ) %>%
    add_row(variable = "Accuracy",
            larger_model_desc = as.character(round(fit$larger_model_acc,4)),
            reduced_model_desc = as.character(round(fit$reduced_model_acc,4)),
            check_confounding_desc = as.character(round(fit$check_coef_variation_acc , 4)),
            check_uni_exclusions_desc = as.character(round(fit$model_fit_acc , 4)),
            final_desc = as.character(round(fit$model_fit_acc , 4)))

  #add multivariate testss
  history_fit <- history_fit %>%
    add_row(variable = "multivariate_wald",
            larger_model_desc = tryCatch(as.character(round(fit$complete_larger_wald_p_value,3)),
                                        error = function(e) {return("NA")}),
            reduced_model_desc = as.character(round(as.numeric(fit$larger_reduced_wald_p_value) , 3)),
            check_confounding_desc = as.character(round(as.numeric(fit$reduced_conf_wald_p_value) , 3))
    ) %>%
    add_row(variable = "multivariate_lr",
            larger_model_desc = as.character(round(fit$complete_larger_lr_p_value , 3)),
            reduced_model_desc = as.character(round(fit$larger_reduced_lr_p_value,3)),
            check_confounding_desc = as.character(round(fit$reduced_conf_wald_p_value , 3))
    )

  #add outcome proportion
  history_fit <- history_fit %>%
    add_row(variable = "outcome_proportion",
            final_desc = as.character(round(max(prop.table(table(fit$data[[fit$outcome_var]]))), 3))
    ) %>%
    add_row(variable = "n_added_back",
            check_confounding_desc = as.character(fit$n_coef_added_back)
    )


  return(list("history" = history_fit,
              "model" = ifelse(return_model, fit$model_fit, "no_model"),
              "predictor" =  fit$model_fit$fitted.values,
              "y" = fit$model_fit$y))


}
