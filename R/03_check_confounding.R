#' Fit all combinations with LR p-values
#'
#' @param multivariable_fit object returned by multivariable_log_fit
#' @param max_term_combinations possible maximum term combinations to be tested
#' @param sig_lr LR significance
#'
#' @return
#' @export
#'
#' @examples
combinations_lr_p_values <- function(multivariable_fit, max_term_combinations = NULL, sig_lr = FALSE){

  if(is.null(max_term_combinations)){
  all_combinations_formulas <- formula_combination(vars = multivariable_fit$excluded_vars_step2,
                                                   term_min = 1,
                                                   term_limit = length(multivariable_fit$excluded_vars_step2),
                                                   term_based = multivariable_fit$term_based,
                                                   lookup_vars_data = multivariable_fit$lookup_vars_data)}else{
    all_combinations_formulas <- formula_combination(multivariable_fit$excluded_vars_step2,
                                                     term_min = 1,
                                                     term_limit = min(c(max_term_combinations,
                                                                      length(multivariable_fit$excluded_vars_step2)), na.rm = TRUE),
                                                     term_based = multivariable_fit$term_based,
                                                     lookup_vars_data = multivariable_fit$lookup_vars_data)
  }

  #store all fits and combinations
  combination_fits <- tibble::tibble( formula = c(sprintf("%s + %s",
                                                      stringr::str_squish(paste(deparse(as.formula(multivariable_fit$reduced_model_fit$formula)), collapse = "")),
                                                      all_combinations_formulas),
                                              stringr::str_squish(paste(deparse(as.formula(multivariable_fit$reduced_model_fit$formula)), collapse = "")))) %>%
    mutate(added_term = c( all_combinations_formulas, NA)) %>%
    mutate(n_added_term = purrr::map_int(.x = added_term, .f = function(x) 1 + sum(stringr::str_count(x, "\\+")))) %>%
    mutate(fit = purrr::map(.x = formula, ~speedglm::speedglm(formula = .x, family = binomial(), data = multivariable_fit$data ) ,
                            .progress = "combinations_fit"))

  if(sig_lr){
    combination_fits = combination_fits%>%
      mutate(results = purrr::map(.x = fit,
                                ~ broom::tidy(.x) %>%
                                  mutate(original_fit = list(.x)) %>%
                                  filter(term != "(Intercept)" ) %>%
                                  left_join(multivariable_fit$lookup_vars_data, by = "term")%>%
                                  mutate(form_without_var = purrr::map2(.x = variable, .y = original_fit, ~update_formula(formula(.y), "-", .x))) %>%
                                  mutate(model_without_var = purrr::map(.x = unlist(form_without_var), ~speedglm::speedglm(formula = .x, family = binomial(), data = multivariable_fit$data ))) %>%
                                  mutate(statistic_lr = purrr::map2_dbl(.x = model_without_var, .y = original_fit,  ~likelihood_ratio( .x, .y )$lr_stat ) ) %>%
                                  mutate(p.value_lr = purrr::map2_dbl(.x = model_without_var,.y = original_fit, ~ likelihood_ratio( .x, .y )$p_value)) %>%
                                  mutate(p.value = p.value_lr),
                                .progress = "getting likelihood ratio p-values fit"))
  }else{
    combination_fits = combination_fits %>%
      mutate(results = purrr::map(.x = fit,
                                  ~ broom::tidy(.x) %>%
                                    filter(term != "(Intercept)" ) %>%
                                    left_join(multivariable_fit$lookup_vars_data, by = "term")%>%
                                    dplyr::rename("p.value_wald" = "p.value") %>%
                                    mutate(p.value = p.value_wald),
                                  .progress = "getting wald p-values fit"))
    }


  multivariable_fit$combination_fits <- combination_fits

  return(multivariable_fit)

}

#' Check confounding (STEP3)
#'
#' @param multivariable_fit object returned from multivariable_log_fit
#' @param sig_level significance level
#' @param delta_beta_threshold Beta delta threshold
#'
#' @return
#' @export
#'
#' @examples
check_confounding_for_excluded_in_step2 <- function(multivariable_fit, sig_level = 0.05, delta_beta_threshold = 20,
                                                    max_term_combinations = NULL) {


  # excluded variables from reduced models
  excluded_vars_step2 <- multivariable_fit$excluded_vars_step2
  if (is.null(max_term_combinations)) {
    max_term_combinations <- length(excluded_vars_step2)
  } else {
    max_term_combinations <- min(c(max_term_combinations, length(multivariable_fit$excluded_vars_step2)), na.rm = TRUE)
  }

  # compare reduced and larger model coefficients
  multivariable_fit$reduced_model_results <- multivariable_fit$reduced_model_results %>%
    left_join(multivariable_fit$larger_model_results[, c("term", "estimate")],
      by = "term", suffix = c("", "_larger_model")
    ) %>%
    mutate(coef_var = abs((estimate - estimate_larger_model) / estimate_larger_model) * 100) %>%
    mutate(being_confounded = case_when(
      coef_var > delta_beta_threshold ~ "yes",
      TRUE ~ "no"
    )) %>%
    dplyr::relocate(c("estimate_larger_model", "coef_var", "being_confounded"), .after = "estimate")


  # all possible formula combinations
  if (length(excluded_vars_step2) == 0 | all(multivariable_fit$reduced_model_results$being_confounded == "no")) {
    model_fit <- multivariable_fit$reduced_model_fit
    check_coef_variation <- NA
  } else {

    #get all combinations formula
    all_combinations_formulas <- formula_combination(vars = excluded_vars_step2,
      term_min = 1,
      term_limit = max_term_combinations,
      term_based = multivariable_fit$term_based,
      lookup_vars_data = multivariable_fit$lookup_vars_data
    )

    # store all fits and combinations
      combinations <- tibble::tibble(formula = c(
        sprintf(
          "%s + %s",
          stringr::str_squish(paste(deparse(as.formula(multivariable_fit$reduced_model_fit$formula)), collapse = "")),
          all_combinations_formulas
        ),
        stringr::str_squish(paste(deparse(as.formula(multivariable_fit$reduced_model_fit$formula)), collapse = ""))
      )) %>%
        mutate(added_term = c(all_combinations_formulas, NA)) %>%
        mutate(n_added_term = purrr::map_int(.x = added_term, .f = function(x) 1 + sum(stringr::str_count(x, "\\+")))) %>%
        mutate(fit = purrr::map(
          .x = formula, ~ try(speedglm::speedglm(formula = .x, family = binomial(), data = multivariable_fit$data)),
          .progress = "combinations_fit"
        )) %>%
        filter(!is.na(is.character(fit)))

      # store lr_pvalues
      if (multivariable_fit$sig_lr) {
        combinations <- combinations %>% mutate(results = purrr::map(
          .x = fit,
          ~ broom::tidy(.x) %>%
            mutate(original_fit = list(.x)) %>%
            filter(term != "(Intercept)") %>%
            left_join(multivariable_fit$lookup_vars_data, by = "term") %>%
            mutate(form_without_var = purrr::map2(.x = twin_terms, .y = original_fit, ~ update_formula(formula(.y), "-", paste0(.x, collapse = "- ") ))) %>%
            mutate(model_without_var = purrr::map(.x = unlist(form_without_var), ~ speedglm::speedglm(formula = .x, family = binomial(), data = multivariable_fit$data))) %>%
            mutate(statistic_lr = purrr::map2_dbl(.x = model_without_var, .y = original_fit, ~ likelihood_ratio(.x, .y)$lr_stat)) %>%
            mutate(p.value_lr = purrr::map2_dbl(.x = model_without_var, .y = original_fit, ~ likelihood_ratio(.x, .y)$p_value)),
          .progress = "getting likelihood ratio p-values fit"
        ))
      } else {
        # store wald_pvalues
        combinations <- combinations %>% mutate(results = purrr::map(
          .x = fit,
          ~ broom::tidy(.x) %>%
            mutate(original_fit = list(.x)) %>%
            filter(term != "(Intercept)") %>%
            left_join(multivariable_fit$lookup_vars_data, by = "term") %>%
            dplyr::rename("p.value_wald" = "p.value") %>%
            mutate(p.value = p.value_wald),
          .progress = "getting wald p-values fit"
        ))
      }

    #calculate the implications/metrics of adding a confounder
    check_coef_variation <- combinations %>%
      mutate(n_added_term_sig = purrr::map_int(
        .x = results, ~ nrow(.x %>% filter(p.value < sig_level & term %in% excluded_vars_step2)),
        .progress = "n_added_term_sig"
      )) %>%
      mutate(results_compare_reduced_model = purrr::map(
        .x = results,
        ~ multivariable_fit$reduced_model_results[, c("term", "estimate", "estimate_larger_model", "being_confounded")] %>%
          right_join(.x[, c("term", "estimate")],
            suffix = c("_reduced", "_reduced_plus_potential_confounder"),
            by = "term"
          ) %>%
          mutate(delta_beta = abs(estimate_reduced_plus_potential_confounder - estimate_larger_model) / estimate_larger_model * 100) %>%
          mutate(controlled = case_when(
            abs(delta_beta) < delta_beta_threshold & term != "(Intercept)" ~ "yes",
            TRUE ~ "no"
          )),
        .progress = "comparing larger to new fits"
      )) %>%
      mutate(controlled_vars = purrr::map_int(
        .x = results_compare_reduced_model, ~ sum(.x$controlled == "yes"),
        .progress = "set conrolled or no"
      )) %>%
      mutate(
        aic = purrr::map_dbl(.x = fit, ~ AIC(.x)),
        bic = purrr::map_dbl(.x = fit, ~ BIC(.x))
      ) %>%
      mutate_at("n_added_term", ~ case_when(
        is.na(.) ~ 0,
        TRUE ~ .
      )) %>%
      arrange(desc(controlled_vars), desc(n_added_term_sig), n_added_term, aic, bic) %>%
      mutate(selected = c("Yes", rep("No", nrow(.) - 1)))

    #choose the first model
    model_fit <- glm(
      formula = as.formula(check_coef_variation$formula[[1]]),
      family = "binomial",
      data = multivariable_fit$data
    )
  }

  #model performance
  model_fit_auc <- suppressMessages(pROC::auc(model_fit$y, model_fit$fitted.values))
  model_fit_acc <- suppressMessages(accuracy(model_fit))

  #final model formula
  model_fit_form <- stringr::str_remove_all(paste(deparse(model_fit$formula), collapse = ""), "\"") %>%
    stringr::str_squish()

  #store final results model

  if (multivariable_fit$sig_lr) {
    model_results <- model_fit %>%
      broom::tidy() %>%
      filter(term != "(Intercept)") %>%
      left_join(multivariable_fit$lookup_vars_data, by = "term") %>%
      mutate(form_without_var = purrr::map(.x = twin_terms, ~ as.formula(paste(model_fit_form, "- ", paste0(.x, collapse = " - "))))) %>%
      mutate(model_without_var = purrr::map(.x = unlist(form_without_var), ~ glm(formula = .x, family = "binomial", data = model_fit$data))) %>%
      mutate(statistic_lr = purrr::map_dbl(.x = model_without_var, ~ likelihood_ratio(.x, model_fit)$lr_stat)) %>%
      mutate(p.value_lr = purrr::map_dbl(.x = model_without_var, ~ likelihood_ratio(.x, model_fit)$p_value)) %>%
      dplyr::rename("p.value_wald" = "p.value") %>%
        dplyr::rename("p.value" = "p.value_lr")
  } else {
    model_results <- model_fit %>%
      broom::tidy() %>%
      filter(term != "(Intercept)") %>%
      left_join(multivariable_fit$lookup_vars_data, by = "term") %>%
      dplyr::relocate("variable", .after = "term") %>%
      mutate(p.value_wald = p.value)
  }


  #sotre added back coefficients and wald/lr comparisons for the current and past models
  if (any(excluded_vars_step2 %in% model_results$term)) {
    message(paste("STEP 3: check confouding from exluded vars in STEP 2: added",
                  excluded_vars_step2[which(excluded_vars_step2 %in% model_results$variable)],
                  "to the model \n"))

    added_back_coefs <- setdiff(names(coef(model_fit)),
                                names(coef(multivariable_fit$reduced_model_fit)))

    reduced_conf_wald_p_value <- multivariate_wald(
      fitted_model = model_fit,
      coefs_to_test = added_back_coefs
    )
    reduced_conf_lr_p_value <- likelihood_ratio(model_fit, multivariable_fit$reduced_model_fit)$p_value
  } else {
    message(paste("STEP 3: check confouding from exluded vars in STEP 2: none were added back to the model \n"))

    added_back_coefs <- NULL

    reduced_conf_wald_p_value <- NA
    reduced_conf_lr_p_value <- NA
  }


  return(list(
    "data" = multivariable_fit$data,
    "outcome_var" = multivariable_fit$outcome_var,
    "sig_lr" = multivariable_fit$sig_lr,
    "lookup" = multivariable_fit$lookup,
    "lookup_vars_data" = multivariable_fit$lookup_vars_data,
    "term_based" = multivariable_fit$term_based,
    "force_entry" = multivariable_fit$force_entry,
    "univariate_results" = multivariable_fit$univariate_results,
    "num_terms" = multivariable_fit$num_terms,
    "cat_terms" = multivariable_fit$cat_terms,

    # reduced
    "reduced_model_results" = multivariable_fit$reduced_model_results,
    "reduced_model_auc" = multivariable_fit$reduced_model_auc,
    "reduced_model_acc" = multivariable_fit$reduced_model_acc,

    # larger
    "larger_model_results" = multivariable_fit$larger_model_results,
    "larger_model_auc" = multivariable_fit$larger_model_auc,
    "larger_model_acc" = multivariable_fit$larger_model_acc,
    "excluded_vars_step1" = multivariable_fit$excluded_vars_step1,
    "excluded_vars_step2" = excluded_vars_step2,

    # confoudning
    "check_coef_variation" = check_coef_variation,
    "model_fit" = model_fit,
    "model_results" = model_results,
    "model_fit_auc" = model_fit_auc,
    "model_fit_acc" = model_fit_acc,
    "n_coef_added_back" =  length(added_back_coefs),

    # multivariate tests
    "complete_larger_lr_p_value" = multivariable_fit$complete_larger_lr_p_value,
    "complete_larger_wald_p_value" = multivariable_fit$complete_larger_wald_p_value,
    "larger_reduced_lr_p_value" = multivariable_fit$larger_reduced_lr_p_value,
    "larger_reduced_wald_p_value" = multivariable_fit$larger_reduced_wald_p_value,
    "reduced_conf_wald_p_value" = reduced_conf_wald_p_value,
    "reduced_conf_lr_p_value" = reduced_conf_lr_p_value
  ))
}
