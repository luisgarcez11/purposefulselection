#' Check Interactions
#'
#' @param fit_to_test
#' @param term_limit
#' @param sig_level
#'
#' @return
#' @export
#'
#' @examples
check_interactions <- function(fit_to_test, term_limit = 2, sig_level = 0.10) {
  # vars to be checked for interactions
  vars <- names(fit_to_test$model_fit$model)
  message(paste("STEP 6: check interactions for", paste(vars, collapse = ", "), "\n"))

  all_combinations_formulas <- formula_combination(setdiff(vars, fit_to_test$outcome_var),
    term_min = 2,
    term_limit = term_limit,
    operator = ":"
  )

  # store all fits and combinations
  combined_fits_interactions <- tibble::tibble(formula = c(
    fit_to_test$model_fit$formula,
    sprintf(
      "%s + %s",
      stringr::str_squish(paste(deparse(as.formula(fit_to_test$model_fit$formula)), collapse = "")),
      all_combinations_formulas
    )
  )) %>%
    mutate(added_variable = c("", all_combinations_formulas)) %>%
    mutate(fit = purrr::map(.x = formula, ~ glm(formula = .x, family = "binomial", data = fit_to_test$data))) %>%
    mutate(results = purrr::map(.x = fit, ~ broom::tidy(.x))) %>%
    mutate(log_likelihood = purrr::map_dbl(.x = fit, ~ stats::logLik(.x))) %>%
    mutate(lr_test = purrr::map(
      .x = fit,
      ~ likelihood_ratio(fit_to_test$model_fit, .x)
    )) %>%
    tidyr::unnest_wider(lr_test)

  # fit multivariate interactions
  added_interaction_terms <- combined_fits_interactions %>%
    filter(p_value < sig_level) %>%
    pull(added_variable)

  message(paste("STEP 6: ADDED INTERACTIONS: interactions", paste(added_interaction_terms, collapse = ", "), "\n"))

  # if any combinations...
  if (length(added_interaction_terms) > 0) {
    interaction_formula <- paste(stringr::str_squish(paste(deparse(as.formula(fit_to_test$model_fit$formula)), collapse = "")),
                                 paste(added_interaction_terms, collapse = " + "), sep = " + ")
    multiinteraction_model <- glm(
      formula = interaction_formula,
      family = "binomial", data = fit_to_test$data
    )

    multiinteraction_model_results <- broom::tidy(multiinteraction_model)

    lookup_interactions <- tail(multiinteraction_model_results$term, length(added_interaction_terms))
    names(lookup_interactions) <- added_interaction_terms
    index_interaction <- which(multiinteraction_model_results$term %in% lookup_interactions[added_interaction_terms])

    stepwise_fits <- tibble::tibble(
      formula = interaction_formula,
      fit = list(multiinteraction_model),
      results = list(multiinteraction_model_results)
    ) %>%
      mutate(lr_test = purrr::map(
        .x = fit,
        ~ likelihood_ratio(fit_to_test$model_fit, .x)
      )) %>%
      tidyr::unnest_wider(lr_test)

    while (any(multiinteraction_model_results$p.value[index_interaction] > sig_level, na.rm = TRUE)) {

      max_p_value <- max(multiinteraction_model_results$p.value[index_interaction] )
      interaction_to_remove_index <- which(multiinteraction_model_results$p.value[index_interaction] == max_p_value)
      interaction_to_remove <- added_interaction_terms[interaction_to_remove_index]

      multiinteraction_model <- update(multiinteraction_model,
                                       drop.terms(multiinteraction_model$terms,
                                                  grep(interaction_to_remove, attr(multiinteraction_model$terms, "term.labels")),
                                                  keep.response = TRUE))

      multiinteraction_model_results <- broom::tidy(multiinteraction_model)
      stepwise_fits <- stepwise_fits %>%
        bind_rows(
          tibble(
            formula = stringr::str_squish(paste(deparse(multiinteraction_model$formula), collapse = "")),
            fit = list(multiinteraction_model),
            results = list(multiinteraction_model_results)
          ) %>%
            mutate(lr_test = purrr::map(.x = fit, ~ likelihood_ratio(fit_to_test$model_fit, .x))) %>%
            tidyr::unnest_wider(lr_test)
        )
      message(paste("interaction", paste(interaction_to_remove, collapse = ", "), "removed" ))
    }
  } else {
    message("no interactions added \n")
  }

  multiinteraction_model_form <- stringr::str_squish(paste(deparse(multiinteraction_model$formula), collapse = ""))


  #results for reduced model
  multiinteraction_model_results <- multiinteraction_model_results %>%
    mutate(or_estimate = exp(estimate),
           or_estimate_lower = exp(estimate - std.error),
           or_estimate_upper = exp(estimate + std.error),
           .after = estimate) %>%
    mutate(variable = purrr::map_chr(.x = term, ~fit_to_test$lookup[.x] ),
           .after = "term") %>%
    mutate(variable = case_when(is.na(variable) ~ stringr::str_remove_all(term,  paste(c(multiinteraction_model$xlevels %>% unlist()), collapse = "|") ) ,
                                TRUE ~ variable)) %>%
    filter(term != "(Intercept)") %>%
    mutate(form_without_var = purrr::map(.x = variable, ~as.formula(paste(multiinteraction_model_form, "- ", .x)))) %>%
    mutate(model_without_var = purrr::map(.x = unlist(form_without_var), ~glm(formula = .x, family = "binomial", data = fit_to_test$data ))) %>%
    mutate(statistic_lr = purrr::map_dbl(.x = model_without_var, ~ likelihood_ratio( .x, multiinteraction_model )$lr_stat)) %>%
    mutate(p.value_lr = purrr::map_dbl(.x = model_without_var, ~ likelihood_ratio( .x, multiinteraction_model )$p_value))


  return(list(
    "data" = fit_to_test$data,
    "outcome_var" = fit_to_test$outcome_var,
    "lookup" = fit_to_test$lookup,
    "univariate_results" = fit_to_test$univariate_results,
    "reduced_model_results" = fit_to_test$reduced_model_results,
    "larger_model_results" = fit_to_test$larger_model_results,
    "check_coef_variation" = fit_to_test$check_coef_variation,
    "check_vars_excluded_step1" = fit_to_test$check_vars_excluded_step1,
    "excluded_vars_step1" = fit_to_test$excluded_vars_step1,
    "excluded_vars_step2" = fit_to_test$excluded_vars_step2,
    "added_interactions" = combined_fits_interactions,
    "added_interactions_history" = stepwise_fits,
    "model_fit" = multiinteraction_model,
    "model_results" = multiinteraction_model_results
  ))
}
