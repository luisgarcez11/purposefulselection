
#' Multivariable fit (STEP 2)
#'
#' @param univariate_fit object returned by univariate_log_fit
#' @param all_variables
#' @param p_value_entrance entrance p-value
#' @param force_entry variables to be forced to entry
#' @param force_exclusion variables to be forced to be excluded
#' @param sig_level significance level
#' @param lookup_vars_data_interactions
#' @param sig_lr
#'
#' @return
#' @export
#'
#' @examples
multivariable_log_fit <- function( univariate_fit,
                                   all_variables = FALSE,
                                   p_value_entrance = 0.25,
                                   force_entry = NULL,
                                   force_exclusion = NULL,
                                   sig_level = 0.051,
                                   lookup_vars_data_interactions = FALSE,
                                   sig_lr = FALSE) {

  #force variables out
  if(!is.null(force_exclusion) | length(force_exclusion) == 0){
    if(!all(force_exclusion %in% names(univariate_fit$data))){stop("all force exclusion must be in the data")}
    univariate_fit$data <- univariate_fit$data %>%
      select(-all_of(force_exclusion))
    univariate_fit$univariate_results <- univariate_fit$univariate_results %>%
      filter(!term %in% force_exclusion)

  }

  #force entry certain vars
  if( length(force_entry) == 1 & !is.null(force_entry)){
    if(force_entry == "categorical" | force_entry == "numeric" ){
      if(force_entry == "categorical"){force_entry = univariate_fit$univariate_results %>% filter(type == "categorical") %>% pull(term)}else{
      if(force_entry == "numeric"){force_entry = univariate_fit$univariate_results %>% filter(type == "numeric") %>% pull(term)}
  }}}

  #if term-based is false, likelihood ratio is needed
  if(univariate_fit$term_based == FALSE){sig_lr = TRUE}

  #set p.value for decisions
  if(sig_lr){
    univariate_fit$univariate_results <- univariate_fit$univariate_results %>% mutate(p.value = p.value_lr)}else{
      univariate_fit$univariate_results <- univariate_fit$univariate_results %>% mutate(p.value = p.value_wald)
    }

  #included variables, filter by 0.25 (wald and lik likelihood might be different!)
  included_terms <- univariate_fit$univariate_results %>%
    filter(term != "(Intercept)") %>%
    filter(p.value < p_value_entrance | term %in% force_entry) %>%
    pull(term) %>% unique()

  #if all variables argument, join all variables
  if(all_variables){included_terms <- univariate_fit$dep_vars}

  #complete models
  all_vars <- unique(univariate_fit$univariate_results$term)

  if(any(is.na(all_vars))){all_vars <- all_vars[-which(is.na(all_vars))]}
  complete_model_form <- paste(univariate_fit$outcome_var, "~", paste(all_vars, collapse = " + "))
  suppressWarnings({complete_model_fit <- glm(formula = complete_model_form, family = "binomial", data = univariate_fit$data )})
  lookup_vars_data <- univariate_fit$lookup_vars_data


  #excluded variables
  excluded_vars_step1 <- univariate_fit$univariate_results %>%
    filter(!term %in% force_exclusion) %>%
    filter(p.value >= p_value_entrance & (!term %in% force_entry)) %>%
    pull(term) %>% unique()

  if(all_variables){excluded_vars_step1 <- c()}

  message(paste("STEP 1: UNIVARIATE ANALYSIS - excluded variables:", paste(excluded_vars_step1, collapse = " | "), "\n"))

  #force entry if needed
  if(!is.null(force_entry)){
    if(any(!force_entry %in% names(univariate_fit$data))){stop("force_entry is not in the data")}
    included_terms <- unique(c(included_terms, force_entry))
  }

  #formula
  larger_model_form <- paste(univariate_fit$outcome_var, "~", paste(included_terms, collapse = " + "))
  message(paste("STEP 2: LARGER MODEL - fit:", larger_model_form, "\n" ))


  #multivariable fit
  larger_model_fit <- suppressMessages(glm(formula = larger_model_form, family = "binomial", data = univariate_fit$data ))
  larger_model_auc <- suppressMessages(pROC::auc(larger_model_fit$y, larger_model_fit$fitted.values ))
  larger_model_acc <-  accuracy(larger_model_fit)

  #store larger model fit model
  results_tib <- broom::tidy(larger_model_fit) %>%
    dplyr::rename("p.value_wald" = "p.value") %>%
    mutate(or_estimate = exp(estimate),
           or_estimate_lower = exp(estimate - std.error),
           or_estimate_upper = exp(estimate + std.error),
           .after = estimate) %>%
    left_join(lookup_vars_data, by = "term") %>%
    relocate(variable, .after = "term")
  larger_model_results <- results_tib

  #likelihood_ratio test for all variables
  if (sig_lr) {
  larger_model_results <- larger_model_results %>%
    filter(term != "(Intercept)") %>%
    mutate(form_without_var = purrr::map(.x = twin_terms, ~ as.formula(paste(larger_model_form, "- ", paste(.x, collapse = " - ") )))) %>%
    mutate(model_without_var = purrr::map(.x = unlist(form_without_var), ~ glm(formula = .x, family = "binomial", data = univariate_fit$data))) %>%
    mutate(statistic_lr = purrr::map_dbl(.x = model_without_var, ~ likelihood_ratio(.x, larger_model_fit)$lr_stat)) %>%
    mutate(p.value_lr = purrr::map_dbl(.x = model_without_var, ~ likelihood_ratio(.x, larger_model_fit)$p_value)) %>%
    mutate(p.value = p.value_lr)
  }else{
    larger_model_results <- larger_model_results %>%
      mutate(p.value = p.value_wald)
  }

  #multivariate test
  excluded_coefs <- setdiff(names(coef(complete_model_fit)), c(names(coef(larger_model_fit)), "(Intercept)" ))

  if(length(excluded_coefs) > 0){

      complete_larger_wald_p_value = try(as.numeric(multivariate_wald(fitted_model = complete_model_fit,
                                                                   coefs_to_test = excluded_coefs)))

      complete_larger_lr_p_value = lmtest::lrtest(complete_model_fit, larger_model_fit)["Pr(>Chisq)"][2,]

    }else{

    complete_larger_wald_p_value = NA
    complete_larger_lr_p_value = NA


  }


  #excluded_variables
  excluded_vars_step2 <- larger_model_results %>%
    filter(term  != "(Intercept)" ) %>%
    filter(p.value > sig_level & (!term %in% force_entry) ) %>%
    pull(term) %>%
    unique()

  message(paste("removing", paste(excluded_vars_step2, collapse = " & "),
                "from larger model. Significance level higher than:",
                sig_level, "\n"))

  #reduced model fit
  if(length(excluded_vars_step2) > 0){
    reduced_model_form <- as.formula(paste(larger_model_form, "-", paste(excluded_vars_step2, collapse = "-")))
    reduced_model_fit <- update(larger_model_fit, reduced_model_form)

    reduced_model_form_string <- stringr::str_squish(paste(deparse(reduced_model_fit$formula), collapse = ""))

  }else{
    reduced_model_form <- larger_model_form
    reduced_model_fit <- larger_model_fit

    reduced_model_form_string <- reduced_model_form
  }

  #performance results for reduced model
  reduced_model_auc <-  suppressMessages(pROC::auc(reduced_model_fit$y, reduced_model_fit$fitted.values ))
  reduced_model_acc <-  suppressMessages(accuracy(reduced_model_fit))

  #results for reduced model
  reduced_model_results <- broom::tidy(reduced_model_fit) %>%
    dplyr::rename("p.value_wald" = "p.value") %>%
    mutate(or_estimate = exp(estimate),
           or_estimate_lower = exp(estimate - std.error),
           or_estimate_upper = exp(estimate + std.error),
           .after = estimate) %>%
    left_join(lookup_vars_data, by = "term") %>%
    relocate(variable, .after = "term")%>%
    filter(term != "(Intercept)")

  #likelihood p-values for reduced model
  if(sig_lr){
    reduced_model_results = reduced_model_results %>%
      mutate(form_without_var = purrr::map(.x = twin_terms, ~as.formula(paste(reduced_model_form_string, "- ", paste0(.x, collapse = "- ") )))) %>%
      mutate(model_without_var = purrr::map(.x = unlist(form_without_var), ~glm(formula = .x, family = "binomial", data = univariate_fit$data ))) %>%
      mutate(statistic_lr = purrr::map_dbl(.x = model_without_var, ~ likelihood_ratio( .x, reduced_model_fit )$lr_stat)) %>%
      mutate(p.value_lr = purrr::map_dbl(.x = model_without_var, ~ likelihood_ratio( .x, reduced_model_fit )$p_value)) %>%
      mutate(p.value = p.value_lr)
    }else{
      reduced_model_results = reduced_model_results %>%
        mutate(p.value = p.value_wald)
      }

  #multivariate test
  excluded_coefs <- setdiff(names(coef(larger_model_fit)), c(names(coef(reduced_model_fit)), "(Intercept)" ))

  if(length(excluded_coefs) > 0){
    larger_reduced_wald_p_value <- try(as.numeric(multivariate_wald(fitted_model = larger_model_fit,
                                                               coefs_to_test = excluded_coefs)))
    larger_reduced_lr_p_value <- lmtest::lrtest(larger_model_fit, reduced_model_fit)["Pr(>Chisq)"][2,]}else{
    larger_reduced_wald_p_value = NA
    larger_reduced_lr_p_value = NA
    }


  message(paste("STEP 2: REDUCED MODEL: ", paste(deparse(reduced_model_form, width.cutoff = 500), collapse=""), "\n"))

  return(list("data" = univariate_fit$data,
              "outcome_var" = univariate_fit$outcome_var,
              "term_based" = univariate_fit$term_based ,
              "force_entry" = force_entry,
              "sig_lr" = sig_lr,
              "univariate_results" = univariate_fit$univariate_results,
              "lookup" = univariate_fit$lookup,
              "num_terms" = univariate_fit$num_terms,
              "cat_terms" = univariate_fit$cat_terms,
              "included_terms" = included_terms,

              "excluded_vars_step1" = excluded_vars_step1,
              "larger_model_fit" = larger_model_fit,
              "larger_model_results" = larger_model_results,
              "larger_model_auc" = larger_model_auc,
              "larger_model_acc" = larger_model_acc,
              "excluded_vars_step2" = excluded_vars_step2,
              "reduced_model_fit" = reduced_model_fit,
              "reduced_model_results" = reduced_model_results,
              "reduced_model_auc" = reduced_model_auc,
              "reduced_model_acc" = reduced_model_acc,
              "reduced_model_form" = stringr::str_remove_all(paste(deparse(reduced_model_form), collapse = ""), "\"" ) %>%
                stringr::str_squish(),

              #utils
              "lookup_vars_data" = lookup_vars_data,

              #multivariate tests
              "complete_larger_lr_p_value" = complete_larger_lr_p_value,
              "complete_larger_wald_p_value" = complete_larger_wald_p_value,
              "larger_reduced_lr_p_value" = larger_reduced_lr_p_value,
              "larger_reduced_wald_p_value" = larger_reduced_wald_p_value))

}
