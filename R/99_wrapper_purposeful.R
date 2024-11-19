
#' Wrap the all workflow of the selection process
#'
#' @param data data to be analyzed
#' @param dep_vars dependent variables
#' @param outcome_var outcome variables
#' @param multivariate_uni multivariate for STEP1
#' @param pe entrance p-value
#' @param sig_level significance level
#' @param delta_beta_threshold dela beta threshold
#'
#' @return
#' @export
#'
#' @examples
wrapper_purposeful <- function(data, dep_vars, outcome_var,
                               multivariate_uni = FALSE,
                               pe = 0.25,
                               sig_level = 0.05,
                               delta_beta_threshold = 20,
                               force_entry_ = NULL,
                               max_term_combinations = NULL,
                               term_based = FALSE,
                               sig_lr = FALSE,
                               return_desc = FALSE
                               ){

  #remove all missing data
  data <- data[,c(dep_vars, outcome_var)] %>% stats::na.omit()

  history_object <- univariate_log_fit(.data = data, outcome_var = outcome_var, dep_vars = dep_vars, term_based = term_based, multivariate = multivariate_uni) %>%
    multivariable_log_fit(p_value_entrance = pe, sig_level = sig_level, sig_lr = sig_lr, force_entry = force_entry_) %>% #reduced model
    check_confounding_for_excluded_in_step2( delta_beta_threshold = delta_beta_threshold, max_term_combinations = max_term_combinations) %>% #preliminary main effects model
    check_significance_for_excluded_in_step1( sig_level = sig_level) %>%
    history()

  if(!return_desc){
    return(history_object)}else{
      history_object[["history_desc"]] = history_object[["history"]] %>% describe_history()
    return(history_object)}


}


#' Wrap 2-step workflow of the selection process
#'
#' @param data
#' @param dep_vars
#' @param outcome_var
#' @param force_entry_step1
#' @param force_exclusion_step1
#' @param force_entry_step2
#' @param force_exclusion_step2
#' @param p_value_entrance
#' @param sig_level
#' @param delta_beta_threshold
#'
#' @return
#' @export
#'
#' @examples
wrapper_purposeful_2step <- function(data, dep_vars, outcome_var,
                               p_value_entrance = 0.25,
                               sig_level = 0.05,
                               force_entry_step1 = NULL,
                               force_exclusion_step1 = NULL,
                               force_entry_step2 = NULL,
                               force_exclusion_step2 = NULL,
                               delta_beta_threshold = 20,
                               max_term_combinations = NULL,
                               term_based = FALSE,
                               sig_lr = FALSE,
                               return_history = FALSE
){

  #remove all missing data
  data <- data[,c(dep_vars, outcome_var)] %>% stats::na.omit()

  #step1
  step1 <- univariate_log_fit(.data = data, outcome_var = outcome_var, dep_vars = dep_vars, term_based = term_based) %>%
    multivariable_log_fit(p_value_entrance = p_value_entrance,
                          sig_level = sig_level,
                          sig_lr = sig_lr,
                          force_entry = force_entry_step1,
                          force_exclusion = force_exclusion_step1) %>% #reduced model
    check_confounding_for_excluded_in_step2(delta_beta_threshold = delta_beta_threshold,
                                            max_term_combinations = max_term_combinations) %>% #preliminary main effects model
    check_significance_for_excluded_in_step1(sig_level = sig_level)

  #history step1
  history_step1 <- history(step1)


  #step2
  step2 <- step1 %>%
    multivariable_log_fit(p_value_entrance = p_value_entrance,
                          sig_level = sig_level,
                          sig_lr = sig_lr,
                          force_entry = step1$num_terms_included,
                          force_exclusion = step1$num_terms_excluded) %>% #reduced model
    check_confounding_for_excluded_in_step2(delta_beta_threshold = delta_beta_threshold, max_term_combinations = max_term_combinations) %>% #preliminary main effects model
    check_significance_for_excluded_in_step1(sig_level = sig_level)

  #history step2
  history_step2 <- history(step2)

  #whole history
  by_vars <- c("term", "variable", "skew" ,"kurtosis", "se" , "type", "largest_share", "univariate" )
  by_vars <- by_vars[which(by_vars %in% intersect(names(history_step1), names(history_step2)))]

  history_object <- full_join(history_step1,
                           history_step2,
                           by = by_vars,
                           suffix = c("_run1", "_run2"))


  if(return_history){return(history_object)}else{return(history_object %>% describe_history())}


}
