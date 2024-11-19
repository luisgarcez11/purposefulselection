#' Describe data
#'
#' @param history_object dato to be described
#'
#' @return list with percentage of categorical-numeric, %variables selected,
#' %variables selected per category, number of variables in/out
#' @export
#'
#' @examples
describe_history <- function(history_object){

  larger_desc = case_when("larger_model_desc_run1" %in% names(history_object) ~ "larger_model_desc_run1",
                         "larger_model_desc" %in% names(history_object) ~ "larger_model_desc")

  check_confounding_desc = case_when("check_confounding_desc_run2" %in% names(history_object) ~ "check_confounding_desc_run2",
                          "check_confounding_desc" %in% names(history_object) ~ "check_confounding_desc")

  final_desc = case_when("final_desc" %in% names(history_object) ~ "final_desc",
                         "final_desc_run2" %in% names(history_object) ~ "final_desc_run2")

  #imbalance categorical info
  if(any(history_object$type == "categorical", na.rm = TRUE)){
    if("n_occurr_run1" %in% names(history_object)){history_object <- history_object %>% rename("n_occurr" = "n_occurr_run1")}

  results_by_variable <- history_object %>%
    filter(type == "categorical") %>%
    group_by(variable) %>%
    summarise(gini_index = 1-sum(prop.table( n_occurr)^2),
              entropy = entropy::entropy(prop.table( n_occurr)),
              n_terms_initial = n(),
              n_terms_retained = sum(final_desc == "in" & type == "categorical", na.rm = TRUE)) %>%
    ungroup()
  }else{
    results_by_variable <- NULL
}

  results_global = tibble::tibble(
    "n_terms_initial" = nrow(history_object) - sum(is.na(history_object$type)),
    "n_terms_total_final" = sum(history_object[[final_desc]] == "in" , na.rm = TRUE),
    "in_terms_final_pct" = sum(history_object[[final_desc]] == "in", na.rm = TRUE)/ n_terms_total_final,

    "n_terms_num_initial" = sum(history_object$type == "numeric", na.rm = TRUE),
    "n_terms_num_final" = sum(history_object[[final_desc]] == "in" & history_object$type == "numeric", na.rm = TRUE),
    "num_in_terms_final_pct" = sum(history_object[[final_desc]] == "in" & history_object$type == "numeric", na.rm = TRUE)/n_terms_num_initial,

    "n_terms_cat_initial" = sum(history_object$type == "categorical", na.rm = TRUE),
    "n_terms_cat_final" = sum(history_object[[final_desc]] == "in" & history_object$type == "categorical", na.rm = TRUE),
    "cat_in_terms_pct" = sum(history_object[[final_desc]] == "in" & history_object$type == "categorical", na.rm = TRUE)/n_terms_cat_initial,

    "num_in_skewness_mean" = mean(abs(history_object$skew)[history_object[[final_desc]] == "in" & history_object$type == "numeric"], na.rm = TRUE),
    "num_out_skewness_mean" = mean(abs(history_object$skew)[history_object[[final_desc]] == "out" & history_object$type == "numeric"], na.rm = TRUE),
    "num_in_kurtosis_mean" = mean(abs(history_object$kurtosis)[history_object[[final_desc]] == "in" & history_object$type == "numeric"], na.rm = TRUE),
    "num_out_kurtosis_mean" = mean(abs(history_object$kurtosis)[history_object[[final_desc]] == "out" & history_object$type == "numeric"], na.rm = TRUE),
    "num_in_se_mean" = mean(abs(history_object$se)[history_object[[final_desc]] == "in" & history_object$type == "numeric"], na.rm = TRUE),
    "num_out_se_mean" = mean(abs(history_object$se)[history_object[[final_desc]] == "out" & history_object$type == "numeric"], na.rm = TRUE),
    # "cat_in_largshare_mean" = mean(abs(history_object$largest_share)[history_object[[final_desc]] == "in" & history_object$type == "categorical"], na.rm = TRUE),
    # "cat_out_largshare_mean" = mean(abs(history_object$largest_share)[history_object[[final_desc]] == "out" & history_object$type == "categorical"], na.rm = TRUE),

    #performace
    "auc_initial" = ifelse(!is.na(larger_desc), history_object[[larger_desc]][which(history_object$variable ==  "AUC")], NA),
    "accuracy_initial" = ifelse(!is.na(larger_desc),history_object[[larger_desc]][which(history_object$variable ==  "Accuracy")],NA),
    "auc" = history_object[[final_desc]][which(history_object$variable ==  "AUC")],
    "accuracy" = history_object[[final_desc]][which(history_object$variable ==  "Accuracy")],
    "outcome_proportion" = history_object[[final_desc]][which(history_object$variable ==  "outcome_proportion")],
    "n_coef_added_back_conf" = ifelse(!is.na(check_confounding_desc),history_object[[check_confounding_desc]][which(history_object$variable ==  "n_added_back")], NA)


  ) %>%
    mutate_at(vars(ends_with("pct")), ~round(.*100,1))



  return(list(results_global = results_global,
              results_by_variable = results_by_variable))
}
