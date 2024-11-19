
#' Comparing History term based and var based history objects
#'
#' @param history_A
#' @param history_B
#' @param var_based_comparison
#'
#' @return
#' @export
#'
#' @examples
compare_history_objects <- function(history_A, history_B,
                                    var_based_comparison = FALSE,
                                    suffix = c("_A", "_B")){

  #acess history
  history_A <- history_A[["history"]]
  history_B <- history_B[["history"]]

  #final desc name
  final_desc_name_A <- case_when("final_desc" %in% names(history_A) ~ "final_desc",
                                 "final_desc_run2" %in% names(history_A) ~ "final_desc_run2" )

  final_desc_name_B <- case_when("final_desc" %in% names(history_B) ~ "final_desc",
                                 "final_desc_run2" %in% names(history_B) ~ "final_desc_run2"  )


  #joined history objects
  if(var_based_comparison){
  vars_different <- history_A[,c("term", "variable", "type", final_desc_name_A )] %>%
    mutate(term = case_when(is.na(term) ~ variable,
                            TRUE ~ term)) %>%
    dplyr::rename("final_desc" = final_desc_name_A)  %>%
    mutate(final_desc = case_when(is.na(final_desc) & type %in% c("numeric", "categorical") ~ "out",
                                       TRUE ~ final_desc))%>%
    rename_all( ~paste0(., "_A")) %>%
    left_join(history_B[,c("variable", final_desc_name_B )] %>% distinct() %>%
                                 dplyr::rename("final_desc" = final_desc_name_B) %>%
                                 dplyr::rename_all( ~paste0(., "_B")),
                               by = c("variable_A" = "variable_B"))
  }else{

    vars_different <- history_A[,c("term", "variable", "type", final_desc_name_A )] %>%
      mutate(term = case_when(is.na(term) ~ variable,
                              TRUE ~ term)) %>%
      dplyr::rename("final_desc" = final_desc_name_A)  %>%
      mutate(final_desc = case_when(is.na(final_desc) & type %in% c("numeric", "categorical") ~ "out",
                                    TRUE ~ final_desc))%>%
      rename_all( ~paste0(., "_A")) %>%
      left_join(history_B[,c("term","variable", final_desc_name_B )] %>%
                  mutate(term = case_when(is.na(term) ~ variable,
                                          TRUE ~ term)) %>%
                  dplyr::rename("final_desc" = final_desc_name_B) %>%
                  dplyr::rename_all( ~paste0(., "_B")),
                by = c("term_A" = "term_B", "variable_A" = "variable_B"))

}



  #performance measure differences
  auc_term_based = as.numeric(vars_different[["final_desc_A"]][vars_different$variable_A == "AUC"])
  acc_term_based = as.numeric(vars_different[["final_desc_A"]][vars_different$variable_A == "Accuracy"])

  auc_var_based = as.numeric(vars_different[["final_desc_B"]][vars_different$variable_A == "AUC"])
  acc_var_based = as.numeric(vars_different[["final_desc_B"]][vars_different$variable_A == "Accuracy"])

  auc_diff = as.numeric(auc_term_based) - as.numeric(auc_var_based)
  acc_diff = as.numeric(acc_term_based) - as.numeric(acc_var_based)


  #measure dissimilarities between variables
  levels_between_history <- vars_different %>%
    filter(!variable_A %in% c("AUC", "Accuracy",
                                      "multivariate_lr", "multivariate_wald", "n_added_back",
                                      "outcome_proportion")) %>%
  mutate_at(c("final_desc_A", "final_desc_B"), ~case_when(is.na(.) ~ "out",
                                     TRUE ~ .)) %>%
    group_by(variable_A, type_A) %>%
    summarise(n_initial_params = n(),
              k_initial_cat_params = sum(type_A == "categorical"),
              k_initial_num_params = sum(type_A == "numeric"),
              k_in_cat_params_A = sum(final_desc_A == "in" & type_A == "categorical"),
              k_in_cat_params_B = sum(final_desc_B == "in" & type_A == "categorical"),
              k_in_num_params_A = sum(final_desc_A == "in" & type_A == "numeric"),
              k_in_num_params_B = sum(final_desc_B == "in" & type_A == "numeric")
              ) %>%
    ungroup() %>%
    mutate(same_decision_cat =  (k_in_cat_params_A > 0 & k_in_cat_params_B > 0) |
             (k_in_cat_params_A == k_in_cat_params_B),
           different_retained_levels_cat = k_in_cat_params_A != k_in_cat_params_B,
           level_differences_cat = abs(k_in_cat_params_A - k_in_cat_params_B))


  #summary of these measure
  levels_between_history_summary = levels_between_history %>%
    filter(type_A == "categorical") %>%
    group_by() %>%
    summarise(k_initial_params_cat = sum(k_initial_cat_params),
              n_initial_vars_cat = sum(type_A == "categorical"),
              same_decision_cat = sum(same_decision_cat),
              different_retained_levels_cat = sum(different_retained_levels_cat),
              level_differences_cat = sum (level_differences_cat),
              ) %>%
    ungroup() %>%
    mutate(same_decision_pct = same_decision_cat / n_initial_vars_cat * 100,
           different_retained_levels_pct = different_retained_levels_cat/n_initial_vars_cat * 100,
           level_differences_pct = level_differences_cat/ k_initial_params_cat * 100) %>%
  bind_cols(tibble("auc_diff" = auc_diff, "acc_diff" = acc_diff)) %>%
    mutate(dissimilarity = (100-same_decision_pct)*3 + (different_retained_levels_pct)*2 + (level_differences_pct))

  return(list("vars_different" = vars_different ,
              "levels_between_history" = levels_between_history,
              "levels_between_history_summary" = levels_between_history_summary))

}
