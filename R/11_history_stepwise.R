
#' History fit
#'
#' @param fit object returned from stepwise
#'
#' @return
#' @export
#'
#' @examples
history_stepwise <- function(fit,
                             return_model = FALSE){



  message("compiling history...")

  variables_included <- unique(fit$final_model_results$term[-1])
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



  if(fit$direction == "forward"){


    stepwise_history <- as_tibble(fit$p_value_matrix_to_include, rownames = "variable") %>%
      left_join(as_tibble(fit$p_value_matrix_to_exclude, rownames = "variable"), by = "variable", suffix = c("_in", "_out")) %>%
      janitor::remove_empty(which = "cols") %>%
      select(variable, order(colnames(.)))


    steps_main <- unique(stringr::str_extract_all(names(stepwise_history)[-which(names(stepwise_history) == "variable")], pattern = "\\d")) %>% unlist()
    number_steps <- max(as.numeric(steps_main))

    if(!is.null(fit$p_value_matrix_to_exclude_interactions)){
    #add interactions
    stepwise_history_int <- as_tibble(fit$p_value_matrix_to_exclude_interactions, rownames = "variable") %>%
      left_join(as_tibble(fit$p_value_matrix_to_include_interactions, rownames = "variable"), by = "variable", suffix = c("_out", "_in")) %>%
      janitor::remove_empty(which = "cols") %>%
      select(variable, order(colnames(.))) %>%
      rename_if(is.double, ~paste0(as.character(as.numeric(stringr::str_extract(pattern = "\\d+", .)) + number_steps), "_",stringr::str_extract(pattern = "[:alpha:]+", .), ""))

    stepwise_history <- stepwise_history %>% full_join(stepwise_history_int)}


    steps <- unique(stringr::str_extract_all(names(stepwise_history)[-which(names(stepwise_history) == "variable")], pattern = "\\d+")) %>% unlist()
    steps <- sort(as.numeric(steps))

    for(step in steps){

      col_in <- paste0(step, "_in")
      col_out <- paste0(step, "_out")
      colname <- paste0(step, "_desc")

      if(!col_in %in% names(stepwise_history)){stepwise_history[[col_in]] <- NA; stepwise_history <-  stepwise_history %>%
        dplyr::relocate(all_of(col_in), .before = col_out)}
      if(!col_out %in% names(stepwise_history)){stepwise_history[[col_out]] <- NA; stepwise_history <-  stepwise_history %>%
        dplyr::relocate(all_of(col_out), .after = col_in)}

      stepwise_history[[colname]] <- case_when( ((stepwise_history[[col_out]] == suppressWarnings(max(stepwise_history[[col_in]], na.rm=TRUE))) &
                                                   (stepwise_history[[col_out]]  > fit$p_value_removal) ) ~ "out",
                 ((stepwise_history[[col_in]] == suppressWarnings(min(stepwise_history[[col_in]], na.rm=TRUE)) ) &
                                                   (stepwise_history[[col_in]]  < fit$p_value_entrance) & (stepwise_history[[col_out]]  < fit$p_value_removal)) ~ "in",
                (is.na(stepwise_history[[col_in]]) & (stepwise_history[[col_out]]  < fit$p_value_removal)) ~ "in",
                (is.na(stepwise_history[[col_in]]) & (stepwise_history[[col_out]]  > fit$p_value_removal)) ~ "out",
                (is.na(stepwise_history[[col_in]]) & is.na(stepwise_history[[col_out]]) ) ~ NA_character_,
                TRUE ~ "out")


      stepwise_history <- stepwise_history %>%
        dplyr::relocate(all_of(colname), .after = col_out)


    }

  }

  if(fit$direction == "backward"){


    stepwise_history <- as_tibble(fit$p_value_matrix_to_exclude, rownames = "variable") %>%
      left_join(as_tibble(fit$p_value_matrix_to_include, rownames = "variable"), by = "variable", suffix = c("_out", "_in")) %>%
      janitor::remove_empty(which = "cols") %>%
      select(variable, order(colnames(.)))

    steps_main <- unique(stringr::str_extract_all(names(stepwise_history)[-which(names(stepwise_history) == "variable")], pattern = "\\d+")) %>% unlist()
    number_steps <- max(as.numeric(steps_main))

    if(!is.null(fit$p_value_matrix_to_exclude_interactions)){
    #add interactions
    stepwise_history_int <- as_tibble(fit$p_value_matrix_to_exclude_interactions, rownames = "variable") %>%
      left_join(as_tibble(fit$p_value_matrix_to_include_interactions, rownames = "variable"), by = "variable", suffix = c("_out", "_in")) %>%
      janitor::remove_empty(which = "cols") %>%
      select(variable, order(colnames(.))) %>%
      rename_if(is.double, ~paste0(as.character(as.numeric(stringr::str_extract(pattern = "\\d+", .)) + number_steps), "_",stringr::str_extract(pattern = "[:alpha:]+", .), ""))

    stepwise_history <- stepwise_history %>% full_join(stepwise_history_int)
}

    steps <- sort(as.numeric(unique(stringr::str_extract_all(names(stepwise_history)[-which(names(stepwise_history) == "variable")], pattern = "\\d+")) %>% unlist()))

    for(step in steps){

      colname_former <- paste0(as.numeric(step) - 1, "_desc")

      col_in <- paste0(step, "_in")
      col_out <- paste0(step, "_out")
      colname <- paste0(step, "_desc")


      if(!col_in %in% names(stepwise_history)){stepwise_history[[col_in]] <- NA}
      if(!col_out %in% names(stepwise_history)){stepwise_history[[col_in]] <- NA}


      stepwise_history[[colname]] <- case_when( ((stepwise_history[[col_in]] == suppressWarnings(min(stepwise_history[[col_in]], na.rm=TRUE))) & (stepwise_history[[col_in]]  < fit$p_value_entrance) ) ~ "in",
                                                ((stepwise_history[[col_out]] == max(stepwise_history[[col_out]], na.rm=TRUE)) & (stepwise_history[[col_out]]  > fit$p_value_removal) & (stepwise_history[[col_out]]  > fit$p_value_removal)) ~ "out",
                                                (!is.na(stepwise_history[[col_out]]) & (stepwise_history[[col_in]]  < fit$p_value_entrance)) ~ "in",
                                                (is.na(stepwise_history[[col_out]]) & (stepwise_history[[col_in]]  > fit$p_value_removal)) ~ "out",
                                                (is.na(stepwise_history[[col_in]]) & is.na(stepwise_history[[col_out]]) ) ~ NA_character_,
                                                TRUE ~ "in")

      if(colname_former != "0_desc"){

        stepwise_history <- stepwise_history %>%
          dplyr::relocate(all_of(col_out), .after = colname_former)

      }

      stepwise_history <- stepwise_history %>%
        dplyr::relocate(all_of(col_in), .after = col_out) %>%
        dplyr::relocate(all_of(colname), .after = col_in)


    }
  }


  auc_history <- fit$auc_history %>%
    as_tibble() %>%
    dplyr::rename_all(~paste0(., "_desc")) %>%
    janitor::remove_empty(which = "cols") %>%
    mutate_all(~as.character(round(., 2)))

  acc_history <- fit$acc_history %>%
    as_tibble() %>%
    dplyr::rename_all(~paste0(., "_desc")) %>%
    janitor::remove_empty(which = "cols") %>%
    mutate_all(~as.character(round(., 2)))

  #add AUC and accuracy
  stepwise_history <- stepwise_history %>%
    mutate_if(is.double, ~round(., 2)) %>%
    add_row(variable = "AUC", auc_history) %>%
    add_row(variable = "Accuracy", acc_history)

  final = stepwise_history[,ncol(stepwise_history)]
  names(final) = "final_desc"
  stepwise_history <- bind_cols(stepwise_history, final)

  #add outcome proportion
  stepwise_history <- stepwise_history %>%
    add_row(variable = "outcome_proportion",
            final_desc = as.character(round(max(prop.table(table(fit$data[[fit$outcome_var]]))), 2))
    )

  #add term variable & description
  stepwise_history <- stepwise_history %>%
    rename("term" = "variable") %>%
    mutate(type = case_when(term %in% num_terms ~ "numeric",
                            term %in% cat_terms ~ "categorical")) %>%
    relocate("type", .after = "term") %>%
    left_join(fit$final_model_results[,c("variable", "term")], by = "term") %>%
    mutate(variable = case_when(term %in% c("AUC", "Accuracy", "outcome_proportion") ~ term,
                                TRUE ~ variable)) %>%
    mutate(n_occurr = purrr::map2_int(.x = term, .y = type, .f = ~case_when(.y == "categorical" ~ sum(fit$data[[.x]])))) %>%
    mutate(share = purrr::map2_dbl(.x = term,
                                   .y = type,
                                   .f = ~case_when(.y == "categorical" ~ sum(fit$data[[.x]])/nrow(fit$data)  ) )) %>%
    left_join(describe_numeric_data, by = c("term" = "variable")) %>%
    relocate(c(names(describe_numeric_data), "n_occurr", "share"), .after = "term")

  return(list("history" = stepwise_history,
              "model" = ifelse(return_model, fit$final_model, "no_model"),
              "predictor" =  fit$final_model$fitted.values,
              "y" = fit$final_model$y))


}
