#' Univariate fit (STEP 1)
#'
#' @param .data data to be analyzed
#' @param dep_vars Dependent variables
#' @param outcome_var Outcome variable
#' @param term_based Term based or not
#' @param multivariate Multivariate insted of Univariate as STEP1
#'
#' @return
#' @export
#'
#' @examples
univariate_log_fit <- function(.data, dep_vars = NULL, outcome_var = NULL,
                               term_based = FALSE, multivariate = FALSE) {

  #store original dep_vars and .data
  original_dep_vars <- if((!is.null(dep_vars))){dep_vars}else{names(.data[,(1:ncol(.data)-1)])}
  original_data <- .data

  #get outcome variable as factor
  original_data <- .data %>%
    mutate_at(outcome_var, ~ as.factor(.)) %>%
    mutate_if(is.logical, ~case_when(. == TRUE ~ 1,
                                     . == FALSE ~ 0))

  # if stated, define and order dep and outcome var
  if (!is.null(dep_vars) & !is.null(outcome_var)) {
    .data <- .data %>%
      select(all_of(c(dep_vars, outcome_var)))
  }

  # dummyfy the dataset
  if (any(sapply(.data, FUN = function(x) (is.character(x) | is.factor(x))))) {
    .data <- .data %>% fastDummies::dummy_columns(
      remove_selected_columns = TRUE,
      remove_most_frequent_dummy = TRUE
    )
    .data <- .data %>% dplyr::relocate(outcome_var, .after = ncol(.))
  }

  #removing NAs
  .data <- .data %>%
    stats::na.omit()

  #setting outcome var
  names(.data)[ncol(.data)] <- outcome_var
  outcome_var <- tail(names(.data), 1)

  #cleaning
  .data <- .data %>%
    mutate_at(outcome_var, ~ as.factor(.)) %>%
    janitor::clean_names() %>%
    rename_all(~ stringr::str_replace_all(., pattern = "-", replacement = "_"))

  # defining outcome and dep_vars
  dep_vars <- names(.data)[1:(ncol(.data) - 1)]
  outcome_var <- tail(names(.data), 1)

  # lookup
  suppressWarnings({
    original_full_model <- glm(
      formula = paste(outcome_var, "~", paste(original_dep_vars, collapse = " + ")),
                      family = "binomial", data = original_data
      )
  })

  data_temp = .data
  var_lookup <- lapply(original_full_model$xlevels, as.data.frame, stringsAsFactors = FALSE)
  if(length(var_lookup) > 0){

    var_lookup <- var_lookup%>%
      bind_rows(.id = "term")%>%
      dplyr::rename_all(~ c("variable", "term")) %>%
      mutate(term = purrr::map2(.x = variable, .y = term,
                              .f = ~ paste0(c(.x, .y), collapse = "_"))) %>%
      mutate_at("term", ~tolower(.)) %>%
      filter( term %in% names(data_temp)) %>%
      mutate_at(c("variable", "term"), ~as.character(.))

    var_lookup <- var_lookup %>% rbind(tibble(variable = setdiff(names(data_temp), c(var_lookup$term, outcome_var)),
                                            term = setdiff(names(data_temp), c(var_lookup$term, outcome_var))))
  }else{

    var_lookup <- tibble::tibble(variable = names(data_temp),
                         term = names(data_temp))

  }

  terms_by_var <- var_lookup %>%
    group_by(variable) %>%
    mutate(twin_terms1 = list(term)) %>%
    ungroup() %>%
    distinct() %>%
    mutate(twin_terms = purrr::map(.x = twin_terms1, ~unlist(.x))) %>%
    select(variable, twin_terms) %>%
    distinct()
  lookup_vars_data <- var_lookup %>% left_join(terms_by_var, by = "variable")
  lookup <- var_lookup$variable
  names(lookup) <- var_lookup$term

  #full model
  suppressWarnings({
    full_model <- glm(
      formula = paste(outcome_var, "~", paste(dep_vars, collapse = " + ")),
      family = "binomial", data = .data
  )
  })

  # defining variable types
  num_vars <- setdiff(variable_types(.data)[["num_vars"]], outcome_var)
  cat_vars <- setdiff(variable_types(.data)[["cat_vars"]], outcome_var)
  term_type <- tibble::tibble(term = dep_vars) %>%
    mutate(type = case_when(
      term %in% num_vars ~ "numeric",
      term %in% cat_vars ~ "categorical"
    ))

  # stop if outcome is not viable
  if (length(outcome_var) != 1) {
    stop("outcome_var must be length == 1")
  }

  # store results tibble
  results_tib <- tibble::tibble(
    term = character(),
    estimate = numeric(),
    std.error = numeric(),
    statistic = numeric(),
    p.value_wald = numeric(),
    statistic_lr = numeric(),
    p.value_lr = numeric()
  )

  # null model
  null_model <- glm(
    formula = as.formula(paste(outcome_var, "~ 1")),
    family = "binomial", data = .data
  )
  df.null_model <- 0

  #multivariate or univariate
  if (!multivariate) {

    # get univariate models
    for (var in dep_vars) {

      #if term already in results, skip it
      if(var %in% results_tib$term){next}

      #build the models for selection
      if(term_based){

        # fit model
        form <- as.formula(paste(outcome_var, "~", var))
        m <- glm(formula = form, family = "binomial", data = .data)


      }else{

        if(var %in% var_lookup$term){
          #terms by variable
          var_of_interest = var_lookup %>%
            filter(term == var) %>%
            pull(variable)
          terms_twins = var_lookup %>%
            filter(variable == var_of_interest ) %>%
            pull(term) %>% unlist()
          }else{

            terms_twins = var

            }

        # fit model
        form <- as.formula(paste(outcome_var, "~", paste(terms_twins, collapse = " + ")))
        m <- glm(formula = form, family = "binomial", data = .data)

      }

        # Likelihood ratio test
        lik_ratio <- likelihood_ratio(null_model, m)

        # broom coeficients
        m_coef <- broom::tidy(m) %>%
          dplyr::rename("p.value_wald" = p.value) %>%
          mutate(variable = lookup[term]) %>%
          mutate(
            statistic_lr = lik_ratio$lr_stat,
            p.value_lr = lik_ratio$p_value
          )


      # results
      results_tib <- results_tib %>%
        dplyr::bind_rows(m_coef[2:nrow(m_coef), ]) %>%
        mutate(
          or_estimate = exp(estimate),
          or_estimate_lower = exp(estimate - std.error),
          or_estimate_upper = exp(estimate + std.error),
          .after = estimate
        ) %>%
        dplyr::relocate("variable", .after = "term")
    }
    }else {

    data1 = .data

    #results
    results_tib <- broom::tidy(full_model) %>%
      filter(term != "(Intercept)" )%>%
      mutate(
        or_estimate = exp(estimate),
        or_estimate_lower = exp(estimate - std.error),
        or_estimate_upper = exp(estimate + std.error),
        .after = estimate
      ) %>%
      dplyr::rename("p.value_wald" = p.value) %>%
      mutate(variable = lookup[term]) %>%
      dplyr::relocate("variable", .after = "term")

    #terms by variable
    terms_by_var = results_tib %>% group_by(variable) %>%
      summarise(terms = list(term)) %>% ungroup()

    #add likelihood raatio
    results_tib <- results_tib%>%
      left_join(terms_by_var) %>%
      mutate(p_value_lr = purrr::map(.x = terms,
                                     .f = ~likelihood_ratio(glm(formula = update_formula(full_model$formula, " - ", paste(.x[[1]], collapse = " - ")),
                                               family = "binomial", data = data1 ), full_model)$p_value)) %>%
      select(-terms)




  }

  # add varuable types
  data2 = .data
  results_tib <- results_tib %>%
    left_join(term_type, by = "term") %>%
    dplyr::relocate("type", .after = "variable")



  # report message with missing values
  # account for categorical vars

  return(list(
    "data" = .data,
    "outcome_var" = outcome_var,
    "dep_vars" = dep_vars,
    "lookup" = lookup,
    "lookup_vars_data" = lookup_vars_data,
    "univariate_results" = results_tib,
    "term_based" = term_based,
    "num_terms" = num_vars,
    "cat_terms" = cat_vars
  ))
}
