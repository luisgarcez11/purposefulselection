#' lookup_var_mapping
#'
#' @param data data to be analyzed
#' @param dep_vars dependent variables
#' @param outcome_var outcome variable
#'
#' @return list with lookup information
#' @export
#'
#' @examples
lookup_var_mapping <- function(.data, dep_vars, outcome_var) {
  # store original dep_vars and .data
  original_dep_vars <- if ((!is.null(dep_vars))) {
    dep_vars
  } else {
    names(.data[, (1:ncol(.data) - 1)])
  }
  original_data <- .data

  # get outcome variable as factor
  original_data <- .data %>%
    mutate_at(outcome_var, ~ as.factor(.)) %>%
    mutate_if(is.logical, ~ case_when(
      . == TRUE ~ 1,
      . == FALSE ~ 0
    ))

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

  # removing NAs
  .data <- .data %>%
    stats::na.omit()

  # setting outcome var
  names(.data)[ncol(.data)] <- outcome_var
  outcome_var <- tail(names(.data), 1)

  # cleaning
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

  data_temp <- .data
  var_lookup <- lapply(original_full_model$xlevels, as.data.frame, stringsAsFactors = FALSE)
  if (length(var_lookup) > 0) {
    var_lookup <- var_lookup %>%
      bind_rows(.id = "term") %>%
      dplyr::rename_all(~ c("variable", "term")) %>%
      mutate(term = purrr::map2(
        .x = variable, .y = term,
        .f = ~ paste0(c(.x, .y), collapse = "_")
      )) %>%
      mutate_at("term", ~ tolower(.)) %>%
      filter(term %in% names(data_temp)) %>%
      mutate_at(c("variable", "term"), ~ as.character(.))

    var_lookup <- var_lookup %>% rbind(tibble(
      variable = setdiff(names(data_temp), c(var_lookup$term, outcome_var)),
      term = setdiff(names(data_temp), c(var_lookup$term, outcome_var))
    ))
  } else {
    var_lookup <- tibble::tibble(
      variable = names(data_temp),
      term = names(data_temp)
    )
  }

  terms_by_var <- var_lookup %>%
    group_by(variable) %>%
    mutate(twin_terms1 = list(term)) %>%
    ungroup() %>%
    distinct() %>%
    mutate(twin_terms = purrr::map(.x = twin_terms1, ~ unlist(.x))) %>%
    select(variable, twin_terms) %>%
    distinct()
  lookup_vars_data <- var_lookup %>% left_join(terms_by_var, by = "variable")
  lookup <- var_lookup$variable
  names(lookup) <- var_lookup$term

  return(list(
    "terms_by_var" = terms_by_var,
    "lookup_vars_data" = lookup_vars_data,
    "lookup" = lookup,
    "data" = .data,
    "dep_vars" = names(.data)[-length(.data)]
  ))
}
