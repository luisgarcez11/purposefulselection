
#' Forward stepwise
#'
#' @param data data
#' @param outcome_var outcome variable
#' @param dep_vars dependent variables
#' @param pe entry p-value
#' @param pr removal p-value
#'
#' @return
#' @export
#'
#' @examples
stepwise <- function(data, outcome_var,
                     dep_vars = c("age", "weight", "height", "bmi", "priorfrac",
                                  "premeno", "momfrac", "armassist",
                                  "smoke", "raterisk" ),
                     pe = 0.25, pr = 0.50,
                     term_based = TRUE){

  #establish relationship between vars and get data straight
  lookup_list <- lookup_var_mapping(.data = data,
                     dep_vars = dep_vars,
                     outcome_var = outcome_var)
  data <- lookup_list$data
  dep_vars <- lookup_list$dep_vars

  testthat::expect_true(pr >= pe, info = "p value for removal must be higher or equal than p value for entrance")

  #AUC history
  auc_history <- matrix(NA, nrow= 1, ncol = 30,
                        dimnames = list("auc", 1:30))

  #ACC history
  acc_history <- matrix(NA, nrow= 1, ncol = 30,
                        dimnames = list("acc", 1:30))

  #matrix with p_values
  p_value_matrix_to_include <- matrix(NA,
                             nrow = length(dep_vars),
                             ncol = 30,
                             dimnames = list(dep_vars, 1:30))

  #to exclude
  p_value_matrix_to_exclude <- matrix(NA,
                           nrow = length(dep_vars),
                           ncol = 30,
                           dimnames = list(dep_vars, 1:30))

  #vectors including the vars included in the model
  included_vars <- c()
  excluded_vars <- dep_vars

  #null model
  form_0 <- as.formula(paste(outcome_var, "~ 1" ))
  m_0 <- glm(formula = form_0, family = "binomial", data = data )
  m_0_results <- broom::tidy(m_0)
  final_model <- m_0
  final_model_results <- m_0_results
  final_form <- form_0

  #p_values for entrance
  p_values_for_entrance <- c()

  #p_values for removal
  p_values_for_removal <- c()

  #models at each step
  stepwise_models <- tibble::tibble(
    step = 0,
    type = "initial",
    formula = as.character(deparse(form_0)),
    fit = list(m_0),
    results = list(m_0_results)
  )

  no_variables_to_exclude <- FALSE
  no_variables_to_include <- FALSE
  # condition_to_proceed_1 <- length(excluded_vars) > 0
  condition_to_proceed_2 <- !all(no_variables_to_exclude & no_variables_to_include)

  step <- 1

  while(condition_to_proceed_2){

    message(paste("STEP", step, "\n"))

    #STEP increment

    #check for variables to include
    for( var in excluded_vars){

      if(term_based){
        vars = var
        }else{
          vars = lookup_list$lookup_vars_data %>% filter(term == var) %>%
            pull(twin_terms) %>% unlist()
      }

        m_increment <- glm(formula = update(final_form, paste("~ . +", paste(vars, collapse = " + "))),
                           family = "binomial", data = data )
        lr_test <- likelihood_ratio(null_model = final_model, model = m_increment)
        p_value_matrix_to_include[var, step] <-  lr_test$p_value



    }

    min_p_value <- suppressWarnings(min(p_value_matrix_to_include[,step], na.rm = TRUE))
    var_to_increment <- names(which(p_value_matrix_to_include[,step] == min_p_value))

    #add variable to the model
    if(min_p_value < pe){

      message(paste(var_to_increment, "added to the final model. p-value:", round(min_p_value,3), "\n"))
      last_model <- final_model
      final_model <- update(final_model, update(final_form, paste("~ . +", paste(var_to_increment, collapse = " + "))))
      final_model_results <- broom::tidy(final_model)
      final_form <- final_model$formula

      stepwise_models <- stepwise_models %>%
        add_row(step = step,
                type = "increment",
                formula = stringr::str_squish(paste(deparse(final_form), collapse = "")),
                fit = list(final_model),
                results = list(final_model_results))

      included_vars <- setdiff(names(final_model$model), outcome_var)
      excluded_vars <- setdiff(dep_vars, included_vars)

    }else{

      no_variables_to_include = TRUE
    }

    #STEP removal

    #check for variables to exclude
    for( var in included_vars){

      if(term_based){
        vars = var
      }else{
        vars = lookup_list$lookup_vars_data %>% filter(term == var) %>%
          pull(twin_terms) %>% unlist()
      }

      m_candidate <- glm(formula = update(final_form, paste("~ . -", paste(vars, collapse = " - ") )),
                         family = "binomial", data = data )

      lr_test <- likelihood_ratio(null_model = m_candidate, model = final_model)

      p_value_matrix_to_exclude[var, step] <-  lr_test$p_value

    }

    max_p_value <- max(p_value_matrix_to_exclude[,step], na.rm = TRUE)
    var_to_remove <- names(which(p_value_matrix_to_exclude[,step] == max_p_value))

    #remove variable from the model
    if(max_p_value > pr){

      message(paste(var_to_remove, "removed from the final model. p-value:", round(max_p_value, 3), "\n"))
      last_model <- final_model
      final_model <- update(final_model, update(final_form, paste("~ . -", paste(var_to_remove, collapse = " - "))))
      final_model_results <- broom::tidy(final_model)
      final_form <- final_model$formula

      stepwise_models <- stepwise_models %>%
        add_row(step = step,
                type = "removal",
                formula = deparse(final_form),
                fit = list(final_model),
                results = list(final_model_results))

      included_vars <- setdiff(names(final_model$model), outcome_var)
      excluded_vars <- setdiff(dep_vars, included_vars)

    }else{

      no_variables_to_exclude = TRUE
    }

    #update conditions to procedd
    # condition_to_proceed_1 <- (length(excluded_vars) > 0)

    p_values_for_entrance <- p_value_matrix_to_include[,step]
    p_values_for_removal <- p_value_matrix_to_exclude[,step]
    condition_to_proceed_2 <- !(no_variables_to_exclude & no_variables_to_include  )

    #save AUC
    auc_history["auc", step] <- pROC::auc(final_model$y, final_model$fitted.values)

    #acc
    acc_history["acc", step] <-  accuracy(final_model)


    if(condition_to_proceed_2 == FALSE ){message("STOP! \n")}




    step <- step + 1



  }

  #add variable names
  final_model_results <- final_model_results %>%
    left_join(lookup_list$lookup_vars_data[,c("term", "variable")], by = "term") %>%
    relocate("variable", 1)

  return(list("stepwise_models" = stepwise_models,
              "outcome_var" = outcome_var,
              "data" = data,
              "dep_vars" = dep_vars,
              "lookup_list" = lookup_list,
              "final_model" = final_model,
              "final_model_results" = final_model_results,
              "p_value_matrix_to_exclude" = p_value_matrix_to_exclude,
              "p_value_matrix_to_include" = p_value_matrix_to_include,
              "direction" = "forward",
              "p_value_entrance" = pe,
              "p_value_removal" = pr,
              "auc_history" = auc_history,
              "acc_history" = acc_history))



}
