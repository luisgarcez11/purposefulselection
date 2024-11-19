#' Check stepwise interactions
#'
#' @param fit
#' @param pe
#' @param pr
#'
#' @return
#' @export
#'
#' @examples
stepwise_interactions <- function(fit, pe = 0.25, pr = 0.50){

  message("INTERACTIONS TESTED \n")

  testthat::expect_true(pr >= pe, info = "p value for removal must be higher or equal than p value for entrance")


  dep_vars <- setdiff(names(fit$final_model$model), fit$outcome_var)
  possible_interactions <- formula_combination(dep_vars, term_min = 2, term_limit = 2, operator = ":")


  #matrix with p_values
  p_value_matrix_to_include <- matrix(NA,
                           nrow = length(possible_interactions),
                           ncol = 30,
                           dimnames = list(possible_interactions, 1:30))

  #to exclude
  p_value_matrix_to_exclude <- matrix(NA,
                                      nrow = length(possible_interactions),
                                      ncol = 30,
                                      dimnames = list(possible_interactions, 1:30))

  if(fit$direction == "forward"){

  #vectors including the vars included in the model
  included_interactions <- c()
  excluded_interactions <- possible_interactions

  #null model
  form_0 <- fit$final_model$formula
  m_0 <- fit$final_model
  m_0_results <- broom::tidy(m_0)
  final_model <- m_0
  final_form <- form_0

  #p_values for entrance
  p_values_for_entrance <- c()

  #p_values for removal
  p_values_for_removal <- c()

  #models at each step
  stepwise_models <- fit$stepwise_models

  no_interactions_to_exclude <- FALSE
  no_interactions_to_include <- FALSE
  # condition_to_proceed_1 <- length(excluded_interactions) == 0
  condition_to_proceed_2 <- !all(no_interactions_to_exclude & no_interactions_to_include)

  step <- 1

  while( condition_to_proceed_2){

    message(paste("STEP", step, "\n"))

    #STEP increment

    #check for interactions to include
    for( interaction in excluded_interactions){

      m_increment <- glm(formula = update(final_form, paste("~ . +", interaction)),
                         family = "binomial", data = fit$data )

      lr_test <- likelihood_ratio(null_model = final_model, model = m_increment)

      p_value_matrix_to_include[interaction, step] <-  lr_test$p_value
    }

    min_p_value <- min(p_value_matrix_to_include[,step], na.rm = TRUE)
    interaction_to_increment <- names(which(p_value_matrix_to_include[,step] == min_p_value))

    #add variable to the model
    if(min_p_value < pe){

      message(paste(interaction_to_increment, "added to the final model. p-value:", round(min_p_value, 3), "\n"))
      last_model <- final_model
      final_model <- update(final_model, update(final_form, paste("~ . +", paste(interaction_to_increment, collapse = " + "))), data = fit$data)
      final_model_results <- broom::tidy(final_model)
      final_form <- final_model$formula

      stepwise_models <- stepwise_models %>%
        add_row(step = step,
                type = "increment",
                formula = paste(stringr::str_squish(deparse(final_form)), collapse = " "),
                fit = list(final_model),
                results = list(final_model_results))

      included_interactions <- setdiff(attr(terms(final_form), "term.labels"), c(fit$outcome_var, dep_vars))
      excluded_interactions <- setdiff(possible_interactions, included_interactions)

    }else{

      no_interactions_to_include = TRUE
    }

    #STEP removal

    #check for variables to exclude
    for( interaction in included_interactions){

      m_candidate <- glm(formula = update(final_form, paste("~ . -", interaction)),
                         family = "binomial", data = fit$data )

      lr_test <- likelihood_ratio(null_model = m_candidate, model = final_model)

      p_value_matrix_to_exclude[interaction, step] <-  lr_test$p_value

    }

    max_p_value <- max(p_value_matrix_to_exclude[,step], na.rm = TRUE)
    interaction_to_remove <- names(which(p_value_matrix_to_exclude[,step] == max_p_value))

    #remove variable from the model
    if(max_p_value > pr){

      message(paste(interaction_to_remove, "removed from the final model. p-value:", round(max_p_value, 3), "\n"))
      last_model <- final_model
      final_model <- update(final_model, update(final_form, paste("~ . -", paste(interaction_to_remove, collapse = " - "))))
      final_model_results <- broom::tidy(final_model)
      final_form <- final_model$formula

      stepwise_models <- stepwise_models %>%
        add_row(step = step,
                type = "removal",
                formula = deparse(final_form),
                fit = list(final_model),
                results = list(final_model_results))

      included_interactions <- setdiff(setdiff(attr(terms(final_form), "term.labels"), c(fit$outcome_var, dep_vars)), fit$outcome_var)
      excluded_interactions <- setdiff(possible_interactions, included_interactions)

    }else{

      no_interactions_to_exclude = TRUE
    }

    #update conditions to procedd
    # condition_to_proceed_1 <- (length(excluded_interactions) == 0)

    p_values_for_entrance <- p_value_matrix_to_include[,step]
    p_values_for_removal <- p_value_matrix_to_exclude[,step]
    condition_to_proceed_2 <- !(no_interactions_to_exclude & no_interactions_to_include  )

    if(condition_to_proceed_2 == FALSE ){message("STOP!")}

    step <- step + 1



  }}

  if(fit$direction == "backward"){

    #vectors including the vars included in the model
    excluded_interactions  <- c()
    included_interactions <- possible_interactions

    #null model
    form_0 <- fit$final_model$formula
    m_0 <- fit$final_model
    m_0_results <- broom::tidy(m_0)

    final_model <- update(m_0, update(form_0, paste("~ . +", paste(included_interactions, collapse = " + "))), data = fit$data)
    final_form <- final_model$formula



    #p_values for entrance
    p_values_for_entrance <- c()

    #p_values for removal
    p_values_for_removal <- c()

    #models at each step
    stepwise_models <- fit$stepwise_models

    no_interactions_to_exclude <- FALSE
    no_interactions_to_include <- FALSE
    # condition_to_proceed_1 <- length(excluded_interactions) == 0
    condition_to_proceed_2 <- !all(no_interactions_to_exclude & no_interactions_to_include)

    step <- 1

    while( condition_to_proceed_2){

      message(paste("STEP", step, "\n"))


      #STEP removal

      #check for variables to exclude
      for( interaction in included_interactions){

        m_candidate <- glm(formula = update(final_form, paste("~ . -", interaction)),
                           family = "binomial", data = fit$data )

        lr_test <- likelihood_ratio(null_model = m_candidate, model = final_model)

        p_value_matrix_to_exclude[interaction, step] <-  lr_test$p_value

      }

      max_p_value <- max(p_value_matrix_to_exclude[,step], na.rm = TRUE)
      interaction_to_remove <- names(which(p_value_matrix_to_exclude[,step] == max_p_value))

      #remove variable from the model
      if(max_p_value > pr){

        message(paste(interaction_to_remove, "removed from the final model. p-value:", round(max_p_value, 3), "\n"))
        last_model <- final_model
        final_model <- update(final_model, update(final_form, paste("~ . -", paste(interaction_to_remove, collapse = " - "))))
        final_model_results <- broom::tidy(final_model)
        final_form <- final_model$formula

        stepwise_models <- stepwise_models %>%
          add_row(step = step,
                  type = "removal",
                  formula = deparse(final_form),
                  fit = list(final_model),
                  results = list(final_model_results))

        included_interactions <- setdiff(setdiff(attr(terms(final_form), "term.labels"), c(fit$outcome_var, dep_vars)), fit$outcome_var)
        excluded_interactions <- setdiff(possible_interactions, included_interactions)

      }else{

        no_interactions_to_exclude = TRUE
      }


      #STEP increment

      #check for interactions to include
      for( interaction in excluded_interactions){

        m_increment <- glm(formula = update(final_form, paste("~ . +", interaction)),
                           family = "binomial", data = fit$data )

        lr_test <- likelihood_ratio(null_model = final_model, model = m_increment)

        p_value_matrix_to_include[interaction, step] <-  lr_test$p_value
      }

      min_p_value <- min(p_value_matrix_to_include[,step], na.rm = TRUE)
      interaction_to_increment <- names(which(p_value_matrix_to_include[,step] == min_p_value))

      #add variable to the model
      if(min_p_value < pe){

        message(paste(interaction_to_increment, "added to the final model. p-value:", round(min_p_value, 3), "\n"))
        last_model <- final_model
        final_model <- update(final_model, update(final_form, paste("~ . +", paste(interaction_to_increment, collapse = " + "))), data = fit$data)
        final_model_results <- broom::tidy(final_model)
        final_form <- final_model$formula

        stepwise_models <- stepwise_models %>%
          add_row(step = step,
                  type = "increment",
                  formula = paste(stringr::str_squish(deparse(final_form)), collapse = " "),
                  fit = list(final_model),
                  results = list(final_model_results))

        included_interactions <- setdiff(attr(terms(final_form), "term.labels"), c(fit$outcome_var, dep_vars))
        excluded_interactions <- setdiff(possible_interactions, included_interactions)

      }else{

        no_interactions_to_include = TRUE
      }


      #update conditions to procedd
      # condition_to_proceed_1 <- (length(excluded_interactions) == 0)

      p_values_for_entrance <- p_value_matrix_to_include[,step]
      p_values_for_removal <- p_value_matrix_to_exclude[,step]
      condition_to_proceed_2 <- !(no_interactions_to_exclude & no_interactions_to_include  )

      if(condition_to_proceed_2 == FALSE ){message("STOP!")}

      step <- step + 1



    }}



  return(list("stepwise_models" = stepwise_models,
              "outcome_var" = fit$outcome_var,
              "data" = data,
              "possible_interactions" = possible_interactions,
              "final_model" = final_model,
              "final_model_results" = final_model_results,
              "p_value_matrix_to_exclude" = fit$p_value_matrix_to_exclude,
              "p_value_matrix_to_include" = fit$p_value_matrix_to_include,
              "p_value_matrix_to_exclude_interactions" = p_value_matrix_to_exclude,
              "p_value_matrix_to_include_interactions" = p_value_matrix_to_include,
              "direction" = fit$direction,
              "p_value_entrance" = fit$p_value_entrance,
              "p_value_removal" = fit$p_value_removal))






}
