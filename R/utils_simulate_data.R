
#' create categorical variable 2
#'
#' @param n_cats number of categories
#' @param n_obs number of observations
#' @param cat_comp_slope lambda from the exponential variable to define categorical values imbalance
#' @param cramer_goal cramer_goal
#'
#' @return
#' @export
#'
#' @examples
simulate_cat_var2 <- function(n_obs,
                              n_cats,
                              y_proportion,
                              cramer_goal){


  #categorical variable compostions weights (uniform)
  cat_comp_slope = 0
  model <- function(cat_comp_slope, x){-cat_comp_slope*(x-0.5)+0.5}
  scale_values <- function(x){(x-min(x))/(max(x)-min(x))}
  partitions <- scale_values(1:n_cats)
  weights_cat_comp <- model(cat_comp_slope, partitions)
  cat_distribution <- round(n_obs * weights_cat_comp/sum(weights_cat_comp), 0)

  #categorical values allocation to outcome variable (1 or 0)
  outcome_alloc_slope = 0.5 #initial
  model <- function(outcome_alloc_slope, x){
    #
    y_pred <- -outcome_alloc_slope*(x-0.5)+0.5
    if(any(y_pred>1)){y_pred[which(y_pred > 1)] = 1}
    if(any(y_pred<0)){y_pred[which(y_pred < 0)] = 0}
    return(y_pred)
  }
  scale_values <- function(x){(x-min(x))/(max(x)-min(x))}
  partitions <- scale_values(partitions)
  p_X_given_Y1 <- model(outcome_alloc_slope, partitions)
  p_X_given_Y0 <- 1 -  p_X_given_Y1

  # Step 1: Generate the binary variable Y
  outcome_values <- factor( rep(c(1,0),
                                times= c(y_proportion*n_obs, (1-y_proportion)*n_obs)) )
  data <- tibble::tibble(y = outcome_values, cat_var = NA)

  # Step 2: Iteratively adjust probabilities
  tolerance <- 0.01  # Tolerance for the Cramér's V difference
  max_iter <- 1000  # Maximum iterations
  iter <- 0
  cramer_stats <- -1
  history <- tibble(cramer_stats_ = numeric(),
                    p_X_given_Y1_ = double(),
                    p_X_given_Y0_ = double()
  )

  i = 1
  while(!between(cramer_stats, cramer_goal - tolerance, right = cramer_goal + tolerance) ){

    #categorical values
    cat_values <- factor( rep( toupper(letters)[1:n_cats],
                               times= cat_distribution))

    #set categorical variable
    data[which(data$y == 1), "cat_var"] <- c(sample(x = levels(cat_values),
                                                    prob = p_X_given_Y1,
                                                    size = sum(data$y == 1)-length(levels(cat_values)),
                                                    replace = TRUE),
                                             levels(cat_values))
    data[which(data$y == 0), "cat_var"] <- c(sample(x = levels(cat_values),
                                                    prob = p_X_given_Y0,
                                                    size = sum(data$y == 0)-length(levels(cat_values)),
                                                    replace = TRUE),
                                             levels(cat_values))


    #contingency
    cont_table <- table(data$y, data$cat_var)
    # print(cont_table)

    #cramer V
    cramer_stats <- vcd::assocstats(cont_table)$cramer
    # print(cramer_stats)
    # print(p_X_given_Y1)
    # print(p_X_given_Y0)

    #adjust slope
    if(cramer_stats > cramer_goal){
      partitions <- scale_values(partitions)
      outcome_alloc_slope = outcome_alloc_slope-0.01
      p_X_given_Y1 <- model(outcome_alloc_slope, partitions)
      p_X_given_Y0 <- 1 -  p_X_given_Y1
    }
    if(cramer_stats < cramer_goal){
      partitions <- scale_values(partitions)
      outcome_alloc_slope = outcome_alloc_slope+0.01
      p_X_given_Y1 <- model(outcome_alloc_slope, partitions)
      p_X_given_Y0 <- 1 -  p_X_given_Y1
    }

    #update history
    history <- history %>% tibble::add_row(cramer_stats_ = cramer_stats,
                                           p_X_given_Y1_ = p_X_given_Y1[1],
                                           p_X_given_Y0_ = p_X_given_Y0[1])

    i = i + 1
    if(i > 10000){stop("maximum number of iterations reached")}

  }

  return(list("cat_values" = cat_values,
              "n_obs" = n_obs,
              "n_cats" = n_cats,
              "cont_table" = cont_table,
              "data" = data,
              "p_X_given_Y1" = p_X_given_Y1,
              "p_X_given_Y0" = p_X_given_Y0,
              "history" = history))

}


#' Generate data
#'
#' @param y_proportion Proportion of the outcome
#' @param n_obs  Number of observations
#' @param x_number_of_variables Total number of variables
#' @param x_cat_vars_prop Proportion of categorical variables
#' @param x_number_of_cats_per_var Vector with number of categories for each categorical variable
#' @param x_comp_imbalance Exponential distribution Lambda for attributing the number of occurrences per category
#' @param xy_imbalance Slope to estimate the allocation of categorical variables to outcome
#'
#' @return
#' @export
#'
#' @examples
simulate_data <- function(y_proportion = 0.1,
                          n_obs = 1000,
                          x_number_of_variables = 10,
                          x_number_of_cats_per_var = rep(3, x_number_of_variables),
                          cramer_goal = rep(0.6, x_number_of_variables)){


  #y outcome
  inverse_y_proportion = 1-y_proportion
  n_controls <- as.integer((inverse_y_proportion * as.integer(n_obs)))
  n_cases <- y_proportion * n_obs
  outcome_values <- c(rep(1, times = n_cases),
                      rep(0, times =  n_controls ))
  data <- tibble::tibble(y =outcome_values)

  #add cat variables
  for(i in 1:x_number_of_variables){

    #cat composition
    categorical_values <- simulate_cat_var2(n_obs = n_obs,
                                            n_cats = x_number_of_cats_per_var[i],
                                            y_proportion = y_proportion,
                                            cramer_goal = cramer_goal[i])

    #set categorical variable
    data[which(data$y == 1), paste0("cat_var_", i)] <- categorical_values$data$cat_var[which(categorical_values$data$y == 1)]
    data[which(data$y == 0), paste0("cat_var_", i)] <- categorical_values$data$cat_var[which(categorical_values$data$y == 0)]

    # message(paste0("creating categorical variable ", i))
  }

  return(data)


}


#' Simualte data based on frequency
#'
#' @param n_obs
#' @param y_prop
#' @param n_cats
#' @param sample_sizes
#' @param prob_Y_given_X
#'
#' @return
#' @export
#'
#' @examples
simulate_data_frequency_based <- function(n_cats, sample_sizes,
                                          prob_Y_given_X = list(c(0.5,0.5), c(0.3,0.7), c(0.3,0.7),
                                                                c(0.3,0.7), c(0.1,0.9)),
                                          add_random = FALSE){

  testthat::expect_equal(length(prob_Y_given_X), n_cats)

  #set names
  names_vars <- LETTERS[1:n_cats]

  #set sample variable
  sample_var <- rep(names_vars, sample_sizes)
  simulated_data <- tibble::tibble(var = sample_var) %>%
    mutate(y =  NA_integer_)

  #simulate data
  for (i in 1:n_cats){

    var = names_vars[i]
    simulated_data$y[which(simulated_data$var == var)] <- sample(c(0,1),
                                                                 size = sample_sizes[i],
                                                                 prob =  prob_Y_given_X[i][[1]],
                                                                 replace = TRUE)
  }

  if(add_random){
    simulated_data <- simulated_data %>%
      mutate(random = rnorm(nrow(.)))
  }

  return(list(data = simulated_data))

}

