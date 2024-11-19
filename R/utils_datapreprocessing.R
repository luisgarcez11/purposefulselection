#' Data scan and preprocessing
#'
#' @param directory directory with files
#'
#' @return
#' @export
#'
#' @examples
data_scan_and_preprocess <- function(directory){

  datasets = tibble::tribble(~data ,  ~dep_vars , ~outcome_var ,~max_term_combinations)

  data_folder_files = list.files(directory, full.names = TRUE)

  for (filename in data_folder_files){

    message(paste("reading", filename))

    #readinf
    datafile = tryCatch(expr = { read.csv(filename) %>% janitor::clean_names() },
                        error = function(e){return(read.csv2( filename) %>% janitor::clean_names())},
                        error = function(e){return(xlsx::read.xlsx(filename) %>% janitor::clean_names())},
                        error = function(e){return(xlsx::read.xlsx2(filename) %>% janitor::clean_names())}
    )
    if(ncol(datafile) == 1){
      datafile <- read.csv(filename, sep = ";") %>% janitor::clean_names()
    }

    #selecting vars
    many_unique = sapply(datafile, FUN = function(x) length(unique(x))> 5)
    one_unique = sapply(datafile, FUN = function(x) length(unique(x))==1)
    categorical = sapply(datafile, FUN = function(x) all(is.na(as.numeric(x))))
    to_remove = names(datafile)[which((many_unique == TRUE& categorical == TRUE) | one_unique)]


    #to turn numeric
    many_unique = sapply(datafile, FUN = function(x) length(unique(x))> 15)
    numeric = sapply(datafile, FUN = function(x) sum(!is.na(as.numeric(x)), na.rm = TRUE)> nrow(datafile)*0.10)
    to_numeric = names(datafile)[which(many_unique == TRUE & numeric == TRUE)]
    datafile = datafile %>% mutate_at(to_numeric, ~as.numeric(.))


    skim_object = as.data.frame(skimr::skim(datafile))
    dep_vars = skim_object[,c("skim_variable", "complete_rate")] %>%
      tibble::tibble() %>%
      dplyr::filter(complete_rate > 0.8) %>%
      dplyr::pull(skim_variable) %>%
      setdiff(to_remove)


    #define outcome var
    outcome_var = tail(names(datafile),1)
    dep_vars = setdiff(dep_vars, outcome_var)


    if(length(unique(datafile[[outcome_var]])) > 2){
      if(is.numeric(datafile[[outcome_var]])){
        median_outcome = median(datafile[[outcome_var]], na.rm = TRUE)
        datafile <- datafile %>% mutate_at(outcome_var, ~case_when(. > median_outcome ~ 1,
                                                                   TRUE ~ 0))
      }
    }
    datafile = datafile %>% mutate_at(outcome_var, ~as.factor(.))
    if(length(unique(datafile[[outcome_var]])) !=  2){outcome_var = NA_character_}

    #storing
    datasets <- datasets %>%
      dplyr::add_row("data" = list(datafile), "dep_vars" = list(dep_vars), "outcome_var" = outcome_var,
                     max_term_combinations = 5)

  }

  return(datasets)

}
