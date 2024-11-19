

#' write history
#'
#' @param history
#' @param write_path
#' @param sheet_name_write
#'
#' @return
#' @export
#'
#' @examples
write_history <- function(history_fit, write_path , sheet_name_write = ""){

  message("writing history...")

  vars_desc <- names(history_fit[["history"]])[which(stringr::str_detect(names(history_fit[["history"]]), pattern = "_desc"))]

  cond_table <- condformat::condformat(history_fit[["history"]])

   cond_table <- cond_table   %>%
     condformat::rule_fill_discrete(columns = "variable",expression = TRUE,
                        colours =  c("TRUE" = "yellow")) %>%
     condformat::rule_fill_discrete(columns = vars_desc[1],
                       colours =c("in" = "lightgreen", "out" = "red"), lockcells = TRUE) %>%
     condformat::rule_fill_discrete(columns = vars_desc[2],
                        colours =c("in" = "lightgreen", "out" = "red"), lockcells = TRUE) %>%
     condformat::rule_fill_discrete(columns = vars_desc[3],
                        colours =c("in" = "lightgreen", "out" = "red"), lockcells = TRUE) %>%
     condformat::rule_fill_discrete(columns = vars_desc[4],
                        colours =c("in" = "lightgreen", "out" = "red"), lockcells = TRUE) %>%
     condformat::rule_fill_discrete(columns = vars_desc[5],
                        colours =c("in" = "lightgreen", "out" = "red"), lockcells = TRUE) %>%
     condformat::rule_fill_discrete(columns = vars_desc[6],
                        colours =c("in" = "lightgreen", "out" = "red"), lockcells = TRUE) %>%
     condformat::rule_fill_discrete(columns = vars_desc[7],
                        colours =c("in" = "lightgreen", "out" = "red"), lockcells = TRUE) %>%
     condformat::rule_fill_discrete(columns = vars_desc[8],
                       colours =c("in" = "lightgreen", "out" = "red"), lockcells = TRUE) %>%
     condformat::rule_fill_discrete(columns = vars_desc[9],
                        colours =c("in" = "lightgreen", "out" = "red"), lockcells = TRUE) %>%
     condformat::rule_fill_discrete(columns = vars_desc[10],
                        colours =c("in" = "lightgreen", "out" = "red"), lockcells = TRUE) %>%
     condformat::rule_fill_discrete(columns = vars_desc[11],
                        colours =c("in" = "lightgreen", "out" = "red"), lockcells = TRUE) %>%
     condformat::rule_fill_discrete(columns = vars_desc[12],
                        colours =c("in" = "lightgreen", "out" = "red"), lockcells = TRUE)%>%
     condformat::rule_fill_discrete(columns = vars_desc[13],
                                    colours =c("in" = "lightgreen", "out" = "red"), lockcells = TRUE)%>%
     condformat::rule_fill_discrete(columns = vars_desc[14],
                                    colours =c("in" = "lightgreen", "out" = "red"), lockcells = TRUE)%>%
     condformat::rule_fill_discrete(columns = vars_desc[15],
                                    colours =c("in" = "lightgreen", "out" = "red"), lockcells = TRUE)%>%
     condformat::rule_fill_discrete(columns = vars_desc[16],
                                    colours =c("in" = "lightgreen", "out" = "red"), lockcells = TRUE)%>%
     condformat::rule_fill_discrete(columns = vars_desc[17],
                                    colours =c("in" = "lightgreen", "out" = "red"), lockcells = TRUE) %>%
     condformat::rule_fill_discrete(columns = vars_desc[18],
                                    colours =c("in" = "lightgreen", "out" = "red"), lockcells = TRUE) %>%
     condformat::rule_fill_discrete(columns = vars_desc[19],
                                    colours =c("in" = "lightgreen", "out" = "red"), lockcells = TRUE)


   attr(cond_table, "condformat")[["rules"]][(length(vars_desc)+2):length(attr(cond_table, "condformat")[["rules"]])] <- NULL

   condformat::condformat2excel(cond_table,
                                filename = write_path,
                                sheet_name = sheet_name_write)
}


