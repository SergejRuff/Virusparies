#' Preprocess virus data for generating a bar plot of the sum of hits for each virus family
#'
#' This function preprocesses virus data for generating a bar plot showing
#' the sum of hits for each virus family from the input dataset.
#'
#' @param vh_file A data frame containing VirusHunter hittables results
#'
#' @return A processed data frame with columns "best_query", "sum", "perc", "res", and "cyl"
#'
#' @details This function preprocesses virus data for generating a bar plot showing the
#' sum of hits for each virus family from the input dataset. It calculates the sum of hits
#' for each virus family, computes the percentage of the sum relative to the total hits,
#' creates a string representation for visualization, and converts the "best_query" column
#' to a factor.
#'
#'
#' @importFrom dplyr group_by summarize mutate
#' @importFrom stringr str_c
#' @importFrom rlang .data
#' @keywords internal
vh_sumhitbar_preprocessing <- function(vh_file){

  vh_group <- vh_file %>%
    group_by(.data$best_query) %>%
    summarize(sum=sum(.data$num_hits))%>%
    mutate(
      perc = round(proportions(.data$sum) * 100, 2),
      res = str_c(.data$sum, " (", .data$perc, "%)"),
      cyl = as.factor(.data$best_query)
    )

  return(vh_group)
}
