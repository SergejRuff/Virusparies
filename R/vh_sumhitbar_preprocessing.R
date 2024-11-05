#' Preprocess virus data for generating a bar plot of the sum of hits for each virus family
#'
#' This function preprocesses virus data for generating a bar plot showing
#' the sum of hits for each virus family from the input dataset.
#'
#' @param vh_file A data frame containing VirusHunter hittables results
#' @param groupby string indicating best_query or ViralRefSeq_taxonomy column
#' @param y_column y_column
#'
#' @return A processed data frame with columns "best_query", "sum", "perc", "res", and "cyl"
#' @author Sergej Ruff
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
#' @noRd
vh_sumhitbar_preprocessing <- function(vh_file,groupby="best_query",y_column ="num_hits"){


  if (y_column == "num_hits" || y_column == "ViralRefSeq_contigs") {

    vh_group <- vh_file %>%
      group_by(.data[[groupby]]) %>%
      summarize(sum=sum(.data[[y_column]]))%>%
      mutate(
        perc = round(proportions(.data$sum) * 100, 2),
        res = str_c(.data$sum, " (", .data$perc, "%)"),
        cyl = as.factor(.data[[groupby]])
      )
  } else {
    row_num <- nrow(vh_file)
    vh_group <- vh_file %>% group_by(.data[[groupby]]) %>% count()%>%
      mutate(
        perc = round(n / row_num *100, 2),
        res = str_c(.data$n, " (", .data$perc, "%)"),
        cyl = as.factor(.data[[groupby]])
      )%>%
      rename(
        sum = n  # Renaming the 'n' column to 'sum'
      )
  }





  return(vh_group)
}
