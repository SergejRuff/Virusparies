#' @title boxp_summary_stats:  Calculate summary statistics for boxplot
#'
#' @description This function calculates summary statistics for boxplots created by vhgBoxplot function.
#'
#' @param vh_file A data frame containing VirusHunter hittables.
#' @param group column to group stats by
#' @param ycol column, which specifies identity or evalue.
#'
#' @return A data frame containing summary statistics including median,
#' Q1 (1st quartile), Q3 (3rd quartile), minimum, and maximum values for each virus group
#' @author Sergej Ruff
#' @details This function is an internal utility function used within the package.
#' It calculates summary statistics such as median, quartiles, minimum,
#' and maximum identity values for each virus group from the input dataset.
#' @importFrom dplyr group_by summarise arrange desc
#' @importFrom rlang .data
#' @importFrom stats median quantile median sd
#' @keywords internal
boxp_summary_stats <- function(vh_file,group="best_query",ycol ="ViralRefSeq_ident"){

  if(ycol == "ViralRefSeq_E"){

    vh_file[[ycol]] <- -log10(vh_file[[ycol]])
  }

  ## calculate summary stats for identity

  summary_stats_identity <- vh_file %>%
    group_by(.data[[group]]) %>%
    summarise(
      median = median(.data[[ycol]]),
      Q1 = quantile(.data[[ycol]], 0.25),
      Q3 = quantile(.data[[ycol]], 0.75),
      mean = mean(.data[[ycol]]),       # Calculate mean
      sd = sd(.data[[ycol]]),           # Calculate standard deviation
      min= min(.data[[ycol]]),
      max = max(.data[[ycol]])
    ) %>%
    arrange(desc(median))

  return(summary_stats_identity)

}




