#' @title boxp_summary_stats:  Calculate summary statistics for boxplot
#'
#' @description This function calculates summary statistics for boxplots created by vhgBoxplot function.
#'
#' @param vh_file A data frame containing VirusHunter hittables.
#' @param group column to group stats by
#' @param ycol column, which specifies sequence identity or evalue.
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
#' @noRd
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

  # Check for NA values in the sd column and provide a specific message
  if (any(is.na(summary_stats_identity$sd))) {
    warning("Standard deviation values are NA for some groups. This typically occurs when there is only one observation in the group.")
  }



  if (ycol == "ViralRefSeq_ident") {
    summary_stats_identity <- summary_stats_identity %>%
      mutate(
        median = round(.data$median, 1),
        Q1 = round(.data$Q1, 1),
        Q3 = round(.data$Q3, 1),
        mean = round(.data$mean, 1),
        sd = round(.data$sd, 1),
        min = round(.data$min, 1),
        max = round(.data$max, 1)
      )
  } else if (ycol == "contig_len") {
    summary_stats_identity <- summary_stats_identity %>%
      mutate(
        median = round(.data$median),
        Q1 = round(.data$Q1),
        Q3 = round(.data$Q3),
        mean = round(.data$mean),
        sd = round(.data$sd),
        min = round(.data$min),
        max = round(.data$max)
      )
  } else if (ycol == "ViralRefSeq_E") {
    summary_stats_identity <- summary_stats_identity %>%
      mutate(
        median = round(.data$median, 2),
        Q1 = round(.data$Q1, 2),
        Q3 = round(.data$Q3, 2),
        mean = round(.data$mean, 2),
        sd = round(.data$sd, 2),
        min = round(.data$min, 2),
        max = round(.data$max, 2)
      )
  }

  return(summary_stats_identity)

}




