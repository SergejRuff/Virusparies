#' Calculate summary statistics and identify hits below a specified E-value threshold for boxplots
#'
#' This function calculates summary statistics, including median, quartiles, minimum, and maximum
#' -log10-transformed E-values, for each virus group from the VirusHunter hittables.
#' Additionally, it identifies hits below a specified E-value threshold and calculates
#' the percentage of hits below the threshold relative to the total hits for each virus group.
#'
#' @param vh_file A data frame containing VirusHunter hittables results.
#' @param cutoff The significance cutoff value for E-values.
#'
#' @return A data frame containing summary statistics and the percentage of hits below the threshold for each virus group
#'
#' @details This function calculates summary statistics such as median, quartiles, minimum,
#' and maximum -log10-transformed E-values for each virus group from the VirusHunter hittables results.
#' It also identifies hits below a specified E-value threshold and calculates the percentage
#' of hits below the threshold relative to the total hits for each virus group.
#'
#' @importFrom dplyr group_by summarise left_join arrange select cur_data
#' @importFrom stats IQR filter median quantile reorder
#' @importFrom rlang .data
#' @keywords internal
vh_sum_stat_evavlue_boxplot <- function(vh_file,cutoff){


  ## calculate median, 25 and 75 quartile

  summary_stats <- vh_file %>%
    group_by(.data$best_query) %>%
    summarise(
      median = median(-log10(.data$ViralRefSeq_E)),
      Q1 = quantile(-log10(.data$ViralRefSeq_E), 0.25),
      Q3 = quantile(-log10(.data$ViralRefSeq_E), 0.75),
      min = min(-log10(.data$ViralRefSeq_E)),
      max= max(-log10(.data$ViralRefSeq_E))
    ) %>%
    arrange(desc(median))

  below_threshold <- vh_file %>%
    filter(cur_data()[["ViralRefSeq_E"]] < 10^(-cutoff)) %>%
    group_by(.data$best_query) %>%
    summarise(below_threshold = sum(.data$num_hits))

  # Calculate the total hits for each best_query group
  total_hits <- vh_file %>%
    group_by(.data$best_query) %>%
    summarise(total_hits = sum(.data$num_hits))

  below_threshold <- below_threshold %>%
    left_join(total_hits, by = "best_query")%>%
    mutate(percentage_below_thr = (below_threshold / total_hits) * 100)




  # Join the summarized statistics with the below_threshold information
  summary_stats<- summary_stats %>%
    left_join(below_threshold, by = "best_query")%>%replace(is.na(.data$.), 0)%>%
    select (-total_hits)

  return(summary_stats)
}
