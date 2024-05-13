#' @title Calculate summary statistics for identity
#'
#' @description This function calculates summary statistics for identity based on the input virus data.
#'
#' @param vh_file A data frame containing VirusHunter hittables.
#'
#' @return A data frame containing summary statistics including median,
#' Q1 (1st quartile), Q3 (3rd quartile), minimum, and maximum identity values for each virus group
#'
#' @details This function is an internal utility function used within the package.
#' It calculates summary statistics such as median, quartiles, minimum,
#' and maximum identity values for each virus group from the input dataset.
#' @importFrom dplyr group_by summarise arrange desc
#' @importFrom rlang .data
#' @keywords internal
summary_stats_identity <- function(vh_file){

  ## calculate summary stats for identity

  summary_stats_identity <- vh_file %>%
    group_by(.data$best_query) %>%
    summarise(
      median = median(.data$ViralRefSeq_ident),
      Q1 = quantile(.data$ViralRefSeq_ident, 0.25),
      Q3 = quantile(.data$ViralRefSeq_ident, 0.75),
      min= min(.data$ViralRefSeq_ident),
      max = max(.data$ViralRefSeq_ident)
    ) %>%
    arrange(desc(median))

  return(summary_stats_identity)

}




