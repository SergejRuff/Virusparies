#' @title Convert summary statistics to a gt table
#'
#' @description This function converts summary statistics to a gt table for visualization.
#'
#' @param summary_stats_identity A data frame containing summary statistics
#'
#' @return A gt table representing the input summary statistics
#'
#' @details This function is an internal utility function used within the package.
#' It converts summary statistics stored in a data frame to a formatted gt table.
#'
#' @import gt
#' @keywords internal
summary_stats_identity_gt <- function(summary_stats_identity){

  summary_stats_identity %>% gt()

  summary_stats_identity

  return(summary_stats_identity)

}



