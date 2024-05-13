#' @title vh_refeval_iden_boxplot: Create a formatted table using the gt package
#'
#' @description This function takes summary statistics and creates a formatted table using
#' the gt package.
#'
#' @param summary_stats A data frame containing summary statistics
#'
#' @return A formatted gt table
#'
#' @details This function is for internal use within the package.
#' It formats summary statistics into a table
#' using the gt package for further evaluation.
#'
#' @seealso \code{\link{vh_refeval_iden_boxplot}}
#' @importFrom gt gt fmt_number
#' @importFrom rlang .data
#' @keywords internal
creat_table_eval_box <- function(summary_stats){

  # Create gt table
  gt_table <- summary_stats %>%
    gt()%>%
    fmt_number(
      columns = !.data$below_threshold,
      decimals = 2
    )

  gt_table

  return(gt_table)
}
