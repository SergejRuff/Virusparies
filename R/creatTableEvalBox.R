#' @title vh_refeval_iden_boxplot: Create a formatted table using the gt package for vhEvalIdenBoxplot
#' with eval_vs_iden="evalue"
#'
#' @description This function takes the summary statistics output from vhEvalIdenBoxplot
#' and creates a formatted table using the gt package.Should only be used if the output stats were
#' generated with eval_vs_iden="evalue".
#'
#' @param summary_stats A data frame containing summary statistics
#'
#' @return A formatted gt table
#'
#' @details This function generates a table for the summary statistics output from the vhEvalIdenBoxplot
#' function.It formats summary statistics into a table using the gt package.
#'
#' @seealso \code{\link{vhEvalIdenBoxplot}}
#' @importFrom gt gt fmt_number
#' @importFrom rlang .data
#' @export
creatTableEvalBox <- function(summary_stats){

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
