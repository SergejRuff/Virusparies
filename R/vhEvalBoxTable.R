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
#' @examples
#' path <- system.file("extdata", "virushunter.tsv", package = "Virusparies")
#' vh_file <- importVirusTable(path)
#'
#' # plot 1: plot boxplot for "evalue"
#' eval <- vhEvalIdenBoxplot(vh_file,eval_vs_iden="evalue",cut = 1e-5)
#'
#' # generate table
#' eval_table <- vhEvalBoxTable(eval$summary_stats)
#'
#' @seealso \code{\link{vhEvalIdenBoxplot}}
#' @import dplyr
#' @importFrom gt gt fmt_number
#' @importFrom rlang .data
#' @export
vhEvalBoxTable <- function(summary_stats){

  # Create gt table
  gt_table <- summary_stats %>%
    gt()%>%
    fmt_number(

      decimals = 2
    )



  return(gt_table)
}
