#' @title creatTableIdentityBox: Create a formatted table using the gt package for vhEvalIdenBoxplot
#' with eval_vs_iden="identity"
#'
#' @description This function takes the summary statistics output from vhEvalIdenBoxplot
#' and creates a formatted table using the gt package.Should only be used if the output stats were
#' generated with eval_vs_iden="identity".
#'
#' @param summary_stats_identity A data frame containing summary statistics
#'
#' @return A gt table representing the summary statistics
#'
#' @details This function generates a table for the summary statistics output from the vhEvalIdenBoxplot
#' function.It formats summary statistics into a table using the gt package.
#'
#' @examples
#' path <- system.file("extdata", "virushunter.tsv", package = "Virusparies")
#' vh_file <- importVirusTable(path)
#'
#'# plot 1: plot boxplot for "identity"
#' identity <- vhEvalIdenBoxplot(vh_file,eval_vs_iden="identity",cut = 1e-5)
#'
#' # generate table
#' iden_table <- vhEvalBoxTable(identity$summary_stats)
#'
#' iden_table
#'
#' @seealso \code{\link{vhEvalIdenBoxplot}}
#' @import gt
#' @export
vhIdentityBoxTable <- function(summary_stats_identity){
  # Create gt table
  gt_table <- summary_stats_identity %>%
    gt()

  # Print the gt table
  print(gt_table)

  # Return the original summary statistics data frame
  return(gt_table)
}


