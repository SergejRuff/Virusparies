#' @title CombineHittables: Combine hittables.
#'
#' @description
#' CombineHittables combines multiple hittables by row-binding them, provided they have the same column names.
#'
#' @param ... Hittables to be combined.
#'
#' @return A single hittable resulting from the row-binding of all input Hittables.
#'
#' @author Sergej Ruff
#' @seealso
#' VirusHunterGatherer is available here: \url{https://github.com/lauberlab/VirusHunterGatherer}.
#'
#' @examples
#' path <- system.file("extdata", "virushunter.tsv", package = "Virusparies")
#' file <- ImportVirusTable(path)
#' file2 <- ImportVirusTable(path)  # both files have 180 observations
#'
#' combined_file <- CombineHittables(file,file2)
#'
#' print(nrow(combined_file))
#' @export
CombineHittables <- function(...) {
  # Capture the list of hittables
  dfs <- list(...)

  # Check if there are any hittables provided
  if (length(dfs) == 0) {
    stop("No hittables provided.")
  }

  # Ensure all elements in the list are dataframes
  if (!all(sapply(dfs, is.data.frame))) {
    stop("All inputs must be dataframes.")
  }

  # Get the column names of the first dataframe
  column_names <- colnames(dfs[[1]])

  # Check that all dataframes have the same column names
  if (!all(sapply(dfs, function(df) identical(colnames(df), column_names)))) {
    stop("All hittables must have the same column names.")
  }

  # Combine the dataframes by row-binding
  combined_df <- do.call(rbind, dfs)

  return(combined_df)
}
