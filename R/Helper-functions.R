#' Check if hittable is empty
#'
#' @param file a Virushuntegatherer hittable
#'
#' @details
#' Internal helper function which checks if input is empty.objects with 0 observations (hittable,summarystats)
#' returns an error.
#'
#'
#' @keywords internal
is_file_empty <- function(file) {
  if (nrow(file) == 0) {
    stop("Error: Input has zero rows.")
  }
}


#' internal function: message to be printed in boxplot depending on chosen y-column
#'
#' @param y_column y_column
#' @param x_column x_column
#' @param cutoff cutoff
#'
#'
#' @keywords internal
plot_boxplot_message <- function(y_column, x_column, cutoff) {
  messages <- list(
    "ViralRefSeq_E" = paste0("boxplot plotting RefSeq evalues for specified groups (", x_column, ")\nusing the following cut off: ", cutoff),
    "ViralRefSeq_ident" = paste0("boxplot plotting RefSeq Identity for specified groups (", x_column, ")"),
    "contig_len" = paste0("boxplot plotting contig length for specified groups (", x_column, ")")
  )

  if (y_column %in% names(messages)) {
    message(messages[[y_column]])
  } else {
    warning("y_column does not match any expected values.")
  }
}


#' Internal function: adds colour blind support
#'
#' @param plot plot obj
#' @param color color
#'
#' @return plot
#' @importFrom colorBlindness cvdPlot
#'
#' @keywords internal
colorblind_support <- function(plot,color){

  return(cvdPlot(plot=plot,layout=color))

}


