#' Check if hittable is empty
#'
#' @param vh_file a Virushuntegatherer hittable
#'
#' @details
#' Internal helper function which checks if hittables are empty. Returns an error message if they are.
#'
#'
#' @keywords internal
is_vh_file_empty <- function(vh_file) {
  if (nrow(vh_file) == 0) {
    stop("Error: Hittables has zero rows.")
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
