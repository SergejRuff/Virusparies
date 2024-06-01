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





#' internal: check if all columns exist
#'
#' @param file file
#' @param columns single column or list of columns
#'
#'
#' @keywords internal
check_columns <- function(file,columns){

  # Get all column names from file
  all_names <- names(file)

  # Define the required columns
  required_columns <- columns

  # Check if each required column exists in file
  for (col in required_columns) {
    if (!(col %in% all_names)) {
      stop("Column '", col, "' not found in file. Available column names: ", paste(all_names, collapse = ", "))
    }
  }

}


#' internal: check if datatype is correct
#'
#' @param vh_file file
#' @param columns object single or list specifying column
#' @param option 1 for chacter,2 for numeric
#'
#'
#' @keywords internal
check_input_type <- function(vh_file, columns, option) {
  # Check if the option is valid (either 1 or 2)
  if (!option %in% c(1, 2)) {
    stop("Invalid option. Please choose 1 for character or 2 for numeric.")
  }

  # Define the expected class based on the option
  expected_class <- if (option == 1) "character" else c("numeric", "integer")

  # Check if the specified columns exist in vh_file
  all_names <- names(vh_file)
  for (col in columns) {
    if (!(col %in% all_names)) {
      stop("Column '", col, "' not found in vh_file. Available column names: ", paste(all_names, collapse = ", "))
    }

    # Check the class of the column
    if (!inherits(vh_file[[col]], expected_class)) {
      stop("Error: Column '", col, "' must be of type ", expected_class, ".")
    }
  }
}

#' internal: add colorblind support by using viridis color palettes
#'
#' @param plot plot obj.
#' @param option A character string indicating the color map option to use. Eight options are available:
#'   - "magma" (or "A")
#'   - "inferno" (or "B")
#'   - "plasma" (or "C")
#'   - "viridis" (or "D")
#'   - "cividis" (or "E")
#'   - "rocket" (or "F")
#'   - "mako" (or "G")
#'   - "turbo" (or "H")
#'
#' @return plot obj
#'
#' @importFrom viridis scale_fill_viridis
#' @keywords internal
colorbildsupport <- function(plot,option){
  return(plot+
           scale_fill_viridis(option=option,discrete = TRUE) )
}

#' internal: add colorblind support by using viridis color palettes for identity scatterplot
#'
#' @param plot plot obj.
#' @param option A character string indicating the color map option to use. Eight options are available:
#'   - "magma" (or "A")
#'   - "inferno" (or "B")
#'   - "plasma" (or "C")
#'   - "viridis" (or "D")
#'   - "cividis" (or "E")
#'   - "rocket" (or "F")
#'   - "mako" (or "G")
#'   - "turbo" (or "H")
#'
#' @return plot obj
#'
#' @importFrom viridis scale_color_viridis
#' @keywords internal
colorbildsupport_ <- function(plot,option){
  return(plot+scale_color_viridis(option=option,discrete = TRUE) )
}


