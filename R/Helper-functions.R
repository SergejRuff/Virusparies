#' Check if hittable or other files are empty
#'
#' @param file a Virushuntegatherer hittable or another file
#'
#' @details
#' Internal helper function which checks if input is empty.objects with 0 observations (hittable,summarystats)
#' returns an error.
#'
#' @author Sergej Ruff
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
#' @author Sergej Ruff
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
#' @author Sergej Ruff
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
#' @author Sergej Ruff
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
#' @author Sergej Ruff
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
#' @author Sergej Ruff
#' @importFrom viridis scale_color_viridis
#' @keywords internal
colorbildsupport_ <- function(plot,option){
  return(plot+scale_color_viridis(option=option,discrete = TRUE) )
}


#' checks if correct datatype is provided to function (character check)
#'
#' @param name argument from function.
#'
#' @author Sergej Ruff
#' @keywords internal
arg_character <- function(name) {
  if (!is.character(name)) {
    errorMessage <- paste("Input argument", deparse(substitute(name)), "must be a character")
    stop(errorMessage)
  }
}

#' checks if correct datatype is provided to function (numeric check)
#'
#' @param name argument from function.
#'
#' @details
#' currently not in use. Might be in a later version.
#'
#' @author Sergej Ruff
#'
#' @keywords internal
arg_numeric_or_integer <- function(name) {
  if (!(is.numeric(name) || is.integer(name))) {
    errorMessage <- paste("Input argument", deparse(substitute(name)), "must be numeric or integer")
    stop(errorMessage)
  }
}

#' checks if correct datatype is provided to function (logical check)
#'
#' @param name argument from function.
#'
#' @author Sergej Ruff
#' @keywords internal
arg_logical <- function(name) {
  if (!is.logical(name)) {
    errorMessage <- paste("Input argument", deparse(substitute(name)), "must be logical")
    stop(errorMessage)
  }
}



#' check if obj is a dataframe
#'
#' @param obj obj
#'
#' @return An error if obj is not a dataframe
#' @author Sergej Ruff
#' @keywords internal
check_is_dataframe <- function(obj) {
  if (!is.data.frame(obj)) {
    stop("The provided object is not a data.frame.")
  }

}



#' check if dataframe has column names.
#'
#' @param df a dataframe obj
#'
#' @return error message, if df has no column names
#' @author Sergej Ruff
#' @keywords internal
has_columnnames <- function(df) {


  if (is.null(colnames(df))) {
    stop("The data frame has no column names.")
  }

  if (any(colnames(df) == "")) {
    stop("The data frame has empty column names.")
  }


}


#' internal function that extracts number of rows below threshold for vhgBoxplot function
#'
#' @param vh_file vh_file
#' @param cut set cutoff in vhGBoxplot function
#'
#' @return filtered dataframe
#'
#' @keywords internal
vhg_filter_belowthresholdboxplot <- function(vh_file,cut){

  return(vh_file[vh_file$ViralRefSeq_E<cut,])
}



#' Internal function to set default title and cutoff in boxplot function
#'
#' @param y_column y-column
#' @param cut cutoff value
#'
#' @return list containing cutoff and title
#'
#' @keywords interlal
get_plot_parameters <- function(y_column, cut) {
  params <- list(
    ViralRefSeq_E = list(cutoff = -log10(cut), ylabel = "-log10 of viral Reference E-values"),
    ViralRefSeq_ident = list(cutoff = NULL, ylabel = "Viral Reference Identity in %"),
    contig_len = list(cutoff = NULL, ylabel = "Contig Length")
  )

  if (!y_column %in% names(params)) {
    stop("Invalid y_column value provided.")
  }

  return(params[[y_column]])
}



#' internal: extract viridae element
#'
#' @param sublist taxonomy column
#'
#'
#' @keywords internal
extract_viridae <- function(sublist) {
  element <- subset(sublist, grepl("viridae", sublist))
  if (length(element) > 0) {
    return(element)
  } else {
    return(NA)  # Return NA if "viridae" is not found in the sublist
  }
}

#' preprocess ViralRefSeq_taxonomy elements
#'
#' @param vh_file hittables file
#'
#' @details
#' Besides best_query the user can utilize the ViralRefSeq_taxonomy column as x_column or groupby
#' in plots. That columns needs to be preprocessed as it it too long and has too many unique elements
#' to be used for grouping. The element containing "viridae" is used because the first one is the tax ID.
#'  NA are replaced by "unclassified" and existing "unclassified" are removed.
#'
#'
#' @return vh_file with preprocessed ViralRefSeq_taxonomy elements
#'
#' @keywords internal
taxonomy_group_preprocess <- function(vh_file){

  # split vh_file.
  my_list <- strsplit(vh_file$ViralRefSeq_taxonomy,split = "|",fixed = TRUE)

  # Apply the function to each sublist
  viridae_elements <- lapply(my_list, extract_viridae)

  # Initialize an empty vector to store the filtered names
  filtered_names <- vector("list", length = length(viridae_elements))

  # Iterate over each element in viridae_elements
  for (i in seq_along(viridae_elements)) {
    # Check if the element contains two or more strings
    if (length(viridae_elements[[i]]) >= 2) {
      # Extract only the virus family names without additional text
      filtered_names[[i]] <- grep("^[[:alnum:]]+viridae$", viridae_elements[[i]], value = TRUE)
    } else {
      # If the element contains less than two strings, keep it unchanged
      filtered_names[[i]] <- viridae_elements[[i]]
    }
  }



  # Remove elements containing "unclassified" outside the function
  viridae_elements <- lapply(filtered_names, function(x) x[!grepl("unclassified", x)])

  viridae_elements[is.na(viridae_elements)] <- "unclassified"

  vh_file$ViralRefSeq_taxonomy <- unlist(viridae_elements)

  return(vh_file)

}

