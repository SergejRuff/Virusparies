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
#' @noRd
is_file_empty <- function(file) {
  if (nrow(file) == 0) {
    message("Warning: Input has zero rows. Skipping further processing.")
    return(TRUE)  # Return TRUE to indicate empty data
  }
  return(FALSE)  # Return FALSE if data is not empty
}


#' internal function: message to be printed in boxplot depending on chosen y-column
#'
#' @param y_column y_column
#' @param x_column x_column
#' @param cutoff cutoff
#'
#' @author Sergej Ruff
#' @noRd
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
#' @noRd
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
#' @noRd
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





#' checks if correct datatype is provided to function (character check)
#'
#' @param name argument from function.
#'
#' @author Sergej Ruff
#' @noRd
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
#' @noRd
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
#' @noRd
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
#' @noRd
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
#' @noRd
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
#' @noRd
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
    ViralRefSeq_E = list(cutoff = -log10(cut), ylabel = "-log10 of viral reference E-values"),
    ViralRefSeq_ident = list(cutoff = NULL, ylabel = "Viral reference identity (%)"),
    contig_len = list(cutoff = NULL, ylabel = "Contig length (nt)")
  )

  if (!y_column %in% names(params)) {
    stop("Invalid y_column value provided.")
  }

  return(params[[y_column]])
}






#' remove group text
#'
#' @param plot plot
#' @param remove_x_axis_labels remove_x_axis_labels
#' @param flip_coords flip_coords
#'
#'
#'
#' @noRd
remove_group_text <- function(plot,remove_x_axis_labels,flip_coords){

  if(remove_x_axis_labels && flip_coords){
    plot <- plot + theme(axis.text.y = element_blank())

  } else if(remove_x_axis_labels && !flip_coords){

    plot <- plot + theme(axis.text.x = element_blank())
  }

  return(plot)

}


#' internal function filter_specific_group
#'
#' @param file a hittable file
#' @param groupby groupby column
#' @param filter_group_criteria string,number or vector of character or number to filte by.
#'
#' @return  preprocessed hittable
#'
#' @noRd
filter_specific_group <- function(file,groupby, filter_group_criteria) {

  # Check if filter_group_criteria is NULL
  if (is.null(filter_group_criteria)) {
    return(file)
  }


  if (!is.atomic(filter_group_criteria)) {
    stop("Error: filter_group_criteria must be a vector, single character, or single numeric value.")
  }



  # Check if filter_group_criteria is character
  if (is.character(filter_group_criteria)) {
    file <- file[file[[groupby]] %in% filter_group_criteria, ]
    return(file)
  }

  # Check if filter_group_criteria is numeric
  if (is.numeric(filter_group_criteria)) {
    # Get unique values of the grouping variable
    unique_groups <- unique(file[[groupby]])

    # Check if any number in filter_group_criteria exceeds the number of groups
    if (any(filter_group_criteria > length(unique_groups))) {
      stop("Error: filter_group_criteria contains numbers larger than the amount of groups.")
    }

    # Extract the names of the specified groups based on numeric indices
    taxonomies_of_interest <- unique_groups[filter_group_criteria]

    # Filter the dataframe for rows where groupby matches the specified categories
    file <- file[file[[groupby]] %in% taxonomies_of_interest, ]

    return(file)
  }

  # Return an error if filter_group_criteria is not character or numeric
  stop("Error: filter_group_criteria must be either character or numeric.")
}


#' facet_plot
#'
#' @param plot plot object
#' @param facet_ncol number of columns to facet by. Default NULL
#' @param flip_coords flip_coords
#'
#' @return plot
#'
#' @noRd
facet_plot <- function(plot,facet_ncol=FALSE,flip_coords=TRUE){

  if(!is.null(facet_ncol)){

    plot <- plot +  # Define colors for TRUE and FALSE
      facet_wrap(~.data$phyl, ncol = facet_ncol,scales = ifelse(flip_coords, "free_y", "free_x"))
  }

  return(plot)
}



#' pivort ICVT Data
#'
#' @param taxa_rank taxa_rank
#'
#' @return pivotted ICTV_data
#'
#' @importFrom tidyr pivot_longer
#' @noRd
format_ICTV <- function(taxa_rank){



  return(ICTV_data %>%
            select(.data$Phylum:.data$Subgenus) %>%
            pivot_longer(.data$Subphylum:.data$Subgenus, names_to = "level", values_to = "name") %>%
            filter(str_detect(.data$name, paste0("\\w+", taxa_rank), negate = TRUE)) %>%
            na.omit() %>%
            distinct())



}

#' change angle of labels
#'
#' @param plot plot obj
#' @param x_angle x angle
#' @param y_angle y angle
#'
#' @return plot
#'
#' @import ggplot2
#' @noRd
adjust_plot_angles <- function(plot, x_angle = NULL, y_angle = NULL) {
  # Check if x_angle is provided, and if so, adjust the x-axis text angle
  if (!is.null(x_angle)) {
    plot <- plot + theme(axis.text.x = element_text(angle = x_angle, vjust = 0.5, hjust = 1))
  }

  # Check if y_angle is provided, and if so, adjust the y-axis text angle
  if (!is.null(y_angle)) {
    plot <- plot + theme(axis.text.y = element_text(angle = y_angle, vjust = 0.5, hjust = 1))
  }

  # Return the modified plot
  return(plot)
}




#' remove_non_group
#'
#' @param file file
#' @param groupby groupby
#' @param chosen_group chosen_group
#' @param label_vector label_vector
#' @param taxa_rank taxa_rank
#'
#' @return list with vector and df
#' @import dplyr
#' @importFrom rlang .data
#' @noRd
remove_non_group <- function(file,groupby,chosen_group,label_vector,taxa_rank){

  valid_phyla_rna <-  c("Ambiviricota","Duplornaviricota","Kitrinoviricota",
                        "Lenarviricota","Negarnaviricota ","Pisuviricota","Artverviricota")


  valid_phyla_smalldna <-  c("Ambiviricota","Duplornaviricota","Kitrinoviricota",
                             "Lenarviricota","Negarnaviricota ","Pisuviricota","Artverviricota")

  valid_phyla_largedna <-  c("Ambiviricota","Duplornaviricota","Kitrinoviricota",
                             "Lenarviricota","Negarnaviricota ","Pisuviricota","Artverviricota")

  chosen_list <- switch(chosen_group,
                        "rna" = valid_phyla_rna,
                        "smalldna" = valid_phyla_smalldna,
                        "largedna" = valid_phyla_largedna,
                        stop("Invalid chosen_group value. Use 'rna', 'smalldna', or 'largedna'."))

  change_label <- switch(chosen_group,
                         "rna" = "Non-RNA-virus",
                         "smalldna" = "Non-Small-DNA-Virus",
                         "largedna" = "Non-Large-DNA-Virus")
  # if (all(grepl("^taxid:", file$ViralRefSeq_taxonomy))) {
  #   file <- VhgPreprocessTaxa(file,taxa_rank)
  # }
  #
  #
  # file <- VhgAddPhylum(file,"ViralRefSeq_taxonomy")

  # Modify the dataframe
  file <- file %>%
    mutate(
      !!sym(groupby) := case_when(
        phyl == "unclassified" ~ !!sym(groupby),
        grepl(paste0(chosen_list, collapse = "|"), .data$phyl, ignore.case = TRUE) ~ !!sym(groupby),
        TRUE ~ change_label
      )
    )

  file <- file %>%
    mutate(
      phyl = case_when(
        phyl == "unclassified" ~ "unclassified",  # Keep as unclassified
        grepl(paste0(chosen_list, collapse = "|"), phyl, ignore.case = TRUE) ~ phyl,  # Keep original value for valid phyla
        TRUE ~ change_label  # Change to label if not in chosen_list
      )
    )

  # Filter out entries that are not in the valid phyla list, excluding "unclassified"
  label_vector <- label_vector[names(label_vector) %in% c(chosen_list, "unclassified")]

  # Add the "Non-RNA-virus" entry with the color black
  label_vector[change_label] <- "#000000"

  return(list(file =file,label=label_vector))


}
