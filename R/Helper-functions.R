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
    ViralRefSeq_E = list(cutoff = -log10(cut), ylabel = "-log10 of viral reference e-values"),
    ViralRefSeq_ident = list(cutoff = NULL, ylabel = "Viral reference identity (%)"),
    contig_len = list(cutoff = NULL, ylabel = "Contig length (nt)")
  )

  if (!y_column %in% names(params)) {
    stop("Invalid y_column value provided.")
  }

  return(params[[y_column]])
}




#' define taxa string
#'
#' @param taxa_rank a character indicating the taxa rank
#'
#' @return a character with the taxa suffix
#'
#' @keywords internal
taxonomy_rank_hierarchy<- function(taxa_rank){


  # Define valid taxa ranks
  valid_ranks <- c("Subphylum", "Class", "Subclass", "Order",
                   "Suborder", "Family", "Subfamily", "Genus")

  # Check if the provided taxa_rank is valid
  if (!(taxa_rank %in% valid_ranks)) {
    stop("Error: Invalid taxa rank provided. Please provide one of: Subphylum, Class, Subclass, Order, Suborder, Family, Subfamily, Genus")
  }


  if(taxa_rank=="Subphylum"){
    return("viricotina")
  }

  if(taxa_rank=="Class"){
    return("viricetes")
  }

  if(taxa_rank=="Subclass"){
    return("viricetidae")
  }

  if(taxa_rank=="Order"){
    return("virales")
  }

  if(taxa_rank=="Suborder"){
    return("virineae")
  }

  if(taxa_rank=="Family"){
    return("viridae")
  }

  if(taxa_rank=="Subfamily"){
    return("virinae")
  }

  if(taxa_rank=="Genus"){
    return("virus")
  }


}

#' internal: extract viridae element
#'
#' @param sublist taxonomy column
#'
#'
#' @keywords internal
extract_taxa <- function(sublist, taxa_rank) {
  element <- subset(sublist, grepl(taxa_rank, sublist))
  if (length(element) > 0) {
    return(element)
  } else {
    return(NA)  # Return NA if "viridae" is not found in the sublist
  }
}

#' @title VhgPreprocessTaxa: preprocess ViralRefSeq_taxonomy elements
#'
#' @details
#' Process the `ViralRefSeq_taxonomy` column.
#'
#'
#' @param vh_file VirusHunter or VirusGatherer hittable
#' @param taxa_rank Specify the taxonomic rank to group your data by.
#' Supported ranks are:
#' - "Subphylum"
#' - "Class"
#' - "Subclass"
#' - "Order"
#' - "Suborder"
#' - "Family" (default)
#' - "Subfamily"
#' - "Genus" (including Subgenus)
#'
#' @details
#' Besides `best_query`, the user can utilize the `ViralRefSeq_taxonomy` column as `x_column` or `groupby` in plots.
#' This column needs preprocessing because it is too long and has too many unique elements for effective grouping.
#' The element containing the taxa suffix specified by the `taxa_rank` argument is used. NA values are replaced by "unclassified".
#'
#' This function is used internally by every function that can use the `ViralRefSeq_taxonomy` column as input.
#' The user can also apply this function independently to process the taxonomy column and filter for the selected taxa rank in their data.
#'
#'
#' @return vh_file with preprocessed ViralRefSeq_taxonomy elements
#' @author Sergej Ruff
#' @examples
#' path <- system.file("extdata", "virushunter.tsv", package = "Virusparies")
#' vh_file <- ImportVirusTable(path)
#'
#' vh_file_filtered <- VhgPreprocessTaxa(vh_file,"Family")
#'
#' print("ViralRefSeq_taxonomy before processing:\n")
#' print(head(vh_file$ViralRefSeq_taxonomy,5))
#'
#' print("ViralRefSeq_taxonomy after processing:\n")
#' print(head(vh_file_filtered$ViralRefSeq_taxonomy,5))
#'
#'
#'
#' @seealso
#' VirusHunterGatherer is available here: \url{https://github.com/lauberlab/VirusHunterGatherer}.
#'
#' @import dplyr
#' @importFrom rlang .data
#' @export
VhgPreprocessTaxa <- function(vh_file,taxa_rank) {

  if (!all(grepl("^taxid:", vh_file$ViralRefSeq_taxonomy))) {

    message("The 'ViralRefSeq_taxonomy' column is expected to start with 'taxid:' followed by taxa ranks separated by '|'.\n",
            "However, Your column has a different structure, possibly indicating that the data has already been processed by 'VhgPreprocessTaxa'.\n",
            "Skipping Taxonomy processing ...")

    return(vh_file)

  }

  taxa_rank <- taxonomy_rank_hierarchy(taxa_rank)

  # Split vh_file.
  my_list <- strsplit(vh_file$ViralRefSeq_taxonomy, split = "|", fixed = TRUE)

  # Apply the function to each sublist
  taxa_elements <- lapply(my_list, function(x) extract_taxa(x, taxa_rank))

  # Initialize an empty vector to store the filtered names
  filtered_names <- vector("list", length = length(taxa_elements))

  # Iterate over each element in taxa_elements
  for (i in seq_along(taxa_elements)) {
    # Check if the element contains two or more strings
    if (length(taxa_elements[[i]]) >= 2) {
      # Extract only the virus family names without additional text
      term <- paste0("^[[:alnum:]]+",taxa_rank,"$")
      filtered_names[[i]] <- grep(term, taxa_elements[[i]], value = TRUE)
    } else {
      # If the element contains less than two strings, keep it unchanged
      filtered_names[[i]] <- taxa_elements[[i]]
    }
  }

  # Load ICTV data (assuming the ICTV data is properly loaded and formatted)
  ICTV_data <- ICTV_data


  # Initialize a vector to store the processed taxonomy names
  processed_taxonomy <- vector("character", length = length(unlist(filtered_names)))

  for (i in seq_along(filtered_names)) {
    entry <- filtered_names[[i]]

    if (is.na(entry) || length(entry) == 0) {
      # Find matching rows in ICTV_data based on various taxonomic levels
      match_row <- ICTV_data %>%
        filter(.data$Subphylum %in% my_list[[i]] |
                 .data$Class %in% my_list[[i]] |
                 .data$Subclass %in% my_list[[i]] |
                 .data$Order %in% my_list[[i]] |
                 .data$Suborder %in% my_list[[i]] |
                 .data$Family %in% my_list[[i]])

      if (nrow(match_row) > 0) {
        phylum <- match_row$Phylum[1]
        processed_taxonomy[i] <- paste("unclassified", phylum)
      } else {
        processed_taxonomy[i] <- "unclassified"
      }
    } else {
      processed_taxonomy[i] <- entry
    }
  }



  # Update the ViralRefSeq_taxonomy column in vh_file
  vh_file$ViralRefSeq_taxonomy <- unlist(processed_taxonomy)

  # Using base R
  vh_file$ViralRefSeq_taxonomy <- ifelse(vh_file$ViralRefSeq_taxonomy ==
                                           "unclassified unclassified", "unclassified", vh_file$ViralRefSeq_taxonomy)

  vh_file$ViralRefSeq_taxonomy <- ifelse(vh_file$ViralRefSeq_taxonomy ==
                                           "unclassified NA", "unclassified", vh_file$ViralRefSeq_taxonomy)


  return(vh_file)
}


#' remove group text
#'
#' @param plot plot
#' @param remove_x_axis_labels remove_x_axis_labels
#' @param flip_coords flip_coords
#'
#'
#'
#' @keywords internal
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
#' @keywords internal
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
