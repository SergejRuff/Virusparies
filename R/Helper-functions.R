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
    "ViralRefSeq_E" = paste0("Boxplot plotting reference sequence E-values for specified groups (", x_column, ")\nusing the following cut off: ", cutoff),
    "ViralRefSeq_ident" = paste0("Boxplot plotting reference sequence identity for specified groups (", x_column, ")"),
    "contig_len" = paste0("Boxplot plotting contig length for specified groups (", x_column, ")")
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
    ViralRefSeq_ident = list(cutoff = NULL, ylabel = "Reference sequence identity (%)"),
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

  # Check if ICTV_env exists and is an environment
  if (exists("ICTV_env") && is.environment(ICTV_env) && exists("new_ICTV_data", envir = ICTV_env)) {

    # Assign new_ICTV_data to df from ICTV_env
    df <- get("new_ICTV_data", envir = ICTV_env)

    message("Warning: User-defined ICTV data is being used in place of the internal default ICTV data.")

  } else {

    # Assign ICTV_data to df
    df <- ICTV_data
  }



  return(df %>%
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
                        "Lenarviricota","Negarnaviricota","Pisuviricota")


  valid_phyla_smalldna <-  c("Hofneiviricota", "Phixviricota", "Cossaviricota",
                             "Cressdnaviricota", "Saleviricota")

  valid_phyla_largedna <-  c("Taleaviricota", "Nucleocytoviricota", "Preplasmiviricota",
                              "Dividoviricota", "Peploviricota", "Uroviricota")

  all_phyla <- c(valid_phyla_rna, valid_phyla_smalldna, valid_phyla_largedna)

  # Add families to exclude for each group

  if(groupby == "best_query"){
    families_rna <- c("^Birna", "^Permutotetra")
    families_smalldna <- c("^Anello")
    families_largedna <- c("^Yara")
    other <- c("^Yara","^Anello","^Birna","^Permutotetra")

  }else{
    families_rna <- c("Birnaviridae", "Permutotetraviridae")
    families_smalldna <- c("Anelloviridae")
    families_largedna <- c("Yaraviridae")
    other <- c("Yaraviridae","Anelloviridae","Birnaviridae","Permutotetraviridae")

  }


  chosen_list <- switch(chosen_group,
                        "rna" = valid_phyla_rna,
                        "smalldna" = valid_phyla_smalldna,
                        "largedna" = valid_phyla_largedna,
                        "others" = NULL,
                        stop("Invalid chosen_group value. Use 'rna', 'smalldna', 'largedna', or 'others'."))


  change_label <- switch(chosen_group,
                         "rna" = "Non-RNA-viruses",
                         "smalldna" = "Non-Small-DNA-Viruses",
                         "largedna" = "Non-Large-DNA-Viruses",
                         "others" = "Other Viruses")




  non_group_families <- switch(chosen_group,
                               "rna" = families_rna,
                               "smalldna" = families_smalldna,
                               "largedna" = families_largedna)

  file <- file %>%
    mutate(
      phyl = case_when(
        grepl(paste0(non_group_families, collapse = "|"), !!sym(groupby), ignore.case = TRUE) ~ phyl,
        phyl == "unclassified" & !!sym(groupby) != "unclassified" ~ change_label,
        TRUE ~ phyl
      ),
      !!sym(groupby) := case_when(
        grepl(paste0(non_group_families, collapse = "|"), !!sym(groupby), ignore.case = TRUE) ~ !!sym(groupby),
        phyl == "unclassified" & !!sym(groupby) != "unclassified" ~ change_label,
        TRUE ~ !!sym(groupby)
      )
    )


  # Step 2: Handle the 'others' category
  if (is.null(chosen_list)) {
    file <- file %>%
      mutate(!!sym(groupby) := case_when(
        !phyl %in% c(valid_phyla_rna, valid_phyla_smalldna, valid_phyla_largedna) ~ change_label,
        TRUE ~ !!sym(groupby)
      ),
        phyl = case_when(
          !phyl %in% c(valid_phyla_rna, valid_phyla_smalldna, valid_phyla_largedna) ~ change_label,
          TRUE ~ phyl
        )

      )

  } else {
    file <- file %>%
      mutate(
        !!sym(groupby) := case_when(
          phyl == "unclassified" ~ !!sym(groupby),
          grepl(paste0(chosen_list, collapse = "|"), .data$phyl, ignore.case = TRUE) ~ !!sym(groupby),
          TRUE ~ change_label
        ),
        phyl = case_when(
          phyl == "unclassified" ~ "unclassified",
          grepl(paste0(chosen_list, collapse = "|"), phyl, ignore.case = TRUE) ~ phyl,
          TRUE ~ change_label
        )
      )
  }




 if(!is.null(label_vector)){

   # Filter label_vector based on chosen_group
   selected_phyla <- if (chosen_group == "others") {
     all_phyla
   } else {
     chosen_list
   }

   label_vector <- label_vector[names(label_vector) %in% c(selected_phyla, "unclassified")]

   # print(label_vector)

   # Add the "Non-RNA-virus" entry with the color black
   label_vector[change_label] <- "#000000"

   return(list(file =file,label=label_vector))


 }else{

   return(list(file =file))
 }





}


#' @title New_ICTV: Assign Custom ICTV Data for Use in Virusparies
#'
#' @description
#' New_ICTV enables users to define their own ICTV data for use in Virusparies.
#'
#'
#' @param new_ICTV_data data frame containing new ICTV data.
#'
#' @return An `ICTV_env` object containing the user-defined ICTV data
#'
#' @details
#' Virusparies utilizes the ICTV taxonomy data set for two main purposes:
#' 1. Assigning the appropriate taxonomy rank via VhgPreprocessTaxa and related functions.
#' 2. In plotting functions, the ICTV data set is used to assign colors based on phylum.
#'
#' `New_ICTV` enables users to define and load their own ICTV data for use in Virusparies.
#' This feature is useful for those who wish to use a subset of ICTV data to speed up functions
#' such as \code{\link{VhgPreprocessTaxa}}, or for users who want to incorporate a newer or older version of
#' ICTV data within Virusparies.
#' The custom ICTV data is stored in a new environment, which means it does not overwrite the
#' internal data. Instead, users need to apply this function each time they start a new R session
#' to load their custom dataset.
#'
#' @examples
#' # Define example data
#'
#' # Sample data
#' Example_ICTV <- data.frame(
#'Phylum = c("Taleaviricota", "Taleaviricota", "Taleaviricota", "Taleaviricota",
#'"Taleaviricota"),Subphylum = c(NA, NA, NA, NA, NA),
#' Class = c("Tokiviricetes", "Tokiviricetes", "Tokiviricetes", "Tokiviricetes",
#' "Tokiviricetes"),
#' Subclass = c(NA, NA, NA, NA, NA),
#' Order = c("Ligamenvirales", "Ligamenvirales", "Ligamenvirales", "Ligamenvirales",
#' "Ligamenvirales"),
#' Suborder = c(NA, NA, NA, NA, NA),
#' Family = c("Lipothrixviridae", "Lipothrixviridae", "Lipothrixviridae",
#' "Lipothrixviridae", "Lipothrixviridae"),
#' Subfamily = c(NA, NA, NA, NA, NA),
#' Genus = c("Alphalipothrixvirus", "Alphalipothrixvirus",
#' "Betalipothrixvirus", "Betalipothrixvirus", "Betalipothrixvirus"),
#' Subgenus = c(NA, NA, NA, NA, NA),
#' stringsAsFactors = FALSE)
#'
#' \donttest{
#' # Assign new example ICTV data for use in Virusparies
#' New_ICTV(Example_ICTV)
#' # check currently used ICTV data in Virusparies
#' Current_ICTV()
#'
#' }
#'
#'
#'
#'
#' @author Sergej Ruff
#' @seealso
#' VirusHunterGatherer is available here: \url{https://github.com/lauberlab/VirusHunterGatherer}.
#' @export
New_ICTV <- function(new_ICTV_data=NULL){

  pos <- 0



  check_is_dataframe(new_ICTV_data)

  # Define the required columns
  required_columns <- c("Phylum", "Subphylum", "Class", "Subclass",
                        "Order", "Suborder", "Family", "Subfamily",
                        "Genus", "Subgenus")

  # Check if all required columns are present
  missing_columns <- setdiff(required_columns, colnames(new_ICTV_data))

  if (length(missing_columns) > 0) {
    stop("new_ICTV_data is missing the following required columns: ",
         paste(missing_columns, collapse = ", "))
  }

  ICTV_env <- new.env()

  assign("new_ICTV_data",new_ICTV_data,envir = ICTV_env)

  if(!is.null(new_ICTV_data)){

    pos <- 1

    assign("ICTV_env", ICTV_env, envir = as.environment(pos))

  }


}


#' @title Current_ICTV: Verify the version of ICTV data currently in use.
#'
#' @description
#' Displays the entire ICTV data object currently used by Virusparies.
#'
#'
#' @return data frame with the currently used ICTV data
#'
#' @details
#' Virusparies utilizes the ICTV taxonomy data set for two main purposes:
#' 1. Assigning the appropriate taxonomy rank via VhgPreprocessTaxa and related functions.
#' 2. In plotting functions, the ICTV data set is used to assign colors based on phylum.
#'
#' Current_ICTV allows users to inspect the current ICTV data object used by Virusparies.
#' It also enables users to check if their custom version of the ICTV data,
#' provided by the \code{\link{New_ICTV}}, is loaded.
#'
#' @examples
#'
#' internal_ICTV <- Current_ICTV()
#'
#' str(internal_ICTV)
#'
#'
#' @author Sergej Ruff
#' @seealso
#' VirusHunterGatherer is available here: \url{https://github.com/lauberlab/VirusHunterGatherer}.
#' @export
Current_ICTV <- function(){

  # Check if ICTV_env exists and is an environment
  if (exists("ICTV_env") && is.environment(ICTV_env) && exists("new_ICTV_data", envir = ICTV_env)) {

    message("Currently, the User-defined ICTV data is being used.")

    # Assign new_ICTV_data to df from ICTV_env
    df <- get("new_ICTV_data", envir = ICTV_env)

  } else {

    message("Currently, the internal ICTV data is being used.")

    # Assign ICTV_data to df
    df <- ICTV_data
  }

  return(df)
}

