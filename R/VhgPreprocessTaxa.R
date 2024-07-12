#' define taxa string
#'
#' @param taxa_rank a character indicating the taxa rank
#'
#' @return a character with the taxa suffix
#'
#' @noRd
taxonomy_rank_hierarchy <- function(taxa_rank) {

  # Define the list with values
  taxa_list <- list(
    "viricotina",
    "viricetes",
    "viricetidae",
    "virales",
    "virineae",
    "viridae",
    "virinae",
    "virus"
  )

  # Define the names for each element
  taxa_names <- c(
    "Subphylum",
    "Class",
    "Subclass",
    "Order",
    "Suborder",
    "Family",
    "Subfamily",
    "Genus"
  )

  # Assign names to the list elements using setNames()
  taxa_list <- setNames(taxa_list, taxa_names)

  valid_ranks <- names(taxa_list)

  if (!(taxa_rank %in% valid_ranks)) {
    stop("Error: Invalid taxa rank provided. Please provide one of: Subphylum, Class, Subclass, Order, Suborder, Family, Subfamily, Genus")
  }

  return(taxa_list[[taxa_rank]])
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
#' @importFrom tidyr unnest pivot_longer
#' @importFrom stringr str_split str_extract
#' @export
VhgPreprocessTaxa <- function(vh_file,taxa_rank) {

  if (!all(grepl("^taxid:", vh_file$ViralRefSeq_taxonomy))) {

    message("The 'ViralRefSeq_taxonomy' column is expected to start with 'taxid:' followed by taxa ranks separated by '|'.\n",
            "However, Your column has a different structure, possibly indicating that the data has already been processed by 'VhgPreprocessTaxa'.\n",
            "Skipping Taxonomy processing ...")

    return(vh_file)

  }

  taxa_rank <- taxonomy_rank_hierarchy(taxa_rank)



  ictv_formatted <- format_ICTV()



  vh_file <- vh_file %>%
    mutate(
      ViralRefSeq_taxonomy_split = str_split(.data$ViralRefSeq_taxonomy, "\\|"),
    ) %>%
    unnest(.data$ViralRefSeq_taxonomy_split) %>%
    left_join(ictv_formatted, join_by("ViralRefSeq_taxonomy_split" == "name")) %>%
    mutate(
      ViralRefSeq_taxonomy = case_when(
        !is.na(str_extract(.data$ViralRefSeq_taxonomy, paste0("\\w+", taxa_rank))) ~ str_extract(.data$ViralRefSeq_taxonomy, paste0("\\w+", taxa_rank)),
        !is.na(.data$Phylum) ~ paste("unclassified", .data$Phylum),
        .default = "unclassified"
      )
    ) %>%
    select(!c(.data$ViralRefSeq_taxonomy_split:.data$level))%>%
    group_by(across(-.data$ViralRefSeq_taxonomy)) %>%  # Specify your grouping columns here
    slice(which.min(.data$ViralRefSeq_taxonomy == "unclassified")) %>%
    ungroup() %>%
    distinct()







  vh_file <- vh_file %>%
    mutate(ViralRefSeq_taxonomy = if_else(.data$ViralRefSeq_taxonomy == "unclassified unclassified",
                                          "unclassified", .data$ViralRefSeq_taxonomy)) %>%
    mutate(ViralRefSeq_taxonomy = if_else(.data$ViralRefSeq_taxonomy == "unclassified NA",
                                          "unclassified", .data$ViralRefSeq_taxonomy))


  return(vh_file)
}
