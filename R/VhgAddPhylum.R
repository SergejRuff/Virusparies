# Function to get the most common taxonomic rank based on suffix
get_most_common_taxonomic_rank <- function(virus_names) {

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

  ranks <- c()
  for (virus_name in virus_names) {
    if (!grepl("unclassified \\+ phylum|unclassified", virus_name, ignore.case = TRUE)) {
      for (j in seq_along(taxa_list)) {
        if (grepl(taxa_list[j], virus_name, ignore.case = TRUE)) {
          ranks <- c(ranks, taxa_names[j])
          break
        }
      }
    }
  }

  if (length(ranks) == 0) {
    return("unclassified")
  }

  most_common_rank <- names(sort(table(ranks), decreasing = TRUE))[1]
  return(most_common_rank)
}



#' @title VhgAddPhylum: extract Phylum information
#'
#' @description
#' VhgAddPhylum adds a `Phylum` column to the provided VirusHunter or VirusGatherer hittable,
#' with each entry in the column reflecting the phylum name derived from the `groupby` column for each observation.
#'
#'
#' @param file A data frame containing VirusHunter or VirusGatherer hittable results.
#' @param groupby (optional): A character specifying the column containing the groups (default: "best_query").
#' Note: Gatherer hittables do not have a "best_query" column. Please provide an appropriate column for grouping.
#'
#' @return Hittable with Phylum column
#'
#' @author Sergej Ruff
#'
#' @seealso
#' VirusHunterGatherer is available here: \url{https://github.com/lauberlab/VirusHunterGatherer}.
#'
#' @examples
#'
#' path <- system.file("extdata", "virushunter.tsv", package = "Virusparies")
#' vh_file <- ImportVirusTable(path)
#'
#' vh_file_filtered <- VhgPreprocessTaxa(vh_file,"Family")
#'
#' processed_taxa <- VhgAddPhylum(vh_file_filtered,"ViralRefSeq_taxonomy")
#'
#' print(unique(processed_taxa$Phylum))
#'
#' @importFrom rlang .data as_string ensym
#' @export
VhgAddPhylum <- function(file,
                         groupby = "best_query"
                         ){

  #is_file_empty(file)
  if (is_file_empty(file)) {
    #message("Skipping VhgBoxplot generation due to empty data.")
    return(invisible(NULL))  # Return invisible(NULL) to stop further execution
  }

  if(rlang::as_string(rlang::ensym(groupby)) != "groupby"){

    groupby <- rlang::as_string(rlang::ensym(groupby))

  }


  if (!(groupby %in% c("best_query", "ViralRefSeq_taxonomy"))) {
    stop('Invalid value for groupby. Please use either "best_query" or "ViralRefSeq_taxonomy".')
  }

  check_input_type(file,groupby,1)


  if (groupby == "ViralRefSeq_taxonomy" && all(grepl("^taxid:", file[[groupby]]))) {

    stop("Error: The 'ViralRefSeq_taxonomy' column doesn't seem to be processed.
Please preprocess this column using the VhgPreprocessTaxa function before using this VhgAddPhylum.")


  }


  if(groupby == "best_query"){
    taxa_rank <- "Family"
  }else{

    taxa_rank <- get_most_common_taxonomic_rank(file[[groupby]])

  }



  if (taxa_rank == "unclassified" && groupby == "ViralRefSeq_taxonomy") {
    # Create the Phylum column based on ViralRefSeq_taxonomy
    file$Phylum <- ifelse(grepl("unclassified", file[[groupby]], ignore.case = TRUE),
                          ifelse(grepl("unclassified\\s+(\\S+)", file[[groupby]], ignore.case = TRUE),
                                 gsub("unclassified\\s+", "", file[[groupby]]),
                                 "unclassified"),
                          NA)

    return(file)
  }




  color_data <- consistentColourPalette(file, groupby = groupby,taxa_rank=taxa_rank)
  legend_labels <- color_data$legend_labels
  labels_manualcol <- color_data$labels
  # matched_vector <- color_data$matched_vector

  # Extract unique values from file$best_query
  unique_queries <- unique(file[[groupby]])

  # Create a vector of corresponding names from legend_labels
  pyhlum_names <- legend_labels[unique_queries]

  # Match names to file$best_query
  file$Phylum <- pyhlum_names[match(file[[groupby]], unique_queries)]

  return(file)


}
