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

# Define the function to process each chunk
process_chunk <- function(chunk, ictv_formatted, taxa_rank) {
  taxon_filter <- paste(unique(ictv_formatted$name), collapse = "|")

  chunk_processed <- chunk %>%
    mutate(
      ViralRefSeq_taxonomy = str_remove_all(.data$ViralRefSeq_taxonomy, "taxid:\\d+\\||\\w+\\s\\w+\\|"),
      name = str_extract(.data$ViralRefSeq_taxonomy, taxon_filter),
      ViralRefSeq_taxonomy = str_extract(.data$ViralRefSeq_taxonomy, paste0("\\w+", taxa_rank))
    ) %>%
    left_join(ictv_formatted, join_by("name" == "name")) %>%
    mutate(
      ViralRefSeq_taxonomy = case_when(
        is.na(.data$ViralRefSeq_taxonomy) & is.na(.data$Phylum) ~ "unclassified",
        is.na(.data$ViralRefSeq_taxonomy) ~ paste("unclassified", .data$Phylum),
        .default = .data$ViralRefSeq_taxonomy
      )
    ) %>% select(-c(.data$name:.data$level)) %>%
    mutate(
      ViralRefSeq_taxonomy = if_else(.data$ViralRefSeq_taxonomy == "unclassified unclassified", "unclassified", .data$ViralRefSeq_taxonomy),
      ViralRefSeq_taxonomy = if_else(.data$ViralRefSeq_taxonomy == "unclassified NA", "unclassified", .data$ViralRefSeq_taxonomy)
    )

  return(chunk_processed)
}


#' @title VhgPreprocessTaxa: preprocess ViralRefSeq_taxonomy elements
#'
#' @details
#' Process the `ViralRefSeq_taxonomy` column.
#'
#'
#' @param file A data frame containing VirusHunter or VirusGatherer hittable results.
#' @param taxa_rank (optional): Specify the taxonomic rank to group your data by.
#' Supported ranks are:
#' - "Subphylum"
#' - "Class"
#' - "Subclass"
#' - "Order"
#' - "Suborder"
#' - "Family" (default)
#' - "Subfamily"
#' - "Genus" (including Subgenus)
#' @param num_cores (optional): Number of cores used (default: 1)
#'
#'
#' @details
#' Besides `best_query`, the user can utilize the `ViralRefSeq_taxonomy` column as `x_column` or `groupby` in plots.
#' This column needs preprocessing because it is too long and has too many unique elements for effective grouping.
#' The element containing the taxa suffix specified by the `taxa_rank` argument is used. NA values are replaced by "unclassified".
#'
#' This function is used internally by every function that can use the `ViralRefSeq_taxonomy` column as input.
#' The user can also apply this function independently to process the taxonomy column and filter for the selected taxa rank in their data.
#'
#' For datasets with significantly more than 100,000 observations, it is recommended to use this function to ensure it is skipped in the plot functions.
#'
#' The `num_cores` parameter allows you to divide the dataset into multiple parts, corresponding to the number of cores available.
#' This enables parallel processing across multiple threads, thereby speeding up the overall processing time.
#'
#' @return file with preprocessed ViralRefSeq_taxonomy elements
#' @author Sergej Ruff
#' @examples
#' path <- system.file("extdata", "virushunter.tsv", package = "Virusparies")
#' file <- ImportVirusTable(path)
#'
#' file_filtered <- VhgPreprocessTaxa(file,"Family")
#'
#' cat("ViralRefSeq_taxonomy before processing:\n")
#' print(head(file$ViralRefSeq_taxonomy,5))
#'
#' cat("ViralRefSeq_taxonomy after processing:\n")
#' print(head(file_filtered$ViralRefSeq_taxonomy,5))
#'
#'
#'
#' @seealso
#' VirusHunterGatherer is available here: \url{https://github.com/lauberlab/VirusHunterGatherer}.
#'
#' @import dplyr
#' @importFrom rlang .data as_string ensym
#' @importFrom tidyr unnest pivot_longer
#' @importFrom stringr  str_extract str_remove_all str_detect
#' @importFrom parallel detectCores mclapply
#' @export
VhgPreprocessTaxa <- function(file,taxa_rank, num_cores = 1) {

  #is_file_empty(file)
  if (is_file_empty(file)) {
    #message("Skipping VhgBoxplot generation due to empty data.")
    return(invisible(NULL))  # Return invisible(NULL) to stop further execution
  }

  if (!all(grepl("^taxid:", file$ViralRefSeq_taxonomy))) {

    message("The 'ViralRefSeq_taxonomy' column is expected to start with 'taxid:' followed by taxa ranks separated by '|'.\n",
            "However, Your column has a different structure, possibly indicating that the data has already been processed by 'VhgPreprocessTaxa'.\n",
            "Skipping Taxonomy processing ...")

    return(file)

  }

  if(rlang::as_string(rlang::ensym(taxa_rank)) != "taxa_rank"){

    taxa_rank <- rlang::as_string(rlang::ensym(taxa_rank))


  }



  taxa_rank <- taxonomy_rank_hierarchy(taxa_rank)



  ictv_formatted <- format_ICTV(taxa_rank)


  # Validate num_cores
  if (num_cores < 1) {
    stop("Number of cores must be at least 1")
  }

  # Use default number of cores if specified is greater than available
  num_cores <- min(num_cores, detectCores() - 1)

  # Split the data into chunks
  num_chunks <- num_cores
  chunks <- split(file, rep(1:num_chunks, length.out = nrow(file)))

  # Process chunks in parallel
  processed_chunks <- mclapply(chunks, process_chunk, ictv_formatted = ictv_formatted, taxa_rank = taxa_rank, mc.cores = num_cores)

  # Combine processed chunks
  file_processed <- bind_rows(processed_chunks)

  return(file_processed)
}
