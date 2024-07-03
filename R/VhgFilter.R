#' @title VhgSubsetHittable: Filter VirusHunter and VirusGatherer hittables
#'
#' @description
#' This function filters a virusHunter or VirusGatherer table based on specified criteria,
#' including specific virus groups, minimum number of hits, and observations below certain
#' E-value or identity percentage criteria.
#'
#' @param file A dataframe containing VirusHunter or VirusGatherer hittable data.
#' @param group_column A string indicating the column containing the virus groups specified in the virus_groups argument.
#' Note: Gatherer hittables do not have a "best_query" column. Please provide an appropriate column for grouping.
#' @param virus_groups A character vector specifying virus groups to filter by.
#' @param num_hits_min Minimum number of hits required. Default is NULL, which means no filter based on num_hits.
#' @param ViralRefSeq_E_criteria Maximum E-value threshold for ViralRefSeq_E criteria. Default is NULL, which means no filter based on ViralRefSeq_E.
#' @param ViralRefSeq_ident_criteria Maximum or minimum sequence identity percentage threshold for ViralRefSeq_ident criteria. Default is NULL, which means no filter based on ViralRefSeq_ident.
#'   If positive, filters where ViralRefSeq_ident is above the threshold. If negative, filters where ViralRefSeq_ident is below the absolute value of the threshold.
#' @param contig_len_criteria description
#'
#' @details
#' The function filters the input VirusHunter or VirusGatherer data (`file`) based on specified criteria:
#'
#' - `group_column`: Specifies the column to filter by, which must be either "ViralRefSeq_taxonomy" or "best_query".
#' - `virus_groups`: Allows filtering by specific virus groups. If NULL, all virus groups are included.
#' - `num_hits_min`: Filters rows where the number of hits ("num_hits") is equal to or larger than the specified minimum.
#' - `ViralRefSeq_E_criteria`: Filters rows where the E-value ("ViralRefSeq_E") is below the specified maximum threshold.
#' - `ViralRefSeq_ident_criteria`: Filters rows where the sequence identity percentage ("ViralRefSeq_ident") is above or below the specified threshold.
#'   Use a positive value to filter where ViralRefSeq_ident is above the threshold, and a negative value to filter where ViralRefSeq_ident is below the absolute value of the threshold.
#' -  `contig_len_criteria`: (Gatherer only) Filters rows where the contig length ("contig_len") is above the specified threshold.
#'
#' @return A filtered dataframe based on the specified criteria.
#' @author Sergej Ruff
#'
#'
#' @examples
#'
#' # import file
#' path <- system.file("extdata", "virushunter.tsv", package = "Virusparies")
#' file <- ImportVirusTable(path)
#'
#' # check dimensions before filtering
#' cat("The dimensions of the VirusHunter Hittable before filtering are: \n");dim(file)
#' # we only want the Anello_ORF1core from "best_query" col in this example
#' cat("\nCounts of 'best_query' values before filtering:\n")
#' table(file$best_query)
#' # summary stats
#' # Summary statistics before filtering
#' cat("\nSummary statistics of 'num_hits' column before filtering:\n")
#' summary(file$num_hits)
#' cat("\nSummary statistics of 'ViralRefSeq_ident' column before filtering:\n")
#' summary(file$ViralRefSeq_ident)
#' cat("\nSummary statistics of 'ViralRefSeq_E' column before filtering:\n")
#' summary(file$ViralRefSeq_E)
#'
#' file_filtered <- VhgSubsetHittable(file,group_column = "best_query",
#' virus_groups = "Anello_ORF1core",
#' num_hits_min = 4,ViralRefSeq_ident_criteria = -90,ViralRefSeq_E_criteria = 0.00001)
#'
#' # check dimensions after filtering
#' cat("The dimensions of the VirusHunter Hittable after filtering are: \n");dim(file_filtered)
#' # Summary statistics after filtering
#' cat("\nCounts of 'best_query' values after filtering:\n")
#' table(file_filtered$best_query)
#' cat("\nSummary statistics of 'num_hits' column after filtering:\n")
#' summary(file_filtered$num_hits)
#' cat("\nSummary statistics of 'ViralRefSeq_ident' column after filtering:\n")
#' summary(file_filtered$ViralRefSeq_ident)
#' cat("\nSummary statistics of 'ViralRefSeq_E' column after filtering:\n")
#' summary(file_filtered$ViralRefSeq_E)
#'
#' @seealso
#' VirusHunterGatherer is available here: \url{https://github.com/lauberlab/VirusHunterGatherer}.
#' @export
VhgSubsetHittable <- function(file,
                      group_column="best_query",
                      virus_groups = NULL,
                      num_hits_min = NULL,
                      ViralRefSeq_E_criteria = NULL,
                      ViralRefSeq_ident_criteria = NULL,
                      contig_len_criteria = NULL) {


  # Apply user-defined virus_groups criteria if provided
  if (!is.null(virus_groups)) {

    valid_columns <- c("ViralRefSeq_taxonomy", "best_query")
    if (!(group_column %in% valid_columns)) {
      stop("Error: 'group_column' must be either 'ViralRefSeq_taxonomy' or 'best_query'.")
    }

    unique_groups <- unique(file[[group_column]])
    if (!all(virus_groups %in% unique_groups)) {
      stop("Error: virus_groups contain entries that do not match unique values in group_column.")
    }

    if (!is.character(virus_groups)) {
      stop("Error: 'virus_groups' must be a character vector.")
    }



    file <- file[file[[group_column]] %in% virus_groups, ]
  }

  # Apply user-defined filtering criteria if provided
  if (!is.null(num_hits_min)) {

    if (!"num_hits" %in% colnames(file)) {
      stop("Error: 'num_hits' column does not exist in the dataframe.")
    }
    file <- file[file$num_hits >= num_hits_min, ]
  }
  if (!is.null(ViralRefSeq_E_criteria)) {

    if (!"ViralRefSeq_E" %in% colnames(file)) {
      stop("Error: 'ViralRefSeq_E' column does not exist in the dataframe.")
    }

    file <- file[file$ViralRefSeq_E <= ViralRefSeq_E_criteria, ]
  }
  if (!is.null(ViralRefSeq_ident_criteria)) {
    if (!"ViralRefSeq_ident" %in% colnames(file)) {
      stop("Error: 'ViralRefSeq_ident' column does not exist in the dataframe.")
    }

    if (ViralRefSeq_ident_criteria > 0) {
      # Filter where ViralRefSeq_ident is above the threshold
      file <- file[file$ViralRefSeq_ident > ViralRefSeq_ident_criteria, ]
    } else {
      # Filter where ViralRefSeq_ident is below the absolute value of the threshold
      file <- file[file$ViralRefSeq_ident < abs(ViralRefSeq_ident_criteria), ]
    }
  }

  if (!is.null(contig_len_criteria)) {
    if (!"contig_len" %in% colnames(file)) {
      stop("Error: 'contig_len' column does not exist in the dataframe.")
    }
    file <- file[file$contig_len > contig_len_criteria, ]
  }

  if (nrow(file) == 0) {
    warning("Warning: After filtering, the dataframe is empty.")
  }
  return(file)

}
