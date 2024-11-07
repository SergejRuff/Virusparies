#' @title VhgSubsetHittable: Filter VirusHunter and VirusGatherer hittables
#'
#' @description
#' VhgSubsetHittable filters a VirusHunter or VirusGatherer hittable based on specified criteria,
#' including specific virus groups, minimum number of hits, and observations below certain
#' E-value or identity percentage criteria.
#'
#' @param file A data frame containing VirusHunter or VirusGatherer hittable results.
#' @param group_column A string indicating the column containing the virus groups specified in the virus_groups argument.
#' Note: Gatherer hittables do not have a "best_query" column. Please provide an appropriate column for grouping.
#' @param virus_groups A character vector specifying virus groups to filter by.
#' @param num_hits_min Minimum number of hits required. Default is NULL, which means no filter based on num_hits.
#' @param ViralRefSeq_E_criteria Maximum E-value threshold for ViralRefSeq_E criteria. Default is NULL, which means no filter based on ViralRefSeq_E.
#' @param ViralRefSeq_ident_criteria Maximum or minimum sequence identity percentage threshold for ViralRefSeq_ident criteria. Default is NULL, which means no filter based on ViralRefSeq_ident.
#'   If positive, filters where ViralRefSeq_ident is above the threshold. If negative, filters where ViralRefSeq_ident is below the absolute value of the threshold.
#' @param contig_len_criteria (Gatherer only): Minimum contig length required.
#'
#' @details
#' The function filters the input VirusHunter or VirusGatherer data (`file`) based on specified criteria:
#'
#' - `group_column`: Specifies the column to filter by, which must be either "ViralRefSeq_taxonomy" or "best_query".
#' - `virus_groups`: Allows filtering by specific virus groups. If NULL, all virus groups are included.
#' - `num_hits_min`: Filters rows where the number of hits ("num_hits") is greater than or equal to the specified minimum.
#' - `ViralRefSeq_E_criteria`: Filters rows where the E-value ("ViralRefSeq_E") is below the specified maximum threshold.
#' - `ViralRefSeq_ident_criteria`: Filters rows where the sequence identity percentage ("ViralRefSeq_ident") is above or below the specified threshold.
#'   Use a positive value to filter where ViralRefSeq_ident is above the threshold, and a negative value to filter where ViralRefSeq_ident is below the absolute value of the threshold.
#' -  `contig_len_criteria`: (Gatherer only) Filters rows where the contig length ("contig_len") is greater than or equal to the specified threshold.
#'
#' @return A filtered dataframe based on the specified criteria.
#' @author Sergej Ruff
#'
#'
#' @examples
#'
#' path <- system.file("extdata", "virushunter.tsv", package = "Virusparies")
#' file <- ImportVirusTable(path)
#'
#' cat("The dimensions of the VirusHunter hittable before filtering are: \n");dim(file)
#'
#' file_filtered <- VhgSubsetHittable(file,group_column = "best_query",
#' virus_groups = "Anello_ORF1core",
#' num_hits_min = 4,ViralRefSeq_ident_criteria = -90,ViralRefSeq_E_criteria = 0.00001)
#'
#' cat("The dimensions of the VirusHunter Hittable after filtering are: \n");dim(file_filtered)
#'
#' # other examples for viral_group
#' \donttest{
#'# Include a single group:
#' result1 <- VhgSubsetHittable(file, virus_groups = "Hepadna-Nackedna_TP")
#' # Include multiple groups:
#' result2 <- VhgSubsetHittable(file, virus_groups = c("Hepadna-Nackedna_TP", "Gemini_Rep"))
#' # Exclude a single group:
#' result3 <- VhgSubsetHittable(file, virus_groups = list(exclude = "Hepadna-Nackedna_TP"))
#' # Exclude multiple groups:
#' result4 <- VhgSubsetHittable(file, virus_groups = list(exclude =
#'  c("Hepadna-Nackedna_TP", "Anello_ORF1core")))
#' }
#'
#' @seealso
#' VirusHunterGatherer is available here: \url{https://github.com/lauberlab/VirusHunterGatherer}.
#' @importFrom rlang .data as_string ensym
#' @export
VhgSubsetHittable <- function(file,
                      group_column="best_query",
                      virus_groups = NULL,
                      num_hits_min = NULL,
                      ViralRefSeq_E_criteria = NULL,
                      ViralRefSeq_ident_criteria = NULL,
                      contig_len_criteria = NULL) {


  # Helper function to validate and process group criteria
  process_groups <- function(groups, include_exclude) {
    if (!is.null(groups) && !is.character(groups)) {
      stop(paste("Error: 'virus_groups$", include_exclude, "' must be a character vector."))
    }
    if (is.character(groups)) {
      unique_groups <- unique(file[[group_column]])
      if (!all(groups %in% unique_groups)) {
        stop(paste("Error:", include_exclude, "groups contain entries that do not match unique values in group_column."))
      }
      if (include_exclude == "include") {
        return(file[file[[group_column]] %in% groups, ])
      } else if (include_exclude == "exclude") {
        return(file[!file[[group_column]] %in% groups, ])
      }
    }
    return(file)
  }

  # Apply user-defined virus_groups criteria if provided
  if (!is.null(virus_groups)) {

    group_column <- rlang::as_string(rlang::ensym(group_column))

    valid_columns <- c("ViralRefSeq_taxonomy", "best_query")
    if (!(group_column %in% valid_columns)) {
      stop("Error: 'group_column' must be either 'ViralRefSeq_taxonomy' or 'best_query'.")
    }

    if (is.list(virus_groups)) {
      file <- process_groups(virus_groups$include, "include")
      file <- process_groups(virus_groups$exclude, "exclude")
    } else {
      file <- process_groups(virus_groups, "include")
    }
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

    file <- file[file$ViralRefSeq_E < ViralRefSeq_E_criteria, ]
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
    file <- file[file$contig_len >= contig_len_criteria, ]
  }

  if (nrow(file) == 0) {
    warning("Warning: After filtering, the dataframe is empty.")
  }
  return(file)

}
