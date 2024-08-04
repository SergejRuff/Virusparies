#' VhgGetSubject: Process and Count Viral Subjects within Groups
#'
#' @param file A data frame containing VirusHunter or VirusGatherer hittable results.
#' @param groupby (optional): A character specifying the column containing the groups (default: "best_query").
#' Note: Gatherer hittables do not have a "best_query" column. Please provide an appropriate column for grouping.
#' @param remove_identifiers (optional): if `TRUE` (default), removes the identifiers in the ViralRefSeq_subject cells.
#'
#' @details
#' The function `VhgGetSubject` counts the number of viral subjects in the `ViralRefSeq_subject` column
#' for each group specified by the `groupby` argument.
#' It returns a tibble with three columns: the first column contains the viral group specified by the `groupby` argument,
#' the second column lists the viral subjects found in that group, and the third column shows how many times each viral subject appears in that group.
#'
#'
#' @return a processed tibble object.
#' @author Sergej Ruff
#'
#' @seealso
#' VirusHunterGatherer is available here: \url{https://github.com/lauberlab/VirusHunterGatherer}.
#'
#' @examples
#'
#' # import data
#' path <- system.file("extdata", "virushunter.tsv", package = "Virusparies")
#' file <- ImportVirusTable(path)
#'
#' # process column and filter for significant groups
#' file <- VhgPreprocessTaxa(file,taxa_rank = "Family")
#' file_filtered <- VhgSubsetHittable(file,ViralRefSeq_E_criteria = 1e-5)
#'
#' subject_df <- VhgGetSubject(file_filtered,groupby = "ViralRefSeq_taxonomy",
#' remove_identifiers = TRUE)
#'
#' print(subject_df)
#'
#' @import dplyr
#' @import stringr
#' @importFrom rlang .data
#' @export
VhgGetSubject <- function(file,
                          groupby = "best_query",
                          remove_identifiers = TRUE){

  #is_file_empty(file)
  if (is_file_empty(file)) {
    #message("Skipping VhgBoxplot generation due to empty data.")
    return(invisible(NULL))  # Return invisible(NULL) to stop further execution
  }

  if (!(groupby %in% c("best_query", "ViralRefSeq_taxonomy"))) {
    stop('Invalid value for groupby. Please use either "best_query" or "ViralRefSeq_taxonomy".')
  }

  arg_logical(remove_identifiers)

  # Process the data based on the remove_identifiers argument
  result <- file %>%
    mutate(Processed_ViralRefSeq_subject = if (remove_identifiers) {
      str_extract(.data$ViralRefSeq_subject, "(?<=\\|).*")
    } else {
      .data$ViralRefSeq_subject
    })

  # Process and count matching subjects
  result <- result %>%
    group_by(across(all_of(groupby)), .data$Processed_ViralRefSeq_subject) %>%
    summarise(subject_count = n(), .groups = 'drop')
  return(result)



}
