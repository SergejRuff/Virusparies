#' VhgGetSubject: Process and Count Viral Subjects within Groups
#'
#' @param file A data frame containing VirusHunter or VirusGatherer hittable results.
#' @param groupby (optional): A character specifying the column containing the groups (default: "best_query").
#' Note: Gatherer hittables do not have a "best_query" column. Please provide an appropriate column for grouping.
#' @param remove_identifiers (optional): if `TRUE` (default), removes the identifiers in the ViralRefSeq_subject cells.
#' @param include_run_ids (optional): If `TRUE` (default is `TRUE`), adds a fourth column named `run_ids` to the output.
#' This column contains a comma-separated list of unique identifiers from either the `SRA_run` or `run_id` column,
#' aggregated for each combination of group and subject.
#' @param extract_brackets (optional): extract content within square brackets [].
#' @param group_unwanted_phyla (optional): A character string specifying which group of viral phyla to retain in the analysis.
#' Valid values are:
#' \describe{
#'   \item{"rna"}{Retain only the phyla specified for RNA viruses (`valid_phyla_rna`).}
#'   \item{"smalldna"}{Retain only the phyla specified for small DNA viruses (`valid_phyla_smalldna`).}
#'   \item{"largedna"}{Retain only the phyla specified for large DNA viruses (`valid_phyla_largedna`).}
#' }
#' All other phyla not in the specified group will be grouped into a single category:
#' "Non-RNA-virus" for `"rna"`, "Non-Small-DNA-Virus" for `"smalldna"`, or "Non-Large-DNA-Virus" for `"largedna"`.
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
#' @importFrom rlang .data as_string ensym
#' @export
VhgGetSubject <- function(file,
                          groupby = "best_query",
                          remove_identifiers = TRUE,
                          include_run_ids = FALSE,
                          extract_brackets = FALSE,
                          group_unwanted_phyla =NULL){

  #is_file_empty(file)
  if (is_file_empty(file)) {
    #message("Skipping VhgBoxplot generation due to empty data.")
    return(invisible(NULL))  # Return invisible(NULL) to stop further execution
  }

  groupby <- rlang::as_string(rlang::ensym(groupby))

  if (!(groupby %in% c("best_query", "ViralRefSeq_taxonomy"))) {
    stop('Invalid value for groupby. Please use either "best_query" or "ViralRefSeq_taxonomy".')
  }

  arg_logical(remove_identifiers)
  arg_logical(include_run_ids)


  # Determine the column to use for run_id_values
  if ("SRA_run" %in% colnames(file)) {
    run_id_column <- "SRA_run"
  } else if ("run_id" %in% colnames(file)) {
    run_id_column <- "run_id"
  } else {
    stop('Neither "SRA_run" nor "run_id" column found in the file.')
  }


  if(!is.null(group_unwanted_phyla)){

    valid_phyla_rna <-  c("Ambiviricota","Duplornaviricota","Kitrinoviricota",
                          "Lenarviricota","Negarnaviricota ","Pisuviricota","Artverviricota")


    valid_phyla_smalldna <-  c("Ambiviricota","Duplornaviricota","Kitrinoviricota",
                               "Lenarviricota","Negarnaviricota ","Pisuviricota","Artverviricota")

    valid_phyla_largedna <-  c("Ambiviricota","Duplornaviricota","Kitrinoviricota",
                               "Lenarviricota","Negarnaviricota ","Pisuviricota","Artverviricota")

    chosen_list <- switch(group_unwanted_phyla,
                          "rna" = valid_phyla_rna,
                          "smalldna" = valid_phyla_smalldna,
                          "largedna" = valid_phyla_largedna,
                          stop("Invalid group_unwanted_phyla value. Use 'rna', 'smalldna', or 'largedna'."))

    change_label <- switch(group_unwanted_phyla,
                           "rna" = "Non-RNA-virus",
                           "smalldna" = "Non-Small-DNA-Virus",
                           "largedna" = "Non-Large-DNA-Virus")



    file <- VhgAddPhylum(file,"ViralRefSeq_taxonomy")

    # Modify the dataframe
    file <- file %>%
      mutate(
        !!sym(groupby) := case_when(
          Phylum == "unclassified" ~ !!sym(groupby),
          grepl(paste0(chosen_list, collapse = "|"), .data$Phylum, ignore.case = TRUE) ~ !!sym(groupby),
          TRUE ~ change_label
        )
      )

    file <- file %>%
      mutate(
        Phylum = case_when(
          Phylum == "unclassified" ~ "unclassified",  # Keep as unclassified
          grepl(paste0(chosen_list, collapse = "|"), Phylum, ignore.case = TRUE) ~ Phylum,  # Keep original value for valid phyla
          TRUE ~ change_label  # Change to label if not in chosen_list
        )
      )

  }

  file <- file %>%
    mutate(
      Processed_ViralRefSeq_subject = if (remove_identifiers) {
        str_extract(.data$ViralRefSeq_subject, "(?<=\\|).*")
      } else {
        .data$ViralRefSeq_subject
      }
    ) %>%
    mutate(
      Processed_ViralRefSeq_subject = if (extract_brackets) {
        str_extract(.data$Processed_ViralRefSeq_subject, "(?<=\\[)[^\\]]*(?=\\])")
      } else {
        .data$Processed_ViralRefSeq_subject
      }
    )



  # Group by the specified column and the processed viral subjects
  result <- file %>%
    group_by(across(all_of(groupby)), .data$Processed_ViralRefSeq_subject) %>%
    summarise(
      subject_count = n(),
      .groups = 'drop'
    )

  # Optionally include run_id_values
  if (include_run_ids) {
    run_id_values <- file %>%
      group_by(across(all_of(groupby)), .data$Processed_ViralRefSeq_subject) %>%
      summarise(
        run_ids = paste(unique(.data[[run_id_column]]), collapse = ", "),
        .groups = 'drop'
      )
    result <- result %>%
      left_join(run_id_values, by = c(groupby, "Processed_ViralRefSeq_subject"))
  }





  return(result)



}
