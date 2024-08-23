#' @title Preprocess virus data for generating a bar plot
#'
#' @description This function preprocesses virus data for generating a bar plot
#' by calculating the number of unique SRA runs for each virus group and
#' creating additional variables for visualization.
#'
#' @param vh_file A data frame containing VirusHunter hittables results.
#' @param groupby string indicating best_query or ViralRefSeq_taxonomy column
#'
#' @return A processed data frame for run_bar plot.
#' @author Sergej Ruff
#' @details This function is an internal utility function used within the package.
#' It calculates the number of unique runs for each virus group from the VirusHunter hittable.
#'
#' @importFrom dplyr group_by summarise mutate n_distinct across any_of coalesce
#' @importFrom stringr str_c
#' @importFrom rlang .data
#'
#' @noRd
preprocess_runs_bar <- function(vh_file,groupby="best_query"){

  all_names <- names(vh_file)
  # Check if either "SRA_run" or "run_id" exists in vh_file
  if (!("SRA_run" %in% all_names) && !("run_id" %in% all_names)) {
    stop("Neither 'SRA_run' nor 'run_id' found in vh_file. Available column names: ", paste(all_names, collapse = ", "))
  }

  # Determine which column to use
  run_column <- if ("SRA_run" %in% all_names) "SRA_run" else "run_id"

  check_columns(vh_file,groupby)

  total_unique_SRA_run <- n_distinct(vh_file[[run_column]])


  sample_run <- vh_file %>%
    group_by(.data[[groupby]]) %>%
    summarise(unique_SRA_run = n_distinct(across(any_of(c("SRA_run", "run_id"))))) %>%
    mutate(
      perc = round(.data$unique_SRA_run / total_unique_SRA_run * 100, 2),
      res = paste(.data$unique_SRA_run, " (", .data$perc, "%)"),
      cyl = as.factor(.data[[groupby]])
    )





  return(sample_run)
}
