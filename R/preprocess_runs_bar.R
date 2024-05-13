#' @title Preprocess virus data for generating a bar plot
#'
#' @description This function preprocesses virus data for generating a bar plot
#' by calculating the number of unique SRA runs for each virus group and
#' creating additional variables for visualization.
#'
#' @param vh_file A data frame containing VirusHunter Hittables results.
#'
#' @return A processed data frame for run_bar plot.
#'
#' @details This function is an internal utility function used within the package.
#' It calculates the number of unique runs for each virus group from the VirusHunters Hittables.
#'
#' @importFrom dplyr group_by summarise mutate n_distinct across any_of
#' @importFrom stringr str_c
#' @importFrom rlang .data
#'
#' @keywords internal
preprocess_runs_bar <- function(vh_file){


  sample_run <- vh_file %>%
    group_by(.data$best_query) %>%
    summarise(unique_SRA_run = n_distinct(across(any_of(c("SRA_run", "run_id"))))) %>%
    mutate(
      perc = round(proportions(.data$unique_SRA_run) * 100, 2),
      res = str_c(.data$unique_SRA_run, " (", .data$perc, "%)"),
      cyl = as.factor(.data$best_query)
    )

  return(sample_run)
}
