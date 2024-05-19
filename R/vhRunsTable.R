#' @title Generate a gt summary table of unique SRA runs for each virus group
#'
#' @description This function generates a summary table of unique SRA runs for each virus group
#' based on the input dataset.
#'
#' @param vh_file A data frame containing the Virushunter Hittables Results
#' @param cut The significance cutoff value for E-values (default: 1e-5). Removes rows in vh_file
#' with values larger than cutoff value in ViralRefSeq_E column.
#'
#' @return A formatted gt table summarizing unique runs for each virus group
#'
#' @details This function is an internal utility function used within the package.
#' It calculates the number of unique runs for each virus group from the input dataset.
#'
#' @seealso \code{\link{vhRunsBarplot}}
#'
#' @importFrom dplyr group_by summarise n_distinct across any_of
#' @importFrom gt gt
#' @importFrom rlang .data
#' @export
vhRunsTable <- function(vh_file,cut = 1e-5){

  vh_file <- vh_file[vh_file$ViralRefSeq_E < cut,]

  which_runs <- vh_file %>%
    group_by(.data$best_query) %>%
    summarise(unique_SRA_run = n_distinct(across(any_of(c("SRA_run", "run_id")))),
              SRAs_found = paste(unique(across(any_of(c("SRA_run", "run_id")))), collapse = ", "))

  # Remove "c()", "\", and quotation marks
  which_runs$SRAs_found <- gsub("[c()\\\\\"]", "", which_runs$SRAs_found)

  # creat a table
  which_runs_table <- which_runs %>% gt()

  return(which_runs_table)

}
