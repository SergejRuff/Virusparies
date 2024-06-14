#' @title ImportVirusTable: Import VirusHunterGatherer hittables into R
#'
#' @description This function imports VirusHunterGatherer hittables into R.
#'
#' @param path A character string specifying the path to the VirusHunter or VirusGatherer hittable file.
#'
#' @details This function reads the VirusHunter or VirusGatherer hittable file specified by 'path'
#' into a dataframe in R. The file should be in tab-separated values (TSV) format.
#' The first row should contain column headers.
#'
#' @return A dataframe containing the data from the VirusHunter or VirusGatherer hittable file.
#' @author Sergej Ruff
#' @examples
#' path <- system.file("extdata", "virushunter.tsv", package = "Virusparies")
#' vh_file <- ImportVirusTable(path)
#'
#' print(head(vh_file))
#'
#' # import gatherer files
#' path2 <- system.file("extdata", "virusgatherer.tsv", package = "Virusparies")
#' vg_file <- ImportVirusTable(path2)
#'
#' print(head(vg_file))
#'
#' @seealso
#' VirusHunterGatherer is available here: \url{https://github.com/lauberlab/VirusHunterGatherer}.
#'
#' @importFrom utils read.table
#' @export
ImportVirusTable <- function(path) {
  if (!is.character(path)) {
    stop("Error: 'path' argument must be a character string.")
  }

  # Try to read the file and handle any potential errors
  tryCatch({
    vh_file <- read.table(file = path, sep = "\t", header = TRUE)

    # Print success message
    message("File successfully imported.")

    return(vh_file)
  }, error = function(e) {
    stop("Error: Failed to import file. ", e$message)
  })
}
