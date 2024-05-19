#' @title importVirusTable: Import VirusHunterGatherer hittables into R
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
#'
#'
#' @importFrom utils read.table
#' @export
importVirusTable <- function(path){

  if (!is.character(path)) {
    stop("Error: 'path' argument must be a character string.")
  }


  vh_file <- read.table(file=path,sep = "\t",header = TRUE)

  return(path)

}
