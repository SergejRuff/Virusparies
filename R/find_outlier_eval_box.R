#' @title find_outlier_eval_box: Detect outliers in a dataset using the interquartile range (IQR) method
#'
#' @description This function detects outliers in a dataset using the interquartile range (IQR) method.
#'
#' @param vh_file  A numeric vector or column of a data frame
#'
#' @return A dataframe indicating outliers in the input data
#'
#' @details This function is an internal utility function used within the package.
#' It identifies outliers in a numeric vector or column of a data frame by comparing the values
#' to the lower and upper bounds determined by the interquartile range (IQR) method.
#'
#' @importFrom dplyr group_by mutate
#' @importFrom rlang .data
#' @keywords internal
find_outlier_eval_box <- function(vh_file){

  ## detect outlier in boxplot

  find_outlier <- function(x) {
    return(x < quantile(x, .25) - 1.5*IQR(x) | x > quantile(x, .75) + 1.5*IQR(x))
  }

  outlier <- vh_file %>%
    group_by(.data$best_query) %>%
    mutate(outlier = ifelse(find_outlier(-log10(.data$ViralRefSeq_E)), -log10(.data$ViralRefSeq_E), NA))

  outlier <- outlier[!is.na(outlier$outlier),]

  return(outlier)
}
