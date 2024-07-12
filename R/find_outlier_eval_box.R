#' @title find_outlier_eval_box: Detect outliers in a dataset using the interquartile range (IQR) method
#'
#' @description This function detects outliers in a dataset using the interquartile range (IQR) method.
#'
#' @param vh_file  A numeric vector or column of a data frame
#' @param group column to group stats by
#' @param y_column column with metric values.
#'
#' @return A dataframe indicating outliers in the input data
#'
#' @details This function is an internal utility function used within the package.
#' It identifies outliers in a numeric vector or column of a data frame by comparing the values
#' to the lower and upper bounds determined by the interquartile range (IQR) method.
#' @author Sergej Ruff
#' @import  dplyr
#' @importFrom stats na.omit IQR
#' @noRd
find_outlier_eval_box <- function(vh_file, group = "best_query", y_column = "ViralRefSeq_E") {

  ## detect outlier in boxplot
  find_outlier <- function(x) {
    return(x < quantile(x, .25) - 1.5 * IQR(x) | x > quantile(x, .75) + 1.5 * IQR(x))
  }

  # Apply transformation only if y_column is "ViralRefSeq_E"
  if(y_column == "ViralRefSeq_E") {
    vh_file$transformed_value <- -log10(vh_file[[y_column]])
  } else {
    vh_file$transformed_value <- vh_file[[y_column]]
  }


  # Identify outliers within each group
  outliers <- lapply(split(vh_file, vh_file[[group]]), function(subgroup) {
    subgroup$outlier <- ifelse(find_outlier(subgroup$transformed_value), subgroup$transformed_value, NA)
    na.omit(subgroup)
  })

  # Combine outliers from all groups
  result <- do.call(rbind, outliers)

  # Remove the transformed_value column
  result$transformed_value <- NULL

  # Reset row names
  rownames(result) <- NULL


  return(as_tibble(result))
}
