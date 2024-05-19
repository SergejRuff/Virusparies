#' Calculate summary statistics and identify hits below a specified E-value threshold for boxplots
#'
#' This function calculates summary statistics, including median, quartiles, minimum, and maximum
#' -log10-transformed E-values, for each virus group from the VirusHunter hittables.
#' Additionally, it identifies hits below a specified E-value threshold and calculates
#' the percentage of hits below the threshold relative to the total hits for each virus group.
#'
#' @param vh_file A data frame containing VirusHunter hittables results.
#' @param cutoff The significance cutoff value for E-values.
#'
#' @return A data frame containing summary statistics and the percentage of hits below the threshold for each virus group
#'
#' @details This function calculates summary statistics such as median, quartiles, minimum,
#' and maximum -log10-transformed E-values for each virus group from the VirusHunter hittables results.
#' It also identifies hits below a specified E-value threshold and calculates the percentage
#' of hits below the threshold relative to the total hits for each virus group.
#'
#' @importFrom dplyr tibble
#' @importFrom stats IQR median quantile reorder aggregate
#' @importFrom rlang .data
#' @keywords internal
vh_sum_stat_evavlue_boxplot <- function(vh_file, cutoff) {
  # Calculate median, quartiles, and other summary statistics
  summary_stats <- aggregate(-log10(vh_file$ViralRefSeq_E) ~ best_query, data = vh_file, FUN = function(x) {
    c(median = median(x), Q1 = quantile(x, 0.25), Q3 = quantile(x, 0.75), min = min(x), max = max(x))
  })

  print(colnames(summary_stats))

  summary_stats <- data.frame(
    best_query = summary_stats$best_query,
    median = summary_stats[, "-log10(vh_file$ViralRefSeq_E)"][, 1],
    Q1 = summary_stats[, "-log10(vh_file$ViralRefSeq_E)"][, 2],
    Q3 = summary_stats[, "-log10(vh_file$ViralRefSeq_E)"][, 3],
    min = summary_stats[, "-log10(vh_file$ViralRefSeq_E)"][, 4],
    max = summary_stats[, "-log10(vh_file$ViralRefSeq_E)"][, 5]
  )

  # Rename columns
  colnames(summary_stats) <- c("best_query", "median", "Q1", "Q3", "min", "max")

  # Identify hits below the cutoff threshold
  below_threshold <- aggregate(num_hits ~ best_query, data = subset(vh_file, .data$ViralRefSeq_E < 10^(-cutoff)), sum)

  # Calculate the total hits for each best_query group
  total_hits <- aggregate(num_hits ~ best_query, data = vh_file, sum)

  # Calculate the percentage of hits below the threshold
  below_threshold <- merge(below_threshold, total_hits, by = "best_query")
  below_threshold$percentage_below_thr <- with(below_threshold, (.data$num_hits.x /.data$num_hits.y) * 100)

  # Rename and keep only the relevant num_hits column
  below_threshold <- within(below_threshold, {
    num_hits <- .data$num_hits.x  # Keeping num_hits.x and renaming it to num_hits
    total_hits <- NULL
    names(percentage_below_thr) <- "percentage_below_cutoff"
  })

  # Join the summarized statistics with the below_threshold information
  summary_stats <- merge(summary_stats, below_threshold, by = "best_query", all.x = TRUE)
  summary_stats[is.na(summary_stats)] <- 0  # Replace NA values with 0

  # Sort by median
  summary_stats <- summary_stats[order(summary_stats$median, decreasing = TRUE), ]

  summary_stats <- summary_stats[,-c(7,8)]

  return(tibble(summary_stats))
}
