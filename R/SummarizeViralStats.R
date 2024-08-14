#' @title SummarizeViralStats: Generate summary stats outside of plot functions
#'
#' @description
#' Summarizes data by grouping it according to a specified metric (contig length, E-value or Identity).
#' SummarizeViralStats generates a summary table that includes counts of observations based on a specified metric cutoff.
#' It computes relevant summary statistics depending on the selected metric, with options to filter rows
#' based on a cutoff value.
#'
#'
#' @param file VirusHunterGatherer hittable.
#'
#' @param groupby (optional): A character specifying the column containing the groups (default: "best_query").
#' Note: Gatherer hittables do not have a "best_query" column. Please provide an appropriate column for grouping.
#'
#' @param metric A character string specifying the name of the metric column to be used for calculations.
#'   This column must be present in \code{file}. Supported metric columns include:
#'   \itemize{
#'     \item "contig_len"
#'     \item "ViralRefSeq_E"
#'     \item "ViralRefSeq_ident"
#'   }
#
#'
#' @param metric_cutoff A numeric value used to classify the metric into two categories: below cutoff and above or equal to cutoff.
#'
#' @param filter_cutoff A numeric value for optional filtering of the data based on E-value. Rows where the specified filtering column has a value
#'   less than this cutoff are retained. If \code{NULL}, no filtering is applied. Default is \code{NULL}.
#'
#' @param show_total A logical value indicating whether to include a row with the total sums for each column in the summary table.
#'   Default is \code{FALSE}.
#'
#' @param extra_stats A character vector specifying additional summary statistics to include in the output. Options include:
#'   \itemize{
#'     \item \code{"mean"}
#'     \item \code{"median"}
#'     \item \code{"Q1"}
#'     \item \code{"Q3"}
#'     \item \code{"sd"}
#'     \item \code{"min"}
#'     \item \code{"max"}
#'   }
#'   If \code{NULL} (the default), only the basic counts are included.
#'
#' @return A data frame summarizing the viral stats. The output includes:
#'   \itemize{
#'     \item The count of observations below and above or equal to the \code{metric_cutoff}.
#'     \item Optional additional summary statistics as specified by \code{extra_stats}.
#'     \item An optional total row if \code{show_total} is \code{TRUE}.
#'   }
#'
#' @author Sergej Ruff
#' @seealso
#' VirusHunterGatherer is available here: \url{https://github.com/lauberlab/VirusHunterGatherer}.
#'
#' @examples
#' path <- system.file("extdata", "virushunter.tsv", package = "Virusparies")
#' file <- ImportVirusTable(path)
#'
#' stats <- SummarizeViralStats(file=file,
#' groupby = "best_query",
#' metric = "ViralRefSeq_ident",
#' metric_cutoff = 90,
#' show_total = TRUE,
#' filter_cutoff = 1e-5,
#' extra_stats = c("median","Q1","Q3"))
#'
#' print(stats)
#'
#'
#'
#' @importFrom rlang .data :=
#' @export
SummarizeViralStats <- function(file,
                                groupby = "best_query",
                                metric,
                                metric_cutoff,
                                filter_cutoff = NULL,
                                show_total = FALSE,
                                extra_stats = NULL) {

  #is_file_empty(file)
  if (is_file_empty(file)) {
    #message("Skipping VhgBoxplot generation due to empty data.")
    return(invisible(NULL))  # Return invisible(NULL) to stop further exefilter_cutoffion
  }

  if (!(groupby %in% c("best_query", "ViralRefSeq_taxonomy"))) {
    stop('Invalid value for groupby. Please use either "best_query" or "ViralRefSeq_taxonomy".')
  }

  required_columns <- c("ViralRefSeq_E",groupby,"ViralRefSeq_ident")
  if (metric == "contig_len") {
    required_columns <- c(required_columns, "contig_len")
  }
  check_columns(file,required_columns)

  # Apply the filter if `filter_cutoff` is not NULL
  if (!is.null(filter_cutoff)) {
    file <- file[file$ViralRefSeq_E < filter_cutoff, ]

    message(paste0("after removing rows based on evalue the hittable has ",nrow(file)," rows left."))
    #is_file_empty(file)
    if (is_file_empty(file)) {
      #message("Skipping VhgBoxplot generation due to empty data.")
      return(invisible(NULL))  # Return invisible(NULL) to stop further exefilter_cutoffion
    }

  }

  # Create dynamic variable names
  less_than_label <- paste0("less_than_", metric_cutoff)
  more_or_equal_label <- paste0("equal_or_more_than_",metric_cutoff)



  # Group by .data[[groupby]] and summarize the data
  summary_table <- file %>%
    group_by(.data[[groupby]]) %>%
    summarize(
      !!less_than_label := sum(.data[[metric]] < metric_cutoff),
      !!more_or_equal_label := sum(.data[[metric]] >= metric_cutoff),
      total = n(),
      mean = if ("mean" %in% extra_stats) mean(.data[[metric]]) else NA_real_,
      median = if ("median" %in% extra_stats) median(.data[[metric]]) else NA_real_,
      Q1 = if ("Q1" %in% extra_stats) quantile(.data[[metric]], probs = 0.25) else NA_real_,
      Q3 = if ("Q3" %in% extra_stats) quantile(.data[[metric]], probs = 0.75) else NA_real_,
      sd = if ("sd" %in% extra_stats) sd(.data[[metric]]) else NA_real_,
      min = if ("min" %in% extra_stats) min(.data[[metric]]) else NA_real_,
      max = if ("max" %in% extra_stats) max(.data[[metric]]) else NA_real_
    ) %>%
    ungroup()

  # Remove columns not in `extra_stats` if `extra_stats` is not NULL
  if (!is.null(extra_stats)) {
    # Define the mandatory columns
    mandatory_columns <- c(less_than_label, more_or_equal_label, "total")
    # Create the list of columns to keep
    selected_columns <- c(mandatory_columns, extra_stats)
    # Keep only the selected columns
    summary_table <- summary_table %>%
      select(.data[[groupby]], all_of(selected_columns))
  } else {
    # If extra_stats is NULL, remove all additional stats columns
    summary_table <- summary_table %>%
      select(.data[[groupby]], !!sym(less_than_label), !!sym(more_or_equal_label), data$total)
  }

  # Optionally add a row with the total for each column
  if (show_total) {
    total_row <- summary_table %>%
      summarize(
        across(all_of(groupby), ~ "Total"),
        across(c(!!sym(less_than_label), !!sym(more_or_equal_label), "total"), sum)
      )

    # Bind the total row to the summary table
    summary_table <- bind_rows(summary_table, total_row)
  }

  return(summary_table)
}
