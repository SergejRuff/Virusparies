#' @title VhgIdenFacetedScatterPlot: Create a Scatterplot of Viral RefSeq Identity vs. -log10 of Viral Reference e-value.
#'
#' @description
#' This function generates a scatterplot of viral RefSeq identity versus -log10 of viral reference e-value
#' for each virus group in the `best_query` column (or another column). The points are colored based on whether the
#' e-value meets a specified cutoff and are faceted by the `best_query` column (or another column).
#'
#' @param file A data frame containing at least the columns `ViralRefSeq_E`, `ViralRefSeq_ident`, and the specified grouping column.
#' @param groupby (optional) A string indicating the column used for grouping the data points in the scatter plot.
#' The values in this column determine the color of each point in the scatter plot. Default is "best_query".
#' Note: Gatherer hittables do not have a "best_query" column. Please provide an appropriate column for grouping.
#' @param taxa_rank (optional) When `groupby` is set to "ViralRefSeq_taxonomy", specify the taxonomic rank to group your data by.
#' Supported ranks are:
#' - "Subphylum"
#' - "Class"
#' - "Subclass"
#' - "Order"
#' - "Suborder"
#' - "Family" (default)
#' - "Subfamily"
#' - "Genus" (including Subgenus)
#' @param cutoff (optional) A numeric value representing the cutoff for the viral reference e-value. Points with `ViralRefSeq_E`
#' less than or equal to this value will be colored blue; otherwise, they will be colored red (default: 1e-5).
#' @param conlen_bubble_plot (optional) Logical value indicating whether the `contig_len` column
#'  should be used to size the bubbles in the plot. Applicable only to VirusGatherer hittables input. Default is FALSE.
#' @param contiglen_breaks (optional) number of breaks for the bubble plot (for `conlen_bubble_plot`=TRUE). Default is 5.
#' @param theme_choice (optional) A character indicating the ggplot2 theme to apply. Options include "minimal",
#'  "classic", "light", "dark", "void", "grey" (or "gray"), "bw", "linedraw", and "test".
#'  Default is "linedraw".
#' @param title (optional) The title of the plot. Default is "Faceted scatterplot of viral reference e-values and identity".
#' @param title_size (optional) The size of the title text. Default is 16.
#' @param title_face (optional) The face (bold, italic, etc.) of the title text. Default is "bold".
#' @param title_colour (optional) The color of the title text. Default is "#2a475e".
#' @param subtitle (optional) The subtitle of the plot. Default is NULL.
#' @param subtitle_size (optional) The size of the subtitle text. Default is 12.
#' @param subtitle_face (optional) The face (bold, italic, etc.) of the subtitle text. Default is "bold".
#' @param subtitle_colour (optional) The color of the subtitle text. Default is "#1b2838".
#' @param xlabel (optional) The label for the x-axis. Default is "Viral reference identity (%)".
#' @param ylabel (optional) The label for the y-axis. Default is "-log10 of viral reference e-values".
#' @param axis_title_size (optional) The size of the axis titles. Default is 12.
#' @param xtext_size (optional) The size of the x-axis text. Default is 10.
#' @param ytext_size (optional) The size of the y-axis text. Default is 10.
#' @param legend_position (optional) The position of the legend. Default is "bottom".
#' @param legend_title_size (optional) The size of the legend title text. Default is 12.
#' @param legend_title_face (optional) The face (bold, italic, etc.) of the legend title text. Default is "bold".
#' @param legend_text_size (optional) The size of the legend text. Default is 10.
#' @param true_colour (optional) The color for points that meet the cutoff condition. Default is "blue".
#' @param false_colour (optional) The color for points that do not meet the cutoff condition. Default is "red".
#' @param wrap_ncol (optional) The number of columns for faceting. Default is 2.
#'
#' @return A list containing the following components:
#' - plot: A plot object representing the faceted scatterplot.
#' - evalue_stats: A tibble data frame with summary statistics for "ViralRefSeq_E" values.
#' - identity_stats: A tibble data frame with summary statistics for "ViralRefSeq_ident" values.
#' - contig_stats (optional): A tibble data frame with summary statistics for "contig_len" values, included only if VirusGatherer is used with `conlen_bubble_plot=TRUE`.
#' @author Sergej Ruff
#' @details
#' 'VhgIdenFacetedScatterPlot' takes a VirusHunter or VirusGatherer Hittable and a cutoff value as inputs.
#' The plot includes:
#' - Points colored based on whether they meet the cutoff condition.
#' - Faceting by the `best_query` column as the default column. The user can provide their own column
#' for grouping.
#' - The option `conlen_bubble_plot` = TRUE generates a bubble plot where the size of points corresponds to "contig_len" (exclusive to VirusGatherer).
#'
#' Tibble data frames containing summary statistics (median, Q1, Q3, mean, sd, min, and max) for 'ViralRefSeq_E'
#' and 'ViralRefSeq_ident' values are generated. Optionally, summary statistics for 'contig_len' values
#' are also included if applicable. These summary statistics, along with the plot object, are returned within a list object.
#'
#' @examples
#' path <- system.file("extdata", "virushunter.tsv", package = "Virusparies")
#' file <- ImportVirusTable(path)
#'
#' # plot 1
#' plot <- VhgIdenFacetedScatterPlot(file,cutoff = 1e-5)
#'
#' plot
#'
#' # plot 2 with custom data
#' custom_plot <- VhgIdenFacetedScatterPlot(file,
#'                                          cutoff = 1e-4,
#'                                          theme_choice = "dark",
#'                                          title = "Custom Scatterplot",
#'                                          title_size = 18,
#'                                          title_face = "italic",
#'                                          title_colour = "orange",
#'                                          xlabel = "Custom X Label",
#'                                          ylabel = "Custom Y Label",
#'                                          axis_title_size = 14,
#'                                          legend_position = "right",
#'                                          true_colour = "green",
#'                                          false_colour = "purple")
#'
#' custom_plot
#'
#' # import gatherer files
#' path2 <- system.file("extdata", "virusgatherer.tsv", package = "Virusparies")
#' vg_file <- ImportVirusTable(path2)
#'
#' # vgplot: virusgatherer plot with ViralRefSeq_taxonomy as custom grouping
#' vgplot <- VhgIdenFacetedScatterPlot(vg_file,groupby = "ViralRefSeq_taxonomy")
#' vgplot
#'
#'
#' @seealso
#' VirusHunterGatherer is available here: \url{https://github.com/lauberlab/VirusHunterGatherer}.
#'
#'
#' @import ggplot2
#' @importFrom stats runif
#' @importFrom rlang .data
#' @export
VhgIdenFacetedScatterPlot <- function(file,
                                     groupby = "best_query",
                                     taxa_rank = "Family",
                                     cutoff = 1e-5,
                                     conlen_bubble_plot = FALSE,
                                     contiglen_breaks = 5,
                                     theme_choice = "linedraw",
                                     title="Faceted scatterplot of viral reference e-values and identity",
                                     title_size = 16,
                                     title_face = "bold",
                                     title_colour = "#2a475e",
                                     subtitle = NULL,
                                     subtitle_size = 12,
                                     subtitle_face = "bold",
                                     subtitle_colour = "#1b2838",
                                     xlabel = "Viral reference identity (%)",
                                     ylabel = "-log10 of viral reference e-values",
                                     axis_title_size = 12,
                                     xtext_size = 10,
                                     ytext_size = 10,
                                     legend_position = "bottom",
                                     legend_title_size = 12,
                                     legend_title_face = "bold",
                                     legend_text_size = 10,
                                     true_colour = "blue",
                                     false_colour = "red",
                                     wrap_ncol = 2
                                     ){




  # check if hittable is complete
  #is_file_empty(file)
  if (is_file_empty(file)) {
    #message("Skipping VhgBoxplot generation due to empty data.")
    return(invisible(NULL))  # Return invisible(NULL) to stop further execution
  }

  if (!(groupby %in% c("best_query", "ViralRefSeq_taxonomy"))) {
    stop('Invalid value for groupby. Please use either "best_query" or "ViralRefSeq_taxonomy".')
  }

  required_columns <- c("ViralRefSeq_E",groupby,"ViralRefSeq_ident")
  if (conlen_bubble_plot) {
    required_columns <- c(required_columns, "contig_len")
  }
  check_columns(file,required_columns)
  check_input_type(file,c("ViralRefSeq_E","ViralRefSeq_ident"),2)
  check_input_type(file,groupby,1)

  if(groupby == "ViralRefSeq_taxonomy"){

    file <- VhgPreprocessTaxa(file,taxa_rank)

  }

  # check arguments
  arg_character(theme_choice)
  arg_character(legend_position)



  file$cutoff_met <- file$ViralRefSeq_E <= cutoff

  # if (!groupby %in% names(file)) {
  #   stop(paste("Error: Column", groupby, "does not exist in file."))
  # }

  # Find the smallest value greater than 0 in ViralRefSeq_E
  min_positive_value <- min(file$ViralRefSeq_E[file$ViralRefSeq_E > 0])

  # Replace all 0 values with the smallest positive value
  file$ViralRefSeq_E[file$ViralRefSeq_E == 0] <- min_positive_value


  # Apply the selected theme
  theme_selected <- select_theme(theme_choice)


  max_y <- max(-log10(file$ViralRefSeq_E)) + 5

  # Calculate the minimum y-value
  min_y <- min(-log10(file$ViralRefSeq_E))
  min_y <- ifelse(min_y > 0, 0, min_y)



  if (conlen_bubble_plot) {
    dynamic_max <- max(file$contig_len)
    dynamic_min <- min(file$contig_len)

    breaks <- round(seq(dynamic_min, dynamic_max, length.out = contiglen_breaks))
    labels <- breaks

    # Reorder file by contig_len in descending order
    file <- file[order(-file$contig_len), ]

    iden_refevalue_seperate <- ggplot(file, aes(x = .data$ViralRefSeq_ident, y = -log10(.data$ViralRefSeq_E))) +
      geom_point(aes(color = .data$cutoff_met, size = .data$contig_len), alpha = 0.5,stroke = 0.5)+
      scale_size(name = 'Contig Length', range = c(.1, 15),breaks = breaks, labels = labels)

  } else {
    iden_refevalue_seperate <- ggplot(file, aes(x = .data$ViralRefSeq_ident, y = -log10(.data$ViralRefSeq_E))) +
      geom_point(aes(color = .data$cutoff_met), alpha = 0.5,stroke = 0.5)
  }


  # Plot the data and color points based on the cutoff condition, faceted by 'best_query'
  iden_refevalue_seperate <- iden_refevalue_seperate +  # Map color to the cutoff condition
    scale_color_manual(name = "Cutoff Met",values = c("TRUE" = true_colour, "FALSE" = false_colour)) +  # Define colors for TRUE and FALSE
    facet_wrap(~.data[[groupby]], ncol = wrap_ncol)+
    labs(x=xlabel,
         y=ylabel,
         title=title,
         subtitle=subtitle)+
    theme_selected+
    theme(legend.position = legend_position)+
    theme(
      plot.title = element_text(
        size = title_size,
        face = title_face ,
        color = title_colour),
      plot.subtitle = element_text(

        size = subtitle_size,
        face = subtitle_face,
        color= subtitle_colour
      ),
      axis.text.y = element_text(size = ytext_size),
      axis.text.x = element_text(size = xtext_size),
      axis.title = element_text(size = axis_title_size),
      legend.text = element_text(size = legend_text_size),
      legend.key.size = unit(1.5, "lines"),
      legend.title = element_text(size = legend_title_size, face = legend_title_face)
    )+ scale_y_continuous(expand = c(0, 0), limits = c(min_y, max_y))+
    scale_x_continuous(
      expand = c(0, 0),
      limits = c(0, 105),
      breaks = seq(0, 100, by = 10)
    )

  # reuse summary stats from boxplot to calculate stats.
  evalue_stats <- boxp_summary_stats(file, group = groupby,ycol ="ViralRefSeq_E")
  identity_stats <- boxp_summary_stats(file, group = groupby,ycol ="ViralRefSeq_ident")


  # Prepare the results list
  results <- list(plot = iden_refevalue_seperate, evalue_stats = evalue_stats,
                  identity_stats=identity_stats)

  if (conlen_bubble_plot){

    contig_stats <- boxp_summary_stats(file, group = groupby,ycol ="contig_len")

    results$contig_stats <- contig_stats
  }








  #plot(iden_refevalue_seperate)
  message("Scatterplot generation completed.")
  return(results)


}
