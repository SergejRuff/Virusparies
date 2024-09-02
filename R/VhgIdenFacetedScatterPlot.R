#' @title VhgIdenFacetedScatterPlot: Create a scatter plot of Viral refseq identity vs. -log10 of viral refseq E-value.
#'
#' @description
#' VhgIdenFacetedScatterPlot generates a scatter plot of viral refseq identity versus -log10 of refseq E-value
#' for each virus group in the `best_query` or `ViralRefSeq_taxonomy` column . The points are colored based on whether the
#' E-value meets a specified cutoff and are faceted by the viral groups in the `best_query` or `ViralRefSeq_taxonomy` column.
#'
#' @param file VirusHunterGatherer hittable.
#' @param groupby (optional): A character specifying the column containing the groups (default: "best_query").
#' Note: Gatherer hittables do not have a "best_query" column. Please provide an appropriate column for grouping.
#' @param taxa_rank (optional): When `groupby` is set to "ViralRefSeq_taxonomy", specify the taxonomic rank to group your data by.
#' Supported ranks are:
#' - "Subphylum"
#' - "Class"
#' - "Subclass"
#' - "Order"
#' - "Suborder"
#' - "Family" (default)
#' - "Subfamily"
#' - "Genus" (including Subgenus)
#' @param cutoff (optional):  A numeric value representing the cutoff for the refseq E-value. Points with `ViralRefSeq_E`
#' less than or equal to this value will be colored blue; otherwise, they will be colored red (default: 1e-5).
#' @param conlen_bubble_plot (optional):  Logical value indicating whether the `contig_len` column
#'  should be used to size the bubbles in the plot. Applicable only to VirusGatherer hittables input (default: FALSE).
#' @param contiglen_breaks (optional): Number of breaks (default: 5) for the bubble plot (for `conlen_bubble_plot`=TRUE).
#' @param theme_choice (optional): A character indicating the ggplot2 theme to apply. Options include "minimal",
#'  "classic", "light", "dark", "void", "grey" (or "gray"), "bw", "linedraw" (default), and "test".
#'  Append "_dotted" to any theme to add custom dotted grid lines (e.g., "classic_dotted").
#' @param title (optional):  The title of the plot (default: "Faceted scatter plot of viral reference E-values and sequence identity").
#' @param title_size (optional):  The size of the title text (default: 16).
#' @param title_face (optional):  The face (bold, italic, etc.) of the title text (default: "bold").
#' @param title_colour (optional):  The color of the title text (default: "#2a475e").
#' @param subtitle (optional):  The subtitle of the plot (default: NULL).
#' @param subtitle_size (optional):  The size of the subtitle text (default: 12).
#' @param subtitle_face (optional):  The face (bold, italic, etc.) of the subtitle text (default: "bold").
#' @param subtitle_colour (optional):  The color of the subtitle text (default: "#1b2838").
#' @param xlabel (optional):  The label for the x-axis (default: "Viral reference sequence identity (%)").
#' @param ylabel (optional):  The label for the y-axis (default: "-log10 of viral reference E-values").
#' @param axis_title_size (optional):  The size of the axis titles (default: 12).
#' @param xtext_size (optional):  The size of the x-axis text (default: 10).
#' @param x_angle (optional): An integer specifying the angle (in degrees) for the x-axis text labels. Default is NULL, meaning no change.
#' @param ytext_size (optional):  The size of the y-axis text (default: 10).
#' @param y_angle (optional): An integer specifying the angle (in degrees) for the y-axis text labels. Default is NULL, meaning no change.
#' @param legend_position (optional):  The position of the legend (default: "bottom).
#' @param legend_title_size (optional):  The size of the legend title text (default: 12).
#' @param legend_title_face (optional):  The face (bold, italic, etc.) of the legend title text (default: "bold").
#' @param legend_text_size (optional):  The size of the legend text (default: 10).
#' @param true_colour (optional):  The color for points that meet the cutoff condition (default: "blue").
#' @param false_colour (optional):  The color for points that do not meet the cutoff condition (default: "red").
#' @param wrap_ncol (optional):  The number of columns for faceting (default: 12).
#' @param filter_group_criteria (optional):  Character vector, numeric vector, or single character/numeric value.
#'   - Character vector: Names of viral groups to filter.
#'   - Numeric vector: Indices of viral groups to filter.
#'   - Single character or numeric value: Filter a single viral group.
#'   - NULL: No filtering is performed (default).
#'
#' @return A list containing the following components:
#' - Plot: A plot object representing the faceted scatter plot.
#' - evalue_stats: A tibble data frame with summary statistics for "ViralRefSeq_E" values.
#' - identity_stats: A tibble data frame with summary statistics for "ViralRefSeq_ident" values.
#' - contig_stats (optional): A tibble data frame with summary statistics for "contig_len" values, included only if VirusGatherer is used with `conlen_bubble_plot=TRUE`.
#' @author Sergej Ruff
#' @details
#' 'VhgIdenFacetedScatterPlot' takes a VirusHunter or VirusGatherer hittable and a cutoff value as inputs.
#' The plot includes:
#' - Points colored based on whether they meet the cutoff condition.
#' - Faceting by the `best_query` column as the default column. The user can provide their own column
#' for grouping.
#' - The option `conlen_bubble_plot` = TRUE generates a bubble plot where the size of points corresponds to "contig_len" (exclusive to VirusGatherer).
#'
#' `filter_group_criteria`: Allows filtering of viral groups by specifying either a single character
#' string or a vector of character strings that match unique entries in `groupby`.
#' Alternatively, a single numeric value, a range, or a vector of numeric values can be used to filter groups.
#'
#' For example, if `groupby` is "best_query" with the following unique groups:
#' - Anello_ORF1core
#' - Gemini_Rep
#' - Genomo_Rep
#' - Hepadna-Nackedna_TP
#'
#' Setting `filter_group_criteria` to `c("Anello_ORF1core", "Genomo_Rep")` will filter the data to only
#' include observations where the "best_query" column has 'Anello_ORF1core' or 'Genomo_Rep'.
#' Alternatively, setting `filter_group_criteria` to `2:3` will return only the second and third
#' alphabetically ordered viral groups from "best_query".
#' The order also matches the order of the viral groups in the faceted scatter plot.
#'
#' This is particularly useful when there are too many viral groups to be plotted in a single plot,
#' allowing for separation into different groups.
#' It also enables the user to focus on specific groups of interest for more detailed analysis.
#'
#'
#'
#' Tibble data frames containing summary statistics (median, Q1, Q3, mean, sd, min, and max) for 'ViralRefSeq_E'
#' and 'ViralRefSeq_ident' values are generated. Optionally, summary statistics for 'contig_len' values
#' are also included if applicable. These summary statistics, along with the plot object, are returned within a list object.
#'
#' Warning: In some cases, E-values might be exactly 0. When these values are transformed using -log10, R
#' returns "inf" as the output. To avoid this issue, we replace all E-values that are 0 with the smallest E-value that is greater than 0.
#' If the smallest E-value is above the user-defined cutoff, we use a value of `cutoff * 10^-10` to replace the zeros.
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
#' @importFrom rlang .data as_string ensym
#' @export
VhgIdenFacetedScatterPlot <- function(file,
                                     groupby = "best_query",
                                     taxa_rank = "Family",
                                     cutoff = 1e-5,
                                     conlen_bubble_plot = FALSE,
                                     contiglen_breaks = 5,
                                     theme_choice = "linedraw",
                                     title="Faceted scatterplot of viral reference E-values and sequence identity",
                                     title_size = 16,
                                     title_face = "bold",
                                     title_colour = "#2a475e",
                                     subtitle = NULL,
                                     subtitle_size = 12,
                                     subtitle_face = "bold",
                                     subtitle_colour = "#1b2838",
                                     xlabel = "Viral reference sequence identity (%)",
                                     ylabel = "-log10 of viral reference E-values",
                                     axis_title_size = 12,
                                     xtext_size = 10,
                                     x_angle = NULL,
                                     ytext_size = 10,
                                     y_angle = NULL,
                                     legend_position = "bottom",
                                     legend_title_size = 12,
                                     legend_title_face = "bold",
                                     legend_text_size = 10,
                                     true_colour = "blue",
                                     false_colour = "red",
                                     wrap_ncol = 2,
                                     filter_group_criteria = NULL
                                     ){




  # check if hittable is complete
  #is_file_empty(file)
  if (is_file_empty(file)) {
    #message("Skipping VhgBoxplot generation due to empty data.")
    return(invisible(NULL))  # Return invisible(NULL) to stop further execution
  }

  groupby <- rlang::as_string(rlang::ensym(groupby))
  taxa_rank <- rlang::as_string(rlang::ensym(taxa_rank))

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

  # Check if the minimum positive value is below the cutoff
  if (min_positive_value < cutoff) {
    # Replace all 0 values with the smallest positive value
    file$ViralRefSeq_E[file$ViralRefSeq_E == 0] <- min_positive_value
  } else {

    # Add 10^-10 to the cutoff value and use that for the 0 values
    file$ViralRefSeq_E[file$ViralRefSeq_E == 0] <- cutoff * 10^-10
  }


  # Apply the selected theme
  theme_selected <- select_theme(theme_choice)


  max_y <- max(-log10(file$ViralRefSeq_E)) + 5

  # Calculate the minimum y-value
  min_y <- min(-log10(file$ViralRefSeq_E))
  min_y <- ifelse(min_y > 0, 0, min_y)


  file <- filter_specific_group(file,groupby,filter_group_criteria)



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
    )+
    theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))


  iden_refevalue_seperate  <- adjust_plot_angles(iden_refevalue_seperate ,x_angle = x_angle,y_angle = y_angle)

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
  message("Scatter plot generation completed.")
  return(results)


}
