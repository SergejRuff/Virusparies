#' @title VhgIdentityScatterPlot: Scatter Plot for Reference Identity vs -log10 of Reference E-value
#'
#' @description
#' This function creates a scatter plot of viral RefSeq identity versus
#' the negative logarithm (base 10) of viral RefSeq E-value. It colors the points
#' based on  best query  and adds a horizontal line representing the cutoff value.
#'
#'
#' @param file A data frame containing VirusHunter or VirusGatherer Hittable results.
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
#' @param cutoff (optional) The significance cutoff value for E-values (default: 1e-5). Removes rows in file
#' with values larger than the cutoff value in the ViralRefSeq_E column.
#' @param conlen_bubble_plot (optional) Logical value indicating whether the `contig_len` column
#'  should be used to size the bubbles in the plot. Applicable only to VirusGatherer hittables input. Default is FALSE.
#' @param contiglen_breaks (optional) number of breaks for the bubble plot (for `conlen_bubble_plot`=TRUE). Default is 5.
#' @param theme_choice (optional) A character indicating the ggplot2 theme to apply. Options include "minimal",
#' "classic", "light", "dark", "void", "grey" (or "gray"), "bw", "linedraw", and "test". Default is "linedraw".
#' @param cut_colour (optional) The color for the horizontal cutoff line. Default is "#990000".
#' @param title (optional) The title of the plot. Default is "Scatterplot of viral reference e-values and identity".
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
#' @param legend_title (optional) The title of the legend. Default is "Phylum".
#' @param legend_position (optional) The position of the legend. Default is "bottom".
#' @param legend_title_size (optional) The size of the legend title text. Default is 12.
#' @param legend_title_face (optional) The face (bold, italic, etc.) of the legend title text. Default is "bold".
#' @param legend_text_size (optional) The size of the legend text. Default is 10.

#'
#' @details
#' VhgIdentityScatterPlot generates a scatter plot for Reference Identity versus -log10 of Reference E-value.
#' It accepts both VirusHunter and VirusGatherer Hittables as input.
#' The plot includes:
#' - A line indicates whether the observed values are above or below the cutoff specified by the 'cutoff' argument (default: 1e-5).
#' - The option `conlen_bubble_plot` = TRUE generates a bubble plot where the size of points corresponds to "contig_len" (exclusive to VirusGatherer).
#'
#' Tibble data frames containing summary statistics (median, Q1, Q3, mean, sd, min, and max) for 'ViralRefSeq_E'
#' and 'ViralRefSeq_ident' values are generated. Optionally, summary statistics for 'contig_len' values
#' are also included if applicable. These summary statistics, along with the plot object, are returned within a list object.
#'
#' @return A list containing the following components:
#' - plot: A plot object representing the faceted scatterplot.
#' - evalue_stats: A tibble data frame with summary statistics for "ViralRefSeq_E" values.
#' - identity_stats: A tibble data frame with summary statistics for "ViralRefSeq_ident" values.
#' - contig_stats (optional): A tibble data frame with summary statistics for "contig_len" values, included only if VirusGatherer is used with `conlen_bubble_plot=TRUE`.
#' @author Sergej Ruff
#' @examples
#' path <- system.file("extdata", "virushunter.tsv", package = "Virusparies")
#' file <- ImportVirusTable(path)
#'
#'
#' # Basic plot
#' plot <- VhgIdentityScatterPlot(file,cutoff = 1e-5)
#'
#' plot(plot$plot)
#'
#'# Custom plot with additional arguments
#' custom_plot <- VhgIdentityScatterPlot(file,
#'                                      cutoff = 1e-5,
#'                                      theme_choice = "dark",
#'                                      cut_colour = "blue",
#'                                      title = "Custom Scatter Plot",
#'                                      title_size = 18,
#'                                      title_face = "italic",
#'                                      title_colour = "darkred",
#'                                      subtitle = "Custom Subtitle",
#'                                      subtitle_size = 14,
#'                                      subtitle_face = "italic",
#'                                      subtitle_colour = "purple",
#'                                      xlabel = "Custom X Label",
#'                                      ylabel = "Custom Y Label",
#'                                      axis_title_size = 14,
#'                                      xtext_size = 12,
#'                                      ytext_size = 12,
#'                                      legend_title = "Custom Legend",
#'                                      legend_position = "top",
#'                                      legend_title_size = 14,
#'                                      legend_title_face = "italic",
#'                                      legend_text_size = 12)
#'
#'plot(custom_plot$plot)
#'
#'# import gatherer files
#' path2 <- system.file("extdata", "virusgatherer.tsv", package = "Virusparies")
#' vg_file <- ImportVirusTable(path2)
#'
#' # vgplot: virusgatherer plot with ViralRefSeq_taxonomy as custom grouping
#' vgplot <- VhgIdentityScatterPlot(vg_file,groupby = "ViralRefSeq_taxonomy")
#' vgplot$plot
#'
#' # plot as bubble plot with contig length as size
#' vgplot_con <- VhgIdentityScatterPlot(vg_file,groupby = "ViralRefSeq_taxonomy",
#' conlen_bubble_plot = TRUE,contiglen_breaks = 4,legend_position = "right")
#'
#' vgplot_con
#'
#' @seealso
#' VirusHunterGatherer is available here: \url{https://github.com/lauberlab/VirusHunterGatherer}.
#'
#' @import ggplot2
#' @importFrom rlang .data
#' @export
VhgIdentityScatterPlot <- function(file,
                                  groupby = "best_query",
                                  taxa_rank = "Family",
                                  cutoff = 1e-5,
                                  conlen_bubble_plot = FALSE,
                                  contiglen_breaks = 5,
                                  theme_choice = "linedraw",
                                  cut_colour = "#990000",
                                  title = "Scatterplot of viral reference e-values and identity",
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
                                  legend_title = "Phylum",
                                  legend_position = "bottom",
                                  legend_title_size = 12,
                                  legend_title_face = "bold",
                                  legend_text_size = 10){

  cutoff <- -log10(cutoff)



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


  # Find the smallest value greater than 0 in ViralRefSeq_E
  min_positive_value <- min(file$ViralRefSeq_E[file$ViralRefSeq_E > 0])

  # Replace all 0 values with the smallest positive value
  file$ViralRefSeq_E[file$ViralRefSeq_E == 0] <- min_positive_value

  # Find the smallest value greater than 0 in ViralRefSeq_E
  min_positive_value <- min(file$ViralRefSeq_E[file$ViralRefSeq_E > 0])

  # Replace all 0 values with the smallest positive value
  file$ViralRefSeq_E[file$ViralRefSeq_E == 0] <- min_positive_value

  # Apply the selected theme
  theme_selected <- select_theme(theme_choice)

  # Calculate the maximum y-value
  max_y <- max(-log10(file$ViralRefSeq_E)) + 5

  # Calculate the minimum y-value
  min_y <- min(-log10(file$ViralRefSeq_E))
  min_y <- ifelse(min_y > 0, 0, min_y)



  color_data <- consistentColourPalette(file, groupby = groupby,taxa_rank=taxa_rank)
  #legend_labels <- color_data$legend_labels
  #labels <- color_data$labels
  matched_vector <- color_data$matched_vector

  # Extract unique values from file$best_query
  #unique_queries <- unique(file[[groupby]])

  # Create a vector of corresponding names from legend_labels
  #pyhlum_names <- legend_labels[unique_queries]

  # Match names to file$best_query
  #file$phylum <- pyhlum_names[match(file[[groupby]], unique_queries)]



  if (conlen_bubble_plot) {
    dynamic_max <- max(file$contig_len)
    dynamic_min <- min(file$contig_len)

    breaks <- round(seq(dynamic_min, dynamic_max, length.out = contiglen_breaks))
    labels <- breaks

    # Reorder file by contig_len in descending order
    file <- file[order(-file$contig_len), ]

    iden_refevalue <- ggplot(file, aes(x = .data$ViralRefSeq_ident, y = -log10(.data$ViralRefSeq_E), size = .data$contig_len))+
      geom_point(aes(color = .data[[groupby]]), alpha = 0.5, shape = 21, fill = "white", colour = "black")+
      geom_point(aes(color=.data[[groupby]]),alpha=0.5)+ geom_hline(aes(yintercept=cutoff), colour=cut_colour)+
      scale_size(name = 'Contig Length', range = c(.1, 15),breaks = breaks, labels = labels)

  } else {
    iden_refevalue <- ggplot(file, aes(x = .data$ViralRefSeq_ident, y = -log10(.data$ViralRefSeq_E)))+
      geom_point(aes(color = .data[[groupby]]), alpha = 0.8, shape = 21, fill = "white", colour = "black")+
      geom_point(aes(color=.data[[groupby]]), alpha = 0.8)+ geom_hline(aes(yintercept=cutoff), colour=cut_colour)
  }



  # Plot the data and color points based on the cutoff condition
  iden_refevalue <- iden_refevalue+
    labs(x= xlabel,
         y= ylabel,
         title=title,
         subtitle = subtitle)+
    theme_selected+
    theme(legend.position = legend_position)+
    guides(color=guide_legend(title=legend_title))+
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



  iden_refevalue  <- iden_refevalue  + scale_color_manual(values = matched_vector)






  # reuse summary stats from boxplot to calculate stats.
  evalue_stats <- boxp_summary_stats(file, group = groupby,ycol ="ViralRefSeq_E")
  identity_stats <- boxp_summary_stats(file, group = groupby,ycol ="ViralRefSeq_ident")




  # Prepare the results list
  results <- list(plot = iden_refevalue, evalue_stats = evalue_stats,
                  identity_stats=identity_stats)

  if (conlen_bubble_plot){

    contig_stats <- boxp_summary_stats(file, group = groupby,ycol ="contig_len")

    results$contig_stats <- contig_stats
  }









  #plot(iden_refevalue)
  message("Scatterplot generation completed.")
  return(results)

}
