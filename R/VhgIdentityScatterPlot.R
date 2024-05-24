#' @title VhgIdentityScatterPlot: Scatter Plot for Reference Identity vs -log10 of Reference E-value
#'
#' @description
#' This function creates a scatter plot of viral RefSeq identity versus
#' the negative logarithm (base 10) of viral RefSeq E-value. It colors the points
#' based on  best query  and adds a horizontal line representing the cutoff value.
#'
#'
#' @param vh_file A data frame containing VirusHunter Hittable results.
#' @param groupby (optional) A string indicating the column used for grouping the data points in the scatter plot.
#' The values in this column determine the color of each point in the scatter plot. Default is "best_query".
#' Note: Gatherer hittables do not have a "best_query" column. Please provide an appropriate column for grouping.
#' @param cutoff (optional) The significance cutoff value for E-values (default: 1e-5). Removes rows in vh_file
#' with values larger than the cutoff value in the ViralRefSeq_E column.
#' @param theme_choice (optional) A character indicating the ggplot2 theme to apply. Options include "minimal",
#' "classic", "light", "dark", "void", "grey" (or "gray"), "bw", "linedraw", and "test". Default is "minimal".
#' @param cut_colour (optional) The color for the horizontal cutoff line. Default is "#990000".
#' @param title (optional) The title of the plot. Default is "scatterplot for reference identity vs -log10 of reference e-value".
#' @param title_size (optional) The size of the title text. Default is 16.
#' @param title_face (optional) The face (bold, italic, etc.) of the title text. Default is "bold".
#' @param title_colour (optional) The color of the title text. Default is "#2a475e".
#' @param subtitle (optional) The subtitle of the plot. Default is NULL.
#' @param subtitle_size (optional) The size of the subtitle text. Default is 12.
#' @param subtitle_face (optional) The face (bold, italic, etc.) of the subtitle text. Default is "bold".
#' @param subtitle_colour (optional) The color of the subtitle text. Default is "#1b2838".
#' @param xlabel (optional) The label for the x-axis. Default is "viral RefSeq Identity in %".
#' @param ylabel (optional) The label for the y-axis. Default is "-log10 of viral Reference e-value".
#' @param axis_title_size (optional) The size of the axis titles. Default is 12.
#' @param xtext_size (optional) The size of the x-axis text. Default is 10.
#' @param ytext_size (optional) The size of the y-axis text. Default is 10.
#' @param legend_title (optional) The title of the legend. Default is "virus family".
#' @param legend_position (optional) The position of the legend. Default is "bottom".
#' @param legend_title_size (optional) The size of the legend title text. Default is 12.
#' @param legend_title_face (optional) The face (bold, italic, etc.) of the legend title text. Default is "bold".
#' @param legend_text_size (optional) The size of the legend text. Default is 10.
#'
#' @return A ggplot object representing the scatter plot.
#'
#' @examples
#' path <- system.file("extdata", "virushunter.tsv", package = "Virusparies")
#' vh_file <- importVirusTable(path)
#'
#'
#' # Basic plot
#' plot <- VhgIdentityScatterPlot(vh_file,cutoff = 1e-5)
#'
#' plot(plot)
#'
#'# Custom plot with additional arguments
#' custom_plot <- VhgIdentityScatterPlot(vh_file,
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
#'plot(custom_plot)
#'
#' @import ggplot2
#' @importFrom rlang .data
#' @export
VhgIdentityScatterPlot <- function(vh_file,
                                  groupby = "best_query",
                                  cutoff = 1e-5,
                                  theme_choice = "minimal",
                                  cut_colour = "#990000",
                                  title = "scatterplot for refrence identity vs -log10 of refrence e-value",
                                  title_size = 16,
                                  title_face = "bold",
                                  title_colour = "#2a475e",
                                  subtitle = NULL,
                                  subtitle_size = 12,
                                  subtitle_face = "bold",
                                  subtitle_colour = "#1b2838",
                                  xlabel = "viral Refseq Identity in %",
                                  ylabel = "-log10 of viral Reference e-value",
                                  axis_title_size = 12,
                                  xtext_size = 10,
                                  ytext_size = 10,
                                  legend_title = "virus family",
                                  legend_position = "bottom",
                                  legend_title_size = 12,
                                  legend_title_face = "bold",
                                  legend_text_size = 10){

  cutoff <- -log10(cutoff)

  if (!groupby %in% names(vh_file)) {
    stop(paste("Error: Column", groupby, "does not exist in vh_file."))
  }

  # Apply the selected theme
  theme_selected <- select_theme(theme_choice)

  # Calculate the maximum y-value
  max_y <- max(-log10(vh_file$ViralRefSeq_E)) + 5

  # Plot the data and color points based on the cutoff condition
  iden_refevalue <- ggplot(vh_file, aes(x = .data$ViralRefSeq_ident, y = -log10(.data$ViralRefSeq_E))) +
    geom_point(aes(color=.data[[groupby]]))+ geom_hline(aes(yintercept=cutoff), colour=cut_colour)+
    labs(x= xlabel,
         y= ylabel,
         title=title,
         subtitle = subtitle)+
    theme_selected+
    theme(legend.position = legend_position)+
    guides(fill=guide_legend(title=legend_title))+
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
    )+ scale_y_continuous(expand = c(0, 0), limits = c(0, max_y))+
    scale_x_continuous(
      expand = c(0, 0),
      limits = c(0, 105),
      breaks = seq(0, 100, by = 10)
    )

  #plot(iden_refevalue)

  return(iden_refevalue)
}
