#' @title VhgIdenFacetedScatterPlot: Create a Scatterplot of Viral RefSeq Identity vs. -log10 of Viral Reference e-value.
#'
#' @description
#' This function generates a scatterplot of viral RefSeq identity versus -log10 of viral reference e-value
#' for each virus group in the `best_query` column (or another column). The points are colored based on whether the
#' e-value meets a specified cutoff and are faceted by the `best_query` column (or another column).
#'
#' @param vh_file A data frame containing at least the columns `ViralRefSeq_E`, `ViralRefSeq_ident`, and the specified grouping column.
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
#' @param colorblind_support (optional): Logical (TRUE or FALSE). If set to TRUE,
#' the function will use color scales that are more accessible for people with color vision deficiencies.
#' The default value is FALSE.
#' @param colormap (optional) Applies if `colorblind_support = TRUE`. A character string indicating the colormap option to use.
#' Default is "viridis". Eight options are available, derived from the Viridis package:
#'   - "magma" (or "A")
#'   - "inferno" (or "B")
#'   - "plasma" (or "C")
#'   - "viridis" (or "D")
#'   - "cividis" (or "E")
#'   - "rocket" (or "F")
#'   - "mako" (or "G")
#'   - "turbo" (or "H")
#'
#' @return A ggplot object representing the faceted scatterplot.
#' @author Sergej Ruff
#' @details
#' 'VhgIdenFacetedScatterPlot' takes a VirusHunter or VirusGatherer Hittable and a cutoff value as inputs..
#' The plot includes:
#' - Points colored based on whether they meet the cutoff condition.
#' - Faceting by the `best_query` column as the default column. The user can provide their own column
#' for grouping.
#'
#' @examples
#' path <- system.file("extdata", "virushunter.tsv", package = "Virusparies")
#' vh_file <- ImportVirusTable(path)
#'
#' # plot 1
#' plot <- VhgIdenFacetedScatterPlot(vh_file,cutoff = 1e-5)
#'
#' plot
#'
#' # plot 2 with custom data
#' custom_plot <- VhgIdenFacetedScatterPlot(vh_file,
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
VhgIdenFacetedScatterPlot <- function(vh_file,
                                     groupby = "best_query",
                                     taxa_rank = "Family",
                                     cutoff = 1e-5,
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
                                     wrap_ncol = 2,
                                     colorblind_support = FALSE,
                                     colormap = "viridis"
                                     ){




  # check if hittable is complete
  #is_file_empty(vh_file)
  if (is_file_empty(vh_file)) {
    #message("Skipping VhgBoxplot generation due to empty data.")
    return(invisible(NULL))  # Return invisible(NULL) to stop further execution
  }

  if (!(groupby %in% c("best_query", "ViralRefSeq_taxonomy"))) {
    stop('Invalid value for groupby. Please use either "best_query" or "ViralRefSeq_taxonomy".')
  }

  required_columns <- c("ViralRefSeq_E",groupby,"ViralRefSeq_ident")
  check_columns(vh_file,required_columns)
  check_input_type(vh_file,c("ViralRefSeq_E","ViralRefSeq_ident"),2)
  check_input_type(vh_file,groupby,1)

  if(groupby == "ViralRefSeq_taxonomy"){

    vh_file <- VhgPreprocessTaxa(vh_file,taxa_rank)

  }

  # check arguments
  arg_character(theme_choice)
  arg_character(legend_position)
  arg_character(colormap)
  arg_logical(colorblind_support)


  vh_file$cutoff_met <- vh_file$ViralRefSeq_E <= cutoff

  # if (!groupby %in% names(vh_file)) {
  #   stop(paste("Error: Column", groupby, "does not exist in vh_file."))
  # }

  # Find the smallest value greater than 0 in ViralRefSeq_E
  min_positive_value <- min(vh_file$ViralRefSeq_E[vh_file$ViralRefSeq_E > 0])

  # Replace all 0 values with the smallest positive value
  vh_file$ViralRefSeq_E[vh_file$ViralRefSeq_E == 0] <- min_positive_value


  # Apply the selected theme
  theme_selected <- select_theme(theme_choice)


  max_y <- max(-log10(vh_file$ViralRefSeq_E)) + 5

  # Calculate the minimum y-value
  min_y <- min(-log10(vh_file$ViralRefSeq_E))
  min_y <- ifelse(min_y > 0, 0, min_y)


  # Plot the data and color points based on the cutoff condition, faceted by 'best_query'
  iden_refevalue_seperate <- ggplot(vh_file, aes(x = .data$ViralRefSeq_ident, y = -log10(.data$ViralRefSeq_E))) +
    geom_point(aes(color = .data$cutoff_met)) +  # Map color to the cutoff condition
    scale_color_manual(values = c("TRUE" = true_colour, "FALSE" = false_colour)) +  # Define colors for TRUE and FALSE
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

  # add colorblind support
  if(colorblind_support){
    iden_refevalue_seperate<- colorbildsupport(iden_refevalue_seperate,colormap)
  }





  #plot(iden_refevalue_seperate)
  message("Scatterplot generation completed.")
  return(iden_refevalue_seperate)


}
