#' @title VhgIdentityScatterPlot: Scatter plot for refseq identity vs -log10 of refseq E-value
#'
#' @description
#' VhgIdentityScatterPlot generates a scatter plot of viral refSeq identity vs. -log10 of viral refseq E-value.
#' It colors the points based on  phylum  and adds a horizontal line representing the cutoff value.
#'
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
#' @param cutoff (optional): A numeric value representing the cutoff for the refseq E-value (default: 1e-5).
#' @param conlen_bubble_plot (optional):  Logical value indicating whether the `contig_len` column
#'  should be used to size the bubbles in the plot. Applicable only to VirusGatherer hittables input (default: FALSE).
#' @param contiglen_breaks (optional): Number of breaks (default: 5) for the bubble plot (for `conlen_bubble_plot`=TRUE).
#' @param theme_choice (optional): A character indicating the ggplot2 theme to apply. Options include "minimal",
#'  "classic", "light", "dark", "void", "grey" (or "gray"), "bw", "linedraw" (default), and "test".
#'  Append "_dotted" to any theme to add custom dotted grid lines (e.g., "classic_dotted").
#' @param cut_colour (optional): The color for the horizontal cutoff line (default: "#990000").
#' @param title (optional): The title of the plot (default: "Scatterplot of viral reference E-values and sequence identity").
#' @param title_size (optional): The size of the title text (default: 16).
#' @param title_face (optional): The face (bold, italic, etc.) of the title text (default: "bold").
#' @param title_colour (optional): The color of the title text (default: "#2a475e").
#' @param subtitle (optional): The subtitle of the plot (default: NULL).
#' @param subtitle_size (optional): The size of the subtitle text (default: 12).
#' @param subtitle_face (optional): The face (bold, italic, etc.) of the subtitle text (default: "bold").
#' @param subtitle_colour (optional): The color of the subtitle text (default: "#1b2838").
#' @param xlabel (optional): The label for the x-axis (default: "Viral reference sequence identity (%)").
#' @param ylabel (optional): The label for the y-axis (default: "-log10 of viral reference E-values").
#' @param axis_title_size (optional): The size of the axis titles (default: 12).
#' @param xtext_size (optional): The size of the x-axis text (default: 10).
#' @param x_angle (optional): An integer specifying the angle (in degrees) for the x-axis text labels. Default is NULL, meaning no change.
#' @param ytext_size (optional): The size of the y-axis text (default: 10).
#' @param y_angle (optional): An integer specifying the angle (in degrees) for the y-axis text labels. Default is NULL, meaning no change.
#' @param legend_title (optional): The title of the legend (default: "Group").
#' @param legend_position (optional): The position of the legend (default: "bottom).
#' @param legend_title_size (optional):  The size of the legend title text (default: 12).
#' @param legend_title_face (optional): The face (bold, italic, etc.) of the legend title text (default: "bold").
#' @param legend_text_size (optional): The size of the legend text (default: 10).
#' @param highlight_groups (optional): A character vector specifying the names of viral groups to be highlighted in the plot (Default:NULL).
#' @param group_unwanted_phyla (optional): A character string specifying which group of viral phyla to retain in the analysis.
#' Valid values are:
#' \describe{
#'   \item{"rna"}{Retain only the phyla specified for RNA viruses.}
#'   \item{"smalldna"}{Retain only the phyla specified for small DNA viruses.}
#'   \item{"largedna"}{Retain only the phyla specified for large DNA viruses.}
#'   \item{"others"}{Retain only the phyla that match small DNA, Large DNA and RNA viruses.}
#' }
#' All other phyla not in the specified group will be grouped into a single category:
#' "Non-RNA-virus" for `"rna"`, "Non-Small-DNA-Virus" for `"smalldna"`,"Non-Large-DNA-Virus" for `"largedna"`,or "Other Viruses" for `"others"`.
#'
#'
#' @details
#' VhgIdentityScatterPlot generates a scatter plot for refseq sequence identity vs -log10 of refseq E-value.
#' It accepts both VirusHunter and VirusGatherer hittables as input.
#' The plot includes:
#' - A line indicates whether the observed values are above or below the cutoff specified by the 'cutoff' argument (default: 1e-5).
#' - The option `conlen_bubble_plot` = TRUE generates a bubble plot where the size of points corresponds to "contig_len" (exclusive to VirusGatherer).
#'
#' Tibble data frames containing summary statistics (median, Q1, Q3, mean, sd, min, and max) for 'ViralRefSeq_E'
#' and 'ViralRefSeq_ident' values are generated. Optionally, summary statistics for 'contig_len' values
#' are also included if applicable. These summary statistics, along with the plot object, are returned within a list object.
#'
#' `highlight_groups` enables the user to specify one or more viral groups from the column indicated in the `groupby` argument. These groups will be highlighted in the plot.
#'
#' Warning: In some cases, E-values might be exactly 0. When these values are transformed using -log10, R
#' returns "inf" as the output. To avoid this issue, we replace all E-values that are 0 with the smallest e-value that is greater than 0.
#' If the smallest E-value is above the user-defined cutoff, we use a value of `cutoff * 10^-10` to replace the zeros.
#'
#' @return A list containing the following components:
#' - Plot: A plot object representing the faceted scatterplot.
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
#' @importFrom rlang .data as_string ensym
#' @importFrom scales hue_pal
#' @export
VhgIdentityScatterPlot <- function(file,
                                  groupby = "best_query",
                                  taxa_rank = "Family",
                                  cutoff = 1e-5,
                                  conlen_bubble_plot = FALSE,
                                  contiglen_breaks = 5,
                                  theme_choice = "linedraw",
                                  cut_colour = "#990000",
                                  title = "Scatterplot of viral reference E-values and sequence identity",
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
                                  legend_title = "Group",
                                  legend_position = "bottom",
                                  legend_title_size = 12,
                                  legend_title_face = "bold",
                                  legend_text_size = 10,
                                  highlight_groups = NULL,
                                  group_unwanted_phyla =NULL){

  cutoff <- -log10(cutoff)



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

  # Calculate the maximum y-value
  max_y <- max(-log10(file$ViralRefSeq_E)) + 5

  # Calculate the minimum y-value
  min_y <- min(-log10(file$ViralRefSeq_E))
  min_y <- ifelse(min_y > 0, 0, min_y)



  color_data <- consistentColourPalette(file, groupby = groupby,taxa_rank=taxa_rank)
  legend_labels <- color_data$legend_labels
  labels_manualcol <- color_data$labels
  # matched_vector <- color_data$matched_vector

  # Extract unique values from file$best_query
  unique_queries <- unique(file[[groupby]])

  # Create a vector of corresponding names from legend_labels
  pyhlum_names <- legend_labels[unique_queries]

  # Match names to file$best_query
  file$phyl <- pyhlum_names[match(file[[groupby]], unique_queries)]


  if(!is.null(group_unwanted_phyla)){

    group_unphyla <- remove_non_group(file=file,groupby=groupby,chosen_group=group_unwanted_phyla,label_vector=labels_manualcol,
                                      taxa_rank=taxa_rank)

    file <- group_unphyla$file

    labels_manualcol <- group_unphyla$label


  }




  if (conlen_bubble_plot) {
    dynamic_max <- max(file$contig_len)
    dynamic_min <- min(file$contig_len)

    breaks <- round(seq(dynamic_min, dynamic_max, length.out = contiglen_breaks))
    labels <- breaks

    # Reorder file by contig_len in descending order
    file <- file[order(-file$contig_len), ]

    iden_refevalue <- ggplot(file, aes(x = .data$ViralRefSeq_ident, y = -log10(.data$ViralRefSeq_E), size = .data$contig_len))+
      geom_point(aes(color = .data$phyl), alpha = 0.5, shape = 21, fill = "white", colour = "black")+
      geom_point(aes(color=.data$phyl),alpha=0.5)+ geom_hline(aes(yintercept=cutoff), colour=cut_colour)+
      scale_size(name = 'Contig Length', range = c(.1, 15),breaks = breaks, labels = labels)

  } else {
    iden_refevalue <- ggplot(file, aes(x = .data$ViralRefSeq_ident, y = -log10(.data$ViralRefSeq_E)))+
      geom_point(aes(color = .data$phyl), alpha = 0.8, shape = 21, fill = "white", colour = "black")+
      geom_point(aes(color=.data$phyl), alpha = 0.8)+ geom_hline(aes(yintercept=cutoff), colour=cut_colour)
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
    )+
    theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

  if (!is.null(highlight_groups)) {

    if (!all(highlight_groups %in% unique(file[[groupby]]))) {
      stop("Error: Some names in 'highlight_groups' are not present in the 'groupby' column.")
    }

    # Generate colors for the highlighted groups
    highlight_colors <- hue_pal()(length(highlight_groups))
    names(highlight_colors) <- highlight_groups

    # Filter the data for highlighting
    highlight_data <- file[file[[groupby]] %in% highlight_groups, ]

    # Create a base plot with all groups and highlight the selected groups
    if (conlen_bubble_plot) {
      iden_refevalue <- iden_refevalue +
        geom_point(data = highlight_data,
                   aes(color = .data[[groupby]], size = .data$contig_len),
                   shape = 21, alpha = 1, stroke = 2)+
        guides(color = guide_legend(override.aes = list(size = 5, stroke = 2)))
    } else {
      iden_refevalue <- iden_refevalue +
        geom_point(data = highlight_data,
                   aes(color = .data[[groupby]]),
                   shape = 21, size = 5, alpha = 1, stroke = 2)
    }

    iden_refevalue <- iden_refevalue +
      scale_color_manual(values = c(labels_manualcol, highlight_colors),
                         breaks = c(names(labels_manualcol), names(highlight_colors))) +
      guides(color = guide_legend(override.aes = list(size = 5, stroke = 2)))
  } else {
    # No highlights; just apply the color scale for base groups
    iden_refevalue <- iden_refevalue +
      scale_color_manual(values = labels_manualcol)
  }


  iden_refevalue  <- adjust_plot_angles(iden_refevalue ,x_angle = x_angle,y_angle = y_angle)









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
