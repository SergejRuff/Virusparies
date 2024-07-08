#' @title VgConLenViolin: Generate a Violinplot of contig length for each group (Gatherer only)
#'
#' @description
#' This function creates a violin plot to visualize the distribution of contig lengths
#' for each group in VirusGatherer hittables.
#'
#'
#' @param vg_file A data frame containing VirusGatherer Hittable results.
#' @param taxa_rank (optional): Specify the taxonomic rank to group your data by.
#' Supported ranks are:
#' - "Subphylum"
#' - "Class"
#' - "Subclass"
#' - "Order"
#' - "Suborder"
#' - "Family" (default)
#' - "Subfamily"
#' - "Genus" (including Subgenus)
#' @param cut (optional): The significance cutoff value for E-values (default: 1e-5).
#' @param log10_scale (optinal) transform y-axis to log10 scale (default: TRUE)
#' @param reorder_criteria Character string specifying the criteria for reordering the x-axis ('max', 'min', 'median'(Default),'mean').
#' NULL sorts alphabetically.
#' @param adjust_bw (optional): control the bandwidth of the kernel density estimator used to create the violin plot.
#' A higher value results in a smoother plot by increasing the bandwidth, while a lower value can make the plot more detailed but potentially noisier.
#' Default is 1.
#' @param jitter_point (optional): logical: TRUE to show all observations, FALSE to show only groups with less than 2 observations.
#' @param jitter_point_colour (optional): colour of jitter points. Default is "blue"
#' @param theme_choice (optional): A character indicating the ggplot2 theme to apply. Options include "minimal",
#'  "classic", "light", "dark", "void", "grey" (or "gray"), "bw", "linedraw", and "test".
#'  Default is "linedraw".
#' @param flip_coords (optional): Logical indicating whether to flip the coordinates of the plot. Default is TRUE.
#' @param title (optional): A character specifying the title of the plot. Default title is set based on y_column.
#' @param title_size (optional): Numeric specifying the size of the title text. Default is 16.
#' @param title_face (optional): A character specifying the font face for the title text. Default is "bold".
#' @param title_colour (optional): A character specifying the color for the title text. Default is "#2a475e".
#' @param subtitle (optional): A character specifying the subtitle of the plot. Default subtitle is set based on y_column.
#' @param subtitle_size (optional): Numeric specifying the size of the subtitle text. Default is 12.
#' @param subtitle_face (optional): A character specifying the font face for the subtitle text. Default is "bold".
#' @param subtitle_colour (optional): A character specifying the color for the subtitle text. Default is "#1b2838".
#' @param xlabel (optional): A character specifying the label for the x-axis. Default is "Virus found in query".
#' @param ylabel (optional): A character specifying the label for the y-axis. Default is set based on y_column.
#' @param axis_title_size (optional): Numeric specifying the size of the axis title text. Default is 12.
#' @param xtext_size (optional): Numeric specifying the size of the x-axis tick labels. Default is 10.
#' @param ytext_size (optional): Numeric specifying the size of the y-axis tick labels. Default is 10.
#' @param remove_group_labels (optional): If `TRUE`, the group labels will be removed; if `FALSE` or omitted, the labels will be displayed.
#' @param legend_title (optional): A character specifying the title for the legend. Default is "Phylum".
#' @param legend_position (optional): A character specifying the position of the legend. Default is "bottom".
#' @param legend_title_size (optional): Numeric specifying the size of the legend title text. Default is 12.
#' @param legend_title_face (optional): A character specifying the font face for the legend title text. Default is "bold".
#' @param legend_text_size (optional): Numeric specifying the size of the legend text. Default is 10.
#' @param min_observations (optional): Minimum number of observations required per group to be included in the plot (default: 1).
#' @param facet_ncol (optional):  The number of columns for faceting. Default is NULL.
#' It is recommended to specify this when the number of viral groups is high, to ensure they fit well in one plot.
#'
#'
#' @details
#' This function creates a violin plot to visualize the distribution of contig lengths
#' for each group in the "ViralRefSeq_taxonomy" column of the VirusGatherer hittable.
#' The x-axis represents the groups as defined by the "ViralRefSeq_taxonomy" column,
#' and the y-axis represents the contig lengths.
#'
#' By default, the y-axis is transformed to a log10 scale to better visualize differences
#' in contig lengths across groups. This transformation can be disabled by setting the
#' `log10_scale` argument to FALSE.
#'
#' `min_observations` filters the datasets to include only groups with at least the specified number of observations
#'  before plotting them. This feature allows users to exclude groups with insufficient data.
#'  By default, every group is plotted, as the minimum requirement is set to at least one observation per group.
#'
#'
#' @examples
#'
#' # import gatherer files
#' path2 <- system.file("extdata", "virusgatherer.tsv", package = "Virusparies")
#' vg_file <- ImportVirusTable(path2)
#'
#' # create a violinplot.
#' violinplot <- VgConLenViolin(vg_file=vg_file,cut = 1e-5,log10_scale = TRUE)
#'
#' violinplot
#'
#'
#' @return A list containing the following components:
#' - plot: A plot object representing the violin plot.
#' - contiglen_stats: A tibble data frame with summary statistics for "contig_len" values.
#'
#' @author Sergej Ruff
#'
#' @seealso
#' VirusHunterGatherer is available here: \url{https://github.com/lauberlab/VirusHunterGatherer}.
#'
#' @import ggplot2
#' @importFrom rlang .data
#' @export
VgConLenViolin <- function(vg_file=vg_file,
                           taxa_rank = "Family",
                           cut = 1e-5,
                           log10_scale = TRUE,
                           reorder_criteria = "median",
                           adjust_bw = 1,
                           jitter_point=FALSE,
                           jitter_point_colour = "blue",
                           theme_choice = "linedraw",
                           flip_coords = TRUE,
                           title = "Violinplot of contig length for each group",
                           title_size = 16,
                           title_face = "bold",
                           title_colour = "#2a475e",
                           subtitle = NULL,
                           subtitle_size = 12,
                           subtitle_face = "bold",
                           subtitle_colour = "#1b2838",
                           xlabel = "Viral group",
                           ylabel = "Contig length (nt)",
                           axis_title_size = 12,
                           xtext_size = 10,
                           ytext_size = 10,
                           remove_group_labels = FALSE,
                           legend_title = "Phylum",
                           legend_position = "bottom",
                           legend_title_size = 12,
                           legend_title_face = "bold",
                           legend_text_size = 10,
                           min_observations = 1,
                           facet_ncol = NULL) {

  #check if table is empty
  #is_file_empty(vg_file)
  if (is_file_empty(vg_file)) {
    #message("Skipping VhgBoxplot generation due to empty data.")
    return(invisible(NULL))  # Return invisible(NULL) to stop further execution
  }

  # Define the required columns
  required_columns <- c("contig_len", "ViralRefSeq_E", "ViralRefSeq_taxonomy")



  check_columns(vg_file,required_columns)
  check_input_type(vg_file,c("contig_len", "ViralRefSeq_E"),2)
  check_input_type(vg_file,"ViralRefSeq_taxonomy",1)

  # check arguments
  arg_character(theme_choice)
  arg_character(legend_position)
  arg_logical(flip_coords)


  vg_file <- VhgPreprocessTaxa(vg_file,taxa_rank)

  ## preprocess data for plotting
  vg_file <- vg_file[vg_file$ViralRefSeq_E < cut,]
  #is_file_empty(vg_file)
  if (is_file_empty(vg_file)) {
    #message("Skipping VhgBoxplot generation due to empty data.")
    return(invisible(NULL))  # Return invisible(NULL) to stop further execution
  }

  # Apply the selected theme
  theme_selected <- select_theme(theme_choice)


  # Set the subtitle based on the input
  if (is.null(subtitle) || subtitle == "") {
    subtitle_text <- NULL
  } else {
    subtitle_text <- subtitle
  }


  color_data <- consistentColourPalette(vg_file, groupby = "ViralRefSeq_taxonomy",taxa_rank=taxa_rank)
  legend_labels <- color_data$legend_labels
  labels <- color_data$labels

  # Extract unique values from vg_file$best_query
  unique_queries <- unique(vg_file[["ViralRefSeq_taxonomy"]])

  # Create a vector of corresponding names from legend_labels
  pyhlum_names <- legend_labels[unique_queries]

  # Match names to vg_file$best_query
  vg_file$phylum <- pyhlum_names[match(vg_file[["ViralRefSeq_taxonomy"]], unique_queries)]



  # Filter out groups with fewer than the specified minimum number of observations
  vg_file <- vg_file %>%
    group_by(.data$ViralRefSeq_taxonomy) %>%
    filter(n() >= min_observations) %>%
    ungroup()


  # Calculate number of observations per group
  obs_per_group <- vg_file %>%
    count(.data$ViralRefSeq_taxonomy) %>%
    mutate(observations = n) %>%
    select(.data$ViralRefSeq_taxonomy, .data$observations)

  # Merge original data with observations count
  vg_file <- merge(vg_file, obs_per_group, by = "ViralRefSeq_taxonomy")



  # Filter data for groups with 2 or fewer observations
  vg_file_filtered <- vg_file %>%
    filter(.data$observations < 2)


  min_y <- 0
  max_y <- max(vg_file$contig_len) + 10
  y_limits <- c(min_y, max_y)




  # Check for valid reorder_criteria
  valid_criteria <- c("max", "min", "median", "mean")
  if (!is.null(reorder_criteria) && !(reorder_criteria %in% valid_criteria)) {
    stop("Invalid reorder_criteria. Please choose one of: max, min, median, mean.")
  }


  # Plotting with tryCatch to handle warning
  p <- ggplot(vg_file,aes(x= if (!is.null(reorder_criteria)) {
    reorder(.data$ViralRefSeq_taxonomy, .data$contig_len,
            FUN = switch(reorder_criteria,
                         "max" = max,
                         "min" = min,
                         "median" = median,
                         "mean" = mean))
  } else {
    factor(.data$ViralRefSeq_taxonomy, levels = rev(unique(sort(.data$ViralRefSeq_taxonomy))))
  },
                          y=.data$contig_len,fill=.data$phylum))+
    geom_violin(drop=FALSE,adjust=adjust_bw) +  # Create violin plot
    labs(x=xlabel,
         y=ylabel,
         title=title,
         subtitle = subtitle_text)+
    theme_selected+
    theme(legend.position = legend_position)+
    guides(fill=guide_legend(title=legend_title))+
    theme(

      plot.title = element_text(

        size = title_size,
        face = title_face,
        color = title_colour),
      axis.text.y = element_text(size = ytext_size),
      axis.text.x = element_text(size = xtext_size),
      axis.title = element_text(size = axis_title_size),
      legend.text = element_text(size = legend_text_size),
      legend.key.size = unit(1.5, "lines"),
      legend.title = element_text(size = legend_title_size, face = legend_title_face),
      plot.subtitle = element_text(

        size = subtitle_size,
        face = subtitle_face,
        color= subtitle_colour
      )
    )+
    theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))+
    scale_y_continuous(expand = c(0, 0), limits = y_limits)

  if(jitter_point){
    p <- p + geom_jitter(width = 0.1, size = 1, color = jitter_point_colour)
  }else{
    p <- p + geom_jitter(data = vg_file_filtered, width = 0.1, size = 1, color = jitter_point_colour)

  }


  if(log10_scale){
    p <- p + scale_y_log10()
  }




  if (flip_coords) {
    p <- p + coord_flip()
  }

  p <- p + scale_fill_manual(values = labels)


  p <- remove_group_text(p,remove_group_labels,flip_coords)

  p <- facet_plot(p,facet_ncol,flip_coords)




  # Check if there are groups with 2 or fewer observations
  if (any(vg_file$observations <= 2)) {
    # Print custom warning message
    message("Warning: Some groups have fewer than two observations. These groups will be represented as points in the violin plot.")
  }


  contiglen_stats <- boxp_summary_stats(vg_file, group = "ViralRefSeq_taxonomy",ycol ="contig_len")


  # Prepare the results list
  results <- list(plot = p, contiglen_stats=contiglen_stats)


  #plot(sum_plot)
  message("Violinplot generation completed.")

  return(results)



}

