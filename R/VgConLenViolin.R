#' @title VgConLenViolin: Generate a Violinplot of contig length for each group (Gatherer only)
#'
#' @description
#' This function creates a violin plot to visualize the distribution of contig lengths
#' for each group in VirusGatherer hittables.
#'
#'
#' @param vg_file A data frame containing VirusGatherer Hittable results.
#' @param cut (optional) The significance cutoff value for E-values (default: 1e-5).
#' @param log10_scale (optinal) transform y-axis to log10 scale (default: TRUE)
#' @param reorder_criteria Character string specifying the criteria for reordering the x-axis ('max', 'min', 'median'(Default),'mean').
#' @param jitter_point (optional) logical: TRUE to show all observations, FALSE to show only groups with less than 2 observations.
#' @param jitter_point_colour (optional) colour of jitter points. Default is "blue"
#' @param theme_choice (optional): A character indicating the ggplot2 theme to apply. Options include "minimal",
#'  "classic", "light", "dark", "void", "grey" (or "gray"), "bw", "linedraw", and "test".
#'  Default is "linedraw".
#' @param flip_coords (optional) Logical indicating whether to flip the coordinates of the plot. Default is TRUE.
#' @param title (optional) A character specifying the title of the plot. Default title is set based on y_column.
#' @param title_size (optional) Numeric specifying the size of the title text. Default is 16.
#' @param title_face (optional) A character specifying the font face for the title text. Default is "bold".
#' @param title_colour (optional) A character specifying the color for the title text. Default is "#2a475e".
#' @param subtitle (optional) A character specifying the subtitle of the plot. Default subtitle is set based on y_column.
#' @param subtitle_size (optional) Numeric specifying the size of the subtitle text. Default is 12.
#' @param subtitle_face (optional) A character specifying the font face for the subtitle text. Default is "bold".
#' @param subtitle_colour (optional) A character specifying the color for the subtitle text. Default is "#1b2838".
#' @param xlabel (optional) A character specifying the label for the x-axis. Default is "Virus found in query".
#' @param ylabel (optional) A character specifying the label for the y-axis. Default is set based on y_column.
#' @param axis_title_size (optional) Numeric specifying the size of the axis title text. Default is 12.
#' @param xtext_size (optional) Numeric specifying the size of the x-axis tick labels. Default is 10.
#' @param ytext_size (optional) Numeric specifying the size of the y-axis tick labels. Default is 10.
#' @param legend_title (optional) A character specifying the title for the legend. Default is "Phylum".
#' @param legend_position (optional) A character specifying the position of the legend. Default is "bottom".
#' @param legend_title_size (optional) Numeric specifying the size of the legend title text. Default is 12.
#' @param legend_title_face (optional) A character specifying the font face for the legend title text. Default is "bold".
#' @param legend_text_size (optional) Numeric specifying the size of the legend text. Default is 10.
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
#' @author Sergej Ruff
#'
#' @return a plot object
#'
#' @import ggplot2
#' @importFrom rlang .data
#' @export
VgConLenViolin <- function(vg_file=vg_file,
                           cut = 1e-5,
                           log10_scale = TRUE,
                           reorder_criteria = "median",
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
                           legend_title = "Phylum",
                           legend_position = "bottom",
                           legend_title_size = 12,
                           legend_title_face = "bold",
                           legend_text_size = 10,
                           colorblind_support = FALSE,
                           colormap = "viridis") {

  #check if table is empty
  is_file_empty(vg_file)

  # Define the required columns
  required_columns <- c("contig_len", "ViralRefSeq_E", "ViralRefSeq_taxonomy")



  check_columns(vg_file,required_columns)
  check_input_type(vg_file,c("contig_len", "ViralRefSeq_E"),2)
  check_input_type(vg_file,"ViralRefSeq_taxonomy",1)

  # check arguments
  arg_character(theme_choice)
  arg_character(legend_position)
  arg_character(colormap)
  arg_logical(flip_coords)
  arg_logical(colorblind_support)

  vg_file <- taxonomy_group_preprocess(vg_file)

  ## preprocess data for plotting
  vg_file <- vg_file[vg_file$ViralRefSeq_E < cut,]
  is_file_empty(vg_file)

  # Apply the selected theme
  theme_selected <- select_theme(theme_choice)


  # Set the subtitle based on the input
  if (is.null(subtitle) || subtitle == "") {
    subtitle_text <- NULL
  } else {
    subtitle_text <- subtitle
  }


  color_data <- consistentColourPalette(vg_file, groupby = "ViralRefSeq_taxonomy")
  legend_labels <- color_data$legend_labels
  labels <- color_data$labels

  # Extract unique values from vg_file$best_query
  unique_queries <- unique(vg_file[["ViralRefSeq_taxonomy"]])

  # Create a vector of corresponding names from legend_labels
  pyhlum_names <- legend_labels[unique_queries]

  # Match names to vg_file$best_query
  vg_file$phylum <- pyhlum_names[match(vg_file[["ViralRefSeq_taxonomy"]], unique_queries)]







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
    geom_violin(drop=FALSE) +  # Create violin plot
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
    )

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

  # add colorblind support
  if(colorblind_support){
    p<- colorbildsupport(p,colormap)
  }


  # Check if there are groups with 2 or fewer observations
  if (any(vg_file$observations <= 2)) {
    # Print custom warning message
    message("Warning: Some groups have fewer than two observations. These groups will be represented as points in the violin plot.")
  }

  #plot(sum_plot)
  message("Violinplot generation completed.")

  return(p)



}

