#' @title VgConLenViolin: Generate a Violinplot of contig length for each group (Gatherer only)
#'
#' @description
#' VgConLenViolin creates a violin plot to visualize the distribution of contig lengths
#' for each group in VirusGatherer hittables.
#'
#'
#' @param vg_file A data frame containing VirusGatherer hittable results.
#'
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
#'
#' @param cut (optional): A numeric value representing the cutoff for the refseq E-value (default: 1e-5).
#'
#' @param log10_scale (optinal): transform y-axis to log10 scale (default: TRUE).
#'
#' @param reorder_criteria Character string specifying the criteria for reordering the x-axis ('max', 'min', 'median'(Default),'mean','phylum').
#' NULL sorts alphabetically. You can also specify criteria with 'phylum_' prefix (e.g., 'phylum_median') to sort by phylum first and then by the specified statistic within each phylum.
#'
#' @param adjust_bw (optional): control the bandwidth of the kernel density estimator used to create the violin plot.
#' A higher value results in a smoother plot by increasing the bandwidth, while a lower value can make the plot more detailed but potentially noisier (default: 1).
#'
#' @param jitter_point (optional): logical: TRUE to show all observations, FALSE to show only groups with less than 2 observations (default: FALSE).
#'
#' @param theme_choice (optional): A character indicating the ggplot2 theme to apply. Options include "minimal",
#'  "classic", "light", "dark", "void", "grey" (or "gray"), "bw", "linedraw" (default), and "test".
#'  Append "_dotted" to any theme to add custom dotted grid lines (e.g., "classic_dotted").
#'
#' @param flip_coords (optional): Logical indicating whether to flip the coordinates of the plot (default: TRUE).
#'
#' @param title (optional): The title of the plot (default: "Violinplot of contig length for each group").
#' @param title_size (optional): The size of the title text (default: 16).
#' @param title_face (optional): The face (bold, italic, etc.) of the title text (default: "bold").
#' @param title_colour (optional): The color of the title text (default: "#2a475e").
#'
#' @param subtitle (optional): A character specifying the subtitle of the plot (default: NULL).
#' @param subtitle_size (optional): Numeric specifying the size of the subtitle text(default: 12).
#' @param subtitle_face (optional): A character specifying the font face for the subtitle text (default: "bold").
#' @param subtitle_colour (optional): A character specifying the color for the subtitle text (default: "#1b2838").
#'
#' @param xlabel (optional): The label for the x-axis (default: "Viral group").
#' @param ylabel (optional): The label for the y-axis (default: "Contig length (nt)").
#' @param axis_title_size (optional): The size of the axis titles (default: 12).
#' @param xtext_size (optional): The size of the x-axis text (default: 10).
#' @param x_angle (optional): An integer specifying the angle (in degrees) for the x-axis text labels. Default is NULL, meaning no change.
#' @param ytext_size (optional): The size of the y-axis text (default: 10).
#' @param y_angle (optional): An integer specifying the angle (in degrees) for the y-axis text labels. Default is NULL, meaning no change.
#'
#' @param remove_group_labels (optional): If `TRUE`, the group labels will be removed; if `FALSE` or omitted, the labels will be displayed.
#'
#' @param legend_title (optional): A character specifying the title for the legend (default: "Phylum").
#' @param legend_position (optional): The position of the legend (default: "bottom).
#' @param legend_title_size (optional): Numeric specifying the size of the legend title text (default: 12).
#' @param legend_title_face (optional): A character specifying the font face for the legend title text (default: "bold").
#' @param legend_text_size (optional): Numeric specifying the size of the legend text (default: 10).
#'
#' @param min_observations (optional): Minimum number of observations required per group to be included in the plot (default: 1).
#' @param facet_ncol (optional):  The number of columns for faceting (default: NULL).
#' It is recommended to specify this when the number of viral groups is high, to ensure they fit well in one plot.
#' @param add_boxplot (optional): Add a boxplot to the violin plot (default: FALSE).
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
#' @details
#' VgConLenViolin creates a violin plot to visualize the distribution of contig lengths
#' for each group in the "ViralRefSeq_taxonomy" column of the VirusGatherer hittable.
#' The x-axis represents the groups as defined by the "ViralRefSeq_taxonomy" column,
#' and the y-axis represents the contig lengths.
#'
#' By default, the y-axis is transformed to a log10 scale to better visualize differences
#' in contig lengths across groups. This transformation can be disabled by setting the
#' `log10_scale` argument to FALSE.
#'
#' `min_observations` filters the data sets to include only groups with at least the specified number of observations
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
#' @importFrom rlang .data as_string ensym
#' @export
VgConLenViolin <- function(vg_file=vg_file,
                           taxa_rank = "Family",
                           cut = 1e-5,
                           log10_scale = TRUE,
                           reorder_criteria = "median",
                           adjust_bw = 1,
                           jitter_point=FALSE,
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
                           x_angle = NULL,
                           ytext_size = 10,
                           y_angle = NULL,
                           remove_group_labels = FALSE,
                           legend_title = "Phylum",
                           legend_position = "bottom",
                           legend_title_size = 12,
                           legend_title_face = "bold",
                           legend_text_size = 10,
                           min_observations = 1,
                           facet_ncol = NULL,
                           add_boxplot =FALSE,
                           group_unwanted_phyla =NULL) {

  #check if table is empty
  #is_file_empty(vg_file)
  if (is_file_empty(vg_file)) {
    #message("Skipping VhgBoxplot generation due to empty data.")
    return(invisible(NULL))  # Return invisible(NULL) to stop further execution
  }

  taxa_rank <- rlang::as_string(rlang::ensym(taxa_rank))

  # Define the required columns
  required_columns <- c("contig_len", "ViralRefSeq_E", "ViralRefSeq_taxonomy")



  check_columns(vg_file,required_columns)
  check_input_type(vg_file,c("contig_len", "ViralRefSeq_E"),2)
  check_input_type(vg_file,"ViralRefSeq_taxonomy",1)

  # check arguments
  arg_character(theme_choice)
  arg_character(legend_position)
  arg_logical(flip_coords)
  arg_logical(jitter_point)
  arg_logical(add_boxplot)


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
  vg_file$phyl <- pyhlum_names[match(vg_file[["ViralRefSeq_taxonomy"]], unique_queries)]


  if(!is.null(group_unwanted_phyla)){

    group_unphyla <- remove_non_group(file=vg_file,groupby="ViralRefSeq_taxonomy",chosen_group=group_unwanted_phyla,label_vector=labels,
                                      taxa_rank=taxa_rank)

    vg_file <- group_unphyla$file

    labels <- group_unphyla$label


  }



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

  # Data for adding boxplots (groups with 2 or more observations)
  if (add_boxplot) {
    vg_file_boxplot <- vg_file %>%
      filter(.data$observations >= 2)
  } else {
    vg_file_boxplot <- NULL
  }


  min_y <- 0
  max_y <- max(vg_file$contig_len) + 10
  y_limits <- c(min_y, max_y)



  # # Check for valid reorder_criteria
  valid_criteria <- c("max", "min", "median", "mean", "phylum", "phylum_max", "phylum_min", "phylum_median", "phylum_mean")

  # Check if reorder_criteria is valid
  if (!is.null(reorder_criteria) && !(reorder_criteria %in% valid_criteria)) {
    stop("Invalid reorder_criteria. Please choose one of: max, min, median, mean, phylum, phylum_max, phylum_min, phylum_median, phylum_mean.")
  }







  p <- ggplot(vg_file, aes(x = {
    if (!is.null(reorder_criteria)) {
      if (reorder_criteria == "phylum") {
        ordered_levels <- rev(unique(.data$ViralRefSeq_taxonomy[order(.data$phyl)]))
        factor(.data$ViralRefSeq_taxonomy, levels = ordered_levels)
      }else if (grepl("^phylum_", reorder_criteria)) {
        # Extract the secondary criterion (e.g., "median" from "phylum_median")
        secondary_criteria <- sub("^phylum_", "", reorder_criteria)

        # Use aggregate to calculate the secondary criteria within each phylum
        agg_data <- aggregate(.data$contig_len, list(phylum = .data$phyl, x_val = .data$ViralRefSeq_taxonomy),
                              FUN = switch(secondary_criteria,
                                           "max" = max,
                                           "min" = min,
                                           "median" = median,
                                           "mean" = mean))




        agg_data$phylum <- factor(agg_data$phylum, levels = rev(sort(unique(agg_data$phylum))))


        agg_data <- agg_data[order(agg_data$phylum, agg_data$x), ]


        ordered_levels <- agg_data$x_val[order(agg_data$phylum, agg_data$x)]
        factor(.data$ViralRefSeq_taxonomy, levels = ordered_levels)
      } else {
        reorder(.data$ViralRefSeq_taxonomy, .data$contig_len,
                FUN = switch(reorder_criteria,
                             "max" = max,
                             "min" = min,
                             "median" = median,
                             "mean" = mean))
      }
    } else {
      factor(.data$ViralRefSeq_taxonomy, levels = rev(unique(sort(.data$ViralRefSeq_taxonomy))))
    }
  }, y = .data$contig_len, fill = .data$phyl))+
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

  if (!is.null(vg_file_boxplot) && nrow(vg_file_boxplot) > 0) {
    alpha_val <- ifelse(jitter_point == 1,0.8,0.5)
    p <- p + stat_boxplot(data = vg_file_boxplot,geom = "errorbar", width = 0.2,alpha = alpha_val) +
      geom_boxplot(
        data = vg_file_boxplot,
        width = 0.2,
        alpha = alpha_val,
        outlier.shape = NA,  # Do not draw outliers
        coef = 1.5,  # Control the length of the whiskers (default is 1.5)
        show.legend = FALSE
      )
  }

  # Assuming jitter_point is a logical variable you have defined earlier
  if (jitter_point) {
    p <- p + geom_point(
      position = position_jitter(width = 0.1),
      size = 1.5,
      aes(color = .data$phyl, fill = .data$phyl),  # Use fill for color by phyl
      show.legend = FALSE,
      shape = 21,
      color = "black",  # Outline color
      stroke = 1  # Width of the outline
    )
  } else {
    p <- p + geom_point(
      data = vg_file_filtered,
      position = position_jitter(width = 0.1),
      size = 1.5,
      aes(color = .data$phyl, fill = .data$phyl),
      show.legend = FALSE,
      shape = 21,
      color = "black",
      stroke = 1
    )
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

  p <- adjust_plot_angles(p,x_angle = x_angle,y_angle = y_angle)





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

