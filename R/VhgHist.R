#' @title VhgHist: Histogram and ridge density plots
#'
#' @description
#' VhgHist generates histogram and ridge density plots, enabling users to visualize the distribution and spread of their data.
#'
#' @param file A data frame containing VirusHunter or VirusGatherer hittable results.
#' @param plot_type  Specifies the type of plot to generate:
#'   - `simple_histogram`: A basic histogram with no color differentiation.
#'   - `taxa_colored_histogram`: A histogram where bars are colored based on taxa/groups.
#'   - `ridge_plot`: A ridge plot showing density distributions for different taxa.
#'   If `NULL`, defaults to `simple_histogram`.
#' @param groupby (optional): A character specifying the column containing the groups (default:"best_query").
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
#' @param y_column A character specifying the column containing the values to be compared. Currently "ViralRefSeq_ident",
#' "contig_len" (column in Gatherer hittable) and "ViralRefSeq_E" are supported columns (default:"ViralRefSeq_E").
#' @param contiglen_log10_scale (optional): When `y_column` is set to "contig_len", this parameter enables logarithmic scaling (log10) of the y-axis (TRUE). By default, this feature is disabled (FALSE).
#' @param cut (optional): The significance cutoff value for E-values (default: 10000).
#' @param alpha (optional): Alpha value for elements in plot (default: 1)
#' @param bins (optional): Number of bins. Overridden by binwidth. Defaults to 30.
#' @param binwidth (optional): The width of the bins. Can be specified as a numeric value or as a function that calculates width from unscaled x.
#' @param scale_ridge_plot (optional):  Determines the vertical scaling of the density curves in the ridge plot.
#' @param ridge_bandwidth (optional): Controls the smoothness of the density estimates in the ridge plot.
#'   This parameter specifies the bandwidth of the kernel used for density estimation.
#' @param add_jittered_points (optional): Add jitter points to ridge plot (default: FALSE).
#' @param rel_min_height (optional): sets a percent cutoff relative to the highest point of any of the density curves to remove trailing tails
#' in ridge plot (default: 0)
#' @param color_simple_hist (optional): fill for `taxa_colored_histogram` (default: "gold").
#' @param theme_choice (optional): A character indicating the ggplot2 theme to apply. Options include "minimal",
#'  "classic", "light", "dark", "void", "grey" (or "gray"), "bw", "linedraw" (default), and "test".
#'  Append "_dotted" to any theme to add custom dotted grid lines (e.g., "classic_dotted").
#' @param title (optional): A character specifying the title of the plot.
#' @param title_size (optional): Numeric specifying the size of the title text (default: 16).
#' @param title_face (optional): A character specifying the font face for the title text (default: "bold").
#' @param title_colour (optional): A character specifying the color for the title text (default: "#2a475e").
#' @param subtitle (optional): A character specifying the subtitle of the plot.
#' @param subtitle_size (optional): Numeric specifying the size of the subtitle text(default: 12).
#' @param subtitle_face (optional): A character specifying the font face for the subtitle text (default: "bold").
#' @param subtitle_colour (optional): A character specifying the color for the subtitle text (default: "#1b2838").
#' @param xlabel (optional): A character specifying the label for the x-axis.
#' @param ylabel (optional): A character specifying the label for the y-axis.
#' @param axis_title_size (optional): Numeric specifying the size of the axis title text (default: 12).
#' @param xtext_size (optional): Numeric specifying the size of the x-axis tick labels (default: 10).
#' @param x_angle (optional): An integer specifying the angle (in degrees) for the x-axis text labels. Default is NULL, meaning no change.
#' @param ytext_size (optional): Numeric specifying the size of the y-axis tick labels (default: 10).
#' @param y_angle (optional): An integer specifying the angle (in degrees) for the y-axis text labels. Default is NULL, meaning no change.
#' @param remove_group_labels (optional): If `TRUE`, the group labels will be removed; if `FALSE` or omitted, the labels will be displayed.
#' @param legend_title (optional): A character specifying the title for the legend (default: "Legend").
#' @param legend_position (optional): A character specifying the position of the legend (default: "bottom").
#' @param legend_title_size (optional): Numeric specifying the size of the legend title text (default: 12).
#' @param legend_title_face (optional): A character specifying the font face for the legend title text (default: "bold").
#' @param legend_text_size (optional): Numeric specifying the size of the legend text (default: 10).
#' @param group_unwanted_phyla (optional): A character string specifying which group of viral phyla to retain in the analysis.
#' Valid values are:
#' \describe{
#'   \item{"rna"}{Retain only the phyla specified for RNA viruses (`valid_phyla_rna`).}
#'   \item{"smalldna"}{Retain only the phyla specified for small DNA viruses (`valid_phyla_smalldna`).}
#'   \item{"largedna"}{Retain only the phyla specified for large DNA viruses (`valid_phyla_largedna`).}
#' }
#' All other phyla not in the specified group will be grouped into a single category:
#' "Non-RNA-virus" for `"rna"`, "Non-Small-DNA-Virus" for `"smalldna"`, or "Non-Large-DNA-Virus" for `"largedna"`.
#'
#'
#'
#'
#' @return A plot obj.
#'
#' @author Sergej Ruff
#'
#' @seealso
#' VirusHunterGatherer is available here: \url{https://github.com/lauberlab/VirusHunterGatherer}.
#'
#' @examples
#' path <- system.file("extdata", "virushunter.tsv", package = "Virusparies")
#' file <- ImportVirusTable(path)
#'
#' # simple plot
#' VhgHist(file,plot_type = "simple_histogram",
#' y_column = ViralRefSeq_E,bins = 50,title = "Hist for -log10 of E-values",
#' ylabel = "Frequency",xlabel = "-log10 of viral reference E-values")
#'
#' # plot filled according to  family
#' VhgHist(file,plot_type = "taxa_colored_histogram",y_column = ViralRefSeq_E,
#' bins = 50,title = "Taxa colored Hist for -log10 of E-values",
#' ylabel = "Frequency",xlabel = "-log10 of viral reference E-values",
#' groupby = ViralRefSeq_taxonomy,alpha = 0.7)
#'
#' # ridge plot
#' VhgHist(file,plot_type = "ridge_plot",y_column = ViralRefSeq_E,
#' title = "Ridge plot for -log10 of E-values",ylabel = "Frequency",
#' xlabel = "-log10 of viral reference E-values",groupby = ViralRefSeq_taxonomy,
#' alpha = 0.7,scale_ridge_plot = 0.6)
#'
#' @import ggplot2
#' @importFrom ggridges geom_density_ridges
#' @importFrom rlang .data as_string ensym
#' @export
VhgHist <- function(file,
                    plot_type=NULL,
                    groupby ="best_query",
                    taxa_rank = "Family",
                    y_column = "ViralRefSeq_E",
                    contiglen_log10_scale = FALSE,
                    cut = 10000,
                    alpha = 1,
                    bins = NULL,
                    binwidth = NULL,
                    scale_ridge_plot = 1,
                    ridge_bandwidth = NULL,
                    add_jittered_points = FALSE,
                    rel_min_height = 0,
                    color_simple_hist = "gold",
                    theme_choice = "linedraw",
                    title = NULL,
                    title_size = 16,
                    title_face = "bold",
                    title_colour = "#2a475e",
                    subtitle = NULL,
                    subtitle_size = 12,
                    subtitle_face = "bold",
                    subtitle_colour = "#1b2838",
                    xlabel = NULL,
                    ylabel = NULL,
                    axis_title_size = 12,
                    xtext_size = 10,
                    x_angle = NULL,
                    ytext_size = 10,
                    y_angle = NULL,
                    remove_group_labels = FALSE,
                    legend_title = "Legend",
                    legend_position = "bottom",
                    legend_title_size = 12,
                    legend_title_face = "bold",
                    legend_text_size = 10,
                    group_unwanted_phyla =NULL
){


  #is_file_empty(file)
  if (is_file_empty(file)) {

    return(invisible(NULL))  # Return invisible(NULL) to stop further execution
  }

  # Convert non-quoted column names to strings
  groupby <- rlang::as_string(rlang::ensym(groupby))
  y_column <- rlang::as_string(rlang::ensym(y_column))
  taxa_rank <- rlang::as_string(rlang::ensym(taxa_rank))

  if (!is.null(plot_type)){
    plot_type <- rlang::as_string(rlang::ensym(plot_type))

  }


  if (!(groupby %in% c("best_query", "ViralRefSeq_taxonomy"))) {
    stop('Invalid value for groupby. Please use either "best_query" or "ViralRefSeq_taxonomy".')
  }

  check_columns(file,groupby)
  check_columns(file,y_column)
  check_input_type(file,y_column,2)
  check_input_type(file,groupby,1)

  # check arguments
  arg_character(theme_choice)
  arg_character(legend_position)

  # Apply the selected theme
  theme_selected <- select_theme(theme_choice)

  # set cutoff for contig_len and identity
  if(y_column=="contig_len" | y_column== "ViralRefSeq_ident"){

    message(paste(y_column," column was selected."))
    message(paste("Removing values with E-values higher than: ",cut,"."))

    file <- file[file$ViralRefSeq_E < cut,]

    message(paste0("after removing rows based on evalue the hittable has ",nrow(file)," rows left."))
    #is_file_empty(file)
    if (is_file_empty(file)) {

      return(invisible(NULL))  # Return invisible(NULL) to stop further execution
    }

  }

  if(groupby == "ViralRefSeq_taxonomy"){

    tagroupby <- file$ViralRefSeq_taxonomy

    file <- VhgPreprocessTaxa(file,taxa_rank)

  }


  # change values of E-values to -log10
  if(y_column == "ViralRefSeq_E"){

    y_aes <- -log10(file[[y_column]])
  }else{
    y_aes <- file[[y_column]]
  }


  color_data <- consistentColourPalette(file, groupby = groupby,taxa_rank=taxa_rank)
  matched_vector <- color_data$matched_vector
  legend_labels <- color_data$legend_labels
  labels <- color_data$labels

  # Extract unique values from file$best_query
  unique_queries <- unique(file[[groupby]])

  # Create a vector of corresponding names from legend_labels
  pyhlum_names <- legend_labels[unique_queries]

  # Match names to file$best_query
  file$phyl <- pyhlum_names[match(file[[groupby]], unique_queries)]

  if(!is.null(group_unwanted_phyla)){


    valid_phyla_rna <-  c("Ambiviricota","Duplornaviricota","Kitrinoviricota",
                          "Lenarviricota","Negarnaviricota ","Pisuviricota","Artverviricota")


    valid_phyla_smalldna <-  c("Ambiviricota","Duplornaviricota","Kitrinoviricota",
                               "Lenarviricota","Negarnaviricota ","Pisuviricota","Artverviricota")

    valid_phyla_largedna <-  c("Ambiviricota","Duplornaviricota","Kitrinoviricota",
                               "Lenarviricota","Negarnaviricota ","Pisuviricota","Artverviricota")

    chosen_list <- switch(group_unwanted_phyla,
                          "rna" = valid_phyla_rna,
                          "smalldna" = valid_phyla_smalldna,
                          "largedna" = valid_phyla_largedna,
                          stop("Invalid group_unwanted_phyla value. Use 'rna', 'smalldna', or 'largedna'."))

    change_label <- switch(group_unwanted_phyla,
                           "rna" = "Non-RNA-viruses",
                           "smalldna" = "Non-Small-DNA-Viruses",
                           "largedna" = "Non-Large-DNA-Viruses")

    file <- file %>%
      mutate(
        !!sym(groupby) := case_when(
          phyl == "unclassified" ~ !!sym(groupby),
          grepl(paste0(chosen_list, collapse = "|"), .data$phyl, ignore.case = TRUE) ~ !!sym(groupby),
          TRUE ~ change_label
        )
      )

    # Filter out entries that are not in the valid phyla list, excluding "unclassified"
    matched_vector <- matched_vector[names(matched_vector) %in% c(unique(file[[groupby]]), "unclassified")]

    # Add the "Non-RNA-virus" entry with the color black
    matched_vector[change_label] <- "#000000"

   if(plot_type == "ridge_plot"){

     group_unphyla <- remove_non_group(file=file,groupby=groupby,chosen_group=group_unwanted_phyla,label_vector=labels,
                                       taxa_rank=taxa_rank)

     file <- group_unphyla$file

     labels <- group_unphyla$label

   }



  }



  if (y_column == "ViralRefSeq_ident") {
    max_y <- 105
    min_y <- 0
    y_limits <- c(min_y, max_y)
  } else if (y_column == "ViralRefSeq_E") {
    max_y <- max(y_aes) + 5
    min_y <- min(y_aes)
    min_y <- ifelse(min_y > 0, 0, min_y)
    y_limits <- c(min_y, max_y)
  } else if (y_column == "contig_len") {
    min_y <- 0
    max_y <- max(y_aes) + 10
    y_limits <- c(min_y, max_y)
  }







  ####################################
  # Generate Plot ####################
  ####################################

  p <- file %>% ggplot()

  if (is.null(plot_type) || plot_type == "simple_histogram") {
    # Default case: Simple histogram without color differentiation
    p <- p +
      geom_histogram(mapping = aes(y_aes), bins = bins, alpha = alpha,binwidth=binwidth,
                     fill = color_simple_hist)
  } else if (plot_type == "taxa_colored_histogram") {

    # Taxa-colored histogram: Bars are filled based on .data[[groupby]]
    p <- p +
      geom_histogram(mapping = aes(y_aes, fill = .data[[groupby]]),
                     bins = bins, position = "identity", alpha = alpha,binwidth=binwidth)
  } else if (plot_type == "ridge_plot") {

    # Create the ggplot
    p <- p +
      geom_density_ridges(mapping = aes(x = y_aes, y = .data[[groupby]],fill = .data$phyl),
                          scale = scale_ridge_plot, alpha = alpha, bandwidth = ridge_bandwidth,
                          jittered_points = add_jittered_points, point_alpha = 1, point_shape = 21,
                          rel_min_height = rel_min_height
      )


    # Calculate number of observations per group
    obs_per_group <- file %>%
      count(.data[[groupby]]) %>%
      mutate(observations = n) %>%
      select(.data[[groupby]], .data$observations)

    # Merge original data with observations count
    file <- merge(file, obs_per_group, by = groupby)

    # Check if y_column is "ViralrefSeq_E"
    if (y_column == "ViralRefSeq_E") {
      # Apply -log10 transformation if y_column is "ViralrefSeq_E"
      file$ViralRefSeq_E <- -log10(file$ViralRefSeq_E)
    }



    p <- p + geom_point(
      data = file[file$observations <=2,],
      size = 1.5,
      aes(x = .data[[y_column]], y = .data[[groupby]],color = .data$phyl, fill = .data$phyl),
      #show.legend = FALSE,
      shape = 21,
      color = "black",
      stroke = 1,
      alpha = 0.7
    )


  } else {
    stop("Invalid plot_type. Choose either 'simple_histogram', 'taxa_colored_histogram', or 'ridge_plot'.")
  }

  p <- p + labs(x=xlabel,
                y=ylabel,
                title=title,
                subtitle = subtitle)+
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
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)))

  if(!is.null(plot_type) && plot_type != "ridge_plot"){

    p <- p+
      scale_y_continuous(expand = expansion(mult = c(0, 0.10)), limits = c(0, NA))
  }

  if(y_column == "contig_len" && contiglen_log10_scale){
    p <- p + scale_y_log10()
  }

  if (!is.null(plot_type) && plot_type != "simple_histogram") {

    colors <- if (plot_type == "taxa_colored_histogram") {
      matched_vector
    } else {
      labels
    }

    p <- p + scale_fill_manual(values = colors)
  }

  p <- remove_group_text(p,remove_group_labels,TRUE)



  p <- adjust_plot_angles(p,x_angle = x_angle,y_angle = y_angle)


  if (plot_type == "ridge_plot") {
    message <- "Ridge density plot generation successfully completed."
  } else {
    message <- "Histogram generation successfully completed."
  }


  message(message)
  return(p)
}
