#' @title VhgRunsBarplot: Generate a bar plot showing the number of data sets with unique runs found for
#' each Virus group.
#'
#' @param file Data frame containing Virushunter or VirusGatherer hittables results.
#' @param groupby (optional) A string indicating the column used for grouping the data points in the plot.
#' "best_query" and "ViralRefSeq_taxonomy" can be used. Default is "best_query".
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
#' @param cut (optional) The significance cutoff value for E-values (default: 1e-5).
#' Removes rows in file with values larger than cutoff value in ViralRefSeq_E column.
#' @param reorder_criteria Character string specifying the criteria for reordering the x-axis ('max' (default), 'min').
#' NULL sorts alphabetically.
#' @param theme_choice (optional) A character indicating the ggplot2 theme to apply.
#' Options include "minimal", "classic", "light", "dark", "void", "grey" (or "gray"), "bw", "linedraw", and "test". Default is "linedraw".
#' @param flip_coords (optional) Logical indicating whether to flip the coordinates of the plot. Default is TRUE.
#' @param title (optional) A character specifying the title of the plot. Default is "Distribution of viral groups detected across query sequences".
#' @param title_size (optional) The size of the plot title. Default is 16.
#' @param title_face (optional) The font face of the plot title. Default is "bold".
#' @param title_colour (optional) The color of the plot title. Default is "#2a475e".
#' @param subtitle (optional) A character specifying the subtitle of the plot.
#' Default is "default", which calculates the total number of data sets with hits and returns it as
#' "total number of data sets with hits: " followed by the calculated number. an empty string ("") removes the subtitle.
#' @param subtitle_size (optional) The size of the plot subtitle. Default is 12.
#' @param subtitle_face (optional) The font face of the plot subtitle. Default is "bold".
#' @param subtitle_colour (optional) The color of the plot subtitle. Default is "#1b2838".
#' @param xlabel (optional) A character specifying the label for the x-axis. Default is "Viral group".
#' @param ylabel (optional) A character specifying the label for the y-axis. Default is "Number of data sets with hits for group".
#' @param axis_title_size (optional) The size of axis titles. Default is 12.
#' @param xtext_size (optional) The size of x-axis text labels. Default is 10.
#' @param ytext_size (optional) The size of y-axis text labels. Default is 10.
#' @param remove_group_labels (optional): If `TRUE`, the group labels will be removed; if `FALSE` or omitted, the labels will be displayed.
#' @param legend_title (optional) A character specifying the title of the legend. Default is "Phylum".
#' @param legend_position (optional) The position of the legend. Default is "bottom"."none" removes legend.
#' @param legend_title_size (optional) The size of the legend title. Default is 12.
#' @param legend_title_face (optional) The font face of the legend title. Default is "bold".
#' @param legend_text_size (optional) The size of the legend text. Default is 10.
#' @param plot_text An index (0-3) to select the variable for text labels.
#' - 0 = None.
#' - 1 = Number of viral groups detected across query sequences.
#' - 2 = Only the percentage.
#' - 3 = Both (Default).
#' @param plot_text_size (optional) The size of the text labels added to the plot. Default is 3.5.
#' @param plot_text_position_dodge (optional) The degree of dodging for positioning text labels. Default is 0.9.
#' @param plot_text_hjust (optional) The horizontal justification of text labels. Default is -0.1.
#' @param plot_text_vjust (optional) The vertical justification of text labels. Default is 0.5.
#' It is recommended to change `vjust` when setting `flip_coords = FALSE`.
#' @param plot_text_colour (optional) The color of the text labels added to the plot. Default is "black".
#'
#'
#' @details
#' 'VhgRunsBarplot' generates a bar plot showing the number of data sets with unique runs found for
#' each Virus group. It takes VirusHunter and VirusGatherer hittables as Input.
#'
#' Only significant values below the threshold specified by the 'cut' argument (default: 1e-5) are included in the plot.
#'
#'
#'
#' @return A list containing the bar plot and tabular data with information from the plot.
#' @author Sergej Ruff
#' @examples
#' path <- system.file("extdata", "virushunter.tsv", package = "Virusparies")
#' file <- ImportVirusTable(path)
#'
#' # plot 1: plot boxplot for "identity"
#' plot <- VhgRunsBarplot(file,cut = 1e-5)
#' plot
#'
#' # return sample_run inside plot object
#' print(plot$sample_run)
#'
#' # Plot 2: Customized plot with modified settings
#' plot_custom <- VhgRunsBarplot(
#'   file,
#'   cut = 1e-6, # Lower cutoff value
#'   theme_choice = "grey", # Classic theme
#'   title = "Customized Bar Plot", # Custom title
#'   subtitle = "test_subtitle", # No subtitle
#'   xlabel = "Custom X Label", # Custom x-axis label
#'   ylabel = "Custom Y Label", # Custom y-axis label
#'   legend_position = "top", # Legend position on top
#'   plot_text_size = 5, # Larger text label size
#'   plot_text_colour = "red" # Red text labels
#' )
#' plot_custom
#'
#'
#' # import gatherer files
#' path2 <- system.file("extdata", "virusgatherer.tsv", package = "Virusparies")
#' vg_file <- ImportVirusTable(path2)
#'
#'
#  plot 3: Virusgatherer plot
#  plot3 <- VhgRunsBarplot(vg_file,groupby = "ViralRefSeq_taxonomy",cut = 1e-5)
#  plot3
#'
#' @seealso
#' VirusHunterGatherer is available here: \url{https://github.com/lauberlab/VirusHunterGatherer}.
#'
#' @import ggplot2
#' @importFrom  dplyr n_distinct coalesce
#' @importFrom rlang .data
#' @export
VhgRunsBarplot <- function(file,
                          groupby = "best_query",
                          taxa_rank = "Family",
                          cut = 1e-5,
                          reorder_criteria = "max",
                          theme_choice = "linedraw",
                          flip_coords = TRUE,
                          title = "Distribution of viral groups detected across query sequences",
                          title_size = 16,
                          title_face = "bold",
                          title_colour = "#2a475e",
                          subtitle = "default",
                          subtitle_size = 12,
                          subtitle_face = "bold",
                          subtitle_colour = "#1b2838",
                          xlabel = "Viral group",
                          ylabel = "Number of datasets with hits for group",
                          axis_title_size = 12,
                          xtext_size = 10,
                          ytext_size = 10,
                          remove_group_labels = FALSE,
                          legend_title = "Phylum",
                          legend_position = "bottom",
                          legend_title_size = 12,
                          legend_title_face = "bold",
                          legend_text_size = 10,
                          plot_text = 3,
                          plot_text_size = 3.5,
                          plot_text_position_dodge = 0.9,
                          plot_text_hjust = -0.1,
                          plot_text_vjust = 0.5,
                          plot_text_colour = "black"
                          ){




  # check if hittable is empty
  #is_file_empty(file)
  if (is_file_empty(file)) {
    #message("Skipping VhgBoxplot generation due to empty data.")
    return(invisible(NULL))  # Return invisible(NULL) to stop further execution
  }

  if (!(groupby %in% c("best_query", "ViralRefSeq_taxonomy"))) {
    stop('Invalid value for groupby. Please use either "best_query" or "ViralRefSeq_taxonomy".')
  }

  # check arguments
  arg_character(theme_choice)
  arg_character(legend_position)

  arg_logical(flip_coords)





  required_columns <- c("ViralRefSeq_E",groupby)
  all_names <- names(file)


  check_columns(file,required_columns)
  # Check if either "SRA_run" or "run_id" exists in file
  if (!("SRA_run" %in% all_names) && !("run_id" %in% all_names)) {
    stop("Neither 'SRA_run' nor 'run_id' found in file. Available column names: ", paste(all_names, collapse = ", "))
  }

  check_input_type(file,"ViralRefSeq_E",2)
  check_input_type(file,groupby,1)

  if(groupby == "ViralRefSeq_taxonomy"){

    file <- VhgPreprocessTaxa(file,taxa_rank)

  }

  # Filter obj
  file <- file[file$ViralRefSeq_E < cut,]
  message(paste0("after removing rows based on evalue the hittable has ",nrow(file)," rows left."))
  #check if obj has 0 ob after filtering
  #is_file_empty(file)
  if (is_file_empty(file)) {
    #message("Skipping VhgBoxplot generation due to empty data.")
    return(invisible(NULL))  # Return invisible(NULL) to stop further execution
  }

  # Apply the selected theme
  theme_selected <- select_theme(theme_choice)


  # preprocess data for plot
  sample_run <- preprocess_runs_bar(file,groupby)

  # Set the subtitle based on the input
  if (subtitle == "default") {
    subtitle_text <- paste0("total number of datasets: ", n_distinct(coalesce(file$SRA_run, file$run_id)))
  } else {
    subtitle_text <- ifelse(is.null(subtitle) || subtitle == "", NULL, subtitle)
  }

  color_data <- consistentColourPalette(file, groupby = groupby,taxa_rank=taxa_rank)
  legend_labels <- color_data$legend_labels
  labels <- color_data$labels

  # Extract unique values from file$best_query
  unique_queries <- unique(sample_run[[groupby]])

  # Create a vector of corresponding names from legend_labels
  pyhlum_names <- legend_labels[unique_queries]

  # Match names to sample_run$best_query
  sample_run$phylum <- pyhlum_names[match(sample_run[[groupby]], unique_queries)]

  # Check for valid reorder_criteria
  valid_criteria <- c("max", "min")
  if (!is.null(reorder_criteria) && !(reorder_criteria %in% valid_criteria)) {
    stop("Invalid reorder_criteria. Please choose one of: max, min.")
  }

  text_var_names <- c("none", "unique_SRA_run", "perc", "res")

  # Check if the provided index is valid
  if (plot_text < 0 || plot_text >= length(text_var_names)) {
    stop("Invalid plot_text. Choose an index from 0 to 3.")
  }

  text_var <- text_var_names[plot_text + 1]

  sample_run <- sample_run %>% mutate(perc=paste(.data$perc,"%"))




  #####################
  ### Generate Plot ###
  ### #################



  run_bar <- ggplot(data = sample_run,aes(x=if (!is.null(reorder_criteria)) {
    reorder(.data[[groupby]], if (reorder_criteria == "max") .data$unique_SRA_run else -.data$unique_SRA_run)
  } else {
    factor(.data[[groupby]], levels = rev(unique(sort(.data[[groupby]]))))
  },y=.data$unique_SRA_run,fill= .data$phylum))+
    geom_bar(stat = "identity")+
    labs(title = title,
         x= xlabel,
         y= ylabel,
         subtitle = subtitle_text)+
    theme_selected+
    theme(legend.position = legend_position,
          axis.text.y = element_text(size = ytext_size),
          axis.text.x = element_text(size = xtext_size),
          axis.title = element_text(size = axis_title_size ),
          legend.text = element_text(size = legend_text_size),
          legend.key.size = unit(1.5, "lines"),
          legend.title = element_text(size = legend_title_size, face = legend_title_face))+
    guides(fill=guide_legend(title=legend_title))+
    scale_y_continuous(expand = c(0, 0), limits = c(0, max(sample_run$unique_SRA_run)+max(sample_run$unique_SRA_run)/10))+
    theme(
      # This is the new default font in the plot
      plot.title = element_text(

        size = title_size,
        face = title_face,
        color = title_colour),
      # Statistical annotations below the main title
      plot.subtitle = element_text(

        size = subtitle_size,
        face = subtitle_face,
        color= subtitle_colour
      ))+
    theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))




  if (text_var != "none") {


    run_bar <- run_bar + geom_text(aes(label = .data[[text_var]]),
                       hjust = plot_text_hjust,
                       vjust = plot_text_vjust,
                       color = plot_text_colour,
                       position = position_dodge(plot_text_position_dodge),
                       size = plot_text_size)
  }



  if (flip_coords) {
    run_bar <- run_bar + coord_flip()
  }

  run_bar <- remove_group_text(run_bar,remove_group_labels,flip_coords)

  if(groupby != "SRA_run"){


    run_bar  <- run_bar  + scale_fill_manual(values = labels)


  }










  message("Bar chart generation completed.")
  return(list(plot=run_bar,
              sample_run=sample_run))




}





