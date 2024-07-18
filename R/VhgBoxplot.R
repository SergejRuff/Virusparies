#' @title VhgBoxplot: Generate box plots comparing E-values,identity or contig length (Gatherer only) for each virus group
#'
#' @description
#'  This function generates box plots comparing either E-values,identity or contig length (Gatherer only)
#' for each group from VirusHunter or VirusGatherer hittables results.
#'
#' @param file A data frame containing VirusHunter or VirusGatherer Hittable results.
#' @param x_column A character specifying the column containing the groups (Default:"best_query").
#' @param taxa_rank (optional): When `x_column` is set to "ViralRefSeq_taxonomy", specify the taxonomic rank to group your data by.
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
#' "contig_len" (Column in Gatherer hittable) and "ViralRefSeq_E" are supported columns (Default:"ViralRefSeq_E").
#' @param contiglen_log10_scale (optional): When `y_column` is set to "contig_len", this parameter enables logarithmic scaling (log10) of the y-axis (TRUE). By default, this feature is disabled (FALSE).
#' @param cut (optional): The significance cutoff value for E-values (default: 1e-5).
#' @param add_cutoff_line (optional): whether to add a horizontal line based on `cut` for `"ViralRefSeq_E"` column (default: TRUE).
#' @param cut_colour (optional): The color for the significance cutoff line (default: "#990000").
#' @param reorder_criteria Character string specifying the criteria for reordering the x-axis ('max', 'min', 'median'(Default),'mean').
#' NULL sorts alphabetically.
#' @param theme_choice (optional):: A character indicating the ggplot2 theme to apply. Options include "minimal",
#'  "classic", "light", "dark", "void", "grey" (or "gray"), "bw", "linedraw", and "test".
#'  Default is "linedraw".
#' @param flip_coords (optional): Logical indicating whether to flip the coordinates of the plot. Default is TRUE.
#' @param add_mean_point (optional): Logical indicating whether to add mean points to the boxplot. Default is FALSE.
#' @param mean_color (optional): change color of point indicating mean value in box plot (default: "white").
#' @param mean_point_size (optional): change size of point indicating mean value in box plot (default: 2)
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
#' @param facet_ncol (optional):  The number of columns for faceting. Default is NULL.
#' It is recommended to specify this when the number of viral groups is high, to ensure they fit well in one plot.
#'
#'
#' @details
#' VhgBoxplot generates box plots comparing either E-values, identity, or contig length (Gatherer only) for each virus group from the VirusHunter or Gatherer Hittable.
#'
#' The user can specify whether to generate boxplots for E-values, identity, or contig length (Gatherer only) by specifying the 'y_column'.
#' This means that 'VhgBoxplot' can generate three different types of boxplots.
#' By default, 'y_column' is set to "ViralRefSeq_E" and will plot the Reference E-Value on the y-axis.
#' Grouping on the x-axis is done by the 'x_column' argument. By default, the "best_query" will be used.
#'
#' Additionally, the function calculates summary statistics and identifies outliers for further analysis ("ViralRefSeq_E" and "contig_len" only).
#' When 'y_column' is set to "ViralRefSeq_E", the output also includes 'rows_belowthres', which contains the hittable filtered for the rows below the threshold specified in the 'cut' argument.
#'
#' The 'cut' argument is used differently depending on the 'y_column' value:
#' - For 'y_column' set to "contig_len" or "ViralRefSeq_ident", the 'cut' argument filters the data to plot only the values with a "ViralRefSeq_E" below the specified threshold (default is 0.0001).
#' - For 'y_column' set to "ViralRefSeq_E", the rows are not filtered. Instead, a horizontal line (h_line) is shown in the plot to indicate the cutoff value.
#'
#' This allows the user to plot only the significant contig lengths and identities while also visualizing the number of non-significant and significant values for comparison.
#'
#' Warning: In some cases, E-values might be exactly 0. When these values are transformed using -log10, R
#' returns "inf" as the output. To avoid this issue, we replace all E-values that are 0 with the smallest E-value that is greater than 0.
#' If the smallest E-value is above the user-defined cutoff, we use a value of `cutoff * 10^-10` to replace the zeros.
#'
#' @return A list containing:
#' - the generated box plot.
#' - summary statistics.
#' - outliers ("ViralRefSeq_E" and "contig_len" only).
#' - and rows_belowthres ("ViralRefSeq_E" only).
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
#' # plot 1 for evalues
#' plot1 <- VhgBoxplot(file, x_column = "best_query", y_column = "ViralRefSeq_E")
#' plot1
#'
#' # plot 2 for identity
#' plot2 <- VhgBoxplot(file, x_column = "best_query", y_column = "ViralRefSeq_ident")
#' plot2
#'
#' # plot 3 custom arguments used
#' plot3 <- VhgBoxplot(file,
#'                   x_column = "best_query",
#'                   y_column = "ViralRefSeq_E",
#'                   theme_choice = "grey",
#'                   subtitle = "Custom subtitle: Identity for custom query",
#'                   xlabel = "Custom x-axis label: Custom query",
#'                   ylabel = "Custom y-axis label: Viral Reference Evalue in -log10 scale",
#'                   legend_position = "right")
#' plot3
#'
#' # import gatherer files
#' path2 <- system.file("extdata", "virusgatherer.tsv", package = "Virusparies")
#' vg_file <- ImportVirusTable(path2)
#'
#'
#' # plot 4: Virusgatherer plot for ViralRefSeq_taxonomy agains contig length
#' plot5 <- VhgBoxplot(vg_file,x_column = "ViralRefSeq_taxonomy",y_column = "contig_len")
#' plot5
#'
#'
#'
#' @import ggplot2
#' @importFrom rlang .data
#' @export
VhgBoxplot <- function(file,
                              x_column ="best_query",
                              taxa_rank = "Family",
                              y_column = "ViralRefSeq_E",
                              contiglen_log10_scale = FALSE,
                              cut = 1e-5,
                              add_cutoff_line = TRUE,
                              cut_colour = "#990000",
                              reorder_criteria = "median",
                              theme_choice = "linedraw",
                              flip_coords = TRUE,
                              add_mean_point = FALSE,
                              mean_color = "white",
                              mean_point_size = 2,
                              title = "default",
                              title_size = 16,
                              title_face = "bold",
                              title_colour = "#2a475e",
                              subtitle = "default",
                              subtitle_size = 12,
                              subtitle_face = "bold",
                              subtitle_colour = "#1b2838",
                              xlabel = NULL,
                              ylabel = NULL,
                              axis_title_size = 12,
                              xtext_size = 10,
                              ytext_size = 10,
                              remove_group_labels = FALSE,
                              legend_title = "Phylum",
                              legend_position = "bottom",
                              legend_title_size = 12,
                              legend_title_face = "bold",
                              legend_text_size = 10,
                              facet_ncol = NULL
                              ){



  #is_file_empty(file)
  if (is_file_empty(file)) {
    #message("Skipping VhgBoxplot generation due to empty data.")
    return(invisible(NULL))  # Return invisible(NULL) to stop further execution
  }


  if (!(x_column %in% c("best_query", "ViralRefSeq_taxonomy"))) {
    stop('Invalid value for x_column. Please use either "best_query" or "ViralRefSeq_taxonomy".')
  }


  check_columns(file,x_column)
  check_columns(file,y_column)
  check_input_type(file,y_column,2)
  check_input_type(file,x_column,1)

  # check arguments
  arg_character(theme_choice)
  arg_character(legend_position)

  arg_logical(flip_coords)





  # Find the smallest value greater than 0 in ViralRefSeq_E
  min_positive_value <- min(file$ViralRefSeq_E[file$ViralRefSeq_E > 0])

  # Check if the minimum positive value is below the cutoff
  if (min_positive_value < cut) {
    # Replace all 0 values with the smallest positive value
    file$ViralRefSeq_E[file$ViralRefSeq_E == 0] <- min_positive_value
  } else {

    # Add 10^-10 to the cutoff value and use that for the 0 values
    file$ViralRefSeq_E[file$ViralRefSeq_E == 0] <- cut * 10^-10
  }





  # Get the parameters based on y_column
  params <- get_plot_parameters(y_column, cut)

  # Access the cutoff and ylabel from the returned list
  cutoff <- params$cutoff
  default_ylabel <- params$ylabel


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
      #message("Skipping VhgBoxplot generation due to empty data.")
      return(invisible(NULL))  # Return invisible(NULL) to stop further execution
    }

  }

  if(x_column == "ViralRefSeq_taxonomy"){

    tax_column <- file$ViralRefSeq_taxonomy

    file <- VhgPreprocessTaxa(file,taxa_rank)

  }


  # change values of evalues to -log10
  if(y_column == "ViralRefSeq_E"){

    y_aes <- -log10(file[[y_column]])
  }else{
    y_aes <- file[[y_column]]
  }


  # add default subtitle for E-values
  if(y_column=="ViralRefSeq_E"){
    # define a cut off fot evalue significance
    default_sub <- paste0("Red line shows viral reference e-values under user-defined threshold: ",10^(-cutoff)," (-log10 scale: ",cutoff,")")

  }else{

    default_sub <- NULL
  }

  # set default titles
  default_titl <- switch(
    y_column,
    "ViralRefSeq_E" = "Boxplot of viral reference e-values for each group",
    "contig_len" = "Boxplot of contig length for each group",
    "ViralRefSeq_ident" = "Boxplot of viral reference identity for each group"
  )


  # Set the title
  title_text <- if (title == "default") default_titl else if (is.null(title) || title == "") NULL else title

  # Set the subtitle
  subtitle_text <- if (subtitle == "default") default_sub else if (is.null(subtitle) || subtitle == "") NULL else subtitle

  # Update xlabel to use user-provided label if provided
  xlabel <- ifelse(!is.null(xlabel), xlabel, "Viral Group")
  # Determine the y-axis label based on user-provided ylabel or default label
  ylabel <- ifelse(!is.null(ylabel), ylabel, default_ylabel)


  color_data <- consistentColourPalette(file, groupby = x_column,taxa_rank=taxa_rank)
  legend_labels <- color_data$legend_labels
  labels <- color_data$labels

  # Extract unique values from file$best_query
  unique_queries <- unique(file[[x_column]])

  # Create a vector of corresponding names from legend_labels
  pyhlum_names <- legend_labels[unique_queries]

  # Match names to file$best_query
  file$phylum <- pyhlum_names[match(file[[x_column]], unique_queries)]


  # Check for valid reorder_criteria
  valid_criteria <- c("max", "min", "median", "mean")
  if (!is.null(reorder_criteria) && !(reorder_criteria %in% valid_criteria)) {
    stop("Invalid reorder_criteria. Please choose one of: max, min, median, mean.")
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




  ########################
  ### generate boxplot ###
  ########################

  # print message to console.
  plot_boxplot_message(y_column=y_column,x_column=x_column,cutoff)







  boxp <- ggplot(file,aes(x= if (!is.null(reorder_criteria)) {
    reorder(.data[[x_column]], y_aes,
            FUN = switch(reorder_criteria,
                         "max" = max,
                         "min" = min,
                         "median" = median,
                         "mean" = mean))
  } else {
    factor(.data[[x_column]], levels = rev(unique(sort(.data[[x_column]]))))
  },
                             y=y_aes,fill=.data$phylum))+
    geom_boxplot(staplewidth = 0.4)+
    labs(x=xlabel,
         y=ylabel,
         title=title_text,
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
      axis.text.x =  element_text(size = xtext_size),
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
    scale_y_continuous(expand = c(0, 0), limits = y_limits)+
    theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

  if (flip_coords) {
    boxp <- boxp + coord_flip()
  }

  if(y_column=="ViralRefSeq_E" && add_cutoff_line){

    boxp <- boxp+geom_hline(aes(yintercept=cutoff), colour=cut_colour)
  }

  if(y_column == "contig_len" && contiglen_log10_scale){
    boxp <- boxp + scale_y_log10()
  }

  if(x_column != "SRA_run"){

    # matched_vector <- consistentColourPalette(file = file, groupby = x_column)
    boxp <- boxp + scale_fill_manual(values = labels)


  }


  if (add_mean_point) {

    # Update the aes() inside geom_point to correctly reference mean_values
    boxp <- boxp + stat_summary(fun=mean, geom="point", size=mean_point_size,
                                shape = 21,      # Use shape 21 for filled circle with outline
                                fill = mean_color,  # Set fill color to white
                                color = "black",show.legend = FALSE)


  }

  boxp <- remove_group_text(boxp,remove_group_labels,flip_coords)

  boxp <- facet_plot(boxp,facet_ncol,flip_coords)












  # delete new column before sum statistics and export.
  file$phylum <- NULL




  summary_stats <- boxp_summary_stats(file, group = x_column,ycol =y_column)

  if(x_column == "ViralRefSeq_taxonomy"){


    file$ViralRefSeq_taxonomy <- tax_column

  }

  if (y_column == "ViralRefSeq_E") {
    outlier <- find_outlier_eval_box(file, group = x_column)
  } else if (y_column == "contig_len") {
    outlier <- find_outlier_eval_box(file, group = x_column, y_column = y_column)
  }







  #plot(boxp)

  # Prepare the results list
  results <- list(boxp = boxp, summary_stats = summary_stats)

  # Add outliers if applicable
  if (y_column == "ViralRefSeq_E" || y_column == "contig_len") {
    results$outlier <- outlier
  }

  # extract rows below threshold
  if (y_column == "ViralRefSeq_E") {
    belowthres_boxp <- vhg_filter_belowthresholdboxplot(file,cut)
    results$rows_belowthres <- belowthres_boxp
  }

  #results$matched_vector <- matched_vector


  # Return the results list
  message("Box plot generation completed.")
  return(results)






}







