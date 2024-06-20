#' @title VhgBoxplot: Generate boxplots comparing E-values,identity or contig length (Gatherer only) for each virus group
#'
#' @description
#'  This function generates boxplots comparing either E-values,identity or contig length (Gatherer only)
#' for each group from VirusHunter or VirusGatherer Hittables results.
#'
#' @param vh_file A data frame containing VirusHunter or VirusGatherer Hittable results.
#' @param x_column A character specifying the column containing the groups (Default:"best_query").
#' @param y_column A character specifying the column containing the values to be compared. Currently "ViralRefSeq_ident",
#' "contig_len" (Column in Gatherer hittable) and "ViralRefSeq_E" are supported columns (Default:"ViralRefSeq_E").
#' @param cut (optional) The significance cutoff value for E-values (default: 1e-5).
#' @param cut_colour (optional) The color for the significance cutoff line (default: "#990000").
#' @param reorder_criteria Character string specifying the criteria for reordering the x-axis ('max', 'min', 'median'(Default),'mean').
#' NULL sorts alphabetically.
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
#' @author Sergej Ruff
#' @details
#' VhgBoxplot generates boxplots comparing either E-values, identity, or contig length (Gatherer only) for each virus group from the VirusHunter or Gatherer Hittable.
#'
#' The user can specify whether to generate boxplots for E-values, identity, or contig length (Gatherer only) by specifying the 'y_column'.
#' This means that 'VhgBoxplot' can generate three different types of boxplots.
#' By default, 'y_column' is set to "ViralRefSeq_E" and will plot the Reference E-Value on the y-axis.
#' Grouping on the x-axis is done by the 'x_column' argument. By default, the "best_query" will be used.
#'
#' Additionally, the function calculates summary statistics and identifies outliers for further analysis ("ViralRefSeq_E" and "contig_len" only).
#' When 'y_column' is set to "ViralRefSeq_E", the output also includes 'rows_belowthres', which contains the Hittable filtered for the rows below the threshold specified in the 'cut' argument.
#'
#' The 'cut' argument is used differently depending on the 'y_column' value:
#' - For 'y_column' set to "contig_len" or "ViralRefSeq_ident", the 'cut' argument filters the data to plot only the values with a "ViralRefSeq_E" below the specified threshold (default is 0.0001).
#' - For 'y_column' set to "ViralRefSeq_E", the rows are not filtered. Instead, a horizontal line (h_line) is shown in the plot to indicate the cutoff value.
#'
#' This allows the user to plot only the significant contig lengths and identities while also visualizing the number of non-significant and significant values for comparison.
#'
#'
#'
#' @return A list containing the generated boxplot, summary statistics,outliers ("ViralRefSeq_E" and "contig_len" only) and rows_belowthres ("ViralRefSeq_E" only).
#'
#' @examples
#' path <- system.file("extdata", "virushunter.tsv", package = "Virusparies")
#' vh_file <- ImportVirusTable(path)
#'
#' # plot 1 for evalues
#' plot1 <- VhgBoxplot(vh_file, x_column = "best_query", y_column = "ViralRefSeq_E")
#' plot1
#'
#' # plot 2 for identity
#' plot2 <- VhgBoxplot(vh_file, x_column = "best_query", y_column = "ViralRefSeq_ident")
#' plot2
#'
#' # plot 3 custom arguments used
#' plot3 <- VhgBoxplot(vh_file,
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
#' @seealso
#' VirusHunterGatherer is available here: \url{https://github.com/lauberlab/VirusHunterGatherer}.
#'
#' @import ggplot2
#' @importFrom rlang .data
#' @export
VhgBoxplot <- function(vh_file,
                              x_column ="best_query",
                              y_column = "ViralRefSeq_E",
                              cut = 1e-5,
                              cut_colour = "#990000",
                              reorder_criteria = "median",
                              theme_choice = "linedraw",
                              flip_coords = TRUE,
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
                              legend_title = "Phylum",
                              legend_position = "bottom",
                              legend_title_size = 12,
                              legend_title_face = "bold",
                              legend_text_size = 10,
                              colorblind_support = FALSE,
                              colormap = "viridis"
                              ){



  is_file_empty(vh_file)


  if (!(x_column %in% c("best_query", "ViralRefSeq_taxonomy"))) {
    stop('Invalid value for x_column. Please use either "best_query" or "ViralRefSeq_taxonomy".')
  }


  check_columns(vh_file,x_column)
  check_columns(vh_file,y_column)
  check_input_type(vh_file,y_column,2)
  check_input_type(vh_file,x_column,1)

  # check arguments
  arg_character(theme_choice)
  arg_character(legend_position)
  arg_character(colormap)
  arg_logical(flip_coords)
  arg_logical(colorblind_support)




  # Find the smallest value greater than 0 in ViralRefSeq_E
  min_positive_value <- min(vh_file$ViralRefSeq_E[vh_file$ViralRefSeq_E > 0])

  # Replace all 0 values with the smallest positive value
  vh_file$ViralRefSeq_E[vh_file$ViralRefSeq_E == 0] <- min_positive_value





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

    vh_file <- vh_file[vh_file$ViralRefSeq_E < cut,]

    message(paste0("after removing rows based on evalue the hittable has ",nrow(vh_file)," rows left."))
    is_file_empty(vh_file)

  }

  if(x_column == "ViralRefSeq_taxonomy"){

    tax_column <- vh_file$ViralRefSeq_taxonomy

    vh_file <- taxonomy_group_preprocess(vh_file)

  }


  # change values of evalues to -log10
  if(y_column == "ViralRefSeq_E"){

    y_aes <- -log10(vh_file[[y_column]])
  }else{
    y_aes <- vh_file[[y_column]]
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


  color_data <- consistentColourPalette(vh_file, groupby = x_column)
  legend_labels <- color_data$legend_labels
  labels <- color_data$labels

  # Extract unique values from vh_file$best_query
  unique_queries <- unique(vh_file[[x_column]])

  # Create a vector of corresponding names from legend_labels
  pyhlum_names <- legend_labels[unique_queries]

  # Match names to vh_file$best_query
  vh_file$phylum <- pyhlum_names[match(vh_file[[x_column]], unique_queries)]


  # Check for valid reorder_criteria
  valid_criteria <- c("max", "min", "median", "mean")
  if (!is.null(reorder_criteria) && !(reorder_criteria %in% valid_criteria)) {
    stop("Invalid reorder_criteria. Please choose one of: max, min, median, mean.")
  }




  ########################
  ### generate boxplot ###
  ########################

  # print message to console.
  plot_boxplot_message(y_column=y_column,x_column=x_column,cutoff)







  boxp <- ggplot(vh_file,aes(x= if (!is.null(reorder_criteria)) {
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

  if (flip_coords) {
    boxp <- boxp + coord_flip()
  }

  if(y_column=="ViralRefSeq_E"){

    boxp <- boxp+geom_hline(aes(yintercept=cutoff), colour=cut_colour)
  }

  if(x_column != "SRA_run"){

    # matched_vector <- consistentColourPalette(vh_file = vh_file, groupby = x_column)
    boxp <- boxp + scale_fill_manual(values = labels)


  }




  # add colorblind support
  if(colorblind_support){
    boxp<- colorbildsupport(boxp,colormap)
  }

  # delete new column before sum statistics and export.
  vh_file$phylum <- NULL




  summary_stats <- boxp_summary_stats(vh_file, group = x_column,ycol =y_column)

  if(x_column == "ViralRefSeq_taxonomy"){


    vh_file$ViralRefSeq_taxonomy <- tax_column

  }

  if (y_column == "ViralRefSeq_E") {
    outlier <- find_outlier_eval_box(vh_file, group = x_column)
  } else if (y_column == "contig_len") {
    outlier <- find_outlier_eval_box(vh_file, group = x_column, y_column = y_column)
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
    belowthres_boxp <- vhg_filter_belowthresholdboxplot(vh_file,cut)
    results$rows_belowthres <- belowthres_boxp
  }

  #results$matched_vector <- matched_vector


  # Return the results list
  message("Boxplot generation completed.")
  return(results)






}







