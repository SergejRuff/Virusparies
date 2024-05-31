#' @title VhgBoxplot: Generate boxplots comparing E-values or identity for each virus group
#'
#' @description
#'  This function generates boxplots comparing either E-values or identity
#' for each group from VirusHunter or VirusGatherer Hittables results.
#'
#' @param vh_file A data frame containing VirusHunter or VirusGatherer Hittable results.
#' @param x_column A character specifying the column containing the groups (Default:"best_query").
#' @param y_column A character specifying the column containing the values to be compared. Currently "ViralRefSeq_ident",
#' "contig_len" (Column in Gatherer hittable) and "ViralRefSeq_E" are supported columns (Default:"ViralRefSeq_E").
#' @param cut (optional) The significance cutoff value for E-values (default: 1e-5).
#' @param cut_colour (optional) The color for the significance cutoff line (default: "#990000").
#' @param theme_choice (optional): A character indicating the ggplot2 theme to apply. Options include "minimal",
#'  "classic", "light", "dark", "void", "grey" (or "gray"), "bw", "linedraw", and "test".
#'  Default is "minimal".
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
#' @param legend_title (optional) A character specifying the title for the legend. Default is "virus family".
#' @param legend_position (optional) A character specifying the position of the legend. Default is "bottom".
#' @param legend_title_size (optional) Numeric specifying the size of the legend title text. Default is 12.
#' @param legend_title_face (optional) A character specifying the font face for the legend title text. Default is "bold".
#' @param legend_text_size (optional) Numeric specifying the size of the legend text. Default is 10.
#'
#'
#' @return A list containing the generated boxplot, summary statistics, and outliers
#'
#' @details This function generates boxplots comparing either E-values or identity for each virus
#' group from the VirusHunter or Gatherer Hittable.
#' It also calculates summary statistics and  identifies outliers for further analysis.
#' The user can specify whether to generate boxplots for E-values or identity by setting
#' the parameter .
#' @examples
#' path <- system.file("extdata", "virushunter.tsv", package = "Virusparies")
#' vh_file <- importVirusTable(path)
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
#' vg_file <- importVirusTable(path2)
#'
#' # plot 4: virusgatherer plot with SRA_run as custom grouping
#' plot4 <- VhgBoxplot(vg_file,x_column = "SRA_run",y_column = "ViralRefSeq_E")
#' plot4
#'
#' # plot 5: Virusgatherer plot for SRA_runs agains contig length
#' plot5 <- VhgBoxplot(vg_file,x_column = "SRA_run",y_column = "contig_len")
#' plot5
#'
#'
#' @import ggplot2
#' @importFrom rlang .data
#' @export
VhgBoxplot <- function(vh_file,
                              x_column ="best_query",
                              y_column = "ViralRefSeq_E",
                              cut = 1e-5,
                              cut_colour = "#990000",
                              theme_choice = "minimal",
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
                              legend_title = "virus family",
                              legend_position = "bottom",
                              legend_title_size = 12,
                              legend_title_face = "bold",
                              legend_text_size = 10
                              ){

  all_names <- names(vh_file)
  # Check if x_column and y_column exist in vh_file
  if (!(x_column %in% all_names)) {
    stop("x_column '", x_column, "' not found in vh_file. Available column names: ", paste(all_names, collapse = ", "))
  }
  if (!(y_column %in% all_names)) {
    stop("y_column '", y_column, "' not found in vh_file. Available column names: ", paste(all_names, collapse = ", "))
  }

  is_vh_file_empty(vh_file)

  # Find the smallest value greater than 0 in ViralRefSeq_E
  min_positive_value <- min(vh_file$ViralRefSeq_E[vh_file$ViralRefSeq_E > 0])

  # Replace all 0 values with the smallest positive value
  vh_file$ViralRefSeq_E[vh_file$ViralRefSeq_E == 0] <- min_positive_value






  if (y_column == "ViralRefSeq_E") {
    # define a cut off for evalue significance
    cutoff <- -log10(cut)
    default_ylabel <- "-log10 of viral Reference E-values"  # Default y-axis label
  } else if (y_column == "ViralRefSeq_ident") {
    cutoff <- NULL
    default_ylabel <- "Viral Reference Identity in %"  # Default y-axis label
  } else if (y_column == "contig_len") {
    cutoff <- NULL
    default_ylabel <- "Contig Length"  # Default y-axis label
  } else {
    stop("Invalid y_column value provided.")
  }


  # Apply the selected theme
  theme_selected <- select_theme(theme_choice)

  # set cutoff for contig_len and identity
  if(y_column=="contig_len" | y_column== "ViralRefSeq_ident"){

    message(paste(y_column," column was selected."))
    message(paste("Removing values with E-values higher than: ",cut,"."))

    vh_file <- vh_file[vh_file$ViralRefSeq_E < cut,]

    message(paste0("after removing rows based on evalue the hittable has: ",nrow(vh_file)," rows left."))
    is_vh_file_empty(vh_file)

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
    default_sub <- paste0("red line shows viral Refrence E-values under user-defined threshold: ",10^(-cutoff)," (-log10 scale: ",cutoff,")")

  }else{

    default_sub <- NULL
  }

  # set default titles
  default_titl <- switch(
    y_column,
    "ViralRefSeq_E" = "Boxplot plotting viral Reference E-Values for each virus family",
    "contig_len" = "Boxplot plotting contig length for each virus family",
    "ViralRefSeq_ident" = "Boxplot plotting viral Reference Identity for each virus family"
  )


  # Set the title
  title_text <- if (title == "default") default_titl else if (is.null(title) || title == "") NULL else title

  # Set the subtitle
  subtitle_text <- if (subtitle == "default") default_sub else if (is.null(subtitle) || subtitle == "") NULL else subtitle

  # Update xlabel to use user-provided label if provided
  xlabel <- ifelse(!is.null(xlabel), xlabel, "Group in query")
  # Determine the y-axis label based on user-provided ylabel or default label
  ylabel <- ifelse(!is.null(ylabel), ylabel, default_ylabel)



  ########################
  ### generate boxplot ###
  ########################

  # print message to console.
  plot_boxplot_message(y_column=y_column,x_column=x_column,cutoff)





  boxp <- ggplot(vh_file,aes(x=reorder(.data[[x_column]],y_aes,FUN=median),
                             y=y_aes,fill=.data[[x_column]]))+
    geom_boxplot(staplewidth = 0.4)+
    labs(x=xlabel,
         y=ylabel,
         title=title_text,
         subtitle = subtitle_text)+
    theme_selected+
    theme(legend.position = legend_position)+
    guides(fill=guide_legend(title=legend_title))+
    coord_flip()+
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






  if(y_column=="ViralRefSeq_E"& x_column=="best_query"){


    summary_stats <- vh_sum_stat_evavlue_boxplot(vh_file,cutoff)
    outlier <- find_outlier_eval_box(vh_file)


  }



  if(y_column=="ViralRefSeq_ident"& x_column=="best_query"){



    summary_stats <- summary_stats_identity(vh_file)


  }

  if(y_column=="ViralRefSeq_E"& x_column != "best_query"){
    summary_stats <- summary_stats_identity(vh_file,group=x_column,ycol =y_column)
    outlier <- find_outlier_eval_box(vh_file,group=x_column)
  }

  if(y_column=="ViralRefSeq_ident"& x_column != "best_query"){
    summary_stats <- summary_stats_identity(vh_file,group=x_column)
  }

  if(y_column=="ViralRefSeq_ident"& x_column != "best_query"){
    summary_stats <- summary_stats_identity(vh_file,group=x_column)
  }

  if(y_column=="contig_len"& x_column != "best_query"){
    summary_stats <- summary_stats_identity(vh_file,group=x_column,ycol =y_column)
    outlier <- find_outlier_eval_box(vh_file,group=x_column,y_column=y_column)
  }




  #plot(boxp)

  # Prepare the results list
  results <- list(boxp = boxp, summary_stats = summary_stats)

  # Add outliers if applicable
  if (y_column == "ViralRefSeq_E" || y_column == "contig_len") {
    results$outlier <- outlier
  }

  # Return the results list
  return(results)




}







