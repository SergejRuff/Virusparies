#' vhRunsBarplot: Generate a bar plot showing the number of datasets with unique runs found for
#' each Virus group.
#'
#' @param vh_file Dataframe containing Virushunter hittables results.
#' @param cut (optional) The significance cutoff value for E-values (default: 1e-5).
#' Removes rows in vh_file with values larger than cutoff value in ViralRefSeq_E column.
#' @param theme_choice (optional) A character indicating the ggplot2 theme to apply.
#' Options include "minimal", "classic", "light", "dark", "void", "grey" (or "gray"), "bw", "linedraw", and "test". Default is "minimal".
#' @param flip_coords (optional) Logical indicating whether to flip the coordinates of the plot. Default is TRUE.
#' @param title (optional) A character specifying the title of the plot. Default is "Number of datasets with hits found for each Virus group".
#' @param title_size (optional) The size of the plot title. Default is 16.
#' @param title_face (optional) The font face of the plot title. Default is "bold".
#' @param title_colour (optional) The color of the plot title. Default is "#2a475e".
#' @param subtitle (optional) A character specifying the subtitle of the plot.
#' Default is "default", which calculates the total number of datasets with hits and returns it as
#' "total number of datasets with hits: " followed by the calculated number. an empty string ("") removes the subtitle.
#' @param subtitle_size (optional) The size of the plot subtitle. Default is 12.
#' @param subtitle_face (optional) The font face of the plot subtitle. Default is "bold".
#' @param subtitle_colour (optional) The color of the plot subtitle. Default is "#1b2838".
#' @param xlabel (optional) A character specifying the label for the x-axis. Default is "Virus family found in query".
#' @param ylabel (optional) A character specifying the label for the y-axis. Default is "Number of datasets with hits".
#' @param axis_title_size (optional) The size of axis titles. Default is 12.
#' @param xtext_size (optional) The size of x-axis text labels. Default is 10.
#' @param ytext_size (optional) The size of y-axis text labels. Default is 10.
#' @param legend_title (optional) A character specifying the title of the legend. Default is "virus family".
#' @param legend_position (optional) The position of the legend. Default is "bottom".
#' @param legend_title_size (optional) The size of the legend title. Default is 12.
#' @param legend_title_face (optional) The font face of the legend title. Default is "bold".
#' @param legend_text_size (optional) The size of the legend text. Default is 10.
#' @param plot_text_size (optional) The size of the text labels added to the plot. Default is 3.5.
#' @param plot_text_position_dodge (optional) The degree of dodging for positioning text labels. Default is 0.9.
#' @param plot_text_hjust (optional) The horizontal justification of text labels. Default is -0.1.
#' @param plot_text_vjust (optional) The vertical justification of text labels. Default is 0.5.
#' It is recommended to change `vjust` when setting `flip_coords = FALSE`.
#' @param plot_text_colour (optional) The color of the text labels added to the plot. Default is "black".
#'
#' @return A list containing the bar plot and optionally the generated table and processed data
#'
#' @examples
#' path <- system.file("extdata", "virushunter.tsv", package = "Virusparies")
#' vh_file <- importVirusTable(path)
#'
#' # plot 1: plot boxplot for "identity"
#' plot <- vhRunsBarplot(vh_file,cut = 1e-5)
#'
#' # return sample_run inside plot object
#' print(plot$sample_run)
#'
#' # Plot 2: Customized plot with modified settings
#' plot_custom <- vhRunsBarplot(
#'   vh_file,
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
#'
#' @import ggplot2
#' @importFrom rlang .data
#' @export
vhRunsBarplot <- function(vh_file,cut = 1e-5,
                          theme_choice = "minimal",
                          flip_coords = TRUE,
                          title = "Number of datasets with hits found for each Virus group",
                          title_size = 16,
                          title_face = "bold",
                          title_colour = "#2a475e",
                          subtitle = "default",
                          subtitle_size = 12,
                          subtitle_face = "bold",
                          subtitle_colour = "#1b2838",
                          xlabel = "Virus family found in query",
                          ylabel = "Number of datesets with hits",
                          axis_title_size = 12,
                          xtext_size = 10,
                          ytext_size = 10,
                          legend_title = "virus family",
                          legend_position = "bottom",
                          legend_title_size = 12,
                          legend_title_face = "bold",
                          legend_text_size = 10,
                          plot_text_size = 3.5,
                          plot_text_position_dodge = 0.9,
                          plot_text_hjust = -0.1,
                          plot_text_vjust = 0.5,
                          plot_text_colour = "black"
                          ){

  vh_file <- vh_file[vh_file$ViralRefSeq_E < cut,]

  # Apply the selected theme
  theme_selected <- select_theme(theme_choice)


  # preprocess data for plot
  sample_run <- preprocess_runs_bar(vh_file)

  # Set the subtitle based on the input
  if (subtitle == "default") {
    subtitle_text <- paste0("total number of datasets with hits: ", sum(sample_run$unique_SRA_run))
  } else if (is.null(subtitle) || subtitle == "") {
    subtitle_text <- NULL
  } else {
    subtitle_text <- subtitle
  }



  #####################
  ### Generate Plot ###
  ### #################

  run_bar <- ggplot(data = sample_run,aes(x=reorder(.data$best_query,.data$unique_SRA_run,FUN=max),y=.data$unique_SRA_run,fill= .data$best_query))+
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
    geom_text(aes(label = .data$res),
              hjust = plot_text_hjust,
              vjust= plot_text_vjust,
              color = plot_text_colour,
              position = position_dodge(plot_text_position_dodge),
              size = plot_text_size)+
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
      ))

  if (flip_coords) {
    run_bar <- run_bar + coord_flip()
  }


  #plot(run_bar)



  return(list(run_bar=run_bar,
              sample_run=sample_run))


}





