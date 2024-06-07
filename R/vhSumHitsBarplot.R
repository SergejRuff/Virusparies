#' @title vhSumHitsBarplot: Generate a bar plot showing the sum of hits for each virus family
#'
#' @description
#'  This function preprocesses virus data for plotting and generates a bar plot showing the sum of hits/ contigs
#' for each virus family from the input dataset.
#'
#' @param vh_file A data frame containing VirusHunters hittables results.
#' @param cut (optional) The significance cutoff value for E-values (default: 1e-5). Removes rows in vh_file
#' with values larger than cutoff value in ViralRefSeq_E column.
#' @param theme_choice (optional) A character indicating the ggplot2 theme to apply. Options include "minimal",
#'  "classic", "light", "dark", "void", "grey" (or "gray"), "bw", "linedraw", and "test".
#'  Default is "linedraw".
#'
#' @param flip_coords (optional) Logical indicating whether to flip the coordinates of the plot. Default is TRUE.
#' @param title (optional) A character specifying the title of the plot. Default is "sum of hits for each family".
#' @param title_size (optional) The size of the plot title. Default is 16.
#' @param title_face (optional) The font face of the plot title. Default is "bold".
#' @param title_colour (optional) The color of the plot title. Default is "#2a475e".
#' @param subtitle (optional) A character specifying the subtitle of the plot.
#' Default is "total number of hits: " followed by the calculated number.empty string ("") removes subttitle.
#' @param subtitle_size (optional) The size of the plot subtitle. Default is 12.
#' @param subtitle_face (optional) The font face of the plot subtitle. Default is "bold".
#' @param subtitle_colour (optional) The color of the plot subtitle. Default is "#1b2838".
#' @param xlabel (optional) A character specifying the label for the x-axis. Default is "Virus family found in query".
#' @param ylabel (optional) A character specifying the label for the y-axis. Default is "sum of hits".
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
#'
#' @details
#' vhSumHitsBarplot preprocesses virus data for plotting by calculating the sum of
#' hits for each virus family from the input dataset (accepts only VirusHunter Hittables).
#' It then generates a bar plot showing the sum of hits/contigs for each virus family.
#' Additionally, it returns the processed data for further analysis.
#'
#' @return A list containing the generated bar plot and processed data
#'
#' @author Sergej Ruff
#' @examples
#' path <- system.file("extdata", "virushunter.tsv", package = "Virusparies")
#' vh_file <- importVirusTable(path)
#'
#' # plot 1: plot boxplot for "identity"
#' plot <- vhSumHitsBarplot(vh_file,cut = 1e-5)
#' plot
#'
#' # return vh_group inside plot object
#' print(plot$vh_group)
#'
#'
#'# plot 2: Customized plot
#' plot2 <- vhSumHitsBarplot(
#'   vh_file,
#'   cut = 1e-5,
#'   title = "Customized Hits Plot",
#'   subtitle = "Total Hits with Significance Cutoff",
#'   xlabel = "Virus Family",
#'   ylabel = "Total Hits",
#'   legend_title = "Virus Families",
#'   plot_text_size = 4,
#'   plot_text_colour = "blue"
#' )
#'
#' plot2
#'
#' @seealso
#' VirusHunterGatherer is available here: \url{https://github.com/lauberlab/VirusHunterGatherer}.
#'
#' @import ggplot2
#' @importFrom rlang .data
#' @importFrom stats reorder
#' @export
vhSumHitsBarplot <- function(vh_file,cut = 1e-5,
                             theme_choice = "linedraw",
                             flip_coords = TRUE,
                             title = "sum of hits for each family",
                             title_size = 16,
                             title_face = "bold",
                             title_colour = "#2a475e",
                             subtitle = "default",
                             subtitle_size = 12,
                             subtitle_face = "bold",
                             subtitle_colour = "#1b2838",
                             xlabel = "Virus family found in query",
                             ylabel = "sum of hits",
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
                             plot_text_colour = "black",
                             colorblind_support = FALSE,
                             colormap = "viridis"){


 #check if table is empty
 is_file_empty(vh_file)


 # Define the required columns
 required_columns <- c("num_hits", "ViralRefSeq_E", "best_query")

 check_columns(vh_file,required_columns)
 check_input_type(vh_file,c("num_hits", "ViralRefSeq_E"),2)
 check_input_type(vh_file,"best_query",1)

 # check arguments
 arg_character(theme_choice)
 arg_character(legend_position)
 arg_character(colormap)
 arg_logical(flip_coords)
 arg_logical(colorblind_support)

 ## preprocess data for plotting
 vh_file <- vh_file[vh_file$ViralRefSeq_E < cut,]
 is_file_empty(vh_file)
 vh_group <- vh_sumhitbar_preprocessing(vh_file)

 # Apply the selected theme
 theme_selected <- select_theme(theme_choice)



 # Set the subtitle based on the input
 if (subtitle == "default") {
    subtitle_text <- paste0("total number of hits: ",sum(vh_group$sum))
 } else if (is.null(subtitle) || subtitle == "") {
    subtitle_text <- NULL
 } else {
    subtitle_text <- subtitle
 }



 ########################
 ### plot sum of hits ###
 ########################
 ########################

 sum_plot <- ggplot(data = vh_group,aes(x=reorder(.data$best_query,.data$sum,FUN=max),y= .data$sum,fill=.data$best_query))+
   geom_bar(stat = "identity")+
   labs(title = title,
        x = xlabel,
        y = ylabel,
        subtitle = subtitle_text)+
   theme_selected+
   theme(legend.position = legend_position,
         axis.text.y = element_text(size = ytext_size),
         axis.text.x = element_text(size = xtext_size),
         axis.title = element_text(size = axis_title_size),
         legend.text = element_text(size = legend_text_size),
         legend.key.size = unit(1.5, "lines"),
         legend.title = element_text(size = legend_title_size, face = legend_title_face))+
   guides(fill=guide_legend(title=legend_title)) +
   geom_text(aes(label = ifelse(sum == 0, "",.data$res)),
             hjust = plot_text_hjust,
             vjust = plot_text_vjust,
             color = plot_text_colour,
             position = position_dodge(plot_text_position_dodge),
             size = plot_text_size)+
   scale_y_continuous(expand = c(0, 0), limits = c(0, max(vh_group$sum)+max(vh_group$sum)/8))+
   theme(
     plot.title = element_text(
        size = title_size,
        face = title_face,
        color = title_colour),
     # Statistical annotations below the main title
     plot.subtitle = element_text(
       size = subtitle_size,
        face = subtitle_face,
        color = subtitle_colour
     ))

 if (flip_coords) {
    sum_plot <- sum_plot + coord_flip()
 }

 # add colorblind support
 if(colorblind_support){
     sum_plot<- colorbildsupport(sum_plot,colormap)
 }



 #plot(sum_plot)
 message("Barplot generation completed.")

 return(list(plot=sum_plot,vh_group=vh_group))



}



