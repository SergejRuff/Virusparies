#' @title VhSumHitsBarplot: Generate a bar plot showing the sum of hits for each virus family
#'
#' @description
#' This function preprocesses virus data for plotting and generates a bar plot showing the sum of hits/ contigs
#' for each virus family from the input dataset.
#'
#' @param vh_file A data frame containing VirusHunters hittables results.
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
#' @param cut (optional) The significance cutoff value for E-values (default: 1e-5). Removes rows in vh_file
#' with values larger than cutoff value in ViralRefSeq_E column.
#' @param reorder_criteria Character string specifying the criteria for reordering the x-axis ('max' (default), 'min').
#' NULL sorts alphabetically.
#' @param theme_choice (optional) A character indicating the ggplot2 theme to apply. Options include "minimal",
#'  "classic", "light", "dark", "void", "grey" (or "gray"), "bw", "linedraw", and "test".
#'  Default is "linedraw".
#'
#' @param flip_coords (optional) Logical indicating whether to flip the coordinates of the plot. Default is TRUE.
#' @param title (optional) A character specifying the title of the plot. Default is "Distribution of hits for each virus group".
#' @param title_size (optional) The size of the plot title. Default is 16.
#' @param title_face (optional) The font face of the plot title. Default is "bold".
#' @param title_colour (optional) The color of the plot title. Default is "#2a475e".
#' @param subtitle (optional) A character specifying the subtitle of the plot.
#' Default is "total number of hits: " followed by the calculated number.empty string ("") removes subttitle.
#' @param subtitle_size (optional) The size of the plot subtitle. Default is 12.
#' @param subtitle_face (optional) The font face of the plot subtitle. Default is "bold".
#' @param subtitle_colour (optional) The color of the plot subtitle. Default is "#1b2838".
#' @param xlabel (optional) A character specifying the label for the x-axis. Default is "Viral group".
#' @param ylabel (optional) A character specifying the label for the y-axis. Default is "Total number of hits".
#' @param axis_title_size (optional) The size of axis titles. Default is 12.
#' @param xtext_size (optional) The size of x-axis text labels. Default is 10.
#' @param ytext_size (optional) The size of y-axis text labels. Default is 10.
#' @param remove_group_labels (optional): If `TRUE`, the group labels will be removed; if `FALSE` or omitted, the labels will be displayed.
#' @param legend_title (optional) A character specifying the title of the legend. Default is "Phylum".
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
#'
#'
#' @details
#' VhSumHitsBarplot preprocesses virus data for plotting by calculating the sum of
#' hits for each virus family from the input dataset (accepts only VirusHunter Hittables).
#' It then generates a bar plot showing the sum of hits/contigs for each virus family.
#' Additionally, it returns the processed data for further analysis.
#'
#' @return A list containing the generated bar plot and processed data.
#'
#' @author Sergej Ruff
#' @examples
#' path <- system.file("extdata", "virushunter.tsv", package = "Virusparies")
#' vh_file <- ImportVirusTable(path)
#'
#' # plot 1: plot boxplot for "identity"
#' plot <- VhSumHitsBarplot(vh_file,cut = 1e-5)
#' plot
#'
#' # return vh_group inside plot object
#' print(plot$vh_group)
#'
#'
#'# plot 2: Customized plot
#' plot2 <- VhSumHitsBarplot(
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
VhSumHitsBarplot <- function(vh_file,
                             groupby = "best_query",
                             taxa_rank = "Family",
                             cut = 1e-5,
                             reorder_criteria = "max",
                             theme_choice = "linedraw",
                             flip_coords = TRUE,
                             title = "Distribution of hits for each virus group",
                             title_size = 16,
                             title_face = "bold",
                             title_colour = "#2a475e",
                             subtitle = "default",
                             subtitle_size = 12,
                             subtitle_face = "bold",
                             subtitle_colour = "#1b2838",
                             xlabel = "Viral group",
                             ylabel = "Total number of hits",
                             axis_title_size = 12,
                             xtext_size = 10,
                             ytext_size = 10,
                             remove_group_labels = FALSE,
                             legend_title = "Phylum",
                             legend_position = "bottom",
                             legend_title_size = 12,
                             legend_title_face = "bold",
                             legend_text_size = 10,
                             plot_text_size = 3.5,
                             plot_text_position_dodge = 0.9,
                             plot_text_hjust = -0.1,
                             plot_text_vjust = 0.5,
                             plot_text_colour = "black"){


 #check if table is empty
 #is_file_empty(vh_file)
 if (is_file_empty(vh_file)) {
    #message("Skipping VhgBoxplot generation due to empty data.")
    return(invisible(NULL))  # Return invisible(NULL) to stop further execution
  }




 # Define the required columns
 required_columns <- c("num_hits", "ViralRefSeq_E", groupby)

 if (!(groupby %in% c("best_query", "ViralRefSeq_taxonomy"))) {
     stop('Invalid value for groupby. Please use either "best_query" or "ViralRefSeq_taxonomy".')
 }

 check_columns(vh_file,required_columns)
 check_input_type(vh_file,c("num_hits", "ViralRefSeq_E"),2)
 check_input_type(vh_file,groupby,1)


 # check arguments
 arg_character(theme_choice)
 arg_character(legend_position)

 arg_logical(flip_coords)


 if(groupby == "ViralRefSeq_taxonomy"){

     vh_file <- VhgPreprocessTaxa(vh_file,taxa_rank)

 }


 ## preprocess data for plotting
 vh_file <- vh_file[vh_file$ViralRefSeq_E < cut,]
 #is_file_empty(vh_file)
 if (is_file_empty(vh_file)) {
   #message("Skipping VhgBoxplot generation due to empty data.")
   return(invisible(NULL))  # Return invisible(NULL) to stop further execution
 }
 vh_group <- vh_sumhitbar_preprocessing(vh_file,groupby)

 # Apply the selected theme
 theme_selected <- select_theme(theme_choice)



 # Set the subtitle based on the input
 if (subtitle == "default") {
    subtitle_text <- paste0("Total number of hits: ",sum(vh_group$sum))
 } else if (is.null(subtitle) || subtitle == "") {
    subtitle_text <- NULL
 } else {
    subtitle_text <- subtitle
 }


 color_data <- consistentColourPalette(vh_file, groupby = groupby,taxa_rank=taxa_rank)
 legend_labels <- color_data$legend_labels
 labels <- color_data$labels

 # Extract unique values from vh_file$best_query
 unique_queries <- unique(vh_group[[groupby]])

 # Create a vector of corresponding names from legend_labels
 pyhlum_names <- legend_labels[unique_queries]

 # Match names to vh_group$best_query
 vh_group$phylum <- pyhlum_names[match(vh_group[[groupby]], unique_queries)]

 # Check for valid reorder_criteria
  valid_criteria <- c("max", "min")
  if (!is.null(reorder_criteria) && !(reorder_criteria %in% valid_criteria)) {
    stop("Invalid reorder_criteria. Please choose one of: max, min.")
  }




 ########################
 ### plot sum of hits ###
 ########################
 ########################



 sum_plot <- ggplot(data = vh_group,aes(x=if (!is.null(reorder_criteria)) {
   reorder(.data[[groupby]], if (reorder_criteria == "max") .data$sum else -.data$sum)
 } else {
     factor(.data[[groupby]], levels = rev(unique(sort(.data[[groupby]]))))
 },y= .data$sum,fill=.data$phylum))+
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

  sum_plot <- remove_group_text(sum_plot,remove_group_labels,flip_coords)

 if(groupby != "SRA_run"){

     sum_plot <- sum_plot + scale_fill_manual(values = labels)


 }








 #plot(sum_plot)
 message("Bar chart generation completed.")

 return(list(plot=sum_plot,vh_group=vh_group))



}



