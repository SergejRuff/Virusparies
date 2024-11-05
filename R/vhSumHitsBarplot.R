#' @title VhgSumHitsBarplot: Generate a bar plot showing the sum of reads/contigs for each virus family
#'
#' @description
#' VhgSumHitsBarplot preprocesses virus data for plotting and generates a bar plot showing the sum of reads/contigs
#' for each virus family from the input data set.
#'
#' @param file A data frame containing VirusHunters hittable results.
#'
#' @param groupby (optional): A character specifying the column containing the groups (default: "best_query").
#'
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
#' @param y_column A character specifying the column containing the values to be compared. Currently "ViralRefSeq_contigs" (micro-contigs),"contigs",
#'  and "num_hits" (reads) are supported columns (default:"num_hits").
#' @param cut (optional): A numeric value representing the cutoff for the refseq E-value (default: 1e-5). Removes rows in file
#' with values larger than cutoff value in "ViralRefSeq_E" column.
#'
#' @param reorder_criteria (optional): Character string specifying the criteria for reordering the x-axis ('max' (default), 'min','phylum',phylum_max,phylum_min).
#' NULL sorts alphabetically.
#'
#' @param theme_choice (optional): A character indicating the ggplot2 theme to apply. Options include "minimal",
#'  "classic", "light", "dark", "void", "grey" (or "gray"), "bw", "linedraw" (default), and "test".
#'  Append "_dotted" to any theme to add custom dotted grid lines (e.g., "classic_dotted").
#'
#' @param flip_coords (optional): Logical indicating whether to flip the coordinates of the plot (default: TRUE).
#'
#' @param title (optional): The title of the plot (default: "Distribution of hits for each virus group").
#'
#' @param title_size (optional): The size of the title text (default: 16).
#' @param title_face (optional): The face (bold, italic, etc.) of the title text (default: "bold").
#' @param title_colour (optional): The color of the title text (default: "#2a475e").
#'
#' @param subtitle (optional): A character specifying the subtitle of the plot.
#' Default is "total number of hits/micro-contigs: " followed by the calculated number.empty string ("") removes subtitle.
#'
#' @param subtitle_size (optional): Numeric specifying the size of the subtitle text(default: 12).
#' @param subtitle_face (optional): A character specifying the font face for the subtitle text (default: "bold").
#' @param subtitle_colour (optional): A character specifying the color for the subtitle text (default: "#1b2838").
#' .
#' @param xlabel (optional): The label for the x-axis (default: "Viral group").
#' @param ylabel (optional): The label for the y-axis (default: "Total number of hits").
#'
#' @param axis_title_size (optional): The size of the axis titles (default: 12).
#' @param xtext_size (optional): The size of the x-axis text (default: 10).
#' @param x_angle (optional): An integer specifying the angle (in degrees) for the x-axis text labels. Default is NULL, meaning no change.
#' @param ytext_size (optional): The size of the y-axis text (default: 10).
#' @param y_angle (optional): An integer specifying the angle (in degrees) for the y-axis text labels. Default is NULL, meaning no change.
#' @param remove_group_labels (optional): If `TRUE`, the group labels will be removed; if `FALSE` or omitted, the labels will be displayed.
#' @param legend_title (optional): A character specifying the title for the legend (default: "Phylum").
#' @param legend_position (optional): A character specifying the position of the legend (default: "bottom").
#' @param legend_title_size (optional): Numeric specifying the size of the legend title text (default: 12).
#' @param legend_title_face (optional): A character specifying the font face for the legend title text (default: "bold").
#' @param legend_text_size (optional): Numeric specifying the size of the legend text (default: 10).
#' @param plot_text (optional): An index (0-3) to select the variable for text labels.
#' - 0 = None.
#' - 1 = Number of hits for each viral group.
#' - 2 = Only the percentage.
#' - 3 = Both (Default).
#' @param plot_text_size (optional): The size of the text labels added to the plot (default: 3.5).
#' @param plot_text_position_dodge (optional): The degree of dodging for positioning text labels (default: 0.9).
#' @param plot_text_hjust (optional): The horizontal justification of text labels (default: -0.1).
#' @param plot_text_vjust (optional): The vertical justification of text labels (default: 0.5).
#' It is recommended to change `vjust` when setting `flip_coords = FALSE`.
#' @param plot_text_colour (optional): The color of the text labels added to the plot (default: "black").
#' @param facet_ncol (optional):  The number of columns for faceting (default: NULL).
#' It is recommended to specify this when the number of viral groups is high, to ensure they fit well in one plot.
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
#'
#' @details
#' VhgSumHitsBarplot preprocesses virus data for plotting by calculating the sum of
#' hits for each virus family from the input data set (accepts only VirusHunter hittables).
#' It then generates a bar plot showing the sum of hits for each virus family.
#' Additionally, it returns the processed data for further analysis.
#'
#' @return A list containing the generated bar plot and processed data.
#'
#' @author Sergej Ruff
#' @examples
#' path <- system.file("extdata", "virushunter.tsv", package = "Virusparies")
#' file <- ImportVirusTable(path)
#'
#' # plot 1: plot boxplot for reads
#' plot <- VhgSumHitsBarplot(file,cut = 1e-5)
#' plot
#'
#' # plot 2: plot boxplot for micro_reads
#' plot_reads <- VhgSumHitsBarplot(file,cut = 1e-5,
#' y_column = "ViralRefSeq_contigs")
#' plot_reads
#'
#' # import gatherer files
#' path2 <- system.file("extdata", "virusgatherer.tsv", package = "Virusparies")
#' vg_file <- ImportVirusTable(path2)
#'
#'
#' # plot 3: contigs in Gatherer
#' contig_plot <- VhgSumHitsBarplot(vg_file,groupby = "ViralRefSeq_taxonomy",
#' y_column = "contig")
#' contig_plot
#'
#' @seealso
#' VirusHunterGatherer is available here: \url{https://github.com/lauberlab/VirusHunterGatherer}.
#'
#' @import ggplot2
#' @importFrom rlang .data as_string ensym
#' @importFrom stats reorder
#' @export
VhgSumHitsBarplot <- function(file,
                             groupby = "best_query",
                             taxa_rank = "Family",
                             y_column = "num_hits",
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
                             x_angle = NULL,
                             ytext_size = 10,
                             y_angle = NULL,
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
                             plot_text_colour = "black",
                             facet_ncol = NULL,
                             group_unwanted_phyla =NULL){


 #check if table is empty
 #is_file_empty(file)
 if (is_file_empty(file)) {
    #message("Skipping VhgBoxplot generation due to empty data.")
    return(invisible(NULL))  # Return invisible(NULL) to stop further execution
 }

 # Convert non-quoted column names to strings
 groupby <- rlang::as_string(rlang::ensym(groupby))
 taxa_rank <- rlang::as_string(rlang::ensym(taxa_rank))
 y_column <- rlang::as_string(rlang::ensym(y_column))




 # Define the required columns
 required_columns <- c("ViralRefSeq_E", groupby)

 if (!(groupby %in% c("best_query", "ViralRefSeq_taxonomy"))) {
     stop('Invalid value for groupby. Please use either "best_query" or "ViralRefSeq_taxonomy".')
 }

 check_columns(file,required_columns)
 check_input_type(file,c("ViralRefSeq_E"),2)
 check_input_type(file,groupby,1)


 # check arguments
 arg_character(theme_choice)
 arg_character(legend_position)

 arg_logical(flip_coords)


 if(groupby == "ViralRefSeq_taxonomy"){

     file <- VhgPreprocessTaxa(file,taxa_rank)

 }


 ## preprocess data for plotting
 file <- file[file$ViralRefSeq_E < cut,]
 #is_file_empty(file)
 if (is_file_empty(file)) {
   #message("Skipping VhgBoxplot generation due to empty data.")
   return(invisible(NULL))  # Return invisible(NULL) to stop further execution
 }


 vh_group <- vh_sumhitbar_preprocessing(file,groupby,y_column)

 # Apply the selected theme
 theme_selected <- select_theme(theme_choice)


 if (y_column == "ViralRefSeq_contigs") {
   micro_contig_or_reads <- "micro-contigs"
 } else if ("contig_len" %in% colnames(file) && y_column == "contig") {
   micro_contig_or_reads <- "contigs"
 } else if (y_column == "num_hits") {
   micro_contig_or_reads <- "reads"
 } else {
   stop("Error: y_column must be 'ViralRefSeq_contigs', 'contig', or 'num_hits'.")
 }



 # Set the subtitle based on the input
 if (!is.null(subtitle) && subtitle == "default") {
   subtitle_text <- paste0("Total number of ",micro_contig_or_reads,": ", sum(vh_group$sum))
 } else if (is.null(subtitle) || subtitle == "") {
   subtitle_text <- NULL
 } else {
   subtitle_text <- subtitle
 }


 color_data <- consistentColourPalette(file, groupby = groupby,taxa_rank=taxa_rank)
 legend_labels <- color_data$legend_labels
 labels <- color_data$labels

 # Extract unique values from file$best_query
 unique_queries <- unique(vh_group[[groupby]])

 # Create a vector of corresponding names from legend_labels
 pyhlum_names <- legend_labels[unique_queries]

 # Match names to vh_group$best_query
 vh_group$phyl <- pyhlum_names[match(vh_group[[groupby]], unique_queries)]

 if(!is.null(group_unwanted_phyla)){

   group_unphyla <- remove_non_group(file=vh_group,groupby=groupby,chosen_group=group_unwanted_phyla,label_vector=labels,
                                     taxa_rank=taxa_rank)

   vh_group <- group_unphyla$file

   labels <- group_unphyla$label

   best_query_col <- if (groupby == "ViralRefSeq_taxonomy") {
     "ViralRefSeq_taxonomy"
   } else {
     "best_query"
   }

   phyl <- "phyl"
   perc <- "perc"
   cyl <- "cyl"
   res <- "res"

   merge_phyl_values <- c("Non-RNA-viruseses", "Non-Small-DNA-Viruseses", "Non-Large-DNA-Viruseses","Other Viruses")

   vh_group <- vh_group %>%
     # Separate the data into those that should be aggregated and those that should not
     mutate(merge_flag = if_else(!!sym(phyl) %in% merge_phyl_values, "Merge", "Keep")) %>%

     # Aggregate the rows where phyl matches specified values
     group_by(!!sym(best_query_col), phyl, .data$merge_flag) %>%
     summarize(
       !!sym(best_query_col) := paste(unique(!!sym(best_query_col)), collapse = ", "),
       sum = sum(sum),  # Aggregate the 'sum' column
       perc = paste0(
         round(
           sum(as.numeric(gsub("%", "", !!sym(perc))) * sum) / sum(sum),
           2
         )
       ),
       res = paste(sum, " (", perc, "%", ")"),
       cyl = paste(unique(!!sym(cyl)), collapse = ", "),
       .groups = 'drop'
     ) %>%
     ungroup()

   # Process non-merged data separately
   vh_group_non_merge <- vh_group %>%
     filter(.data$merge_flag == "Keep") %>%
     select(-.data$merge_flag)  # Remove the merge_flag after processing

   # Combine the merged and non-merged data
   vh_group <- bind_rows(
     vh_group %>%
       filter(.data$merge_flag == "Merge") %>%
       select(-.data$merge_flag),  # Remove the merge_flag after processing
     vh_group_non_merge
   )

   # Select final columns
   vh_group <- vh_group %>%
     select(!!sym(best_query_col), sum, perc, res, cyl, phyl)



 }

 # Check for valid reorder_criteria
 valid_criteria <- c("max", "min","phylum","phylum_max","phylum_min")
 if (!is.null(reorder_criteria) && !(reorder_criteria %in% valid_criteria)) {
   stop("Invalid reorder_criteria. Please choose one of: max, min,phylum,phylum_max,phylum_min.")
 }

  text_var_names <- c("none", "sum", "perc", "res")

  # Check if the provided index is valid
  if (plot_text < 0 || plot_text >= length(text_var_names)) {
    stop("Invalid plot_text. Choose an index from 0 to 3.")
  }

  text_var <- text_var_names[plot_text + 1]

  vh_group <- vh_group %>% mutate(perc=paste(.data$perc,"%"))





 ########################
 ### plot sum of hits ###
 ########################
 ########################





 sum_plot <-   ggplot(data = vh_group, aes(
   x = if (!is.null(reorder_criteria)) {
     if (reorder_criteria == "phylum") {
       # Sort `groupby` levels alphabetically by `phyl`, and reverse the order
       factor(.data[[groupby]], levels = rev(unique(.data[[groupby]][order(.data$phyl)])))
     }else if (grepl("^phylum_", reorder_criteria)) {
       # Extract the secondary criterion (e.g., "median" from "phylum_median")
       secondary_criteria <- sub("^phylum_", "", reorder_criteria)

       # Use aggregate to calculate the secondary criteria within each phylum
       agg_data <- aggregate(.data$sum, list(phylum = .data$phyl, x_val = .data[[groupby]]),
                             FUN = switch(secondary_criteria,
                                          "max" = max,
                                          "min" = min))




       agg_data$phylum <- factor(agg_data$phylum, levels = rev(sort(unique(agg_data$phylum))))


       agg_data <- agg_data[order(agg_data$phylum, agg_data$x), ]



       ordered_levels <- agg_data$x_val[order(agg_data$phylum, if (secondary_criteria == "max") agg_data$x else -agg_data$x)]
       factor(.data[[groupby]], levels = ordered_levels)
     } else {
       # Reorder `groupby` based on `unique_SRA_run` using the specified criteria
       reorder(.data[[groupby]], if (reorder_criteria == "max") .data$sum else -.data$sum)
     }
   } else {
     # Default ordering if no `reorder_criteria` is specified
     factor(.data[[groupby]], levels = rev(unique(sort(.data[[groupby]]))))
   },
   y = .data$sum,
   fill = .data$phyl
 ))+
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
   guides(fill=guide_legend(title=legend_title))+
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
     ))+
    theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

  if (text_var != "none") {


    sum_plot <- sum_plot + geom_text(aes(label = ifelse(sum == 0, "",.data[[text_var]])),
                                     hjust = plot_text_hjust,
                                     vjust = plot_text_vjust,
                                     color = plot_text_colour,
                                     position = position_dodge(plot_text_position_dodge),
                                     size = plot_text_size)


  }






 if (flip_coords) {
    sum_plot <- sum_plot + coord_flip()
 }

  sum_plot <- remove_group_text(sum_plot,remove_group_labels,flip_coords)

 if(groupby != "SRA_run"){

     sum_plot <- sum_plot + scale_fill_manual(values = labels)


 }

  sum_plot <- facet_plot(sum_plot,facet_ncol,flip_coords)

  sum_plot <- adjust_plot_angles(sum_plot,x_angle = x_angle,y_angle = y_angle)








 #plot(sum_plot)
 message("Bar chart generation completed.")

 return(list(plot=sum_plot,vh_group=vh_group))



}



