#' Generate a bar plot showing the sum of hits for each virus family
#'
#' This function preprocesses virus data for plotting and generates a bar plot showing the sum of hits
#' for each virus family from the input dataset.
#'
#' @param vh_file A data frame containing VirusHunters hittables results.
#' @param cut The significance cutoff value for E-values (default: 1e-5). Removes rows in vh_file
#' with values larger than cutoff value in ViralRefSeq_E column.
#'
#' @return A list containing the generated bar plot and processed data
#'
#' @details This function preprocesses virus data for plotting by calculating the sum of
#' hits for each virus family from the input dataset. It then generates a bar plot showing
#' the sum of hits for each virus family. Additionally, it returns the processed data
#' for further analysis.
#'
#' @examples
#' path <- system.file("extdata", "virushunter.tsv", package = "Virusparies")
#' vh_file <- importVirusTable(path)
#'
#' # plot 1: plot boxplot for "identity"
#' plot <- vhSumHitsBarplot(vh_file,cut = 1e-5)
#'
#' # return vh_group inside plot object
#' print(plot$vh_group)
#'
#' @import ggplot2
#' @importFrom rlang .data
#'
#' @export
vhSumHitsBarplot <- function(vh_file,cut = 1e-5){



 ## preprocess data for plotting
 vh_file <- vh_file[vh_file$ViralRefSeq_E < cut,]
 vh_group <- vh_sumhitbar_preprocessing(vh_file)

 ## need to genrate a table function that will also be returned.


 ########################
 ### plot sum of hits ###
 ########################
 ########################

 sum_plot <- ggplot(data = vh_group,aes(x=reorder(.data$best_query,.data$sum,FUN=max),y= .data$sum,fill=.data$best_query))+
   geom_bar(stat = "identity")+
   labs(title = "sum of hits for each family",
        x="Virus family found in query",
        y="sum of hits",
        subtitle = paste0("total number of hits: ",sum(vh_group$sum)))+
   theme_minimal()+
   theme(legend.position = "bottom",
         axis.text.y = element_text(size = 10),
         axis.text.x = element_text(size = 10),
         axis.title = element_text(size = 12),
         legend.text = element_text(size = 10),
         legend.key.size = unit(1.5, "lines"),
         legend.title = element_text(size = 12, face = "bold"))+
   guides(fill=guide_legend(title="virus family")) +
   geom_text(aes(label = ifelse(sum == 0, "",.data$res)),
             hjust = -0.1,
             color = "black",
             position = position_dodge(0.9),
             size = 3.5)+
   coord_flip()+
   scale_y_continuous(expand = c(0, 0), limits = c(0, max(vh_group$sum)+max(vh_group$sum)/8))+
   theme(
     # This is the new default font in the plot
     text = element_text( size = 8, color = "black"),
     plot.title = element_text(
       size = 16,
       face = "bold",
       color = "#2a475e"),
     # Statistical annotations below the main title
     plot.subtitle = element_text(
       size = 12,
       face = "bold",
       color="#1b2838"
     ))


 plot(sum_plot)

 return(list(plot=sum_plot,vh_group=vh_group))

}



