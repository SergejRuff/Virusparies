#' Create a Scatterplot of Viral RefSeq Identity vs. -log10 of Viral Reference e-value faceted by
#' the `best_query` column.
#'
#' This function generates a scatterplot of viral RefSeq identity versus -log10 of viral reference e-value
#' for each virus group in the `best_query` column. The points are colored based on whether the
#' e-value meets a specified cutoff and are faceted by the `best_query` column.
#'
#' @param vh_file A data frame containing at least the columns `ViralRefSeq_E`, `ViralRefSeq_ident`, and `best_query`.
#' @param cutoff A numeric value representing the cutoff for the viral reference e-value. Points with `ViralRefSeq_E`
#' less than or equal to this value will be colored blue; otherwise, they will be colored red (default: 1e-5).
#'
#' @return A ggplot object representing the scatterplot.
#'
#' @details This function takes a data frame and a cutoff value as inputs, adds a new column to the
#' data frame indicating whether each `ViralRefSeq_E` value meets the cutoff, and then plots the data.
#' The plot includes:
#' - Points colored based on whether they meet the cutoff condition.
#' - Faceting by the `best_query` column.
#'
#'
#' @import ggplot2
#' @importFrom stats runif
#' @importFrom rlang .data
#' @export
vhIdenFacetedScatterPlot <- function(vh_file,cutoff = 1e-5){

  vh_file$cutoff_met <- vh_file$ViralRefSeq_E <= cutoff

  max_y <- max(-log10(vh_file$ViralRefSeq_E)) + 5

  # Plot the data and color points based on the cutoff condition, faceted by 'best_query'
  iden_refevalue_seperate <- ggplot(vh_file, aes(x = .data$ViralRefSeq_ident, y = -log10(.data$ViralRefSeq_E))) +
    geom_point(aes(color = .data$cutoff_met)) +  # Map color to the cutoff condition
    scale_color_manual(values = c("TRUE" = "blue", "FALSE" = "red")) +  # Define colors for TRUE and FALSE
    facet_wrap(~.data$best_query, ncol = 2)+
    labs(x="viral Refseq Identity in %",
         y="-log10 of viral Reference e-value",
         title="scatterplot for refrence identity vs -log10 of refrence e-value for each virus seperatly")+
    theme_linedraw()+
    theme(legend.position = "bottom")+
    guides(fill=guide_legend(title="virus family"))+
    theme(
      # This is the new default font in the plot
      text = element_text(size = 8, color = "black"),
      plot.title = element_text(
        size = 16,
        face = "bold",
        color = "#2a475e"),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      axis.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      legend.key.size = unit(1.5, "lines"),
      legend.title = element_text(size = 12, face = "bold")
    )+ scale_y_continuous(expand = c(0, 0), limits = c(0, max_y))+
    scale_x_continuous(
      expand = c(0, 0),
      limits = c(0, 105),
      breaks = seq(0, 100, by = 10)
    )

  plot(iden_refevalue_seperate)

  return(iden_refevalue_seperate)
}
