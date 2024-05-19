#' @title vhIdentityScatterPlot: Scatter Plot for Reference Identity vs -log10 of Reference E-value
#'
#' @description
#' This function creates a scatter plot of viral RefSeq identity versus
#' the negative logarithm (base 10) of viral RefSeq E-value. It colors the points
#' based on  best query  and adds a horizontal line representing the cutoff value.
#'
#' @param vh_file A data frame containing VirusHunter Hittable results.
#' @param cutoff The significance cutoff value for E-values (default: 1e-5). Removes rows in vh_file
#' with values larger than cutoff value in ViralRefSeq_E column.
#'
#' @return A ggplot object representing the scatter plot.
#'
#'
#' @import ggplot2
#' @importFrom rlang .data
#' @export
vhIdentityScatterPlot <- function(vh_file,cutoff = 1e-5){

  cutoff <- -log10(cutoff)

  # Calculate the maximum y-value
  max_y <- max(-log10(vh_file$ViralRefSeq_E)) + 5

  # Plot the data and color points based on the cutoff condition
  iden_refevalue <- ggplot(vh_file, aes(x = .data$ViralRefSeq_ident, y = -log10(.data$ViralRefSeq_E))) +
    geom_point(aes(color=.data$best_query))+ geom_hline(aes(yintercept=cutoff), colour="#990000")+
    labs(x="viral Refseq Identity in %",
         y="-log10 of viral Reference e-value",
         title="scatterplot for refrence identity vs -log10 of refrence e-value")+
    theme_linedraw()+
    theme(legend.position = "bottom")+
    guides(fill=guide_legend(title="virus family"))+
    theme(
      # This is the new default font in the plot
      text = element_text( size = 8, color = "black"),
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

  plot(iden_refevalue)

  return(iden_refevalue)
}
