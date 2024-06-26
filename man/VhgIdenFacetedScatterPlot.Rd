% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VhgIdenFacetedScatterPlot.R
\name{VhgIdenFacetedScatterPlot}
\alias{VhgIdenFacetedScatterPlot}
\title{VhgIdenFacetedScatterPlot: Create a Scatter plot of Viral RefSeq Identity vs. -log10 of Viral Reference e-value.}
\usage{
VhgIdenFacetedScatterPlot(
  file,
  groupby = "best_query",
  taxa_rank = "Family",
  cutoff = 1e-05,
  conlen_bubble_plot = FALSE,
  contiglen_breaks = 5,
  theme_choice = "linedraw",
  title = "Faceted scatterplot of viral reference e-values and identity",
  title_size = 16,
  title_face = "bold",
  title_colour = "#2a475e",
  subtitle = NULL,
  subtitle_size = 12,
  subtitle_face = "bold",
  subtitle_colour = "#1b2838",
  xlabel = "Viral reference identity (\%)",
  ylabel = "-log10 of viral reference e-values",
  axis_title_size = 12,
  xtext_size = 10,
  ytext_size = 10,
  legend_position = "bottom",
  legend_title_size = 12,
  legend_title_face = "bold",
  legend_text_size = 10,
  true_colour = "blue",
  false_colour = "red",
  wrap_ncol = 2
)
}
\arguments{
\item{file}{A data frame containing at least the columns \code{ViralRefSeq_E}, \code{ViralRefSeq_ident}, and the specified grouping column.}

\item{groupby}{(optional) A string indicating the column used for grouping the data points in the scatter plot.
The values in this column determine the color of each point in the scatter plot. Default is "best_query".
Note: Gatherer hittables do not have a "best_query" column. Please provide an appropriate column for grouping.}

\item{taxa_rank}{(optional) When \code{groupby} is set to "ViralRefSeq_taxonomy", specify the taxonomic rank to group your data by.
Supported ranks are:
\itemize{
\item "Subphylum"
\item "Class"
\item "Subclass"
\item "Order"
\item "Suborder"
\item "Family" (default)
\item "Subfamily"
\item "Genus" (including Subgenus)
}}

\item{cutoff}{(optional) A numeric value representing the cutoff for the viral reference e-value. Points with \code{ViralRefSeq_E}
less than or equal to this value will be colored blue; otherwise, they will be colored red (default: 1e-5).}

\item{conlen_bubble_plot}{(optional) Logical value indicating whether the \code{contig_len} column
should be used to size the bubbles in the plot. Applicable only to VirusGatherer hittables input. Default is FALSE.}

\item{contiglen_breaks}{(optional) number of breaks for the bubble plot (for \code{conlen_bubble_plot}=TRUE). Default is 5.}

\item{theme_choice}{(optional) A character indicating the ggplot2 theme to apply. Options include "minimal",
"classic", "light", "dark", "void", "grey" (or "gray"), "bw", "linedraw", and "test".
Default is "linedraw".}

\item{title}{(optional) The title of the plot. Default is "Faceted scatter plot of viral reference e-values and identity".}

\item{title_size}{(optional) The size of the title text. Default is 16.}

\item{title_face}{(optional) The face (bold, italic, etc.) of the title text. Default is "bold".}

\item{title_colour}{(optional) The color of the title text. Default is "#2a475e".}

\item{subtitle}{(optional) The subtitle of the plot. Default is NULL.}

\item{subtitle_size}{(optional) The size of the subtitle text. Default is 12.}

\item{subtitle_face}{(optional) The face (bold, italic, etc.) of the subtitle text. Default is "bold".}

\item{subtitle_colour}{(optional) The color of the subtitle text. Default is "#1b2838".}

\item{xlabel}{(optional) The label for the x-axis. Default is "Viral reference identity (\%)".}

\item{ylabel}{(optional) The label for the y-axis. Default is "-log10 of viral reference e-values".}

\item{axis_title_size}{(optional) The size of the axis titles. Default is 12.}

\item{xtext_size}{(optional) The size of the x-axis text. Default is 10.}

\item{ytext_size}{(optional) The size of the y-axis text. Default is 10.}

\item{legend_position}{(optional) The position of the legend. Default is "bottom".}

\item{legend_title_size}{(optional) The size of the legend title text. Default is 12.}

\item{legend_title_face}{(optional) The face (bold, italic, etc.) of the legend title text. Default is "bold".}

\item{legend_text_size}{(optional) The size of the legend text. Default is 10.}

\item{true_colour}{(optional) The color for points that meet the cutoff condition. Default is "blue".}

\item{false_colour}{(optional) The color for points that do not meet the cutoff condition. Default is "red".}

\item{wrap_ncol}{(optional) The number of columns for faceting. Default is 2.}
}
\value{
A list containing the following components:
\itemize{
\item plot: A plot object representing the faceted scatter plot.
\item evalue_stats: A tibble data frame with summary statistics for "ViralRefSeq_E" values.
\item identity_stats: A tibble data frame with summary statistics for "ViralRefSeq_ident" values.
\item contig_stats (optional): A tibble data frame with summary statistics for "contig_len" values, included only if VirusGatherer is used with \code{conlen_bubble_plot=TRUE}.
}
}
\description{
This function generates a scatter plot of viral RefSeq identity versus -log10 of viral reference e-value
for each virus group in the \code{best_query} column (or another column). The points are colored based on whether the
e-value meets a specified cutoff and are faceted by the \code{best_query} column (or another column).
}
\details{
'VhgIdenFacetedScatterPlot' takes a VirusHunter or VirusGatherer Hittable and a cutoff value as inputs.
The plot includes:
\itemize{
\item Points colored based on whether they meet the cutoff condition.
\item Faceting by the \code{best_query} column as the default column. The user can provide their own column
for grouping.
\item The option \code{conlen_bubble_plot} = TRUE generates a bubble plot where the size of points corresponds to "contig_len" (exclusive to VirusGatherer).
}

Tibble data frames containing summary statistics (median, Q1, Q3, mean, sd, min, and max) for 'ViralRefSeq_E'
and 'ViralRefSeq_ident' values are generated. Optionally, summary statistics for 'contig_len' values
are also included if applicable. These summary statistics, along with the plot object, are returned within a list object.
}
\examples{
path <- system.file("extdata", "virushunter.tsv", package = "Virusparies")
file <- ImportVirusTable(path)

# plot 1
plot <- VhgIdenFacetedScatterPlot(file,cutoff = 1e-5)

plot

# plot 2 with custom data
custom_plot <- VhgIdenFacetedScatterPlot(file,
                                         cutoff = 1e-4,
                                         theme_choice = "dark",
                                         title = "Custom Scatterplot",
                                         title_size = 18,
                                         title_face = "italic",
                                         title_colour = "orange",
                                         xlabel = "Custom X Label",
                                         ylabel = "Custom Y Label",
                                         axis_title_size = 14,
                                         legend_position = "right",
                                         true_colour = "green",
                                         false_colour = "purple")

custom_plot

# import gatherer files
path2 <- system.file("extdata", "virusgatherer.tsv", package = "Virusparies")
vg_file <- ImportVirusTable(path2)

# vgplot: virusgatherer plot with ViralRefSeq_taxonomy as custom grouping
vgplot <- VhgIdenFacetedScatterPlot(vg_file,groupby = "ViralRefSeq_taxonomy")
vgplot


}
\seealso{
VirusHunterGatherer is available here: \url{https://github.com/lauberlab/VirusHunterGatherer}.
}
\author{
Sergej Ruff
}
