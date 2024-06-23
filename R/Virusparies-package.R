#' Virusparies Package
#'
#'@description
#' Virusparies contains a collection of Data Visualisations for VirusHunterGatherer output.
#'
#'@details
#' VirusHunterGatherer is a data-driven tool designed for high-throughput virus discovery.
#' The process involves Virushunter conducting sensitive homology-based detection of
#' viral sequence reads in unprocessed data. This approach identifies the most conserved regions
#' of a virus, which serve as seeds in the Virusgatherer step for a seed-based assembly
#' of full-length viral genome sequences.
#'
#' The Virusparies package provides a set of plotting functions tailored for visualizing
#' VirushunterGatherer Hittables. The name draws inspiration from the hunter-gatherer metaphor,
#' with "paries" derived from Latin meaning "wall". It symbolizes the parietal art left by
#' ancient hunters and gatherers on walls, summarizing their stories and beliefs.
#'
#' @section Functions:
#'
#'
#' This package includes the following functions:
#'
#' \strong{Import}:
#' \itemize{
#'   \item \code{\link{ImportVirusTable}}: Import VirusHunterGatherer hittables into R.
#' }
#'
#' \strong{VirushunterGatherer Plots}:
#' \itemize{
#'   \item \code{\link{VhgBoxplot}}: boxplot plotting RefSeq identity, evalues or contig length for each group.
#'   \item \code{\link{VhgIdenFacetedScatterPlot}}: Faceted Scatter Plot for Reference Identity vs -log10 of Reference E-value.
#'   \item \code{\link{VhgIdentityScatterPlot}}: Scatter Plot for Reference Identity vs -log10 of Reference E-value.
#'   \item \code{\link{VhgRunsBarplot}}: barplot showing how many unique Runs map against each virus.
#' }
#'
#' \strong{VirusGatherer only Plots}:
#' \itemize{
#'   \item \code{\link{VgConLenViolin}}: violin plot to visualize the distribution of contig lengths.
#' }
#'
#'
#' \strong{VirusHunter only Plots}:
#' \itemize{
#'   \item \code{\link{VhSumHitsBarplot}}: barplot for the sum of hits for each virus found in group.
#' }
#'
#'
#' \strong{Graphical Tables(gt)}:
#' \itemize{
#'   \item \code{\link{VhgRunsTable}}: Generate a gt-table for VhgRunsBarplot.
#'   \item \code{\link{VhgTabularRasa}}: Generate custome gt-tables.
#' }
#'
#'\strong{Export}:
#' \itemize{
#'   \item \code{\link{ExportVirusGt}}: Export Graphical Tables.
#'   \item \code{\link{ExportVirusPlot}}: Export plots.
#' }
#'
#' \strong{Utils}:
#' \itemize{
#'   \item \code{\link{VhgPreprocessTaxa}}: process ViralRefSeq_taxonomy column.
#' }
#'
#'
#' @name Virusparies
#' @aliases Virusparies-package
#'
#'
#'
#'
#' @seealso
#' VirusHunterGatherer is available on: \url{https://github.com/lauberlab/VirusHunterGatherer}.
#'
#' If you encounter any errors or issues, please report them at: \url{https://github.com/SergejRuff/Virusparies/issues}.
#'
#'
#'
#'@author
#'
#'
#' \strong{Maintainer}: Sergej Ruff (\email{serijnh@gmail.com})
#'
#' \strong{Other contributors}:
#'
#' \itemize{
#'   \item Chris Lauber (\email{chris.lauber@twincore.de})
#'   \item Li Chuin, Chong (\email{lichuin.chong@twincore.de})
#' }
#'
#' If you have any questions, suggestions, or issues, please feel free to contact the maintainer,
#' Sergej Ruff (\email{serijnh@gmail.com}).
#'
#'
"_PACKAGE"
