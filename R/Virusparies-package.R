#' @title Virusparies Package
#'
#' @description
#' Virusparies designed for visualizing output from VirusHunterGatherer.
#'
#' @details
#' VirusHunterGatherer is a pipeline designed for data-driven virus discovery (DDVD).
#' It involves two steps: (i) Virushunter conducts sensitive homology-based detection of viral sequence reads in unprocessed data,
#' identifying the most conserved regions of a virus, which serve as seeds for the
#' (ii) Virusgatherer step that assembles full-length viral genome sequences.
#'
#' The Virusparies package provides a set of plotting functions tailored for visualizing
#' VirushunterGatherer hittables. The name draws inspiration from the hunter-gatherer metaphor,
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
#' \strong{VirusHunterGatherer Plots}:
#' \itemize{
#'   \item \code{\link{VhgBoxplot}}: Box plot plotting refseq sequence identity, E-values or contig length for each group.
#'   \item \code{\link{VhgIdenFacetedScatterPlot}}: Faceted scatter plot for reference sequence identity vs -log10 of reference E-value.
#'   \item \code{\link{VhgIdentityScatterPlot}}: Scatter plot for reference sequence identity vs -log10 of reference E-value.
#'   \item \code{\link{VhgRunsBarplot}}: Bar plot showing how many unique runs map against each virus.
#'   \item \code{\link{VhgSumHitsBarplot}}: Bar plot for the sum of hits for each virus found in group.
#' }
#'
#'
#' \strong{VirusGatherer only plots}:
#' \itemize{
#'   \item \code{\link{VgConLenViolin}}: Violin plot to visualize the distribution of contig lengths.
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
#'   \item \code{\link{ExportVirusDataFrame}}: Export data frames.
#'   \item \code{\link{ExportVirusGt}}: Export graphical tables.
#'   \item \code{\link{ExportVirusPlot}}: Export plots.
#' }
#'
#' \strong{Utils}:
#' \itemize{
#'   \item \code{\link{CombineHittables}}: Combine hittables.
#'   \item \code{\link{Current_ICTV}}: Verify the version of ICTV data currently in use.
#'   \item \code{\link{New_ICTV}}: Assign Custom ICTV Data for Use in Virusparies.
#'   \item \code{\link{SummarizeViralStats}}: Generate summary stats outside of plot functions.
#'   \item \code{\link{VhgAddPhylum}}: Extract phylum information.
#'   \item \code{\link{VhgGetSubject}}: Process and Count Viral Subjects within Groups.
#'   \item \code{\link{VhgPreprocessTaxa}}: Process ViralRefSeq_taxonomy column.
#'   \item \code{\link{VhgSubsetHittable}}: Filter VirusHunterGatherer data based on user’s own criteria.
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
#' If you have any questions, suggestions, or issues, please feel free to contact Sergej Ruff (\email{serijnh@gmail.com}).
#'
#'
"_PACKAGE"
