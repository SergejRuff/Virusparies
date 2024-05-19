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
#'
#' \strong{Virushunter Plots}:
#' \itemize{
#'   \item \code{\link{vhSumHitsBarplot}}: barplot for the sum of hits for each virus found in pHMM query.
#'   \item \code{\link{vhEvalIdenBoxplot}}: boxplot plotting identity or RefSeq evalues for each virus query.
#'   \item \code{\link{vhRunsBarplot}}: barplot showing how many unique Runs map agains each virus.
#' }
#'
#' \strong{Virusgatherer Plots}:
#' \itemize{
#'   \item \code{\link{vhSumHitsBarplot}}: Placeholder
#' }
#'
#' \strong{Tables}:
#' \itemize{
#'   \item \code{\link{vhRunsTable}}: Generate a gt-table for vhRunsBarplot.
#'   \item \code{\link{creatTableEvalBox}}: Generate a gt-table for vhEvalIdenBoxplot (eval_vs_iden="evalue") summary stats.
#'   \item \code{\link{creatTableIdentityBox}}: Generate a gt-table for vhEvalIdenBoxplot (eval_vs_iden="identity") summary stats.
#' }
#'
#'
#'
#' @name Virusparies
#' @aliases Virusparies-package
#'
#'
#'
#'
#'@seealso
#'VirusHunterGatherer is available on: \url{https://github.com/lauberlab/VirusHunterGatherer}.
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
