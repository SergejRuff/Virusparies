#' Virusparies Package
#'
#'@description
#'A comprehensive toolkit for repeated high-dimensional analysis.
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
#'   \item \code{\link{vh_sum_hits_bar}}: barplot for the sum of hits for each virus found in pHMM query.
#'   \item \code{\link{vh_refeval_iden_boxplot}}: boxplot plotting identity or RefSeq evalues for each virus query.
#'   \item \code{\link{vh_runs_bar}}: barplot showing how many unique Runs map agains each virus
#' }
#'
#' \strong{Virusgatherer Plots}:
#' \itemize{
#'   \item \code{\link{vh_sum_hits_bar}}: Placeholder
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
#'   \item Chong, Li Chuin (\email{lichuin.chong@twincore.de})
#' }
#'
#' If you have any questions, suggestions, or issues, please feel free to contact the maintainer,
#' Sergej Ruff (\email{serijnh@gmail.com}).
#'
#'
"_PACKAGE"
