#' Spatial transcriptomics data of Selaginella
#'
#' Single-molecule spatial transcriptomics smFISH data of Selaginella moellendorffii roots of a
#' replicated experiment by \insertCite{Yang2023}{spatrans}. Molecule locations, gene identity and design
#' variables are included. Only a subset of the data, consisting of sections 1-3 of roots 1-3 is included for computational reasons.
#'
#' @format A data matrix
#' \describe{
#'   \item{x,y}{Molecule coordinates}
#'   \item{gene}{Character vector with gene identities}
#'   \item{root,section,day}{Design variables}
#' }
#' @source \url{https://doi.org/10.1016/j.cub.2023.08.030}
#' @references
#' \insertAllCited{}
"Yang"