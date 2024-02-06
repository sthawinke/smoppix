#' Spatial transcriptomics data of Selaginella moellendorffii
#'
#' Single-molecule spatial transcriptomics smFISH data of Selaginella moellendorffii roots of a
#' replicated experiment by \insertCite{Yang2023}{spatrans}. Molecule locations, gene identity and design
#' variables are included. Only a subset of the data, consisting of roots 1-3 and sections 1-5 is included for computational and memory reasons.
#'
#' @format A data matrix
#' \describe{
#'   \item{x,y}{Molecule coordinates}
#'   \item{gene}{Character vector with gene identities}
#'  \item{root,section,day}{Design variables}
#' }
#' @source \doi{10.1016/j.cub.2023.08.030}
#' @references
#' \insertAllCited{}
#' @usage data(Yang)
"Yang"
#' Spatial transcriptomics data of mouse fibroblast cells
#'
#' Single-molecule spatial transcriptomics seqFISH+ data containing measurements of 10,000 genes in NIH/3T3 mouse fibroblast cells by
#' \insertCite{Eng2019}{spatrans}. Molecule locations, gene identity and design
#' variables are included. In addition, a list of regions of interest (rois) is given describing the cell boundaries
#' @format A data matrix
#' \describe{
#'   \item{x,y}{Molecule coordinates}
#'   \item{gene}{Character vector with gene identities}
#'  \item{experiment,fov}{Design variables}
#' }
#' @source \doi{10.1016/j.cub.2023.08.030}
#' @references
#' \insertAllCited{}
#' @usage data(Eng)
"Eng"
