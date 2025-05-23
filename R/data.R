#' Spatial transcriptomics data of Selaginella moellendorffii roots
#'
#' Single-molecule spatial transcriptomics smFISH data of Selaginella moellendorffii roots of a
#' replicated experiment by \insertCite{Yang2023}{smoppix}. Molecule locations, gene identity and design
#' variables are included. Only a subset of the data, consisting of roots 1-3 and sections 1-5 is included in the package for computational and memory reasons.
#' The data are in table format to illustrate conversion to \link[spatstat.geom]{hyperframe} using \link{buildHyperFrame}.
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
#' \insertCite{Eng2019}{smoppix}. Molecule locations, gene identity and design
#' variables are included, a subset of eight most expressed genes is included in the package,
#' and the dataset was subsampled to 100,000 observations for memory reasons.
#' In addition, a list of regions of interest (rois) is given describing the cell boundaries.
#' @format
#' \enumerate{
#' \item \strong{Eng} A data frame with variables
#' \describe{
#'   \item{x,y}{Molecule coordinates}
#'   \item{gene}{Character vector with gene identities}
#'  \item{experiment,fov}{Design variables}
#' }
#' \item \strong{EngRois} A list of lists of regions of interest (ROIs): the cell boundaries
#' }
#' @source \doi{10.1038/s41586-019-1049-y}
#' @references
#' \insertAllCited{}
#' @usage data(Eng)
#' @aliases Eng EngRois
"Eng"
