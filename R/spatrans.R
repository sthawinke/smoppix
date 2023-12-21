#' Spatrans: Analyze Spatial Transcriptomics and Other Spatial Omics Data on the Single Molecule Scale
#' @details
#' Test for univariate and multivariate spatial patterns in
#' spatial omics data with single-molecule resolution. The tests implemented
#' allow for analysis of nested designs and are automatically calibrated to different
#' biological specimens. Tests for aggregation, colocalization and
#' vicinity to cell edge are provided.
#' @useDynLib spatrans, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @author Stijn Hawinkel stijn.hawinkel@ugent.be
#' @name spatrans
#' @examples
#' data(Yang)
#' hypYang <- buildHyperFrame(Yang,
#'   coordVars = c('x', 'y'),
#'   imageVars = c('day', 'root', 'section')
#' )
#' # Fit a subset of features to limit computation time
#' yangPims <- estPims(hypYang[c(seq_len(5), seq(25, 29)),], pis = 'nn')
#' # Univariate nearest neighbour distances
NULL
