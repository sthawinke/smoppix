#' Helper function to spit gene pairs
#'
#' @param x character string
#' @param sep The character used to split
sund <- function(x, sep = "--") {
    strsplit(x, sep)[[1]]
}
#' Make design variable by combining different design variables
#'
#' @param x the design matrix
#' @param designVars the design variables to be combined
#' @param sep The string to separate the components
#'
#' @return a vector of design levels
makeDesignVar <- function(x, designVars, sep = "_") {
    apply(x[, designVars], 1, collapse = sep)
}
#' An aux function to build gene pairs
#'
#' @param genes The genes to be combined
#'
#' @return A chcracter vector of gene pairs
makePairs <- function(genes) {
    apply(combn(genes, 2), 2, paste, collapse = "--")
}
getN <- function(ecdf) {
    environment(ecdf)$nobs
}
#' Subsample a point pattern when it is too large
#'
#' @param p The point pattern
#' @param nSims The maximum size
#'
#' @return A point pattern, subsample if necessary
subSampleP <- function(p, nSims) {
    if ((NP <- npoints(p)) > nSims) p[sample(NP, nSims), ] else p
}
#' Extract elements from a distance matrix by index
#'
#' @param distMat the distance matrix of class 'dist'
#' @param i,j Indices for which entries are extracted
#' @details See documentation of the stats::dist function
getElDist <- function(distMat, i, j) {
    id <- outer(attr(distMat, "Size") * (i - 1) - i * (i - 1) / 2 - i, j, FUN = "+")
    # Following dist help
    distMat[c(id)]
}
#' Remove the diagonal from a square matrix
#' @details The diagonal elements are removed, and the elements right from the
#' diagonal are shifted left to form a new matrix with one column less
#'
#' @param x The matrix from which to remove the diagonal
#'
#' @return A matrix with the same number of rows as x and one column less
dropDiagonal <- function(x) {
    diagId <- seq(1, length(x), by = (NR <- nrow(x)) + 1)
    t(matrix(t(x)[-diagId], ncol = NR))
}
