#' Helper function to spit gene pairs
#'
#' @param x character string
#' @param sep The character used to split
#' @return The split string
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
#' Return the number of observations used to fit the ecdf
#'
#' @param ecdf The ecdf function
#'
#' @return The number of observations
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
#' @return The required distances
getElDist <- function(distMat, i, j) {
    id <- outer(attr(distMat, "Size") * (i-1) - i*(i-1)/2 - i, j, FUN = "+")
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
#' A version of contr.sum that retains names, a bit controversial but also
#' clearer
#'
#' @param x,... passed on to contr.sum
#' @importFrom stats contr.sum
#' @return The matrix of contrasts
named.contr.sum <- function(x, ...) {
    lev <- x
    x <- contr.sum(x, ...)
    colnames(x) <- lev[-length(lev)]
    x
}
# After https://stackoverflow.com/questions/24515892/r-how-to-contrast-code-factors-and-retain-meaningful-labels-in-output-summary
#' Add tables with gene counts to the hyperframe
#'
#' @param hypFrame The hyperframe
#'
#' @return The hyperframe with tabObs added
addTabObs = function(hypFrame){
    hypFrame$tabObs = lapply(hypFrame$ppp, function(x) {
        table(marks(x, drop = FALSE)$gene)
        })
    hypFrame
}
#' Nest random effects within point patterns, by pasting the design factor
#' onto it
#'
#' @param df The dataframe
#' @param randomVars The random variables
#'
#' @return The dataframe with adapted randomVars
nestRandom = function(df, randomVars){
    stopifnot(any("design" == names(df)))
    for(i in randomVars){
        df[, i] = apply(df[, c("design", i)], 1, paste, collapse = "_")
    }
    df
}
#' A faster way to find the first element in a vector larger than a fixed value,
#'  by stopping early. This function is only faster than which.max is the match
#'  indeed comes early, as for nearest neighbour in a set of all distances.
#'
#' @param x An ordered vector
#' @param a A value to be queried
#'
#' @return The index
which.max.early = function(x, a){
    for(i in seq_along(x)){
        if(a<x[i])
            return(i)
    }
}


