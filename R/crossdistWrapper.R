#' A wrapper for C-functions calculating cross-distance matrix fast
#' @param x,y the matrices or point patterns between which to calculate the
#' cross distances
#' @return a matrix of cross distances
crossdistWrapper <- function(x, y) {
    crossdistFastLocal(getCoordsMat(x), getCoordsMat(y))
}
