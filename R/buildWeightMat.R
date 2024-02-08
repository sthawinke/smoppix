#' Build a weight matrix for Moran's I calculations
#' @param coordList A list of coordinate matrices
#' @return A matrix of weights for Moran's I statistic, inversely proportional to the distance between cell centroids
#' @importFrom stats dist
buildWeightMat <- function(coordList) {
    coordMat <- t(vapply(coordList, FUN.VALUE = double(2), getCoordsMat))
    as.matrix(1 / dist(coordMat))
}
#Need to work with sparse, symmetric Matrix class. Threshold on 8 closest centroids or so
