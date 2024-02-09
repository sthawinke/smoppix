#' Build a weight matrix for Moran's I calculations
#' @param coordList A list of coordinate matrices
#' @return A sparse matrix of weights for Moran's I statistic, of equal value
#' for the numNNs nearest neighbours and zero otherwise
#' @importFrom Matrix sparseMatrix
buildMoransIWeightMat <- function(coordList, numNNs) {
    coordMat <- t(vapply(coordList, FUN.VALUE = double(2), getCoordsMat))
    nns = nnwhich(X = coordMat[, "x"], Y = coordMat[, "y"], k = seq.int(numNNs))
    # Assign non-zero slots to all numNNs nearest neighbours
    sparseMatrix(i = rep(seq_len(nrow(coordMat)), numNNs), j = nns, x = 1/numNNs, symmetric = FALSE)
    # Sparse, but not symmetric
}

