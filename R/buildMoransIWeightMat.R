#' Build a weight matrix for Moran's I calculations
#' @param coordMat A matrix of centroid coordinates
#' @inheritParams fitLMMs
#' @return A sparse matrix of weights for Moran's I statistic, of equal value
#' for the numNNs nearest neighbours and zero otherwise
#' @importFrom Matrix sparseMatrix
#' @importFrom spatstat.geom nnwhich
buildMoransIWeightMat <- function(coordMat, numNNs) {
    numNNs <- min(numNNs, nrow(coordMat) - 3)
    if (numNNs < 0) {
        matrix(0, nrow(coordMat), nrow(coordMat))
    }
    # Maximum as many NN as there are neighbouring cells
    nns <- nnwhich(X = coordMat[, "x"], Y = coordMat[, "y"], k = seq.int(numNNs))
    # Assign non-zero slots to all numNNs nearest neighbours
    sparseMatrix(i = rep(seq_len(nrow(coordMat)), numNNs), j = nns, x = 1/numNNs,
        symmetric = FALSE, dimnames = list(rownames(coordMat), rownames(coordMat)))
    # Sparse, but not symmetric, as nearest neighbours are not always
    # reciprocal
}

