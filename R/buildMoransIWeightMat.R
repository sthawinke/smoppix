#' Build a weight matrix for Moran's I calculations
#' @param coordList A list of coordinate matrices
#' @inheritParams fitLMMs
#' @return A sparse matrix of weights for Moran's I statistic, of equal value
#' for the numNNs nearest neighbours and zero otherwise
#' @importFrom Matrix sparseMatrix
#' @importFrom spatstat.geom nnwhich
buildMoransIWeightMat <- function(coordList, numNNs) {
    numNNs <- min(numNNs, length(coordList)-3)
    if(numNNs < 0){
        matrix(0, nrow(coordList), nrow(coordList))
    }
    #Maximum as many NN as there are neighbouring cells
    coordMat <- t(vapply(coordList, FUN.VALUE = double(2), getCoordsMat))
    nns = nnwhich(X = coordMat[, 1], Y = coordMat[, 2], k = seq.int(numNNs))
    # Assign non-zero slots to all numNNs nearest neighbours
    sparseMatrix(i = rep(seq_len(nrow(coordMat)), numNNs), j = nns, x = 1/numNNs,
                 symmetric = FALSE, dimnames = list(names(coordList), names(coordList)))
    # Sparse, but not symmetric, as nearest neighbours are not always reciprocal
}

