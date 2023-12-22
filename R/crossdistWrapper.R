#' A wrapper for C-functions calculating cross-distance matrix fast and with
#' low memory usage. Big matrices are also supported
#'@param returnBigMatrix a boolean, should the result be returned as
#'a big matrix (low RAM use) or a regular r-matrix (faster)
#'@param x,y the matrices or point patterns between which to calculate the
#'cross distances
#'@return a matrix of cross distances
#'@importFrom bigmemory big.matrix
crossdistWrapper = function(x, y, returnBigMatrix = FALSE, tmpFile = NULL){
    if(returnBigMatrix){
        x = getCoordsMat(x);y = getCoordsMat(y)
        descFile = if(!is.null(tmpFile)) paste0(tmpFile, ".desc")
        bigMat = big.matrix(nrow = nrow(x), ncol = nrow(y), type = "double",
                           backingfile = tmpFile, descriptorfile = descFile)
        crossdistFast(x, y, bigMat@address)
        return(bigMat)
    } else {
        crossdistFastLocal(getCoordsMat(x), getCoordsMat(y))
    }

}
