#' A wrapper for C-functions calculating rowSums(x <y) with y a vector and x a
#' (possible big) matrix
#'@param x,y the matrices or point patterns between which to calculate the
#'cross distances
#'@return a matrix of cross distances
#'@importFrom bigmemory is.big.matrix
rowSumsWrapper = function(x, y){
    if(is.big.matrix(x)){
        rowSumsLarger(y, x@address)
    } else {
        rowSums(x < y)
    }

}
