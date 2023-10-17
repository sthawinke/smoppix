#'Extract elements from a distance matrix by index
#'
#' @param distMat the distance matrix of class 'dist'
#' @param i,j Indices for which entries are extracted
#' @details See documentation of the stats::dist function
getElDist = function(distMat, i, j){
    id = outer(attr(distMat, "Size")*(i-1) - i*(i-1)/2-i, j, FUN = "+")
    #Following dist help
    distMat[c(id)]
}
