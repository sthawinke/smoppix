#' Extract all unique and estimated features from an object
#'
#' @param x A hyperframe or a results list containing a hyperframe
#'
#' @return A vector of features
#' @export
#' @examples
#' data(Yang)
#' hypYang <- buildHyperFrame(Yang,
#'     coordVars = c("x", "y"),
#'     imageVars = c("day", "root", "section")
#' )
#' head(getFeatures(hypYang))
getFeatures <- function(x) {
    unique(unlist(lapply(getHypFrame(x)$tabObs, names)))
}
getEstFeatures <- function(x){
    x$features 
}
