#' Extract all unique features from an object, or the ones for which PIs were estimated
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
