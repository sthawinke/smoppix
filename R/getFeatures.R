#' Extract all unique features from an object
#'
#' @param x A hyperframe or a results list containing a hyperframe
#'
#' @return A vector of features.
#' @export
#' @details If x is a hyperFrame, all of its features are returned. If it is a result from estPis,
#' only features with estimated PIs are returned.
#' @examples
#' data(Yang)
#' hypYang <- buildHyperFrame(Yang,
#'     coordVars = c("x", "y"),
#'     imageVars = c("day", "root", "section")
#' )
#' head(getFeatures(hypYang))
getFeatures <- function(x) {
    if (is.hyperframe(x)) {
        unique(unlist(lapply(x$tabObs, names)))
    } else if (is.hyperframe(x$hypFrame)) {
        x$features
    }
}
