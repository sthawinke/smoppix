#' Extract all unique features from an object
#'
#' @param x A hyperframe or a results list containing a hyperframe
#'
#' @return A vector of features
#' @export
#' @examples
#' data(Yang)
#' hypYang <- buildHyperFrame(Yang,
#'     coordVars = c('x', 'y'),
#'     imageVars = c('day', 'root', 'section')
#' )
#' head(getFeatures(hypYang))
getFeatures <- function(x) {
    x <- if (is.hyperframe(x)) {
        x
    } else if (is.hyperframe(x$hypFrame)) {
        x$hypFrame
    }
    unique(unlist(lapply(x$tabObs, names)))
}
