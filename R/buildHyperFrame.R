#' Built the hyperframe containing all separate point patterns that will
#'  be used throughout the analysis. Both matrices, dataframe and SpatialExperiment inputs are accepted.
#' @rdname buildHyperFrame
#' @export
setGeneric("buildHyperFrame", function(x, ...) standardGeneric("buildHyperFrame"))
#' @rdname buildHyperFrame
#' @export
setMethod("buildHyperFrame", "data.frame", function(x, ...) {
})
#' @rdname buildHyperFrame
#' @export
setMethod("buildHyperFrame", "SpatialExperiment", function(x, ...) {
})
