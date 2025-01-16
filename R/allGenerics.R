#' Build a hyperframe containing all point patterns of an experiment.
#' @description
#' Build a spatstat hyperframe with point patterns and metadata.
#' Matrices, dataframe, lists and SpatialExperiment inputs are accepted.
#' @rdname buildHyperFrame
#' @importFrom spatstat.geom ppp hyperframe
#' @importFrom methods setGeneric setMethod
#' @export
#' @param x the input object, see methods('buildHyperFrame')
#' @param covariates A matrix or dataframe of covariates
#' @param featureName The name of the feature identifier for the molecules.
#' @param ... additional constructor arguments
#' @return An object of class 'hyperframe' from the 'spatstat.geom' package
#' @seealso \link[spatstat.geom]{hyperframe}
#' @examples
#' data(Yang)
#' hypYang <- buildHyperFrame(Yang,
#'     coordVars = c("x", "y"),
#'     imageVars = c("day", "root", "section")
#' )
setGeneric("buildHyperFrame", function(x, ...) standardGeneric("buildHyperFrame"))