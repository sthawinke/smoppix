#' Evaluate a variance weighting function
#' @description
#' Evaluate the variance weighting function to  return unnormalized weights
#'
#' @param wf The weighting function
#' @param newdata A data frame with new data
#' @return A vector of weights, so the inverse of predicted variances, unnormalized
#' @importFrom scam predict.scam
#' @export
#' @examples
#' data(Yang)
#' hypYang <- buildHyperFrame(Yang, coordVars = c("x", "y"),
#'     imageVars = c("day", "root", "section")
#' )
#' yangPims <- estPis(hypYang, pis = "nn", features = getFeatures(hypYang)[12:19], nPointsAll = 1e3)
#' # First Build the weighting function
#' yangObj <- addWeightFunction(yangPims, designVars = c("day", "root"))
#' evalWeightFunction(yangObj$Wfs$nn, newdata = data.frame("NP" = 2))
#' @seealso \link[scam]{predict.scam}, \link{addWeightFunction}
evalWeightFunction <- function(wf, newdata) {
    1 / exp(predict.scam(wf, newdata = newdata))
}
