#' Evaluate a variance weighting function
#' @description
#' Evaluate the variance weighting function to  return unnormalized weights
#'
#' @param wf The weighting function
#' @param newdata Optional, a data frame with new data
#' @return A vector of predictions
#' @importFrom scam predict.scam
#' @export
#' @examples
#' data(Yang)
#' hypYang <- buildHyperFrame(Yang,
#'     coordVars = c("x", "y"),
#'     imageVars = c("day", "root", "section")
#' )
#' yangPims <- estPims(hypYang[c(seq_len(5), seq(25, 29)), ],
#'     nPointsAll = 1e3,
#'     pis = "nn", features = getFeatures(hypYang)[1:10]
#' )
#' # First Build the weighting function
#' yangObj <- addWeightFunction(yangPims, designVars = c("day", "root"))
#' evalWeightFunction(yangObj$Wfs$nn, newdata = data.frame("NP" = 2))
evalWeightFunction <- function(wf, newdata) {
    1 / exp(predict.scam(wf, newdata = newdata))
}
