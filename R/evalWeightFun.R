#' Evaluate the weight function and return (unnormalized) weights
#'
#' @param wf The weight function
#' @param newdata Optional, a data frame with new data
#' @return A vector of predictions
#' @export
#' @importFrom scam predict.scam
#' @examples
#' data(Yang)
#' hypYang = buildHyperFrame(Yang, coordVars = c("x", "y"),
#' designVar = c("day", "root", "section"))
#' #Fit a subset of features to limit computation time
#' yangPims = estPims(hypYang, pis = c("nn", "nnPair"),
#' features = attr(hypYang, "features")[1:20])
#' #First build the weight function
#' wf <- buildWeightFunction(yangPims, pi = "nn", hypFrame = hypYang,
#' designVars = c("day", "root"))
#' evalWeightFunction(wf, newdata = data.frame("NP" = 1))
evalWeightFunction = function(wf, newdata){
    1/exp(predict.scam(wf, newdata = newdata))
}
