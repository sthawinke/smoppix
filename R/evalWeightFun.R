#' Evaluate the weighting function and return (unnormalized) weights
#'
#' @param wf The weighting function
#' @param newdata Optional, a data frame with new data
#' @return A vector of predictions
#' @importFrom scam predict.scam
#' @export
#' @examples
#' data(Yang)
#' hypYang <- buildHyperFrame(Yang, coordVars = c('x', 'y'),
#'     imageVars = c('day', 'root', 'section')
#' )
#' yangPims <- estPims(hypYang[c(seq_len(5), seq(25, 29)),],
#'     pis = 'nn', features = getFeatures(hypYang)[1:10])
#' # First Build the weighting function
#' yangObj <- addWeightFunction(yangPims, designVars = c('day', 'root'))
#' evalWeightFunction(yangObj$Wfs$nn, newdata = data.frame('NP' = 2))
evalWeightFunction <- function(wf, newdata) {
    1/exp(predict.scam(wf, newdata = newdata))
}
#' Add weights to PI matrix
#'
#' @param piMat The exisitng matrix with PI and design variables
#' @inheritParams buildDfMM
#'
#' @return The matrix with the weights column added
addWeights = function(piMat, pi, obj){
    weight <- evalWeightFunction(obj$Wfs[[pi]],
                                 newdata = piMat[, if(grepl("Pair", pi)) c("minP", "maxP") else "NP", drop = FALSE])
    weight <- weight/sum(weight, na.rm = TRUE)
    piMat <- cbind(piMat, weight)
}
