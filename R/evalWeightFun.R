#' Evaluate the weight function and return (unnormalized) weights
#'
#' @param wf The weight function
#' @param newdata Optional, a data frame with new data
#' @return A vector of predictions
#' @importFrom scam predict.scam
evalWeightFunction <- function(wf, newdata) {
    1 / exp(predict.scam(wf, newdata = newdata))
}
