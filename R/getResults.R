#' Extract test results from the linear model
#' @description Extract effect size estimates and p-values for a certain parameter
#' from a linear model.
#'
#' @param obj The result object
#' @param parameter The desired parameter
#' @param moransI A boolean, should results for the Moran's I be returned?
#'
#' @return The matrix with results, with p-values in ascending order
#' \item{Estimate}{The estimated PI}
#' \item{se}{The corresponding standard error}
#' \item{pVal}{The p-value}
#' \item{pAdj}{The Benjamini-Hochberg adjusted p-value}
#' @export
#' @examples
#' data(Yang)
#' hypYang <- buildHyperFrame(Yang,
#'     coordVars = c('x', 'y'),
#'     imageVars = c('day', 'root', 'section')
#' )
#' yangPims <- estPis(hypYang, pis = 'nn')
#' # First build the weighting function
#' yangObj <- addWeightFunction(yangPims, designVars = c('day', 'root'))
#' fittedModels <- fitLMMs(yangObj,
#'     features = getFeatures(yangObj)[seq_len(10)],
#'     fixedVars = 'day', randomVars = 'root',
#'     pi = 'nn'
#' )
#' res <- getResults(fittedModels, 'Intercept')
#' head(res)
getResults <- function(obj, parameter, moransI = FALSE) {
    resChar <- if (moransI) {
        "resultsMoran"
    } else {
        "results"
    }
    if (parameter == "Intercept") {
        obj[[resChar]][[parameter]]
    } else {
        obj[[resChar]]$fixedEffects[[parameter]]
    }
}
