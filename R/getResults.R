#' Extract results for a certain parameter
#'
#' @param obj The result object
#' @param parameter The desired parameter
#'
#' @return The matrix with result, with p-values in ascending order
#' @export
#'
#' @examples
#' data(Yang)
#' hypYang <- buildHyperFrame(Yang,
#'     coordVars = c("x", "y"),
#'     imageVars = c("day", "root", "section")
#' )
#' # Fit a subset of features to limit computation time
#' yangPims <- estPims(hypYang, pis = "nn", features = attr(hypYang, "features")[seq_len(20)])
#' # First build the weight function
#' yangObj <- addWeightFunction(yangPims, designVars = c("day", "root"))
#' fittedModels <- fitLMMs(yangObj, fixedVars = "day", randomVars = "root")
#' res <- getResults(fittedModels, "Intercept")
#' head(res)
getResults <- function(obj, parameter) {
    if (parameter == "Intercept") {
        obj$results[[parameter]]
    } else {
        obj$results$fixedEffects[[parameter]]
    }
}
