#' Extract results for a certain parameter
#'
#' @param obj The result object
#' @param parameter The desired parameter
#'
#' @return The matrix with result, with p-values in ascending order
#' @export
#'
#' @examples
getResults = function(obj, parameter){
    if(parameter == "Intercept")
        obj$results[[parameter]]
    else
        obj$results$fixedEffects[[parameter]]
}
