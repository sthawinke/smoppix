#' @description getResults() Extracts effect size estimates, standard errors and adjusted p-values for a certain parameter
#' from a linear model.
#'
#' @param obj The result object
#' @param pi The desired PI
#' @param parameter The desired parameter
#' @param moransI A boolean, should results for the Moran's I be returned?
#'
#' @return For getResults(), the matrix with results, with p-values in ascending order
#' \item{Estimate}{The estimated PI}
#' \item{se}{The corresponding standard error}
#' \item{pVal}{The p-value}
#' \item{pAdj}{The Benjamini-Hochberg adjusted p-value}
#' @export
#' @rdname fitLMMs
#' @order 3
getResults <- function(obj, pi, parameter, moransI = FALSE) {
    if (!(pi %in% names(obj))) {
        stop("PI ", pi, " not estimated! Run estPis() again with desired pi as argument")
    }
    resChar <- if (moransI) {
        "resultsMoran"
    } else {
        "results"
    }
    Obj <- if (parameter == "Intercept") {
        obj[[pi]][[resChar]]
    } else {
        obj[[pi]][[resChar]]$fixedEffects
    }
    if (!(parameter %in% names(Obj))) {
        stop("Parameter ", parameter, " not estimated in linear model!
             Rerun fitLMMs with the correct variables or formula supplied.")
    } else {
        Obj[[parameter]]
    }
}
