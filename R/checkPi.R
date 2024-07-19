#' Check if the required PI's are present in the object
#'
#' @param x The result of the PI calculation, or a weighting function
#' @param pi A character string indicating the desired PI
#'
#' @return Throws an error when the PIs are not found, otherwise returns invisible
checkPi <- function(x, pi) {
    if (!((pi %in% x$pis) || pi %in% names(x))) {
        stop("Required PI not present in object.", "Rerun estPis with correct 'pis' argument")
    } else {
        invisible()
    }
}
#' Check if features are present in hyperframe
#'
#' @inheritParams estPis
#'
#' @return Throws error when features not found
checkFeatures <- function(hypFrame, features) {
    if (any(id <- !(features %in% getFeatures(hypFrame)))) {
        stop("Features ", features[id], " not found in hyperframe")
    } else {
        invisible
    }
}
