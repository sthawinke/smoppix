#' Check if the required pi's are present in the object
#'
#' @param x The result of the PI calculation, or a weighting function
#' @param pi A character string indicating the desired PI
#'
#' @return Throws an error when the PIs are not found, otherwise returns invisible
checkPi <- function(x, pi) {
    if (!(pi %in% x$pis)) {
        stop("Required PI not present in object.",
        "Rerun estPis with correct 'pis' argument")
    } else {
        invisible()
    }
}
