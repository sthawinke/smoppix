#' Helper function to get gene pair or order reversed from a vector or list
#'
#' @param x The object in which to look
#' @param gp A character string describing the gene pair
#' @param Collapse The character separating the gene pair
#' @param drop A boolean, should matrix attributes be dropped in [] subsetting
#'
#' @return The element sought
#' @export
#' @examples
#' mat <- t(cbind(
#'     "gene1--gene2" = c(1, 2),
#'     "gene1--gene3" = c(2, 3)
#' ))
#' getGp(mat, "gene3--gene1")
getGp <- function(x, gp, drop = TRUE, Collapse = "--") {
    if (!length(x)) {
        return(NULL)
    }
    if (isVec <- (is.vector(x) || is.list(x))) {
        Names <- names(x)
    } else if (isMat <- is.matrix(x)) {
        Names <- rownames(x)
    }
    if (gp %in% Names) {
        if (isVec) {
            x[[gp]]
        } else {
            x[gp, , drop = drop]
        }
    } else {
        geneSplit <- sund(gp)
        gp <- paste(rev(geneSplit), collapse = Collapse)
        if (gp %in% Names) {
            if (isVec) {
                x[[gp]]
            } else {
                x[gp, , drop = drop]
            }
        } else {
            NULL
        }
    }
}
