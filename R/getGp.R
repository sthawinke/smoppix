#' Helper function to get gene pair from a vector or list
#'
#'When provided with argument "geneA--geneB", looks for this gene pair as well as for
#'"geneB--geneA" in the provided object.
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
    if (isVec <- (is.vector(x) || is.list(x) || is.table(x))) {
        Names <- names(x)
    } else if (isMat <- is.matrix(x)) {
        Names <- rownames(x)
    }
    if (any(gp == Names)) {
        if (isVec) {
            x[[gp]]
        } else {
            x[gp, , drop = drop]
        }
    } else {
        gp <- paste(rev(sund(gp)), collapse = Collapse)
        if (any(gp == Names)) {
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
