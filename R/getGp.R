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
    isVec <- (is.vector(x) || is.list(x) || is.table(x))
    if(inherits(attempt <- tryGetEl(x, gp, isVec), "try-error")){
      gp <- paste(rev(sund(gp)), collapse = Collapse)
      attempt <- tryGetEl(x, gp, isVec)
    }
    attempt <- if(inherits(attempt, "try-error")){
      NULL
    } else {attempt}
    return(attempt)
}
getEl = function(x, gp, isVec){
  if (isVec) {
    x[[gp]]
  } else {
    x[gp, , drop = drop]
  }
}
tryGetEl = function(...){
  try(getEl(...), silent = TRUE)
}