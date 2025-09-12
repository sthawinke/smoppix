#' Helper function to get gene pair from a vector or list
#'
#'When provided with argument "geneA--geneB", looks for this gene pair as well as for
#'"geneB--geneA" in the provided object.
#'
#' @param x The object in which to look
#' @param gp A character string describing the gene pair
#' @param Collapse The character separating the gene pair
#' @param drop A boolean, should matrix attributes be dropped in [] subsetting
#' @param notFoundReturn value to return if element is not found
#'
#' @return The element sought
#' @export
#' @examples
#' mat <- t(cbind(
#'     "gene1--gene2" = c(1, 2),
#'     "gene1--gene3" = c(2, 3)
#' ))
#' getGp(mat, "gene3--gene1")
getGp <- function(x, gp, drop = TRUE, Collapse = "--", notFoundReturn = NULL) {
    if (!length(x)) {
        return(notFoundReturn)
    }
    isVec <- (is.vector(x) || is.list(x) || (is.table(x) && is.na(ncol(x))))

    if(inherits(attempt <- tryGetEl(x, gp, isVec, drop), "try-error")){
      gp <- paste(rev(sund(gp)), collapse = Collapse)
      attempt <- tryGetEl(x, gp, isVec, drop)
    }
    attempt <- if(inherits(attempt, "try-error")){
      if(!isVec){
         rep(notFoundReturn, ncol(x))
      } else {
        notFoundReturn
      }
    } else {attempt}
    return(attempt)
}
getEl <- function(x, gp, isVec, drop){
  if (isVec) {
    x[[gp]]
  } else {
    x[, gp, drop = drop]
  }
}
tryGetEl <- function(...){
  try(getEl(...), silent = TRUE)
}