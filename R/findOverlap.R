#' Try to find overlap between list of windows, and throw error when found
#'
#' @param owins the list of windows
#' @param returnIds A boolean, should the indices of the overlap be returned?
#' If FALSE an error is thrown at the first overlap
#'
#' @return Throws an error when overlap found, otherwise returns invisible
#' @importFrom utils combn
#' @importFrom spatstat.geom overlap.owin is.owin
#' @importFrom BiocParallel bplapply
#' @export
#' @examples
#' library(spatstat.geom)
#' owins = replicate(10, owin(xrange = runif(1)+c(0,0.2),
#' yrange = runif(1)+c(0,0.1)), simplify = FALSE)
#' \dontrun{
#' findOverlap(owins)
#' }
#' idOverlap = findOverlap(owins, returnIds = TRUE)
findOverlap <- function(owins, returnIds = FALSE) {
    stopifnot(all(vapply(owins, is.owin, FUN.VALUE = TRUE)))
    combs <- combn(length(owins), 2)
    Ids = bplapply(seq_len(ncol(combs)), function(i) {
        Overlap <- overlap.owin(owins[[combs[1, i]]], owins[[combs[2, i]]]) >0
        if(returnIds){
            return(Overlap)
        } else if(Overlap){
                stop("Overlap detected between windows ", combs[1, i], " and ",
                     combs[2, i], "!\nWindows must be non-overlapping.")
        }
    })
    return( if (returnIds) {combs[, unlist(Ids)]} else {invisible()})
}
