#' Try to find overlap between list of windows, and throw error when found
#'
#' @param owins the list of windows
#' @param returnIds A boolean, should the indices of the overlap be returned?
#' If FALSE an error is thrown at the first overlap
#'
#' @return Throws an error when overlap found, otherwise returns invisible
#' @importFrom utils combn
#' @importFrom spatstat.geom overlap.owin
findOverlap <- function(owins, returnIds = FALSE) {
    combs <- combn(length(owins), 2)
    Ids = vapply(seq_len(ncol(combs)), FUN.VALUE = TRUE, function(i) {
        if(Overlap <- overlap.owin(owins[[combs[1, i]]], owins[[combs[2, i]]]) > 0){
            if (returnIds) {
                return(Overlap)
            } else {
                stop("Overlap detected between windows ", combs[1, i], " and ",
                     combs[2, i], "!\nWindows must be non-overlapping.")
            }
        }
    })
    return(combs[, Ids])
}
