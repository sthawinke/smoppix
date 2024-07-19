#' Find overlap between list of windows
#' @description The function seeks overlap between the list of windows supplied,
#' and throws an error when found or returns the id's when found.
#'
#' @param owins the list of windows
#' @param returnIds A boolean, should the indices of the overlap be returned?
#' If FALSE an error is thrown at the first overlap
#' @param centroids The centroids of the windows
#' @param numCentroids An integer, the number of cells with closest centroids to consider looking for overlap
#'
#' @return Throws an error when overlap found, otherwise returns invisible.
#' When returnIds=TRUE, the indices of overlapping windows are returned.
#' @importFrom utils combn
#' @importFrom spatstat.geom overlap.owin is.owin is.ppp
#' @export
#' @examples
#' library(spatstat.geom)
#' owins <- replicate(10, owin(
#'     xrange = runif(1) + c(0, 0.2),
#'     yrange = runif(1) + c(0, 0.1)
#' ), simplify = FALSE)
#' idOverlap <- findOverlap(owins, returnIds = TRUE)
findOverlap <- function(owins, centroids = NULL, returnIds = FALSE, numCentroids = 30) {
    stopifnot(all(vapply(owins, is.owin, FUN.VALUE = TRUE)), is.null(centroids) ||
        all(vapply(centroids, is.ppp, FUN.VALUE = TRUE)))
    if (is.null(centroids)) {
        centroids <- lapply(owins, centroid.owin, as.ppp = TRUE)
    }
    combs <- combn(length(owins), 2)
    # First find distances between centroids
    distCentroids <- vapply(seq_len(ncol(combs)), FUN.VALUE = double(1), function(i) {
        crossdistWrapper(centroids[[combs[1, i]]], centroids[[combs[2, i]]])
    })
    # Find x-th closest centroid for every window, and forget about the other
    # windows
    for (winNum in seq_along(owins)) {
        idWin <- which(colSums(combs == winNum) == 1)
        idNotClosest <- -idWin[which(distCentroids[idWin] > sort(distCentroids[idWin])[numCentroids])]
        if (length(idNotClosest)) {
            combs <- combs[, idNotClosest, drop = FALSE]
            distCentroids <- distCentroids[idNotClosest]
        }
    }
    Ids <- loadBalanceBplapply(seq_len(ncol(combs)), function(ss) {
        vapply(ss, FUN.VALUE = TRUE, function(i) {
            Overlap <- overlap.owin(owins[[combs[1, i]]], owins[[combs[2, i]]]) >
                0
            if (returnIds) {
                return(Overlap)
            } else if (Overlap) {
                stop("Overlap detected between windows ", combs[1, i], " and ", combs[
                    2,
                    i
                ], "!\nWindows must be non-overlapping.")
            }
        })
    })
    return(if (returnIds) {
        combs[, unlist(Ids)]
    } else {
        invisible()
    })
}
