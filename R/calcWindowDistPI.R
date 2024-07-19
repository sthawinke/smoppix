#' Estimate the PI for the distance to a fixed object of interest, such as a cell wall or centroid
#'
#' @param pSub The subset point pattern containing only a single gene
#' @param ecdfAll the cumulative distribution function under the null
#' @inheritParams estPisSingle
#' @param pi The type of PI to calculate
#'
#' @return A list of vectors of estimated probabilistic indeces per event
#' @importFrom spatstat.geom nncross edges split.ppp
#' @details Analysis of the distance to the border was introduced by \insertCite{Joyner2013}{smoppix} in the form of the B-function.
#' The independent evaluations of the B-functions under the null hypothesis represented by \textit{ecdfAll} 
#' per cell are here returned as realizations of the probabilistic index.
#' @seealso \link{addCell}, \link{estPis}
#' @references
#' \insertAllCited{}
calcWindowDistPI <- function(pSub, owins, centroids, ecdfAll, pi) {
    splitPPP <- split.ppp(pSub, f = "cell")
    splitPPP <- splitPPP[names(splitPPP) != "NA"]
    obsDistEdge <- lapply(names(splitPPP), function(x) {
        Dist <- switch(pi,
            centroid = crossdistFastLocal(getCoordsMat(splitPPP[[x]]), centroids[x, , drop = FALSE]),
            edge = nncross(splitPPP[[x]], edges(owins[[x]]), what = "dist")
        )
        ecdfAll[[x]][[pi]](Dist) # Do not average here, independent observations
    })
    names(obsDistEdge) <- names(splitPPP)
    obsDistEdge
}
