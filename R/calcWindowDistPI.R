#' Estimate the PI for the distance to a fixed region of interest
#'
#' @param pSub The subset point pattern containing only a single gene
#' @param ecdfAll the cumulative distribution function under the null
#' @inheritParams estPisSingle
#' @param pi The type of PI to calculate
#' @param feat The feature
#' @param ecdfs The ecdfs of the distances of every single point in the cell
#' @param cellAndGene The marks dataframe with cell and gene, used to match the
#' observations with the ecdfs
#'
#' @return A list of vectors of estimated probabilistic indeces per event
#' @importFrom spatstat.geom nncross edges split.ppp
#' @importFrom stats dist
#' @details Analysis of the distance to the border was introduced by \insertCite{Joyner2013}{smoppix} in the form of the B-function.
#' The independent evaluations of the B-functions per cell are here returned as realizations of the probabilistic index.
#' @seealso \link{addCell}, \link{estPis}
#' @references
#' \insertAllCited{}
calcWindowDistPI <- function(pSub, owins, centroids, ecdfAll, pi, null, ecdfs,
                             cellAndGene, feat) {
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
