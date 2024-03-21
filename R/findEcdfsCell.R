#' Construct empirical cumulative distribution functions (ecdfs) for distances within the cell
#' @description The distance distribution under the null hypothesis of complete spatial randomness (CSR)
#' within the cell is the same for all genes. This function precalculates this distribution using Monte-Carlo simulation under CSR,
#' and summarizes it in a ecdf object
#'
#' @inheritParams estPisSingle
#' @importFrom spatstat.random runifpoint
#' @importFrom spatstat.geom nncross
#' @return The list of ecdf functions
#' @seealso \link{ecdf}
findEcdfsCell <- function(p, owins, nPointsAllWin, centroids, null, pis, loopFun) {
    ecdfsCell <- loadBalanceBplapply(names(owins), function(nam) {
        pSub <- switch(null,
            CSR = runifpoint(nPointsAllWin, win = owins[[nam]]),
            background = subSampleP(pCell <- p[marks(p, drop = FALSE)$cell == nam, ], nPointsAllWin)
        )
        edge <- if (any(pis == "edge")) {
            ecdf(nncross(pSub, edges(owins[[nam]]), what = "dist"))
        }
        centroid <- if (any(pis == "centroid")) {
            ecdf(crossdistFastLocal(getCoordsMat(pSub), centroids[nam, , drop = FALSE]))
        }
        list(edge = edge, centroid = centroid)
    }, loopFun = loopFun)
    names(ecdfsCell) <- names(owins)
    ecdfsCell
}
