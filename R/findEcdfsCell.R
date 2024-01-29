#' Construct ecdfs for cellwise measures, such as distance to edge or centroid
#'
#' @inheritParams estPisSingle
#' @importFrom spatstat.random runifpoint
#' @importFrom spatstat.geom crossdist nncross
#' @return The list of ecdf functions
findEcdfsCell <- function(p, owins, nPointsAllWin, centroids, null, pis, loopFun) {
    ecdfsCell <- loadBalanceBplapply(names(owins), function(nam) {
        pSub <- switch(null,
            "CSR" = runifpoint(nPointsAllWin, win = owins[[nam]]),
            "background" = subSampleP(pCell <- p[marks(p,
                drop = FALSE
            )$cell == nam, ], nPointsAllWin)
        )
        edge <- if (any(pis == "edge")) {
            ecdf(nncross(pSub, edges(owins[[nam]]), what = "dist"))
        }
        centroid <- if (any(pis == "centroid")) {
            ecdf(crossdist(pSub, centroids[[nam]]))
        }
        list(edge = edge, centroid = centroid)
    }, loopFun = loopFun)
    names(ecdfsCell) <- names(owins)
    ecdfsCell
}
