#' Construct ecdfs for cellwise measures, such as distance to edge or centroid
#'
#' @inheritParams estPimsSingle
#'
#' @return The list of ecdf functions
findEcdfsCell = function(p, owins, nPointsAllWin, centroids, null){
    ecdfsCell <- lapply(names(owins), function(nam) {
        pSub <- switch(null,
                       "CSR" = runifpoint(nPointsAllWin, win = owins[[nam]]),
                       "background" = subSampleP(pCell <- p[marks(p, drop = FALSE)$cell == nam, ], nPointsAllWin)
        )
        edge <- if (any(pis == "edge")) {
            ecdf(nncross(pSub, edges(owins[[nam]]), what = "dist"))
        }
        midpoint <- if (any(pis == "midpoint")) {
            ecdf(crossdist(pSub, centroids[[nam]]))
        }
        list("edge" = edge, "midpoint" = midpoint)
    })
    names(ecdfsCell) <- names(owins)
    ecdfsCell
}
