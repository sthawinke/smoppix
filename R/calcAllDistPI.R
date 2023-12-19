#' Estimate the probabilistic index for the collection of all distances
#'
#' @param pSub The subset point pattern with only one feature
#' @param ecdfAll The empirical cumulative distribution function of
#' all distances under the null
#'
#' @return The probabilistic index for all distances
#' @inheritParams estPimsSingle
#' @inheritParams calcUniPIs
#' @importFrom stats dist
calcAllDistPI <- function(pSub, p, ecdfAll, null, cd) {
    obsDist <- stats::dist(coords(pSub))
    if (null == "CSR") {
        mean(ecdfAll(obsDist)) # Need to condition on point too?
    } else if (null == "background") {
        obsDist <- dropDiagonal(as.matrix(obsDist))
        mean(cd < obsDist)
    }
}
calcAllDistPIpair <- function(id1, id2, ecdfAll, null, crossDistSub, cd) {
    if (null == "CSR") {
        mean(ecdfAll(crossDistSub))
    } else if (null == "background") {
        # #Keep observed points fixed
        piEsts1 <- rowMeans(cd[id1,] < crossDistSub[id1, ])
        piEsts2 <- rowMeans(cd[id2,] < crossDistSub[, id2])
        mean(c(piEsts1, piEsts2))
    }
}
