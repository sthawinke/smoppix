#' Estimate the probabilistic index for the collection of all distances
#'
#' @param pSub The subset point pattern with only one feature
#' @param ecdfAll The empirical cumulative distribution function of
#' all distances under the null
#'
#' @return The probabilistic index for all distances
#' @inheritParams estPimsSingle
#' @importFrom stats dist
calcAllDistPI <- function(pSub, p, ecdfAll, null, cd) {
    obsDist <- dist(coords(pSub))
    if (null == "CSR") {
        mean(ecdfAll(obsDist)) # Need to condition on point too?
    } else if (null == "background") {
        obsDist <- as.matrix(obsDist)
        piEsts <- vapply(seq_len(npoints(pSub)), FUN.VALUE = double(1),
                         function(i) {
            mean(ecdfPreSort(cd[i,])(obsDist[i, -i]))
        })
        mean(piEsts)
    }
}
calcAllDistPIpair <- function(id1, id2, ecdfAll, null, crossDistSub, cd) {
    if (null == "CSR") {
        mean(ecdfAll(crossDistSub))
    } else if (null == "background") {
        # #Keep observed points fixed
        piEsts1 <- vapply(seq_along(id1), FUN.VALUE = double(length(id2)),
                          function(i) {
            ecdfPreSort(cd[i,])(crossDistSub[i, ])
        })
        piEsts2 <- vapply(seq_along(id2), FUN.VALUE = double(length(id1)),
                          function(i) {
            ecdfPreSort(cd[i+ length(id1),])(crossDistSub[, i])
        })
        mean(c(piEsts1, piEsts2))
    }
}
