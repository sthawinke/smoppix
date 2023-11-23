#' Estimate the probabilistic index for the collection of all distances
#'
#' @param pSub The subset point pattern with only one feature
#' @param ecdfAll The empirical cumulative distribution function of all distances under the null
#' @param ecdfs A list of empirical distribution functions, one for each event
#'
#' @return The probabilistic index for all distances
#' @inheritParams estPimsSingle
#' @importFrom stats dist
calcAllDistPI = function(pSub, p, ecdfAll, null, ecdfs){
    obsDist = as.matrix(dist(coords(pSub)))
    if(null == "CSR"){
        mean(ecdfAll(obsDist)) #Need to condition on point too?
    } else if(null == "background"){
        piEsts = vapply(seq_along(ecdfs), FUN.VALUE = double(1), function(i){
            mean(ecdfs[[i]](obsDist[i, -i]))
        })
        mean(piEsts)
    }
}
calcAllDistPIpair = function(id1, id2, ecdfAll, null, crossDist, ecdfs){
    if(null == "CSR"){
        mean(ecdfAll(crossDist))
    } else if(null == "background"){
        # #Keep observed points fixed
        piEsts1 = vapply(seq_along(id1), FUN.VALUE = double(length(id2)), function(i){
            ecdfs[[i]](crossDist[i, ])
        })
        piEsts2 = vapply(seq_along(id2), FUN.VALUE = double(length(id1)), function(i){
            ecdfs[[i + length(id1)]](crossDist[,i])
        })
        mean(c(piEsts1, piEsts2))
    }
}

