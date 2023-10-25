#' Estimate the PI for the distance to a fixed region of interest
#'
#' @param pSub The subset point pattern containing only a single gene
#' @param edge The psp object to which the distances are calculated
#' @param ecdfAll the cumulative distribution function under the null
#'
#' @importFrom matrixStats rowMins
#' @return The estimated probabilistic index
#' @importFrom spatstat.geom nncross centroid.owin edges
calcWindowDistPI = function(pSub, owins, ecdfAll, midPoint = FALSE){
    splitPPP = split.ppp(pSub, f = "gene")
    obsDistEdge = vapply(FUN.VALUE = double(1), names(splitPPP), function(x){
        Dist = if(midPoint){
            crossdist(splitPPP[[x]], centroid.owin(owins[[x]], as.ppp = TRUE))
        } else {
            nncross(splitPPP[[x]], edges(owins[[x]]), what = "dist")
        }
        mean(ecdfAll[[x]](Dist))
    })

}

