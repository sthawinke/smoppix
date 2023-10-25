#' Estimate the PI for the distance to a fixed region of interest
#'
#' @param pSub The subset point pattern containing only a single gene
#' @param edge The psp object to which the distances are calculated
#' @param ecdfAll the cumulative distribution function under the null
#'
#' @importFrom matrixStats rowMins
#' @return The estimated probabilistic index
#' @importFrom spatstat.geom nncross
calcWindowDistPI = function(pSub, owins, ecdfAll, midPoint = FALSE){
    splitPPP = split.ppp(pSub, f = "gene")
    obsDistEdge = lapply(names(splitPPP), function(x){
        nncross(splitPPP[[x]], (if(midPoint) centroid.owin else edge)(owins[[x]]), what = "dist")
    })
    mean(ecdfAll(unlist(obsDistEdge)))
}

