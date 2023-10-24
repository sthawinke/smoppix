#' Estimate the PI for the distance to a fixed region of interest
#'
#' @param pSub The subset point pattern containing only a single gene
#' @param edge The psp object to which the distances are calculated
#' @param ecdfAll the cumulative distribution function under the null
#'
#' @importFrom matrixStats rowMins
#' @return The estimated probabilistic index
#' @importFrom spatstat.geom nncross
calcEdgeDistPI = function(pSub, edge, ecdfAll){
    obsDistEdge = nncross(pSub, edge, what = "dist")
    mean(ecdfAll(obsDistEdge))
    #Convert to psp, plus vectorize. FIX ME!
}

