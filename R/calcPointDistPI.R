#' Estimate the PI for the distance to a fixed region of interest
#'
#' @param point The fixed point to which the distance is to be calculated
#' @inheritParams calcEdgeDistPI
#' @importFrom spatstat.geom crossdist
calcPointDistPI = function(pSub, point, ecdfAll){
    obsDistPoint = crossdist(pSub, point)
    mean(ecdfAll(obsDistPoint))
}
