#' Estimate the PI for the distance to a fixed region of interest
#'
#' @param point The fixed point to which the distance is to be calculated
#' @inheritParams calcWindowDistPI
#' @importFrom spatstat.geom crossdist
calcFixedPointDistPI = function(pSub, pointPPP, ecdfAll){
    obsDistPoint = crossdist(pSub, pointPPP)
    mean(ecdfAll(obsDistPoint))
}
