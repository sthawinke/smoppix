#' Estimate the PI for the distance to a fixed region of interest
#'
#' @param pSub The subset point pattern containing only a single gene
#' @param ecdfAll the cumulative distribution function under the null
#' @param owins the list of windows describing the cells
#' @param towhat a character vector, either "edge" or "midpoint", indicating with respect to what the distance should be calculated
#'
#' @return The estimated probabilistic index
#' @importFrom spatstat.geom nncross centroid.owin edges
#' @details Analysis of the distance to the border was introduced by \insertCite{Joyner2013}{spatrans} in the form of the B-function.
#' The independent evaluations of the B-functions per cell are here returned as realizations of the probabilistic index.
#' @seealso \link{addCell}
#' @references
#' \insertAllCited{}
calcWindowDistPI = function(pSub, owins, ecdfAll, towhat){
    splitPPP = split.ppp(pSub, f = "cell")
    obsDistEdge = lapply(names(splitPPP), function(x){
        Dist = switch(towhat,
            "edge" = crossdist(splitPPP[[x]], centroid.owin(owins[[x]], as.ppp = TRUE)),
            "midpoint" = nncross(splitPPP[[x]], edges(owins[[x]]), what = "dist"))
        ecdfAll[[x]](Dist) #Do not average here, independent observations
    })
    names(obsDistEdge) = names(splitPPP)
    obsDistEdge
}

