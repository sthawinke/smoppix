#' Return all design variables, both at the level of the point pattern and the
#' level of the event
#'
#' @param x The results list
#' @param exclude variables to exclude
#'
#' @return A vector of likely design variables
getDesignVars = function(x, exclude = c("x", "y", "z", "tabObs", "centroids",
                                        "owins", "ppp")){
    pppVars = getPPPvars
    eventVars = getEventVars(x, exclude)
    c(pppVars, eventVars)
}
getPPPvars = function(x, exclude){
    setdiff(names(x$hypFrame), exclude)
}
getEventVars = function(x, exclude){
    setdiff(unique(lapply(x$hypFrame$ppp, function(ppp){
        names(marks(ppp, drop = FALSE))
    })), exclude)
}
