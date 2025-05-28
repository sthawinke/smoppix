#' getDesignVars() returns all design variables, both at the level of the point pattern and the
#' level of the event
#'
#' @param x The results list, output from estPis
#'
#' @return A vector of design variables
#' @details getDesignVars() returns all design variables, \link{getPPPvars} returns
#' design variables related to the different images and \link{getEventVars} returns
#' design variables related to the individual events
getDesignVars <- function(x) {
    c(getPPPvars(x), getEventVars(x))
}
#' getPPPvars() extract variables related to point patterns
#'
#' @param exclude variables to exclude
#' @rdname getDesignVars
getPPPvars <- function(x, exclude = c(
                           "tabObs", "centroids", "owins", "ppp", "pimRes",
                           "image", "nuclei"
                       )) {
    setdiff(names(getHypFrame(x)), exclude)
}
#' getEventVars() extracts variables related to events, e.g. the cell (the marks)
#' @rdname getDesignVars
getEventVars <- function(x, exclude = c("x", "y", "z")) {
    setdiff(unique(unlist(lapply(getHypFrame(x)$ppp, function(ppp) {
        names(marks(ppp, drop = FALSE))
    }))), exclude)
}
getDiscreteVars <- function(x){
    eventVars <- unique(unlist(lapply(getHypFrame(x)$ppp, function(ppp) {
        names(marks(ppp, drop = FALSE))[!vapply(marks(ppp, drop = FALSE), FUN.VALUE = TRUE, is.numeric)]
    })))
    pppVars <- names(getHypFrame(x))[!vapply(getHypFrame(x), FUN.VALUE = TRUE, is.numeric)]
    c(eventVars, pppVars)
}
