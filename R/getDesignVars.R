#' Return all design variables, both at the level of the point pattern and the
#' level of the event
#'
#' @param x The results list
#'
#' @return A vector of likely design variables
#' @details getDesignVars() returns all design variables, getPPPvars returns
#' design variables related to the different images and getEventVars returs
#' design variables related to the individual events
getDesignVars <- function(x) {
    pppVars <- getPPPvars(x)
    eventVars <- getEventVars(x)
    c(pppVars, eventVars)
}
#' Extract variables from point patterns
#'
#' @param exclude variables to exclude
#' @inheritParams getDesignVars
#' @return A vector of variables
getPPPvars <- function(x, exclude = c("tabObs", "centroids", "owins", "ppp", "pimRes", "image", "inSeveralCells")) {
    setdiff(names(x$hypFrame), exclude)
}
#' Extract variables from events (the marks)
#' @return A vector of variables
#' @inheritParams getDesignVars
#' @inheritParams getPPPvars
getEventVars <- function(x, exclude = c("x", "y", "z")) {
    setdiff(unique(unlist(lapply((if (is.hyperframe(x)) {
        x
    } else {
        x$hypFrame
    })$ppp, function(ppp) {
        names(marks(ppp, drop = FALSE))
    }))), exclude)
}
