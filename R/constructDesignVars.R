#' Run checks on design variables, or construct them if missing
#'
#' @param designVars The initial design variables
#' @param lowestLevelVar Variable indicating the lowest level of nesting
#' @param allCell A boolean, are all PIs cell-related?
#' @param resList The results list
#'
#' @return A vector of design variables
constructDesignVars <- function(designVars, lowestLevelVar, allCell, resList) {
    allVars <- getDesignVars(resList)
    if (missing(designVars) && missing(lowestLevelVar) && !allCell) {
        stop("Provide either designVars or lowestLevelVar for measures ", "not on cell level!")
    }
    if (!missing(designVars) && !missing(lowestLevelVar)) {
        stop("Provide either designVars or lowestLevelVar, both not both!")
    }
    if (missing(designVars)) {
        if (allCell) {
            lowestLevelVar <- "cell"
        }
        if (length(lowestLevelVar) != 1) {
            stop("lowestLevelVar must have length 1")
        }
        if (!(lowestLevelVar %in% allVars)) {
            stop("lowestLevelVar must be part of the design,
                    see getFeatures(resList)")
        }
        designVars <- setdiff(getPPPvars(resList), lowestLevelVar)
        # If missing take everything but the ones that are obviously no design
        # variables
    }
    if (any(idMissing <- !(designVars %in% allVars))) {
        stop("Design variables\n", designVars[idMissing], "\nnot found in hypFrame object")
    }
    return(designVars)
}
