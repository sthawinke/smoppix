#' Check in which cell each event lies and add a cell marker, plus include the list of regions
#'
#' @param hypFrame The hyperframe
#' @param owins the list of list of owins per hyperframe
#' @param checkOverlap a boolean, should windows be checked for overlap
#' @param warnOut a boolean, should warning be issued when points are not contained in window
#' @return the modified hyperframe
#' @importFrom matrixStats rowAnys
#' @importFrom spatstat.geom is.owin
#' @export
addCell = function(hypFrame, owins, checkOverlap = TRUE, warnOut = TRUE){
    stopifnot(nrow(hypFrame) == length(owins), all(unlist(lapply(owins, function(x) sapply(x, is.owin)))),
              all(rownames(hypFrame) %in% names(owins)))
    for(nn in rownames(hypFrame)){
        ppp = hypFrame[[nn, "ppp"]]
        if(any("cell" == names(marks(ppp, drop = FALSE)))){
            stop("Cell markers already present in point pattern ", nn)
        }
        if(checkOverlap){
            foo = findOverlap(owins[[nn]])
        }
        idWindow = vapply(owins[[nn]], FUN.VALUE = logical(NP <- npoints(ppp)), function(rr){
                inside.owin(ppp, w = rr)
        })
        idOut = which(!rowAnys(idWindow))
        if(warnOut){
            if(length(idOut) == npoints(ppp)){
                stop("All points lie outside all windows for point pattern ", nn, " check your input!")
            }
            warning(length(idOut), " points lie outside all windows for point pattern ", nn,
                    " and were not assigned to a cell.\n")
        }
        cellOut = character(NP)
        cellOut[idOut] = NA
        cellOut[-idOut] = rownames(which(t(idWindow[-idOut, ]), arr.ind = TRUE))
        hypFrame[[nn, "ppp"]] = setmarks(hypFrame[[nn, "ppp"]], cbind(marks(hypFrame[[nn, "ppp"]], drop = FALSE), cell = cellOut))
    }
    hypFrame$owins = owins
    return(hypFrame)
}

