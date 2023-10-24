#' Check in which cell each event lies and add a cell marker, plus include the list of regions
#'
#' @param hypFrame The hyperframe
#' @param rois the list of list of rois per hyperframe
#' @param checkOverlap a boolean, should windows be checked for overlap
#' @param warnout a boolean, should warning be issued when points are not contained in window
#' @return the modified hyperframe
#' @importFrom matrixStats rowAnys
#' @importFrom spatstat.geom is.owin
#' @export
addCell = function(hypFrame, rois, checkOverlap = TRUE, warnOut = TRUE){
    stopifnot(nrow(hypFrame) == length(rois), all(unlist(lapply(rois, function(x) sapply(x, is.owin)))),
              all(rownames(hypFrame) %in% names(rois)))
    for(nn in rownames(hypFrame)){
        ppp = hypFrame[[nn, "ppp"]]
        if(any("cell" == names(marks(ppp)))){
            stop("Cell markers already present in point pattern ", nn)
        }
        if(checkOverlap){
            foo = findOverlap(rois[[nn]])
        }
        idWindow = vapply(rois[[nn]], FUN.VALUE = logical(NP <- npoints(ppp)), function(rr){
                inside.owin(ppp, w = rr)
        })
        idOut = which(!rowAnys(idWindow))
        if(warnOut){
            if(length(idOut) == npoints(ppp)){
                stop("All points lie outside all windows for point pattern ", nn, " check your input!")
            }
            warning(length(idOut), " points lie outside all windows for point pattern ", nn,
                    "and were not assigned to a cell.\n")
        }
        cellOut = character(NP)
        cellOut[idOut] = NA
        cellOut[-idOut] = names(rois[[nn]])[which(idWindow, arr.ind = TRUE)[, "col"]]
        hypFrame[[nn, "ppp"]] = setmarks(hypFrame[[nn, "ppp"]], cbind(marks(hypFrame[[nn, "ppp"]]), cell = cellOut))
    }
    hypFrame$rois = rois
    return(hypFrame)
}
findOverlap = function(rois){
    combs = combn(length(rois), 2)
    lapply(seq_len(ncol(combs)), function(i){
        if(overlap.owin(rois[[combs[1, i]]], rois[[combs[2, i]]])>0){
            stop("Overlap detected between windows ", combs[1, i], " and ", combs[2, i], " !")
        }
    })
}
