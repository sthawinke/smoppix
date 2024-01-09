#' Plot the n cells with highest abundance of a feature
#' @description After testing for within-cell patterns, it may be useful to
#' look at the cells with the most events for a certain gene. These are plotted
#' here, but the spatial location of the cells is lost
#'
#' @param obj A hyperframe, or an object containing one
#' @param feature The feature to be plotted, a character string
#' @param nCells An integer, the number of cells to be plotted
#'
#' @return Plots to the plotting window, returns invisible
#' @export
#' @importFrom spatstat.geom plot.owin
#'
#' @examples
#' example(addCell, "spatrans")
#' plotCells(hypFrame2, "gene1")
plotCells = function(obj, feature, nCells = 1e2){
    if(!is.hyperframe(obj)){
        obj = obj$hypFrame
    }
    stopifnot(length(feature)==1, is.hyperframe(obj), !is.null(obj$owins),
              length(nCells) == 1)
    tablesCell = lapply(obj$ppp, function(p){
        table(marks(p[marks(p, drop = FALSE)$gene == feature, ], drop = FALSE)$cell)
    })
    nthcell = sort(unlist(tablesCell), decreasing = TRUE)[nCells]
    tablesCell = lapply(tablesCell, function(x){
        x[x >= nthcell]
    })
    counter = 0L;limVec = c(0, ceiling(sqrt(nCells)))
    plot(type = "n", x = limVec, y = limVec, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
    for(i in seq_along(tablesCell)){
        nam = names(tablesCell)[i]
        for(j in seq_along(obj$owins[[nam]])){
            namIn = names(obj$owins[[nam]])[j]
            shiftedOwin = owinToUnitSquare(obj$owins[[nam]][[namIn]], Shift = c(counter %% nCells, counter %/% nCells))
            plot.owin(shiftedOwin, main = paste0("Point pattern ", nam, "\nCell", namIn), add = TRUE)
            counter = counter + 1
        }
    }
    invisible()
}
#' @importFrom spatstat.geom affine.owin
owinToUnitSquare = function(owin, Shift){
    vec = -c(owin$xrange[1], owin$yrange[1]) + Shift
    shrink = 1/rep(max(diff(owin$xrange), diff(owin$yrange)), 2)
    affine.owin(owin, diag(shrink), vec*shrink)
}
