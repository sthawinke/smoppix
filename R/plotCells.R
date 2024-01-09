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
#'
#' @examples
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
        x[x>=nthcell]
    })
    on.exit(par(par(no.readonly = TRUE)))
    nrowAndCol = ceiling(sqrt(nCells))
    par(mfrow = rep(nrowAndCol, 2))
    foo = lapply(names(tablesCell), function(nam){
        lapply(names(obj$owins[[nam]]), function(namIn){
            plot(obj$owins[[nam]][[namIn]], main = paste0("Point pattern ", nam, "\nCell", namIn))
        })
    })
    invisible()
}
