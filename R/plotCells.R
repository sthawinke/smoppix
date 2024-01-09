#' Plot the n cells with highest abundance of a feature
#' @description After testing for within-cell patterns, it may be useful to
#' look at the cells with the most events for certain genes. These are plotted
#' here, but the spatial location of the cells is lost
#'
#' @param obj A hyperframe, or an object containing one
#' @param feature The feature to be plotted, a character string
#' @param nCells An integer, the number of cells to be plotted
#'
#' @return Plots to the plotting window, returns invisible
#' @export
#' @importFrom spatstat.geom plot.owin subset.ppp
#'
#' @examples
#' example(addCell, "spatrans")
#' plotCells(hypFrame2, "gene1")
plotCells = function(obj, features = getFeatures(obj)[seq_len(3)], nCells = 1e2, Cex =1){
    if(!is.hyperframe(obj)){
        obj = obj$hypFrame
    }
    Cols <- makeCols(features, hypFrame)
    stopifnot(is.hyperframe(obj), !is.null(obj$owins),
              length(nCells) == 1)
    tablesCell = lapply(obj$ppp, function(p){
        table(marks(p[marks(p, drop = FALSE)$gene %in% features, ], drop = FALSE)$cell)
    })
    nCells = min(nCells - 1, length(ul <- unlist(tablesCell)))
    nthcell = sort(ul, decreasing = TRUE)[nCells]
    tablesCell = lapply(tablesCell, function(x){
        x[x >= nthcell]
    })
    counter = 0L;limVec = c(0, 1 + (Ceil <- ceiling(sqrt(nCells))))
    plot(type = "n", x = limVec, y = limVec, xlab = "", ylab = "", xaxt = "n",
         yaxt = "n", frame.plot = FALSE)
    for(i in seq_along(tablesCell)){
        nam = names(tablesCell)[i]
        ppp = subset.ppp(obj$ppp[[nam]], gene %in% features)
        for(j in seq_along(tablesCell[[nam]])){
            namIn = names(tablesCell[[nam]])[j]
            shifted = toUnitSquare(obj$owins[[nam]][[namIn]], ppp = subset.ppp(ppp, cell == namIn),
                                       Shift = Shift <- c(counter %% Ceil, counter %/% Ceil))
            plot.owin(shifted$owin, add = TRUE)
            #text(x = Shift[1], y = Shift[2], labels = paste0("Point pattern ", nam, "\nCell", namIn))
            points(coords(shifted$ppp), col = Cols[marks(shifted$ppp, drop = FALSE)$gene], pch = ".", cex = Cex)
            counter = counter + 1
        }
    }
    addLegend(Cols, Shift+ c(1,0))
    invisible()
}
#' @importFrom spatstat.geom affine.owin affine.ppp
toUnitSquare = function(owin, ppp, Shift){
    shrink = 1/rep(max(diff(owin$xrange), diff(owin$yrange)), 2)
    vec = -c(owin$xrange[1], owin$yrange[1])*shrink+ Shift
    list("owin" = affine.owin(owin, diag(shrink), vec),
         "ppp" = affine.ppp(ppp, diag(shrink), vec))
}
