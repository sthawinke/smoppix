#' Plot the hyperframe for exploratory purposes
#'
#' @param hypFrame The hyperframa
#' @param features A small number of features to be fitted. Defaults to the first 5
#' @param ppps The rownames or indices of the point patterns to be plotted. Defaults to the first 9.
#'
#' @return
#' @importFrom spatstat.geom is.hyperframe
#' @export
#'
#' @examples
#' data(Yang)
#' hypYang = buildHyperFrame(Yang, coordVars = c("x", "y"),
#' designVar = c("day", "root", "section"))
#' plotExplore(hypYang)
plotExplore = function(hypFrame, features = attr(hypFrame, "features"), ppps){
    stopifnot(is.hyperframe(hypFrame), is.character(features))
    npp = nrow(hypFrame)
    if(missing(ppps)){
        ppps = rownames(hypFrame)[seq_len(min(8, npp))]
    }
    Cols = palette(); Cols = c(Cols, rep("grey", length(features) - length(Cols)))
    names(Cols) = features
    old.par = par(no.readonly = TRUE)
    on.exit(par(old.par))
    par(mfrow = if((LL <- length(ppps)) <= 3) c(2,2) else c(3,3))
    baa = lapply(ppps, function(i){
        id = marks(hypFrame$ppp[[i]], drop = FALSE)$gene %in% features
        if(any(id))
            plot(coords(hypFrame$ppp[[i]][id,]), main = ppps[i], pch = ".", asp = 1,
                 col = Cols[marks(hypFrame$ppp[[i]], drop = FALSE)$gene[id]])
    })
    plot(0,0, type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n")
    legend("center", legend = names(Cols), pch = 20, col = Cols)
}
