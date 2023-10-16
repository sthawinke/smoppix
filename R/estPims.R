#' A wrapper fucntion for the different pims
#'
#' @param p The point patterns
#' @param pis the Probabilitstic indices to be estimated
#' @param ... Additional arguments passed on to the appropriate functions
#'
#' @return Data frames with estimated quantities per gene and/or gene pair
#' @export
#' @importFrom matrixStats rowMins
#' @importFrom BiocParallel bplapply
#' @examples
estPims = function(p, pis = c("nn", "allDist", "nnPair", "allDistPair", "edge", "centroid"),
                   null = c("background", "CSR"), nSims = 1e2, ...){
    pis = match.arg(pis, several.ok = TRUE)
    null = match.arg(null)
    tabObs = table(marks(p)$gene)
    unFeatures = names(tabObs); names(unFeatures) = unFeatures
    dfFeat = data.frame("feature" = unFeatures)
    if(any(pis == "allDist")){
        simDists = vapply(integer(nSims), FUN.VALUE = double(npoints(pSub)), function(i){

        })
    }
    uniPIs = bplapply(unFeatures, function(feat){
        pSub = subset(p, gene == feat)
        NNdistPI = if(any(pis == "nn")){calcNNPI(pSub, null, nSims)} else NULL
        allDistPI = if(any(pis == "allDist")){} else NULL
        c("NNdistPI" = NNdistPI, "allDistPI" = allDistPI)
    })
    }

}
