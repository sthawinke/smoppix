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
    if(any(pis == "nn")){
        #obsNNdists = nndist(pp, by = "gene")
        NNdists = bplapply(unFeatures, function(feat){
            pSub = subset(p, gene == feat)
            obsDist = nndist(pSub)
            simDists = vapply(integer(nSims), FUN.VALUE = double(npoints(pSub)), function(i){
                distMat = crossdist(pSub, p[sample(npoints(p), npoints(pSub)-1),])
                #Keep observed points fixed
                matrixStats::rowMins(distMat)
            })
            piEsts = vapply(seq(npoints(pSub)), FUN.VALUE = double(1), function(i){
                ecdf(simDists[i,])(obsDist[i])
            })
            mean(piEsts)#Average over all points, single ecdf per point
            #Resampling amounts to permutation
        })
    } else if(any(pis == "allDist")){

    }

}
