#' Build a data frame for a single gene and measure to prepare for mixed model building
#'
#' @param pimRes The result of a call to estPims()
#' @param gene A character string indicating the desired gene or gene pair
#' @param pi character string indicating the desired pim as outcome in the linear mixed model
#'
#' @return A dataframe
#' @export
#'
#' @examples
buildDfMM = function(pimRes, gene, pi = c("nn", "allDist", "nnPair", "allDistPair", "edge", "midpoint", "fixedpoint")){
    stopifnot(length(gene)==1)
    pi = match.arg(pi)
    piListName = if(pairId <- grepl("Pair", pi)) "biPIs" else "uniPIs"
    piListNameInner = if(winId <- any(pi == c("edge", "midpoint"))) "windowDists" else "pointDists"
    nameVec = if(winId) pi else if(pairId) c(pi, "minNP", "maxNP") else c(pi, "NP")
    Design = names(pimRes)
    piEsts = lapply(pimRes, function(x){
        vec = (if(pairId) getGp(x[[piListName]], gene) else x[[piListName]][[gene]])[[piListNameInner]][nameVec]
        if(winId) {
            cbind("pi" = vec, "cell" = names(vec))
        } else
            cbind(vec)
    })
    piMat = Reduce(piEsts, f = cbind)
    if(!winId){
        piMat = cbind(piMat, "weight" = if(pairId){
                        wVec = switch(pi,
                                      "allDistPair" = piMat[, "minNP"],
                                      "nnPair" = rowSums(piMat[, c("minNP", "maxNP")], na.rm = TRUE))
                        wVec/sum(wVec, na.rm = TRUE)
                      } else {
                          piMat[, "NP"]/sum(piMat[, "NP"], na.rm = TRUE)
                      })
    }
    piDf = data.frame("design" = Design, piMat)
    return(piDf)
}
