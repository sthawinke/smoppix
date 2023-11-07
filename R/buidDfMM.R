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
    stopifnot(length(gene) %in% c(1, 2))
    if(length(gene) ==2){
        if(!any(grepl(pattern = "Pair", pi)))
            stop("Provide a single gene for univariate PIs!")
        gene = paste(gene, collapse = "_")
    }
    pi = match.arg(pi)
    foo = checkAttr(pimRes, pim)
    piListName = if(pairId <- grepl("Pair", pi)) "biPIs" else "uniPIs"
    piListNameInner = if(fixedId <- any(pi == c("edge", "midpoint", "fixedpoint"))) "windowDists" else "pointDists"
    nameVec = if(fixedId) pi else if(pairId) c(pi, "minNP", "maxNP") else c(pi, "NP")
    winId = any(pi == c("edge", "midpoint"))
    piEsts = lapply(pimRes, function(x){
        vec = (if(pairId) getGp(x[[piListName]], gene)[nameVec] else x[[piListName]][[gene]][[piListNameInner]][nameVec])
        if(winId) {
            data.frame("pi" = unlist(vec), "cell" = rep(names(vec[[pi]]), times = vapply(vec[[pi]], FUN.VALUE = integer(1), length)))
        } else if (fixedId){
            cbind(vec[[pi]])
        } else
            c(vec)
    })
    Design = if(fixedId) rep(names(pimRes), times = vapply(piEsts, FUN.VALUE = integer(1), if(winId) NROW else length)) else names(pimRes)
    piMat = Reduce(piEsts, f = rbind)
    if(is.null(piMat))
        stop("Gene or gene pair not found!\n")
    if(!fixedId){
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

