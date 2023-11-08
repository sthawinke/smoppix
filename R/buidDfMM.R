#' Build a data frame for a single gene and measure to prepare for mixed model building
#'
#' @param pimRes The result of a call to estPims()
#' @param gene A character string indicating the desired gene or gene pair
#' @param pi character string indicating the desired pim as outcome in the linear mixed model
#' @param weightFunction The weighing function, obtained by the buildWeightFunction() function
#'
#' @return A dataframe
#' @export
#' @importFrom scam predict.scam
buildDfMM = function(pimRes, gene, hypFrame, pi = c("nn", "allDist", "nnPair", "allDistPair", "edge", "midpoint", "fixedpoint"), weightFunction){
    stopifnot(length(gene) %in% c(1, 2))
    if(length(gene) ==2){
        if(!any(grepl(pattern = "Pair", pi)))
            stop("Provide a single gene for univariate PIs!")
        gene = paste(gene, collapse = "_")
    }
    if(!all(id <- (sund(gene) %in% attr(hypFrame, "features")))){
        stop("Features", sund(gene)[id], "not found in hyperframe")
    }
    pi = match.arg(pi)
    foo = checkAttr(pimRes, pi)
    piListName = if(pairId <- grepl("Pair", pi)) "biPIs" else "uniPIs"
    piListNameInner = if(fixedId <- any(pi == c("edge", "midpoint", "fixedpoint"))) "windowDists" else "pointDists"
    winId = any(pi == c("edge", "midpoint"))
    piEsts = lapply(seq_along(pimRes), function(n){
        x = pimRes[[n]]
        vec = (if(pairId) getGp(x[[piListName]], gene)[pi] else x[[piListName]][[gene]][[piListNameInner]][pi])
        if(idNull <- is.null(vec))
            vec = rep(NA, if(pairId) 2 else 1)
        piOut = if(winId) {
            if(idNull) data.frame("pi" = NA, "cell" = "NA") else data.frame("pi" = unlist(vec), "cell" = rep(names(vec[[pi]]), times = vapply(vec[[pi]], FUN.VALUE = integer(1), length)))
        } else if (fixedId){
            cbind("pi" = if(idNull) vec else vec[[pi]])
        } else
            c("pi" = unname(vec))
        if(!fixedId){
            npVec = hypFrame[n, "table", drop = TRUE][if(pairId) sund(gene) else gene]
            if(anyNA(npVec))
                npVec = rep(NA, if(pairId) 2 else 1)
            else if(pairId)
                npVec = sort(npVec)
                names(npVec) = if(pairId) c("minP", "maxP") else "NP"
        }
        if(fixedId) {piOut} else c(piOut, npVec)
    })
    Design = if(fixedId) rep(names(pimRes), times = vapply(piEsts, FUN.VALUE = integer(1), if(winId) NROW else length)) else names(pimRes)
    piMat = Reduce(piEsts, f = rbind)
    if(is.null(piMat))
        stop("Gene or gene pair not found!\n")
    if(!fixedId){
        weight = if(missing(weightFunction)){
            if(pairId){
                switch(pi, "allDistPair" = piMat[, "minP"],
                       "nnPair" = rowSums(piMat[, c("minP", "maxP")], na.rm = TRUE))

            } else {
                piMat[, "NP"]/sum(piMat[, "NP"], na.rm = TRUE)
            }
        } else {
            foo = checkAttr(weightFunction, pi)
            weights = evalWeightFunction(weightFunction, newdata = data.frame(piMat[, -1]))
        }
        weight = weight/sum(weight, na.rm = TRUE)
        piMat = cbind(piMat, weight)
    }
    piDf = data.frame("design" = Design, piMat)
    if(!fixedId)
        rownames(piDf) = names(pimRes)
    return(piDf)
}

