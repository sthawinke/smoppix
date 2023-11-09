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
buildDfMM = function(pimRes, gene, pi = c("nn", "allDist", "nnPair", "allDistPair", "edge", "midpoint", "fixedpoint"),
                     weightFunction, ...){
    #Establish whether pi and gene match, and call separate functions
    stopifnot((lg <- length(gene)) %in% c(1, 2))
    #Check whether genes are present
    if(!all(id <- (sund(gene) %in% attr(hypFrame, "features")))){
        stop("Features", sund(gene)[id], "not found in hyperframe")
    }
    pi = match.arg(pi)
    foo = checkAttr(pimRes, pi)
    if(!missing(weightFunction))
        foo = checkAttr(weightFunction, pi)
    #To do: switch to quixotic gene pair separation tag
    if(pairId <- grepl(pattern = "Pair", pi)){
        if(lg == 2){
            gene = paste(gene, collapse = "&")
        }
        buildDfMMBi(gene = gene, pi = pi, pimRes = pimRes, weightFunction = weightFunction, ...)
    } else {
        if(lg !=1){
            stop("Provide a one gene for univariate PIs!")
        }
        buildDfMMUni(gene = gene, pi = pi, pimRes = pimRes, weightFunction = weightFunction,...)
    }
}
buildDfMMBi = function(pimRes, gene, hypFrame, pi, weightFunction){
    piEsts = t(vapply(seq_along(pimRes), FUN.VALUE = double(3), function(n){
        if(idNull <- is.null(vec <- getGp(pimRes[[n]][["biPIs"]], gene)[pi]))
            vec = NA
        npVec = sort(hypFrame[n, "table", drop = TRUE][sund(gene)])
        names(npVec) = c("minP", "maxP")
        c("pi" = unname(vec), npVec)
    }))
    if(all(is.na(piEsts[,"pi"])))
        stop("Gene pair not found!\n")
    weight = if(missing(weightFunction)){
            piEsts[,"minP"]
        } else {
            evalWeightFunction(weightFunction, newdata = data.frame(piEsts[, c("minP", "maxP")]))
        }
    weight = weight/sum(weight, na.rm = TRUE)
    piMat = cbind(piMat, weight)
    piDf = data.frame("design" = names(pimRes), piMat)
    rownames(piDf) = names(pimRes)
    return(piDf)
    }
buildDfMMUni = function(pimRes, gene, hypFrame, pi, weightFunction){
    piListNameInner = if(fixedId <- any(pi == c("edge", "midpoint", "fixedpoint"))) "windowDists" else "pointDists"
    #fixedId: Is the distance to a fixed point?
    winId = any(pi == c("edge", "midpoint")) #winId: is a window involved?
    piEsts = lapply(seq_along(pimRes), function(n){
        x = pimRes[[n]]
        if(idNull <- is.null(vec <- x[["uniPIs"]][[gene]][[piListNameInner]][pi]))
            vec = NA
        piOut = if(winId) {
            if(idNull) data.frame("pi" = NA, "cell" = "NA") else data.frame("pi" = unlist(vec), "cell" = rep(names(vec[[pi]]), times = vapply(vec[[pi]], FUN.VALUE = integer(1), length)))
        } else if (fixedId){
            cbind("pi" = if(idNull) vec else vec[[pi]])
        } else
            c("pi" = unname(vec))
        if(!fixedId){
            npVec = hypFrame[n, "table", drop = TRUE][gene]
            names(npVec) = "NP"
            return(c(piOut, npVec))
        } else {
            return(piOut)
        }
    })
    Design = if(fixedId) {
            rep(names(pimRes), times = vapply(piEsts, FUN.VALUE = integer(1), if(winId) NROW else length))
        } else {
            names(pimRes)
        }
    piMat = Reduce(piEsts, f = rbind)
    if(is.null(piMat))
        stop("Gene  not found!\n")
    if(!fixedId){
        weight = if(missing(weightFunction)){
                piMat[, "NP"]
        } else {
            weights = evalWeightFunction(weightFunction, newdata = data.frame(piMat[, "NP"]))
        }
        weight = weight/sum(weight, na.rm = TRUE)
        piMat = cbind(piMat, weight)
    }
    piDf = data.frame("design" = Design, piMat)
    if(!fixedId)
        rownames(piDf) = names(pimRes)
    return(piDf)
    }

