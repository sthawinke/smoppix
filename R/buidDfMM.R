#' Build a data frame for a single gene and measure to prepare for mixed model building
#'
#' @param pimRes The result of a call to estPims()
#' @param gene A character string indicating the desired gene or gene pair
#' @param pi character string indicating the desired pim as outcome in the linear mixed model
#' @param weightFunction The weighing function, obtained by the buildWeightFunction() function
#' @param hypFrame The hyperframe with the data
#' @param designVars Optionally, the design variables to be included in the dataframe. Defaults to all variables in the hyperframe.
#'
#' @return A dataframe
#' @export
#' @importFrom scam predict.scam
#' @examples
#' data(Yang)
#' hypYang = suppressWarnings(buildHyperFrame(Yang, coordVars = c("x", "y"),
#' designVar = c("day", "root", "section")))
#' yangPims = estPims(hypYang, pis = c("nn", "nnPair"))
#' #First build the weight function
#' wf <- buildWeightFunction(yangPims, pi = "nn", hypFrame = hypYang, designVars = c("day", "root"))
#' #Now build the data frame for mixed model analysis
#' dfUniNN = buildDfMM(yangPims, gene = "SmAHK4e", pi  = "nn", hypFrame = hypYang)
#' dfBiNN = buildDfMM(yangPims, gene = "SmAHK4e--SmAHP6b", pi  = "nnPair", hypFrame = hypYang)
#' #Example analysis
buildDfMM = function(pimRes, gene, pi = c("nn", "allDist", "nnPair", "allDistPair", "edge", "midpoint", "fixedpoint"),
                     weightFunction, hypFrame, designVars){
    #Establish whether pi and gene match, and call separate functions
    stopifnot((lg <- length(gene)) %in% c(1, 2))
    #Check whether genes are present
    if(any(id <- !(sund(gene) %in% attr(hypFrame, "features")))){
        stop("Features\n", sund(gene)[id], "\nnot found in hyperframe")
    }
    pi = match.arg(pi)
    foo = checkAttr(pimRes, pi)
    if(!missing(weightFunction))
        foo = checkAttr(weightFunction, pi)
    df = if(pairId <- grepl(pattern = "Pair", pi)){
        if(lg == 2){
            gene = paste(gene, collapse = "--")
        }
        buildDfMMBi(gene = gene, pi = pi, pimRes = pimRes, hypFrame = hypFrame, weightFunction = weightFunction)
    } else {
        if(lg !=1){
            stop("Provide a one gene for univariate PIs!")
        }
        buildDfMMUni(gene = gene, pi = pi, pimRes = pimRes, hypFrame = hypFrame, weightFunction = weightFunction)
    }
    if(missing(designVars))
        designVars = names(hypFrame)[!names(hypFrame) %in% c("ppp", "design", "tabObs", "owins")]
        #If missing take everything but the ones that are obviously no design variables
    data.frame(df, as.data.frame(hypFrame[match(df$design, names(pimRes)), designVars]))
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
            npVec = hypFrame[n, "tabObs", drop = TRUE][gene]
            names(npVec) = "NP"
            return(c(piOut, npVec))
        } else {
            return(piOut)
        }
    })
    design = if(fixedId) {
            rep(names(pimRes), times = vapply(piEsts, FUN.VALUE = integer(1), if(winId) NROW else length))
        } else {
            names(pimRes)
        }
    piMat = data.frame(Reduce(piEsts, f = rbind), "design" = design)
    if(is.null(piMat))
        stop("Gene  not found!\n")
    if(!fixedId){
        weight = if(missing(weightFunction)){
            piMat[, "NP"]
        } else {
            evalWeightFunction(weightFunction, newdata = data.frame(piMat[, "NP"]))
        }
        weight = weight/sum(weight, na.rm = TRUE)
        piMat = cbind(piMat, weight)
    }
    if(!fixedId)
        rownames(piMat) = names(pimRes)
    return(piMat)
    }
buildDfMMBi = function(pimRes, gene, hypFrame, pi, weightFunction){
    piEsts = t(vapply(seq_along(pimRes), FUN.VALUE = double(3), function(n){
        if(idNull <- is.null(vec <- getGp(pimRes[[n]][["biPIs"]], gene)[pi]))
            vec = NA
        npVec = sort(hypFrame[n, "tabObs", drop = TRUE][sund(gene)])
        names(npVec) = c("minP", "maxP")
        c("pi" = unname(vec), npVec)
    }))
    rownames(piEsts) = names(pimRes)
    if(all(is.na(piEsts[,"pi"])))
        stop("Gene pair not found!\n")
    weight = if(missing(weightFunction)){
        piEsts[,"minP"]
    } else {
        evalWeightFunction(weightFunction, newdata = data.frame(piEsts[, c("minP", "maxP")]))
    }
    weight = weight/sum(weight, na.rm = TRUE)
    piMat = data.frame(piEsts, weight, design = names(pimRes))
    return(piMat)
}
