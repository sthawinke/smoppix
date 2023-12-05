#' Build a data frame for a single gene and measure to prepare
#' for mixed model building
#'
#' @param resList The result of a call to estPims()
#' @param gene A character string indicating the desired gene or gene pair
#' @param pi character string indicating the desired pim as outcome
#'  in the linear mixed model
#'
#' @return A dataframe
#' @export
#' @importFrom scam predict.scam
#' @examples
#' data(Yang)
#' hypYang = buildHyperFrame(Yang, coordVars = c("x", "y"),
#' imageVars = c("day", "root", "section"))
#' #Fit a subset of features to limit computation time
#' yangPims = estPims(hypYang, pis = c("nn", "nnPair"),
#' features = c(attr(hypYang, "features")[1:12], "SmVND2", "SmPINR"))
#' #First build the weight function
#' yangObj <- addWeightFunction(yangPims, designVars = c("day", "root"))
#' #Now build the data frame for mixed model analysis
#' dfUniNN = buildDfMM(yangObj, gene = "SmAHK4e", pi  = "nn")
#' #A bivariate weight function for the pairs
#' dfBiNN = buildDfMM(yangObj, gene = "SmVND2--SmPINR", pi  = "nnPair")
#' #Example analysis with linear mixed model
#' library(lmerTest)
#' #Use sum coding for day to maintain interpretability of the intercept
#' mixedMod = lmer(pi - 0.5 ~ day + (1|root), weight = weight, data = dfBiNN,
#' contrasts = list("day" = "contr.sum"))
#' summary(mixedMod)
#'#Evidence for anticorrelation
buildDfMM = function(resList, gene,
                     pi = c("nn", "allDist", "nnPair", "allDistPair", "edge",
                            "midpoint", "fixedpoint", "nnCell", "allDistCell",
                            "nnPairCell", "allDistPairCell")){
    pi = match.arg(pi)
    stopifnot((lg <- length(gene)) %in% c(1, 2))
    #Check whether genes are present
    if(any(id <- !(sund(gene) %in% attr(resList$hypFrame, "featuresEst")))){
        stop("PIs for features\n", sund(gene)[id], "\nnot found in object")
    }
    #Establish whether pi and gene match, and call separate functions
    foo = checkAttr(resList$hypFrame, pi)
    # if(!(pi %in% c("edge", "midpoint", "fixedpoint")) && is.null(resList$Wfs[[pi]]))
    #     message("No weight function supplied, we will use approximative weigths.\n",
    #             "Consider fitting a weight function with addWeightFunction() " ,
    #             "and supplying it in the 'weightFunction' argument for improved inference.")
    df = if(pairId <- grepl(pattern = "Pair", pi)){
        if(lg == 2){
            gene = paste(gene, collapse = "--")
        } else if(!grepl("--", gene)){
            stop("Provide gene pair as character vector of length 2 or separated by '--'.")
        }
        buildDfMMBi(gene = gene, pi = pi, resList = resList)
    } else {
        if(lg !=1){
            stop("Provide a one gene for univariate PIs!")
        }
        buildDfMMUni(gene = gene, pi = pi, resList = resList)
    }
    out = data.frame(df, as.data.frame(resList$hypFrame[match(df$design,
                            rownames(resList$hypFrame)), resList$designVars]))
    attr(out, "pi") = pi
    out
}
buildDfMMUni = function(resList, gene, pi){
    piListNameInner = if(fixedId <- any(pi == c("edge", "midpoint", "fixedpoint"))) "windowDists" else "pointDists"
    #fixedId: Is the distance to a fixed point?
    winId = any(pi == c("edge", "midpoint")) #winId: is a window involved?
    piEsts = lapply(seq_along(resList$hypFrame$pimRes), function(n){
        x = resList$hypFrame$pimRes[[n]]
        if(idNull <- is.null(vec <- x[["uniPIs"]][[gene]][[piListNameInner]][pi]))
            vec = NA
        piOut = if(winId) {
            if(idNull) {
                data.frame("pi" = NA, "cell" = "NA")
                }else {
                    data.frame("pi" = unlist(vec), "cell" = rep(names(vec[[pi]]),
                    times = vapply(vec[[pi]], FUN.VALUE = integer(1), length)))}
        } else if (fixedId){
            cbind("pi" = if(idNull) vec else vec[[pi]])
        } else
            c("pi" = unname(vec))
        if(!fixedId){
            npVec = resList$hypFrame[n, "tabObs", drop = TRUE][gene]
            names(npVec) = "NP"
            return(c(piOut, npVec))
        } else {
            return(piOut)
        }
    })
    design = if(fixedId) {
            rep(rownames(resList$hypFrame),
                times = vapply(piEsts, FUN.VALUE = integer(1),
                               if(winId) NROW else length))
        } else {
            rownames(resList$hypFrame)
        }
    piMat = data.frame(Reduce(piEsts, f = rbind), "design" = design)
    if(is.null(piMat))
        stop("Gene  not found!\n")
    if(!fixedId){
        weight = evalWeightFunction(resList$Wfs[[pi]],
                                    newdata = piMat[, "NP", drop = FALSE])
        weight = weight/sum(weight, na.rm = TRUE)
        piMat = cbind(piMat, weight)
    }
    if(!fixedId)
        rownames(piMat) = rownames(resList$hypFrame)
    return(piMat)
    }
buildDfMMBi = function(resList, gene, pi){
    piEsts = t(vapply(seq_along(resList$hypFrame$pimRes),
                      FUN.VALUE = double(3), function(n){
        if(idNull <- is.null(vec <- getGp(resList$hypFrame$pimRes[[n]][["biPIs"]], gene, drop = FALSE)[,pi])){
            vec = NA
            npVec = rep(NA, 2)
        } else {
            npVec = sort(resList$hypFrame[n, "tabObs", drop = TRUE][sund(gene)])
        }
        names(npVec) = c("minP", "maxP")
        c("pi" = unname(vec), npVec)
    }))
    rownames(piEsts) = rownames(resList$hypFrame)
    if(all(is.na(piEsts[,"pi"])))
        stop("Gene pair not found!\n")
    weight = evalWeightFunction(resList$Wfs[[pi]],
                                newdata = data.frame(piEsts[, c("minP", "maxP")]))
    weight = weight/sum(weight, na.rm = TRUE)
    piMat = data.frame(piEsts, weight, design = rownames(piEsts))
    return(piMat)
}
