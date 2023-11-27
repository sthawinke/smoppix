#' Build a weight function based on the data
#'
#' @param pimRes The object containing all
#' @param pi A character string, the PI required
#' @param hypFrame The hyperframe with all the data
#' @param designVars A character vector containing all design factors (both fixed and random), that are also present as variables in hypFrame
#' @param ... Additional arguments passed on to the scam::scam function
#' @param maxObs,maxFeatures The maximum number of observations respectively features for fitting the weight function. See details
#' @details For computational and memory reasons, for large datasets the trend fitting
#' is restricted to a subset of the data through the maxObs and maxFeatures parameters.
#'
#' @return A fitted scam object, which can be used to
#' @details The scam functions fits a decreasing spline of the variance as a function of the number of observations
#' @importFrom scam scam
#' @importFrom stats formula
#' @export
#' @seealso \link{buildDfMM}
#' @examples
#' data(Yang)
#' hypYang = buildHyperFrame(Yang, coordVars = c("x", "y"),
#' designVar = c("day", "root", "section"))
#' #Fit a subset of features to limit computation time
#' yangPims = estPims(hypYang, pis = c("nn", "nnPair"), features = attr(hypYang, "features")[1:12])
#' #First Build the weight function
#' wf <- adddWeightFunction(yangPims, pi = "nn", hypFrame = hypYang,
#' designVars = c("day", "root"))
adddWeightFunction = function(pimRes, pi = c("nn", "allDist", "nnPair", "allDistPair"),
                               hypFrame, designVars, maxObs = 1e5, maxFeatures = 5e2,...){
    if(any(pi==c("edge", "midpoint", "fixedpoint")))
        stop("Calculating weight matrices for distances to fixed points is unnecessary as they are independent.
             Simply proceed with fitting the model on the indiviual evaluations of the B-function.")
    if(any(idMissing <- !(designVars %in% colnames(hypFrame)))){
        stop("Design variables\n", designVars[idMissing], "\nnot found in hypFrame object")
    }
    if(is.null(hypFrame$pimRes)){
        stop("No pims found in the hyperframe. First estimate them using the estPims() function.")
    }
    pi = match.arg(pi)
    foo = checkAttr(pimRes, pi)
    piListName = if(pairId <- grepl("Pair", pi)) "biPIs" else "uniPIs"
    piList = lapply(pimRes, function(x){
        if(pairId){
            x[["biPIs"]][,pi ,drop = FALSE]
        } else {
            vapply(x[["uniPIs"]], FUN.VALUE = double(1), function(y){
                y$pointDists[pi]
            })
        }
    })
    designVec = apply(as.data.frame(hypFrame[, designVars, drop = FALSE]), 1, paste, collapse = "_")
    ordDesign = order(designVec) #Ensure correct ordering for tapply
    features = attr(pimRes, "features")
    if(pairId){
        features = apply(combn(features, 2), 2, paste, collapse = "--")
    }
    varEls = lapply(features, function(gene){
        tmp = vapply(piList[ordDesign], FUN.VALUE = double(1), function(x){
            if(is.null(baa <- getGp(x,gene))) NA else baa
        })
        quadDeps = unlist(tapply(tmp, designVec, function(x){
            if(sum(!is.na(x)) >= 2) (x-mean(x, na.rm = TRUE))^2 else rep_len(NA, length(x))
        })) #The quadratic departures from the conditional mean
        if(pairId){
            gene = sund(gene)
        }
        tabEntries = vapply(hypFrame$tabObs[ordDesign], FUN.VALUE = double(if(pairId) 2 else 1), function(x){
            if(pairId) {
                if(all(gene %in% names(x))) sort(x[gene]) else rep(NA, 2)
            } else x[gene]
        })
        out = rbind("quadDeps" = quadDeps, tabEntries)
        rownames(out)[-1] = if(pairId) c("minP", "maxP") else "NP"
        out
    })
    varElMat = matrix(unlist(varEls), ncol = if(pairId) 3 else 2, byrow = TRUE,
                      dimnames = list(NULL, c("quadDeps", if(pairId) c("minP", "maxP") else "NP"))) #Faster than cbind
    varElMat = varElMat[!is.na(varElMat[,"quadDeps"]) & varElMat[,"quadDeps"] != 0,]
    if(nrow(varElMat) > maxObs){
        varElMat = varElMat[sample(nrow(varElMat), maxObs),]
    }
    scamForm = formula(paste("log(quadDeps) ~", if(pairId) {
        "s(log(maxP), bs = 'mpd') + s(log(minP), bs = 'mpd')"
    } else {
        "s(log(NP), bs = 'mpd')"
    }))
    scamMod = scam(scamForm, data = data.frame(varElMat), ...)
    attr(scamMod, "pis") = pi
    return(scamMod)
}

