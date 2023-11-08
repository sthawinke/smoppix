#' Build a weight function based on the data
#'
#' @param pimRes The object containing all
#' @param pi A character string, the PI required
#' @param hypFrame The hyperframe with all the data
#' @param designVars A character vector containing all design factors (both fixed and random), that are also present as variables in hypFrame
#' @param ... Additional arguments passed on to the scam::scam function
#'
#' @return A fitted scam object, which can be used to
#' @details The scam functions fits a decreasing spline of the variance as a function of the number of observations
#' @importFrom scam scam
#' @importFrom stats formula
#' @export
#' @examples
buildWeightFunction = function(pimRes, pi = c("nn", "allDist", "nnPair", "allDistPair"), hypFrame, designVars, ...){
    if(any(pi==c("edge", "midpoint", "fixedpoint")))
        stop("Calculating weight matrices for distances to fixed points is unnecessary as they are independent.
             Simply proceed with fitting the model on the indiviual evaluations of the B-function.")
    if(!all(idMissing <- designVars %in% colnames(hypFrame))){
        stop("Design variables", designVars[idMissing], "not found in hypFrame object")
    }
    pi = match.arg(pi)
    foo = checkAttr(pimRes, pi)
    piListName = if(pairId <- grepl("Pair", pi)) "biPIs" else "uniPIs"
    piList = lapply(pimRes, function(x){
        if(pairId){
            x[["biPIs"]][pi, ]
        } else {
            vapply(x[["uniPIs"]], FUN.VALUE = double(1), function(y){
                y$pointDists[pi]
            })
        }
    })
    designVec = apply(as.data.frame(hypFrame[, designVars, drop = FALSE]), 1, paste, collapse = "_")
    ordDesign = order(designVec) #Ensure correct ordering for tapply
    features = unique(unlist(lapply(hypFrame$ppp, function(x) unique(marks(x, drop = FALSE)$gene))))
    if(pairId){
        features = apply(combn(features, 2), 2, paste, collapse = "_")
    }
    varEls = lapply(features, function(gene){
        tmp = vapply(piList[ordDesign], FUN.VALUE = double(1), function(x){
            x[gene]
        })
        quadDeps = unlist(tapply(tmp, designVec, function(x){
            if(sum(!is.na(x)) >= 2) (x-mean(x, na.rm = TRUE))^2 else rep_len(NA, length(x))
        })) #The quadratic departures from the conditional mean
        tabEntries = vapply(hypFrame$table[ordDesign], FUN.VALUE = double(if(pairId) 2 else 1), function(x){
            if(pairId) sort(x[gene]) else x[gene]
        })
        out = rbind( "quadDeps" = quadDeps, tabEntries)
        rownames(out)[-1] = if(pairId) c("minP", "maxP") else "NP"
        out
    })
    ncolMat = if(pairId) 3 else 2
    varElMat = matrix(unlist(varEls), ncol = ncolMat, byrow = TRUE,
                      dimnames = list(NULL, c("quadDeps", if(pairId) c("minP", "maxP") else "NP"))) #Faster than cbind
    scamForm = formula(paste("log(quadDeps) ~", if(pairId) {
        "s(log(maxNP), bs = 'mpd') + s(log(minNP), bs = 'mpd')"
    } else {
        "s(log(NP), bs = 'mpd')"
    }))
    scamMod = scam(scamForm, data = data.frame(varElMat), ...)
    attr(scamMod, "pi") = pi
    return(scamMod)
}
#' Evaluate the weight function and return (unnormalized) weights
#'
#' @param wf The weight function
#' @param newdata Optional, a data frame with new data
#' @return A vector of predictions
#' @export
evalWeightFunction = function(wf, newdata = NULL){
    1/exp(predict(wf, newdata = newdata))
}
