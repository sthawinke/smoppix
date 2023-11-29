#' Fit linear mixed nodels for all features of a pimRes object
#'
#' @inheritParams buildDfMM
#' @param fixedVars,randomVars Names of fixed and random variables
#' @param verbose A boolean, should the formula be printed?
#' @param resultsOnly a boolean: should only the results be returned?
#' If false, also the linear mixed models are returned
#' @param Formula A formula; if not supplied it will be constructed
#' from the fixed and random variables
#'
#' @return A fitted linear mixed model of class 'lmerTest'
#' @export
#' @import lmerTest
#' @importFrom stats formula na.omit
#'
#' @examples
#' data(Yang)
#' hypYang = buildHyperFrame(Yang, coordVars = c("x", "y"),
#' imageVars = c("day", "root", "section"))
#' #Fit a subset of features to limit computation time
#' yangPims = estPims(hypYang, pis = "nn", features = attr(hypYang, "features")[seq_len(20)])
#' #First build the weight function
#' yangObj <- addWeightFunction(yangPims, designVars = c("day", "root"))
#' fittedModels = fitLMMs(yangObj, fixedVars = "day", randomVars = "root")
fitLMMs = function(resList, pi = c("nn", "allDist", "nnPair", "allDistPair", "edge", "midpoint", "fixedpoint"),
                   fixedVars, randomVars, verbose = TRUE, resultsOnly = TRUE, Formula = NULL){
    pi = match.arg(pi)
    designVars = c(fixedVars, randomVars)
    stopifnot(all(designVars %in% resList$designVars))
    if(is.null(Formula))
        Formula = formula(formChar <- paste("pi - 0.5 ~", paste(fixedVars, sep = "+"), "+",
                                            paste0("(1|", randomVars, ")", collapse = "+")))
    if(verbose)
        cat("Fitted formula:\n", formChar)
    Control = lmerControl(## convergence checking options
                          check.conv.grad     = .makeCC("ignore", tol = 2e-3, relTol = NULL),
                          check.conv.singular = .makeCC(action = "ignore", tol = formals(isSingular)$tol),
                          check.conv.hess     = .makeCC(action = "ignore", tol = 1e-6))
    contrasts = lapply(fixedVars, function(x) "contr.sum");names(contrasts) = fixedVars
    Features = if(grepl("Pair", pi)) makePairs(attr(resList$hypFrame, "featuresEst")) else attr(resList$hypFrame, "featuresEst")
    models = lapply(Features, function(gene){
        df = buildDfMM(resList, gene = gene, pi  = pi)
        if(sum(!is.na(df$pi))<3)
            return(NULL)
        try(lmer(Formula, data = df, na.action = na.omit, weights = weight,
                 contrasts = contrasts, control = Control), silent = TRUE)
    })
    results = extractResults(models, fixedVars)
    #Effect size, standard error, p-value and adjusted p-value per mixed effect
    if(resultsOnly){
        return(list("results" = results, "resList" = resList))
    } else{
        return(list("results" = results, "models" = models, "resList" = resList))
    }
}

