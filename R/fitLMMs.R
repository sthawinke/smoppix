#' Fit linear mixed nodels for all features of a pimRes object
#'
#' @inheritParams buildDfMM
#' @param fixedVars,randomVars Names of fixed and random variables
#'
#' @return A fitted linear mixed model of class 'lmerTest'
#' @export
#' @importFrom lmerTest lmer
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
                   fixedVars, randomVars, verbose = TRUE){
    pi = match.arg(pi)
    designVars = c(fixedVars, randomVars)
    stopifnot(all(designVars %in% resList$designVars))
    Formula = formula(paste("pi - 0.5 ~", paste(fixedVars, sep = "+"), "+", paste("(1|", T, ")", collapse = "+")))
    if(verbose)
        cat("Fitted formula:\n", Formula)
    contrasts = lapply(fixedVars, function(x) "contr.sum");names(contrasts) = fixedVars
    Features = if(grepl("Pair", pi)) makePairs(attr(resList$hypFrame, "featuresEst")) else attr(resList$hypFrame, "featuresEst")
    models = lapply(Features, function(gene){
        df = buildDfMM(resList, gene = gene, pi  = pi)
        try(lmer(Formula, data = df, na.action = na.omit, weights = weight,
                 contrasts = contrasts), silent = TRUE)
    })
}
