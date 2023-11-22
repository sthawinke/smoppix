#' Title
#'
#' @inheritParams buildDfMM
#'
#' @return A fitted linear mixed model of class 'lmerTest'
#' @export
#'
#' @examples
#' #' data(Yang)
#' hypYang = buildHyperFrame(Yang, coordVars = c("x", "y"),
#' designVar = c("day", "root", "section"))
#' yangPims = estPims(hypYang, pis = "nn")
#' #First build the weight function
#' wf <- buildWeightFunction(yangPims, pi = "nn", hypFrame = hypYang,
#' designVars = c("day", "root"))
#' fittedModels = fitLMMs(yangPims, pi = "nn", weightFunction = wf, fiexedVars = "time", randomVars = "root")
fitLMMs = function(pimRes, pi = c("nn", "allDist", "nnPair", "allDistPair", "edge", "midpoint", "fixedpoint"),
                   weightFunction, hypFrame, fixedVars, randomVars){
    designVars = c(fixedVars, randomVars)
    Fromula = formula(paste("pi ~", paste(fixedVars, sep = "+"), "+", paste("(1|", randomVars, ")", sep = "+")))
    models = lapply(attr(pimRes, "features"), function(gene){
        df = buildDfMM(nnObj, gene = gene, pi  = pi, hypFrame = hypFrame,
                       weightFunction = weightFunction, designVars = designVars)
        lmer(pi ~ day + (1|root), data = df, na.action = na.omit, offset = Offset,
             weights = weights, contrasts = list("day" = "contr.sum"))
    })
}
