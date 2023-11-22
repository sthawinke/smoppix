#' Fti linear mixed nodels for all features of a pimRes object
#'
#' @inheritParams buildDfMM
#'
#' @return A fitted linear mixed model of class 'lmerTest'
#' @export
#' @importFrom lmerTest lmer
#' @importFrom stats formula
#'
#' @examples
#' #' data(Yang)
#' hypYang = buildHyperFrame(Yang, coordVars = c("x", "y"),
#' designVar = c("day", "root", "section"))
#' yangPims = estPims(hypYang, pis = "nn", features = attr(hypYang, "features)[seq_len(20)])
#' #First build the weight function
#' wf <- buildWeightFunction(yangPims, pi = "nn", hypFrame = hypYang,
#' designVars = c("day", "root"))
#' fittedModels = fitLMMs(yangPims, pi = "nn", weightFunction = wf,
#' fixedVars = "time", randomVars = "root")
fitLMMs = function(pimRes, pi = c("nn", "allDist", "nnPair", "allDistPair", "edge", "midpoint", "fixedpoint"),
                   weightFunction, hypFrame, fixedVars, randomVars){
    designVars = c(fixedVars, randomVars)
    stopifnot(all(designVars %in% names(hypFrame)))
    Formula = formula(paste("pi - 0.5 ~", paste(fixedVars, sep = "+"), "+", paste("(1|", randomVars, ")", collapse = "+")))
    contrasts = lapply(fixedVars, function(x) "contr.sum");names(contrasts) = contrasts
    models = lapply(attr(pimRes, "features"), function(gene){
        df = buildDfMM(nnObj, gene = gene, pi  = pi, hypFrame = hypFrame,
                       weightFunction = weightFunction, designVars = designVars)
        lmer(Formula, data = df, na.action = na.omit,
             weights = weight, contrasts = contrasts)
    })
}
