#' Title
#'
#' @inheritParams buildDfMM
#'
#' @return A fitted linear mixed model of class 'lmerTest'
#' @export
#'
#' @examples
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
