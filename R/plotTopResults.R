#' Plot the most significant findings for a certain PI
#'
#' @param hypFrame The hyperframe
#' @param results The results frame
#' @param pi The probabilistic index
#' @param smallPI A boolean, should features with PI smaller than piThreshold be shown? See details
#' @param sigLevel The significance level
#' @param numFeats The number of features to plot
#' @param effect The name of the effect
#' @param piThreshold The threshold for PI, a minimum effect size
#' @param ... passed onto plotting functions plotCells or plotExplore
#'
#' @details If smallPI is set to TRUE, features with aggregation, colocalization and vicinity to cell boundary or centroid are shown.
#' If smallPI is set to FALSE, features with regularity, antilocalization and remoteness from cell boundary or centroid are shown.
#'
#' @return A plot from plotCells or plotExplore, throws a warning when no features meet the criteria
#' @export
#'
#' @examples
#' example(fitLMMs, "spatrans")
#' plotTopResults(hypYang, lmmModels, "nn")
#' plotTopResults(hypYang, lmmModels, "nn", effect = "Intercept")
plotTopResults = function(hypFrame, results, pi, effect = "Intercept", smallPI = TRUE, sigLevel = 0.05,
                         numFeats = 2, piThreshold = switch(effect, "Intercept" = 0.5, 0), ...){
    pi = match.arg(pi, choices = c("nn", "nnPair", "edge", "centroid", "nnCell", "nnPairCell"))
    if(is.null(Res <- results[[pi]]$results)){
        stop("PI not present in results object!")
    }
    if(is.null(subRes <- if(effectId <- effect == "Intercept") Res$Intercept else Res$fixedEffects[[effect]])){
        stop("Effect ", effect, " was not estimated in the linear model!")
    }
    Fun = match.fun(if(smallPI) "<" else ">")
    estId = if(effectId) {
            !is.na(subRes[, "Estimate"]) & Fun(subRes[, "Estimate"], piThreshold)
        } else if(ncol(subRes)==4){
            Fun(apply(subRes[, seq_len(2)], 1, max, na.rm = TRUE), piThreshold)
        } else TRUE
    Feats = rownames(subRes)[estId & subRes[, "pAdj"] < sigLevel ][seq_len(numFeats)]
    if(all(is.na(Feats))){
        stop("No significant features found for PI ", pi, " significance level ", sigLevel, " and smallPI ", smallPI, "!")
    } else {
        Feats = Feats[!is.na(Feats)]
    }
    plotFun = if(nnId <- pi %in% c("nn", "nnPair")){
        plotExplore
    } else {
        plotCells
    }
    plotFun(hypFrame, Feats, titleVar = if(!effectId && !nnId) effect, ...)
}