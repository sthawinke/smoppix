#' Plot the most significant findings for a certain PI
#'
#' @param hypFrame The hyperframe
#' @param results The results frame
#' @param pi The probabilistic index
#' @param smallPI A boolean, should features with PI smaller than piThreshold be shown? See details
#' @param sigLevel The significance level
#' @param numFeats The number of features to plot
#' @param ... passed onto plotting functions plotCells or plotExplore
#' @details If smallPI is set to TRUE, features with aggregation, colocalization and vicinity to cell boundary or centroid are shown.
#' If smallPI is set to FALSE, features with regularity, antilocalization and remoteness from cell boundary or centroid are shown.
#'
#' @return A plot from plotCells or plotExplore, throws a warning when no features meet the criteria
#' @export
#'
#' @examples
#' example(fitLMMs, "spatrans")
#' plotTopResult(hypYang, yangObj, "nn")
plotTopResult = function(hypFrame, results, pi, smallPI = TRUE, sigLevel = 0.05,
                         numFeats = 2, piThreshold = 0.5, ...){
    Res = results[[pi]]
    pi = match.arg(pi, choices = c("nn", "nnPair", "edge", "midpoint", "nnCell", "nnPairCell"))
    Fun = if(smallPI) "<" else ">"
    Feats = rownames(Res)[Fun(Res[, "Estimate"], piThreshold) & Res[, "pAdj"] < sigLevel ][seq_len(numFeats)]
    if(length(Feats) == 0){
        stop("No feature match these criteria")
    }
    plotFun = if(pi %in% c("nn", "nnPair")){
        plotExplore
    } else {
        plotCells
    }
    plotFun(hypFrame, features, ...)
}
