#' Plot the most significant findings for a certain PI
#'
#' Extract the most significant features for a certain PI and direction of effect,
#' and plot them using an appropriate function.
#'
#' @param hypFrame The hyperframe with the data
#' @param results The results frame
#' @param pi A character string, specifying the probabilistic index
#' @param what Which features should be detected? See details.
#' @param sigLevel The significance level
#' @param numFeats The number of features to plot
#' @param effect The name of the effect
#' @param piThreshold The threshold for PI, a minimum effect size
#' @param ... passed onto plotting functions plotCells or plotExplore
#'
#' @details The "what" argument indicates if features far from or close to
#' cell wall or centroid should be shown for pi "edge" or "centroid",
#' aggregated or regular features for "nn" and "nnCell" and colocalized or
#' antilocalized features for "nnPair" and "nnPairCell". Partial matching is allowed.
#' Defaults to small probabilistic indices: proximity, aggregation and colocalization
#'
#' @return A plot from plotCells or plotExplore, throws an error when no features meet the criteria
#' @export
#' @seealso \link{plotCells},\link{plotExplore},\link{fitLMMs}
#' @examples
#' example(fitLMMs, 'spatrans')
#' plotTopResults(hypYang, lmmModels, 'nn')
#' plotTopResults(hypYang, lmmModels, 'nn', effect = 'Intercept', what = "reg")
plotTopResults <- function(hypFrame, results, pi, effect = "Intercept",
                           what = if(pi %in% c("nn", "nnCell")){"aggregated"} else
                               if(pi %in% c("nnPair", "nnPairCell")){"colocalized"
                          } else if(pi %in% c("edge", "centroid")){"close"},
    sigLevel = 0.05, numFeats = 2, piThreshold = switch(effect, Intercept = 0.5, 0), ...) {
    stopifnot(is.hyperframe(hypFrame), is.character(what), sigLevel > 0,
              sigLevel < 1)
    what = match.arg(what, choices = c("close", "far", "regular", "aggregated", "antilocalized", "colocalized"))
    smallPI = if(what %in% c("far", "regular", "antilocalized")){
        FALSE
    } else if (what %in% c("close", "aggregated", "colocalized")){
        TRUE
        }
    pi <- match.arg(pi, choices = c("nn", "nnPair", "edge", "centroid", "nnCell",
        "nnPairCell"))
    if (is.null(Res <- results[[pi]]$results)) {
        stop("PI not present in results object!")
    }
    if (is.null(subRes <- if (effectId <- effect == "Intercept") Res$Intercept else Res$fixedEffects[[effect]])) {
        stop("Effect ", effect, " was not estimated in the linear model!")
    }
    Fun <- match.fun(if (smallPI) {"<"} else {">"})
    estId <- if (effectId) {
        !is.na(subRes[, "Estimate"]) & Fun(subRes[, "Estimate"], piThreshold)
    } else if (ncol(subRes) == 4) {
        Fun(apply(subRes[, seq_len(2)], 1, max, na.rm = TRUE), piThreshold)
    } else {
        TRUE
    }
    Feats <- rownames(subRes)[estId & subRes[, "pAdj"] < sigLevel][seq_len(numFeats)]
    if (all(is.na(Feats))) {
        stop("No significant features found for PI ", pi, " significance level ",
            sigLevel, " and what ", what, "!")
    } else {
        Feats <- Feats[!is.na(Feats)]
    }
    plotFun <- if (nnId <- pi %in% c("nn", "nnPair")) {
        plotExplore
    } else {
        plotCells
    }
    plotFun(hypFrame, Feats, titleVar = if (!effectId && !nnId) {
        effect
    }, ...)
}
