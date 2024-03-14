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
#' @param effectParameter A character string, indicating which parameter to lool for when
#' effecct is provided
#' @param ... passed onto plotting functions plotCells or plotExplore
#'
#' @details The "what" argument indicates if features far from or close to
#' cell wall or centroid should be shown for pi "edge" or "centroid",
#' aggregated or regular features for "nn" and "nnCell" and colocalized or
#' antilocalized features for "nnPair" and "nnPairCell". Partial matching is allowed.
#' Defaults to small probabilistic indices: proximity, aggregation and colocalization.
#' For fixed effects, provide the name of the parameter, in combination with what. For instance,
#' what = "regular", effect = "Var1" and effectParameter = "level1" will return
#' features more regular at level1 of the variable than at baseline
#'
#' @return A plot from plotCells or plotExplore, throws an error when no features meet the criteria
#' @export
#' @seealso \link{plotCells},\link{plotExplore},\link{fitLMMs}
#' @examples
#' example(fitLMMs, 'smoppix')
#' plotTopResults(hypYang, lmmModels, 'nn')
#' plotTopResults(hypYang, lmmModels, 'nn', effect = 'Intercept', what = "reg")
#' #For the sake of illustration, set high significance level
#' plotTopResults(hypYang, lmmModels, 'nn', effect = 'day', what = "reg",
#' effectParameter ="day0", sigLevel = 0.9)
plotTopResults <- function(hypFrame, results, pi, effect = "Intercept",
                           what = if(pi %in% c("nn", "nnCell")){"aggregated"} else
                               if(pi %in% c("nnPair", "nnPairCell")){"colocalized"
                          } else if(pi %in% c("edge", "centroid")){"close"},
    sigLevel = 0.05, numFeats = 2, piThreshold = switch(effect, Intercept = 0.5, 0), effectParameter = NULL, ...) {
    stopifnot(is.hyperframe(hypFrame), is.character(what), sigLevel > 0,
              sigLevel < 1)
    pi <- match.arg(pi, choices = c("nn", "nnPair", "edge", "centroid", "nnCell",
        "nnPairCell"))
    what = match.arg(what, choices = c("close", "far", "regular", "aggregated", "antilocalized", "colocalized"))
    smallPI = if(what %in% c("far", "regular", "antilocalized")){
        FALSE
    } else if (what %in% c("close", "aggregated", "colocalized")){
        TRUE
    }
    if (is.null(Res <- results[[pi]]$results)) {
        stop("PI not present in results object!")
    }
    if (is.null(subRes <- if (interceptId <- effect == "Intercept") Res$Intercept else Res$fixedEffects[[effect]])) {
        stop("Effect ", effect, " was not estimated in the linear model!")
    }
    if(!interceptId){
        if(is.null(effectParameter)){
            stop("Provide desired parameter for fixed effects!")
        } else if(length(effectParameter)>1){
            stop("Provide effect parameter of length 1")
        } else if(!((effectParameter <- paste0(effect, effectParameter)) %in% colnames(subRes))){
            stop("Effect parameter ", effectParameter, " not found in parameters\n",
                 colnames(subRes)[seq_len(ncol(subRes)-2)])
        }
    }
    Fun <- match.fun(if (smallPI) {"<"} else {">"})
    estId <- if (interceptId) {
        !is.na(subRes[, "Estimate"]) & Fun(subRes[, "Estimate"], piThreshold)
    } else {
        !is.na(subRes[, effectParameter]) & Fun(subRes[, effectParameter], piThreshold)
    }
    Feats <- rownames(subRes)[estId & subRes[, "pAdj"] < sigLevel][seq_len(numFeats)]
    if (all(is.na(Feats))) {
        stop("No significant features found for PI ", pi, " significance level ",
            sigLevel, " and what ", what, if(!interceptId) {c("for", effectParameter)}, "!")
    } else {
        Feats <- Feats[!is.na(Feats)]
    }
    plotFun <- if (nnId <- pi %in% c("nn", "nnPair")) {
        plotExplore
    } else {
        plotCells
    }
    plotFun(hypFrame, Feats, titleVar = if (!interceptId && !nnId) {
        effect
    }, ...)
}
