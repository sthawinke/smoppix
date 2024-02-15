#' Estimate gradients for single-molecule localization patterns
#' @description Estimate gradients on all point patterns of a hyperframe,
#' and integrate the results in the same hyperframe
#' @export
#' @param ... additional arguments, passed on to \link{fitGradient}.
#' @param gradients The gradients indices to be estimated
#' @param features A character vector, for which features should the
#' gradients indices be calculated?
#' @param silent A boolean, should error messages from ppm.ppp be printed?
#' @param fixedEffects,randomEffects Character vectors of fixed and random effects present in the hyperframe,
#' modifying the baseline intensity. See details.
#' @inheritParams estPis
#' @inheritParams estPisSingle
#' @note Fitting Poisson point processes is computation-intensive.
#' @details The test for existence of a gradient revolves around interaction terms
#' between x and y coordinates and image identifiers. If this interactions are
#' significant, this implies existence of gradients in the different point patterns,
#' albeit with different directions. Yet be aware that a gradient that is significant
#' for a computer may look very different from the human perspective; many
#' spatial patterns can be captured by a gradient to some extent.
#' Baseline intensity corrections for every image or cell are included by default.
#' The fixed and random effects modify the baseline intensity of the point pattern, not the gradient!
#' Random effects can lead to problems with fitting and are dissuaded.
#' @return A list with the estimated gradients
#' @examples
#' data(Yang)
#' hypYang <- buildHyperFrame(Yang,
#'     coordVars = c('x', 'y'),
#'     imageVars = c('day', 'root', 'section')
#' )
#' yangGrads <- estGradients(hypYang, features = getFeatures(hypYang)[seq_len(3)])
#' #Overall Gradients
#' @seealso \link{fitGradient}, \link{estGradientsSingle}
estGradients <- function(hypFrame, gradients = c("overall", if(!is.null(hypFrame$owins)) "cell"),
                         fixedEffects = NULL, randomEffects = NULL,
                         verbose = FALSE, features = getFeatures(hypFrame), silent = TRUE, loopFun = "bplapply", ...) {
    gradients <- match.arg(gradients, several.ok = TRUE, choices = c("overall", "cell"))
    stopifnot(is.hyperframe(hypFrame), is.character(features))
    if (any(gradients == "cell") && is.null(hypFrame$owins)) {
        stop("No window provided for gradient calculation within cell. ",
            "Add it using the addCell() function")
    }
    if(!all(id <- ((allEffects <- c(fixedEffects, randomEffects)) %in% getDesignVars(hypFrame)))){
        stop("Variables ", allEffects[!id], " not found in hyperframe.")
    }
    foo = checkFeatures(hypFrame, features)
    if (verbose) {
        message("Calculating gradients")
    }
    fixedForm = buildFormula(outcome = "ppp", fixedVars = c(fixedEffects, "x:id", "y:id", "id"),
                             randomVars = NULL)
    fixedFormSimple = buildFormula(outcome = "ppp", fixedVars = c(fixedEffects, "id"),
                             randomVars = NULL)
    randomForm = if(!is.null(randomEffects))
        buildFormula(outcome = "", fixedVars = NULL, randomVars = randomEffects)
    grads <- loadBalanceBplapply(loopFun = loopFun, features, function(gene) {
        hypFrame$ppp = lapply(hypFrame$ppp, function(x) x[marks(x, drop = FALSE)$gene == gene, ])
        out <- estGradientsSingle(hypFrame, gradients = gradients,
                                  silent = silent, fixedForm = fixedForm,
            randomForm = randomForm, fixedFormSimple = fixedFormSimple, effects = allEffects, ...)
        return(out)
    })
    names(grads) = features
    grads
}
#' A wrapper function for the estimation of the different gradients, applied to individual point patterns
#'
#' @inheritParams estGradients
#' @param fixedForm,randomForm,fixedFormSimple Formulae for fixed effects, random effects and fixed effects without
#' slopes respectively
#' @param ... Passed onto fitGradient
#'
#' @return A list of data frames with estimated gradients per gene
#' \item{overall}{overall gradients}
#' \item{cell}{gradients per cell}
#' @importFrom Rdpack reprompt
#' @importFrom spatstat.geom unmark cbind.hyperframe
#' @seealso \link{estGradients}
estGradientsSingle <- function(hypFrame, gradients, fixedForm, randomForm, fixedFormSimple, effects = NULL, ...) {
    # Within cell: recurse into estGradientsSingle but now per cell
    cell <- if ("cell" %in% gradients) {
        hypSub = hyperframe("ppp" = unlist(recursive = FALSE,
            lapply(rownames(hypFrame), function(rn) {
                sp = split.ppp(hypFrame$ppp[[rn]], f = "cell")
                sp = unmark(sp[names(sp)!="NA"])
                for(nam in names(sp)){
                    sp[[nam]]$window = hypFrame$owins[[rn]][[nam]]
                }
                return(sp)
                })))
        if(!is.null(effects)){
            hypSub = cbind(hypSub, hypFrame[rep(seq_len(nrow(hypFrame)), times = vapply(hypFrame$ppp, function(x) {
                length(unique(marks(x)$cell))
            }, FUN.VALUE = integer(1))), effects, drop = FALSE])
        }
       estGradientsSingle(gradients = "overall", hypFrame = hypSub, fixedForm = fixedForm,
                          randomForm = randomForm, fixedFormSimple = fixedFormSimple, ...)$overall
    }
    #Unmark for ppm.ppp
    overall <- if ("overall" %in% gradients) {
        hypFrame$ppp = lapply(hypFrame$ppp, unmark)
        fitGradient(hypFrame, fixedForm = fixedForm, randomForm = randomForm,
                    fixedFormSimple = fixedFormSimple,...)
    }
    list("overall" = overall, "cell" = cell)
}
