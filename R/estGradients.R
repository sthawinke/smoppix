#' Estimate gradients for single-molecule localization patterns
#' @description Estimate gradients on all point patterns of a hyperframe,
#' and integrate the results in the same hyperframe
#' @export
#' @param ... additional arguments, passed on to \link{fitGradient}.
#' @param gradients The gradients indices to be estimated
#' @param features A character vector, for which features should the
#' gradients indices be calculated?
#' @param silent A bolean, should error messages from ppm.ppp be printed?
#' @inheritParams estPis
#' @return A list with the hyperframe and the estimated gradients
#' @examples
#' data(Yang)
#' hypYang <- buildHyperFrame(Yang,
#'     coordVars = c('x', 'y'),
#'     imageVars = c('day', 'root', 'section')
#' )
#' yangGrads <- estGradients(hypYang)
#' Overall Gradients
#' @seealso \link{fitGradient, estGradientsSingle}
estGradients <- function(hypFrame, gradients = c("overall", if(!is.null(hypFrame$owins)) "cell"), fixedEffects = NULL,
                         verbose = FALSE, features = getFeatures(hypFrame), silent = TRUE, loopFun = "bplapply", ...) {
    gradients <- match.arg(gradients, several.ok = TRUE, choices = c("overall", "cell"))
    stopifnot(is.hyperframe(hypFrame), is.character(features))
    if (any(gradients == "cell") && is.null(hypFrame$owins)) {
        stop("No window provided for gradient calculation within cell. ",
            "Add it using the addCell() function")
    }
    if(!all(id <- (fixedEffects %in% getPPPvars(hypFrame)))){
        stop("Variables ", fixedEffects[!id], " not found in hyperframe.")
    }
    foo = checkFeatures(hypFrame, features)
    if (verbose) {
        message("Calculating gradients")
    }
    grads <- loadBalanceBplapply(loopFun = loopFun, features, function(gene) {
        hypFrameSub = lapply(hypFrame$ppp, function(x) x[marks(x, drop = FALSE)$gene == gene, ])
        out <- estGradientsSingle(hypFrameSub, owins = hypFrame[x, "owins", drop = TRUE],
            gradients = gradients, silent = silent, fixedEffects = fixedEffects, ...)
        return(out)
    })
    list(hypFrame = hypFrame, grads = grads)
}
#' A wrapper function for the estimation of the different gradients, applied to individual point patterns
#'
#' @inheritParams estGradients
#' @param ... Passed onto fitGradient
#'
#' @return A list of data frames with estimated gradients per gene
#' \item{overall}{overall gradients}
#' \item{cell}{gradients per cell}
#' @importFrom Rdpack reprompt
#' @importFrom spatstat.geom unmark
#' @seealso \link{estGradients}
estGradientsSingle <- function(hypFrame, gradients, fixedEffects, owins, ...) {
    #Unmark for ppm.ppp
    overall <- if ("overall" %in% gradients) {
        fitGradient(hypFrame, fixedEffects = fixedEffects,...)
    }
    # Within cell: recurse into estGradientsSingle but now per cell
    cell <- if ("cell" %in% gradients) {
        hypSub = hyperframe("ppp" = unlist(recursive = FALSE,
            lapply(rownames(hypFrame), function(rn) {
                sp = split.ppp(hypFrame$ppp[[rn]], f = marks(hypFrame$ppp[[rn]])$cell)
                sp = sd[names(sp)!="NA"]
                for(nam in names(sp)){
                    sp[[nam]]$window = owins[[rn]][[nam]]
                }
                })))
       estGradientsSingle(gradients = "overall", hypFrame = hypSub, fixedEffects = fixedEffects, ...)
    }
    cbind("overall" = overall, "cell" = cell)
}
