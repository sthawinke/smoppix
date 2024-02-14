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
#' @return The hyperframe with the estimated gradients present in it
#' @examples
#' data(Yang)
#' hypYang <- buildHyperFrame(Yang,
#'     coordVars = c('x', 'y'),
#'     imageVars = c('day', 'root', 'section')
#' )
#' yangGrads <- estGradients(hypYang)
#' Overall Gradients
#' @seealso \link{fitGradient, estGradientsSingle}
estGradients <- function(hypFrame, gradients = c("overall", if(!is.null(hypFrame$owins)) "cell"),
                         verbose = TRUE, features = getFeatures(hypFrame), silent = TRUE, ...) {
    gradients <- match.arg(gradients, several.ok = TRUE, choices = c("overall", "cell"))
    stopifnot(is.hyperframe(hypFrame), is.character(features))
    if (any(gradients == "cell") && is.null(hypFrame$owins)) {
        stop("No window provided for gradient calculation within cell. ",
            "Add it using the addCell() function")
    }
    foo = checkFeatures(hypFrame, features)
    if (verbose) {
        message("Calculating gradients for point pattern ")
    }
    hypFrame$gradients <- lapply(seq_len(nrow(hypFrame)), function(x) {
        if (verbose) {
            message(x, " of ", nrow(hypFrame), "  ")
        }
        out <- estGradientsSingle(hypFrame[[x, "ppp"]], owins = hypFrame[x, "owins", drop = TRUE],
            gradients = gradients, features = features, tabObs = hypFrame[[x, "tabObs"]], silent = silent,...)
        return(out)
    })
    list(hypFrame = hypFrame, gradients = gradients)
}
#' A wrapper function for the estimation of the different gradients, applied to individual point patterns
#'
#' @inheritParams estPisSingle
#' @param ... Passed onto fitGradient
#'
#' @return A list of data frames with estimated gradients per gene
#' \item{overall}{overall gradients}
#' \item{cell}{gradients per cell}
#' @importFrom Rdpack reprompt
#' @importFrom spatstat.geom unmark
#' @seealso \link{estGradients}
estGradientsSingle <- function(p, tabObs, gradients, owins = NULL, loopFun = "bplapply",
                         features, ...) {
    if (!length(features)) {
        return(list(overall = NULL, cell = NULL))
    }
    features <- sample(intersect(features, names(tabObs)))
    # Scramble to ensure equal calculation times in multithreading
    pSplit <- unmark(split.ppp(p, f = factor(marks(p, drop = FALSE)$gene)))
    #Unmark for ppm.ppp
    overall <- if ("overall" %in% gradients) {
        tmp = t(simplify2array(loadBalanceBplapply(features, function(gene){
            fitGradient(pSplit[[gene]], ...)
        }, loopFun = loopFun)))
        rownames(tmp) = features
        tmp
    }
    # Within cell: recurse into estGradientsSingle but now per cell
    cell <- if ("cell" %in% gradients) {
        unCells <- setdiff(unique(marks(p, drop = FALSE)$cell), "NA")
        cellGrads <- loadBalanceBplapply(unCells, function(nam) {
            pSub <- p[marks(p, drop = FALSE)$cell == nam, ]
            pSub$window = owins[[nam]]
            estGradientsSingle(gradients = "cell", p = pSub, features = features,
                               tabObs = table(marks(pSub, drop = FALSE)$gene),
                loopFun = "lapply")
        }, loopFun = loopFun)
        names(cellGrads) <- unCells
        cellGrads
    }
    list("overall" = overall, "cell" = cell)
}
