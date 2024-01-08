#' Build a data frame for a single gene and measure to prepare
#' for mixed model building
#'
#' @param obj The result of a call to buildWeightFunction()
#' @param gene A character string indicating the desired gene or gene pair
#' @param pi character string indicating the desired pim as outcome
#'  in the linear mixed model
#' @return A dataframe
#' @export
#' @importFrom scam predict.scam
#' @examples
#' data(Yang)
#' hypYang <- buildHyperFrame(Yang,
#'     coordVars = c('x', 'y'),
#'     imageVars = c('day', 'root', 'section')
#' )
#' # Fit a subset of features to limit computation time
#' yangPims <- estPims(hypYang[c(seq_len(5), seq(25, 29)),], pis = 'nn',nPointsAll = 1e3,
#' features = c('SmVND2', 'SmPINR'))
#' # First build the weighting function
#' yangPims <- addWeightFunction(yangPims, designVars = c('day', 'root'))
#' # Now build the data frame for mixed model analysis
#' dfUniNN <- buildDataFrame(yangPims, gene = 'SmVND2', pi = 'nn')
#' # Example analysis with linear mixed model
#' library(lmerTest)
#' # Use sum coding for day to maintain interpretability of the intercept
#' mixedMod <- lmer(pi - 0.5 ~ day + (1 | root),
#'     weight = weight, data = dfUniNN,
#'     contrasts = list('day' = 'contr.sum')
#' )
#' summary(mixedMod)
#' # Evidence for aggregation
buildDataFrame <- function(obj, gene, pi = c("nn", "nnPair", "edge", "midpoint",
                                             "nnCell", "nnPairCell")) {
    pi <- match.arg(pi)
    stopifnot((lg <- length(gene)) %in% c(1, 2))
    foo <- checkPi(obj, pi)
    if (!(pi %in% c("edge", "midpoint")) && is.null(obj$Wfs[[pi]])) {
        stop("No weight function added yet, run addWeightFunction first!")
    }
    # Establish whether pi and gene match, and call separate functions
    windowId <- pi %in% c("edge", "midpoint")
    cellId <- any(grepl("Cell", pi))
    if (cellId) {
        piSub <- sub("Cell", "", pi)
    }
    if (pairId <- grepl("Pair", pi)) {
        if (lg == 2) {
            gene <- paste(gene, collapse = "--")
        } else if (!grepl("--", gene)) {
            stop("Provide gene pair as character vector of length 2",
                 " or separated by '--'.")
        }
    } else {
        if (lg != 1 || grepl("--", gene)) {
            stop("Provide a single gene for univariate PIs!")
        }
    }
    geneSplit <- sund(gene)
    # Check whether genes are present
    if (any(id <- !(geneSplit %in% getFeatures(obj)))) {
        stop("PIs for features\n", geneSplit[id], "\nnot found in object")
    }
    piListNameInner <- if (windowId)
        "windowDists" else if (cellId)
        "withinCellDists" else "pointDists"
    piDfs <- lapply(seq_along(obj$hypFrame$pimRes), function(n) {
        x <- obj$hypFrame$pimRes[[n]]
        # Extract pi, and add counts and covariates
        df <- if (windowId) {
            piEst <- if(is.null(vec <- getGp(x[[piListNameInner]], gene)[[pi]]))
                NA else vec
            data.frame(pi = unlist(piEst), cell = rep(names(piEst),
                times = vapply(piEst, FUN.VALUE = double(1), length)))
        } else if (cellId) {
            piEst <- vapply(x[[piListNameInner]], FUN.VALUE = double(1), function(y){
                if(is.null(vec <- getGp(y[[piSub]], gene))) NA else vec
            })
            cellCovars <- marks(obj$hypFrame$ppp[[n]], drop = FALSE)[,
                                            getEventVars(obj), drop = FALSE]
            cellCovars <- cellCovars[match(names(piEst), cellCovars$cell),
                                     setdiff(colnames(cellCovars), "gene")]
            tabCell <- table(marks(obj$hypFrame$ppp[[n]], drop = FALSE)$gene,
                             marks(obj$hypFrame$ppp[[n]], drop = FALSE)$cell)
            npVec <- if (all(geneSplit %in% rownames(tabCell))) {
            tmp <- tabCell[geneSplit, match(names(piEst), colnames(tabCell)), drop = FALSE]
                t(matrix(tmp, ncol = ncol(tmp), dimnames = dimnames(tmp)))
                    } else matrix(NA, ncol = if (pairId) 2 else 1, nrow = length(piEst))
             colnames(npVec) <- if (pairId)
                       c("minP", "maxP") else "NP"
            data.frame(pi = piEst, cellCovars, npVec)
        } else {
            piEst <- if(is.null(vec <- getGp(x[[piListNameInner]][[pi]], gene)))
                NA else vec
            npVec <- if (all(geneSplit %in% names(obj$hypFrame[n, "tabObs", drop = TRUE]))){
                obj$hypFrame[n, "tabObs", drop = TRUE][geneSplit]
            } else NA
            if (pairId) {
                data.frame(pi = piEst, minP = min(npVec), maxP = max(npVec))
            } else {
                data.frame(pi = piEst, NP = npVec)
            }
            data.frame(pi = piEst, NP = npVec)
        }
    })
    if (all((Times <- vapply(piDfs, FUN.VALUE = integer(1), NROW)) == 0)) {
        stop("Gene not found!\n")
    }
    design <- rep(rownames(obj$hypFrame), times = Times)
    piMat <- data.frame(Reduce(piDfs, f = rbind), obj$hypFrame[design,
        getPPPvars(obj)], design = design)
    if (!windowId) {
        # Add weights
        weight <- evalWeightFunction(obj$Wfs[[pi]],
            newdata = piMat[, if(grepl("Pair", pi)) c("minP", "maxP") else "NP",
                            drop = FALSE])
        weight <- weight/sum(weight, na.rm = TRUE)
        piMat <- cbind(piMat, weight = weight)
    }
    rownames(piMat) <- NULL
    return(piMat)
}
