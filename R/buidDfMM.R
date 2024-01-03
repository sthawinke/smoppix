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
#' yangPims <- estPims(hypYang[c(seq_len(5), seq(25, 29)),], pis = 'nn',
#' features = c('SmVND2', 'SmPINR'))
#' # First build the weighting function
#' yangPims <- addWeightFunction(yangPims, designVars = c('day', 'root'))
#' # Now build the data frame for mixed model analysis
#' dfUniNN <- buildDfMM(yangPims, gene = 'SmVND2', pi = 'nn')
#' # Example analysis with linear mixed model
#' library(lmerTest)
#' # Use sum coding for day to maintain interpretability of the intercept
#' mixedMod <- lmer(pi - 0.5 ~ day + (1 | root),
#'     weight = weight, data = dfUniNN,
#'     contrasts = list('day' = 'contr.sum')
#' )
#' summary(mixedMod)
#' # Evidence for aggregation
buildDfMM <- function(obj, gene, pi = c("nn", "nnPair", "edge", "midpoint", "nnCell",
    "nnPairCell")) {
    pi <- match.arg(pi)
    stopifnot((lg <- length(gene)) %in% c(1, 2))
    if (is.null(obj$Wfs)) {
        stop("No weight function added yet, run addWeightFunction first!")
    }
    # Check whether genes are present
    if (any(id <- !(sund(gene) %in% getFeatures(obj)))) {
        stop("PIs for features\n", sund(gene)[id], "\nnot found in object")
    }
    # Establish whether pi and gene match, and call separate functions
    foo <- checkPi(obj, pi)
    df <- if (pairId <- grepl(pattern = "Pair", pi)) {
        if (lg == 2) {
            gene <- paste(gene, collapse = "--")
        } else if (!grepl("--", gene)) {
            stop("Provide gene pair as character vector of length 2 or separated by '--'.")
        }
        buildDfMMBi(gene = gene, pi = pi, obj = obj)
    } else {
        if (lg != 1) {
            stop("Provide a single gene for univariate PIs!")
        }
        buildDfMMUni(gene = gene, pi = pi, obj = obj)
    }
    out <- data.frame(df, as.data.frame(obj$hypFrame[match(df$design, rownames(obj$hypFrame)), obj$designVars]))
    for (i in names(out)) {
        if (is.character(out[[i]])) {
            out[[i]] <- factor(out[[i]])
        }
    }
    attr(out, "pi") <- pi
    out
}
buildDfMMUni <- function(obj, gene, pi) {
    windowId <- pi %in% c("edge", "midpoint")
    cellId <- any(grepl("Cell", pi))
    if(cellId){piSub <- sub("Cell", "", pi)}
    piListNameInner = if(windowId) "windowDists" else if(cellId) "withinCellDists" else "pointDists"
    piDfs <- lapply(seq_along(obj$hypFrame$pimRes), function(n) {
        x <- obj$hypFrame$pimRes[[n]]
        #Extract pi, and add counts and covariates
        df = if (windowId) {
            piEst = if(gene %in% names(x[[piListNameInner]])) x[[piListNameInner]][[gene]][[pi]] else NA
            data.frame("pi" = unlist(piEst),
                "cell" = rep(names(piEst),
                        times = vapply(piEst, FUN.VALUE = double(1), length)))
        } else if (cellId) {
            piEst = vapply(x[[piListNameInner]], FUN.VALUE = double(1), function(y) {
                if(gene %in% names(y[[piSub]])) y[[piSub]][[gene]] else NA
                })
            cellCovars <- marks(obj$hypFrame$ppp[[n]], drop = FALSE)[, getEventVars(obj), drop = FALSE]
            cellCovars = cellCovars[match(names(piEst), cellCovars$cell),
                       setdiff(colnames(cellCovars), "gene")]
            tabCell <- table(marks(obj$hypFrame$ppp[[n]], drop = FALSE)$gene, marks(obj$hypFrame$ppp[[n]], drop = FALSE)$cell)
            npVec <- tabCell[gene, match(names(piEst), colnames(tabCell))]
            data.frame(pi = piEst, cellCovars, "NP" = npVec[cellCovars$cell])
        } else {
            piEst = if(gene %in% names(x[[piListNameInner]][[pi]])) x[[piListNameInner]][[pi]] else NA
            npVec <- obj$hypFrame[n, "tabObs", drop = TRUE][gene]
            data.frame(piEst, "NP" = npVec)
        }
    })
    if (all((Times <- vapply(piEsts, FUN.VALUE = integer(1), NROW))==0)){
        stop("Gene  not found!\n")
    }
    design <- rep(rownames(obj$hypFrame), times = Times)
    piMat <- data.frame(Reduce(piEsts, f = rbind), obj$hypFrame[design, getPPPvars(obj)])
    if (!windowId) {#Add weights
        piMat = addWeights(piMat, pi, ibj)
    }
    return(piMat)
}
#' @importFrom Rfast colSort
buildDfMMBi <- function(obj, gene, pi) {
    piListNameInner <- if (winId <- grepl("Cell", pi))
        "windowDists" else "pointDists"
    piEsts0 <- lapply(seq_along(obj$hypFrame$pimRes), function(n) {
        if (idNull <- is.null(vec <- getGp(obj$hypFrame$pimRes[[n]][["biPIs"]], gene, drop = FALSE)[[piListNameInner]][[pi]])) {
            NULL
        } else {
            if (winId) {
                tabCell <- table(marks(obj$hypFrame$ppp[[n]], drop = FALSE)$gene, marks(obj$hypFrame$ppp[[n]], drop = FALSE)$cell)
                npVec0 <- tabCell[sund(gene), match(names(vec), colnames(tabCell)), drop = FALSE]
                npVec <- colSort(npVec0)
                colnames(npVec) <- colnames(npVec0)
                rownames(npVec) <- c("minP", "maxP")
                cbind(pi = unlist(vec), t(npVec))
            } else {
                npVec <- sort(obj$hypFrame[n, "tabObs", drop = TRUE][sund(gene)])
                names(npVec) <- c("minP", "maxP")
                c(pi = unname(vec), npVec)
            }
        }
    })
    piEsts <- matrix(unlist(piEsts0[id <- !vapply(piEsts0, FUN.VALUE = TRUE, is.null)]), ncol = 3, byrow = TRUE, dimnames = list(NULL,
        c("pi", "minP", "maxP")))
    rownames(piEsts) <- if (winId) {
        rep(rownames(obj$hypFrame), times = vapply(piEsts0, FUN.VALUE = integer(1), NROW))
    } else {
        rownames(obj$hypFrame)[id]
    }
    weight <- evalWeightFunction(obj$Wfs[[pi]], newdata = data.frame(piEsts[, c("minP", "maxP"), drop = FALSE]))
    weight <- weight/sum(weight, na.rm = TRUE)
    piMat <- data.frame(piEsts, weight, design = rownames(piEsts))
    return(piMat)
}
