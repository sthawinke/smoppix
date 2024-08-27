#' Build a data frame for a certain gene and PI, in preparation for mixed model building.
#'
#' @param obj A results object. For distances to fixed objects, the result of a call to \link{estPis}();
#' for nearest neighbour distances, the result of a call to \link{addWeightFunction}
#' @param gene A character string indicating the desired gene or gene pair (gene separated by double hyphens)
#' @param pi character string indicating the desired PI
#' @param piMat A data frame. Will be constructed if not provided, for internal use
#' @param moransI A boolean, should Moran's I be calculated
#' in the linear mixed model
#' @param weightMats List of weight matrices for Moran's I calculation
#' @param pppDf Dataframe of point pattern-wise variables. It is precalculated
#' in fitLMMsSingle for speed, but will be newly constructed when not provided
#' @inheritParams fitLMMs
#' @return A dataframe with estimated PIs and covariates
#' @export
#' @importFrom scam predict.scam
#' @seealso \link{addWeightFunction}, \link{buildMoransIDataFrame}
#' @examples
#' example(addWeightFunction, "smoppix")
#' dfUniNN <- buildDataFrame(yangObj, gene = "SmVND2", pi = "nn")
#' # Example analysis with linear mixed model
#' library(lmerTest)
#' mixedMod <- lmer(pi - 0.5 ~ day + (1 | root),
#'     weight = weight, data = dfUniNN,
#'     contrasts = list("day" = "contr.sum")
#' )
#' summary(mixedMod)
#' # Evidence for aggregation
buildDataFrame <- function(obj, gene, pi = c("nn", "nnPair", "edge", "centroid",
                               "nnCell", "nnPairCell"
                           ), piMat, moransI = FALSE, numNNs = 8, weightMats, pppDf) {
    pi <- match.arg(pi)
    if (missing(pppDf)) {
        pppDf <- as.data.frame(obj$hypFrame[, c("image", getPPPvars(obj)), drop = FALSE])
    }
    if (missing(piMat)) {
        stopifnot((lg <- length(gene)) %in% c(1, 2))
        foo <- checkPi(obj, pi)
        if (!(pi %in% c("edge", "centroid")) && is.null(obj$Wfs[[pi]])) {
            stop("No weight function added yet, run addWeightFunction first!")
        }
        # Establish whether pi and gene match, and call separate functions
        windowId <- pi %in% c("edge", "centroid")
        cellId <- grepl("Cell", pi)
        if (cellId) {
            piSub <- sub("Cell", "", pi)
        }
        if (pairId <- grepl("Pair", pi)) {
            if (lg == 2) {
                gene <- paste(gene, collapse = "--")
            } else if (!grepl("--", gene)) {
                stop("Provide gene pair as character vector of length 2", " or separated by '--'.")
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
        piListNameInner <- if (windowId) {
            "windowDists"
        } else if (cellId) {
            "withinCellDists"
        } else {
            "pointDists"
        }
        if (cellId || windowId) {
            eventVars <- getEventVars(obj)
        }
        piDfs <- lapply(seq_len(nrow(obj$hypFrame)), function(n) {
            nList <- obj$hypFrame[n, , drop = TRUE]
            class(nList) <- "list" # Simplify things
            # Extract pi, and add counts and covariates
            if (cellId || windowId) {
                Marks <- marks(nList$ppp, drop = FALSE)
            }
            df <- if (windowId) {
                piEst <- if (is.null(vec <- getGp(nList$pimRes[[piListNameInner]], gene)[[pi]])) {
                    NULL
                } else {
                    vec
                }
                dfWin <- data.frame(pi = unlist(piEst), cell = rep(names(piEst),
                    times = vapply(piEst, FUN.VALUE = double(1), length)
                ))
                if (length(cellVars <- setdiff(eventVars, c("gene", "cell"))) &&
                    NROW(dfWin)) {
                    mat <- Marks[match(dfWin$cell, Marks$cell), cellVars, drop = FALSE]
                    colnames(mat) <- cellVars
                    dfWin <- cbind(dfWin, mat)
                }
                return(dfWin)
            } else if (cellId) {
                piEst <- vapply(nList$pimRes[[piListNameInner]], FUN.VALUE = double(1), function(y) {
                    if (is.null(vec <- getGp(y[[piSub]], gene))) {
                        NA
                    } else {
                        vec
                    }
                })
                cellCovars <- Marks[, eventVars, drop = FALSE]
                cellCovars <- cellCovars[match(names(piEst), cellCovars$cell), setdiff(
                    colnames(cellCovars),
                    "gene"
                ), drop = FALSE]
                tabCell <- table(Marks$gene, Marks$cell)
                npVec <- if (all(geneSplit %in% rownames(tabCell))) {
                    tmp <- tabCell[geneSplit, match(names(piEst), colnames(tabCell)),
                        drop = FALSE
                    ]
                    t(matrix(tmp, ncol = ncol(tmp), dimnames = dimnames(tmp)))
                } else {
                    matrix(NA, ncol = if (pairId) {
                        2
                    } else {
                        1
                    }, nrow = length(piEst))
                }
                colnames(npVec) <- if (pairId) {
                    c("minP", "maxP")
                } else {
                    "NP"
                }
                data.frame(pi = piEst, cellCovars, npVec)
            } else {
                piEst <- if (is.null(vec <- getGp(nList$pimRes[[piListNameInner]][[pi]], gene))) {
                    NA
                } else {
                    vec
                }
                npVec <- if (all(geneSplit %in% names(to <- nList$tabObs))) {
                    to[geneSplit]
                } else {
                    NA
                }
                if (pairId) {
                    cbind(pi = piEst, minP = min(npVec), maxP = max(npVec))
                } else {
                    cbind(pi = piEst, NP = npVec)
                }
            }
        })
        if (all((Times <- vapply(piDfs, FUN.VALUE = integer(1), NROW)) == 0)) {
            return(NULL)
        }
        image <- rep(rownames(obj$hypFrame), times = Times)
        piDfsMat <- Reduce(piDfs, f = rbind)
        piMat <- data.frame(piDfsMat, pppDf[image, ,drop = FALSE])
        if (!windowId && !moransI) {
            # Add weights
            weight <- evalWeightFunction(obj$Wfs[[pi]], newdata = piMat[, if (grepl(
                "Pair",
                pi
            )) {
                c("minP", "maxP")
            } else {
                "NP"
            }, drop = FALSE])
            weight <- weight / sum(weight, na.rm = TRUE)
            piMat <- cbind(piMat, weight = weight)
        }
        rownames(piMat) <- NULL
    }
    if (moransI) {
        if (missing(weightMats)) {
            weightMats <- lapply(getHypFrame(obj)$centroids, buildMoransIWeightMat,
                numNNs = numNNs
            )
            names(weightMats) <- getHypFrame(obj)$image
        }
        piMat <- buildMoransIDataFrame(piMat = piMat, pi = pi, weightMats = weightMats)
    }
    return(piMat)
}
