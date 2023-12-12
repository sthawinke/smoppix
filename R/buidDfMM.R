#' Build a data frame for a single gene and measure to prepare
#' for mixed model building
#'
#' @param resList The result of a call to estPims()
#' @param gene A character string indicating the desired gene or gene pair
#' @param pi character string indicating the desired pim as outcome
#'  in the linear mixed model
#' @return A dataframe
#' @export
#' @importFrom scam predict.scam
#' @examples
#' data(Yang)
#' hypYang <- buildHyperFrame(Yang,
#'     coordVars = c("x", "y"),
#'     imageVars = c("day", "root", "section")
#' )
#' # Fit a subset of features to limit computation time
#' yangPims <- estPims(hypYang,
#'     pis = c("nn", "nnPair"),
#'     features = c(attr(hypYang, "features")[1:12], "SmVND2", "SmPINR")
#' )
#' # First build the weighting function
#' yangObj <- addWeightFunction(yangPims, designVars = c("day", "root"))
#' # Now build the data frame for mixed model analysis
#' dfUniNN <- buildDfMM(yangObj, gene = "SmAHK4e", pi = "nn")
#' # A bivariate weighting function for the pairs
#' dfBiNN <- buildDfMM(yangObj, gene = "SmVND2--SmPINR", pi = "nnPair")
#' # Example analysis with linear mixed model
#' library(lmerTest)
#' # Use sum coding for day to maintain interpretability of the intercept
#' mixedMod <- lmer(pi - 0.5 ~ day + (1 | root),
#'     weight = weight, data = dfBiNN,
#'     contrasts = list("day" = "contr.sum")
#' )
#' summary(mixedMod)
#' # Evidence for anticorrelation
buildDfMM <- function(resList, gene,
                      pi = c(
                          "nn", "allDist", "nnPair", "allDistPair", "edge",
                          "midpoint", "nnCell", "allDistCell",
                          "nnPairCell", "allDistPairCell"
                      )) {
    pi <- match.arg(pi)
    stopifnot((lg <- length(gene)) %in% c(1, 2))
    # Check whether genes are present
    if (any(id <- !(sund(gene) %in% attr(resList$hypFrame, "featuresEst")))) {
        stop("PIs for features\n", sund(gene)[id], "\nnot found in object")
    }
    # Establish whether pi and gene match, and call separate functions
    foo <- checkAttr(resList$hypFrame, pi)
    df <- if (pairId <- grepl(pattern = "Pair", pi)) {
        if (lg == 2) {
            gene <- paste(gene, collapse = "--")
        } else if (!grepl("--", gene)) {
            stop("Provide gene pair as character vector of length 2 or separated by '--'.")
        }
        buildDfMMBi(gene = gene, pi = pi, resList = resList)
    } else {
        if (lg != 1) {
            stop("Provide a single gene for univariate PIs!")
        }
        buildDfMMUni(gene = gene, pi = pi, resList = resList)
    }
    out <- data.frame(df, as.data.frame(resList$hypFrame[match(
        df$design, rownames(resList$hypFrame)
    ), resList$designVars]))
    for(i in names(out)){
        if(is.character(out[[i]]))
            out[[i]] = factor(out[[i]])
    }
    attr(out, "pi") <- pi
    out
}
buildDfMMUni <- function(resList, gene, pi) {
    piListNameInner <- if (winId <- any(pi == c("edge", "midpoint", "nnCell", "allDistCell"))) "windowDists" else "pointDists"
    # winId: is a window involved?
    fixedId <- any(pi == c("edge", "midpoint"))
    # fixedId: Is the distance to a fixed point? So no weighting needed
    piEsts <- lapply(seq_along(resList$hypFrame$pimRes), function(n) {
        x <- resList$hypFrame$pimRes[[n]]
        if (idNull <- is.null(vec <- x[["uniPIs"]][[gene]][[piListNameInner]][pi])) {
            return(NULL)
        }
        piOut <- if (winId) {
            if (idNull) {
                NULL
            } else {
                Cell = rep(names(vec[[pi]]),
                           times = vapply(vec[[pi]], FUN.VALUE = integer(1), length)
                )
                cellCovars = marks(resList$hypFrame$ppp[[n]], drop = FALSE)[
                    ,attr(resList$hypFrame, "cellVars"), drop = FALSE]
                data.frame("pi" = unlist(vec), cellCovars[match(Cell, cellCovars$cell),])
            }
        } else if (fixedId) {
            cbind("pi" = if (idNull) vec else vec[[pi]])
        } else {
            c("pi" = unname(vec))
        }
        if (!fixedId) {
            # If not fixed, provide counts to allow weighting
            if (winId) {
                tabCell <- with(
                    marks(resList$hypFrame$ppp[[n]], drop = FALSE),
                    table(gene, cell)
                )
                npVec <- tabCell[gene, match(piOut[, "cell"], colnames(tabCell))]
                return(cbind(piOut, "NP" = npVec))
            } else {
                npVec <- resList$hypFrame[n, "tabObs", drop = TRUE][gene]
                names(npVec) <- "NP"
                return(c(piOut, npVec))
            }
        } else {
            return(piOut)
        }
    })
    design <- if (winId) {
        rep(rownames(resList$hypFrame),
            times = vapply(piEsts, FUN.VALUE = integer(1), NROW)
        )
    } else {
        rownames(resList$hypFrame)
    }
    piMat <- data.frame(Reduce(piEsts, f = rbind),
                        "design" = design[id <- vapply(piEsts, FUN.VALUE = TRUE, function(x) !is.null(x))])
    if (is.null(piMat)) {
        stop("Gene  not found!\n")
    }
    if (!fixedId) {
        weight <- evalWeightFunction(resList$Wfs[[pi]],
            newdata = piMat[, "NP", drop = FALSE]
        )
        weight <- weight / sum(weight, na.rm = TRUE)
        piMat <- cbind(piMat, weight)
    }
    return(piMat)
}
#' @importFrom Rfast colSort
buildDfMMBi <- function(resList, gene, pi) {
    piListNameInner <- if (winId <- grepl("Cell", pi)) "windowDists" else "pointDists"
    piEsts0 <- lapply(seq_along(resList$hypFrame$pimRes), function(n) {
        if (idNull <- is.null(vec <- getGp(resList$hypFrame$pimRes[[n]][["biPIs"]],
            gene,
            drop = FALSE
        )[[piListNameInner]][[pi]])) {
            NULL
        } else {
            if (winId) {
                tabCell <- with(
                    marks(resList$hypFrame$ppp[[n]], drop = FALSE),
                    table(gene, cell)
                )
                npVec0 <- tabCell[sund(gene), match(names(vec), colnames(tabCell)), drop = FALSE]
                npVec <- colSort(npVec0)
                colnames(npVec) <- colnames(npVec0)
                rownames(npVec) <- c("minP", "maxP")
                cbind("pi" = unlist(vec), t(npVec))
            } else {
                npVec <- sort(resList$hypFrame[n, "tabObs", drop = TRUE][sund(gene)])
                names(npVec) <- c("minP", "maxP")
                c("pi" = unname(vec), npVec)
            }
        }
    })
    piEsts <- matrix(unlist(piEsts0[id <- !vapply(piEsts0, FUN.VALUE = TRUE, is.null)]),
                     ncol = 3, byrow = TRUE, dimnames = list(NULL, c("pi", "minP", "maxP")))
    if (all(is.na(piEsts[, "pi"]))) {
        stop("Gene pair not found!\n")
    }
    rownames(piEsts) <- if (winId) {
        rep(rownames(resList$hypFrame),
            times = vapply(piEsts0, FUN.VALUE = integer(1), NROW)
        )
    } else {
        rownames(resList$hypFrame)[id]
    }
    weight <- evalWeightFunction(resList$Wfs[[pi]],
        newdata = data.frame(piEsts[, c("minP", "maxP"), drop = FALSE])
    )
    weight <- weight / sum(weight, na.rm = TRUE)
    piMat <- data.frame(piEsts, weight, design = rownames(piEsts))
    return(piMat)
}
