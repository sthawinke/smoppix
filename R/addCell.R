#' Check in which cell each event lies and add a cell marker,
#' plus include the list of regions
#'
#' @param hypFrame The hyperframe
#' @param owins the list containing a list of owins per point pattern.
#' The length of the list must match the length of the hyperframe, and the names must match.
#' Also lists of geojson objects or rois are accepted
#' @param cellTypes A dataframe of cell types and other cell-associated covariates.
#' If supplied, it must contain a variable 'cell' that is matched with the names of the owins
#' @param checkOverlap a boolean, should windows be checked for overlap?
#' @param warnOut a boolean, should warning be issued when points are not
#' contained in window
#' @param coords The names of the coordinates, if the windows are given as sets of coordinates
#' @return the modified hyperframe
#' @importFrom Rfast rowAny
#' @importFrom spatstat.geom inside.owin setmarks centroid.owin
#' @export
#' @seealso \link{buildHyperFrame}
#' @details First the different cells are checked for overlap per point pattern.
#' If no overlap is found, each event is assigned the cell that it falls into.
#' Events not belonging to any cell will trigger a warning and be assigned 'NA'.
#' Cell types and other variables are added to the marks if applicable
#' @examples
#' library(spatstat.random)
#' n <- 1e3 # number of molecules
#' ng <- 25 # number of genes
#' nfov <- 3 # Number of fields of view
#' conditions <- 3
#' # sample xy-coordinates in [0, 1]
#' x <- runif(n)
#' y <- runif(n)
#' # assign each molecule to some gene-cell pair
#' gs <- paste0('gene', seq(ng))
#' gene <- sample(gs, n, TRUE)
#' fov <- sample(nfov, n, TRUE)
#' condition <- sample(conditions, n, TRUE)
#' # construct data.frame of molecule coordinates
#' df <- data.frame(gene, x, y, fov, 'condition' = condition)
#' # A list of point patterns
#' listPPP <- tapply(seq(nrow(df)), df$fov, function(i) {
#'     ppp(x = df$x[i], y = df$y[i], marks = df[i, 'gene', drop = FALSE])
#' }, simplify = FALSE)
#' # Regions of interest (roi): Diamond in the center plus four triangles
#' w1 <- owin(poly = list(x = c(0, .5, 1, .5), y = c(.5, 0, .5, 1)))
#' w2 <- owin(poly = list(x = c(0, 0, .5), y = c(.5, 0, 0)))
#' w3 <- owin(poly = list(x = c(0, 0, .5), y = c(1, 0.5, 1)))
#' w4 <- owin(poly = list(x = c(1, 1, .5), y = c(0.5, 1, 1)))
#' w5 <- owin(poly = list(x = c(1, 1, .5), y = c(0, 0.5, 0)))
#' hypFrame <- buildHyperFrame(df,
#'     coordVars = c('x', 'y'),
#'     imageVars = c('condition', 'fov')
#' )
#' nDesignFactors <- length(unique(hypFrame$image))
#' wList <- lapply(seq_len(nDesignFactors), function(x) {
#'     list('w1' = w1, 'w2' = w2, 'w3' = w3, 'w4' = w4, 'w5' = w5)
#' })
#' names(wList) <- rownames(hypFrame) # Matching names is necessary
#' hypFrame2 <- addCell(hypFrame, wList)
#' # TO DO: accept coordinate matrices, rois geojson
addCell <- function(hypFrame, owins, cellTypes = NULL, checkOverlap = FALSE,
                    warnOut = TRUE, coords = c("x", "y")) {
    stopifnot(nrow(hypFrame) == length(owins),
              all(rownames(hypFrame) %in% names(owins)),
              is.hyperframe(hypFrame),
              is.null(cellTypes) || is.data.frame(cellTypes))
    owins = lapply(owins, convertToOwins, coords = coords)
    if (ct <- is.data.frame(cellTypes)) {
        if (!("cell" %in% names(cellTypes))) {
            stop("Variable names 'cell' must be contained in cell types dataframe")
        }
        otherCellNames <- names(cellTypes)[names(cellTypes) != "cell"]
        # Names of the other covariates
    }
    for (nn in rownames(hypFrame)) {
        ppp <- hypFrame[[nn, "ppp"]]
        NP <- npoints(ppp)
        if (any("cell" == names(marks(ppp, drop = FALSE)))) {
            stop("Cell markers already present in point pattern ", nn)
        }
        if (checkOverlap) {
            foo <- findOverlap(owins[[nn]])
        }
        idWindow <- lapply(owins[[nn]], function(rr) {
            which(inside.owin(ppp, w = rr))
        })
        if (!any((lll <- vapply(idWindow, FUN.VALUE = double(1), length))>0)) {
            stop("All points lie outside all windows for point pattern ", nn, " check your input!")
        }
        if (((nOut <- (NP- length(ul <- unlist(idWindow)))) > 0) && warnOut) {
            warning(nOut, " points lie outside all windows for point pattern ",
                    nn, " and were not assigned to a cell.\n")
        }
        if(anyDuplicated(ul)){
            warning("Some points lie in several overlapping cells! See findOverlap()")
        }
        cellOut <- rep("NA", NP)
        for(i in names(idWindow)){
            cellOut[idWindow[[i]]][cellOut[idWindow[[i]]] == "NA"] = i
            #Don't overwrite, stick to first match
        }
        hypFrame[[nn, "ppp"]] <- setmarks(hypFrame[[nn, "ppp"]], {
            m <- cbind(marks(hypFrame[[nn, "ppp"]], drop = FALSE), cell = cellOut)
            if (ct) {
                m <- cbind(m, cellTypes[match(cellOut, cellTypes$cell), otherCellNames, drop = FALSE])
            }
            m
        })
    }
    hypFrame$owins <- owins
    hypFrame$centroids <- lapply(hypFrame$owins, function(x) {
        lapply(x, function(y) centroid.owin(y, as.ppp = TRUE))
    })
    # Add centroids for all windows
    return(hypFrame)
}
