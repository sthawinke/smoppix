#' Add cell boundaries and event-wise cell identifiers to a hyperframe.
#' @description Add the list of the cells and their centroids in the hyperframe,
#' check in which cell each event lies and add a cell marker.
#'
#' @param hypFrame The hyperframe
#' @param owins the list containing a list of owins per point pattern.
#' The length of the list must match the length of the hyperframe, and the names must match.
#' Also lists of geojson objects, coordinate matrices or rois are accepted, see details.
#' @param cellTypes A dataframe of cell types and other cell-associated covariates.
#' If supplied, it must contain a variable 'cell' that is matched with the names of the owins
#' @param findOverlappingOwins a boolean, should windows be checked for overlap?
#' @param warnOut a boolean, should warning be issued when points are not
#' contained in window?
#' @param coords The names of the coordinates, if the windows are given as sets of coordinates.
#' @param verbose A boolean, should verbose output be printed?
#' @param addCellMarkers A boolean, should cell identities be added? Set this to
#'  FALSE if cell identifiers are already present in the data, and you only want
#'  to add windows and centroids.
#' @param ... Further arguments passed onto \link{convertToOwins}
#' @return The hyperframe with cell labels added in the marks of the point patterns
#' @importFrom spatstat.geom inside.owin marks centroid.owin marks<-
#' @export
#' @seealso \link{buildHyperFrame}, \link{convertToOwins}
#' @details First the different cells are checked for overlap per point pattern.
#' If no overlap is found, each event is assigned the cell that it falls into.
#' Events not belonging to any cell will trigger a warning and be assigned 'NA'.
#' Cell types and other variables are added to the marks if applicable.
#' This function employs multithreading through the BiocParallel package for
#' converting windows to owins as well as for finding overlap between the windows
#' in the findOverlap() function when checkOverlap is TRUE.
#' If this leads to excessive memory usage and crashes, try serial processing by
#' setting register(SerialParam()).
#' Different formats of windows are allowed, if the corresponding packages are installed.
#' A dataframe of coordinates or a list of spatstat.geom owins is always allowed, as necessary packages are required by smoppix.
#' A 'SpatialPolygonsDataFrame' object is allowed if the polycub package is installed,
#' and a lsit of 'ijroi' object or a single 'ijzip' of the 'RImageJROI' package is installed.
#' @note By default, there is no checking for overlap between windows.
#' Events are assigned to the first window in which they fall
#' Do check your input or set checkOverlap to TRUE, even when this make take time.
#' @examples
#' library(spatstat.random)
#' set.seed(54321)
#' n <- 1e3 # number of molecules
#' ng <- 25 # number of genes
#' nfov <- 3 # Number of fields of view
#' conditions <- 3
#' # sample xy-coordinates in [0, 1]
#' x <- runif(n)
#' y <- runif(n)
#' # assign each molecule to some gene-cell pair
#' gs <- paste0("gene", seq(ng))
#' gene <- sample(gs, n, TRUE)
#' fov <- sample(nfov, n, TRUE)
#' condition <- sample(conditions, n, TRUE)
#' # construct data.frame of molecule coordinates
#' df <- data.frame(gene, x, y, fov, "condition" = condition)
#' # A list of point patterns
#' listPPP <- tapply(seq(nrow(df)), df$fov, function(i) {
#'     ppp(x = df$x[i], y = df$y[i], marks = df[i, "gene", drop = FALSE])
#' }, simplify = FALSE)
#' # Regions of interest (roi): Diamond in the center plus four triangles
#' w1 <- owin(poly = list(x = c(0, .5, 1, .5), y = c(.5, 0, .5, 1)))
#' w2 <- owin(poly = list(x = c(0, 0, .5), y = c(.5, 0, 0)))
#' w3 <- owin(poly = list(x = c(0, 0, .5), y = c(1, 0.5, 1)))
#' w4 <- owin(poly = list(x = c(1, 1, .5), y = c(0.5, 1, 1)))
#' w5 <- owin(poly = list(x = c(1, 1, .5), y = c(0, 0.5, 0)))
#' hypFrame <- buildHyperFrame(df,
#'     coordVars = c("x", "y"),
#'     imageVars = c("condition", "fov")
#' )
#' nDesignFactors <- length(unique(hypFrame$image))
#' wList <- lapply(seq_len(nDesignFactors), function(x) {
#'     list("w1" = w1, "w2" = w2, "w3" = w3, "w4" = w4, "w5" = w5)
#' })
#' names(wList) <- rownames(hypFrame) # Matching names is necessary
#' hypFrame2 <- addCell(hypFrame, wList)
addCell <- function(hypFrame,
                    owins,
                    cellTypes = NULL,
                    findOverlappingOwins = FALSE,
                    warnOut = TRUE,
                    coords = c("x", "y"),
                    verbose = TRUE,
                    addCellMarkers = TRUE,
                    ...) {
    stopifnot(
        nrow(hypFrame) == length(owins),
        all(rownames(hypFrame) %in% names(owins)),
        is.hyperframe(hypFrame),
        is.null(cellTypes) || is.data.frame(cellTypes)
    )
    if (verbose) {
        message("Converting windows to spatstat owins")
    }
    # Convert different types of windows to owins
    owins <-
        loadBalanceBplapply(Nam <- names(owins), function(nam) {
            convertToOwins(owins[[nam]], coords = coords, namePPP = nam, ...)
        })
    names(owins) <- Nam
    if (addCellMarkers) {
        if (ct <- is.data.frame(cellTypes)) {
            if (!("cell" %in% names(cellTypes))) {
                stop("Variable names 'cell' must be contained in cell types dataframe")
            }
            otherCellNames <-
                names(cellTypes)[names(cellTypes) != "cell"]
            # Names of the other covariates
        }
        if (verbose) {
            message("Adding cell names for point pattern")
        }
        for (nn in rownames(hypFrame)) {
            if (verbose) {
                message(match(nn, rownames(hypFrame)), " of ", nrow(hypFrame))
            }
            ppp <- hypFrame[[nn, "ppp"]]
            NP <- npoints(ppp)
            if (any("cell" == names(marks(ppp, drop = FALSE)))) {
                stop("Cell markers already present in point pattern ",
                     nn)
            }
            if (findOverlappingOwins) {
                # Find overlap between cells
                foo <-
                    findOverlap(owins[[nn]], hypFrame[nn, "centroids"])
            }
            # Assign events to cells
            cellOut <- rep("NA", NP)
            idLeft <- seq_len(NP)
            for (i in names(owins[[nn]])) {
                idIn <- which(inside.owin(ppp[idLeft,], w = owins[[nn]][[i]]))
                cellOut[idLeft[idIn]] <- i
                # Don't overwrite, stick to first match, and do not detect
                # overlap
                idLeft <- idLeft[-idIn]
            }
            if (NP == (nOut <- length(idLeft))) {
                stop(
                    "All points lie outside all windows for point pattern ",
                    nn,
                    " check your input!"
                )
            }
            if ((nOut > 0) && warnOut) {
                warning(
                    nOut,
                    " points lie outside all windows for point pattern ",
                    nn,
                    " and were not assigned to a cell.\n",
                    immediate. = TRUE
                )
            }
            newmarks <-
                cbind(marks(ppp, drop = FALSE), cell = cellOut)
            if (ct) {
                newmarks <-
                    cbind(newmarks, cellTypes[match(cellOut, cellTypes$cell),
                                              otherCellNames,
                                              drop = FALSE])
            }
            marks(ppp) <- newmarks
            hypFrame[nn, "ppp"] <- ppp
        }
    }
    hypFrame$owins <- owins[rownames(hypFrame)]
    if (is.null(hypFrame$centroids)) {
        hypFrame$centroids <- lapply(hypFrame$owins, function(y) {
            t(vapply(y, FUN.VALUE = double(2), function(x) {
                tmp <- centroid.owin(x, as.ppp = FALSE)
                c("x" = tmp$x, "y" = tmp$y)
            }))
            # Store as simple matrix for memory reasons
        })
    }
    return(hypFrame)
}
