#' @rdname buildHyperFrame
#' @export
#' @param coordVars Names of coordinates
#' @param imageIdentifier A character vector of variables whose unique combinations
#' define the separate point patterns (images)
#' @param imageVars Covariates belonging to the point patterns
#' @param pointVars Names of event-wise covariates such as gene or cell for each single point
setMethod("buildHyperFrame", "data.frame", function(x, coordVars, imageIdentifier = imageVars, imageVars, pointVars = setdiff(
                                                        names(x),
                                                        c(imageVars, imageIdentifier, coordVars, featureName)
                                                    ), featureName = "gene", ...) {
    buildHyperFrame(as.matrix(x[, coordVars]), imageIdentifier = x[,imageIdentifier],
                    imageVars = x[,imageVars],
        covariates = x[, c(pointVars, featureName), drop = FALSE], featureName = featureName,  ...
    )
})
#' @rdname buildHyperFrame
#' @export
setMethod("buildHyperFrame", "matrix", function(x, imageVars, imageIdentifier = imageVars,
                                                covariates, featureName = "gene", ...) {
    if (length(intersect(names(imageVars), names(covariates)))) {
        stop("Overlap detected between 'imageVars' and 'covariates'.
        The combination of the first separates the images,
             wheras the second defines point-specific data such as
             gene identifiers.")
    }
    if (ncol(x) != 2) {
        stop("Coordinates must be 2 dimensional")
    }
    if (nrow(x) != NROW(imageVars)) {
        stop("Number of rows of imageVars and coordinate matrix must match")
    }
    designVec <- makeDesignVar(imageIdentifier)
    covariates <- as.data.frame(covariates)
    if (!any(featureName == colnames(covariates))) {
        stop("Gene or cell identity\n", featureName, "\nmust be supplied in covariate matrix")
    } else {
        covariates$gene <- as.character(covariates[[featureName]])
    }
    stopifnot(is.null(covariates) || nrow(x) == NROW(covariates))
    message("Found ", length(unDesignFactors <- unique(designVec)), " unique images")
    ppps <- tapply(seq_len(nrow(x)), designVec, simplify = FALSE, function(i) {
        i <- i[order(x[i, 1])]
        # Pre-order for nncross function
        spatstat.geom::ppp(
            x = x[i, 1], y = x[i, 2], marks = covariates[i, , drop = FALSE],
            xrange = range(x[i, 1]), yrange = range(x[i, 2]), drop = FALSE
        )
    })
    hypFrame <- spatstat.geom::hyperframe(ppp = ppps, image = names(ppps))
    hypFrame <- addTabObs(hypFrame)
    hypFrame <- addDesign(hypFrame, imageVars, designVec)
    return(hypFrame)
})
#' @param list A list of matrices or of point patterns of class 'spatstat.geom::ppp'
#' @param idVar An optional id variable present in covariates, that is matched
#' with the names of covariates
#'
#' @rdname buildHyperFrame
#' @export
#' @importFrom spatstat.geom is.ppp
setMethod("buildHyperFrame", "list", function(
    x, coordVars = c("x", "y"), covariates = NULL,
    idVar = NULL, featureName = "gene", ...) {
    stopifnot(is.null(covariates) || (nrow(covariates) == length(x)), is.null(idVar) ||
        idVar %in% names(covariates))
    if (!is.null(names(x)) && !is.null(idVar)) {
        covariates <- covariates[match(names(x), covariates[, idVar]), ]
        # Match the ranking
    } else if (!is.null(covariates)) {
        message("No matching information supplied, matching point patterns ", "and covariates based on order.")
    }

    if (all(vapply(x, is.ppp, FUN.VALUE = TRUE))) {
        hypFrame <- spatstat.geom::hyperframe(ppp = x, image = names(x))
    } else if (all(vapply(x, function(y) {
        is.matrix(y) || is.data.frame(y)
    }, FUN.VALUE = TRUE))) {
        hypFrame <- spatstat.geom::hyperframe(ppp = lapply(x, function(z) {
            PPcovariates <- setdiff(colnames(z), coordVars)
            if (!(featureName %in% PPcovariates)) {
                stop("Gene or cell marker\n", featureName, "\nis missing in at least one point pattern")
            }
            i <- order(z[, coordVars[1]])
            spatstat.geom::ppp(x = z[i, coordVars[1]], y = z[i, coordVars[2]], marks = z[i,
                PPcovariates,
                drop = FALSE
            ], xrange = range(z[, coordVars[1]]), yrange = range(z[
                ,
                coordVars[2]
            ]), drop = FALSE)
        }), image = names(x))
    } else {
        stop("Supply a list of point patterns (ppp)", "or of dataframes or matrices")
    }
    hypFrame <- addTabObs(hypFrame)
    hypFrame <- addDesign(hypFrame, covariates, hypFrame$image)
    return(hypFrame)
})
#' @rdname buildHyperFrame
#' @export
#' @importFrom SpatialExperiment SpatialExperiment spatialCoords
#' @importFrom SummarizedExperiment colData
setMethod("buildHyperFrame", "SpatialExperiment", function(
    x, imageVars, pointVars, imageIdentifier = imageVars,
    ...) {
    buildHyperFrame(spatialCoords(x),
                    imageVars = as.data.frame(colData(x)[, imageVars, drop = FALSE]),
                    imageIdentifier = as.data.frame(colData(x)[, imageIdentifier,drop = FALSE]),
                    covariates = cbind("gene" = rownames(colData(x)), as.data.frame(colData(x))[,pointVars,drop = FALSE]))
})
