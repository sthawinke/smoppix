#' Build a hyperframe containing all point patterns of an experiment.
#' @description
#' Build a spatstat hyperframe with point patterns and metadata.
#' Matrices, dataframe, lists and SpatialExperiment inputs are accepted.
#' @rdname buildHyperFrame
#' @importFrom spatstat.geom ppp hyperframe
#' @importFrom methods setGeneric setMethod
#' @export
#' @param x the input object, see methods('buildHyperFrame')
#' @param ... additional constructor arguments
#' @param covariates A vector of covariates
#' @return An object of class 'hyperframe' from the 'spatstat.geom' package
#' @examples
#' data(Yang)
#' hypYang <- buildHyperFrame(Yang,
#'     coordVars = c("x", "y"),
#'     imageVars = c("day", "root", "section")
#' )
setGeneric("buildHyperFrame", function(x, ...) standardGeneric("buildHyperFrame"))
#'
#' @rdname buildHyperFrame
#' @export
#' @param coordVars Names of coordinates
#' @param imageVars A character vector of variables whose unique combinations
#' define the separate point patterns (images)
#' @param coVars Names of event-wise covariates such as gene or cell for each single point

setMethod("buildHyperFrame", "data.frame", function(x, coordVars, imageVars, coVars = setdiff(names(x), c(
                                                        imageVars,
                                                        coordVars
                                                    )), ...) {
    buildHyperFrame(as.matrix(x[, coordVars]),
        image = x[, imageVars, drop = FALSE], covariates = x[, coVars, drop = FALSE],
        ...
    )
})
#' @param matrix The input matrix
#'
#' @rdname buildHyperFrame
#' @export
setMethod("buildHyperFrame", "matrix", function(x, imageVars, covariates, ...) {
    if (length(intersect(names(imageVars), names(covariates)))) {
        stop("Overlap detected between 'imageVars' and 'coVars'.
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
    designVec <- if (NCOL(imageVars) != 1) {
        apply(imageVars, 1, paste, collapse = "_")
    } else {
        c(imageVars)
    }
    if (!any("gene" == colnames(covariates))) {
        stop("Gene identity must be supplied in covariate matrix")
    } else {
        covariates$gene <- as.character(covariates$gene)
    }
    stopifnot(is.null(covariates) || nrow(x) == NROW(covariates))
    message("Found ", length(unDesignFactors <- unique(designVec)), " unique images")
    ppps <- tapply(seq_len(nrow(x)), designVec, simplify = FALSE, function(i) {
        i = i[order(x[i, 1])]
        #Pre-order for nncross function
        spatstat.geom::ppp(
            x = x[i, 1], y = x[i, 2], marks = covariates[i, , drop = FALSE],
            xrange = range(x[i, 1]),
            yrange = range(x[i, 2]), drop = FALSE
        )
    })
    # Replace underscore in gene names => Used to build pairs
    hypFrame <- spatstat.geom::hyperframe(ppp = ppps, image = names(ppps))
    desMat <- matrix(simplify2array(strsplit(names(ppps), "_")),
        nrow = nrow(hypFrame), byrow = TRUE
    )
    colnames(desMat) <- names(imageVars)
    hypFrame <- addTabObs(hypFrame, desMat)
    return(hypFrame)
})
#' @param list A list of matrices or of point patterns of class 'ppp'
#' @param covariatesDf A data frame or matrix of covariates
#' @param idVar An optional id variable present in covariatesDf, that is matched
#' with the names of covariatesDf
#'
#' @rdname buildHyperFrame
#' @export
#' @importFrom spatstat.geom is.ppp
setMethod("buildHyperFrame", "list", function(x, coordVars = c("x", "y"),
                                              covariatesDf = NULL, idVar = NULL, ...) {
    stopifnot(
        is.null(covariatesDf) || (nrow(covariatesDf) == length(x)),
        is.null(idVar) || idVar %in% names(covariatesDf)
    )
    if (!is.null(names(x)) && !is.null(idVar)) {
        covariatesDf <- covariatesDf[match(names(x), covariatesDf[, idVar]), ]
        # Match the ranking
    } else if (!is.null(covariatesDf)) {
        message(
            "No matching information supplied, matching point patterns ",
            "and covariates based on order."
        )
    }

    if (all(vapply(x, is.ppp, FUN.VALUE = TRUE))) {
        hypFrame <- spatstat.geom::hyperframe(ppp = x, image = names(x))
    } else if (all(vapply(x, function(y) {
        is.matrix(y) || is.data.frame(y)
    }, FUN.VALUE = TRUE))) {
        hypFrame <- spatstat.geom::hyperframe(ppp = lapply(x, function(z) {
            PPcovariates <- setdiff(colnames(z), coordVars)
            if (!("gene" %in% PPcovariates)) {
                stop("Gene marker is missing in at least one point pattern")
            }
            i = order(z[, coordVars[1]])
            spatstat.geom::ppp(
                x = z[i, coordVars[1]], y = z[i, coordVars[2]], marks = z[i, PPcovariates, drop = FALSE],
                xrange = range(z[, coordVars[1]]), yrange = range(z[, coordVars[2]]), drop = FALSE
            )
        }), image = names(x))
    } else {
        stop("Supply a list of point patterns (ppp)", "or of dataframes or matrices")
    }
    hypFrame <- addTabObs(hypFrame, covariatesDf)
    return(hypFrame)
})
#' @rdname buildHyperFrame
#' @export
#' @importFrom SpatialExperiment SpatialExperiment spatialCoords
#' @importFrom SummarizedExperiment colData
setMethod("buildHyperFrame", "SpatialExperiment", function(x, imageVars, coVars, ...) {
    buildHyperFrame(spatialCoords(x), image = as.data.frame(colData(x)[, imageVars, drop = FALSE]), covariates = cbind(
        gene = rownames(colData(x)),
        as.data.frame(colData(x))[, coVars, drop = FALSE]
    ))
})
