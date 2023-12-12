#' Built the hyperframe containing all separate point patterns that will
#'  be used throughout the analysis. Both matrices, dataframe and SpatialExperiment inputs are accepted.
#' @rdname buildHyperFrame
#' @importFrom spatstat.geom ppp hyperframe
#' @import methods
#' @export
#' @param x the input object, see methods("buildHyperFrame")
#' @param ... additional constructor arguments
#' @param covariates A vector of covariates
#' @return An object of class 'hyperframe' from the 'spatstat' package
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
#' @param imageVars A character vector of variables whose unique combinations define the separate point patterns (images)
#' @param coVars Names of covariates such as gene or cell for each single point

setMethod(
    "buildHyperFrame", "data.frame",
    function(x, coordVars, imageVars, coVars = setdiff(names(x), c(imageVars, coordVars)), ...) {
        buildHyperFrame(as.matrix(x[, coordVars]),
            image = x[, imageVars, drop = FALSE],
            covariates = x[, coVars, drop = FALSE], ...
        )
    }
)
#' @param matrix The input matrix
#' @param image A set of design variable whose unique combination distinguishes
#' the different point patterns
#'
#' @rdname buildHyperFrame
#' @export
setMethod("buildHyperFrame", "matrix", function(x, image, covariates, ...) {
    if (length(intersect(names(image), names(covariates)))) {
        stop("Overlap detected between 'imageVars' and 'coVars'.
        The combination of the first separates the images,
             wheras the second defines point-specific data such as
             gene identifiers.")
    }
    if (ncol(x) != 2) {
        stop("Coordinates must be 2 dimensional")
    }
    if (nrow(x) != NROW(image)) {
        stop("Number of rows of image and coordinate matrix must match")
    }
    designVec <- if (NCOL(image) != 1) {
        apply(image, 1, paste, collapse = "_")
    } else {
        image
    }

    if (!any("gene" == colnames(covariates))) {
        stop("Gene identity must be supplied in covariate matrix")
    } else {
        covariates$gene <- as.character(covariates$gene)
    }
    stopifnot(is.null(covariates) || nrow(x) == NROW(covariates))
    message("Found ", length(unDesignFactors <- unique(designVec)), " unique images")
    ppps <- tapply(seq_len(nrow(x)), designVec, simplify = FALSE, function(i) {
        spatstat.geom::ppp(
            x = x[i, 1], y = x[i, 2], marks = covariates[i, , drop = FALSE],
            xrange = range(x[i, 1]), yrange = range(x[i, 2]), drop = FALSE
        )
    })
    # Replace underscore in gene names => Used to build pairs
    hypFrame <- spatstat.geom::hyperframe("ppp" = ppps, image = names(ppps))
    rownames(hypFrame) <- hypFrame$image
    desMat <- matrix(simplify2array(strsplit(names(ppps), "_")),
        nrow = nrow(hypFrame), byrow = TRUE
    )
    colnames(desMat) <- names(image)
    for (i in names(image)) {
        hypFrame[, i] <- desMat[, i]
    }
    hypFrame <- addTabObs(hypFrame)
    attr(hypFrame, "features") <- unique(unlist(lapply(hypFrame$tabObs, names)))
    attr(hypFrame, "imageVars") <- colnames(image)
    return(hypFrame)
})
#' @param list A list of matrices or of point patterns of class "ppp"
#'
#' @rdname buildHyperFrame
#' @export
setMethod("buildHyperFrame", "list", function(x, coordVars = c("x", "y"),...) {
    if(all(vapply(x, is.ppp, FUN.VALUE = TRUE))){
        hypFrame <- spatstat.geom::hyperframe(
            "ppp" = x, "image" = names(x),
            "tabObs" = lapply(x, function(y) table(marks(y, drop = FALSE)$gene))
        )
    } else if(all(vapply(x, is.matrix, FUN.VALUE = TRUE))){
        hypFrame <- spatstat.geom::hyperframe(
            "ppp" = lapply(x, function(z) {
                covariates = setDiff(colnames(z), coordVars)
                if(!("gene" %in% covariates)){
                    stop("Gene marker is missing in at least one point pattern")
                }
                spatstat.geom::ppp(
                x = z[, coordVars[1]], y = z[, coordVars[2]],
                marks = z[, covariates, drop = FALSE],
                xrange = range(z[, coordVars[1]]),
                yrange = range(z[, coordVars[2]]), drop = FALSE
            )}),
            "image" = names(x)
        )
        hypFrame <- addTabObs(hypFrame)
        buildHyperFrame(Reduce(f = rbind, x), ...)
    } else{
        stop("Supply a list of point patterns (ppp) or of matrices")
    }
    return(hypFrame)
})
#' @rdname buildHyperFrame
#' @export
#' @importFrom SpatialExperiment SpatialExperiment spatialCoords
#' @importFrom SummarizedExperiment colData
setMethod("buildHyperFrame", "SpatialExperiment", function(x, imageVars, coVars, ...) {
    buildHyperFrame(spatialCoords(x),
        image = as.data.frame(colData(x)[, imageVars, drop = FALSE]),
        covariates = cbind("gene" = rownames(colData(x)), as.data.frame(colData(x))[, coVars, drop = FALSE])
    )
})
