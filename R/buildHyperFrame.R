#' Built the hyperframe containing all separate point patterns that will
#'  be used throughout the analysis. Both matrices, dataframe and SpatialExperiment inputs are accepted.
#' @rdname buildHyperFrame
#' @importFrom spatstat.geom ppp hyperframe
#' @import methods
#' @export
#' @param x the input object, see methods("buildHyperFrame")
#' @param ... additional constructor arguments
#' @param covariates A vector of covariates
#' @examples
#' data(Yang)
#' hypYang = suppressWarnings(buildHyperFrame(Yang, coordVars = c("x", "y"),
#' designVar = c("day", "root", "section")))
setGeneric("buildHyperFrame", function(x, ...) standardGeneric("buildHyperFrame"))
#'
#' @rdname buildHyperFrame
#' @export
#' @param coordVars Names of coordinates
#' @param designVar A single character vector, name of the design variable
#' @param coVars Names of covariates such as gene or cell for each single point

setMethod("buildHyperFrame", "data.frame",
        function(x, coordVars, designVar, coVars = setdiff(names(x), c(designVar, coordVars)),...) {
        buildHyperFrame(as.matrix(x[, coordVars]), design = x[,designVar], covariates = x[, coVars, drop = FALSE],...)
})
#' @param matrix The input matrix
#' @param design A single design variable distinguishing the different point patterns
#'
#' @rdname buildHyperFrame
#' @export
setMethod("buildHyperFrame", "matrix", function(x, design, covariates, ...) {
    if(ncol(x)!=2){
        stop("Coordinates must be 2 dimensional")
    }
    if(nrow(x) != NROW(design)){
        stop("Number of rows of design and coordinate matrix must match")
    }
    designVec = if(NCOL(design)!=1){
        apply(design, 1, paste, collapse = "_")
    } else design

    if(!any("gene" == colnames(covariates))){
        stop("Gene identity must be supplied in covariate matrix")
    } else {
        covariates$gene = as.character(covariates$gene)
    }
    stopifnot(is.null(covariates) || nrow(x) == NROW(covariates))
    message("Found ", length(unDesignFactors <- unique(designVec)), " unique design factors")
    ppps = tapply(seq_len(nrow(x)), designVec, simplify = FALSE, function(i){
        spatstat.geom::ppp(x = x[i, 1], y = x[i, 2], marks = covariates[i,,drop = FALSE],
                           xrange = range(x[i, 1]), yrange = range(x[i, 2]), drop = FALSE)
    })
    #Replace underscore in gene names => Used to build pairs
    hypFrame = spatstat.geom::hyperframe("ppp" = ppps, design = names(ppps))
    desMat = t(simplify2array(strsplit(names(ppps), "_"))); colnames(desMat) = names(design)
    for(i in names(design)){hypFrame[, i] = desMat[, i]}
    hypFrame$tabObs = lapply(hypFrame$ppp, function(x) table(marks(x, drop = FALSE)$gene))
    attr(hypFrame, "features") = unique(unlist(lapply(hypFrame$tabObs, names)))
    return(hypFrame)
})
#' @param list A list of point patterns of class "ppp"
#'
#' @rdname buildHyperFrame
#' @export
setMethod("buildHyperFrame", "list", function(x,...) {
    hypFrame = spatstat.geom::hyperframe("ppp" = x, design = names(x))
    return(hypFrame)
})
#' @rdname buildHyperFrame
#' @export
#' @importFrom SpatialExperiment SpatialExperiment spatialCoords
#' @importFrom SummarizedExperiment colData
setMethod("buildHyperFrame", "SpatialExperiment", function(x, designVar, coVars, ...) {
    buildHyperFrame(spatialCoords(x), design = colData(x)[, designVar],
                    covariates = cbind("gene" = rownames(colData(x)), as.data.frame(colData(x))[, coVars, drop = FALSE]))
})

