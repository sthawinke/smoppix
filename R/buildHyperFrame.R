#' Built the hyperframe containing all separate point patterns that will
#'  be used throughout the analysis. Both matrices, dataframe and SpatialExperiment inputs are accepted.
#' @rdname buildHyperFrame
#' @importFrom spatstat.geom ppp hyperframe
#' @import methods
#' @export
#' @param x the input object, see methods("buildHyperFrame")
#' @param ... additional constructor arguments
#' @param covariates A vector of covariates
setGeneric("buildHyperFrame", function(x, ...) standardGeneric("buildHyperFrame"))
#'
#' @rdname buildHyperFrame
#' @export
#' @param coordVars Names of coordinates
#' @param designVar A single character vector, name of the design variable
#' @param coVars Names of covariates such as gene or cell for each single point

setMethod("buildHyperFrame", "data.frame",
          function(x, coordVars, designVar, coVars = setdiff(names(x), c(designVar, coordVars)), ...) {
    buildHyperFrame(as.matrix(x[, coordVars]), design = x[,designVar], covariates = x[, coVars],...)
})
#' @param matrix The input matrix
#' @param design A single design variable distinguishing the different point patterns
#'
#' @rdname buildHyperFrame
#' @export
setMethod("buildHyperFrame", "matrix", function(x, design, covariates,...) {
    if(ncol(x)!=2){
        stop("Coordinates must be 2 dimensional")
    }
    if(nrow(x) != NROW(design)){
        stop("Number of rows of design and coordinate matrix must match")
    }
    if(NCOL(design)!=1){
        stop("Design variable distinguishes different ppint patterns and should be one dimensionsal")
    }
    if(!any("gene" == colnames(covariates))){
        stop("Gene identity must be supplied in covariate matrix")
    }
    stopifnot(is.null(covariates) || nrow(x) == NROW(covariates))
    message("Found ", length(unDesignFactors <- unique(design)), " unique design factors")
    ppps = tapply(seq_len(nrow(x)), design, function(i){
        spatstat.geom::ppp(x = x[i, 1], y = x[i, 2], marks = covariates[i,])
    })
    hypFrame = spatstat.geom::hyperframe("ppp" = ppps, design = names(ppps))
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
setMethod("buildHyperFrame", "SpatialExperiment", function(x, ...) {

})

