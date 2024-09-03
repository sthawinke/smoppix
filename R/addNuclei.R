#' Add nuclei to a hyperframe already containing cells.
#' @description Add the list of the cells and their centroids in the hyperframe,
#' check in which cell each event lies and add a cell marker.
#'
#' @param hypFrame The hyperframe
#' @param nucleiList the list containing a list of owins per point pattern.
#' The length of the list must match the length of the hyperframe, and the names must match.
#' Also lists of geojson objects, coordinate matrices or rois are accepted, see \link{addCell}
#' @param coords The names of the coordinates, if the windows are given as sets of coordinates.
#' @param verbose A boolean, should verbose output be printed?
#' @param checkSubset A boolean, should be checked whether nuclei are encompassed by cells?
#' @param overwriteNuclei A boolean, should existing nuclei be replaced?
#' @param ... Further arguments passed onto \link{convertToOwins}
#'
#' @return The hyperframe with nuclei added as entry
#' @importFrom spatstat.geom is.subset.owin
#' @export
#' @seealso \link{addCell}, \link{convertToOwins}
#' @details The nuclei names must match the cell names already present, all other nuclei are dropped. 
#' A warning is issued when nuclei are not encompassed by their cell
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
#' # The nuclei
#' n1 <- owin(poly = list(x = c(0.2, .4, 0.8, .4), y = c(.4, .2, .4, .8)))
#' n2 <- owin(poly = list(x = c(0.1, 0.1, .4), y = c(.4, .1, .1)))
#' n3 <- owin(poly = list(x = c(0.1, 0.1, .4), y = c(1, .75, 1)))
#' n4 <- owin(poly = list(x = c(1, 1, .6), y = c(.7, .9, .9)))
#' n5 <- owin(poly = list(x = c(.95, .95, .7), y = c(.1, .4, .1)))
#' nList <- lapply(seq_len(nDesignFactors), function(x) {
#'     list("w1" = n1, "w2" = n2, "w3" = n3, "w4" = n4, "w5" = n5)
#' })
#' names(nList) <- rownames(hypFrame) # Matching names is necessary
#' hypFrame3 <- addNuclei(hypFrame2, nList)
addNuclei <- function(hypFrame, nucleiList, checkSubset = TRUE,
                    verbose = TRUE, coords = c("x", "y"), ...) {
    stopifnot(
        nrow(hypFrame) == length(nucleiList),
        all(rownames(hypFrame) %in% names(nucleiList)),
        is.hyperframe(hypFrame),
        is.list(nucleiList)
    )
    if(!is.null(hypFrame$nuclei) && !overwriteNuclei){
      stop("Cells already present in hyperframe!
             Set overwriteNuclei=TRUE to replcae them")
    }
    if (verbose) {
        message("Converting nuclei to spatstat owins")
    }
    # Convert different types of windows to owins
    nucleiList <- loadBalanceBplapply(Nam <- names(nucleiList), function(nam){
            convertToOwins(nucleiList[[nam]], coords = coords, namePPP = nam, ...)
        })
    names(nucleiList) <- Nam
    hypFrame$nuclei <- nucleiList[rownames(hypFrame)]
    if(checkSubset){
      if(is.null(hypFrame$owins)){
        stop("No cells present in hyperframe so I cannot check if nuclei are encompassed by them!")
      }
      noSubSet <- lapply(rownames(hypFrame), function(i){
        which(!vapply(names(hypFrame$nuclei[[i]]), FUN.VALUE = logical(1), function(x){
          is.subset.owin(hypFrame$nuclei[[i]][[x]], hypFrame$owins[[i]][[x]])
        }))
      })
      if(any(vapply(noSubSet, FUN.VALUE = integer(1), length) > 0)){
        warning("Not all nuclei are completely contained in their cells")
      }
    }
    return(hypFrame)
}
