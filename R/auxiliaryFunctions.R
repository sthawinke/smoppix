#' Helper function to spit gene pairs
#'
#' @param x character string
#' @param sep The character used to split
#' @return The split string
#' @export
#' @examples
#' GenePair <- "gene1--gene2"
#' sund(GenePair)
sund <- function(x, sep = "--") {
    strsplit(x, sep)[[1]]
}
#' Make design variable by combining different design variables
#'
#' @param x the design matrix
#' @param designVars the design variables to be combined
#' @param sep The string to separate the components
#'
#' @return a vector of design levels
makeDesignVar <- function(x, designVars, sep = "_") {
    if(NCOL(x)==1){
        unlist(x)
    } else {
        apply(x[, designVars, drop = FALSE], 1, paste, collapse = sep)
    }
}
#' An aux function to build gene pairs
#'
#' @param genes The genes to be combined
#'
#' @return A character vector of gene pairs
#' @export
#' @examples
#' genes <- paste0("gene", seq_len(4))
#' makePairs(genes)
makePairs <- function(genes) {
    apply(combn(genes, 2), 2, paste, collapse = "--")
}
#' Subsample a point pattern when it is too large
#'
#' @param p The point pattern
#' @param nSims The maximum size
#' @param returnId A boolean, should the id of the sampled elements be returned?
#'
#' @return A point pattern, subsampled if necessary
subSampleP <- function(p, nSims, returnId = FALSE) {
    Pout <- if (tooBig <- (NP <- npoints(p)) > nSims) {
        p[id <- sample(NP, nSims), ]
    } else {
        p
    }
    if (!tooBig && returnId) {
        id <- seq_len(npoints(p))
    }
    if (returnId) {
        return(list(Pout = Pout, id = id))
    } else {
        return(Pout)
    }
}
#' A version of contr.sum that retains names, a bit controversial but also
#' clearer
#' @note After
#' https://stackoverflow.com/questions/24515892/r-how-to-contrast-code-factors-and-retain-meaningful-labels-in-output-summary
#' @param x,... passed on to contr.sum
#' @importFrom stats contr.sum
#' @return The matrix of contrasts
named.contr.sum <- function(x, ...) {
    lev <- x
    x <- contr.sum(x, ...)
    colnames(x) <- lev[-length(lev)]
    x
}
#' Add tables with gene counts to the hyperframe, presort by gene and x-ccordinate
#' and add design varibales
#'
#' @param hypFrame The hyperframe
#'
#' @return The hyperframe with tabObs added
#' @importFrom spatstat.geom marks
addTabObs <- function(hypFrame) {
    hypFrame$tabObs <- lapply(hypFrame$ppp, function(x) {
        table(marks(x, drop = FALSE)$gene)
    })
    hypFrame
}
#' Add design variables to hyperframe
#'
#' @param hypFrame The hyperframe
#' @param desMat The design matrix
#' @param designVec The design vector
#'
#' @return The hyperframe with design variables added
addDesign = function(hypFrame, desMat, designVec){
    # Add design variables
    id <- match(hypFrame$image, designVec)
    for (i in colnames(desMat)) {
        hypFrame[, i] <- desMat[id, i]
    }
    rownames(hypFrame) <- hypFrame$image
    hypFrame
}
#' Nest random effects within fixed variables, in case the names are the same
#'
#' @param df The dataframe
#' @param randomVars The random variables
#' @param fixedVars The fixed variables
#'
#' @return The dataframe with adapted randomVars
nestRandom <- function(df, randomVars, fixedVars) {
    for (i in randomVars) {
        df[, i] <- apply(df[, c(fixedVars, i), drop = FALSE], 1, paste, collapse = "_")
    }
    df
}
#' Extract coordinates from a point pattern or data frame
#' @param x the point pattern, dataframe or matrix
#' @return the matrix of coordinates
getCoordsMat <- function(x) {
    if (is.ppp(x)) {
        cbind(x = x$x, y = x$y)
    } else if (is.data.frame(x)) {
        as.matrix(x)
    } else if (is.matrix(x)) {
        x
    }
}
#' Extract en element from a matrix or vector
#'
#' @param x the matrix or vector
#' @param e The column or element name
#'
#' @return The desired element
getElement <- function(x, e) {
    if (is.matrix(x)) {
        x[, e]
    } else {
        x[e]
    }
}
#' Extract the hyperframe
#' @param x The hyperframe, or list containing one
#' @return the hyperframe
#' @importFrom spatstat.geom is.hyperframe
getHypFrame <- function(x) {
    if (is.hyperframe(x)) {
        x
    } else {
        x$hypFrame
    }
}
#' Split a number of plots into rows and colums
#'
#' @param x The number of plots
#'
#' @return A vector of length 2 with required number of rows and colums
splitWindow <- function(x) {
    Sqrt <- sqrt(x)
    a <- ceiling(Sqrt)
    b <- floor(Sqrt)
    if (a * b < x) {
        b <- b + 1
    }
    c(a, b)
}
#' A wrapper for C-functions calculating cross-distance matrix fast
#' @param x,y the matrices or point patterns between which to calculate the
#' cross distances
#' @return a matrix of cross distances
crossdistWrapper <- function(x, y) {
    crossdistFastLocal(getCoordsMat(x), getCoordsMat(y))
}
#' Center numeric variables
#'
#' @param x The dataframe whose numeric variables are being centered
#'
#' @return The adapeted dataframe
centerNumeric = function(x){
    numId <- vapply(x, FUN.VALUE = TRUE, is.numeric)
    for(i in which(numId)){
        x[,i] <- x[,i] - mean(x[,i])
    }
    x
}
