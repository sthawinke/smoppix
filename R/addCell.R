#' Check in which cell each event lies and add a cell marker,
#' plus include the list of regions
#'
#' @param hypFrame The hyperframe
#' @param owins the list of list of owins per hyperframe
#' @param checkOverlap a boolean, should windows be checked for overlap
#' @param warnOut a boolean, should warning be issued when points are not
#' contained in window
#' @return the modified hyperframe
#' @importFrom Rfast rowAny
#' @importFrom spatstat.geom is.owin
#' @export
#' @examples
#' library(spatstat.random)
#' n <- 1e3 # number of molecules
#' ng <- 25 # number of genes
#' nfov = 3 # Number of fields of view
#' conditions = 3
#' # sample xy-coordinates in [0, 1]
#' x <- runif(n)
#' y <- runif(n)
#' # assign each molecule to some gene-cell pair
#' gs <- paste0("gene", seq(ng))
#' gene <- sample(gs, n, TRUE)
#' fov <- sample(nfov, n, TRUE)
#' condition = sample(conditions, n, TRUE)
#' # construct data.frame of molecule coordinates
#' df <- data.frame(gene, x, y, fov, "condition" = condition)
#' #A list of point patterns
#' listPPP = tapply(seq(nrow(df)), df$fov, function(i){
#'     ppp(x = df$x[i], y = df$y[i], marks = df[i, "gene", drop = FALSE])
#' }, simplify = FALSE)
#' # Regions of interest (roi): Diamond in the center plus four triangles
#' w1 <- owin(poly=list(x=c(0,.5,1,.5),y=c(.5,0,.5,1)))
#' w2 <- owin(poly=list(x=c(0,0,.5),y = c(.5,0,0)))
#' w3 <- owin(poly=list(x=c(0,0,.5),y=c(1,0.5,1)))
#' w4 <- owin(poly=list(x=c(1,1,.5),y=c(0.5,1,1)))
#' w5 <- owin(poly=list(x=c(1,1,.5),y=c(0,0.5,0)))
#' hypFrame <- buildHyperFrame(df, coordVars = c("x", "y"),
#' imageVars = c("condition", "fov"))
#' nDesignFactors = length(unique(hypFrame$image))
#' wList = lapply(seq_len(nDesignFactors), function(x){
#'     list("w1" = w1, "w2" = w2, "w3" = w3, "w4" = w4, "w5" = w5)
#' })
#' names(wList) = rownames(hypFrame) #Matching names is necessary
#' hypFrame2 = addCell(hypFrame, wList)
addCell = function(hypFrame, owins, checkOverlap = TRUE, warnOut = TRUE){
    stopifnot(nrow(hypFrame) == length(owins),
    all(unlist(lapply(owins,
                      function(x) {vapply(x, FUN.VALUE = FALSE, is.owin)}))),
              all(rownames(hypFrame) %in% names(owins)),
    is(hypFrame, "hyperframe"))
    Feat = attr(hypFrame, "features")
    for(nn in rownames(hypFrame)){
        ppp = hypFrame[[nn, "ppp"]]
        if(any("cell" == names(marks(ppp, drop = FALSE)))){
            stop("Cell markers already present in point pattern ", nn)
        }
        if(checkOverlap){
            foo = findOverlap(owins[[nn]])
        }
        idWindow = vapply(owins[[nn]], FUN.VALUE = logical(NP <- npoints(ppp)),
                          function(rr){
                inside.owin(ppp, w = rr)
        })
        idIn = rowAny(idWindow)
        if(all(!idIn)){
            stop("All points lie outside all windows for point pattern ", nn,
                 " check your input!")
        }
        if(warnOut && (LL <- sum(!idIn))){
            warning(LL," points lie outside all windows for point pattern ", nn,
                    " and were not assigned to a cell.\n")
        }
        cellOut = rep("NA", NP)
        cellOut[idIn] = rownames(which(t(idWindow[idIn, ]), arr.ind = TRUE))
        hypFrame[[nn, "ppp"]] = setmarks(hypFrame[[nn, "ppp"]],
            cbind(marks(hypFrame[[nn, "ppp"]], drop = FALSE), cell = cellOut))
    }
    hypFrame$owins = owins
    attr(hypFrame, "features") = Feat #Retain features attribute,
    #which gets lost when setting the marks somehow
    return(hypFrame)
}

