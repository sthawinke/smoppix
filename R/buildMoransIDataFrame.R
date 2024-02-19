#' Build a data frame with Moran's I as outcome variable
#'
#' @inheritParams buildDataFrame
#'
#' @return A modified data frame
buildMoransIDataFrame <- function(pi, piMat, weightMats) {
    dfMoran <- if (pi %in% c("edge", "centroid")) {
        aggregate(pi ~ cell + image, piMat, FUN = mean, na.rm = TRUE)
    } else {
        piMat
    }
    moransIs <- vapply(unIm <- unique(dfMoran$image), FUN.VALUE = double(2), function(im) {
        subDf <- dfMoran[dfMoran$image == im, ]
        moransI(subDf$pi, W = weightMats[[im]][subDf$cell, subDf$cell])
    })
    W = 1/moransIs["VarMoransI",]
    cbind(MoransI = moransIs["MoransI",], "weight" = W/sum(W, na.rm = TRUE),
          piMat[match(unIm, piMat$image), setdiff(colnames(piMat), "weight")])
}
