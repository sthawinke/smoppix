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
        subDf <- dfMoran[dfMoran$image == im & !is.na(dfMoran$pi), ]
        moransI(subDf$pi, W = weightMats[[im]][subDf$cell, subDf$cell, drop = FALSE])
    })
    cbind(MoransI = moransIs["MoransI",], "Variance" = moransIs["VarMoransI",],
          piMat[match(unIm, piMat$image), setdiff(colnames(piMat), "weight")])
}
