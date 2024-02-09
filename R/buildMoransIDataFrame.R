#' Build a data frame with Moran's I as outcome variable
#'
#' @inheritParams buildDataFrame
#'
#' @return A modified data frame
buildMoransIDataFrame = function(pi, piMat, weightMats){
    dfMoran <- if (pi %in% c("edge", "centroid")) {
        aggregate(pi ~ cell + image, piMat, FUN = mean, na.rm = TRUE)
    } else {
        piMat
    }
    moransIs <- vapply(unIm <- unique(dfMoran$image), FUN.VALUE = double(1), function(im) {
        subDf <- dfMoran[dfMoran$image == im, ]
        moransI(subDf$pi, W = weightMats[[im]][subDf$cell, subDf$cell])
    })
    cbind(MoransI = moransIs, piMat[match(unIm, piMat$image), setdiff(colnames(piMat), "weight")])
}
