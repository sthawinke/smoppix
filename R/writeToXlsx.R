#' Write linear modelling results to excel worksheet
#'
#' The results are written to different tabs for every sign (PI smaller than or larger than 0.5) of every PI,
#' sorted by increasing p-value.
#'
#' @param obj The results of linear model fitting
#' @param file The file to write the results to
#'
#' @return Returns invisible with a message when writing operation successful,
#' otherwise throws an error.
#' @export
#'
#' @examples
#' example(fitLMMs, 'spatrans')
#' writeToXslx(yangPims, "tmpFile.xlsx")
#' file.remove("tmpFile.xlsx")
#' @importFrom openxlsx createWorkbook writeData addWorksheet saveWorkbook
writeToXlsx = function(obj, file, overwrite = FALSE){
    if(!grepl("\\.xlsx", file)){
        message("Adding .xlsx extension to file")
        file = paste0(file, ".xlsx")
    }
    if(file.exists(file)){
        if(overwrite){
            message("Overwriting existing file")
        } else {
            stop("File ", file, "already exists! Set overwrite = TRUE to overwrite")
        }
    }
    tmp = PIsignifNN$Intercept
    tmp[, "pi"] = round(tmp[, "pi"], 3)
    tmp[, "pval"] = signif(tmp[, "pval"], 3)
    tmp[, "qval"] = signif(tmp[, "qval"], 3)
    geneMat = t(sapply(rownames(tmp), sund))
    colnames(geneMat) = paste("Gene", 1:2)
    tmp = data.frame(tmp, geneMat)
    idQval = tmp[, "qval"] < sigLevel
    colocDf = tmp[idQval & tmp[, "pi"] < 0.5, ]
    colocDf = colocDf[order(colocDf[, "qval"]), ]
    wb = createWorkbook()
    for(i in c("Colocalization", "Antilocalization")) {addWorksheet(wb, i)}
    writeData(wb, sheet = "Colocalization", x = data.frame(colocDf), colNames = TRUE, rowNames = TRUE)
    antilocDf = tmp[idQval & tmp[, "pi"] > 0.5, ]
    antilocDf = antilocDf[order(antilocDf[, "qval"]), ]
    writeData(wb, sheet = "Antilocalization", x = data.frame(antilocDf), colNames = TRUE, rowNames = TRUE)
    saveWorkbook(wb, file = xlsxName, overwrite = TRUE)
}
