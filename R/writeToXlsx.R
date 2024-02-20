#' Write linear modelling results to excel worksheet
#'
#' The results are written to different tabs for every sign (PI smaller than or larger than 0.5) of every PI,
#' sorted by increasing p-value.
#'
#' @param obj The results of linear model fitting
#' @param file The file to write the results to
#' @param overwrite A boolean, should the file be overwritten if it exists already?
#' @param digits An integer, the number of significant digits to retain for the PI,
#' raw and adjusted p-values
#' @param sigLevel The significance level threshold to use for the adjusted p-values,
#' only features exceeding the threshold are written to the file. Set this parameter to 1 to write all features
#'
#' @details If no feature exceeds the significance threshold for a certain pi and parameter combination,
#' an empty tab is created. For fixed effects, a single tab is written for PI differences of any sign
#' @return Returns invisible with a message when writing operation successful,
#' otherwise throws an error.
#' @export
#'
#' @examples
#' example(fitLMMs, 'spatrans')
#' writeToXlsx(lmmModels, "tmpFile.xlsx")
#' file.remove("tmpFile.xlsx")
#' @importFrom openxlsx createWorkbook writeData addWorksheet saveWorkbook getSheetNames
writeToXlsx = function(obj, file, overwrite = FALSE, digits = 3, sigLevel = 0.05){
    stopifnot(is.logical(overwrite), is.character(file), is.numeric(digits), is.numeric(sigLevel))
    if(!grepl("\\.xlsx", file)){
        message("Adding .xlsx extension to file")
        file = paste0(file, ".xlsx")
    }
    if(file.exists(file)){
        if(overwrite){
            message("Overwriting existing file")
        } else {
            stop("File ", file, " already exists! Set overwrite = TRUE to overwrite")
        }
    }
    pis = names(obj)
    fixedEffects = names(obj[[1]]$results$fixedEffects)
    wb = createWorkbook()
    for(pi in pis){
        for(effect in c("Intercept", fixedEffects)){
            mat = switch(effect, "Intercept" = obj[[pi]]$results$Intercept,
                         obj[[pi]]$results$fixedEffects[[effect]])
            mat = mat[!is.na(mat[, "pAdj"]),,drop = FALSE]
            #Only significant features
            subMat = mat[mat[, "pAdj"] < sigLevel,,drop = FALSE]
            #Rounding
            for(i in c("pVal", "pAdj")){
                subMat[, i] = signif(subMat[, i], digits)
            }
            for(i in setdiff(colnames(subMat), c("pVal", "pAdj"))){
                subMat[, i] = round(subMat[, i], digits)
            }
            for(smallPI in switch(effect, "Intercept" = c(TRUE, FALSE), TRUE)){
                sheetName = makeSheetName(pi, effect, smallPI)
                subMat2 = if(effect == "Intercept") {
                    subMat[match.fun(if (smallPI) "<" else ">")(subMat[, "Estimate"], 0.5),,drop = FALSE]
                    } else subMat
                addWorksheet(wb, sheetName) #Create sheet and write data to it
                writeData(wb, sheet = sheetName, x = data.frame(subMat2), colNames = TRUE, rowNames = TRUE)
            }
        }
    }
    saveWorkbook(wb, file = file, overwrite = overwrite)
    message(length(getSheetNames(file)), " tabs successfully written to ", file)
}
makeSheetName = function(pi, effect, smallPI = TRUE){
    keyword = if(grepl("Pair", pi)){
        if(smallPI) "Colocalization" else "Antilocalization"
    } else if(grepl("nn", pi)){
        if(smallPI) "Aggregation" else "Regularity"
    } else {
        paste(if(smallPI) "Close to" else "Far from", pi)
    }
    paste(keyword, "|", switch(effect, "Intercept" = "overall", effect))
}

