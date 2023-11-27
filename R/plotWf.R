#' Plot the weighting function: The observation weigth as a function of the number of observations
#'
#' @param wf The fitted scam object
#' @param pi
#' @param ... Passed onto the plot.scam function for 1D splines
#'
#' @return For univariate PI, returns a line plot; for bivariate PI a ggplot object
#' @import ggplot2
#' @export
#'
#' @examples
#' data(Yang)
#' hypYang = buildHyperFrame(Yang, coordVars = c("x", "y"),
#' designVar = c("day", "root", "section"))
#' yangPims = estPims(hypYang, pis = c("nn", "nnPair"), features = attr(hypYang, "features")[1:12])
#' #First Build the weight function
#' yangObj <- addWeightFunction(yangPims, designVars = c("day", "root"))
#' plotWf(yangObj, "nn")
#' plotWf(yangObj, "nnPair")
plotWf = function(resList, pi = c("nn", "nnPair", "allDist", "allDistPair"), ...){
    pi = match.arg(pi)
    if(is.null(wf <- resList$Wfs[[pi]])){
        stop("This type of pi (", pi, ") is not available in the object)")
    }
    if(grepl("Pair", pi)){
        tmp = exp(wf$model[, c("log(maxP)", "log(minP)")])
        colnames(tmp) = c("maxP", "minP")
        df = cbind("variance" = exp(predict.scam(wf, newdata = tmp)), tmp)
        ggplot(df, aes(x = log10(minP), y = log10(maxP), col = 1/variance)) +
            geom_point(size = 1) +
            scale_colour_gradient(low = "yellow", high = "blue", name = "Weight") +
            xlab("Log10 number of events for least expressed gene") +
            ylab("Log10 number of events for most expressed gene")
    } else {
        tmp = data.frame("NP" = exp(wf$model[, "log(NP)"]))
        df = cbind("weight" = 1/exp(predict.scam(wf, newdata = tmp)), tmp)
        plot(weight ~ NP, data = df[order(df$NP),], type = "l",
             xlab = "Number of observations", ylab = "Weight",
             main = paste("Weighing function for probabilistic indices of type", attr(wf, "pi")))
    }
}
