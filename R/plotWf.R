#' Plot the weighting function: The observation weigth as a function of the number of observations
#'
#' @param wf The fitted scam object
#' @param se A boolean, should standard errors be plotted?
#' @param ... Passed onto the plot.scam function for 1D splines
#'
#' @return
#' @import ggplot2
#' @export
#'
#' @examples
plotWf = function(wf, ...){
    if(grepl("Pair", attr(wf, "pi"))){
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
        plot(weight ~ NP, data = df, type = "l", xlab = "Number of observations", ylab = "Weight")
    }
}
