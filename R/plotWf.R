#' Title
#'
#' @param wf The fitted scam object
#' @param se A boolean, should standard errors be plotted?
#' @param ... Passed onto the plot.scam function for 1D splines
#'
#' @return
#' @export
#'
#' @examples
plotWf = function(wf, se = FALSE, ...){
    if(grepl(attr(wf, "Pair"))){
        df = cbind("variance" = c(predict.scam(wf)), exp(wf$model[, c("log(maxP)", "log(minP)")]))
        colnames(df) = c("maxP", "minP")
        ggplot(df, aes(x = log10(minP), y = log10(maxP), col = variance)) +
            geom_point(size = 1) +
            scale_colour_gradient(low = "yellow", high = "blue", name = "Weight") +
            xlab("Log10 number of events for least expressed gene") +
            ylab("Log10 number of events for most expressed gene")
    } else {
        plot(wf, xlab = "Log number of observations", ylab = "Log-variance", se = se)
    }
}
