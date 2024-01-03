#' Plot the weighting function: The observation weigth as a function of the number of observations
#'
#' @inheritParams buildDataFrame
#' @param pi The PI for which to plot the weighting function
#'
#' @return For univariate PI, returns a line plot; for bivariate PI a ggplot object
#' @importFrom ggplot2 ggplot geom_point scale_colour_gradient xlab ylab ggtitle aes
#' @export
#'
#' @examples
#' data(Yang)
#' hypYang <- buildHyperFrame(Yang,
#'     coordVars = c('x', 'y'),
#'     imageVars = c('day', 'root', 'section')
#' )
#' yangPims <- estPims(hypYang[c(seq_len(5), seq(25, 29)),],
#'     pis = c('nn', 'nnPair')
#' )
#' # First Build the weighting function
#' yangObj <- addWeightFunction(yangPims, designVars = c('day', 'root'))
#' plotWf(yangObj, 'nn')
#' plotWf(yangObj, 'nnPair')
plotWf <- function(obj, pi = obj$pis[1]) {
    pi <- match.arg(pi, choices = c("nn", "nnPair", "nnCell", "nnPairCell"))
    if (is.null(obj$Wfs)) {
        stop("No weight function added yet, run addWeightFunction first!")
    }
    if (!pi %in% obj$pis) {
        stop("This type of pi (", pi, ") is not available in the object)")
    }
    wf <- obj$Wfs[[pi]]
    if (grepl("Pair", pi)) {
        tmp <- wf$model[, c("maxP", "minP")]
        colnames(tmp) <- c("maxP", "minP")
        tmp <- tmp[tmp[, "maxP"] > 0 & tmp[, "minP"] > 0, ]
        df <- cbind(Weight = evalWeightFunction(wf, newdata = tmp), exp(tmp))
        df[, "Weight"] <- df[, "Weight"]/max(df[, "Weight"])
        ggplot(df, aes(x = log10(minP), y = log10(maxP), col = Weight)) +
            geom_point(size = 1) +
            scale_colour_gradient(low = "yellow",
            high = "blue", name = "Weight") +
            xlab("Log10 number of events for least expressed gene") +
            ylab("Log10 number of events for most expressed gene") +
        ggtitle(paste("Weighting function for probabilistic indices of type", pi))
    } else {
        tmp <- data.frame(NP = wf$model[, "log(NP)"])
        df <- cbind(Weight = evalWeightFunction(wf, newdata = tmp), exp(tmp))
        df[, "Weight"] <- df[, "Weight"]/max(df[, "Weight"])
        plot(Weight ~ NP, data = df[order(df$NP), ], type = "l",
             xlab = "Number of observations", ylab = "Weight",
        main = paste("Weighting function for probabilistic indices of type",
            pi))
    }
}
