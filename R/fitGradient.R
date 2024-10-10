#' Test for presence of gradient in a hyperframe of point patterns
#'
#' A Poisson process is fitted to the data assuming exponential relationship wit intensity of
#' the interaction between x and y variables and image identifier. This is
#' compared to a model without this interaction to test for the significance of the gradient.
#'
#' @param hypFrame the hyperframe
#' @param returnModel A boolean, should the entire model be returned?
#' Otherwise the p-value and coefficient vector are returned
#' @param ... passed onto \link[spatstat.model]{mppm}
#' @inheritParams estGradientsSingle
#' @inheritParams estGradients
#'
#' @return A list contraining
#' \item{pVal}{The p-value for existence of gradients}
#' \item{coef}{The model coefficients}
#' or a mppm model when returnModel is true
#' @importFrom stats formula coef
#' @importFrom spatstat.model mppm anova.mppm
#' @seealso \link{estGradients}
fitGradient <- function(hypFrame, fixedForm, randomForm, fixedFormSimple,
                        returnModel = FALSE, silent, ...) {
    xyModel <- try(mppm(data = hypFrame, fixedForm, random = randomForm, ...),
        silent = silent
    )
    out <- if (returnModel) {
        xyModel
    } else if (is(xyModel, "try-error")) {
        list("pVal" = 1, "coef" = NULL)
    } else {
        intModel <- try(mppm(data = hypFrame, fixedFormSimple, random = randomForm, ...), silent = silent)
        pVal <- if (is(intModel, "try-error")) {
            NA
        } else {
            anovaRes <- anova(xyModel, intModel, test = "Chisq")
            anovaRes$"Pr(>Chi)"[2]
        }
        list("pVal" = pVal, "coef" = coef(xyModel))
    }
    return(out)
}
