#' Extract results from a list of fitted LMMs
#'
#' @param models The models
#' @param hypFrame The original hyperframe
#' @param fixedVars The fixed effects for which the effect is to be reported
#' @param subSet The name of the subset to be extracted, either PI or Moran's I
#' @param method Multiplicity correction method passed onto p.adjust
#'
#' @importFrom stats p.adjust
#' @importFrom methods is
#' @return A list of matrices, all containing estimate, standard error,
#' p-value and ajdusted p-value
#' @seealso \link{fitLMMs}, \link{p.adjust}
extractResults <- function(
    models, hypFrame, subSet = "piMod", fixedVars = NULL,
    method = "BH") {
    id <- vapply(models, FUN.VALUE = TRUE, function(x) {
        is(x[[subSet]], "lmerModLmerTest") || is(x[[subSet]], "lm")
    })
    if (!any(id)) {
        return(list("Intercept" = NULL, "fixedEffects" = NULL))
    }
    ints <- t(vapply(models, FUN.VALUE = double(3), function(x) {
        if (is.null(x[[subSet]]) || is(x[[subSet]], "try-error")) {
            c(Estimate = NA, `Std. Error` = NA, `Pr(>|t|)` = NA)
        } else {
            summary(x[[subSet]])$coef["(Intercept)", c(
                "Estimate", "Std. Error",
                "Pr(>|t|)"
            )]
        } # Identical code for lmerTest or lm
    }))
    colnames(ints) <- c("Estimate", "SE", "pVal")
    ints[, "Estimate"] <- ints[, "Estimate"] + switch(subSet,
        piMod = 0.5,
        moranMod = 0
    )
    intMat <- cbind(ints, pAdj = p.adjust(ints[, "pVal"], method = method))[order(ints[
        ,
        "pVal"
    ]), ]
    # Order by p-value
    AnovaTabs <- lapply(models[id], function(x) anova(x[[subSet]]))
    fixedOut <- lapply(fixedVars, function(Var) {
        if(discr <- Var %in% getDiscreteVars(hypFrame)){
            unVals <- if (Var %in% getEventVars(hypFrame)) {
                unique(unlist(lapply(hypFrame$ppp, function(x) {
                    marks(x, drop = FALSE)[[Var]]
                })))
            } else {
                unique(hypFrame[[Var]])
            }
            emptyCoef <- rep_len(NA, length(unVals))
            names(emptyCoef) <- paste0(Var, unVals) # Prepare empty coefficient
        } else {
            emptyCoef <- double(1)
            names(emptyCoef) <- Var
        }
        pVal <- vapply(AnovaTabs, FUN.VALUE = double(1), function(x) x[Var, "Pr(>F)"])
        coefs <- lapply(models[id], function(x) {
            # Prepare the empty coefficient vector with all levels present. If
            # outcome is NA for all levels, the factor level gets dropped,
            # causing problems downstream.
            coefObj <- summary(x[[subSet]])$coef
            if(discr){
                rn <- intersect(names(emptyCoef), rownames(coefObj))
                emptyCoef[rn] <- coefObj[rn, "Estimate"]
                emptyCoef
            } else {
                coefObj[Var, "Estimate"]
            }
        })
        coefMat <- matrix(unlist(coefs), byrow = TRUE, nrow = length(pVal))
        colnames(coefMat) <- names(emptyCoef)
        if (ncol(coefMat) == 2) {
            # If two levels, the two coefficients are each other's opposite
            # in sum coding, and NA can be replaced by the negative. For more than two
            # levels, it's more complicated as NA may indicate that the parameter was not estimated
            oneNA <- which(rowSums(naMat <- is.na(coefMat)) == 1)
            for (i in oneNA) {
                coefMat[i, naMat[i, ]] <- -coefMat[i, !naMat[i, ]]
            }
        }
        cbind(coefMat, pVal = pVal, pAdj = p.adjust(pVal, method = method))[order(pVal), ]
    })
    names(fixedOut) <- fixedVars
    list(Intercept = intMat, fixedEffects = fixedOut)
}
