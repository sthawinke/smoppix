#' Extract results from a list of fitted LMMs. For internal use mainly.
#'
#' @param models The models
#' @param hypFrame The original hyperframe
#' @param fixedVars The fixed effects for which the effect is to be reported
#' @param method Multiplicity correction method passed onto p.adjust
#'
#' @importFrom stats p.adjust anova
#' @importFrom methods is
#' @return A list of matrices, all containing estimate, standard error,
#' p-value and adjusted p-value
#' @seealso \link{fitLMMs}, \link[stats]{p.adjust}
extractResults <- function(models, hypFrame, fixedVars = NULL, method = "BH") {
    id <- vapply(models, FUN.VALUE = TRUE, function(x) {
        is(x, "lmerModLmerTest") || is.matrix(x)
    })
    out <- if (!any(id)) {
        return(list("Intercept" = NULL, "fixedEffects" = NULL))
    } else {
      ints <- t(vapply(models, FUN.VALUE = double(3), function(x) {
          if (is.null(x) || is(x, "try-error")) {
              c(Estimate = NA, `Std. Error` = NA, `Pr(>|t|)` = NA)
          } else if(is(x, "lmerModLmerTest")){
              summary(x)$coef["(Intercept)", c(
                  "Estimate", "Std. Error",
                  "Pr(>|t|)"
              )]
          } else {
            getLmWfitPvalues(x)["(Intercept)", c("Estimate", "Std. Error","Pr(>|t|)")]
          }
      }))
      colnames(ints) <- c("Estimate", "SE", "pVal")
      ints[, "Estimate"] <- ints[, "Estimate"] + 0.5
      intMat <- cbind(ints, pAdj = p.adjust(ints[, "pVal"], method = method))[order(ints[
          ,"pVal"]), ]
      # Order by p-value
      AnovaTabs <- loadBalanceBplapply(models[id], anova)
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
              coefObj <- if(is(x, "lmerModLmerTest")) summary(x)$coef else x
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
    return(out)
}
