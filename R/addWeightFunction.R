#' Build a weighting function based on the data
#'
#' @param hypFrame The hyperframe with all the data
#' @param designVars A character vector containing all design factors
#' (both fixed and random), that are also present as variables in hypFrame
#' @param ... Additional arguments passed on to the scam::scam function
#' @param maxObs,maxFeatures The maximum number of observations respectively
#' features for fitting the weighting function. See details
#' @param pis The PIs for which weighting functions are constructed
#' @param lowestLevelVar The design variable at the lowest level of nesting,
#' often separating technical replicates. The conditional variance is calculated
#' within the groups of PIs defined by this variable
#' @details For computational and memory reasons,
#'  for large datasets the trend fitting
#' is restricted to a subset of the data through the maxObs and maxFeatures
#' parameters.
#' Provide either 'designVars' or 'lowestLevelVar'. In the latter case, the
#' design variables are set to all imageVars in the hypFrame object except
#' lowestLevelVar. When the PI is calculated on the cell level,
#' this the lowest nesting level so inputs will be ignored for these PIs
#'
#' @return A fitted scam object, which can be used to
#' @details The scam functions fits a decreasing spline of the variance as a
#'  function of the number of observations
#' @importFrom scam scam
#' @importFrom stats formula
#' @importFrom BiocParallel bplapply
#' @importFrom Rfast rowAll
#' @export
#' @seealso \link{buildDfMM}
#' @examples
#' data(Yang)
#' hypYang <- buildHyperFrame(Yang,
#'   coordVars = c("x", "y"),
#'   imageVars = c("day", "root", "section")
#' )
#' # Fit a subset of features to limit computation time
#' yangPims <- estPims(hypYang,
#'   pis = c("nn", "nnPair"),
#'   features = attr(hypYang, "features")[1:10]
#' )
#' # Build the weighting function for all PIs present
#' yangObj <- addWeightFunction(yangPims, designVars = c("day", "root"))
#' # Alternative formulation with 'lowestLevelVar'
#' yangObj2 <- addWeightFunction(yangPims, lowestLevelVar = "section",
#' pi = "nn")
addWeightFunction <- function(hypFrame, pis = attr(hypFrame, "pis")[
                                !attr(hypFrame, "pis") %in% c("edge", "midpoint")
                              ],
                              designVars, lowestLevelVar, maxObs = 1e5,
                              maxFeatures = 5e2,
                              ...) {
  if (any(pis %in% c("edge", "midpoint"))) {
    stop("Calculating weight matrices for distances to fixed points is ", "unnecessary as they are independent.
             Simply proceed with fitting the model on the", " indiviual evaluations of the B-function.")
  }
  if (missing(designVars) && missing(lowestLevelVar) &&
    any(!grepl("Cell", pis))) {
    stop("Provide either designVars or lowestLevelVar for measures not on cell level!")
  }
  if (!missing(designVars) && !missing(lowestLevelVar)) {
    stop("Provide either designVars or lowestLevelVar, both not both!")
  }
  allVars <- c(attr(hypFrame, "imageVars"), attr(hypFrame, "cellVars"))
  if (missing(designVars)) {
    if (length(lowestLevelVar) != 1) {
      stop("lowestLevelVar must have length 1")
    }
    if (!(lowestLevelVar %in% allVars)) {
      stop("lowestLevelVar must be part of the design,
            see attr(hypFrame, 'imageVars') and attr(hypFrame, 'cellVars')")
    }
    designVars <- setdiff(attr(hypFrame, "imageVars"), lowestLevelVar)
    # If missing take everything but the ones that are obviously
    # no design variables
  }
  if (any(idMissing <- !(designVars %in% c(
    colnames(hypFrame),
    attr(hypFrame, "cellVars")
  )))) {
    stop(
      "Design variables\n", designVars[idMissing],
      "\nnot found in hypFrame object"
    )
  }
  pis <- match.arg(pis,
    choices = c(
      "nn", "nnPair", "allDist", "allDistPair",
      "nnCell", "allDistCell", "nnPairCell", "allDistPairCell"
    ),
    several.ok = TRUE
  )
  if (is.null(hypFrame$pimRes)) {
    stop(
      "No pims found in the hyperframe.",
      "First estimate them using the estPims() function."
    )
  }
  Wfs <- bplapply(pis, function(pi) {
    piListName <- if (pairId <- grepl("Pair", pi)) "biPIs" else "uniPIs"
    subListName <- if (cellId <- grepl("Cell", pi)) "windowDists" else "pointDists"
    features <- attr(hypFrame, "featuresEst")
    piList <- lapply(hypFrame$pimRes, function(x) {
      lapply(x[[piListName]], function(y) {
        y[[subListName]][[pi]]
      })
    })
    if (pairId) {
      features <- apply(combn(features, 2), 2, paste, collapse = "--")
    }
    if (cellId) {
      ordDesign <- seq_len(nrow(hypFrame))
    } else {
      designVec <- apply(as.data.frame(hypFrame[, designVars, drop = FALSE]),
        1, paste,
        collapse = "_"
      )
      ordDesign <- order(designVec) # Ensure correct ordering for tapply
    }
    varEls <- lapply(features, function(gene) {
      geneSplit <- if (pairId) {
        sund(gene)
      } else {
        gene
      }
      if (cellId) {
        # If cellId, there is no tapply, cells are the lowest level anyway
        tmp <- lapply(ordDesign, function(x) {
          if (!all((if (pairId) gene else geneSplit) %in% names(piList[[x]]))) {
            return(NULL)
          }
          CellGene <- table(
            marks(hypFrame$ppp[[x]])$cell,
            marks(hypFrame$ppp[[x]])$gene
          )[, geneSplit, drop = FALSE]
          deps <- vapply(piList[[x]][if (pairId) gene else geneSplit],
            FUN.VALUE = double(sum(id <- rowAll(CellGene > 0))), function(y) {
              xx <- unlist(y)
              if (sum(!is.na(xx)) >= 2) {
                (xx - mean(xx, na.rm = TRUE))^2
              } else {
                rep_len(NA, length(xx))
              }
            }
          )
          cbind("quadDeps" = deps, rowSort(CellGene[id, , drop = FALSE]))
        })
        if (all(vapply(tmp, FUN.VALUE = logical(1), is.null))) {
          return(NULL)
        } else {
          out <- t(Reduce(tmp, f = rbind))
        }
      } else {
        tmp <- vapply(piList[ordDesign], FUN.VALUE = double(1), function(x) {
          if (is.null(baa <- getGp(x, gene))) NA else baa
        })
        quadDeps <- unlist(tapply(tmp, designVec, function(x) {
          if (sum(!is.na(x)) >= 2) {
            (x - mean(x, na.rm = TRUE))^2
          } else {
            rep_len(NA, length(x))
          }
        })) # The quadratic departures from the conditional mean
        tabEntries <- vapply(hypFrame$tabObs[ordDesign],
          FUN.VALUE = double(if (pairId) 2 else 1), function(x) {
            if (pairId) {
              if (all(geneSplit %in% names(x))) sort(x[geneSplit]) else rep(NA, 2)
            } else {
              x[geneSplit]
            }
          }
        )
        out <- rbind("quadDeps" = quadDeps, tabEntries)
      }
      rownames(out)[-1] <- if (pairId) c("minP", "maxP") else "NP"
      out
    })
    varElMat <- matrix(unlist(varEls),
      ncol = if (pairId) 3 else 2,
      byrow = TRUE, dimnames = list(
        NULL,
        c("quadDeps", if (pairId) c("minP", "maxP") else "NP")
      )
    )
    # Faster than cbind
    varElMat <- varElMat[!is.na(varElMat[, "quadDeps"]) &
      varElMat[, "quadDeps"] != 0, ]
    if (nrow(varElMat) > maxObs) {
      varElMat <- varElMat[sample(nrow(varElMat), maxObs), ]
    }
    scamForm <- formula(paste("log(quadDeps) ~", if (pairId) {
      "s(log(maxP), bs = 'mpd') + s(log(minP), bs = 'mpd')"
    } else {
      "s(log(NP), bs = 'mpd')"
    }))
    scamMod <- scam(scamForm, data = data.frame(varElMat), ...)
    attr(scamMod, "pis") <- pi
    return(scamMod)
  })
  names(Wfs) <- pis
  return(list("hypFrame" = hypFrame, "Wfs" = Wfs, "designVars" = designVars))
}
