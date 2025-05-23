% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/addCell.R
\name{addCell}
\alias{addCell}
\title{Add cell boundaries and event-wise cell identifiers to a hyperframe.}
\usage{
addCell(
  hypFrame,
  owins,
  cellTypes = NULL,
  findOverlappingOwins = FALSE,
  warnOut = TRUE,
  coords = c("x", "y"),
  verbose = TRUE,
  addCellMarkers = TRUE,
  overwriteCells = FALSE,
  ...
)
}
\arguments{
\item{hypFrame}{A hyperframe}

\item{owins}{A list containing a list of owins per point pattern.
The length of the list must match the length of the hyperframe, and the names must match.
Also lists of geojson objects, coordinate matrices or rois are accepted, see details.}

\item{cellTypes}{A dataframe of cell types and other cell-associated covariates.
If supplied, it must contain a variable 'cell' that is matched with the names of the owins}

\item{findOverlappingOwins}{a boolean, should windows be checked for overlap? Can be computationally intensive.}

\item{warnOut}{a boolean, should warning be issued when points are not
contained in window?}

\item{coords}{The names of the coordinates, if the windows are given as sets of coordinates.}

\item{verbose}{A boolean, should verbose output be printed?}

\item{addCellMarkers}{A boolean, should cell identities be added? Set this to
FALSE if cell identifiers are already present in the data, and you only want
to add windows and centroids.}

\item{overwriteCells}{A boolean, should cells already present in hyperframe be overwritten?}

\item{...}{Further arguments passed onto \link{convertToOwins}}
}
\value{
The hyperframe with cell labels added in the marks of the point patterns
}
\description{
Add the list of the cells and their centroids in the hyperframe,
check in which cell each event lies and add a cell marker.
}
\details{
First the different cells are checked for overlap per point pattern if 'findOverlappingOwins' is TRUE.
If no overlap is found, each event is assigned the cell that it falls into.
Events not belonging to any cell will trigger a warning and be assigned 'NA'.
Cell types and other variables are added to the marks if applicable.
This function employs multithreading through the BiocParallel package. If this leads to excessive memory usage and crashes, try serial processing by
setting register(SerialParam()).
Different formats of windows are allowed, if the corresponding packages are installed.
A dataframe of coordinates or a list of spatstat.geom owins is always allowed, as the necessary packages are required by smoppix.
A 'SpatialPolygonsDataFrame' object is allowed if the 'polycub' package is installed,
and a list of 'ijroi' object or a single 'ijzip' object if the 'RImageJROI' package is installed.
}
\note{
By default, overlap between windows is not checked.
Events are assigned to the first window they fall in. If you are not sure of the quality of the segmentation,
do check your input or set checkOverlap to TRUE, even when this make take time.
}
\examples{
library(spatstat.random)
set.seed(54321)
n <- 1e3 # number of molecules
ng <- 25 # number of genes
nfov <- 3 # Number of fields of view
conditions <- 3
# sample xy-coordinates in [0, 1]
x <- runif(n)
y <- runif(n)
# assign each molecule to some gene-cell pair
gs <- paste0("gene", seq(ng))
gene <- sample(gs, n, TRUE)
fov <- sample(nfov, n, TRUE)
condition <- sample(conditions, n, TRUE)
# construct data.frame of molecule coordinates
df <- data.frame(gene, x, y, fov, "condition" = condition)
# A list of point patterns
listPPP <- tapply(seq(nrow(df)), df$fov, function(i) {
    ppp(x = df$x[i], y = df$y[i], marks = df[i, "gene", drop = FALSE])
}, simplify = FALSE)
# Regions of interest (roi): Diamond in the center plus four triangles
w1 <- owin(poly = list(x = c(0, .5, 1, .5), y = c(.5, 0, .5, 1)))
w2 <- owin(poly = list(x = c(0, 0, .5), y = c(.5, 0, 0)))
w3 <- owin(poly = list(x = c(0, 0, .5), y = c(1, 0.5, 1)))
w4 <- owin(poly = list(x = c(1, 1, .5), y = c(0.5, 1, 1)))
w5 <- owin(poly = list(x = c(1, 1, .5), y = c(0, 0.5, 0)))
hypFrame <- buildHyperFrame(df,
    coordVars = c("x", "y"),
    imageVars = c("condition", "fov")
)
nDesignFactors <- length(unique(hypFrame$image))
wList <- lapply(seq_len(nDesignFactors), function(x) {
    list("w1" = w1, "w2" = w2, "w3" = w3, "w4" = w4, "w5" = w5)
})
names(wList) <- rownames(hypFrame) # Matching names is necessary
hypFrame2 <- addCell(hypFrame, wList)
}
\seealso{
\link{buildHyperFrame}, \link{convertToOwins}
}
