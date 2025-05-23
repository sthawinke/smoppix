% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotExplore.R
\name{plotExplore}
\alias{plotExplore}
\title{Plot a hyperframe with chosen features highlighted}
\usage{
plotExplore(
  hypFrame,
  features = getFeatures(hypFrame)[seq_len(6)],
  ppps,
  numPps,
  maxPlot = 1e+05,
  Cex = 1,
  plotWindows = !is.null(hypFrame$owins),
  plotPoints = TRUE,
  plotNuclei = !is.null(hypFrame$nuclei),
  piEsts = NULL,
  Xlim = NULL,
  Ylim = NULL,
  Cex.main = 1.1,
  Mar = c(0.5, 0.1, 0.9, 0.1),
  titleVar = NULL,
  piColourCell = NULL,
  palCols = c("blue", "yellow"),
  nucCol = "lightblue",
  border = NULL,
  CexLegend = 1.4,
  CexLegendMain = 1.7,
  Nrow,
  Cols
)
}
\arguments{
\item{hypFrame}{The hyperframe}

\item{features}{A small number of features to be highlighted. Defaults to the first 5.}

\item{ppps}{The rownames or indices of the point patterns to be plotted.
Defaults to maximum 99.}

\item{numPps}{The number of point patterns with highest expression to be shown.
Ignored is pps is given.}

\item{maxPlot}{The maximum number of events plotted per point pattern}

\item{Cex, Cex.main}{Point and title expansion factors, repsectively}

\item{plotWindows}{A boolean, should windows be plotted too?}

\item{plotPoints}{A boolean, should the molecules be plotted as points?}

\item{plotNuclei}{A boolean, should the nuclei be plotted?}

\item{piEsts}{Set of PI estimates, returned by estPis}

\item{Xlim, Ylim}{plotting limits}

\item{Mar}{the margins}

\item{titleVar}{Image variable to be added to the title}

\item{piColourCell}{PI by which to colour the cell}

\item{palCols}{Two extremes of the colour palette for colouring the cells}

\item{nucCol}{The colour for the nucleus window}

\item{border}{Passed on to plot.owin, and further to graphics::polygon}

\item{CexLegend, CexLegendMain}{Expansion factor for the legend and its title
respectively}

\item{Nrow}{Number of rows of the facet plot. Will be calculated if missing.}

\item{Cols}{colours vector named by features. If missing a default palette is used}
}
\value{
Plots a facet of point patterns to output
}
\description{
All points of the hyperframe are plotted in grey, with a subset of features highlighted in colour.
A selection of point patterns is plotted that fit in the window, by default the first six.
This function is meant for exploratory purposes as well as for visual confirmation of findings.
}
\details{
When cell-specific PIs are calculated ("nnCell', "nnCellPair", "edge", "centroid"),
the cells can be coloured by them to investigate their spatial distribution,
for instance those discovered through Moran's I statistic.
The colour palette is taken from the output of palette(),
so set that one to change the colour scheme.
}
\note{
palCols sets the pseudo-continuous scale to colour cells.
}
\examples{
example(buildHyperFrame, "smoppix")
plotExplore(hypYang)
plotExplore(hypYang, titleVar = "day")
plotExplore(hypYang, features = c("SmRBRb", "SmTMO5b", "SmWER--SmAHK4f"))
}
