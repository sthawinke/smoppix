% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convertToOwins.R
\name{convertToOwins}
\alias{convertToOwins}
\title{Convert windows in different formats to spatstat owin format}
\usage{
convertToOwins(windows, namePPP, coords, ...)
}
\arguments{
\item{windows}{The list of windows}

\item{namePPP}{the name of the point pattern, will be added to the cell names}

\item{coords}{The names of the coordinates, if the windows are given as sets of coordinates}

\item{...}{Further arguments passed onto read functions}
}
\value{
A list of owins
}
\description{
Convert a list of differently formatted windows to owins, for
addition to a hyperframe
}