% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getDesignVars.R
\name{getDesignVars}
\alias{getDesignVars}
\alias{getPPPvars}
\alias{getEventVars}
\title{Return all design variables, both at the level of the point pattern and the
level of the event}
\usage{
getDesignVars(x)

getPPPvars(
  x,
  exclude = c("tabObs", "centroids", "owins", "ppp", "pimRes", "image", "inSeveralCells",
    "nuclei")
)

getEventVars(x, exclude = c("x", "y", "z"))
}
\arguments{
\item{x}{The results list, output from estPis}

\item{exclude}{variables to exclude}
}
\value{
A vector of design variables

A vector of variables

A vector of variables
}
\description{
Return all design variables, both at the level of the point pattern and the
level of the event

Extract variables from point patterns

Extract variables from events (the marks)
}
\details{
getDesignVars() returns all design variables, \link{getPPPvars} returns
design variables related to the different images and \link{getEventVars} returns
design variables related to the individual events
}
