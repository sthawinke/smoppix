% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{named.contr.sum}
\alias{named.contr.sum}
\title{A version of contr.sum that retains names, a bit controversial but also clearer}
\usage{
named.contr.sum(x, ...)
}
\arguments{
\item{x, ...}{passed on to contr.sum}
}
\value{
The matrix of contrasts
}
\description{
A version of contr.sum that retains names, a bit controversial but also clearer
}
\note{
After
https://stackoverflow.com/questions/24515892/r-how-to-contrast-code-factors-and-retain-meaningful-labels-in-output-summary
}
