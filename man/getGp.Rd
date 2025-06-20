% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getGp.R
\name{getGp}
\alias{getGp}
\title{Helper function to get gene pair from a vector or list}
\usage{
getGp(x, gp, drop = TRUE, Collapse = "--", notFoundReturn = NULL)
}
\arguments{
\item{x}{The object in which to look}

\item{gp}{A character string describing the gene pair}

\item{drop}{A boolean, should matrix attributes be dropped in [] subsetting}

\item{Collapse}{The character separating the gene pair}

\item{notFoundReturn}{value to return if element is not found}
}
\value{
The element sought
}
\description{
When provided with argument "geneA--geneB", looks for this gene pair as well as for
"geneB--geneA" in the provided object.
}
\examples{
mat <- t(cbind(
    "gene1--gene2" = c(1, 2),
    "gene1--gene3" = c(2, 3)
))
getGp(mat, "gene3--gene1")
}
