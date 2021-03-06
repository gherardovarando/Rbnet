% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/2-bnet-fit-inference.R
\name{bnetPar}
\alias{bnetPar}
\title{Set fitting parameters}
\usage{
bnetPar(...)
}
\arguments{
\item{...}{see details.}
}
\description{
Set appropriate parameters to be used by fitting
functions.
}
\details{
This function set parameters used for estimation.
          \code{N}:  The number of knots in every dimensions.
          \code{order}: The order of the B-spline in every dimensions.
          \code{alpha}: The exponent to compute the number of knots.
         \code{knotsMethod}: \code{"uniform"} or \code{"quantiles"} knots.
}

