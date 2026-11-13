\name{kdrcc_single}
\alias{kdrcc_single}
\title{Evaluate Fixed-Background Single-Exposure Curves}
\description{
Evaluate Fixed-Background Single-Exposure Curves.
}
\usage{
kdrcc_single (object, exposure = object$specification$exposures, values = NULL, background = NULL, 
    points = 41L, level = 0.95)
}
\arguments{
\item{object}{A fitted kdrcc object.}
\item{exposure}{Exposure names for single-coordinate curves.}
\item{values}{Concentration vector, or a named list of concentration vectors.}
\item{background}{Complete fixed background exposure combination.}
\item{points}{Number of concentrations per coordinate before adding the anchor.}
\item{level}{Pointwise confidence level.}
}
\details{
The fitted response uses centered, whitened Gaussian features. Complete matched sets are preserved throughout estimation and validation.
}
\value{A kdrcc_effect object with a data component containing coordinates, log OR, standard error, pointwise limits and exponentiated OR results.}
\author{Peizheng Li}
