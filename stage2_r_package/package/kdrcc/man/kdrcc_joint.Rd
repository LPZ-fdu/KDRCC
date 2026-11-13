\name{kdrcc_joint}
\alias{kdrcc_joint}
\title{Evaluate a Joint Exposure Path}
\description{
Evaluate a Joint Exposure Path.
}
\usage{
kdrcc_joint (object, quantiles = seq(0.1, 0.9, length.out = 21L), variables = object$specification$exposures, 
    grid = NULL, reference = NULL, level = 0.95)
}
\arguments{
\item{object}{A fitted kdrcc object.}
\item{quantiles}{Exposure quantiles defining a joint path in the saved evaluation distribution.}
\item{variables}{Fitted exposure names defining the varying path or contribution set.}
\item{grid}{Candidate data frame with M, kappa and lambda columns; or a complete exposure grid for joint effects.}
\item{reference}{Complete reference exposure vector, one-row data frame, or row-specific reference combinations for prediction.}
\item{level}{Pointwise confidence level.}
}
\details{
The fitted response uses centered, whitened Gaussian features. Complete matched sets are preserved throughout estimation and validation.
}
\value{A kdrcc_effect object with a data component containing coordinates, log OR, standard error, pointwise limits and exponentiated OR results.}
\author{Peizheng Li}
