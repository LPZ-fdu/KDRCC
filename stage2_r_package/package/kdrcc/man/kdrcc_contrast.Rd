\name{kdrcc_contrast}
\alias{kdrcc_contrast}
\title{Estimate Quantitative Exposure Contrasts}
\description{
Estimate Quantitative Exposure Contrasts.
}
\usage{
kdrcc_contrast (object, newdata, reference = NULL, level = 0.95)
}
\arguments{
\item{object}{A fitted kdrcc object.}
\item{newdata}{Exposure combinations with the fitted exposure column names.}
\item{reference}{Complete reference exposure vector, one-row data frame, or row-specific reference combinations for prediction.}
\item{level}{Pointwise confidence level.}
}
\details{
The fitted response uses centered, whitened Gaussian features. Complete matched sets are preserved throughout estimation and validation. The criterion is one half of the squared doubly robust moment norm plus lambda/2 times the mixed-derivative quadratic penalty. Pointwise covariance uses the exact criterion curvature and full finite-dimensional auxiliary-estimation contributions. Its population target is the regularized working response conditional on independently learned and selected objects.
}
\value{A kdrcc_effect object with a data component containing coordinates, log OR, standard error, pointwise limits and exponentiated OR results.}
\author{Peizheng Li}
