\name{predict.kdrcc}
\alias{predict.kdrcc}
\title{Methods for KDR-CC Objects}
\description{
Standard extraction, printing and plotting methods for fitted KDR-CC objects.
}
\usage{
\method{predict}{kdrcc} (object, newdata, reference = NULL, type = c("log_or", "or"), interval = TRUE, 
    level = object$specification$control$confidence, ...)
}
\arguments{
\item{object}{A fitted kdrcc object.}
\item{newdata}{Exposure combinations with the fitted exposure column names.}
\item{reference}{Complete reference exposure vector, one-row data frame, or row-specific reference combinations for prediction.}
\item{type}{Reported prediction scale: log_or or or. Both scales are retained in the returned data frame.}
\item{interval}{Logical; calculate coefficient-covariance-based pointwise intervals.}
\item{level}{Pointwise confidence level.}
\item{...}{Additional graphical or generic-method arguments.}
}
\details{
The fitted response uses centered, whitened Gaussian features. Complete matched sets are preserved throughout estimation and validation. The criterion is one half of the squared doubly robust moment norm plus lambda/2 times the mixed-derivative quadratic penalty. Pointwise covariance uses the exact criterion curvature and full finite-dimensional auxiliary-estimation contributions. Its population target is the regularized working response conditional on independently learned and selected objects.
}
\value{The configured object, numerical result, data frame, diagnostic list or graphical object described by the function.}
\author{Peizheng Li}
