\name{vcov.kdrcc}
\alias{vcov.kdrcc}
\title{Methods for KDR-CC Objects}
\description{
Standard extraction, printing and plotting methods for fitted KDR-CC objects.
}
\usage{
\method{vcov}{kdrcc} (object, ...)
}
\arguments{
\item{object}{A fitted kdrcc object.}
\item{...}{Additional graphical or generic-method arguments.}
}
\details{
The fitted response uses centered, whitened Gaussian features. Complete matched sets are preserved throughout estimation and validation. The criterion is one half of the squared doubly robust moment norm plus lambda/2 times the mixed-derivative quadratic penalty. Pointwise covariance uses the exact criterion curvature and full finite-dimensional auxiliary-estimation contributions. Its population target is the regularized working response conditional on independently learned and selected objects.
}
\value{The configured object, numerical result, data frame, diagnostic list or graphical object described by the function.}
\author{Peizheng Li}
