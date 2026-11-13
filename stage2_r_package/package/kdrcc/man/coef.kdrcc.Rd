\name{coef.kdrcc}
\alias{coef.kdrcc}
\title{Methods for KDR-CC Objects}
\description{
Standard extraction, printing and plotting methods for fitted KDR-CC objects.
}
\usage{
\method{coef}{kdrcc} (object, ...)
}
\arguments{
\item{object}{A fitted kdrcc object.}
\item{...}{Additional graphical or generic-method arguments.}
}
\details{
The fitted response uses centered, whitened Gaussian features. Complete matched sets are preserved throughout estimation and validation.
}
\value{The configured object, numerical result, data frame, diagnostic list or graphical object described by the function.}
\author{Peizheng Li}
