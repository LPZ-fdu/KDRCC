\name{print.kdrcc_test}
\alias{print.kdrcc_test}
\title{Methods for KDR-CC Objects}
\description{
Standard extraction, printing and plotting methods for fitted KDR-CC objects.
}
\usage{
\method{print}{kdrcc_test} (x, ...)
}
\arguments{
\item{x}{A fitted object or a returned effect, test, tuning or contribution object.}
\item{...}{Additional graphical or generic-method arguments.}
}
\details{
The fitted response uses centered, whitened Gaussian features. Complete matched sets are preserved throughout estimation and validation.
}
\value{The configured object, numerical result, data frame, diagnostic list or graphical object described by the function.}
\author{Peizheng Li}
