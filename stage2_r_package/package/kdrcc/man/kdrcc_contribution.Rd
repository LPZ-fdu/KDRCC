\name{kdrcc_contribution}
\alias{kdrcc_contribution}
\title{Allocate Baseline Shapley Contributions}
\description{
Allocate Baseline Shapley Contributions.
}
\usage{
kdrcc_contribution (object, data = object$evaluation, variables = object$specification$exposures, 
    reference = NULL, chunk = 128L)
}
\arguments{
\item{object}{A fitted kdrcc object.}
\item{data}{Matched long data, or exposure combinations for effect evaluation. Each estimation stratum contains one event and at least one reference date.}
\item{variables}{Fitted exposure names defining the varying path or contribution set.}
\item{reference}{Complete reference exposure vector, one-row data frame, or row-specific reference combinations for prediction.}
\item{chunk}{Number of evaluation combinations per contribution block.}
}
\details{
The fitted response uses centered, whitened Gaussian features. Complete matched sets are preserved throughout estimation and validation. Only the specified contribution coordinates change from reference to current levels; the remaining coordinates stay at their reference background. Gaussian separability permits exact Gauss-Legendre integration of the Shapley polynomial. Shares normalize mean absolute allocations within the specified set.
}
\value{A kdrcc_contribution object containing allocations, mean absolute importance, signed averages, share_percent and the allocation identity error.}
\author{Peizheng Li}
