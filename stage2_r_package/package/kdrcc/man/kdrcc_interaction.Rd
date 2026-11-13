\name{kdrcc_interaction}
\alias{kdrcc_interaction}
\title{Evaluate and Test Four-Point Interactions}
\description{
Evaluate and Test Four-Point Interactions.
}
\usage{
kdrcc_interaction (object, pairs, values = NULL, background = NULL, points = 9L, alpha = 0.05, 
    adjust = "holm", level = 0.95)
}
\arguments{
\item{object}{A fitted kdrcc object.}
\item{pairs}{Two-column matrix of exposure names specifying interaction pairs.}
\item{values}{Concentration vector, or a named list of concentration vectors.}
\item{background}{Complete fixed background exposure combination.}
\item{points}{Number of concentrations per coordinate before adding the anchor.}
\item{alpha}{Family significance level.}
\item{adjust}{Multiplicity adjustment accepted by stats::p.adjust; the default is holm.}
\item{level}{Pointwise confidence level.}
}
\details{
The fitted response uses centered, whitened Gaussian features. Complete matched sets are preserved throughout estimation and validation. The four-point log-OR contrast is h(x,y)-h(x,b)-h(b,y)+h(b,b). The overall zero-grid null uses a quadratic-form reference; pointwise signs and two-sided Wald p values are returned. Overall pairs and eligible points form separate correction families.
}
\value{A kdrcc_interaction object with overall pair tests and data containing pointwise contrasts, intervals, p values and directions.}
\author{Peizheng Li}
