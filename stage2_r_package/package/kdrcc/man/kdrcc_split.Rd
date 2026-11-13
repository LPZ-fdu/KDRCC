\name{kdrcc_split}
\alias{kdrcc_split}
\title{Partition Complete Matched Sets}
\description{
Partition Complete Matched Sets.
}
\usage{
kdrcc_split (data, strata, learning = 0.25, tuning = 0, seed = 27182818L)
}
\arguments{
\item{data}{Matched long data, or exposure combinations for effect evaluation. Each estimation stratum contains one event and at least one reference date.}
\item{strata}{Name of the matched-set identification column.}
\item{learning}{Learning matched-set count or a fraction between zero and one.}
\item{tuning}{Selection matched-set count or fraction; zero gives a two-way split.}
\item{seed}{Integer seed for reproducible independent random streams.}
}
\details{
The fitted response uses centered, whitened Gaussian features. Complete matched sets are preserved throughout estimation and validation.
}
\value{The configured object, numerical result, data frame, diagnostic list or graphical object described by the function.}
\author{Peizheng Li}
