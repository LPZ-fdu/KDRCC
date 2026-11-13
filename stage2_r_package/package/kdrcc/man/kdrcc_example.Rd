\name{kdrcc_example}
\alias{kdrcc_example}
\title{Generate a Matched Example Dataset}
\description{
Generate a Matched Example Dataset.
}
\usage{
kdrcc_example (n = 800L, seed = 13579L, discrete = FALSE)
}
\arguments{
\item{n}{Number of synthetic matched sets.}
\item{seed}{Integer seed for reproducible independent random streams.}
\item{discrete}{Logical; generate a three-level adjustment covariate.}
}
\details{
The fitted response uses centered, whitened Gaussian features. Complete matched sets are preserved throughout estimation and validation.
}
\value{The configured object, numerical result, data frame, diagnostic list or graphical object described by the function.}
\author{Peizheng Li}
\examples{
d <- kdrcc_example(30, seed = 11)
table(d$case)
}
