\name{kdrcc_variables}
\alias{kdrcc_variables}
\title{Test Derivative-Based Variable Participation}
\description{
Test Derivative-Based Variable Participation.
}
\usage{
kdrcc_variables (object, data = object$evaluation, alpha = 0.05, adjust = "holm")
}
\arguments{
\item{object}{A fitted kdrcc object.}
\item{data}{Matched long data, or exposure combinations for effect evaluation. Each estimation stratum contains one event and at least one reference date.}
\item{alpha}{Family significance level.}
\item{adjust}{Multiplicity adjustment accepted by stats::p.adjust; the default is holm.}
}
\details{
The fitted response uses centered, whitened Gaussian features. Complete matched sets are preserved throughout estimation and validation. For each standardized coordinate the statistic is beta' Q beta, where Q is the evaluation mean of feature-derivative outer products. The working null is Q^(1/2) beta_lambda = 0. Its centered Gaussian reference is a weighted chi-square mixture. Holm adjustment covers the candidate family. For a finite positive-distance Gaussian space, a zero coordinate derivative on an open set implies a zero centered function; a full-rank finite-grid Q likewise gives the zero-coefficient working null. Generating-label recognition and false inclusion are evaluated in the separate simulation archive.
}
\value{A kdrcc_test object with statistic, p_value, p_adjusted, selected, numerical accuracy and algorithm status for each exposure.}
\author{Peizheng Li}
