\name{kdrcc}
\alias{kdrcc}
\title{Fit Kernel-Based Doubly Robust Case-Crossover Regression}
\description{
Fit Kernel-Based Doubly Robust Case-Crossover Regression.
}
\usage{
kdrcc (data, case = "case", strata = "set_id", exposures, adjust = ~1, learning_data = NULL, 
    nuisance = kdrcc_nuisance(), control = kdrcc_control(), reference = NULL, kernel = NULL, 
    fold_id = NULL)
}
\arguments{
\item{data}{Matched long data, or exposure combinations for effect evaluation. Each estimation stratum contains one event and at least one reference date.}
\item{case}{Name of the zero/one event column.}
\item{strata}{Name of the matched-set identification column.}
\item{exposures}{Unique names of numeric exposure columns.}
\item{adjust}{One-sided formula for time-varying adjustment covariates. Transformations and factor levels are learned independently.}
\item{learning_data}{Independent complete matched sets for kernel and design learning.}
\item{nuisance}{Auxiliary-model design specification returned by kdrcc_nuisance().}
\item{control}{Computation and regularization settings returned by kdrcc_control().}
\item{reference}{Complete reference exposure vector, one-row data frame, or row-specific reference combinations for prediction.}
\item{kernel}{A learned kdrcc_kernel object, for example the kernel returned by kdrcc_tune(). Its geometry and auxiliary specification are reused.}
\item{fold_id}{Optional cross-fitting labels in the order of first appearance of estimation strata.}
}
\details{
The fitted response uses centered, whitened Gaussian features. Complete matched sets are preserved throughout estimation and validation. The independent learning sample determines the retained kernel basis by screening features for negligible within-set variation and numerical linear dependence. The criterion is one half of the squared doubly robust moment norm plus lambda/2 times the mixed-derivative quadratic penalty. Pointwise covariance uses the exact criterion curvature and full finite-dimensional auxiliary-estimation contributions. Its population target is the regularized working response in the retained kernel basis, conditional on independently learned and selected objects.
}
\value{A fitted object of class kdrcc containing coefficients, covariance, learned kernel, auxiliary fits, reference, diagnostics and reproducibility records.}
\author{Peizheng Li}
