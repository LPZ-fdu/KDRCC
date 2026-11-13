\name{kdrcc_tune}
\alias{kdrcc_tune}
\title{Select a Configuration by Independent Grouped Validation}
\description{
Select a Configuration by Independent Grouped Validation.
}
\usage{
kdrcc_tune (data, learning_data, case = "case", strata = "set_id", exposures, adjust = ~1, 
    nuisance = kdrcc_nuisance(), grid = expand.grid(M = c(96L, 144L), kappa = c(3e-05, 
        1e-04), lambda = c(0.3, 1, 3)), control = kdrcc_control(), validation_folds = 3L, 
    fallback = TRUE)
}
\arguments{
\item{data}{Matched long data reserved for configuration selection. Each stratum contains one event and at least one reference date.}
\item{learning_data}{Independent complete matched sets for kernel and design learning.}
\item{case}{Name of the zero/one event column.}
\item{strata}{Name of the matched-set identification column.}
\item{exposures}{Unique names of numeric exposure columns.}
\item{adjust}{One-sided formula for time-varying adjustment covariates. Transformations and factor levels are learned independently.}
\item{nuisance}{Auxiliary-model design specification returned by kdrcc_nuisance().}
\item{grid}{Candidate data frame with M, kappa and lambda columns. M is the number of candidate kernel centers.}
\item{control}{Computation and regularization settings returned by kdrcc_control().}
\item{validation_folds}{Number of outer matched-set validation folds.}
\item{fallback}{Logical. If no candidate passes all fold checks, warn and use the supplied control configuration. Set to FALSE to require a validated candidate.}
}
\details{
The fitted response uses centered, whitened Gaussian features. Each kernel is learned on independent matched sets. The learner excludes kernel features with negligible within-set variation or numerical linear dependence and retains the resulting basis for validation and final estimation.

Each outer training fold uses matched-set nuisance cross-fitting and computes the coefficient sandwich covariance. A candidate is eligible when every training fold yields a usable covariance and every held-out fold yields a finite conditional negative log likelihood. Eligible candidates are ranked by held-out loss per matched set. Ties follow ascending M, kappa and lambda. Final effects use a fresh estimation sample, and the covariance is computed again on that sample. The workers control sets the number of independent kernel profiles that can run concurrently.
}
\value{A kdrcc_tuning object containing the candidate grid with retained dimensions and eligibility, fold scores and covariance diagnostics, selected control, and learned kernel.}
\author{Peizheng Li}
