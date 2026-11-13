\name{kdrcc_control}
\alias{kdrcc_control}
\title{Configure Kernel Learning and Estimation}
\description{
Configure Kernel Learning and Estimation.
}
\usage{
kdrcc_control (centers = 144L, kappa = 1e-04, lambda = 1, folds = 3L, coverage = 0.25, 
    seed = 27182818L, workers = 1L, verbose = FALSE, maxit = 70L, tolerance = 1e-07, 
    step_tolerance = 1e-06, rho_bounds = c(0.01, 100), initial_scale = 0.0625, kernel_maxit = 80L, 
    kernel_continuations = 1L, kernel_tolerance = 1e-05, penalty_rows = 256L, evaluation_points = 1000L, 
    confidence = 0.95, tail_accuracy = 1e-08, tail_terms = 100000L, tail_retry_terms = 1000000L, 
    tail_draws = 199999L, keep_data = FALSE)
}
\arguments{
\item{centers}{Number of Gaussian centers M.}
\item{kappa}{Positive ridge coefficient in independent kernel learning.}
\item{lambda}{Nonnegative effective mixed-derivative penalty in final estimation.}
\item{folds}{Number of matched-set auxiliary cross-fitting folds.}
\item{coverage}{Fraction of centers assigned to farthest-point coverage, subject to at least two clustering centers.}
\item{seed}{Integer seed for reproducible independent random streams.}
\item{workers}{Number of independent PSOCK processes used for tuning profiles.}
\item{verbose}{Logical; print timestamped computational progress.}
\item{maxit}{Maximum coefficient iterations.}
\item{tolerance}{Maximum-gradient convergence tolerance.}
\item{step_tolerance}{Relative coefficient-step convergence tolerance.}
\item{rho_bounds}{Positive lower and upper distance multipliers relative to the geometric baseline.}
\item{initial_scale}{Initial distance multiplier relative to the geometric baseline.}
\item{kernel_maxit}{Outer log-distance optimization iterations per block.}
\item{kernel_continuations}{Additional outer optimization blocks.}
\item{kernel_tolerance}{Projected log-distance gradient tolerance.}
\item{penalty_rows}{Number of independently learned exposure rows used in the mixed-derivative matrix.}
\item{evaluation_points}{Maximum learning-sample exposure rows saved for default evaluation grids.}
\item{confidence}{Default pointwise confidence level.}
\item{tail_accuracy}{Requested Davies absolute integration accuracy.}
\item{tail_terms}{Initial Davies integration budget.}
\item{tail_retry_terms}{Retry integration budget.}
\item{tail_draws}{Minimum reference-distribution draws for Monte Carlo fallback.}
\item{keep_data}{Logical; retain estimation data in the fitted object.}
}
\details{
The fitted response uses centered, whitened Gaussian features. Complete matched sets are preserved throughout estimation and validation.
}
\value{The configured object, numerical result, data frame, diagnostic list or graphical object described by the function.}
\author{Peizheng Li}
