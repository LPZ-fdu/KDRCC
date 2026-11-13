\name{kdrcc_nuisance}
\alias{kdrcc_nuisance}
\title{Specify Finite-Dimensional Auxiliary Models}
\description{
Specify Finite-Dimensional Auxiliary Models.
}
\usage{
kdrcc_nuisance (g_exposure = "spline_kernel", f_method = c("basis", "cells"), g_df = 4L, 
    kernel_rank = 24L, f_df = 3L, f_formula = NULL)
}
\arguments{
\item{g_exposure}{Auxiliary exposure dictionary: spline_kernel, spline, kernel, a one-sided formula, or a function returning a finite numeric design matrix.}
\item{f_method}{Conditional feature mean design: basis or cells. Cells use the ordered raw adjustment-covariate combinations.}
\item{g_df}{Spline design size for each auxiliary exposure coordinate.}
\item{kernel_rank}{Maximum whitened kernel columns in the auxiliary exposure dictionary.}
\item{f_df}{Spline design size per ordered conditional-mean coordinate.}
\item{f_formula}{Optional one-sided formula on ordered adjustment columns, named plus.V and minus.V for covariate V.}
}
\details{
The fitted response uses centered, whitened Gaussian features. Complete matched sets are preserved throughout estimation and validation. Both auxiliary models have fixed finite-dimensional designs learned from the independent sample. The g fit is a joint conditional regression of exposure basis and adjustment covariates; its adjustment coefficients define the ordered g function. The f fit is weighted least squares among event-second pairs. The cells option uses a saturated ordered-covariate indicator design. The basis option includes an intercept and finite spline coordinates or f_formula. The covariance includes their full influence contributions.
}
\value{The configured object, numerical result, data frame, diagnostic list or graphical object described by the function.}
\author{Peizheng Li}
