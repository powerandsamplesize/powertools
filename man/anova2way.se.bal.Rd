% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anova2way.se.bal.R
\name{anova2way.se.bal}
\alias{anova2way.se.bal}
\title{Power calculation for test of simple effect for two-way balanced analysis of variance}
\usage{
anova2way.se.bal(
  n = NULL,
  mmatrix = NULL,
  cmatrix = NULL,
  sd = 1,
  Rsq = 0,
  ncov = 0,
  alpha = 0.05,
  power = NULL,
  sides = 2,
  v = FALSE
)
}
\arguments{
\item{n}{The sample size per cell.}

\item{mmatrix}{A matrix of cell means (see example).}

\item{cmatrix}{A matrix of contrast coefficients (see example).}

\item{sd}{The estimated standard deviation within each cell; defaults to 1.}

\item{Rsq}{The estimated R^2 for regressing the outcome on the covariates; defaults to 0.}

\item{ncov}{The number of covariates adjusted for in the model; defaults to 0.}

\item{alpha}{The significance level (type 1 error rate); defaults to 0.05.}

\item{power}{The specified level of power.}

\item{sides}{Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.}

\item{v}{Either TRUE for verbose output or FALSE (default) to output computed argument only.}
}
\value{
A list of the arguments (including the computed one).
}
\description{
Conducts power and sample size calculations for a test of a simple effect in a two-way
balanced (equal cell sizes) ANOVA. A "simple effect" is a contrast among the cell means.
For a test of a contrast in an unbalanced (unequal
cell sizes) two-way ANOVA, see anova2way.se.unbal. For a test of contrast among
factor levels, see anova2way.c.bal.
}
\examples{
mmatrix <- matrix(c(9.3, 8.9, 8.5, 8.7, 8.3, 7.3), nrow = 2, byrow = TRUE)
cmatrix <- matrix(c(-1, 0, 0, 1, 0, 0), nrow = 2, byrow = TRUE)
anova2way.se.bal(n = 30, mmatrix = mmatrix, cmatrix = cmatrix, sd = 2, alpha = 0.025, sides = 1)
}
