% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/es.fsq.overall.R
\name{es.fsq.overall}
\alias{es.fsq.overall}
\title{Cohen f-squared effect size for overall F test in multiple linear regression}
\usage{
es.fsq.overall(Rsq = 0.02)
}
\arguments{
\item{Rsq}{The squared sample multiple correlation coefficient.}
}
\value{
A list of the arguments and the f^2 effect size.
}
\description{
Computes the f-squared (f^2) effect size for an overall F test in a multiple linear regression model
based on the model R^2 (Rsq). Based on Cohen (1988).
}
\details{
Cohen J (1988) Statistical Power Analysis for the Behavioral Sciences, 2nd edition.
Lawrence Erlbaum Associates, Hillsdale, New Jersey
}
\examples{
es.fsq.overall(Rsq = 0.02)
}
