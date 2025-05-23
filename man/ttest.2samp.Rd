% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ttest.2samp.R
\name{ttest.2samp}
\alias{ttest.2samp}
\title{Power calculation for two-sample t test}
\usage{
ttest.2samp(
  n1 = NULL,
  n.ratio = 1,
  delta = NULL,
  sd1 = 1,
  sd.ratio = 1,
  df.method = c("welch", "classical"),
  alpha = 0.05,
  power = NULL,
  sides = 2,
  v = FALSE
)
}
\arguments{
\item{n1}{The sample size for group 1.}

\item{n.ratio}{The ratio n2/n1 between the sample sizes of two groups; defaults to 1 (equal group sizes).}

\item{delta}{DeltaA (true difference mu1 - mu2) - Delta0 (difference under the null). For a noninferiority
or superiority by a margin test,
the margin should be subtracted, that is, delta = DeltaA - Delta0 - margin.}

\item{sd1}{The standard deviation for group 1; defaults to 1 (equal standard deviations in the two groups).}

\item{sd.ratio}{The ratio sd2/sd1 between the standard deviations of the two groups.}

\item{df.method}{Method for calculating the degrees of freedom: "welch" (default) or "classical", which is n1+n2-2.}

\item{alpha}{The significance level (type 1 error rate); defaults to 0.05.}

\item{power}{The specified level of power.}

\item{sides}{Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.}

\item{v}{Either TRUE for verbose output or FALSE (default) to output computed argument only.}
}
\value{
A list of the arguments (including the computed one).
}
\description{
This function computes power and sample size for a two-sample t test.
Unequal sample sizes and/or unequal variances in the two groups are allowed.
Can solve for power, n1, n.ratio, sd1, sd.ratio, delta or alpha.
}
\details{
For a noninferiority or superiority by a margin test, the margin should be subtracted
as part of the specification of delta and a one-sided test should be used. See Crespi (2025)
for specific guidance. For an equivalence test for two means, see the sampleN.TOST function from the
PowerTOST package.
}
\examples{
ttest.2samp(n1 = 50, delta = 2, sd1 = 5, sides = 1)
ttest.2samp(n1 = NULL, n.ratio = 2, delta = 0.5, sd1 = 1, power = 0.8, sides = 2)
ttest.2samp(n1 = 49, n.ratio = 2, delta = 0.5, sd1 = 1, power = NULL, sides = 2)
ttest.2samp(n1 = 25, n.ratio = 3, delta = 3, sd1 = 4, sd.ratio = 1.5, alpha = 0.025, sides = 1)
ttest.2samp(n1 = NULL, delta = 0.5, sd1 = 1, power = 0.8, sides = 2)
}
