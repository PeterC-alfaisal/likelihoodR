% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/L_efficacy.R
\name{L_efficacy}
\alias{L_efficacy}
\title{Likelihood Support for Efficacy}
\usage{
L_efficacy(a, n, null=0, exp.eff=NULL, L.int=2, alpha=0.05, toler=0.0001, verb=TRUE)
}
\arguments{
\item{a}{the number of affected in control group.}

\item{n}{total number of participants.}

\item{null}{the null value for efficacy, if no effect then it would be 0, default = 0.}

\item{exp.eff}{the expected or hypothesized efficacy, default = NULL.}

\item{L.int}{likelihood interval given as support values, e.g. 2 or 3, default = 2.}

\item{alpha}{the significance level used, 1 - alpha interval calculated, default = 0.05.}

\item{toler}{the desired accuracy using optimise, default = 0.0001.}

\item{verb}{show output, default = TRUE.}
}
\value{
$S.val - support for the observed efficacy versus the null value.

$obs.eff - the observed efficacy.

$null - the null efficacy.

$exp.eff - expected efficacy as specified.

$S.exp.vsObs - support for expected efficacy versus observed.

$S.exp.versus.null - support for the expected efficacy versus the null.

$L.int - the likelihood interval for the observed efficacy.

$S_int - the specified likelihood interval.

$observed - observed numbers affected in control and intervention groups.

$expected - expected numbers according to the null.

$chi.sq - chi-squared statistic.

$p.value - p value associated with chi-squared statistic.

$df - degrees of freedom for chi-squared.

$residuals - the Pearson residuals.

$conf.int - likelihood-based confidence interval according to specified alpha.

$alpha - specified alpha for confidence interval.

$all.err.acc - error accuracy for each application of the optimize function.
}
\description{
This function calculates the support for the efficacy, the likelihood interval
and the likelihood-based confidence interval.
It uses the optimize function to locate desired limits and their error.
}
\examples{
# pfizer covid-19 efficacy 2020
m = L_efficacy(a = 86, n = 94, null=0.8, exp.eff=0.95, L.int=2, alpha=0.05, toler=0.0001)
m

}
\references{
Aitkin, M. et al (1989) Statistical Modelling in GLIM, Clarendon Press, ISBN : 978-0198522041

Cahusac, P.M.B. (2020) Evidence-Based Statistics, Wiley, ISBN : 978-1119549802

Royall, R. M. (1997). Statistical evidence: A likelihood paradigm. London: Chapman & Hall, ISBN : 978-0412044113

Edwards, A.W.F. (1992) Likelihood, Johns Hopkins Press, ISBN : 978-0801844430
}
\keyword{confidence}
\keyword{efficacy;}
\keyword{interval}
\keyword{interval;}
\keyword{likelihood}
\keyword{likelihood-based}
\keyword{support;}
