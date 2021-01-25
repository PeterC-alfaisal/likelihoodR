#' Likelihood-based Confidence Interval for the Binomial
#'
#' This function allows you to calculate the usual confidence interval, but is likelihood-based.
#' It uses the optimize function to locate desired limits and their error.
#' @usage binpL(a, n, alpha, toler)
#' @param a the number of successes.
#' @param n total number of trials.
#' @param alpha the significance level used, 1 - alpha interval calculated, default = 0.05.
#' @param toler the desired accuracy using optimise, default = 0.0001.
#'
#' @return
#' $prob.val - MLE for probability from the data.
#'
#' $succ.fail - number of successes and failures.
#'
#' $conf.int - likelihood-based confidence interval.
#'
#' $err.acc - error accuracy for optimize function.
#'
#' @keywords Likelihood-based; binomial, confidence interval
#'
#' @export
#'
#' @importFrom graphics curve
#' @importFrom graphics lines
#' @importFrom graphics segments
#'
#' @examples binpL(86,94,0.05)
#'
#' @references Aitkin, M. et al (1989) Statistical Modelling in GLIM, Clarendon Press, ISBN : 978-0198522041
#'
#' Cahusac, P.M.B. (2020) Evidence-Based Statistics, Wiley, ISBN : 978-1119549802
#'
#' Royall, R. M. (1997). Statistical evidence: A likelihood paradigm. London: Chapman & Hall, ISBN : 978-0412044113
#'
#' Edwards, A.W.F. (1992) Likelihood, Johns Hopkins Press, ISBN : 978-0801844430


binpL <- function(a,n,alpha=0.05, toler=0.0001){
  r = n - a; p = a/n;x = 0; goal = -qchisq(1-alpha,1)/2
  f <- function(x,a,r,p,goal) (a*log(x)+r*log(1-x)-(a*log(p)+r*log(1-p))-goal)^2
  xmin1 <- optimize(f, c(0, p), tol = toler, a, r, p, goal)
  xmin2 <- optimize(f, c(p, 1), tol = toler, a, r, p, goal)
  if (p < .5) { lolim <- 0; hilim <- p + 6*sqrt(p*(1-p)/n)
  }  else {hilim <- 1; lolim <- p - 6*sqrt(p*(1-p)/n)}
  if (lolim < 0) {lolim <- 0}
  if (hilim > 1) {hilim <- 1}
  curve((x^a*(1-x)^r)/(p^a*(1-p)^r), xlim = c(lolim,hilim), xlab = "Probability", ylab = "Likelihood")
  lines(c(p,p),c(0,1),lty=2) # add MLE as dashed line
  segments(xmin1$minimum, exp(goal), xmin2$minimum, exp(goal), lwd = 0.2, col = "red")
  cat("Likelihood-based ", 100*(1-alpha), sep= "", "% confidence interval: ",
      c(round(xmin1$minimum,5), "  ", round(xmin2$minimum,5)))
  invisible(list(prob.val = p, succ.fail = c(a,r), conf.int = c(xmin1$minimum, xmin2$minimum),
       err.acc = c(xmin1$objective, xmin2$objective)))
}

