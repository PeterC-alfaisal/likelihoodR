#' Likelihood Support for Efficacy
#'
#' This function calculates the support for the efficacy, the likelihood interval
#' and the likelihood-based confidence interval.
#' It uses the optimize function to locate desired limits and their error.
#'
#' @usage L_efficacy(a, n, null=0, exp.eff=NULL, L.int=2, alpha=0.05, toler=0.0001)
#' @param a the number of affected in control group.
#' @param n total number of participants.
#' @param null the null value for efficacy, if no effect then it would be 0, default = 0.
#' @param exp.eff the expected or hypothesized efficacy, default = NULL.
#' @param L.int likelihood interval given as support values, e.g. 2 or 3, default = 2.
#' @param alpha the significance level used, 1 - alpha interval calculated, default = 0.05.
#' @param toler the desired accuracy using optimise, default = 0.0001.
#'
#' @return
#' $S.val - support for the observed efficacy versus the null value.
#'
#' $obs.eff - the observed efficacy.
#'
#' $null - the null efficacy.
#'
#' $exp.eff - expected efficacy as specified.
#'
#' $S.exp.vsObs - support for expected efficacy versus observed.
#'
#' $S.exp.versus.null - support for the expected efficacy versus the null.
#'
#' $L.int - the likelihood interval for the observed efficacy.
#'
#' $S_int - the specified likelihood interval.
#'
#' $observed - observed numbers affected in control and intervention groups.
#'
#' $expected - expected numbers according to the null.
#'
#' $chi.sq - chi-squared statistic.
#'
#' $p.value - p value associated with chi-squared statistic.
#'
#' $df - degrees of freedom for chi-squared.
#'
#' $residuals - the Pearson residuals.
#'
#' $conf.int - likelihood-based confidence interval according to specified alpha.
#'
#' $alpha - specified alpha for confidence interval.
#'
#' $all.err.acc - error accuracy for each application of the optimize function.
#'
#' @keywords efficacy; support; likelihood interval; likelihood-based confidence interval
#' @export
#'
#' @importFrom stats optimize
#' @importFrom stats chisq.test
#' @importFrom graphics lines
#' @importFrom graphics segments
#'
#' @examples # pfizer covid-19 efficacy 2020
#' m = L_efficacy(a = 86, n = 94, null=0.8, exp.eff=0.95, L.int=2, alpha=0.05, toler=0.0001)
#' m
#'
#' @references Aitkin, M. et al (1989) Statistical Modelling in GLIM, Clarendon Press, ISBN : 978-0198522041
#'
#' Cahusac, P.M.B. (2020) Evidence-Based Statistics, Wiley, ISBN : 978-1119549802
#'
#' Royall, R. M. (1997). Statistical evidence: A likelihood paradigm. London: Chapman & Hall, ISBN : 978-0412044113
#'
#' Edwards, A.W.F. (1992) Likelihood, Johns Hopkins Press, ISBN : 978-0801844430


L_efficacy <- function(a, n, null=0, exp.eff=NULL, L.int=2, alpha=0.05, toler=0.0001){

  r = n - a; eff = (a - r)/a; goal = -qchisq(1-alpha,1)/2

  n_a <- n - n/(2 - null)  # a for null
  Sexeff_null=-(a*log(a)+r*log(r)-(a*log(n-n_a)+r*log(n_a)))  #NULL when exp.OR not specified
  n_a_h <- exp(Sexeff_null)
  S.eff <- log(1) - Sexeff_null  # support for observed versus null

  S.exp_null = NULL       # support for expected efficacy versus null
  Sexeff_obs = NULL

  if (!is.null(exp.eff)) {
    x_a <- n - n/(2 - exp.eff) # a for exp efficacy
    Sexef <- -(a*log(a)+r*log(r)-(a*log(n-x_a)+r*log(x_a)))
    x_a_h <- exp(Sexef)
    S.exp_null <- Sexef - Sexeff_null
    Sexeff_obs <- S.exp_null - S.eff
  }

  # direct calculation of efficacy support
  obs <- c(a,n-a); exp.p <- c((n-n_a)/n,n_a/n)

  suppressWarnings(lt <- chisq.test(obs,p=exp.p)) # ignore warning message

  f <- function(x,a,r,goal) (-(a*log(a)+r*log(r)-(a*log(n-x)+r*log(x)))-goal)^2
  xmin1 <- optimize(f, c(-1, r), tol = toler, a, r, goal)
  xmin2 <- optimize(f, c(r, n), tol = toler, a, r, goal)
  beg <- (n - 2*xmin2$minimum)/(n-xmin2$minimum)
  end <- (n - 2*xmin1$minimum)/(n-xmin1$minimum)

  # same for likelihood
  goalL <- -L.int
  xmin1L <- optimize(f, c(-1, r), tol = toler, a, r, goalL)
  xmin2L <- optimize(f, c(r, n), tol = toler, a, r, goalL)
  begL <- (n - 2*xmin2L$minimum)/(n-xmin2L$minimum)
  endL <- (n - 2*xmin1L$minimum)/(n-xmin1L$minimum)

  if (eff < 0) { lolim <- -1; hilim <- 0}
  if (eff < .5) { lolim <- 0; hilim <- eff + 6*sqrt(eff*(1-eff)/n)
  } else {hilim <- 1; lolim <- eff - 6*sqrt(eff*(1-eff)/n)}
  if (lolim < -1) {lolim <- -1}
  if (hilim > 1) {hilim <- 1}

  res <- 100
  arrlen <- res*n-1
  xs <- 0; ys <- 0

  for (i in 1:arrlen) {
    dv <- i/res
    ys[i] <- exp(-(a*log(a)+r*log(r)-(a*log(n-dv)+r*log(dv))))
    xs[i] <- (n-2*dv)/(n-dv)
  }

  plot(xs, ys, xlim = c(lolim,hilim), type="l", lwd = 1,
        xlab = "Efficacy", ylab = "Likelihood")
  lines(c(eff,eff),c(0,1),lty=2) # add MLE as dashed line
  segments(begL, exp(goalL), endL, exp(goalL), lwd = 0.2, col = "red")
  lines(c(null,null), c(0,n_a_h), lty=1, col = "black") # add black line for null
  if (!is.null(exp.eff)) {
    lines(c(exp.eff,exp.eff), c(0,x_a_h), lty=1, col = "blue") # add H prob as blue line
  }

  cat("\nSupport for observed efficacy ", round(eff,4), " (dashed line) versus null of ", null,
      " (black line) = ", round(S.eff,3),"\n Support for specified efficacy of ",
      exp.eff, " (blue line) versus observed = ", if (!is.null(exp.eff)) round(Sexeff_obs,3),
      "\n Support for specified efficacy versus null value = ",
      if (!is.null(exp.eff)) round(S.exp_null,3), sep= "", "\n   S-", L.int,
      " likelihood interval (red line) is from ", c(round(begL,5), " to ", round(endL,5)),
      "\n\nChi-square(1) = ", round(lt$statistic,3), ",  p = ", round(lt$p.value,4), ", N = ",
      n, "\n   Likelihood-based ", 100*(1-alpha), "% confidence interval from ",
      c(round(beg,5), " to ", round(end,5)))


  invisible(list(S.val = S.eff,  obs.eff = eff, null = null, exp.eff = exp.eff, S.exp.vsObs = Sexeff_obs,
                 S.exp.versus.null = S.exp_null, L.int = c(begL, endL), S_int = L.int,
                 observed = lt$observed, expected = lt$expected,
                 chi.sq = unname(lt$statistic), p.value = lt$p.value, df = unname(lt$parameter),
                 residuals = lt$residuals, conf.int = c(beg, end), alpha = alpha,
                 all.err.acc = c(xmin1$objective, xmin2$objective,
                                 xmin1L$objective, xmin2L$objective)))

}




