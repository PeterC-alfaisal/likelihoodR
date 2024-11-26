#' Likelihood Support for Correlation
#'
#' This function calculates the support for a correlation from 2 vectors of data.
#' An expected correlation can be specified and the support calculated for this relative to the observed
#' and the null (which is assumed to be 0, but can also be specified) values. A likelihood function
#' is plotted for the obtained correlation with a likelihood interval added and expected correlation,
#' if specified. Conventional p value is also given.
#'
#' @usage L_corr(xv, yv, null=0, exp.r=NULL, L.int=2, alpha=.05,
#' toler=0.0001, logplot=FALSE, supplot=-10, verb=TRUE)
#' @param xv a numeric vector.
#' @param yv a numeric vector the same length as xv.
#' @param null the null value, default = 0.
#' @param exp.r a specified correlation (could be expected value for the study), default = NULL.
#' @param L.int likelihood interval given as support values, e.g. 2 or 3, default = 2.
#' @param alpha the significance level used, 1 - alpha interval calculated, default = 0.05.
#' @param toler the desired accuracy using optimise, default = 0.0001.
#' @param logplot plot vertical axis as log likelihood, default = FALSE
#' @param supplot set minimum likelihood display value in plot, default = -10
#' @param verb show output, default = TRUE.
#'
#' @return
#' $obs.r - observed correlation.
#'
#' $S.0 - support for observed correlation versus the null.
#'
#' $S.1 - support for the specified correlation versus observed correlation.
#'
#' $S.10 - support for the specified correlation versus the null value.
#'
#' $exp.r - the specified correlation.
#'
#' $N - the sample size.
#'
#' $p.value - the p value for significance test versus 0.
#'
#' $like.int - the likelihood interval.
#'
#' $like.int.spec - the specified likelihood interval in terms of support.
#'
#' $conf.int - the % confidence interval for the correlation.
#'
#' $alpha.spec - the specified alpha for the % confidence interval.
#'
#' @keywords Likelihood; support; correlation; likelihood interval
#'
#' @export
#'
#' @importFrom stats cor.test
#' @importFrom graphics segments
#' @importFrom graphics lines
#' @importFrom graphics curve
#'
#' @examples # for heptathlon example, p 104
#' m200 <- c(22.6,	23.7,	23.1,	23.6,	23.6,	23.6,	25.5,
#' 23.9,	24.5,	23.9,	24.9,	24.8,	24.7,
#' 25.0,	24.6,	24.9,	25.0,	25.6,	24.8,
#' 25.5,	25.7,	24.9,	26.6,	25.2,	26.2)
#' m800 <- c(128.5,	126.1,	124.2,	132.5,
#' 134.7,	132.5,	138.5,	127.9,	133.7,	132.2,
#' 136.1,	142.8,	125.8, 131.5,	137.1,	134.9,
#' 146.7,	133.9,	146.4,	144.0,	133.4,
#' 138.0,	139.2,	137.3,	163.4)
#' m=L_corr(m200, m800, null=0, exp.r=.4, L.int=3, alpha=.05,
#' toler=0.0001, logplot=FALSE, supplot=-10, verb=TRUE)
#' m
#' #Note: the support for observed vs 0 is different from book (5.776 vs 5.700)
#' #due to differences in calculation of r by Excel and R
#'
#' @references
#' Cahusac, P.M.B. (2020) Evidence-Based Statistics, Wiley, ISBN : 978-1119549802
#'
#' Royall, R. M. (1997). Statistical evidence: A likelihood paradigm. London: Chapman & Hall, ISBN : 978-0412044113
#'
#' Edwards, A.W.F. (1992) Likelihood, Johns Hopkins Press, ISBN : 978-0801844430
#'

L_corr <- function(xv, yv, null=0, exp.r=NULL, L.int=2, alpha=.05, toler=0.0001, logplot=FALSE, supplot=-10, verb=TRUE) {

  m <- cor.test(xv, yv, conf.level =1-alpha)
  r <- unname(m$estimate)
  N <- unname(m$parameter) + 2

# using z
  rtrans <- function(r_value) (0.5 * log(abs((1 + r_value)/(1 - r_value))))
  se <- 1/(sqrt(N - 3))
  z0 <- (null - rtrans(r))/se

# support for observed versus null
  S0 <- -z0^2/2
  SM <- -S0

  z1 <- (rtrans(exp.r)-rtrans(r))/se
  S1 <- -z1^2/2             # support for exp.r versus observed r

# support for specified r versus null
  S_10 <- S1 - S0

# likelihood-based % confidence interval
  x=0

  f <- function(x, r, goal) {
    (-((rtrans(x)-rtrans(r))/se)^2/2-goal)^2
  }

  goal = -qchisq(1-alpha,1)/2
  xmin1 <- optimize(f, c(-1, r), tol = toler, r, goal)
  xmin2 <- optimize(f, c(r, 1), tol = toler, r, goal)
  beg <- xmin1$minimum
  end <- xmin2$minimum

# to determine x axis space for plot

  goalx = supplot
  xmin1x <- optimize(f, c(-1, r), tol = toler, r, goalx)
  xmin2x <- optimize(f, c(r, 1), tol = toler, r, goalx)
  x_lim <- c(xmin1x$minimum, xmin2x$minimum)

  r_dash_lower <- rtrans(r)-se*sqrt(L.int*2)
  r_dash_upper <- rtrans(r)+se*sqrt(L.int*2)

  lolim <- (exp(2*r_dash_lower)-1)/(exp(2*r_dash_lower)+1)
  hilim <- (exp(2*r_dash_upper)-1)/(exp(2*r_dash_upper)+1)

  if(isFALSE(logplot)) {
  curve(exp(-((rtrans(r)-rtrans(x))/se)^2/2), xlim = x_lim, xlab = "r", ylab = "Likelihood")
  lines(c(r,r),c(0,1),lty=2) # add MLE as dashed line
  lines(c(null,null),c(0,exp(-((rtrans(r)-rtrans(null))/se)^2/2)), lty=1, col = "black") # add H prob as black line
  if (!is.null(exp.r)) {
    lines(c(exp.r,exp.r), c(0,exp(-((rtrans(r)-rtrans(exp.r))/se)^2/2)), lty=1, col = "blue") # add H prob as blue line
  }
  segments(lolim, exp(-L.int), hilim, exp(-L.int), lwd = 0.2, col = "red")
  } else {
    curve(-((rtrans(r)-rtrans(x))/se)^2/2, xlim = x_lim, xlab = "r", ylab = "Log Likelihood")
    lines(c(r,r),c(supplot,0),lty=2) # add MLE as dashed line
    lines(c(null,null),c(supplot,(-((rtrans(r)-rtrans(null))/se)^2/2)), lty=1, col = "black") # add H prob as black line
    if (!is.null(exp.r)) {
      lines(c(exp.r,exp.r), c(supplot,(-((rtrans(r)-rtrans(exp.r))/se)^2/2)), lty=1, col = "blue") # add H prob as blue line
    }
    segments(lolim, -L.int, hilim, -L.int, lwd = 0.2, col = "red")
  }

  if(verb) cat("\nSupport for observed correlation ", round(r,4), " (dashed line) versus null of ", null,
     " (black line) = ", round(SM,3), sep= "",
      "\n Support for specified correlation of ", exp.r,
      " (blue line) versus observed r = ", if (!is.null(exp.r)) round(S1,3),
     "\n Support for specified correlation versus null = ",
      if (!is.null(exp.r)) round(S_10,3),
     "\n S-", L.int," likelihood interval (red line) is from ",round(lolim,5), " to ", round(hilim,5),
     "\n\n P value = ", format.pval(m$p.value,4),
     "\n Likelihood-based ", 100*(1-alpha), "% confidence interval from ",
     c(round(beg,5), " to ", round(end,5)), " \n N = ", N, "\n ")

  invisible(list(obs.r = r, S.0 = SM, S.1 = S1, S.10 = S_10, exp.r = exp.r, N = N, p.value = m$p.value,
               like.int = c(lolim,hilim), like.int.spec = L.int, conf.int = m$conf.int[1:2],
               conf.int = c(beg, end),
               alpha.spec = alpha))
}
