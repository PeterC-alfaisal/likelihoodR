#' Likelihood Support for Relative Risk (RR)
#'
#' This function calculates the support for an RR from a 2 x 2 categorical data table.
#' An expected RR can be specified and the support calculated for this relative to the observed
#' and null (which is assumed to be 1, but can also be specified) values. A likelihood function
#' is plotted for the obtained RR with a likelihood interval, and expected RR,
#' if specified. The log likelihood plot can optionally be given instead.
#' Chi-squared statistics are also provided and a likelihood-based % confidence interval.
#' It uses the optimize function to locate desired limits for both intervals and other
#' support calculations.
#'
#' @usage L_RR(table, null=1, exp.RR=NULL, L.int=2, alpha=0.05,
#' cc=FALSE, toler=0.0001, logplot=FALSE, supplot=-10, verb=TRUE)
#' @param table a 2 x 2 matrix or contingency table containing counts.
#' @param null the value against which the obtained RR is tested, default = 1.
#' @param exp.RR an expected or hypothetical RR.
#' @param L.int likelihood interval given as support values, e.g. 2 or 3, default = 2.
#' @param alpha the significance level used, 1 - alpha interval calculated, default = 0.05.
#' @param cc logical indicating whether to apply continuity correction, default = FALSE.
#' @param toler the desired accuracy using optimise, default = 0.0001.
#' @param logplot plot vertical axis as log likelihood, default = FALSE
#' @param supplot set minimum likelihood display value in plot, default = -10
#' @param verb show output, default = TRUE.
#'
#' @return
#' $S.val - support for observed RR from expected.
#'
#' $df - degrees of freedom.
#'
#' $exp.RR - expected RR.
#'
#' $S.exp.RRvsObs - support for expected RR versus observed.
#'
#' $S.exp.RRvsNull - support for expected RR versus the null.
#'
#' $L.int - likelihood interval of observed RR for specified level of support.
#'
#' $S_int - specified likelihood interval in units of support.
#'
#' $observed - observed frequencies.
#'
#' $expected - the expected values for null hypothesis of no interaction.
#'
#' $chi.sq - chi-squared statistic.
#'
#' $corrected - whether chi-squared was corrected, default = FALSE.
#'
#' $p.value - p value.
#'
#' $residuals - the Pearson residuals.
#'
#' $alpha - specified significance level.
#'
#' $conf.int - likelihood-based confidence interval for observed RR.
#'
#' $all.err.acc - error accuracy for each application of the optimize function.
#'
#' @keywords Likelihood; support; relative risk; likelihood interval; confidence interval
#'
#' @export
#'
#' @importFrom graphics segments
#' @importFrom graphics plot
#' @importFrom graphics lines
#' @importFrom stats optimize
#' @importFrom stats chisq.test
#' @importFrom stats qchisq
#'
#' @examples # for folic acid and neural tube defects example
#' tab <- as.table(rbind(c(6,587),c(21,581)))
#' dimnames(tab) <- list(Treatment=c("Folic acid","None"),Defect=c("Yes","No"))
#' L_RR(tab, exp.RR = 0.5, L.int = 2, alpha=0.05, cc=FALSE,
#' toler=0.0001, logplot=FALSE, supplot=-10, verb=TRUE)
#'
#'# S. Korea COVID-19 patient mortality
#' tabcor <- as.table(rbind(c(41,3095),c(34,4992)))
#' dimnames(tabcor) <- list(Sex=c("Male","Female"),Status=c("Dead","Alive"))
#' L_RR(tabcor, exp.RR = 0.5, L.int = 2, alpha=0.05, cc=FALSE, toler=0.0001,
#' logplot=FALSE, supplot=-10, verb=TRUE)
#'
#' @references Aitkin, M. et al (1989) Statistical Modelling in GLIM, Clarendon Press, ISBN : 978-0198522041
#'
#' Cahusac, P.M.B. (2020) Evidence-Based Statistics, Wiley, ISBN : 978-1119549802
#'
#' Royall, R. M. (1997). Statistical evidence: A likelihood paradigm. London: Chapman & Hall, ISBN : 978-0412044113
#'
#' Edwards, A.W.F. (1992) Likelihood, Johns Hopkins Press, ISBN : 978-0801844430
#'
#' Dienes, Z. (2008) Understanding Psychology as a Science: An Introduction to Scientific and Statistical
#' Inference, Palgrave, MacMillan, ISBN : 978-0230542303
#'

L_RR <- function(table, null=1, exp.RR=NULL, L.int=2, alpha=0.05, cc=FALSE, toler=0.0001, logplot=FALSE, supplot=-10, verb=TRUE) {

  SexRR_null=NULL  #NULL when exp.RR not specified
  SexRR_obs=NULL

  exa=NULL
  exa$objective=NULL

  r1tot <- sum(table[1],table[3]) #sum of 1st row
  r2tot <- sum(table[2],table[4]) #sum of 2nd row
  c1tot <- sum(table[1:2])
  c2tot <- sum(table[3:4])
  grandtot <- c1tot+c2tot
  minmarg <- min(r1tot,r2tot,c1tot,c2tot)
  maxmarg <- max(r1tot,r2tot,c1tot,c2tot)


  a <- table[1]
  b <- table[2]
  c <- table[3]
  d <- table[4]

  rr <- (a*(b+d))/(b*(a+c)) # RR from the contingency table

  res <- 100
  mintot <- min(r1tot,r2tot,c1tot,c2tot)
  arrlen <- res*mintot-1
  xs <- 0; ys <- 0

    for (i in 1:arrlen) {
    dv <- i/res
    ys[i] <- exp(-sum(a*log(a/dv), b*log(b/(c1tot-dv)), c*log(c/(r1tot-dv)), d*log(d/(r2tot-c1tot+dv))))
    xs[i] <- dv*r2tot/((c1tot-dv)*(r1tot))
  }

  # likelihood-based % confidence interval

  arry <- numeric(maxmarg)   # finding endpoints for S and RR values
  for(x in 1:maxmarg) {
    arry[x] <- x*(r2tot-c1tot+x)/((c1tot-x)*(r1tot-x))
  }
  arry[!is.finite(arry)] <- 0
  ind <- which(arry > 0)
  aa <- split(ind, cumsum(c(0, diff(ind) > 1)))
  dvs <- min(aa$'0')-1
  dve <- max(aa$'0')+1

  goal = -qchisq(1-alpha,1)/2
  f <- function(x,a,b,c,d,c1tot,r1tot,r2tot,goal) {
    (-sum(a*log(a/x), b*log(b/(c1tot-x)), c*log(c/(r1tot-x)), d*log(d/(r2tot-c1tot+x)))-goal)^2
  }
  xmin1 <- optimize(f, c(0, a), tol = toler, a, b, c, d, c1tot, r1tot, r2tot, goal)
  xmin2 <- optimize(f, c(a, dve), tol = toler, a, b, c, d, c1tot, r1tot, r2tot, goal)
  beg <- xmin1$minimum*r2tot/((c1tot-xmin1$minimum)*r1tot)
  end <- xmin2$minimum*r2tot/((c1tot-xmin2$minimum)*r1tot)

  # same for likelihood
  goalL <- -L.int
  xmin1L <- optimize(f, c(0, a), tol = toler, a, b, c, d, c1tot, r1tot, r2tot, goalL)
  xmin2L <- optimize(f, c(a, dve), tol = toler, a, b, c, d, c1tot, r1tot, r2tot, goalL)
  begL <- xmin1L$minimum*r2tot/((c1tot-xmin1L$minimum)*r1tot)
  endL <- xmin2L$minimum*r2tot/((c1tot-xmin2L$minimum)*r1tot)

  # to determine x axis space for plot

  goalx <- supplot   # with e^-10 we get x values for when curve is down to 0.00004539
  suppressWarnings(xmin1x <- optimize(f, c(0, a), tol = toler, a, b, c, d, c1tot, r1tot, r2tot, goalx))
  suppressWarnings(xmin2x <- optimize(f, c(a, dve), tol = toler, a, b, c, d, c1tot, r1tot, r2tot, goalx))
  lolim <- xmin1x$minimum*(r2tot-c1tot+xmin1x$minimum)/((c1tot-xmin1x$minimum)*(r1tot-xmin1x$minimum))
  hilim <- xmin2x$minimum*(r2tot-c1tot+xmin2x$minimum)/((c1tot-xmin2x$minimum)*(r1tot-xmin2x$minimum))

  # to determine height of exp.RR on likelihood function
  h <- function(x,c1tot,r1tot,r2tot,goal) {
    (x*(r2tot)/((c1tot-x)*(r1tot))-goal)^2
  }
  if (!is.null(exp.RR)) {
    goal <- exp.RR
    exa <- optimize(h, c(dvs, dve), tol = toler, c1tot, r1tot, r2tot, goal)
    xa <- unname(unlist(exa[1]))
    xah <- exp(-sum(a*log(a/xa), b*log(b/(c1tot-xa)), c*log(c/(r1tot-xa)), d*log(d/(r2tot-c1tot+xa))))
  }

  # and likelihood for 1 (null value)

  goal <- null
  exan <- optimize(h, c(dvs, dve), tol = toler, c1tot, r1tot, r2tot, goal)
  xa <- unname(unlist(exan[1]))
  nullh <- exp(-sum(a*log(a/xa), b*log(b/(c1tot-xa)), c*log(c/(r1tot-xa)), d*log(d/(r2tot-c1tot+xa))))

  S2way <- -log(nullh) # check that this should be the same as S for observed RR

  # do the plot with lines

  if(isFALSE(logplot)) {
  plot(xs,ys,xlim=c(lolim,hilim),type="l", lwd = 1, xlab = "Relative Risk", ylab = "Likelihood")
  lines(c(rr,rr),c(0,1),lty=2) # add MLE as dashed line
  segments(begL, exp(goalL), endL, exp(goalL), lwd = 0.2, col = "red")
  lines(c(null,null),c(0,nullh), lty=1, col = "black") # add H prob as black line
  if (!is.null(exp.RR)) {
    lines(c(exp.RR,exp.RR), c(0,xah), lty=1, col = "blue") # add H prob as blue line
    }
  } else {  plot(xs,log(ys),xlim=c(lolim,hilim), ylim=c(supplot,0), type="l", lwd = 1, xlab = "Relative Risk", ylab = "Log Likelihood")
    lines(c(rr,rr),c(supplot,0),lty=2) # add MLE as dashed line
    segments(begL, goalL, endL, goalL, lwd = 0.2, col = "red")
    lines(c(null,null),c(supplot,log(nullh)), lty=1, col = "black") # add H prob as black line
    if (!is.null(exp.RR)) {
      lines(c(exp.RR,exp.RR), c(supplot,log(xah)), lty=1, col = "blue") # add H prob as blue line
    }
  }

  # direct calculation of RR support
  suppressWarnings(lt <- chisq.test(table,correct=cc)) # ignore warning message
  if (null == 1) S2way <- sum(lt$observed * log(lt$observed/lt$expected))

  if (!is.null(exp.RR)) {
    SexRR_null <- log(xah) - log(nullh)
    SexRR_obs <- SexRR_null - S2way
    }

  if(verb) {
    print(table)
    cat("\nSupport for observed RR ", round(rr,4), " (dashed line) versus null of ", null,
      " (black line) = ", round(S2way,3), "\n Support for specified RR of ", exp.RR,
      " (blue line) versus observed = ", if (!is.null(exp.RR)) round(SexRR_obs,3),
      "\n Support for specified RR versus null value = ",
      if (!is.null(exp.RR)) round(SexRR_null,3), sep= "", "\n   S-", L.int,
      " likelihood interval (red line) is from ", c(round(begL,5), " to ", round(endL,5)),
      "\n\nChi-square(1) = ", round(lt$statistic,3), ",  p = ", format.pval(lt$p.value,4), ", N = ",
      grandtot, "\n   Likelihood-based ", 100*(1-alpha), "% confidence interval from ",
      c(round(beg,5), " to ", round(end,5)), "\n ")
  }

  invisible(list(S.val = S2way, df = unname(lt$parameter), exp.RR = exp.RR,
                 S.exp.RRvsObs = SexRR_obs,
                 S.exp.RRvsNull = SexRR_null, L.int = c(begL, endL), S_int = L.int,
                 observed = lt$observed, expected = lt$expected,
                 chi.sq = lt$statistic, corrected = cc, p.value = lt$p.value,
                 residuals = lt$residuals, alpha = alpha, conf.int = c(beg, end),
                 all.err.acc = c(xmin1$objective, xmin2$objective,
                                 xmin1L$objective, xmin2L$objective,
                                 exa$objective, exan$objective)))
}
