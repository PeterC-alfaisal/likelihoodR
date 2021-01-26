#' Likelihood Support for Odds Ratio (OR)
#'
#' This function calculates the support for an OR from a 2 x 2 categorical data table.
#' An expected OR can be specified and the support calculated for this relative to the observed
#' and null (which is assumed to be 1, but can also be specified) values. A likelihood function
#' is plotted for the obtained OR with a specified likelihood interval, and expected OR,
#' if specified.
#' Chi-squared statistics are also provided and a likelihood-based % confidence interval.
#' It uses the optimize function to locate desired limits for both intervals and other
#' support calculations.
#'
#' @usage L_OR(table, null=1, exp.OR=NULL, L.int=2, alpha=0.05, cc=FALSE, toler=0.0001)
#' @param table a 2 x 2 matrix or contingency table containing counts.
#' @param null the value against which the obtained OR is tested, default = 1.
#' @param exp.OR an expected or hypothetical OR.
#' @param L.int likelihood interval given as support values, e.g. 2 or 3, default = 2.
#' @param alpha the significance level used, 1 - alpha interval calculated, default = 0.05.
#' @param cc logical indicating whether to apply continuity correction, default = FALSE.
#' @param toler the desired accuracy using optimise, default = 0.0001.
#'
#' @return
#' $S.val - support for observed OR from expected.
#'
#' $df - degrees of freedom.
#'
#' $exp.OR - expected OR.
#'
#' $S.exp.ORvsObs - support for expected OR versus observed.
#'
#' $S.exp.ORvsNull - support for expected OR versus the null.
#'
#' $L.int - likelihood interval of observed OR for specified level of support.
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
#' @keywords Likelihood; support; odds ratio; likelihood interval; confidence interval
#'
#' @export
#'
#' @importFrom graphics segments
#' @importFrom graphics lines
#' @importFrom stats optimize
#' @importFrom stats chisq.test
#' @importFrom stats qchisq
#'
#' @examples # for folic acid and neural tube defects example, p 146
#' tab <- as.table(rbind(c(6,587),c(21,581)))
#' dimnames(tab) <- list(Treatment=c("Folic acid","None"),Defect=c("Yes","No"))
#' L_OR(tab, exp.OR = 0.5, L.int = 2)
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


L_OR <- function(table, null=1, exp.OR=NULL, L.int=2, alpha=0.05, cc=FALSE, toler=0.0001) {

  SexOR_null=NULL  #NULL when exp.OR not specified
  SexOR_obs=NULL
  exa=NULL
  exa$objective=NULL

r1tot <- sum(table[1],table[3]) #sum of 1st row
r2tot <- sum(table[2],table[4]) #sum of 2nd row
c1tot <- sum(table[1:2])
c2tot <- sum(table[3:4])
grandtot <- c1tot+c2tot

a <- table[1]
b <- table[2]
c <- table[3]
d <- table[4]
(orv <- (a*d)/(b*c)) # actual odds ratio from the contingency table

res <- 100            # resolution, increase for greater resolution
arrlen <- res*c1tot-1
xs <- 0; ys <- 0

for (i in 1:arrlen) {     # arrays to plot likelihood
  dv <- i/res
  ys[i] <- exp(-sum(a*log(a/dv), b*log(b/(c1tot-dv)), c*log(c/(r1tot-dv)), d*log(d/(r2tot-c1tot+dv))))
  xs[i] <- dv*(r2tot-c1tot+dv)/((c1tot-dv)*(r1tot-dv))
}

# likelihood-based % confidence interval
goal = -qchisq(1-alpha,1)/2
f <- function(x,a,b,c,d,c1tot,r1tot,r2tot,goal) {
  (-sum(a*log(a/x), b*log(b/(c1tot-x)), c*log(c/(r1tot-x)), d*log(d/(r2tot-c1tot+x)))-goal)^2
    }
xmin1 <- optimize(f, c(1, a), tol = toler, a, b, c, d, c1tot, r1tot, r2tot, goal)
xmin2 <- optimize(f, c(a, c1tot), tol = toler, a, b, c, d, c1tot, r1tot, r2tot, goal)
beg <- xmin1$minimum*(r2tot-c1tot+xmin1$minimum)/((c1tot-xmin1$minimum)*(r1tot-xmin1$minimum))
end <- xmin2$minimum*(r2tot-c1tot+xmin2$minimum)/((c1tot-xmin2$minimum)*(r1tot-xmin2$minimum))

# same for likelihood
goalL <- -L.int
xmin1L <- optimize(f, c(1, a), tol = toler, a, b, c, d, c1tot, r1tot, r2tot, goalL)
xmin2L <- optimize(f, c(a, c1tot), tol = toler, a, b, c, d, c1tot, r1tot, r2tot, goalL)
begL <- xmin1L$minimum*(r2tot-c1tot+xmin1L$minimum)/((c1tot-xmin1L$minimum)*(r1tot-xmin1L$minimum))
endL <- xmin2L$minimum*(r2tot-c1tot+xmin2L$minimum)/((c1tot-xmin2L$minimum)*(r1tot-xmin2L$minimum))

# to determine x axis space for plot
dif <- orv-beg
lolim <- orv - 3*dif; hilim <- orv + 4*dif
if (orv < 1 ) { hilim <- orv + 6*dif}
if (lolim < 0) {lolim <- 0}

# to determine height of exp.OR on likelihood function
if (!is.null(exp.OR)) {
  goal <- exp.OR
g <- function(x,c1tot,r1tot,r2tot,goal) {
  ((x*(r2tot-c1tot+x)/((c1tot-x)*(r1tot-x)))-goal)^2
  }
exa <- optimize(g, c(1, c1tot), tol = toler, c1tot, r1tot, r2tot, goal)
xa <- unname(unlist(exa[1]))
xah <- exp(-sum(a*log(a/xa), b*log(b/(c1tot-xa)), c*log(c/(r1tot-xa)), d*log(d/(r2tot-c1tot+xa))))
}

# and likelihood for 1 (null value)

goal <- null
h <- function(x,c1tot,r1tot,r2tot,goal) {
  ((x*(r2tot-c1tot+x)/((c1tot-x)*(r1tot-x)))-goal)^2
}
exan <- optimize(h, c(1, c1tot), tol = toler, c1tot, r1tot, r2tot, goal)
xa <- unname(unlist(exan[1]))
nullh <- exp(-sum(a*log(a/xa), b*log(b/(c1tot-xa)), c*log(c/(r1tot-xa)), d*log(d/(r2tot-c1tot+xa))))

S2way <- log(1) - log(nullh) # check that this should be the same as S for observed OR

# do the plot with lines
plot(xs,ys,xlim=c(lolim,hilim),type="l", lwd = 1, xlab = "Odds Ratio", ylab = "Likelihood")
lines(c(orv,orv),c(0,1),lty=2) # add MLE as dashed line
segments(begL, exp(goalL), endL, exp(goalL), lwd = 0.2, col = "red")
lines(c(null,null),c(0,nullh), lty=1, col = "black") # add H prob as black line
if (!is.null(exp.OR)) {
  lines(c(exp.OR,exp.OR), c(0,xah), lty=1, col = "blue") # add H prob as blue line
}

# direct calculation of OR support
suppressWarnings(lt <- chisq.test(table,correct=cc)) # ignore warning message
if (null == 1) S2way <- sum(lt$observed * log(lt$observed/lt$expected))

if (!is.null(exp.OR)) {
  SexOR_null <- log(xah) - log(nullh)
  SexOR_obs <- SexOR_null - S2way
  }

cat("\nSupport for observed OR ", round(orv,4), " (dashed line) versus null of ", null,
    " (black line) = ", round(S2way,3), "\n Support for specified OR of ", exp.OR,
    " (blue line) versus observed = ", if (!is.null(exp.OR)) round(SexOR_obs,3),
    "\n Support for specified OR versus null value = ",
    if (!is.null(exp.OR)) round(SexOR_null,3),
    sep= "", "\n   S-", L.int," likelihood interval (red line) is from ",
    c(round(begL,5), " to ", round(endL,5)),
    "\n\nChi-square(1) = ", round(lt$statistic,3), ",  p = ", round(lt$p.value,4), ", N = ", grandtot,
    "\n   Likelihood-based ", 100*(1-alpha), "% confidence interval from ",
    c(round(beg,5), " to ", round(end,5)))


invisible(list(S.val = S2way, df = unname(lt$parameter), exp.OR = exp.OR, S.exp.ORvsObs = SexOR_obs,
               S.exp.ORvsNull = SexOR_null, L.int = c(begL, endL), S_int = L.int,
               observed = lt$observed, expected = lt$expected,
               chi.sq = lt$statistic, corrected = cc, p.value = lt$p.value,
               residuals = lt$residuals, alpha = alpha, conf.int = c(beg, end),
               all.err.acc = c(xmin1$objective, xmin2$objective,
                               xmin1L$objective, xmin2L$objective,
                               exa$objective, exan$objective)))
}
