#' Likelihood Supports for the One Sample and Related Samples t Test
#'
#' This function calculates several different supports. Effect size (Cohen's d) and a
#' second alternative hypothesis value can be specified. The maximum support is the support
#' for the observed mean versus the null value. The support for the specified d versus
#' the null is also calculated. If a second hypothesis value is specified (in units of
#' the original measurements) then two further supports are calculated: d versus 2nd
#' alternative hypothesis, and 2nd alternative hypothesis versus the null.
#' The likelihood curve graphic with MLE and specified hypothesis values is produced.
#' The requested likelihood interval is provided and displayed on likelihood curve.
#' The t and p values for the test against the null value are given.
#'
#' @usage L_ttest(data1, data2, null=0, d=0.5, alt.2=NULL, L.int=2)
#' @param data1 a (non-empty) numeric vector of data values.
#' @param data2 a (non-empty) numeric vector of data values for related sample, default = NULL.
#' @param null value for the null hypothesis, default = 0.
#' @param d Cohen's effect size, default = 0.5.
#' @param alt.2 value for an alternative hypothesis, in units used for data, default = NULL.
#' @param L.int likelihood interval given for a given support value, e.g. 2 or 3, default = 2.
#'
#' @return
#' $obs.mean - the observed mean or difference in mean for related samples.
#'
#' $df - degrees of freedom.
#'
#' $alt.H1 - mean value according to specified d.
#'
#' $alt.H2 - specified second hypothesis value.
#'
#' $S_max - maximum support for observed mean against the null.
#'
#' $S_10 - support for d versus null.
#'
#' $S_12 - support for d versus specified second hypothesis.
#'
#' $S_20 - support for second hypothesis versus the null.
#'
#' $like.int - likelihood interval.
#'
#' $L.int.spec - specified likelihood interval in units of support.
#'
#' $null.value - null value.
#'
#' $t.val - t value for test against null.
#'
#' $p.val - p value for test against null.
#'
#'
#' @keywords Likelihood; support; t test; likelihood interval
#'
#' @export
#'
#' @importFrom stats complete.cases
#' @importFrom stats t.test
#' @importFrom graphics curve
#' @importFrom graphics lines
#' @importFrom graphics segments
#'
#' @examples # one sample Gosset's original additional hours of sleep data, p 29
#' mysample <- c(0.7, -1.6, -0.2, -1.2, -0.1, 3.4, 3.7, 0.8, 0.0, 2.0)
#' L_ttest(mysample, d=.5, alt.2=2, L.int=2)
#'
#' # related samples, p 56
#' mysample2 <- c(1.9, 0.8, 1.1, 0.1, -0.1, 4.4, 5.5, 1.6, 4.6, 3.4)
#' L_ttest(mysample, mysample2, d=1, alt.2=2, L.int=2)
#' @references Cahusac, P.M.B. (2020) Evidence-Based Statistics, Wiley, ISBN : 978-1119549802
#'
#' Baguley, T. (2012) Serious Stats, Palgrave Macmillan, ISBN: 978-0230577183
#'
#' Edwards, A.W.F. (1992) Likelihood, Johns Hopkins Press, ISBN : 978-0801844430
#'
#' Royall, R. M. (2000). On the probability of observing misleading statistical evidence. Journal of the American Statistical Association, 95, 760–768.


L_ttest <- function(data1, data2=NULL, null=0, d=0.5, alt.2=NULL, L.int=2) {
  x=0
  if (is.null(d)) d = 0
  adata <- data1
  if (!is.null(data2)) {
    dat <- data.frame(data1,data2)
    ad <- dat[complete.cases(dat), ] # remove missing, NA or NaN, case-wise
    adata <- ad$data2-ad$data1
  }
  t0 <- t.test(adata, mu = null) # t test
  tval <- unname(unlist(abs(t0$statistic))) # t value, abs & remove label
  m.obs <- unname(unlist(t0$estimate)) # observed mean
  df <- unname(unlist(t0$parameter)) # df
  pval <- unname(unlist(t0$p.value))
  N <- df+1
  sed <- t0$stderr # standard error
  alt.1 <- d*sed*sqrt(N)
  t0val <- (null - m.obs)/sed # t value for H0
  t1val <- (alt.1 - m.obs)/sed # t value for H1
  like0 <- (1 + t0val^2/df)^-(N/2) # L0
  like1 <- (1 + t1val^2/df)^-(N/2) # L1
# Supports
  S_10 <- log(like1)-log(like0)
  S_12 = NULL
  S_20 = NULL
  if (!is.null(alt.2 )) {
    t2val <- (alt.2 - m.obs)/sed
    like2 <- (1 + t2val^2/df)^-(N/2) # L2
    S_12 <- log(like1)-log(like2)
    S_20 <- S_10 - S_12
  }
# Maximum Support
  S_m <- log(1) - log(like0)

# Plot the likelihood function
  curve((1 + ((m.obs-x)/sed)^2/df)^-(N/2),
      xlim = c(m.obs-5*sed,m.obs+5*sed), ylab = "Likelihood",
      xlab = "Mean value")
  lines(c(m.obs,m.obs),c(0,1),lty=2) # add mean as dashed line
  lines(c(null, null), c(0,(1 + ((m.obs-null)/sed)^2/df)^-(N/2)),lty=1, col = "black") # for null
  lines(c(alt.1,alt.1),c(0,(1 + ((m.obs-alt.1)/sed)^2/df)^-(N/2)),lty=1, col = "blue") # d
  if (!is.null(alt.2 )) lines(c(alt.2,alt.2),c(0,(1 + ((m.obs-alt.2)/sed)^2/df)^-(N/2)),lty=1, col = "green") # alt.2 value
# Add Likelihood interval
  lolim <- m.obs - sed*sqrt((exp(L.int*2/N)-1)*df)
  hilim <- m.obs + sed*sqrt((exp(L.int*2/N)-1)*df)
  segments(lolim, exp(-L.int), hilim, exp(-L.int), lwd = 0.2, col = "red")

  cat("\nMaximum support for the observed mean ", m.obs, " (dashed line) against the null ",
      null, " (black line) = ", round(S_m,3), sep= "", "\n Support for d of ", d,
      " (", alt.1, ", blue line) versus null = ", round(S_10,3),
      "\n Support for d versus 2nd alt Hypothesis ", alt.2, " (green line) = ", if (!is.null(alt.2 )) round(S_12,3),
      "\n Support for 2nd alt Hypothesis versus null = ", if (!is.null(alt.2 )) round(S_20,3),
      "\n\n S-", L.int, " likelihood interval (red line) is from ", c(round(lolim,5),
                                                                      " to ", round(hilim,5)),
      "\n\nt(", df, ") = ", round(tval,3), ", p = ", pval, "\n ")

  invisible(list(obs.diff = m.obs, df = df, alt.H1 = alt.1, alt.H2 = alt.2, S_max = S_m,
                 S_10 = S_10, S_12 = S_12, S_20 = S_20, like.int = c(lolim, hilim),
                 L.int.spec = L.int, null.value = null, t.val = tval, p.val = pval))
}
