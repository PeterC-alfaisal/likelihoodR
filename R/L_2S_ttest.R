#' Likelihood Supports for Independent Samples t Test
#'
#' This function calculates several different supports for independent samples. Effect
#' size (Cohen's d) and a second alternative hypothesis value can be specified.
#' The maximum support is the support for the observed mean versus the null value.
#' The support for the specified d versus the null is also calculated. If a second
#' hypothesis value is specified (in units of the original measurements) then two
#' further supports are calculated: d versus 2nd alternative hypothesis, and 2nd
#' alternative hypothesis versus the null. The likelihood curve graphic with MLE and
#' specified hypothesis values is produced. Finally, the requested likelihood interval
#' is provided. The t, p and observed d values for the test against the null are given.
#' If variances are specified as unequal then uses Welch's test where
#' homogeneity of variance is not required.
#'
#' @usage L_2S_ttest(data, group, veq=0, null=0, d=0.5, alt.2=NULL, L.int=2, verb=TRUE)
#' @param data a (non-empty) numeric vector of data values.
#' @param group an integer vector the same length as data, coding for 2 groups.
#' @param veq whether variances are equal: 1 = Yes, 0 = No, default = 0.
#' @param null value for the null hypothesis, default = 0.
#' @param d Cohen's effect size, default = 0.5.
#' @param alt.2 value for an alternative hypothesis, in units used for data, default = NULL.
#' @param L.int likelihood interval given as support values, e.g. 2 or 3, default = 2.
#' @param verb show output, default = TRUE.
#'
#'
#' @return
#' $obs.diff - the observed difference in means.
#'
#' $df - degrees of freedom.
#'
#' $var.eq - if not equal (0) then Welch's test used.
#'
#' $alt.H1 - mean value according to specified d.
#'
#' $alt.H2 - specified second hypothesis value.
#'
#' $S_max - maximum support for observed mean difference against the null.
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
#' $null.value - null value
#'
#' $t.val - t value for test against null.
#'
#' $p.val - p value for test against null.
#'
#' $d.obs - observed effect size (from null).
#'
#'
#' @keywords Likelihood; support; independent samples t test; likelihood interval
#'
#' @export
#'
#' @importFrom graphics curve
#' @importFrom graphics segments
#' @importFrom graphics lines
#' @importFrom stats complete.cases
#' @importFrom stats t.test
#' @importFrom stats qt
#' @importFrom stats sd
#'
#' @examples # using a variation on Gosset's original additional hours of sleep data, p 59
#' mysample <- c(0.7, -1.6, -0.2, -1.2, -0.1, 3.4, 3.7, 0.8, 0.0, 2.0)
#' treat <- rep(1:0,each=5)
#' L_2S_ttest(mysample, treat, veq=0, null=0, d=0.5, alt.2=2, L.int=2)
#'
#' @references Cahusac, P.M.B. (2020) Evidence-Based Statistics, Wiley, ISBN : 978-1119549802
#'
#' Baguley, T. (2012) Serious Stats, Palgrave Macmillan, ISBN: 978-0230577183
#'
#' Edwards, A.W.F. (1992) Likelihood, Johns Hopkins Press, ISBN : 978-0801844430
#'
#' Royall, R. M. (1997) Statistical Evidence: A Likelihood Paradigm, Chapman & Hall, ISBN : 978-0412044113
#'
#' Royall, R. M. (2000). On the probability of observing misleading statistical evidence.
#' Journal of the American Statistical Association, 95, 760.


L_2S_ttest <- function(data, group, veq=0, null=0, d=0.5, alt.2=NULL, L.int=2, verb=TRUE) {
  if (is.null(d)) d = 0
  dat <- data.frame(data, group)
  ad <- dat[complete.cases(dat), ] # remove missing, NA or NaN, case-wise
  dat <- ad$data
  gp <- ad$group
  lev <- levels(factor(gp))
  ns <- length(sort(dat[gp==lev[1]])) #sort to remove NA
  nc <- length(sort(dat[gp==lev[2]]))
  sd1 <- sd(dat[gp==lev[1]],na.rm=TRUE)
  sd2 <- sd(dat[gp==lev[2]],na.rm=TRUE)
  SD <- sqrt(sd1^2+sd2^2)
  tres0 <- t.test(dat~gp, mu = null, paired=FALSE,var.equal=veq) #
  m.obs <- unname(tres0$estimate[1]-tres0$estimate[2]) # t test does it this way
  sed <- unname(tres0$stderr)   #SE
  SD <- sqrt(sd1^2 + sd2^2)     # non-pooled SD
  if (veq == 1) {
    pool_var <- ((ns - 1)*sd1^2 + (nc - 1)*sd2^2)/(ns + nc - 2)
    SD <- sqrt(pool_var)   # pooled SD for equal variance assumption
  }
  obs_es <- m.obs/SD      # observed effect size
  alt.1 <- SD*d           # specified effect size
  tres1 <- t.test(dat~gp, mu = alt.1, paired=FALSE, var.equal=veq) # sp. effect size
  N <- ns + nc
  df <- unname(tres0$parameter)
  pval <- unname(unlist(tres0$p.value))
  tval <- unname(tres0$statistic)
  like0 <- (1 + tval^2/df)^-(N/2) #L0
  like1 <- unname((1 + tres1$statistic^2/df)^-(N/2)) #L1

  # Supports
  S_10 <- log(like1)-log(like0)
  S_12 = NULL
  S_20 = NULL
  if (!is.null(alt.2 )) {
    tres2 <- t.test(dat~gp, mu = alt.2, paired=FALSE, var.equal=veq)  #alt H2
    like2 <- unname((1 + tres2$statistic^2/df)^-(N/2)) #L2
    S_12 <- log(like1)-log(like2)
    S_20 <- S_10 - S_12
  }

  # Maximum likelihood ratio and S
  S_m <- log(1) - log(like0)

  # Plot the likelihood function
  x=0
  curve((1 + ((m.obs-x)/sed)^2/df)^-(N/2),
        xlim = c(m.obs-5*sed,m.obs+5*sed), ylab = "Likelihood",
        xlab = "Mean difference")
  lines(c(m.obs,m.obs),c(0,1),lty=2) # add mean as dashed line
  lines(c(null, null), c(0,(1 + ((m.obs-null)/sed)^2/df)^-(N/2)),lty=1, col = "black") # for null
  lines(c(alt.1,alt.1),c(0,(1 + ((m.obs-alt.1)/sed)^2/df)^-(N/2)),lty=1, col = "blue") # d
  if (!is.null(alt.2)) lines(c(alt.2,alt.2),c(0,(1 + ((m.obs-alt.2)/sed)^2/df)^-(N/2)),lty=1, col = "green") # alt.2 value

  # Add Likelihood interval
  lolim <- m.obs - sed*sqrt((exp(L.int*2/N)-1)*df)
  hilim <- m.obs + sed*sqrt((exp(L.int*2/N)-1)*df)
  segments(lolim, exp(-L.int), hilim, exp(-L.int), lwd = 0.2, col = "red")

  if(verb) cat("\nMaximum support for the observed mean ", m.obs, " (dashed line) against the null ",
      null, " (black line) = ", round(S_m,3), sep= "", "\n Support for d of ", d,
      " (", round(alt.1,3), ", blue line) versus null = ", round(S_10,3),
      "\n Support for d versus 2nd alt Hypothesis ", alt.2, " (green line) = ", if (!is.null(alt.2 )) round(S_12,3),
      "\n Support for 2nd alt Hypothesis versus null = ", if (!is.null(alt.2 )) round(S_20,3),
      "\n\n S-", L.int, " likelihood interval (red line) is from ", c(round(lolim,5),
                                                                      " to ", round(hilim,5)),
      "\n\nt(", round(df,1), ") = ", round(tval,3), ", p = ",
      pval, ", d = ", round(obs_es,3), "\n ")


  invisible(list(obs.diff = m.obs, df = df, var.eq = veq, alt.H1 = alt.1, alt.H2 = alt.2, S_max = S_m,
                 S_10 = S_10, S_12 = S_12, S_20 = S_20, like.int = c(lolim, hilim),
                 L.int.spec = L.int, null.value = null, t.val = tval, p.val = pval, d.obs = obs_es))
}
