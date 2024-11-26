#' Likelihood Supports for One-way Repeated Measures ANOVA
#'
#' This function calculates support for the treatment means versus the null model,
#' Type of correction for the parameters in the null versus the means model can be selected.
#' support for the first contrast versus the group means model and the
#' support for the first versus the second contrast.
#' Also gives the F, p and partial eta-squared values for the overall analysis,
#' and F and p for the first contrast. Assumes sphericity.
#'
#'
#' @usage L_1way_RM_ANOVA(dat, group, ID, correct=1, contrast1=NULL, contrast2=NULL, verb=TRUE)
#'
#' @param dat a (non-empty) numeric vector of data values.
#' @param group a vector the same length as data, coding for k groups.
#' @param ID is an identifier for each case.
#' @param correct specifies the correction: 0 = none, 1 = Occam's bonus (default), 2 = AIC
#' @param contrast1 first contrast, default = NULL.
#' @param contrast2 second contrast, default = NULL.
#' @param verb show output, default = TRUE.
#'
#' @return
#' $S.12 - support for groups means hypothesis versus null, uncorrected
#'
#' $S.12c - corrected
#'
#' $capplied - type of correction
#'
#' $S.1m = support for first contrast versus means model.
#'
#' $S.cont.12 - support for first versus second contrast.
#'
#' $contrast1 - first contrast.
#'
#' $contrast2 - second contrast.
#'
#' $gp.means - group means.
#'
#' $df - degrees of freedom for ANOVA.
#'
#' $F.val - F value for overall analysis.
#'
#' $P.val - p value for the overall analysis.
#'
#' Fval.c1 - F value for the first contrast.
#'
#' $df.1 - degrees of freedom for the contrast.
#'
#' $P.val1 - p value for the contrast.
#'
#'
#' @keywords Likelihood; support; repeated measures ANOVA
#'
#' @export
#'
#' @importFrom stats anova
#' @importFrom graphics plot
#' @importFrom stats lm
#' @importFrom stats aggregate
#' @importFrom stats pf
#' @importFrom stats logLik
#' @importFrom lme4 lmer
#'
#' @examples # sleep data example, p 97
#' dat <- c(0.7,	-1.6,	-0.2,	-1.2,	-0.1,	3.4,	3.7,	0.8,	0,	2,
#' 1.9,	0.8,	1.1,	0.1,	-0.1,	4.4,	5.5,	1.6,	4.6,	3.4,
#' 1.5,	1.4,	0.0,	-0.7,	0.5,	5.1,	5.7,	1.5,	4.7,	3.5)
#' treat <- gl(3,10,30)
#' patients <- gl(10,1,30)
#' contrast1 <- c(-1, 0, 1)  # linear
#' contrast2 <- c(-2, 1, 1)  # 1st vs treatments 2 & 3
#'
#' m=L_1way_RM_ANOVA(dat, treat, patients, correct=1, contrast1, contrast2)
#' m
#'
#' @references Cahusac, P.M.B. (2020) Evidence-Based Statistics, Wiley, ISBN : 978-1119549802
#'
#' Dixon P. The effective number of parameters in post hoc models. Behavior
#' Research Methods. 2013; 45(3):604.
#'
#' Dixon P. The p-value fallacy and how to avoid it. Canadian Journal of
#' Experimental Psychology/Revue canadienne de psychologie expérimentale. 2003;
#' 57(3):189.
#'
#' Glover S, Dixon P. Likelihood ratios: a simple and flexible statistic for empirical
#' psychologists. Psychonomic Bulletin and Review. 2004; 11(5):791.
#'
#'

L_1way_RM_ANOVA <- function(dat, group, ID, correct=1, contrast1=NULL, contrast2=NULL, verb=TRUE) {

  m1=anova(lm(dat~ID + group))

  cp <- data.frame(dat,group,ID)

  m1.lmer<-lme4::lmer(dat~group+(1|ID),cp,REML=F)
  m0.lmer<-lme4::lmer(dat~1+(1|ID),cp,REML=F)
  np <- attr(logLik(m0.lmer),"df")  # Null: parameter each for variance, grand mean and ID
  mp <- attr(logLik(m1.lmer),"df")

  within_ss <- sum(m1$`Sum Sq`[2:3])
  eta_sq <- m1$`Sum Sq`[2]/within_ss
  dfv <- m1$Df
  dff <- dfv[2]
  N <- dfv[1]+1
  S_12 <- -0.5 * N * dff * (log(m1$`Sum Sq`[3]) - log(within_ss))

  # Correction
  Ac <- function(c,k1,k2,N) {
    if(c==0) { 0
    } else if(c==1) { 0.5*(k2-k1)
    } else if(c==2) { 1*(k2-k1)
    }
  }
  S_12c <- S_12 - Ac(correct,np,mp,N)

# contrasts
  k <- m1$Df[2]+1
  nulfg <- 0
  if (k > 2) {
  if (is.null(contrast1)) {
    conta <- contr.poly(k, scores = 1:k)
    contrast1 <- conta[,1]               # linear
    contrast2 <- conta[,2]               # quadratic
    nulfg <- NULL
  }}
  gp_means <- tapply(dat, group, mean)


  SS_1 <- N * sum(contrast1*gp_means)^2/(sum(contrast1^2))
  SS_2 <- N * sum(contrast2*gp_means)^2/(sum(contrast2^2))

  r_SS_1 <- within_ss - SS_1
  r_SS_2 <- within_ss - SS_2

  S_cont_12 <- -0.5*N*(log(r_SS_1) - log(r_SS_2))  # support contrast1 vs contrast2
  S_cont1_means <- -0.5*N*(log(r_SS_1) - log(m1$`Sum Sq`[3])) # support contrast1
                                                            # vs means model
  datf <- data.frame(dat,group)
  mean_out <- aggregate(datf[1],datf[2],mean)
  plot(mean_out)

# conventional frequentist statistics
  Fval <- m1$`F value`[2]
  Pval <- m1$`Pr(>F)`[2]

  Fval_c1 <- SS_1/m1$`Mean Sq`[3]
  Pval_c1 <- pf(Fval_c1, 1, m1$Df[3], lower.tail = FALSE)


  if (correct==0) {text="No"
  } else if (correct==1) {text="Occam's bonus"
  } else if (correct==2) {text="AIC"
  } else {text="Undefined"}
  correcttext <- paste(text,"correction was applied")

  if(verb) cat("\nSupport for group means versus null = ", round(S_12c,3), sep= "",
  "\n", correcttext,
  "\n Support for contrast 1 ", if (is.null(nulfg)) "(linear) ",
  "versus group means model = ", round(S_cont1_means,3),
  "\n Support for contrast 1 versus contrast 2 ", if (is.null(nulfg)) "(quadratic) ", "= ", round(S_cont_12,3),
  "\n\nOverall analysis F(",dfv[2],",",dfv[3],") = ", round(Fval,3),
  ", p = ", format.pval(Pval,4), ", partial eta-squared = ", round(eta_sq,3),
  "\nContrast 1 F(1,",dfv[3],") = ", round(Fval_c1,3),
  ", p = ", format.pval(Pval_c1,4), "\n ")

  df1 <- c(1, m1$Df[3])
  invisible(list(S.12 = S_12, S.S_12c = S_12c, capplied = correcttext,
                 S.1m = S_cont1_means, S.cont.12 = S_cont_12,
                 contrast1 = contrast1, contrast2 = contrast2,
                 gp.means = mean_out, df = m1$Df, F.val = Fval, P.val = Pval,
                 eta.sq = eta_sq, Fval.c1 = Fval_c1, df.1 = df1, P.val1 = Pval_c1))
}
