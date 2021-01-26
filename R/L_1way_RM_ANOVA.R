#' Likelihood Supports for One-way Repeated Measures ANOVA
#'
#' This function calculates support for the treatment means versus the null model.
#' No correction is made for treatment degrees of freedom.
#'
#'
#' @usage L_1way_RM_ANOVA(dat, group, ID)
#'
#' @param dat a (non-empty) numeric vector of data values.
#' @param group a vector the same length as data, coding for k groups.
#' @param ID is an identifier for each case.
#'
#' @return
#' $S.12 - support for groups means hypothesis versus null.
#'
#' $df - degrees of freedom for ANOVA.
#'
#' $gp.means - group means.
#'
#' @keywords Likelihood; support; repeated measures ANOVA
#'
#' @export
#'
#' @importFrom stats anova
#' @importFrom stats lm
#' @importFrom stats aggregate
#'
#' @examples # sleep data example, p 97
#' dat <- c(0.7,	-1.6,	-0.2,	-1.2,	-0.1,	3.4,	3.7,	0.8,	0,	2,
#' 1.9,	0.8,	1.1,	0.1,	-0.1,	4.4,	5.5,	1.6,	4.6,	3.4,
#' 1.5,	1.4,	0.0,	-0.7,	0.5,	5.1,	5.7,	1.5,	4.7,	3.5)
#' treat <- gl(3,10,30)
#' patients <- gl(10,1,30)
#'
#' m=L_1way_RM_ANOVA(dat, treat, patients)
#' m
#'
#' @references Cahusac, P.M.B. (2020) Evidence-Based Statistics, Wiley, ISBN : 978-1119549802
#'
#' Dixon P. The effective number of parameters in post hoc models. Behavior
#' Research Methods. 2013; 45(3):604–12.
#'
#' Dixon P. The p-value fallacy and how to avoid it. Canadian Journal of
#' Experimental Psychology/Revue canadienne de psychologie expérimentale. 2003;
#' 57(3):189–202.
#'
#' Glover S, Dixon P. Likelihood ratios: a simple and flexible statistic for empirical
#' psychologists. Psychonomic Bulletin and Review. 2004; 11(5):791–806.
#'
#'

L_1way_RM_ANOVA <- function(dat, group, ID) {

m1=anova(lm(dat~ID + group))

within_ss <- sum(m1$`Sum Sq`[2:3])
N <- m1$Df[1]+1
S_12 <- -0.5 * N * (log(m1$`Sum Sq`[3]) - log(within_ss))

datf <- data.frame(dat,group)
mean_out <- aggregate(datf[1],datf[2],mean)
plot(mean_out)

cat("Support for group means versus null = ", round(S_12,3), sep= "")

invisible(list(S.12 = S_12, df = m1$Df, gp.means = mean_out))
}
