#' Likelihood Supports for One-way Independent Samples ANOVA
#'
#' This function calculates supports for independent samples ANOVA. One support is
#' for the model of group means against the null (no grouping), and the other for 2 contrasts,
#' both of which should be either NULL or specified.
#' If the contrasts use the default of NULL, then it calculates a linear versus a
#' quadratic contrast. The corrected support is given for groups versus null, using
#' Akaike's correction (Hurvich & Tsai (1989)). No correction is necessary for the
#' two contrasts support since they both involve 1 parameter. Unequal group sizes are
#' accommodated.
#'
#' @usage L_1way_ANOVA(data, group, contrast1=NULL, contrast2=NULL)
#'
#' @param data a (non-empty) numeric vector of data values.
#' @param group an integer vector the same length as data, coding for k groups.
#' @param contrast1 first contrast, default = NULL.
#' @param contrast2 second contrast, default = NULL.
#'
#' @return
#' $S.12c - corrected support for groups hypothesis versus null.
#'
#' $S.12 - uncorrected support.
#'
#' $df - degrees of freedom for ANOVA.
#'
#' $S.cont.12 - support for contrast 1 versus contrast 2.
#'
#' $gp.means - group means.
#'
#' @keywords Likelihood; support; independent samples ANOVA; contrasts
#'
#' @export
#'
#' @importFrom stats anova
#' @importFrom stats lm
#' @importFrom stats contr.poly
#' @importFrom stats aggregate
#' @importFrom stats complete.cases
#'
#' @examples # fitness example, p 81
#' dat <- c(7,	5,	9,	8,	3,	12,	10,	8,	7,	9,
#' 5,	7,	7,	6,	4,	8,	12,	9,	7,	8,
#' 3,	2,	7,	6,	8,	6,	5,	3,	4,	3,
#' 4,	3,	3,	1,	2,	5,	7,	6,	8,	7)
#' gp <- as.factor(rep(1:4,each=10))
#' gp = gl(4,10,40, labels=c("6 hr",	"3 hr",	"Sports club",	"Video games"))
#' contrast1 <- c(-3, -1, 1, 3)  # linear
#' contrast2 <- c(1, -1, -1, 1)  # quadratic
#' contrast3 <- c(1, 1, -1, -1)
#' L_1way_ANOVA(dat,gp,contrast3, contrast1)
#'
#' @references Cahusac, P.M.B. (2020) Evidence-Based Statistics, Wiley, ISBN : 978-1119549802
#'
#' Hurvich CM, Tsai C-L. Regression and time series model selection in small
#' samples. Biometrika. 1989; 76(2):297–307.
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

L_1way_ANOVA <- function(data, group, contrast1=NULL, contrast2=NULL) {

  datf <- data.frame(data, group)
  ad <- datf[complete.cases(data), ] # remove missing, NA or NaN, case-wise
  dat <- ad$data
  gp <- ad$group

m1 <- anova(lm(dat~gp))
tss <- sum(m1$`Sum Sq`)
N <- (sum(m1$Df)+1)
k <- m1$Df[1]+1
S_12 <- -0.5 * N * (log(m1$`Sum Sq`[2]) - log(tss))
k2 <- 2       # parameters for variance and grand mean
k1 <- m1$Df[1] + k2

# Akaike's correction
Ac <- function(k1,k2,N) { k2 * N/(N - k2 - 1) - k1 * (N/(N - k1 - 1)) }
S_12c <- S_12 + Ac(k1,k2,N)

# contrasts
if (is.null(contrast1)) {
  conta <- contr.poly(k, scores = 1:k)
  contrast1 <- conta[,1]
  contrast2 <- conta[,2]
}
n <- N/(m1$Df[1]+1)
gp_means <- tapply(dat, gp, mean)
gp_n <- tapply(dat, gp, length)

SS_1 <- sum(contrast1*gp_means)^2/(sum(contrast1^2/(gp_n)))
SS_2 <- sum(contrast2*gp_means)^2/(sum(contrast2^2/(gp_n)))

r_SS_1 <- tss - SS_1
r_SS_2 <- tss - SS_2

S_cont_12 <- -0.5*N*(log(r_SS_1) - log(r_SS_2))

mean_out <- aggregate(datf[1],datf[2],mean)

cat("Support for group means versus null = ", round(S_12c,3), sep= "",
    "\n Support for contrast 1 versus contrast 2 = ", round(S_cont_12,3))

invisible(list(S.12c = S_12c, S.12 = S_12, df = m1$Df, S.cont.12 = S_cont_12,
               gp.means = mean_out))

}
