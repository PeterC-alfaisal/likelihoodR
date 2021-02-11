#' Likelihood Supports for One-way Independent Samples ANOVA
#'
#' This function calculates supports for independent samples ANOVA. One support is
#' for the model of group means against the null (no grouping), for the first contrast
#' versus the group means model, and the other for 2 contrasts.
#' Both contrasts should be either NULL or specified.
#' If the contrasts use the default of NULL, then it calculates a linear versus a
#' quadratic contrast. The corrected support is given for groups versus null, using
#' Akaike's correction (Hurvich & Tsai (1989)). No correction is necessary for the
#' two contrasts' support since they both involve 1 parameter.
#' Conventional frequentist F and p value statistics are given for the overall
#' analysis and for contrast 1. Unequal group sizes are
#' accommodated.
#'
#' @usage L_1way_ANOVA(data, group, contrast1=NULL, contrast2=NULL, verb=TRUE)
#'
#' @param data a (non-empty) numeric vector of data values.
#' @param group an integer vector the same length as data, coding for k groups.
#' @param contrast1 first contrast, default = NULL.
#' @param contrast2 second contrast, default = NULL.
#' @param verb show output, default = TRUE.
#'
#' @return
#' $S.12c - corrected support for groups hypothesis versus null.
#'
#' $S.12 - uncorrected support for groups hypothesis versus null.
#'
#' S.1mc - support for contrast 1 versus the group means model.
#'
#' $S.1m = uncorrected support for contrast 1 vs group means.
#'
#' $S.cont.12 - support for contrast 1 versus contrast 2.
#'
#' $contrast1 - first contrast.
#'
#' $contrast2 - second contrast.
#'
#' $gp.means - group means.
#'
#' $df - degrees of freedom for groups and error.
#'
#' $F.val - F value for overall ANOVA analysis.
#'
#' $P.val - p value for overall analysis.
#'
#' $eta.sq - eta-squared.
#'
#' $Fval.c1 - F value for contrast 1.
#'
#' $df.1 - degrees of freedom for contrast.
#'
#' $P.val1 - p value for contrast 1.
#'
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
#' @importFrom stats pf
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
#' samples. Biometrika. 1989; 76(2):297.
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

L_1way_ANOVA <- function(data, group, contrast1=NULL, contrast2=NULL, verb=TRUE) {

  datf <- data.frame(data, group)
  ad <- datf[complete.cases(data), ] # remove missing, NA or NaN, case-wise
  dat <- ad$data
  gp <- ad$group

  m1 <- anova(lm(dat~gp))
  tss <- sum(m1$`Sum Sq`)  # total SS
  eta_sq <- m1$`Sum Sq`[1]/tss  # eta-squared
  N <- (sum(m1$Df)+1)
  k <- m1$Df[1]+1
  S_12 <- -0.5 * N * (log(m1$`Sum Sq`[2]) - log(tss))
  k2 <- 2       # parameters for variance and grand mean
  k1 <- m1$Df[1] + k2

# Akaike's correction
  Ac <- function(k1,k2,N) { k2 * N/(N - k2 - 1) - k1 * (N/(N - k1 - 1)) }
  S_12c <- S_12 + Ac(k1,k2,N)

# contrasts
  nulfg <- 0
  if (is.null(contrast1)) {
    conta <- contr.poly(k, scores = 1:k)
    contrast1 <- conta[,1]
    contrast2 <- conta[,2]
    nulfg <- NULL
  }
  n <- N/(m1$Df[1]+1)
  gp_means <- tapply(dat, gp, mean)
  gp_n <- tapply(dat, gp, length)

  SS_1 <- sum(contrast1*gp_means)^2/(sum(contrast1^2/(gp_n)))
  SS_2 <- sum(contrast2*gp_means)^2/(sum(contrast2^2/(gp_n)))

  r_SS_1 <- tss - SS_1
  r_SS_2 <- tss - SS_2

  S_cont_12 <- -0.5*N*(log(r_SS_1) - log(r_SS_2))  # support contrast1 vs contrast2
  S_cont1_means <- -0.5*N*(log(r_SS_1) - log(m1$`Sum Sq`[2])) # support contrast1
                                                            # vs means model

  k1 <- 1       # parameters for contrast and means model
  k2 <- m1$Df[1] + 2
  S_cont1_meansc <- S_cont1_means + Ac(k1,k2,N)

  mean_out <- aggregate(datf[1],datf[2],mean)

# conventional frequentist statistics
  Fval <- m1$`F value`[1]
  dfv <- m1$Df
  Pval <- m1$`Pr(>F)`[1]
  Fval_c1 <- SS_1/m1$`Mean Sq`[2]
  Pval_c1 <- pf(Fval_c1, 1, m1$Df[2], lower.tail = FALSE)
# Fval_c12 <- SS_1/SS_2                              #comparing 2 contrasts
# Pval_c12 <- pf(Fval_c12, 1, 1, lower.tail = FALSE)

  plot(mean_out)

  if(verb) cat("\nSupport for group means model versus null = ", round(S_12c,3), sep= "",

    "\n Support for contrast 1 ", if (is.null(nulfg)) "(linear) ",
    "versus group means model = ", round(S_cont1_meansc,3),
    "\n Support for contrast 1 versus contrast 2 ", if (is.null(nulfg)) "(quadratic) ", "= ", round(S_cont_12,3),
     "\n\nOverall analysis F(",dfv[1],",",dfv[2],") = ", round(Fval,3),
    ", p = ", round(Pval,5), ", eta-squared = ", round(eta_sq,3),
    "\nContrast 1 F(1,",dfv[2],") = ", round(Fval_c1,3),
    ", p = ", round(Pval_c1,5), "\n ")
  df1 <- c(1, m1$Df[2])
  invisible(list(S.12c = S_12c, S.12 = S_12, S.1mc = S_cont1_meansc,
               S.1m = S_cont1_means, S.cont.12 = S_cont_12,
               contrast1 = contrast1, contrast2 = contrast2,
               gp.means = mean_out, df = m1$Df, F.val = Fval, P.val = Pval,
               eta.sq = eta_sq, Fval.c1 = Fval_c1, df.1 = df1, P.val1 = Pval_c1))

}
