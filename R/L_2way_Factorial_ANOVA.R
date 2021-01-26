#' Likelihood Supports for Two-way Independent Samples Factorial ANOVA
#'
#' This function calculates supports for independent samples ANOVA. One support is
#' for the full model versus null (no factors), and the second is for full model
#' versus main effects. Two contrasts can be specified. Each should be given
#' as a vector arranged as means for factor1 changing first (see example). If only the first
#' contrast is specified then this is compared to the main effects model. If a second
#' contrast is specified then the first contrast is compared to it.
#' Corrected support is given where appropriate, using
#' Akaike's correction (Hurvich & Tsai (1989)). No correction is necessary for the
#' two contrasts support since they both involve 1 parameter. Balanced designs only.
#'
#' @usage L_2way_Factorial_ANOVA(data, factor1, factor2, contrast1=NULL, contrast2=NULL)
#'
#' @param data a (non-empty) numeric vector of data values.
#' @param factor1 a vector the same length as data, coding the first factor.
#' @param factor2 a vector the same length as data, coding the second factor.
#' @param contrast1 first contrast, with values for factor1 changing first, default = NULL.
#' @param contrast2 second contrast, default = NULL.
#'
#' @return
#' $S.12c - corrected support for full model versus null.
#'
#' $S.12 - uncorrected support for full model versus null.
#'
#' $df - degrees of freedom for the ANOVA.
#'
#' $S_FMc - corrected support for full model versus main effects model.
#'
#' $S.FM - uncorrected support for full versus main effects.
#'
#' $S.c1.Mc - corrected support for first contrast versus main effects model.
#'
#' $S.c1.M - uncorrected support for first contrast versus main effects.
#'
#' $S.c1.c2 - support for first versus second contrast, $Means - 2 way table of means.
#'
#' @keywords Likelihood; support; independent samples factorial ANOVA; contrasts
#'
#' @export
#'
#' @importFrom stats anova
#' @importFrom stats lm
#' @importFrom stats aggregate
#' @importFrom stats cor
#' @importFrom stats interaction.plot
#'
#' @examples # blood clotting times example, p 91
#' time <- c(6.4,	4.6,	6.4,	5.6,	5.9, 6.1,	6.3,	4.5,
#' + 4.8,	6.6, 7,	9.3,	7.9,	9.4,	8.2, 4.4,	4.2,	5,
#' + 6.9,	4.5, 4,	4.3,	6.9,	5.5,	5.8,
#' + 4.4,	4.2,	5.1,	6.9,	4.5)
#' Treatment = gl(3,5,30, labels=c("T1","T2","T3"))
#' Health = gl(2,15,30, labels=c("Hemophiliac","Normal"))
#'
#' L_2way_Factorial_ANOVA(time, Treatment, Health)
#'
#' contrast1 <- c(-1, -1, 5,
#'                -1, -1, -1) # interaction Hemo T3 higher than others
#' L_2way_Factorial_ANOVA(time, Treatment, Health, contrast1)
#'
#' contrast2 <- c(-1, -1, -1,
#'                1, 1, 1) # main effect of health status (Hemo higher than Normal)
#'
#' m=L_2way_Factorial_ANOVA(time, Treatment, Health, contrast1, contrast2)
#' m     #show outputs
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

L_2way_Factorial_ANOVA <- function(data, factor1, factor2, contrast1=NULL, contrast2=NULL) {

  dat <- data
  F1 <- factor1
  F2 <- factor2

  m1=anova(lm(dat~F1*F2))

  tss <- sum(m1$`Sum Sq`)
  N <- (sum(m1$Df)+1)
  ka <- m1$Df[1]+1
  kb <- m1$Df[2]+1

  S_12 <- -0.5 * N * (log(m1$`Sum Sq`[4]) - log(tss))
  k2 <- 2       # parameters for variance and grand mean
  k1 <- sum(m1$Df[1:3]) + k2

  # Akaike's correction
  Ac <- function(k1,k2,N) { k2 * N/(N - k2 - 1) - k1 * (N/(N - k1 - 1)) }

  S_12c <- S_12 + Ac(k1,k2,N)

  # full model versus main effects only
  S_FM <- -0.5 * N * (log(m1$`Sum Sq`[4]) - log(sum(m1$`Sum Sq`[3:4])))
  k2 <- 2 + sum(m1$Df[1:2])
  S_FMc <- S_FM + Ac(k1,k2,N)

  datf <- data.frame(dat, F1, F2)
  mean_out <- aggregate(datf[1], by=c(datf[2], datf[3]), mean)

  S_c1M <- NULL
  S_c1Mc <- NULL
  S_c1c2 <- NULL

  if (!is.null(contrast1)) {
  # compare contrast versus main effects
    means <- unlist(mean_out[3])
    alerting <- (cor(means,contrast1))^2
    SS_cont1 <- alerting * (tss - m1$`Sum Sq`[4])
    SS_c1_unex <- tss - SS_cont1
    S_c1M <- -0.5 * N * (log(SS_c1_unex) - log(sum(m1$`Sum Sq`[3:4])))
    k1 <- 2 + 1    #2 parameters (variance and grand mean) + one for contrast
    k2 <- 2 + sum(m1$Df[1:2])
    S_c1Mc <- S_c1M + Ac(k1,k2,N)

    if (!is.null(contrast2)) {
      alerting <- (cor(means,contrast2))^2
      SS_cont2 <- alerting * (tss - m1$`Sum Sq`[4])
      SS_c2_unex <- tss - SS_cont2
      S_c1c2 <- -0.5 * N * (log(SS_c1_unex) - log(SS_c2_unex))
    # no correction necessary using contrasts
      }
    }


  interaction.plot(F1, F2, dat, ylab="Means")


  cat("Support for full model versus null = ", round(S_12c,3), sep= "",
      "\n Support for full model versus main effects = ", round(S_FMc,3),
      "\n Support for contrast1 versus main effects = ", if (!is.null(contrast1)) round(S_c1Mc,3),
      "\n Support for contrast1 versus contrast2 = ", if (!is.null(contrast2)) round(S_c1c2,3))

  invisible(list(S.12c = S_12c, S.12 = S_12, df = m1$Df, S.FMc = S_FMc, S.FM = S_FM,
                 S.c1.Mc = S_c1Mc, S.c1.M = S_c1M, S.c1.c2 = S_c1c2, Means = mean_out))
}
