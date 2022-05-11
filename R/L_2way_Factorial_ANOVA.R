#' Likelihood Supports for Two-way Independent Samples Factorial ANOVA
#'
#' This function calculates supports for independent samples ANOVA. One support is
#' for the full model versus null (no factors), and the second is for full model
#' versus main effects. Two contrasts can be specified which can be used
#' to explore interactions. Each should be given
#' as a vector arranged as means for factor1 changing first (see example). If only the first
#' contrast is specified then this is compared to the main effects model. If a second
#' contrast is specified then the first contrast is compared to it.
#' Corrected support is given where appropriate, using
#' Akaike's correction (Hurvich & Tsai (1989)). No correction is necessary for the
#' two contrasts support since they both involve 1 parameter. Unequal group sizes are
#' accommodated, using type III sums of squares. F, p and partial eta-squared values
#' are given for the two factors and their interaction.
#'
#' @usage L_2way_Factorial_ANOVA(data, factor1, factor2, contrast1=NULL, contrast2=NULL, verb=TRUE)
#'
#' @param data a (non-empty) numeric vector of data values.
#' @param factor1 a vector the same length as data, coding the first factor.
#' @param factor2 a vector the same length as data, coding the second factor.
#' @param contrast1 first contrast, with values for factor1 changing first, default = NULL.
#' @param contrast2 second contrast, default = NULL.
#' @param verb show output, default = TRUE.
#'
#' @return
#' $S.12c - corrected support for full model versus null.
#'
#' $S.12 - uncorrected support for full model versus null.
#'
#' $S_FMc - corrected support for full model versus main effects model.
#'
#' $S.FM - uncorrected support for full versus main effects.
#'
#' $S.c1.Mc - corrected support for first contrast versus main effects model.
#'
#' $S.c1.M - uncorrected support for first contrast versus main effects.
#'
#' $S.c1.c2 - support for first versus second contrast.
#'
#' $Means - 2 way table of means.
#'
#' $df - degrees of freedom for the ANOVA.
#'
#' $F.f1 - F value for first factor main effect.
#'
#' $Pval.f1 - P value for first factor main effect.
#'
#' $eta.sq.1 - partial eta-squared for first factor main effect.
#'
#' $F.f2 - F value for second factor main effect.
#'
#' $Pval.f2 - P value for second factor main effect.
#'
#' $eta.sq.2 - partial eta-squared for second factor main effect.
#'
#' $F.int - F value for interaction.
#'
#' $Pval.int - P value for interaction.
#'
#' $eta.sq.12 - partial eta-squared for the interaction.
#'
#' $F.val.c1 - F value for first contrast.
#'
#' $P.val.c1 - P value for first contrast.
#'
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
#' @importFrom stats pf
#'
#' @examples # blood clotting times example, p 91
#' time <- c(6.4,	4.6,	6.4,	5.6,	5.9, 6.1,	6.3,	4.5,
#' 4.8,	6.6, 7,	9.3,	7.9,	9.4,	8.2, 4.4,	4.2,	5,
#' 6.9,	4.5, 4,	4.3,	6.9,	5.5,	5.8,
#' 4.4,	4.2,	5.1,	6.9,	4.5)
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

L_2way_Factorial_ANOVA <- function(data, factor1, factor2, contrast1=NULL, contrast2=NULL, verb=TRUE) {

  dat <- data
  F1 <- factor1
  F2 <- factor2

  options(contrasts = c("contr.sum","contr.poly")) # for type III SS
  model <- lm(dat~F1*F2)
  m1 <- drop1(model, .~., test="F")

  tss <- sum(m1$`Sum of Sq`[2:4],m1$RSS[1])
  eta_sq_1 <- m1$`Sum of Sq`[2]/(m1$`Sum of Sq`[2] + m1$RSS[1]) # partial eta-squared
  eta_sq_2 <- m1$`Sum of Sq`[3]/(m1$`Sum of Sq`[3] + m1$RSS[1])
  eta_sq_12 <- m1$`Sum of Sq`[4]/(m1$`Sum of Sq`[4] + m1$RSS[1])
  N <- length(dat)
  ka <- m1$Df[2]+1
  kb <- m1$Df[3]+1

  S_12 <- -0.5 * N * (log(m1$RSS[1]) - log(tss))
  k2 <- 2       # parameters for variance and grand mean
  k1 <- sum(m1$Df[2:4]) + k2

  # Akaike's correction
  Ac <- function(k1,k2,N) { k2 * N/(N - k2 - 1) - k1 * (N/(N - k1 - 1)) }

  S_12c <- S_12 + Ac(k1,k2,N)

  # full model versus main effects only
  S_FM <- -0.5 * N * (log(m1$RSS[1]) - log(sum(m1$`Sum of Sq`[4],m1$RSS[1])))
  k2 <- 2 + sum(m1$Df[2:3])
  S_FMc <- S_FM + Ac(k1,k2,N)

  datf <- data.frame(dat, F1, F2)
  mean_out <- aggregate(datf[1], by=c(datf[2], datf[3]), mean)
  gp_n <- aggregate(datf[1], by=c(datf[2], datf[3]), FUN = length)

  S_c1M <- NULL
  S_c1Mc <- NULL
  S_c1c2 <- 0
  SS_cont1 <- NULL

  if (!is.null(contrast1)) {
  # compare contrast versus main effects
    means <- unlist(mean_out[3])
    numbers <- unlist(gp_n[3])
    SS_cont1 <- sum(contrast1*means)^2/(sum(contrast1^2/(numbers)))
    SS_c1_unex <- tss - SS_cont1

    S_c1M <- -0.5 * N * (log(SS_c1_unex) - log(sum(m1$`Sum of Sq`[4],m1$RSS[1])))
    k1 <- 2 + 1    #2 parameters (variance and grand mean) + one for contrast
    k2 <- 2 + sum(m1$Df[2:3])
    S_c1Mc <- S_c1M + Ac(k1,k2,N)

    if (!is.null(contrast2)) {
      SS_cont2 <- sum(contrast2*means)^2/(sum(contrast2^2/(numbers)))
      SS_c2_unex <- tss - SS_cont2
      S_c1c2 <- -0.5 * N * (log(SS_c1_unex) - log(SS_c2_unex))
    # no correction necessary using contrasts
      }
    }

  interaction.plot(F1, F2, dat, ylab="Means")

  Fval <- m1$`F value`[2:4]
  Pval <- m1$`Pr(>F)`[2:4]
  dfv <- m1$Df[2:4]
  dfres <- (N-sum(dfv)-1)

  res_msq <- m1$RSS[1]/dfres
  Fval_c1 <- unname(SS_cont1/res_msq)
  Pval_c1 <- pf(Fval_c1, 1, dfres, lower.tail = FALSE)

  if(verb) cat("\nSupport for full model (including interaction) versus null = ", round(S_12c,3), sep= "",
      "\n Support for full model versus main effects = ", round(S_FMc,3),
      "\n Support for contrast 1 versus main effects = ", if (!is.null(contrast1)) round(S_c1Mc,3),
      "\n Support for contrast 1 versus contrast 2 = ", if (!is.null(contrast2)) round(S_c1c2,3),
      "\n\nFirst factor main effect F(", dfv[1],",", dfres,") = ", round(Fval[1],3),
      ", p = ", Pval[1], ", partial eta-squared = ", round(eta_sq_1,3),
      "\n Second factor main effect F(", dfv[2],",", dfres,") = ", round(Fval[2],3),
      ", p = ", Pval[2], ", partial eta-squared = ", round(eta_sq_2,3),
      "\n Interaction F(", dfv[3],",", dfres,") = ", round(Fval[3],3),
      ", p = ", Pval[3],  ", partial eta-squared = ", round(eta_sq_12,3),
      "\n Contrast 1 F(1,",dfres,") = ", if (!is.null(contrast1)) round(Fval_c1,3),
      ", p = ", Pval_c1, "\n ")

  invisible(list(S.12c = S_12c, S.12 = S_12, S.FMc = S_FMc, S.FM = S_FM,
                 S.c1.Mc = S_c1Mc, S.c1.M = S_c1M, S.c1.c2 = S_c1c2, Means = mean_out, df = c(m1$Df[2:4],dfres),
                 F.f1 = Fval[1], Pval.f1 = Pval[1], eta.sq.1 = eta_sq_1, F.f2 = Fval[2],
                 Pval.f2 = Pval[2], eta.sq.2 = eta_sq_2,
                 F.int = Fval[3], Pval.int = Pval[3], eta.sq.12 = eta_sq_12,
                 F.val.c1 = Fval_c1, P.val.c1 = Pval_c1))
}

