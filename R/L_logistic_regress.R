#' Likelihood Support for Logistic Regression
#'
#' This function calculates the supports for multiple logistic regression.
#' A binary dependent variable is entered into the function, followed by up to 6 predictor
#' variables (which need to be dummy coded if nominal and more than 2 levels). Outputs
#' give the overall support for the full model versus the null (constant) model, supports
#' for each of the predictor variables. Outputs include the usual chi-squared and p values.
#'
#' @usage L_logistic_regress(yv, p1, p2=NULL, p3=NULL, p4=NULL, p5=NULL, p6=NULL)
#' @param yv a binomial numeric vector for dependent variable.
#' @param p1 vector for predictor variable, same length as yv.
#' @param p2 vector for predictor variable, same length as yv, default = NULL.
#' @param p3 vector for predictor variable, same length as yv, default = NULL.
#' @param p4 vector for predictor variable, same length as yv, default = NULL.
#' @param p5 vector for predictor variable, same length as yv, default = NULL.
#' @param p6 vector for predictor variable, same length as yv, default = NULL.
#'
#' @return
#' $S.FNc - support for the full versus null (constant value) model.
#'
#' $S.each - support for each of the predictors, from first to last.
#'
#' $chi.sq.FN - chi-squared for full versus null model.
#'
#' $df - degrees of freedom for chi-squared.
#'
#' $chi.sq.FN.p - p value for chi-squared.
#'
#' $p.vals - p values for each of the predictors
#'
#'
#' @keywords Likelihood; support; logistic regression
#'
#' @export
#'
#' @importFrom stats as.formula
#' @importFrom stats drop1
#' @importFrom stats pchisq
#'
#' @examples # prescription errors example, p 114
#' p_error <- c(rep(1,6),rep(0,9))
#' score <- c(4,	5,	6,	5,	4,	6,	6,	4,
#' 5,	8,	9,	7,	10,	8,	9)
#' med1 <- c(1,	1,	0,	0,	1,	1,	0,	0,
#' 0,	0,	0,	0,	0,	0,	1)
#' med2 <- c(0,	0,	1,	0,	0,	0,	1,	0,
#' 0,	1,	1,	0,	1,	1,	0)
#' m1 = L_logistic_regress(p_error, score, med1, med2)
#' m1
#'
#' @references
#' Cahusac, P.M.B. (2020) Evidence-Based Statistics, Wiley, ISBN : 978-1119549802
#'
#' 1 Akaike H. A new look at the statistical model identification. IEEE Transactions
#' on Automatic Control. 1974;19(6):716–23.
#'
#' Glover S, Dixon P. Likelihood ratios: a simple and flexible statistic for empirical
#' psychologists. Psychonomic Bulletin & Review. 2004;11(5):791–806.
#'
#' Tabachnick BG, Fidell LS. Using Multivariate Statistics. Boston: Pearson Education; 2007.
#'

L_logistic_regress <- function(yv, p1, p2=NULL, p3=NULL, p4=NULL, p5=NULL, p6=NULL) {

  ya = "yv"
  p1a = "p1"; p2a = "p2"; p3a = "p3"; p4a = "p4"; p5a = "p5"; p6a = "p6"
  if (is.null(p2)) {
    num_preds=1; fo <- as.formula(paste(ya, collapse = "", "~", p1a))
  } else {
    if (is.null(p3)) {
      num_preds=2; fo <- as.formula(paste(ya, collapse = "", "~", p1a, "+", p2a))
    }  else {
      if (is.null(p4)) {
        num_preds=3; fo <- as.formula(paste(ya, collapse = "", "~", p1a, "+", p2a, "+", p3a))
      } else {
        if (is.null(p5)) {
         num_preds=4; fo <- as.formula(paste(ya, collapse = "", "~", p1a, "+", p2a, "+", p3a, "+", p4a))
        } else {
          if (is.null(p6)) {
            num_preds=5; fo <- as.formula(paste(ya, collapse = "", "~", p1a, "+", p2a, "+", p3a, "+", p4a, "+", p5a))
            } else {
              num_preds=6; fo <- as.formula(paste(ya, collapse = "", "~", p1a, "+", p2a, "+", p3a, "+", p4a, "+", p5a, "+", p6a))
  }}}}}

  fit.1 <- do.call("glm", list(fo, quote(binomial)))

  m = drop1(fit.1, test="Chisq")

# full versus null model
  S_FNc <- -fit.1$deviance/2 - (-fit.1$null.deviance/2) - num_preds # correcting for parameters
  X_FN <- -fit.1$deviance - (-fit.1$null.deviance)
  X_FN_p <- 1-pchisq(X_FN, df =  num_preds)

  Pv = m$`Pr(>Chi)`
  S = 0
  for (i in 1:num_preds) {
    S[i] <- (m$AIC[i+1]-m$AIC[1])/2     #comparing without 1st factor
  }
  cat("\nSupport for full over null model ", round(S_FNc,3), sep= "",
      "\n Support for 1st predictor versus the full model = ", round(S[1],3),
      "\n Support for 2nd predictor versus the full model = ", round(S[2],3),
      "\n Support for 3rd predictor versus the full model = ", round(S[3],3),
      "\n Support for 4th predictor versus the full model = ", round(S[4],3),
      "\n Support for 5th predictor versus the full model = ", round(S[5],3),
      "\n Support for 6th predictor versus the full model = ", round(S[6],3),
      "\n\nFor full over null chi-squared(", num_preds,") = ", X_FN, " p = ", round(X_FN_p,5),
      "\nP values for each of the predictors versus the full model = ",
      round(Pv[2], 5), "  ", round(Pv[3], 5), "  ",
      round(Pv[4], 5), "  ", round(Pv[5], 5), "  ",
      round(Pv[6], 5), "  ", round(Pv[5], 5), "\n")

  invisible(list(S.FNc = S_FNc, S.each = S, chi.sq.FN = X_FN, df = num_preds,
               chi.sq.FN.p = X_FN_p, p.vals = Pv[1:num_preds+1]))

}
