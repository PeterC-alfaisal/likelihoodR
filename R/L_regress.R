#' Likelihood Support for Regression
#'
#' This function calculates the supports for different regression fits from 2 vectors of data.
#' Models include linear, quadratic and cubic (given sufficient data). A plot is
#' included showing linear (black), quadratic (red) and cubic (blue dashed) lines. P values for
#' tehe model fits are also given.
#'
#' @usage L_regress(y, x)
#' @param x a numeric vector.
#' @param y a numeric vector the same length as x.
#'
#' @return
#' $S.LNc - corrected support for linear versus null model.
#'
#' $S.LN - uncorrected support for linear versus null model.
#'
#' $S.QLc - corrected support for quadratic versus linear model.
#'
#' $S.QL - uncorrected support for quadratic versus linear model.
#'
#' S.QCc = support for quadratic versus cubic model.
#'
#' $N - sample size.
#'
#' $p.vals - p values for 3 fits.
#'
#' @keywords Likelihood; support; regression; linear; quadratic; cubic
#'
#' @export
#'
#' @importFrom stats anova
#' @importFrom stats lm
#' @importFrom graphics lines
#'
#' @examples # for women's world record times for 1500m event example, p 108
#' years <- c(0.0,	7.1,	8.9,	8.9,	10.1,	12.8,	17.0,	19.1,
#' 25.0, 28.7, 29.7,	29.9,	35.3, 39.8,	40.2,	41.9,	42.1,	44.0,
#' 44.9, 45.0,	45.1, 45.1,	48.9,	52.9,	53.0,	66.1,	87.9)
#' time <- c(5.30,	5.12,	5.03,	4.79,	4.75,	4.70,	4.63,	4.63,
#' 4.62, 4.59,	4.50,	4.50,	4.32,   4.29,	4.26,	4.21,	4.18,
#' 4.16,	4.12, 4.11,	4.09,	4.02,	3.93,	3.92,	3.87,	3.84,	3.83)
#'
#' m=L_regress(time, years)
#' m
#'
#' @references
#' Cahusac, P.M.B. (2020) Evidence-Based Statistics, Wiley, ISBN : 978-1119549802
#'
#' Royall, R. M. (1997). Statistical evidence: A likelihood paradigm. London: Chapman & Hall, ISBN : 978-0412044113
#'
#' Edwards, A.W.F. (1992) Likelihood, Johns Hopkins Press, ISBN : 978-0801844430
#'

L_regress <- function(y, x) {

  m1=anova(lm(y ~ x))

  tss <- sum(m1$`Sum Sq`)
  N <- (sum(m1$Df)+1)
  lin_df <- m1$Df[1]

# support for linear versus null
  S_LN <- -0.5 * N * (log(m1$`Sum Sq`[2]) - log(tss))

# Akaike's correction
  k2 <- 1       # parameters for grand mean
  k1 <- m1$Df[1] + 1  # parameters for slope and intercept
  Ac <- k1 - k2
  S_LNc <- S_LN - Ac

# examining non-linearity
# comparing quadratic fit to linear
  m3 <- anova(lm(y ~ x + I(x^2) + I(x^3)))
  m3

  quad_ss <- m3$`Sum Sq`[2] + m3$`Sum Sq`[1]
  unex_q_ss <- tss - quad_ss
  S_QL <- -0.5 * N * (log(unex_q_ss) - log(m1$`Sum Sq`[2]))
  S_QLc <- S_QL - 1   # additional parameter for the quadratic

# check
#Rsq_l <- summary(lm(y ~ x))$r.squared
#Rsq_q <- summary(lm(y ~ x + I(x^2)))$r.squared
#(LR_ql <- ((1-Rsq_q)/(1-Rsq_l))^-(N/2))
#S <- log(LR_ql)
#(SC <- S-1)
# alternate calculation for quadratic versus linear, gives same answer

# quadratic versus cubic
  unex_c_ss <- tss - sum(m3$`Sum Sq`[1:3])
  S_QC <- -0.5 * N * (log(unex_q_ss) - log(unex_c_ss))
  S_QCc <- S_QC + 1   # quadratic has one less parameter

# plot
  plot(x, y)
  fit1<-lm(y ~ poly(x,1,raw=TRUE))
  linear = fit1$coefficient[2]*x + fit1$coefficient[1]
  fit2<-lm(y ~ poly(x,2,raw=TRUE))
  quadratic = fit2$coefficient[3]*x^2 + fit2$coefficient[2]*x + fit2$coefficient[1]
  fit3<-lm(y ~ poly(x,3,raw=TRUE))
  cubic = fit3$coefficient[4]*x^3 + fit3$coefficient[3]*x^2 +
          fit3$coefficient[2]*x + fit3$coefficient[1]
  lines(x,linear, col="black")          # linear
  lines(x,quadratic, col="red")         # quadratic
  lines(x,cubic, col="blue",lty=2)      # cubic


  cat("\nSupport for linear fit over null model ", round(S_LNc,3), sep= "",
      "\n Support for quadratic versus linear fit = ", round(S_QLc,3),
      "\n Support for quadratic versus cubic = ", round(S_QCc,3),
       "\n N = ", N,
      "\n P values for linear, quadratic and cubic fits = ", m3$`Pr(>F)`[1],
      "  ", m3$`Pr(>F)`[2], "  ", m3$`Pr(>F)`[3], "\n ")

  invisible(list(S.LNc = S_LNc, S.LN = S_LN, S.QLc = S_QLc, S.QL = S_QL,
               S.QCc = S_QCc, N = N, p.vals = m3$`Pr(>F)`[1:3]))

}
