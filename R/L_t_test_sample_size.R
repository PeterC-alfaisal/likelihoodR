#' Sample size calculation using the evidential approach for t tests
#'
#' This function calculates
#' the required sample size for t tests. The standard deviation and effect size are
#' specified. Calculations given for one sample and independent samples t tests. For a related
#' samples test calculation use the sd for paired differences.
#'
#' @usage L_t_test_sample_size(MW = 0.05, sd = 1, d = 1.2, S = 3, paired = FALSE)
#'
#' @param MW set M1 + W1 probability, default = .05.
#' @param sd set standard deviation, default = 1.
#' @param d  set desired effect size, default = 1.2.
#' @param S  set strength of evidence (support), default = 3.
#' @param paired set to TRUE for one sample and FALSE for independent samples, default = FALSE.
#'
#' @return
#' $N - required sample size.
#'
#' $S - specified strength (support) for evidence from the test.
#'
#' $sd - specified standard deviation.
#'
#' $d - Cohen's effect size specified.
#'
#' $m1.w1 - specified probability for combined misleading and weak evidence.
#'
#'
#' @keywords Likelihood; t test; sample size; weak evidence; misleading evidence
#'
#' @export
#'
#' @importFrom stats t.test
#' @importFrom stats qt
#' @importFrom stats sd
#' @importFrom stats optimize
#'
#' @examples # for one sample or related samples (differences)
#' v = L_t_test_sample_size(MW = 0.2, sd = 1, d = 1, S = 3, paired = TRUE)
#' v
#' # for 2 independent samples
#' v = L_t_test_sample_size(MW = 0.05, sd = 1, d = 1.2, S = 3, paired = FALSE)
#' v
#'
#' @references
#' Cahusac, P.M.B. (2020) Evidence-Based Statistics, Wiley, ISBN : 978-1119549802
#'
#' Royall, R. (2000). "On the Probability of Observing Misleading Statistical Evidence." Journal of the
#' American Statistical Association 95(451): 760-768.
#'
#' Royall, R. (2004). The Likelihood paradigm for statistical evidence. The Nature of Scientific Evidence.
#' M. L. Taper and S. R. Lele. Chicago, University of Chicago: 119-152.
#'
#' Royall, R. M. (1997). Statistical evidence: A likelihood paradigm. London: Chapman & Hall, ISBN : 978-0412044113
#'
#' Edwards, A.W.F. (1992) Likelihood, Johns Hopkins Press, ISBN : 978-0801844430
#'
#'

L_t_test_sample_size <- function(MW = 0.05, sd = 1, d = 1.2, S = 3, paired = FALSE) {

  toler <- 0.0001   # for optimize function
  ns <- 0; nul <- 0
  to <- vector(); tmw <- vector(); td <- vector()

  if(paired) {

    f <- function(n, MW, nul, d, sd, S) {
      tmw <- qt(MW,df=n-1, lower.tail =TRUE)
      Xu <- nul+(d*sd*(1-sqrt((exp(-(2*S/n))-(1-exp(-(2*S/n)))^2*(n-1)/
                                 (n*d^2))))/((1-exp(-(2*S/n)))))
      to <- (Xu-d*sd)*sqrt(n)/sd
      td <- (to-tmw)^2
    }
    xmin1 <- optimize(f, c(3, 1000), tol = toler, MW, nul, d, sd, S)
    ns <- round(xmin1$minimum,0)

    cat("\nFor 1 sample t test with M1 + W1 probability of ", MW,
        "\n Strendth of evidence required = ", S,
        ", and effect size = ", d, "\n Required sample size = ", ns, sep ="")
    return(invisible(list(N = ns, S = S, sd = sd, d = d, m1.w1 = MW)))

  } else {

    g <- function(n, MW, nul, d, sd, S) {
      tmw <- qt(MW,df=n-2, lower.tail =TRUE)
      Xu <- nul+(d*sd*(1-sqrt((exp(-(2*S/n))-(1-exp(-(2*S/n)))^2*(n-2)/
                                 (n*d^2))))/((1-exp(-(2*S/n)))))
      to <- (Xu-d*sd)*sqrt(n)/(2*sd)
      td <- (to-tmw)^2
    }
    xmin1 <- optimize(g, c(3, 1000), tol = toler, MW, nul, d, sd, S)
    ns <- round(xmin1$minimum,0)

    if (lapply(ns, "%%", 2) != 0) ns <- ns+1 # make even N for 2 samples
    cat("\nFor independent samples t test with M1 + W1 probability of ", MW,
        "\n Strendth of evidence required = ", S,
        ", and effect size = ", d, "\n Required total sample size (n1 + n2) = ", ns, sep ="")
    return(invisible(list(N = ns, S = S, sd = sd, d = d, m1.w1 = MW)))
  }
}
