#' Likelihood Support for Two-way Categorical Data
#'
#' This function calculates supports for two-way categorical data. This consists of the
#' support for the interaction and the two main effects. Support for the interaction
#' being closer than expected (Edwards p 187, Cahusac p 158) is calculated. The support
#' for trend across the columns is given (assuming the levels for columns are ordered),
#' and conventional p value for trend.
#' Finally, the usual chi-squared statistic and p value are given.
#'
#' @usage L_2way_cat(table)
#' @param table a 2 x 2 matrix or contingency table containing counts.
#' @return
#' $S.int - support for the interaction.
#'
#' $df - the degrees of freedom for the interaction.
#'
#' $S.int.unc - the uncorrected support for the interaction.
#'
#' $S.Main.rows - support for the rows main effect.
#'
#' $S.Main.cols - support for the columns main effect.
#'
#' $S.Mr.uncorr - uncorrected support for rows main effect.
#'
#' $S.Mc.uncorr - uncorrected support for the columns main effect.
#'
#' $df.rows - degrees of freedom for rows.
#'
#' $df.cols - degrees of freedom for columns.
#'
#' $S.total - support for the whole table.
#'
#' $S.trend - support for the trend across columns (if ordered).
#'
#' $too.good - support for the interaction being closer than expected (Edwards p 187, Cahusac p 158).
#'
#' $observed - the observed table frequencies.
#'
#' $expected - the expected values for null hypothesis of no interaction.
#'
#' $residuals - the Pearson residuals.
#'
#' $chi.sq = the chi-squared statistic.
#'
#' $p.value - the p value associated with the chi-squared statistic.
#'
#' $trend.p - p value for trend (from chi-squared dist.).
#'
#' @keywords Likelihood-based; 2-way contingency table
#'
#' @export
#'
#' @importFrom stats prop.trend.test
#' @importFrom stats chisq.test
#'
#' @examples # S. mansoni eggs in stools example, p 151
#' eggs <- as.table(rbind(c(14, 16, 14, 7, 6), c(87, 33, 66, 34, 11)))
#' dimnames(eggs) = list("Infested" = c("Positive","Negative"),
#'                   "Age Group" = c("0-","10-", "20-",
#'                   "30-", "40-"))
#' L_2way_cat(eggs)
#'
#' # or as a matrix
#' eggs <- as.matrix(c(14, 87, 16, 33, 14, 66, 7, 34, 6, 11))
#' dim(eggs) <- c(2,5)
#' L_2way_cat(eggs)
#'
#'
#' @references
#'
#' Cahusac, P.M.B. (2020) Evidence-Based Statistics, Wiley, ISBN : 978-1119549802
#'
#' Royall, R. M. (1997). Statistical evidence: A likelihood paradigm. London: Chapman & Hall, ISBN : 978-0412044113
#'
#' Edwards, A.W.F. (1992) Likelihood, Johns Hopkins Press, ISBN : 978-0801844430


L_2way_cat <- function(table) {

# calculating the interaction
  S2way <- 0
  suppressWarnings(lt <- chisq.test(table)) # ignore warning message
# lt$observed * log(lt$observed/lt$expected) # individual S terms
  S2way <- sum( lt$observed * log(lt$observed/lt$expected) )
  df <- unname(lt$parameter)
  S2wayc <- S2way - (df-1)/2
# main marginal totals

  row_sum <- rowSums(table)
  col_sum <- colSums(table)
  grandtot <- sum(table)
  RowMain <- sum(row_sum*log(row_sum))-grandtot*log(grandtot) + grandtot*log(length(row_sum))
  RowMain_c <- RowMain - ((length(row_sum)-1)-1)/2 # corrected for row df

  ColMain <- sum(col_sum*log(col_sum))-grandtot*log(grandtot) + grandtot*log(length(col_sum))
  ColMain_c <- ColMain - ((length(col_sum)-1)-1)/2 # corrected for column df

# Total S
  Tot_S <- sum(table*log(table))-sum(table)*
    log(sum(table)/length(table))
# same as components added together (without correction for df)
# S2way + RowMain + ColMain

  chi.s <- unname(lt$statistic)
  toogood <- df/2*(log(df/chi.s)) - (df - chi.s)/2

# evidence for trend
  table.pos <- table[1:length(col_sum)]
  trX <- prop.trend.test(table[1,], col_sum)
  tr <- unname(trX$statistic)/2      # S for trend

  cat("\nSupport for interaction corrected for ", df, " df = ", round(S2wayc,3), sep= "",
    "\n Support for rows main effect corrected for ",
    length(row_sum)-1, " df = ", round(RowMain_c,3),
    "\n Support for columns main effect corrected for ", length(col_sum)-1, " df = ",
    round(ColMain_c,3), "\n Total support for whole table = ", round(Tot_S,3),
    "\n Support for trend across columns = ", round(tr,3),
    "\n Support for data fitting too well = ", round(toogood,3),
    "\n\n P value for trend from chi-squared = ", trX$p.value,
    "\n Chi-square(", df, ") = ", round(chi.s,3),
    ",  p = ", round(lt$p.value,5), ", N = ", grandtot)
invisible(list(S.int = S2wayc, df = df, S.int.unc = S2way,
               S.Main.rows = RowMain_c, S.Main.cols = ColMain_c,
               S.Mr.uncorr = RowMain, S.Mc.uncorr = ColMain,
               df.rows = length(row_sum)-1, df.cols = length(col_sum)-1,
               S.total = Tot_S, S.trend = tr, too.good = toogood,
               observed = lt$observed, expected = lt$expected,
               residuals = lt$residuals,
                      chi.sq = chi.s, p.value = lt$p.value, trend.p = trX$p.value))
}
