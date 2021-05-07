#' 
#' @title Create a quantile-quantile plot with ggplot2.
#' @description 
#' Assumptions:
#' 
#'   - Expected P values are uniformly distributed.\cr
#'   - Confidence intervals assume independence between tests.\cr
#'     We expect deviations past the confidence intervals if the tests are
#'     not independent.
#'     For example, in a genome-wide association study, the genotype at any
#'     position is correlated to nearby positions. Tests of nearby genotypes
#'     will result in similar test statistics.
#'
#' @description Creates a QQ-plot to visualize inflation
#' 
#' @param ps Vector of p-values.
#' @param ci Size of the confidence interval, 95\% by default.
#' @return A ggplot2 plot.
#' @export
#' @import ggplot2
#' 

qqplot <- function(ps, ci = 0.95) {
  
  inflation <- function(x) {
    chisq <- qchisq(1 - x, 1)
    lambda <- median(chisq) / qchisq(0.5, 1)
    lambda
  }
  
  n  <- length(ps)
  observed <- -log10(sort(ps))
  expected <- -log10(ppoints(n))
  df <- data.frame(
    observed = observed,
    expected = expected,
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], plain((p.value))))
  log10Po <- expression(paste("Observed -log"[10], plain((p.value))))
  ggplot2::ggplot(df) +
    ggplot2::geom_ribbon(
      mapping = aes(x = expected, ymin = clower, ymax = cupper),
      alpha = 0.1
    ) +
    ggplot2::geom_point(aes(expected, observed), size = 1.2) +
    ggplot2::geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    ggplot2::xlab(log10Pe) + ggplot2::ylab(log10Po) + 
    ggplot2::annotate(geom = "text", x = min(expected), y = max(observed), hjust = 0, 
              label = paste("lambda == ", round(inflation(ps),2)),
              size = 4, parse = TRUE)
}
