#' @title Quantile-quantile plot 
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
  ggplot(df) +
    geom_ribbon(
      mapping = aes(x = expected, ymin = clower, ymax = cupper),
      alpha = 0.1
    ) +
    geom_point(aes(expected, observed), size = 1.2) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    # geom_line(aes(expected, cupper), linetype = 2, size = 0.5) +
    # geom_line(aes(expected, clower), linetype = 2, size = 0.5) +
    xlab(log10Pe) + ylab(log10Po) + 
    annotate(geom = "text", x = min(expected), y = max(observed), hjust = 0, 
             label = sprintf("λ = %.2f", inflation(ps)),
             size = 4)
}
