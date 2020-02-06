##' Creates a QQ-plot from omic data analysis
##' 
##' @title Generate a QQ-plot of p-values output from any omic data analysis
##' @description The function produces a QQ-plot of the observed vs expected p-values (with confidence bands)
##' @param x a "ds.dsLmFeature" type object containing a column of p-values output from a regression analysis
##' @param ... All other standard R plot arguments available
##' 
##' @export
##' 
##' 
##' @examples
##' 

plot.dsLmFeature <- function(x, ...){
  
  concBand <- function(N, alpha=0.01, ...)
  {
    
    # Calculate 95% confidence intervals. The jth order statistic from a uniform(0,1) 
    # sample has a beta(j,n-j+1) distribution (Casella & Berger, 2002, 2nd edition, pg
    # 230, Duxbury)
    
    conf <- 1-(alpha/2)
    c95<- vector()
    c05<- vector()
    for(i in 1:N){
      c95[i] <- qbeta(conf,i,N-i+1)
      c05[i] <- qbeta(1-conf,i,N-i+1)
    }
    ans <- cbind(c05, c95)
    ans
  }
  
  shade <- function(x1, y1, x2, y2, color = "lightblue") {
    n <- length(x2)
    polygon(c(x1, x2[n:1]), c(y1, y2[n:1]), border = NA, 
            col = color)
  }
  
  qqplotBand <- function(x,  main=NULL, lab, pos,  cex1=1, cex2=1, cex3=1, ...)
  {
    o <- -log10(sort(x,decreasing=F))
    e <- -log10( 1:length(o)/length(o) )
    names(e) <- names(o)
    plot(e,o, main=main, 
         xlab=expression(Expected~~-log[10](italic(p))),
         ylab=expression(Observed~~-log[10](italic(p))),
         xlim=c(0,max(e)), ylim=c(0,max(o)), type="n", 
         cex.lab=cex1, cex.axis=cex2)
    
    N <-length(x)
    band <- concBand(N, ...)
    
    shade(e, -log(band[,1],10), e, -log(band[,2],10))
    
    points(e, o, ...)
    
    lines(e, e, col="blue")
    
  }
  
  pValues = x[, "p-value"]
  
  qqplotBand(pValues)
  
}
