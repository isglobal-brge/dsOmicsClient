# a <- merge(results$server1, reference, by.x = "rs", by.y = "Rs")
# ggplot(a) + geom_point(aes(x = freq, y = Freq))
# 
# reference3 <- cbind(reference2, str_split_fixed(reference2$ChrPosID, ":", 2))
# 
# b <- merge(results$server1, reference3, by.x = "pos", by.y = "2")
# ggplot(b) + geom_point(aes(x = freq, y = Freq1))
# 
# 
# 
# a <- merge(results$cohort3, reference, by.x = "rs", by.y = "Rs")
# ggplot(a) + geom_point(aes(x = 1-freq, y = Freq))


eafPlot <- function(x, reference){
  
  # Check columns do exist
  
  # Single cohort
  if(inherits(x, "data.frame")){
    merge(x, reference, by.x = "rs", by.y = "Rs")
    return(.eafPlot(x))
  }
  # Multi-cohort
  if(inherits(x, "list")){
    eafPlots <- lapply(names(x), function(y){
      .eafPlot(x[[y]])
    })
    return(gridExtra::grid.arrange(grobs = eafPlots, ncol = 2))
  }
  
}

.eafPlot <- function(){
  
}