comb <- function(x, ...) {  
  
  mapply(rbind,x,...,SIMPLIFY=FALSE)
  
}