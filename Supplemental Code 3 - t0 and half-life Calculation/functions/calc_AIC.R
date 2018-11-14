calc_AIC <- function(MSE, n, p){
  
  n + n * log(2 * pi) + n * log(MSE) + 2 * p
  
}