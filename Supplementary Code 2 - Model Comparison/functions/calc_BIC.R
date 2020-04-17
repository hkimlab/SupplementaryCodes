calc_BIC <- function(MSE, n, p){
  
  n + n * log(2 * pi) + n * log(MSE) + log(n) * p
  
}
