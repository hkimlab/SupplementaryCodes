exponential <- function(data, ID, Ft_upper=1, Ft_lower=0, time_upper=60, Ft_min=1, weight=F, est_t0=T){
  
  # data : data                                                                                                       
  # ID : cell line ID                                                                                                 
  # Ft_upper : the upper value of the Ft 
  # Ft_lower : the lower value of the Ft                               
  # time_upper : the upper value of the time                       
  # Ft_min : the lower vaule of min(Ft), if min(Ft) is higher than Ft_min, error message appears
  # weight : whether to use nps weight
  # est_t0 : whether to estimate t0
  
  df <- extract_cell_line(data=data, ID=ID, Ft_upper=Ft_upper, Ft_lower=Ft_lower, time_upper=time_upper, Ft_min=Ft_min)
  
  if(weight==F){
    weights <- rep(1,nrow(df))
  }else if(weight==T){
    weights <- df$nps
  }else{
    stop("Check weight")
  }
  
  if(nrow(df)%in%c(0,1)){
    stop("time points are not enough to estimate parameter")
  }
  
  if(est_t0==F){
    
    exp_t0 <- 0
    suppressWarnings(exp_m <- nls(Ft~exp(-lambda*(x-exp_t0)), start=list(lambda=0.001), data=df,
                                  lower=c(0), upper=c(Inf), weights=weights, algorithm="port",
                                  nls.control(maxiter=100, tol=1e-04, minFactor=1/10^100,warnOnly =T)))
    exp_lambda <- exp_m$m$getPars()
    
  }else{
    
    suppressWarnings(exp_m <- nls(Ft~exp(-lambda*(x-exp_t0)), start=list(lambda=0.001,exp_t0=1), data=df,
                                  lower=c(0,0), upper=c(Inf,Inf), weights=weights, algorithm="port",
                                  nls.control(maxiter=100, tol=1e-04, minFactor=1/10^100,warnOnly =T)))
    
    exp_lambda <- exp_m$m$getPars()[1]
    exp_t0 <- exp_m$m$getPars()[2]
  }
  
  
  F_hat <- exp(-exp_lambda*(df$x-exp_t0))
  
  n <- sum(is.na(df$Ft)==F)
  nps_w <- df$nps/sum(df$nps, na.rm=T)
  
  RSS <- sum((df$Ft-F_hat)^2, na.rm=T)
  MSE <- RSS/n
  MAE <- mean(abs(df$Ft-F_hat), na.rm=T)
  w_MSE <- sum((df$Ft-F_hat)^2*nps_w, na.rm=T)
  w_MAE <- sum(abs(df$Ft-F_hat)*nps_w, na.rm=T)
  
  AIC <- calc_AIC(MSE, n, 1)
  BIC <- calc_BIC(MSE, n, 1)
  w_AIC <- calc_AIC(w_MSE, n, 1)
  w_BIC <- calc_BIC(w_MSE, n, 1)
  
  df$Ft[df$Ft==0] <- NA
  t_hat <- -log(df$Ft)/exp_lambda+exp_t0
  
  RRSE_lin <- sqrt(sum(((df$time-t_hat)/df$time)^2, na.rm=T)/(length(na.omit(t_hat))-1))
  RRAE_lin <- sum((abs(df$time-t_hat)/df$time), na.rm=T)/(length(na.omit(t_hat))-1)
  
  param <- data.frame(exp_t0, exp_lambda)
  error <- data.frame(RSS_exp=RSS, MSE_exp=MSE, MAE_exp=MAE, w_MSE_exp=w_MSE, w_MAE_exp=w_MAE)
  AIC_BIC <- data.frame(AIC_exp=AIC, BIC_exp=BIC, w_AIC_exp=w_AIC, w_BIC_exp=w_BIC)
  
  return(list(param=param, error=error, AIC_BIC=AIC_BIC))
}