linear <- function(data, ID, Ft_upper=1, Ft_lower=0, time_upper=60, Ft_min=1, weight=F){
  
  # data : data                                                                                                       
  # ID : cell line ID                                                                                                 
  # Ft_upper : the upper value of the Ft 
  # Ft_lower : the lower value of the Ft                               
  # time_upper : the upper value of the time                       
  # Ft_min : the lower vaule of min(Ft), if min(Ft) is higher than Ft_min, error message appears
  # weight : whether to use nps weight
  
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
  
  lin_m <- lm(Ft~x+0, data=df, weights=weights, offset=rep(1,nrow(df)))
  beta <- lin_m$coefficients
  
  F_hat <- 1 + df$x*beta
  
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
  t_hat <- (df$Ft - 1)/beta
  
  RRSE_lin <- sqrt(sum(((df$time-t_hat)/df$time)^2, na.rm=T)/(length(na.omit(t_hat))-1))
  RRAE_lin <- sum((abs(df$time-t_hat)/df$time), na.rm=T)/(length(na.omit(t_hat))-1)
  
  param <- data.frame(beta)
  error <- data.frame(RSS_lin=RSS, MSE_lin=MSE, MAE_lin=MAE, w_MSE_lin=w_MSE, w_MAE_lin=w_MAE)
  AIC_BIC <- data.frame(AIC_lin=AIC, BIC_lin=BIC, w_AIC_lin=w_AIC, w_BIC_lin=w_BIC)
  
  return(list(param=param,error=error,AIC_BIC=AIC_BIC))
}