gomperz <- function(data, ID, Ft_upper=1, Ft_lower=0, time_upper=60, Ft_min=1, weight=F){
  
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
  
  old_c <- 0.1; old_b <- 1
  diff <- 1
  
  while(diff>10^-3){
    suppressWarnings(gom_m_temp <- nls(Ft ~ exp(-old_b*exp(c*x)+old_b),start=list(c=old_c), data=df,
                                       lower=c(0), upper=c(Inf), weights=weights, algorithm="port",
                                       nls.control(maxiter=100, tol=1e-04, minFactor = 1/10^100, warnOnly=T)))
    new_c <- as.vector(gom_m_temp$m$getPars())
    suppressWarnings(gom_m <- nls(Ft ~ exp(-b*exp(new_c*x)+b), start=list(b=old_b), data=df,
                                  lower=c(0), upper=c(Inf), weights=weights, algorithm="port",
                                  nls.control(maxiter=100, tol = 1e-04, minFactor = 1/10^100, warnOnly=T)))
    new_b <- as.vector(gom_m$m$getPars())
    diff <- gom_m_temp$m$deviance()-gom_m$m$deviance()
    old_c <- new_c;old_b <- new_b
  }
  
  gom_b <- new_b
  gom_c <- new_c
  
  F_hat <- exp(-gom_b*exp(gom_c*df$x)+gom_b)
  
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
  t_hat <- log(1-log(df$Ft)/gom_b)/gom_c
  
  RRSE_lin <- sqrt(sum(((df$time-t_hat)/df$time)^2, na.rm=T)/(length(na.omit(t_hat))-1))
  RRAE_lin <- sum((abs(df$time-t_hat)/df$time), na.rm=T)/(length(na.omit(t_hat))-1)
  
  param <- data.frame(gom_b, gom_c)
  error <- data.frame(RSS_gom=RSS, MSE_gom=MSE, MAE_gom=MAE, w_MSE_gom=w_MSE, w_MAE_gom=w_MAE)
  AIC_BIC <- data.frame(AIC_gom=AIC, BIC_gom=BIC, w_AIC_gom=w_AIC, w_BIC_gom=w_BIC)
  
  return(list(param=param,error=error,AIC_BIC=AIC_BIC))
}