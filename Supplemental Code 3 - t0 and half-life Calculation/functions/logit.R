logit <- function(data, ID, Ft_upper=1, Ft_lower=0, time_upper=60, Ft_min=1, weight=F){
  
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
  
  logit_m_base <- suppressWarnings(nls(Ft ~ (1+exp(-0.5*x0))/(1+exp(0.5*(x-x0))), start=list(x0=0), data=df,
                                       lower=c(-Inf), upper=c(Inf), weights=weights, algorithm="port",
                                       nls.control(maxiter = 100, tol = 1e-04, minFactor = 1/10^10,warnOnly =T)))
  
  old_x0 <- logit_m_base$m$getPars(); old_k <- 0
  diff <- 1; j <- 1
  
  suppressWarnings(while(diff>10^-3&j<100){
    
    try(logit_m_temp <- nls(Ft ~ (1+exp(-k*old_x0))/(1+exp(k*(x-old_x0))), start=list(k=old_k), data=df,
                            lower=c(0), upper=c(Inf), weights=weights, algorithm="port",
                            nls.control(maxiter = 100, tol = 1e-04, minFactor = 1/10^10, warnOnly =T)),silent=T)
    new_k <- as.vector(logit_m_temp$m$getPars())
    
    try(logit_m <- nls(Ft ~ (1+exp(-new_k*x0))/(1+exp(new_k*(x-x0))),start=list(x0=old_x0), data=df,
                       lower=c(-Inf), upper=c(Inf), weights=weights, algorithm="port",
                       nls.control(maxiter = 100, tol = 1e-04, minFactor = 1/10^10, warnOnly =T)),silent=T)
    new_x0 <- as.vector(logit_m$m$getPars())
    diff <- logit_m_temp$m$deviance()-logit_m$m$deviance()
    old_x0 <- new_x0;old_k <- new_k
    
    j <- j + 1
  })
  
  logit_k <- new_k
  logit_x0 <- new_x0
  
  F_hat <- (1+exp(-logit_k*logit_x0))/(1+exp(logit_k*(df$x-logit_x0)))
  
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
  t_hat <- log((1+exp(-logit_k*logit_x0)/df$Ft)-1)/logit_k+logit_x0
  
  RRSE_lin <- sqrt(sum(((df$time-t_hat)/df$time)^2, na.rm=T)/(length(na.omit(t_hat))-1))
  RRAE_lin <- sum((abs(df$time-t_hat)/df$time), na.rm=T)/(length(na.omit(t_hat))-1)
  
  param <- data.frame(logit_k,logit_x0)
  error <- data.frame(RSS_logit=RSS, MSE_logit=MSE, MAE_logit=MAE, w_MSE_logit=w_MSE, w_MAE_logit=w_MAE)
  AIC_BIC <- data.frame(AIC_logit=AIC, BIC_logit=BIC, w_AIC_logit=w_AIC, w_BIC_logit=w_BIC)
  
  return(list(param=param,error=error,AIC_BIC=AIC_BIC))
}