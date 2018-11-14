model_comparison <- function(data, ID, Ft_upper=1, Ft_lower=0, time_upper=60, Ft_min=1, weight=F, est_t0=T){
  
  # data : data
  # ID : cell line ID            
  # Ft_upper : the upper value of the Ft 
  # Ft_lower : the lower value of the Ft                               
  # time_upper : the upper value of the time                       
  # Ft_min : the lower vaule of min(Ft), if min(Ft) is higher than Ft_min, error message appears
  # weight : whether to use nps weight
  # est_t0 : whether to estimate t0 for exponential model
  
  lin_res <- linear(data, ID, Ft_upper=Ft_upper, Ft_lower=Ft_lower, time_upper=time_upper, Ft_min=Ft_min, weight=Ft_min)
  exp_res <- exponential(data, ID, Ft_upper=Ft_upper, Ft_lower=Ft_lower, time_upper=time_upper, Ft_min=Ft_min, weight=Ft_min, est_t0=est_t0)
  gom_res <- gomperz(data, ID, Ft_upper=Ft_upper, Ft_lower=Ft_lower, time_upper=time_upper, Ft_min=Ft_min, weight=Ft_min)
  logit_res <- logit(data, ID, Ft_upper=Ft_upper, Ft_lower=Ft_lower, time_upper=time_upper, Ft_min=Ft_min, weight=Ft_min)
  
  param <- data.frame(ID,lin_res$param,exp_res$param,gom_res$param,logit_res$param)
  error <- data.frame(ID,lin_res$error,exp_res$error,gom_res$error,logit_res$error)
  
  AIC_BIC <- data.frame(ID,lin_res$AIC_BIC,exp_res$AIC_BIC,gom_res$AIC_BIC,logit_res$AIC_BIC)
  
  row.names(param) <- row.names(error) <- row.names(AIC_BIC) <- ID
  
  return(list(param=param,error=error,AIC_BIC=AIC_BIC))
}