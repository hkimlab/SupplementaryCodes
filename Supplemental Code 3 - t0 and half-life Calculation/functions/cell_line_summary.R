cell_line_summary <- function(data, ID, minimize="RSE_Ft", Ft_upper=0.95, Ft_lower=0.02, time_upper=60, Ft_min=0.85, t0=1.021,
                              use_HL=F, HL, weight=T, plot=F, data_name, ...){
  
  # data : data                                                                                                       
  # minimize : object function to minimize
  # Ft_upper : the upper value of the Ft 
  # Ft_lower : the lower value of the Ft                               
  # time_upper : the upper value of the time                       
  # Ft_min : the lower vaule of min(Ft), if min(Ft) is higher than Ft_min, error message appears
  # t0 : waiting time
  # use_HL : whether to use certain helf life
  # HL : helf life
  # weight : whether to use nps weight
  # plot : whether to display the graph
  # data_name : used data name       
  
  df <- extract_cell_line(data=data, ID=ID, Ft_upper=Ft_upper, Ft_lower=Ft_lower, time_upper=time_upper, Ft_min=Ft_min)
  
  RSE_Ft_ind <- which(df$Ft!=1&df$Ft>Ft_lower&df$Ft<Ft_upper)
  RRSE_t_ind <- which(df$Ft!=0&df$Ft>Ft_lower&df$Ft<Ft_upper)
  
  if(weight==F&minimize=="RSE_Ft"){
    
    weights <- rep(1, nrow(df[RSE_Ft_ind,]))
    
  }else if(weight==T&minimize=="RSE_Ft"){
    
    weights <- df$nps[RSE_Ft_ind]
    
  }else if(weight==F&minimize=="RRSE_t"){
    
    weights <- rep(1, nrow(df[RRSE_t_ind,]))
    
  }else if(weight==T&minimize=="RRSE_t"){
    
    weights <- df$nps[RRSE_t_ind]
    
  }else{
    
    stop("Check weight and minimize function")
    
  }
  
  if(length(RRSE_t_ind)==1|length(RRSE_t_ind)==0|length(RSE_Ft_ind)==1|length(RSE_Ft_ind)==0){
    
    stop("time points are not enough to estimate lambda and t0")
    
  }
  
  if(use_HL==F){
    
    if(minimize=="RRSE_t"){
      
      suppressWarnings(model <- nls(x~-log(Ft)/lambda+t0, start=list(lambda=0.001), data=df[RRSE_t_ind,],
                                    lower=c(0,0), upper=c(Inf,Inf), weights=weights, algorithm="port",
                                    nls.control(maxiter=100, tol=1e-04, minFactor = 1/10^100)))
      
    }else if(minimize=="RSE_Ft"){
      
      suppressWarnings(model <- nls(Ft~exp(-lambda*(x-t0)), start=list(lambda=0.001), data=df[RSE_Ft_ind,],
                                    lower=c(0,0), upper=c(Inf,Inf), weights=weights, algorithm="port",
                                    nls.control(maxiter=100, tol=1e-04, minFactor=1/10^100)))
      
    }else{
      
      stop("Check minimize function")
      
    }
    
    lambda <- model$m$getPars()[1]
    
  }else{
    
    lambda <- log(2)/(HL[HL$ID==ID,2]-t0)
    
  }
  
  half_life <- log(2)/lambda+t0
  
  Ft_hat <- exp(-lambda*(df$x-t0))
  RSE_Ft <- sqrt(sum((df$Ft-Ft_hat)^2, na.rm=T)/(length(na.omit(Ft_hat))-1))
  
  plot_x <- df$x
  plot_y <- df$Ft
  
  df$Ft[df$Ft==0] <- NA
  
  t_hat <- -log(df$Ft)/lambda+t0
  RRSE_t <- sqrt(sum(((df$time-t_hat)/df$time)^2, na.rm=T)/(length(na.omit(t_hat))-1))
  RAE_t <- abs(t_hat-df$time)/df$time
  
  if(plot==T){
    
    plot(plot_x, plot_y, ylim=c(0,1), xlim=c(0, time_upper), xlab="Time", ylab="Ft", pch=19,
        main=paste(data_name, " / ID : ", ID, sep=""))
    curve(0*x+1, 0, t0, add=T, lwd=3, lty=2, col="red")
    curve(exp(-lambda*(x-t0)), t0, time_upper, add=T, lwd=3, lty=2, col="red")
    
  }
  
  res <- matrix(c(ID, lambda, t0, half_life, RSE_Ft, RRSE_t, Ft_hat, t_hat, RAE_t), nrow=1)
  colnames(res) <- c("ID", "lambda", "t0", "half_life", "RSE_Ft", "RRSE_t", paste("Ft_hat", round(df$time, 1), sep="_"),
                     paste("t_hat", round(df$time, 1), sep="_"), paste("RAE_t", round(df$time, 1), sep="_"))
  row.names(res) <- ID
  
  return(res)
}
