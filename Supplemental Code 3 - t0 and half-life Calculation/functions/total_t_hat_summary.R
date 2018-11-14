total_t_hat_summary <- function(data, minimize="RSE_Ft", Ft_upper=0.95, Ft_lower=0.02, time_upper=60, Ft_min=0.85, t0=1.021,
                                HL, weight=T, ...){
  
  # data : data                                                                                                       
  # ID : cell line ID      
  # minimize : object function to minimize
  # Ft_upper : the upper value of the Ft 
  # Ft_lower : the lower value of the Ft                               
  # time_upper : the upper value of the time                       
  # Ft_min : the lower vaule of min(Ft), if min(Ft) is higher than Ft_min, error message appears
  # t0 : waiting time
  # HL : helf life
  # weight : whether to use nps weight
  
  pb <- txtProgressBar(max=length(unique(data$ID)), style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  
  len <- length(cell_line_summary(data=data, ID=unique(data$ID)[1], minimize=minimize, Ft_upper=Ft_upper, Ft_lower=Ft_lower, 
                                  time_upper=time_upper, Ft_min=Ft_min, t0=t0, use_HL=T, HL=HL, weight=weight, plot=F))
  
  HL_summary <- foreach(i=1:length(unique(data$ID)), .export=c("cell_line_summary","extract_cell_line"),
                        .options.snow=opts, .combine='rbind') %dopar% {
    
    tryCatch(cell_line_summary(data=data, ID=unique(data$ID)[i], minimize=minimize, Ft_upper=Ft_upper, Ft_lower=Ft_lower, 
                               time_upper=time_upper, Ft_min=Ft_min, t0=t0, use_HL=T, HL=HL, weight=weight),
             error=function(e){c(unique(data$ID)[i],rep(NA,len-1))},
             warning=function(w){c(unique(data$ID)[i],rep(NA,len-1))})
  
  }
 
  ind <- grep("t_hat",colnames(HL_summary))[grep("t_hat",colnames(HL_summary))!=grep("Ft_hat",colnames(HL_summary))]
  t_hat <- HL_summary[,ind]

  ## extract nps
  
  len <- sum(unique(data$Day)<60)
  
  nps_mat <- foreach(i=1:length(unique(data$ID)), .options.snow=opts, .combine='rbind', .export="extract_cell_line") %dopar% {
    
    tryCatch(extract_cell_line(data=data, ID=unique(data$ID)[i], Ft_upper=Ft_upper, Ft_lower=Ft_lower, 
                               time_upper=time_upper, Ft_min=Ft_min)$nps,
             error=function(e){rep(NA,len-1)})
    
  }
  
  cut <- apply(t_hat, 2, quantile, c(0.25,0.75), na.rm=T)
  t_hat[sweep(t_hat, 2, cut[1,], FUN="<")] <- NA; t_hat[sweep(t_hat, 2, cut[2,], FUN=">")] <- NA
  nps_mat[is.na(t_hat)] <- NA
  
  weight_nps <- sweep(nps_mat, 2, apply(nps_mat, 2, sum, na.rm=T), FUN="/")
  nps_w_total_t_hat <- apply(t_hat*weight_nps, 2, sum, na.rm=T)
  
  return(nps_w_total_t_hat)
}
  
