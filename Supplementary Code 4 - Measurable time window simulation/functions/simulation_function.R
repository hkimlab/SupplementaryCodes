HL_extract <- function(data,ID,minimize="RSE_Ft",Ft_upper=0.95,Ft_lower=0.02,time_upper=60,Ft_min=0.85,weight=T){
  
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
  
  sub_data <- data[which(data$ID==ID),]
  
  time <- unique(data$Day)
  
  df <- data.frame(time,x=NA,Ft=NA,nps=NA)
  
  x <- sub_data$Day
  Ft <- sub_data$Ft/100
  nps <- sub_data$nps
  
  for(i in 1:length(x)){
    df[df$time==x[i],2:4] <- c(x[i],Ft[i],nps[i])
  }
  
  df$Ft[df$Ft>1] <- 1
  df$Ft[df$Ft<0] <- 0
  df <- df[-which(df$time==0|df$time>time_upper),]
  
  if(min(df$Ft,na.rm=T)>Ft_min){
    stop("min(Ft) is upper than Ft_min")
  }
  
  RSE_Ft_ind <- which(df$Ft!=1&df$Ft>Ft_lower&df$Ft<Ft_upper)
  RRSE_t_ind <- which(df$Ft!=0&df$Ft>Ft_lower&df$Ft<Ft_upper)
  
  if(weight==F&minimize=="RSE_Ft"){
    weights <- rep(1,nrow(df[RSE_Ft_ind,]))
  }else if(weight==T&minimize=="RSE_Ft"){
    weights <- df$nps[RSE_Ft_ind]
  }else if(weight==F&minimize=="RRSE_t"){
    weights <- rep(1,nrow(df[RRSE_t_ind,]))
  }else if(weight==T&minimize=="RRSE_t"){
    weights <- df$nps[RRSE_t_ind]
  }else{
    stop("Check weight and minimize function")
  }
  
  if(length(RRSE_t_ind)==1|length(RRSE_t_ind)==0|length(RSE_Ft_ind)==1|length(RSE_Ft_ind)==0){
    stop("time points are not enough to estimate lambda and t0")
  }
  t0 <- 1.021
  
  if(minimize=="RRSE_t"){
    suppressWarnings(model <- nls(x~-log(Ft)/lambda+t0,start=list(lambda=0.001),data=df[RRSE_t_ind,],
                                  lower=c(0,0),upper=c(Inf,Inf),weights=weights,algorithm="port",
                                  nls.control(maxiter = 100, tol = 1e-04, minFactor = 1/10^100)))
  }else if(minimize=="RSE_Ft"){
    suppressWarnings(model <- nls(Ft~exp(-lambda*(x-t0)),start=list(lambda=0.001),data=df[RSE_Ft_ind,],
                                  lower=c(0,0),upper=c(Inf,Inf),weights=weights,algorithm="port",
                                  nls.control(maxiter = 100, tol = 1e-04, minFactor = 1/10^100)))
  }else{
    stop("Check minimize function")
  }
  
  lambda <- model$m$getPars()[1]
  
  HL <- log(2)/lambda+t0
  t_hat <- rep(NA,length(time))
  t_hat[time%in%df$time] <- -log(df$Ft)/lambda+t0
  nps <- sub_data$nps
  
  res <- matrix(c(ID,lambda,t0,HL),nrow=1)
  t_hat <- matrix(c(ID,t_hat),nrow=1)
  nps <- matrix(c(ID,nps),nrow=1)
  
  colnames(res) <- c("ID","lambda","t0","HL")
  colnames(t_hat) <- c("ID",paste0("time_",round(time,3)))
  colnames(nps) <- c("ID",paste0("time_",round(time,3)))
  
  row.names(res) <- ID
  return(list(res=res, t_hat=t_hat, nps=nps))
}

simul_graph1 <- function(data, HL_data, HL_lower=0, HL_upper=1000, simulation=F){
  
  # data : data                                                                                                       
  # HL_data : data for half-life(HL)
  # HL_lower : the upper value of the HL
  # HL_upper : the upper value of the HL
  # simulation : whether to use simulated Ft
  
  data <- data[data$ID%in%HL_data$ID,]
  u_id <- unique(data$ID)
  t0 <- 1.021
  HL_data$lambda <- log(2)/(HL_data$HL-t0)
  
  if(simulation==F){
    data <- data[,c("ID","Day","k","nps","Ft")]
  }else{
    data <- data[,c("ID","Day","k")]
    data$nps <- data$Ft <- NA
    
    for(i in 1:length(u_id)){
      Ft <- exp(-HL_data$lambda[HL_data$ID==u_id[i]]*(data$Day[data$ID==u_id[i]]-t0))
      Ft[-1] <- Ft[-1] + rnorm(length(Ft)-1, 0.02739254/100, 1.441617/100)
      Ft[Ft<0] <- 0
      Ft[Ft>1] <- 1
      
      sub_k <- data$k[data$ID==u_id[i]]
      data$nps[data$ID==u_id[i]] <- suppressWarnings(sub_k-rnorm(length(Ft),sub_k*Ft,sqrt(sub_k*Ft*(1-Ft))))
      data$nps[data$ID==u_id[i]&data$Day==0] <- 0
      data$Ft[data$ID==u_id[i]] <- (1-data$nps[data$ID==u_id[i]]/sub_k)*100
    }
    
    data$Ft[which(data$Ft<0)] <- 0
    data$Ft[which(data$Ft>100)] <- 100
  }
  
  
  data$t_hat <- NA
  for(i in 1:length(u_id)){
    temp <- tryCatch(HL_extract(data=data, ID=unique(data$ID)[i]),
                     error=function(e){list(res=c(unique(data$ID)[i],rep(NA,4-1)),
                                            t_hat=c(unique(data$ID)[i],rep(NA,16-1)),
                                            nps=c(unique(data$ID)[i],rep(NA,16-1)))},
                     warning=function(w){list(res=c(unique(data$ID)[i],rep(NA,4-1)),
                                              t_hat=c(unique(data$ID)[i],rep(NA,16-1)),
                                              nps=c(unique(data$ID)[i],rep(NA,16-1)))})
    data$t_hat[data$ID==unique(data$ID)[i]] <- temp$t_hat[-1]
  }
  
  data$RAE <- abs(data$t_hat-data$Day)/data$t_hat*100
  data$Ft_group <- cut(data$Ft, c(0,seq(0.5,99.5,1),100), include=TRUE)
  data$Ft_group_mid <- cut(data$Ft, c(0,seq(0.5,99.5,1),100), include=TRUE, labels=c(0.25,1:99,97.5))
  
  data <- merge(data, HL_data[,c("ID","HL")])
  HL_ind <- which(data$HL>=HL_lower,data$HL<=HL_upper)
  if(length(HL_ind)>0){
    data <- data[HL_ind,]
  }
  
  res <- as.data.frame(data[!is.na(data$Ft_group),] %>% group_by(Ft_group, Ft_group_mid) %>% 
                         summarise(RAE_mean=mean(RAE, na.rm=T), RAE_median=median(RAE, na.rm=T)))
  res$Ft_group_mid <- as.numeric(as.character(res$Ft_group_mid))
  
  return(res=res)
}


simul_graph2 <- function(data, HL_data, simulation=F){
  
  # data : data                                                                                                       
  # HL_data : data for half-life(HL)
  # simulation : whether to use simulated Ft
  
  data <- data[data$ID%in%HL_data$ID,]
  u_id <- unique(data$ID)
  t0 <- 1.021
  HL_data$lambda <- log(2)/(HL_data$HL-t0)
  
  if(simulation==F){
    data <- data[,c("ID","Day","k","nps","Ft")]
  }else{
    data <- data[,c("ID","Day","k")]
    data$nps <- data$Ft <- NA
    
    for(i in 1:length(u_id)){
      Ft <- exp(-HL_data$lambda[HL_data$ID==u_id[i]]*(data$Day[data$ID==u_id[i]]-t0))
      Ft[-1] <- Ft[-1] + rnorm(length(Ft)-1, 0.02739254/100, 1.441617/100)
      Ft[Ft<0] <- 0
      Ft[Ft>1] <- 1
      
      sub_k <- data$k[data$ID==u_id[i]]
      data$nps[data$ID==u_id[i]] <- suppressWarnings(sub_k-rnorm(length(Ft),sub_k*Ft,sqrt(sub_k*Ft*(1-Ft))))
      data$nps[data$ID==u_id[i]&data$Day==0] <- 0
      data$Ft[data$ID==u_id[i]] <- (1-data$nps[data$ID==u_id[i]]/sub_k)*100
    }
    
    data$Ft[which(data$Ft<0)] <- 0
    data$Ft[which(data$Ft>100)] <- 100
  }
  
  for(i in 1:length(u_id)){
    temp <- tryCatch(HL_extract(data=data, ID=unique(data$ID)[i]),
                     error=function(e){list(res=c(unique(data$ID)[i],rep(NA,4-1)),
                                            t_hat=c(unique(data$ID)[i],rep(NA,16-1)),
                                            nps=c(unique(data$ID)[i],rep(NA,16-1)))},
                     warning=function(w){list(res=c(unique(data$ID)[i],rep(NA,4-1)),
                                              t_hat=c(unique(data$ID)[i],rep(NA,16-1)),
                                              nps=c(unique(data$ID)[i],rep(NA,16-1)))})
    if(i==1){
      res <- temp$res
      t_hat <- temp$t_hat
      nps <- temp$nps
    }else{
      res <- rbind(res, temp$res)
      t_hat <- rbind(t_hat, temp$t_hat)
      nps <- rbind(nps, temp$nps)
    }
  }
  
  t_hat <- t_hat[,-1]
  t_hat[which(t_hat==Inf)] <- NA
  nps <- nps[,-1]
  
  mean_total_t_hat <- apply(t_hat,2,mean,na.rm=T)[-1]
  median_total_t_hat <- apply(t_hat,2,median,na.rm=T)[-1]
  
  cut <- apply(t_hat,2,quantile,c(0.25,0.75),na.rm=T)
  t_hat[sweep(t_hat,2,cut[1,],FUN="<")] <- NA; t_hat[sweep(t_hat,2,cut[2,],FUN=">")] <- NA
  IQ_mean_total_t_hat <- apply(t_hat,2,mean,na.rm=T)[-1]
  
  nps[is.na(t_hat)] <- NA
  weight_nps <- sweep(nps,2,apply(nps,2,sum,na.rm=T),FUN="/")
  nps_w_total_t_hat <- apply(t_hat*weight_nps,2,sum,na.rm=T)[-1]
  
  time <- unique(data$Day)[-1]
  RAE_mean <- abs(mean_total_t_hat-time)/mean_total_t_hat*100
  RAE_median <- abs(median_total_t_hat-time)/median_total_t_hat*100
  RAE_IQ_mean <- abs(IQ_mean_total_t_hat-time)/IQ_mean_total_t_hat*100
  RAE_nps_mean <- abs(nps_w_total_t_hat-time)/nps_w_total_t_hat*100
  
  return(data.frame(time, mean=mean_total_t_hat, median=median_total_t_hat, IQ_mean_total_t_hat, nps_mean=nps_w_total_t_hat,
                    RAE_mean, RAE_median, RAE_IQ_mean, RAE_nps_mean))
}

simul_graph3 <- function(data, HL_data, HL_list, simulation=T){
  
  # data : data                                                                                                       
  # HL_data : data for half-life(HL)
  # HL_list : grid for dividing sections of HL
  # simulation : whether to use simulated Ft
  
  data <- data[data$ID%in%HL_data$ID,]
  u_id <- unique(data$ID)
  t0 <- 1.021
  HL_data$lambda <- log(2)/(HL_data$HL-t0)
  
  if(simulation==F){
    data <- data[,c("ID","Day","k","nps","Ft")]
  }else{
    data <- data[,c("ID","Day","k")]
    data$nps <- data$Ft <- NA
    
    for(i in 1:length(u_id)){
      Ft <- exp(-HL_data$lambda[HL_data$ID==u_id[i]]*(data$Day[data$ID==u_id[i]]-t0))
      Ft[-1] <- Ft[-1] + rnorm(length(Ft)-1, 0.02739254/100, 1.441617/100)
      Ft[Ft<0] <- 0
      Ft[Ft>1] <- 1
      
      sub_k <- data$k[data$ID==u_id[i]]
      data$nps[data$ID==u_id[i]] <- suppressWarnings(sub_k-rnorm(length(Ft),sub_k*Ft,sqrt(sub_k*Ft*(1-Ft))))
      data$nps[data$ID==u_id[i]&data$Day==0] <- 0
      data$Ft[data$ID==u_id[i]] <- (1-data$nps[data$ID==u_id[i]]/sub_k)*100
    }
    
    data$Ft[which(data$Ft<0)] <- 0
    data$Ft[which(data$Ft>100)] <- 100
  }
  
  data$t_hat <- NA
  for(i in 1:length(u_id)){
    temp <- tryCatch(HL_extract(data=data, ID=unique(data$ID)[i]),
                     error=function(e){list(res=c(unique(data$ID)[i],rep(NA,4-1)),
                                            t_hat=c(unique(data$ID)[i],rep(NA,16-1)),
                                            nps=c(unique(data$ID)[i],rep(NA,16-1)))},
                     warning=function(w){list(res=c(unique(data$ID)[i],rep(NA,4-1)),
                                              t_hat=c(unique(data$ID)[i],rep(NA,16-1)),
                                              nps=c(unique(data$ID)[i],rep(NA,16-1)))})
    data$t_hat[data$ID==unique(data$ID)[i]] <- temp$t_hat[-1]
  }
  
  data <- merge(data, HL_data[,c("ID","HL")])
  
  data$RAE <- abs(data$t_hat-data$Day)/data$t_hat*100
  data$HL_group <- cut(data$HL, HL_list, include=TRUE)
  data$HL_group_mid <- cut(data$HL, HL_list, include=TRUE, labels=(HL_list[-1]+HL_list[-length(HL_list)])/2)
  data <- data[!is.na(data$HL_group)&data$Day!=0,]
  data$Day <- as.factor(data$Day)
  res <- as.data.frame(data %>% group_by(Day, HL_group) %>% summarise(RAE=median(RAE, na.rm=T)))
  return(res)
}

