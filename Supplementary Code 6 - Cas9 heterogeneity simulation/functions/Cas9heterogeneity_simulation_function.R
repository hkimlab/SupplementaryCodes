HL_extract <- function(data,ID,minimize="RSE_Ft",Ft_upper=0.95,Ft_lower=0.02,time_upper=60,Ft_min=0.85,weight=T){
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


time_predict <- function(data, HL_data, heterogeneity=NA, t0=1.021, Ideal=T, k_scale=1, calc_type, ...){

  HL_data$lambda <- log(2)/(HL_data$HL-t0)
  
  data <- data[data$ID%in%HL_data$ID,]
  u_id <- unique(data$ID)
  
  HL_data <- HL_data[HL_data$ID%in%u_id,]
    
  data <- data[,c("ID","Day","k")]
  data$nps <- data$Ft <- NA
  
  HL_data$s_lambda <- NA

  for(i in 1:length(u_id)){
    HL_list <- c()
    HL <- HL_data$HL[HL_data$ID==u_id[i]]
    temp <- data[data$ID==u_id[i],]
    t_list <- temp$Day
    t_list <- t_list[t_list!=0]
    for(j in 1:length(t_list)){
      k <- temp$k[temp$Day==t_list[j]]
      k <- round(k*k_scale + rnorm(1),0)
      k <- ifelse(k<0,0,k)
      if(Ideal==F){
        s_HL <- rgamma(k,HL/heterogeneity,1/heterogeneity)
        s_lambda <- log(2)/(s_HL-t0)
      }else{
        s_HL <- HL
        s_lambda <- log(2)/(s_HL-t0)
      }
      HL_list <- c(HL_list,s_HL)
      Ft <- exp(-s_lambda*(t_list[j]-t0))
      Ft <- mean(Ft,na.rm=T)
      # diff mean 0.02739254, diff sd 1.441617
      Ft <- Ft + rnorm(1, 0.02739254/100, 1.441617/100)
      Ft[Ft<0] <- 0; Ft[Ft>1] <- 1
      temp$k[j+1] <- k
      temp$nps[j+1] <- sum(rbinom(k,1,1-Ft))
      temp$Ft[j+1] <- (temp$k[j+1]-temp$nps[j+1])/temp$k[j+1]*100
    }
    temp[temp$Day==0,c("Ft","nps")] <- c(100,0)
    data[data$ID==u_id[i],] <- temp
    HL_data$s_HL[HL_data$ID==u_id[i]] <- mean(HL_list)
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
  
  if(calc_type=="mean"){
    RAE <- apply(calc_RAE(t_true=unique(data$Day), t_hat=t_hat), 1, mean, na.rm=T)
  }else{
    RAE <- apply(calc_RAE(t_true=unique(data$Day), t_hat=t_hat), 1, median, na.rm=T)
  }
  return(data.frame(HL=HL_data$s_HL,RAE=RAE))
}

calc_RAE <- function(t_true, t_hat){
  RAE <- t_hat
  for(i in 1:length(t_true)){
    RAE[,i] <- abs(t_hat[,i]-t_true[i])/t_true[i]*100
  }
  return(RAE[,-1])
}

HL_histogram <- function(s_HL, HL_data, ID){
  id_ind <- which(HL_data$ID==ID)
  HL <- HL_data$HL[id_ind]
  hist(s_HL[,id_ind],main=paste0("ID=",ID,", HL=",round(HL,3)), xlab="Simulated HL")
  abline(v=HL, col="red")
}

comb <- function(x, ...) {  
  mapply(rbind,x,...,SIMPLIFY=FALSE)
}

combine_RAE <- function(data, RAE, heterogeneity){
  for(i in 1:ncol(RAE)){
    temp <- data.frame(heterogeneity=heterogeneity,time=round(unique(data$Day)[i+1],1),RAE=RAE[,i])
    if(i==1){
      RAE_data <- temp
    }else{
      RAE_data <- rbind(RAE_data, temp)
    }
  }
  return(RAE_data)
}

