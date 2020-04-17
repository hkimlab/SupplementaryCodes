HL_extract <- function(data,ID,minimize="RSE_Ft",Ft_upper=0.95,Ft_lower=0.02,time_upper=60,Ft_min=0.85,weight=T,
                       use_HL=F,HL_data=NA){
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
  
  if(use_HL==F){
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
    
  }else{
    HL <- HL_data$HL[HL_data$ID==ID]
    lambda <- log(2)/(HL-t0) 
  }
  
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

Ft_simul <- function(data, HL_data, t0=1.021, total_k="origin"){
  
  HL_data$lambda <- log(2)/(HL_data$HL-t0)
  
  data <- data[data$ID%in%HL_data$ID,]
  u_id <- unique(data$ID)
  
  data <- data[,c("ID","Day","k")]
  data$nps <- data$Ft <- NA
  
  tp <- unique(data$Day)
  for(t in 1:length(tp)){
    if(is.numeric(total_k)){
      temp_k <- data$k[data$Day==tp[t]]
      data$k[data$Day==tp[t]] <- as.vector(rmultinom(1,total_k,temp_k*total_k/sum(temp_k)))
    }
  }
  
  HL_data$s_lambda <- NA
  
  for(i in 1:length(u_id)){
    HL_list <- c()
    HL <- HL_data$HL[HL_data$ID==u_id[i]]
    temp <- data[data$ID==u_id[i],]
    t_list <- temp$Day
    t_list <- t_list[t_list!=0]
    for(j in 1:length(t_list)){
      k <- temp$k[temp$Day==t_list[j]]
      s_HL <- HL
      s_lambda <- log(2)/(s_HL-t0)
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
  return(list(data=data, HL_data=HL_data))
}

time_predict <- function(iter, data, HL_data, test_data, test_HL_data, t0=1.021, HL_k="origin", Ft_k="origin",
                         calc_type, ...){
  
  origin_data <- Ft_simul(data=data, HL_data=HL_data, total_k=HL_k)
  data <- origin_data$data
  
  data <- data[data$Day==sample(unique(data$Day)[-1],1),]
  
  u_id <- unique(data$ID)
  new_HL_data <- data.frame(ID=u_id, HL=NA)
  
  for(i in 1:length(u_id)){
    temp <- data[data$ID==u_id[i],]
    lambda <- log(temp$Ft/100)/(t0-temp$Day)
    new_HL_data$HL[i] <- log(2)/lambda+t0
  }
  new_HL_data$HL[is.infinite(new_HL_data$HL)] <- NA
  
  new_data <- Ft_simul(data=test_data, HL_data=test_HL_data, total_k=Ft_k)
  test_data <- new_data$data
  
  t_len <- length(unique(test_data$Day))
  
  u_id <- unique(test_data$ID)
  for(i in 1:length(u_id)){
    temp <- tryCatch(
      HL_extract(data=test_data, ID=unique(test_data$ID)[i], HL_data=new_HL_data, use_HL=T),
      error=function(e){list(res=c(unique(test_data$ID)[i],rep(NA,4-1)),
                             t_hat=c(unique(test_data$ID)[i],rep(NA,t_len)),
                             nps=c(unique(test_data$ID)[i],rep(NA,t_len)))},
      warning=function(w){list(res=c(unique(test_data$ID)[i],rep(NA,4-1)),
                               t_hat=c(unique(test_data$ID)[i],rep(NA,t_len)),
                               nps=c(unique(test_data$ID)[i],rep(NA,t_len)))})
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
  
  cut <- apply(t_hat,2,quantile,c(0.25,0.75),na.rm=T)
  t_hat[sweep(t_hat,2,cut[1,],FUN="<")] <- NA; t_hat[sweep(t_hat,2,cut[2,],FUN=">")] <- NA
  nps[is.na(t_hat)] <- NA
  weight_nps <- sweep(nps,2,apply(nps,2,sum,na.rm=T),FUN="/")
  nps_w_total_t_hat <- apply(t_hat*weight_nps,2,sum,na.rm=T)
  
  if(calc_type=="mean"){
    RAE <- mean(calc_RAE(t_true=unique(test_data$Day), t_hat=t(nps_w_total_t_hat)), na.rm=T)
  }else{
    RAE <- median(calc_RAE(t_true=unique(test_data$Day), t_hat=t(nps_w_total_t_hat)), na.rm=T)
  }
  return(data.frame(iter=iter, HL_k=HL_k, Ft_k=Ft_k, RAE=RAE))
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

