## data import & library

setwd("PathToWorkDir") 		# set working directory 
library(doSNOW) 			# install.packages("doSNOW")
library(tcltk) 				# install.packages("tcltk")
source("function/sourceDir.R")
sourceDir("function/", trace=F)

data 				 <- read.csv("SampleData/TestInput.csv")
data 				 <- data[data$TP_measured>0,]
data$Ft[data$k==0] 	 <- NA
data$Ft[data$mutK<0] <- NA

## parallel computation setting

cl <- makeSOCKcluster(3) 	# input : number of cores
registerDoSNOW(cl)			# stopCluster(cl)

## model comparison result

pb 			<- txtProgressBar(max=length(unique(data$ID)), style=3)
progress 	<- function(n) setTxtProgressBar(pb, n)
opts 		<- list(progress=progress)

result <- foreach(i=1:length(unique(data$ID)), .options.snow=opts, .combine='comb', .multicombine=T) %dopar% {
  
  tryCatch(suppressWarnings(model_comparison(data=data, ID=unique(data$ID)[i], Ft_upper=1, Ft_lower=0, 
                                          time_upper=60, Ft_min=1, weight=T, est_t0=T)),
           error=function(e){list(param=c(unique(data$ID)[i], rep(NA,7)), error=c(unique(data$ID)[i], rep(NA,28)),
                                  AIC_BIC=c(unique(data$ID)[i] ,rep(NA,16)))})

}

write.csv(result$param,"parameter_summary.csv",row.names=F)
write.csv(result$error,"error_summary.csv",row.names=F)
write.csv(result$AIC_BIC,"AIC_BIC_summary.csv",row.names=F)

## extract nps

len 	<- sum(unique(data$Day)<60)
nps_mat <- foreach(i=1:length(unique(data$ID)), .options.snow=opts, .combine='rbind', .multicombine=T) %dopar% {
  
  tryCatch(extract_cell_line(data=data, ID=unique(data$ID)[i], Ft_upper=1, Ft_lower=0, time_upper=60, Ft_min=1)$nps,
           error=function(e){rep(NA,len-1)})
  
}

## calculate total AIC, BIC

nps_weight 		<- apply(nps_mat,1,sum,na.rm=T)/sum(nps_mat,na.rm=T)
total_AIC_BIC 	<- apply(result$AIC_BIC[,-1],2,mean,na.rm=T)
w_total_AIC_BIC <- apply(result$AIC_BIC[,-1]*matrix(rep(nps_weight,16),nrow=length(unique(data$ID))),2,sum,na.rm=T)

write.csv(as.data.frame(total_AIC_BIC),"total_AIC_BIC_summary.csv",row.names=T)
write.csv(as.data.frame(w_total_AIC_BIC),"w_total_AIC_BIC_summary.csv",row.names=T)


