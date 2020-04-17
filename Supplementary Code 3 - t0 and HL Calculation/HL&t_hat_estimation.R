## data import & library

setwd("PathtoWorkDir") 	# set working directory 
library(doSNOW) 		# install.packages("doSNOW")
library(tcltk) 			# install.packages("tcltk")
source("function/sourceDir.R")
sourceDir("function/", trace=F)

data 				 <- read.csv("SampleData/TestInput.csv") # data import
data 				 <- data[data$TP_measured>0,]
data$Ft[data$k==0] 	 <- NA
data$Ft[data$mutK<0] <- NA

## parallel computation setting

cl <- makeSOCKcluster(3) 	# input : number of cores
registerDoSNOW(cl) 			# stopCluster(cl) 

## cell line summary

pb 		 <- txtProgressBar(max=length(unique(data$ID)), style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts 	 <- list(progress=progress)

len 	 <- length(cell_line_summary(data=data, ID=unique(data$ID)[1], minimize="RSE_Ft", Ft_upper=0.95, Ft_lower=0.02, 
                         time_upper=60, Ft_min=0.85, t0=1.021, use_HL=F, weight=T, plot=F))

result 	 <- foreach(i=1:length(unique(data$ID)), .options.snow=opts, .combine='rbind') %dopar% {
  
  tryCatch(cell_line_summary(data=data, ID=unique(data$ID)[i], minimize="RSE_Ft", Ft_upper=0.95, Ft_lower=0.02, 
                             time_upper=60, Ft_min=0.85, t0=1.021, use_HL=F, weight=T, plot=F),
           error=function(e){c(unique(data$ID)[i],rep(NA,len-1))},
           warning=function(w){c(unique(data$ID)[i],rep(NA,len-1))})
  
}

write.csv(as.data.frame(result),"SampleData/TestOutput.csv",row.names=F) # save result
