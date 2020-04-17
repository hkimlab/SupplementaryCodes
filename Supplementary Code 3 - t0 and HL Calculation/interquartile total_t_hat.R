## data import & libaray

setwd("PathtoWorkDir") 	# set working directory 
library(doSNOW) 		# install.packages("doSNOW")
library(tcltk) 			# install.packages("tcltk")
source("function/sourceDir.R")
sourceDir("function/", trace=F)

data 				 <- read.csv("SampleData/TestInput.csv") # data import
data 				 <- data[data$TP_measured>0,]
data$Ft[data$k==0] 	 <- NA
data$Ft[data$mutK<0] <- NA

HL <- read.csv("SampleData/TestInput.csv")

## parallel computation setting

cl <- makeSOCKcluster(3) # input : number of cores
registerDoSNOW(cl) 		 # stopCluster(cl)

## model comparison result
nps_w_total_t_hat <- total_t_hat_summary(data=data, minimize="RSE_Ft", Ft_upper=0.95, Ft_lower=0.02, time_upper=60, Ft_min=0.85, t0=1.021,
                                   HL=HL[,c(1,2)], weight=T)

write.csv(as.data.frame(nps_w_total_t_hat),"SampleData/TestOutput.csv",row.names=T)
