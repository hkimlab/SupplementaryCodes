## data import & library

# set working directory 
setwd("PathtoWorkDir") 
source("simulation_function.R")
library(openxlsx)
library(dplyr)
library(doSNOW)
library(ggplot2)

## parallel computation setting
cl <- makeSOCKcluster(6) # input : number of cores
registerDoSNOW(cl)

## Import data
data <- read.csv("TestInput.csv")

HL_data <- read.csv("HL.csv")
HL_data <- na.omit(HL_data)
names(HL_data) <- c("ID","HL")


## Graph 1
# original data
graph_1 <- simul_graph1(data, HL_data, HL_lower=0, HL_upper=1000, simulation=F)
graph_1

# simulation data (n=300)
iter <- 300
pb <- txtProgressBar(max=iter, style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

graph_1_simul <- foreach(i=1:iter, .options.snow=opts, .combine='rbind', .packages=c('dplyr')) %dopar% {
  simul_graph1(data, HL_data, HL_lower=0, HL_upper=1000, simulation=T)
}

graph_1_simul <- as.data.frame(graph_1_simul %>% group_by(Ft_group, Ft_group_mid) %>% 
                                 summarise(RAE_mean=mean(RAE_mean, na.rm=T), RAE_median=mean(RAE_median, na.rm=T)))
graph_1_simul



## Graph 2
# original data
graph_2 <- simul_graph2(data, HL_data, HL_list=c(0:100), simulation=F)
graph_2


# simulation data (n=300)
iter <- 300
pb <- txtProgressBar(max=iter, style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

graph_2_simul <- foreach(i=1:iter, .options.snow=opts, .combine='rbind', .packages=c('dplyr')) %dopar% {
  simul_graph2(data, HL_data, simulation=T)
}

graph_2_simul <- as.data.frame(graph_2_simul %>% group_by(time) %>% summarise(RAE_mean=mean(RAE_mean, na.rm=T),
                                                              RAE_median=mean(RAE_median, na.rm=T),
                                                              RAE_IQ_mean=mean(RAE_IQ_mean, na.rm=T),
                                                              RAE_nps_mean=mean(RAE_nps_mean, na.rm=T)))
graph_2_simul




## Graph 3
# original data
graph_3 <- simul_graph3(data, HL_data, HL_list=c(0,4,10,18,28,40,54,70,88,108,130,154), simulation=F)
graph_3



# simulation data (n=300)
iter <- 300
pb <- txtProgressBar(max=iter, style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

graph_3_simul <- foreach(i=1:iter, .options.snow=opts, .combine='rbind', .packages=c('dplyr')) %dopar% {
  simul_graph3(data, HL_data, HL_list=c(0,4,10,18,28,40,54,70,88,108,130,154), simulation=T)
}

graph_3_simul <- as.data.frame(graph_3_simul %>% group_by(Day, HL_group) %>% summarise(RAE=median(RAE, na.rm=T)))
graph_3_simul
