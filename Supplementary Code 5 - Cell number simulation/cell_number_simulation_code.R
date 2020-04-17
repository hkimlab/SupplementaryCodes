
# set working directory 
setwd("PathtoWorkDir") 
source("cell_number_simulation_function.R")

# Import package & function
library(doSNOW)
library(data.table)
library(ggplot2)
library(dplyr)  

## parallel computation setting
cl <- makeSOCKcluster(6) # input : number of cores
registerDoSNOW(cl)

## Import data
data <- read.csv("TestInput.csv")
test_data <- read.csv("TestInput.csv") 

HL_data <- read.csv("HL.csv")
HL_data <- na.omit(HL_data)
names(HL_data) <- c("ID","HL")

test_HL_data <- read.csv("HL.csv")
test_HL_data <- na.omit(test_HL_data)
names(test_HL_data) <- c("ID","HL")

## Simulation 

grid <- expand.grid(iter=1:100, HL_k=c(100,215,464,1000,2150,4640,10000,21500,46400,100000,215000,464000,1000000,2150000,4640000,10000000),Ft_k=c(100,215,464,1000,2150,4640,10000,21500,46400,100000,215000,464000,1000000,2150000,4640000,10000000))
# iter : number of simulation
# HL_k : total_k list for HL estimation
# Ft_k : total_k list for t_hat estimation

pb <- txtProgressBar(max=nrow(grid), style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

# 1
result_Ideal <- foreach(i=1:nrow(grid), .options.snow=opts, .combine='rbind', .multicombine=T) %dopar% {
  time_predict(iter=grid$iter[i], data=data, test_data=test_data, HL_data=HL_data, test_HL_data=test_HL_data,
               HL_k=grid$HL_k[i], Ft_k=grid$Ft_k[i], calc_type="median")
}

result_Ideal %>% group_by(HL_k, Ft_k) %>% summarise(mean=mean(RAE))
