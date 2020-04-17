# set working directory 
setwd("PathtoWorkDir") 
source("Cas9heterogeneity_simulation_function.R")

## Import package & function
library(doSNOW)
library(data.table)
library(ggplot2)

## parallel computation setting
cl <- makeSOCKcluster(6) # input : number of cores
registerDoSNOW(cl)

## Import data
data <- read.csv("TestInput.csv")
HL_data <- read.csv("HL.csv")
HL_data <- na.omit(HL_data)
names(HL_data) <- c("ID","HL")

## Simulation 
iter <- 300 # number of simulations
pb <- txtProgressBar(max=iter, style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

# heterogeneity 
# ExtremelyHigh
result_exthigh <- foreach(i=1:iter, .options.snow=opts, .combine='rbind', .multicombine=T) %dopar% {
  time_predict(data, HL_data, heterogeneity=100, k_scale=1, calc_type="mean")
}

# High
result_high <- foreach(i=1:iter, .options.snow=opts, .combine='rbind', .multicombine=T) %dopar% {
  time_predict(data, HL_data, heterogeneity=10, k_scale=1, calc_type="mean")
}

# Low-Medium
result_lowmedium <- foreach(i=1:iter, .options.snow=opts, .combine='rbind', .multicombine=T) %dopar% {
  time_predict(data, HL_data, heterogeneity=1, k_scale=1, calc_type="mean")
}

# Ideal simulation
result_Ideal <- foreach(i=1:iter, .options.snow=opts, .combine='rbind', .multicombine=T) %dopar% {
  time_predict(data, HL_data, Ideal=T, k_scale=1, calc_type="mean")
}


## Plotting
# ExtremelyHigh
ggplot(na.omit(result_exthigh), aes(x=HL, y=RAE)) + geom_point() + geom_smooth() +
  theme_bw() + xlab("Half Life") + scale_y_continuous(name="Relative abolute error") + 
  theme(axis.line = element_line(size=0.5, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.background = element_blank(),
        plot.title = element_text(size = 13, face = "bold", hjust = 0.5, vjust = 1),
        axis.title.x=element_text(colour="black", size = 15),
        axis.title.y=element_text(colour="black", size = 15),
        axis.text.x=element_text(colour="black", size = 13),
        axis.text.y=element_text(colour="black", size = 13),
        legend.position = "bottom",
        legend.text = element_text(colour="black", size = 15),
        legend.title = element_text(colour="black", size = 15))

# High
ggplot(na.omit(result_high), aes(x=HL, y=RAE)) + geom_point() + geom_smooth() +
  theme_bw() + xlab("Half Life") + scale_y_continuous(name="Relative abolute error") + 
  theme(axis.line = element_line(size=0.5, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.background = element_blank(),
        plot.title = element_text(size = 13, face = "bold", hjust = 0.5, vjust = 1),
        axis.title.x=element_text(colour="black", size = 15),
        axis.title.y=element_text(colour="black", size = 15),
        axis.text.x=element_text(colour="black", size = 13),
        axis.text.y=element_text(colour="black", size = 13),
        legend.position = "bottom",
        legend.text = element_text(colour="black", size = 15),
        legend.title = element_text(colour="black", size = 15))

# Low-Medium
ggplot(na.omit(result_lowmedium), aes(x=HL, y=RAE)) + geom_point() + geom_smooth() +
  theme_bw() + xlab("Half Life") + scale_y_continuous(name="Relative abolute error") + 
  theme(axis.line = element_line(size=0.5, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.background = element_blank(),
        plot.title = element_text(size = 13, face = "bold", hjust = 0.5, vjust = 1),
        axis.title.x=element_text(colour="black", size = 15),
        axis.title.y=element_text(colour="black", size = 15),
        axis.text.x=element_text(colour="black", size = 13),
        axis.text.y=element_text(colour="black", size = 13),
        legend.position = "bottom",
        legend.text = element_text(colour="black", size = 15),
        legend.title = element_text(colour="black", size = 15))

# No heterogeneity
ggplot(na.omit(result_Ideal), aes(x=HL, y=RAE)) + geom_point() + geom_smooth() +
  theme_bw() + xlab("Half Life") + scale_y_continuous(name="Relative abolute error") + 
  theme(axis.line = element_line(size=0.5, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.background = element_blank(),
        plot.title = element_text(size = 13, face = "bold", hjust = 0.5, vjust = 1),
        axis.title.x=element_text(colour="black", size = 15),
        axis.title.y=element_text(colour="black", size = 15),
        axis.text.x=element_text(colour="black", size = 13),
        axis.text.y=element_text(colour="black", size = 13),
        legend.position = "bottom",
        legend.text = element_text(colour="black", size = 15),
        legend.title = element_text(colour="black", size = 15))

