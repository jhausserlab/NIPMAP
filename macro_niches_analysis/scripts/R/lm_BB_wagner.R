rm(list=ls())
library(tidyverse)
library(ade4)
library(factoextra)
library(readxl)
library(reshape2)
library(plotly)
#library(igraph)



########******PRELIMINARY SETTINGS******########
#Set working directory
dirName <- dirname(rstudioapi::getSourceEditorContext()$path )
setwd(dirName)#setwd(dir="/srv/mfs/hausserlab/anissa.el/ImmuneStates/wagner_analysis")
MacroMicroData <- readRDS(file="./outputs/MacroMicroData.rds")

##########--------LINEAR REGRESSION--------##########
# VARIABLE TO EXPLAIN : cellular abundance of BC tumors (Wagner dataset)
# EXPLANATORY VARIABLES :
# BUILDING BLOCK 1 (TLS)
# BUILDING BLOCK 2 (INFLAMMATORY BLOCK)
# BUILDING BLOCK 3 (CANCER BLOCK)
# All variables are numerical continuous
# Building lbocks are vectors of size 10(number of cell types)
# The observed data (Y) is the dataset of wagner tumors of size
# 134 tumors x 10 cell types whose row sums are all equal to 100 (percentages)

Y <- MacroMicroData$MacroTumorsW
BB1 <- MacroMicroData$BBs[1,]
BB2 <- MacroMicroData$BBs[2,]
BB3 <- MacroMicroData$BBs[3,]

B <-t(MacroMicroData$BBs)

model1 <- lm(t(Y)~BB1 + BB2 + BB3 )
model2 <- lm(Y[1,]~B+0)
summary(model2)
## in lm(), we ocnsider that the response variable is as vector of observations
# Whereas here, we have just one observation of dimension 10 (we observe 10 
# cell ytpes proportiosn per tumor) 
# So the most ocnvenient way os to compute this as a matrix equation 

OmegaT <- solve(t(B) %*% B) %*% t(B) %*% t(Y)
# Round for easier viewing
Omega <- t(round(OmegaT, 2))
fig1 <-plot_ly(x=Omega[,1],
        y=Omega[,2],
        z=Omega[,3],
        type="scatter3d",
        mode="markers",
        name="Wagner BC tumors",
        showlegend=TRUE,
        marker=list(size=5, color ="red"))
bestcol <- rgb(0, 0, 255, max = 255, alpha = 125, names = "blue20")
fig1 <- fig1%>%add_mesh(x=c(1,0,0),
                        y=c(0,1,0),
                        z=c(0,0,1),
                        name = "BB plane",
                        facecolor= bestcol,
                        opacity=0.3,
                        inherit=FALSE,
                        showlegend=TRUE)

fig1 <- fig1 %>%layout(scene = list(xaxis = list(title = "BB1"), yaxis = list(title = "BB2"), zaxis = list(title = "BB3") ),
                       title = "BC tumors from Wagner dataset fitted in the Building Blocks space")
fig1


BBscores <- as_tibble(Omega,rownames=NA)%>%rownames_to_column(var="sample")%>%pivot_longer(c(BB1,BB2,BB3),names_to="building_block",values_to="proportion")
ggplot(data=BBscores,aes( x=sample,fill=building_block,y=proportion)) + 
  geom_bar(position="fill",stat="identity")+
  labs(title="Bar charts of inferred proportions of building blocks across samples of BC")+
  theme(axis.text.x=element_text(angle=45,hjust=0.8,vjust=0.2))+
  xlab("")
