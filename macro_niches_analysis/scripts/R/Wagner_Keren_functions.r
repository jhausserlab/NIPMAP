##### This script contains all the functions that process data from Wagner
## but also Keren  papers : they contain cell abundance of breast cancer samples
## Here, the aim is to bridge macro (Wagner) and micro-architecture (Keren) with
## the "building blocks theory"

## Prior to that, a reference space common to Wagner and Keren is created
## with the same cell types to allow comparison between them.

## Then, after PCA and archetypal analysis, we can compare the immune states 
## visible in the macro arhcitecture with the building blocks observed in local
## view
## The functions are mainly used in two files :
## scBC_analysis.Rmd
## lm_BB_wagner.Rmd

rm(list=ls())
library(tidyverse)
library(ade4)
library(factoextra)
#library(readxl)
library(reshape2)
library(plotly)
library(igraph)

#######------PROCESSING OF WAGNER DATASET------#######

set_new_cells_proportions <- function(CellsProp.raw){
  CellsProp.raw%>%filter(`LiveCells [%]`>50)
  #CellsPropBC <- left_join(CellsPropBC,MetaData,by="Patient")%>%filter(`Neoadjuvant Therapy Received` !="yes")
  CellsProp.raw <-  CellsProp.raw%>%
    mutate(LiveCells_tot = `EndothelialCells [%]`+`EpithelialCells [%]`+`Fibroblasts [%]`+`ImmuneCells [%]`+`Other [%]`)
  CellsProp.raw <-  CellsProp.raw%>%
    mutate(EndothelialCells = ((`EndothelialCells [%]`/LiveCells_tot))*100)%>%
    mutate(EpithelialCells = ((`EpithelialCells [%]`/LiveCells_tot))*100)%>%
    mutate(Fibroblasts = ((`Fibroblasts [%]`/LiveCells_tot))*100)%>%
    mutate(ImmuneCells = ((`ImmuneCells [%]`/LiveCells_tot))*100)%>%
    mutate(Other = ((`Other [%]`/LiveCells_tot))*100)
  
  #Correction to ImmuneCells 
  CellsProp.raw <-  CellsProp.raw%>%
    mutate(TCells = (`TCells [%]`/100)*ImmuneCells)%>%
    mutate(NKCells = (`NaturalKillerCells [%]`/100)*ImmuneCells)%>%
    mutate(Granulocytes = (`Granulocytes [%]`/100)*ImmuneCells)%>%
    mutate(BCells = (`BCells [%]`/100)*ImmuneCells)%>%
    mutate(PlasmaCells = (`PlasmaCells [%]`/100)*ImmuneCells)%>%
    mutate(DCs = (`plasmacytoidDendriticCells [%]`/100)*ImmuneCells)%>%
    mutate(MyeloidCells = (`MyeloidCells [%]`/100)*ImmuneCells)%>%
    mutate(Basophils = (`Basophils [%]`/100)*ImmuneCells)
  
  #%>%mutate(across(c(`TCells [%]`,`NaturalKillerCells [%]`,`Granulocytes [%]`,`BCells [%]`,`PlasmaCells [%]`,`plasmacytoidDendriticCells [%]`,`MyeloidCells [%]`,`Basophils [%]`),(.x/100)*ImmuneCells))
  
  ## Correction to T cells
  CellsProp.raw <- CellsProp.raw%>%
    mutate(sumTcells = `TCells CD4+ [%]`+`TCells CD8+ [%]`)%>%
    mutate(OtherTcells = ((100 - sumTcells)/100)*TCells)%>%
    mutate(CD4TCells = (`TCells CD4+ [%]`/100)*TCells)%>%
    mutate(CD8TCells = (`TCells CD8+ [%]`/100)*TCells)
  
  ## Correction to myeloid cells
  CellsProp.raw <- CellsProp.raw %>% mutate(Monocytes = (M06 + M15)*MyeloidCells)%>%
    mutate (Macrophages = (M01 + M02 + M03 + M04 + M08 + M09 + M11 + M13 + M14 + M16 + M17)*MyeloidCells)%>%
    mutate(OtherMyeloid = (M05 + M07 + M10 + M12 + M18 + M19)*MyeloidCells)%>%
    mutate(OtherImmune = OtherMyeloid + OtherTcells + Basophils + Granulocytes)
  
  return(CellsProp.raw)
}


#######--CREATION OF REFERENCE CELL SPACE ACROSS WAGNER & KEREN DATASETS--#######

## Creates incidence matrix from excel file
# @param: filename {str} name of the file containing the incidence matrix
# @return : GMatrix {matrix} incidence matrix
##
create_incidence_matrix <- function(filename)
{
  GMatrix<- readxl::read_excel(filename,col_names=TRUE)%>%column_to_rownames(var="...1")
  #Create incidence matrix
  GMatrix <- GMatrix%>%
    replace_na(replace=list("EndothelialCells"=0,"EpithelialCells"=0,
                                               "Fibroblasts"=0,"Other"=0,"NK"=0,
                                               "BCells"=0,"PlasmaCells"=0,"DC"=0,
                                               "CD4-T"=0,"CD8-T"=0,
                                               "Macrophages"=0,"OtherImmune"=0))
  return(GMatrix)
  
}


## Plots a graph from incidence matrix
# @param: IncidenceMatrix {matrix} nCellTypes(Keren)xnCellType(Wagner)
# @return: None
# Plots the bipartite graph of the Incidence Matrix that maps cell types from
# both datasets
##
plot_graph_incidence_matrix <- function(IncidenceMatrix){
  # FIrst we create from our incidence matrix 
  rownames(IncidenceMatrix)[4] <- "CD3+/CD4- T"
  rownames(IncidenceMatrix)[1] <- "CD8+ T"
  rownames(IncidenceMatrix)[9] <- "CD4+ T"
  colnames(IncidenceMatrix) <- c("EndothelialCells","EpithelialCells",
                                 "Fibroblasts","Other","NK",
                                 "BCells","PlasmaCells","DC",
                                 "CD4+T","CD8+T",
                                 "Macrophages","OtherImmune")
  graph <- graph_from_incidence_matrix(IncidenceMatrix, directed = FALSE,multiple=FALSE,
                                       mode = "in")
  from_incidence_matrix(graph)
  LO = layout_as_bipartite(graph) # Set up to bipartite graph (two groups : Keren and Wagner dataset)
  LO = LO[,c(2,1)]
  V(graph)$label.cex <- 0.82
  # Setting figure + other parameters
  pdf("./figs/bipartiteGraphCells.pdf",height=4, width=5)
  fig <- plot(graph,
       layout=LO,
       vertex.color = c("cornflowerblue","darkorange1")[V(graph)$type+1],
       vertex.shape = "circle",
       vertex.size = 8,
       label.font=2,
       vertex.label.dist=7,
       edge.width=2,
       asp = 0.95,
       margin=-0.05,
       vertex.label.degree= c(2*3.14,3.14)[V(graph)$type+1],
       main=" ")#"Graph (incidence matrix) of cell labels from Keren and Wagner data"
  fig
  dev.off()
}


## Computes the transformation matrix for Keren dataset
# @param: G {matrix}
# @return: Gw {matrix} transformation matrix for cellular % Wagner dataset
#
##
compute_Gw_matrix <- function(G){
  sum1 <- colSums(G)
  Gw<-G
  for (i in 1:length(sum1)){
    
    if(sum1[[i]] >1){
      ocurrences <-which(Gw[,i]==1)
      #print(ocurrences)
      #print("okk")
      firstIndex <-ocurrences[[1]] #match(1, GMatrix[,i])
      toDelete <- as.vector(which(Gw[,i]==1))[-1]
      #print(toDelete)
      rownames(Gw)[firstIndex]<-colnames(Gw)[i]
      Gw <- Gw[-toDelete,]
      #rm(ocurrences)
    }
  }
  return(Gw)
}

## Computes the transformation matrix for Wagner dataset
# @param: G {matrix}
# @return: Gk {matrix} transformation matrix for cell abundance dataset (Keren)
#
##
compute_Gk_matrix <- function(G){
  sum2 <- rowSums(G)
  Gk <-G
  for (i in 1:length(sum2)){
    
    if(sum2[[i]] >1){
      ocurrences <-which(Gk[i,]==1)
      #print(ocurrences)
      #print("okk")
      firstIndex <-ocurrences[[1]] #match(1, GMatrix[,i])
      toDelete <- as.vector(which(Gk[i,]==1))[-1]
      #print(toDelete)
      colnames(Gk)[firstIndex] <- rownames(Gk)[i]
      Gk <- Gk[,-toDelete]
    }
  }
  return(Gk)
}


####------DISSECTION OF HEALTHY TISSUE IN WAGNER BC SAMPLES-----####

## Computes Bmatrix which corresponds to the cell proportions of archetypes
# found in PCA from Wagner tumors 
# @param: ArchetypesCoordPCA matrix (nbPCs x nbArchetypes) coordinates of arch
#         in 3D PCA space (3 first Principal Components)
# @param: eigenVectPCA matrix (nbCellTypes x nbPCs) scores of each PC from Wagner data
# @param: meansPCA , vector of means for each variable used to do the PCA
# @return: B matrix (nbCellTYpes x nbArchetypes)
compute_cells_prop_archetypes <- function(ArchetypesCoordPCA,eigenVectPCA,meansPCA){
  B <- t(ArchetypesCoordPCA)%*% t(eigenVectPCA)   #t(scale(t(as.matrix(newPCAw$c1)),center=-(newPCAw$cent),scale=FALSE)) %*% archCoordw 
  ## Decentering the matrix to obtain percentages of cell types ==> the sum for each archetype must be =1
  B0 <- apply(t(B), 2, function(x) x+meansPCA)
  return(as.matrix(B0))
}



# Remove "healthy" archetype from tumor samples in silico
#@param: ArchToKeep {vector} vecotr of struings of archetype names to keep
#
# @return: cellAbundCorrected {matrix}
remove_healthy_archetype <- function(ArchCellsProp,ArchWeights,ArchToKeep, ArchToRemove){#,eigenVectPCA
  Bcorr0 <- ArchCellsProp[,ArchToKeep] #c("arch1", "arch2","arch4")
  ## Correction of the archetypes ==> remove the part of archetype 3 (made essentially of mesenchymal-like cells)
  archWPCAcorr <- matrix(nrow=0,ncol=0) #matrix(nrow=length(ArchToKeep),ncol=ncol(ArchScoresPCA))
  for (i in 1:length(ArchToKeep)) {
    archWPCAcorr <- rbind(archWPCAcorr,ArchWeights[ArchToKeep[i],]/(1-ArchWeights[ArchToRemove,]))
  }
  
  cellAbundCorrected <- t(as.matrix(Bcorr0) %*%as.matrix(archWPCAcorr))
  return(cellAbundCorrected)
}

remove_healthy_archetype2 <- function(eigenVectPCA,meansPCA,ArchetypesCoordPCA,ArchWeights,ArchToKeep, ArchToRemove){#,eigenVectPCA
  #Bcorr0 <- ArchCellsProp[,ArchToKeep] #c("arch1", "arch2","arch4")
  ## Correction of the archetypes ==> remove the part of archetype 3 (made essentially of adipocytes) 
  
  archWPCAcorr <- matrix(nrow=0,ncol=0) #matrix(nrow=length(ArchToKeep),ncol=ncol(ArchScoresPCA))
  for (i in 1:length(ArchToKeep)){
    archWPCAcorr <- rbind(archWPCAcorr,ArchWeights[ArchToKeep[i],]/(1-ArchWeights[ArchToRemove,]))
  }
  print(dim(archWPCAcorr))
  newCellProps <- scale(t(eigenVectPCA%*%(ArchetypesCoordPCA[,ArchToKeep]%*%t(archWPCAcorr))),center=-meansPCA,scale=FALSE)
  #cellAbundCorrected <- t(as.matrix(Bcorr0) %*%as.matrix(archWPCAcorr))
  return(newCellProps)#(cellAbundCorrected)
}


############-----FETCH/PROCESS DATA FROM KEREN ARTICLE-----############


###########-------LINEAR REGRESSION OF MACRO CELL ABUNDANCE ON BBS ------###########


# Returns Omega, the matrix of TMENs weights, solution of the equation:
# Y = Omega.B
# @param: Y {matrix} nSamples x nCellTYpes matrix of cellular abundance in BC
# @param: B {matrix} nCellTypes x nTMENs cellular abundance in TMENs
# @return: results {list} Omega {matrix} nSamples x nTMENs
#
#linear_regression <- function(Y, B){
  
  #OmegaT <- solve(t(B) %*% B) %*% t(B) %*% t(Y)
  # Round for easier viewing
  #Omega <- t(round(OmegaT, 2))
#  Omega <- Y %*% B %*% solve(t(B) %*% B)
#  residual <-(Y - Omega%*% t(B))# t(t(Y)- (B%*%OmegaT))
  #Ymean <-apply(Y, 1,function(x) mean(x))
#  R2 <- matrix(ncol = 1,nrow=nrow(Y)) # r square value for each sample
#  rownames(R2) <- rownames(Y)
#  colnames(R2) <- "R2"
#  R2 <- mapply(function(x, y){cor(x, y)}, as.data.frame(t(B%*%OmegaT)), as.data.frame(t(Y)))
  #apply(R2,1, function(x))
  #for (i in 1:nrow(Y)){
  #  R2[i] <- cor((B%*%OmegaT[,i]),t(Y)[,i],method = "pearson")^2#1- ((apply(Y-t(residual)))^2/(sum(Y-Ymean)^2))
  #}
  #rownames(R2) <- rownames(Y)
#  results = list("OmegaT" = OmegaT,
#                 "residuals" = residual,
#                 "R2_values"= R2)
#  return(results)
#}

# Performs linear regression on TMENs cell % of 
# macroscopic cellular abundance of tumors
# @param: Y{matrix} nbSamples x nbCellTypes, cell % across tumors
# @param: B{matrix} nbCellTypes x nbTMENs, cell % across TMENs
# @return: results{list}
#         OmegaT{matrix} nbTMENs x nbTMENs, proportions of TMENs across tumors
#         R2_values{list} R^2 values from linear regression across tumor samples
linear_regression <- function(Y, B){
  
  OmegaT <- solve(t(B) %*% B) %*% t(B) %*% t(Y)
  # Round for easier viewing
  Omega <- t(round(OmegaT, 3))
  #Omega <- round(Y %*% t(B) %*% solve(t(B) %*% B),3)
  residual <- (Y - Omega%*% t(B))# t(t(Y)- (B%*%OmegaT))
  #Ymean <-apply(Y, 1,function(x) mean(x))
  R2 <- matrix(ncol = 1,nrow=nrow(Y)) # r square value for each sample
  rownames(R2) <- rownames(Y)
  colnames(R2) <- "R2"
  R2 <- mapply(function(x, y){cor(x, y)}, as.data.frame(B%*%OmegaT), as.data.frame(t(Y)))
  #apply(R2,1, function(x))
  #for (i in 1:nrow(Y)){
  #  R2[i] <- cor((B%*%OmegaT[,i]),t(Y)[,i],method = "pearson")^2#1- ((apply(Y-t(residual)))^2/(sum(Y-Ymean)^2))
  #}
  #rownames(R2) <- rownames(Y)
  results = list("OmegaT" = OmegaT,
                 "residual" = residual,
                 "R2_values"= R2)
  return(results)
}

## Violin plot of R-square values from linear regression
# @param: R2{tibble} with first column called patients(str) and second:value
# @param: savePlot{lgl} (default=FALSE) save the file in pdf format in ./figs
# @param: nameToSave {str} (default=NULL)
# @return: vlnplt{ggplot object}
violinPlot_R_square <- function(R2,savePlot=FALSE,nameToSave=NULL){
  #label1 <- as.expression(bquote(italic(r)^2))
  vlnplt <-ggplot(data = R2,mapping = aes(x = patients,y = value,fill="#EBB064"))+ #77CA5D
    geom_violin()+
    theme(legend.position = "none")+
    ylim(0,1)+
    xlab("")+ylab(as.expression(bquote(R^2 ~ predicted ~ vs ~ observed)))
  
  vlnplt
  if(savePlot==TRUE & !is.null(nameToSave)){
    ggsave(paste0("./figs/",nameToSave,".pdf"), height = 3, width = 4) 
  }
  return(vlnplt)
}

# Returns scatter plot of obsrved versus predicted cell % of tumors 
# from lin regression for an average patient and eventually saves it in figs directory
# with following name: namePlot_measure.pdf
# @param: Ywagner{tibble} nbPatients x nbCellTypes
# @param: B{matrix} nbCellTypes x nbTMENs
# @param: R2{tibble} nbPatients x 1, col: value and rownames: patients names
#                   R^2 values from eq t(Y) = B.t(Omega)
#
# @param: error{dbl} (default: 0) error to add in case of log10 transform
# @param: log10Scale{lgl}(default: FALSE) plots in log10 scale if TRUE
# @param: measure{str} "mean" or "median" (default is mean) the patient whose
#                     r^2 value is either median or mean
# @param: savePlot{lgl} TRUE or FALSE (default: FALSE)
# @param: namePlot{str} name of the plot(_measure is added automatically) (default: NULL)
# @return: plot{ggplot object}
plot_obs_VS_pred_patient <- function(Y, B, R2, Omega,error = 0,log10Scale = FALSE,measure = "mean",savePlot = FALSE,namePlot = " "){
  Y <- Y + error
  B <- B + error
  R2 <- R2 %>%rownames_to_column(var = "names")

  colnames(Y) <- colnames(Y)%>%set_cells_names()#change cell labels for plot
  rownames(B) <- rownames(B)%>%set_cells_names()#change cell labels for plot
  plots <- list()
  #names(plots) <- rownames(Y)
  names_patients <- rownames(Y)
  #for (i in rownames(Y)){
  plots <- lapply(names_patients, function(x){
    avPatient <- R2[which(R2$names == x),]
    #print(i)
    namePlot=x
    #Create df of observed cell% (Y) and predicted cell % (B%*%t(Omega))
    ObsVSPred <- data.frame(Observed = t(Y)[,avPatient$names],
                             Predicted=(B%*%t(Omega)[,avPatient$names]))
    #ObsVSPred <- data.frame(Observed = t(Y),
    #                        Predicted=(B%*%t(Omega)))
    ObsVSPred<- as_tibble(ObsVSPred,rownames=NA)%>%
      rownames_to_column(var="cell_type")%>%
      mutate(cell_type = set_cells_names(cell_type))# set same cell labels

    #Scatterplot of cell % (Observed VS predicted)
    
    r2 <- round(as.numeric(avPatient$value),2)
    label1 <- as.expression(bquote(italic(r)^2== .(r2)))
    
    p <- ggplot(data = ObsVSPred,aes(x = Predicted,y = Observed))+
      geom_point()+
      ggrepel::geom_text_repel(aes(label = cell_type),max.overlaps=13)+
      geom_abline(slope = 1, intercept = 0) + 
      xlab("% predicted tumor cellular composition")+
      ylab("% observed tumor cellular composition ")+

    #stat_cor(method = "pearson", label.x = 3, label.y = 35)+
    annotate("text",x =5,y=40,label=avPatient$names)+
    annotate("text", x = 5, y = 33, label = label1)#expression(paste(italic(R)^2, " = ", format(r2, digits = 2)))
    #plot
    
    ## Saving fig in pdf file in figs directory
     if(savePlot==TRUE){
       figName  <- paste0("./figs/",namePlot,".pdf") #add the name of the measure (mean or median) in the name of the file
       ggsave(figName,p, height = 3, width = 4)
     }     
  })
 names(plots) <- names_patients
  return(plots)
}


###########------UTILS: Useful functions for figures------###########

# Returns modified names of cell types for graphs to make it more readable
# @param: x, a vector containing the names that need to be changed (Keren 
#         or Wagner dataset)
# @return: newNames, a vector, same length as the one entered as input
#          containing modified names
set_cells_names <- function(x){
  newNames <- str_replace_all(x,c(EndothelialCells = "endothelial",
                                  EpithelialCells = "cancer",
                                  OtherImmune = "other immune", 
                                  `CD4-T` = "CD4 T",`CD8-T` = "CD8 T",
                                  Macrophages = "macrophage",BCells = "B",
                                  PlasmaCells = "plasma B",
                                  Other = "healthy tissue",
                                  Fibroblasts = "fibroblast",
                                  DC = "DC", NK = "NK",
                                  `Other immune`="other immune",
                                  `DC / Mono`="DC / Mono",`CD3-T`="CD3 T",
                                  `Keratin-positive tumor` = "keratin+ cancer",
                                  Tumor = "cancer",Tregs = "Treg",
                                  `Mesenchymal-like` = "mesenchymal-like",
                                  Endothelial = "endothelial",
                                  Unidentified = "unidentified",
                                  `Mono / Neu` = "Mono / Neu", 
                                  Neutrophils = "neutrophil"))
  return(newNames)
}







