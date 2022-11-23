library(tidyverse)
library(fdrtool)
library(geometry)

set_tmens_names <- function(x){
  names1 <- str_replace_all(x,c(a1 = "TLS",
                                a2 = "inflammatory",
                                a3 = "cancer",
                                a4 = "fibrotic",
                                archetype1 = "TLS",
                                archetype2 = "inflamm",
                                archetype3 = "cancer",
                                archetype4 = "fibrotic",
                                a1a2 = fixed("TLSxinflamm"),
                                a1a3 = fixed("TLSxcancer"),
                                a1a4 = fixed("TLSxfibrotic"),
                                a2a3 = fixed("inflammxcancer"),
                                a2a4 = fixed("inflammxfibrotic"),
                                a3a4 = fixed("cancerxfibrotic")))
  return(names1)
}


## Plot ternary plot: 3d tetrahedron
# Each vertex represents TMENs coordinates (proportion)
# A cell is represented by a point with a neighborhood of TMENs
# Each point is colored by a marker
# Plot markers expressed in a specific cell type
plot_3D_ternary_cell_markers <- function (cellType, marker){
  MarkersCellsTMENS.ct <- MarkersCellsTMENS%>%
    filter(cell_type == cellType)%>%
    pivot_wider(id_cols = c("site_id","patient_id","arch1","arch2","arch3","arch4"),
                names_from = "marker",
                values_from = "value")
  X <- rbind(c(1, 0, 0),
             c(0, 1, 0),
             c(0, 0, 1),
             c(0, 0, 0))
  rownames(X) <- c("TLS", "inflammatory","cancer","fibrotic")
  
  Coords.3D <- bary2cart(X,as.matrix(MarkersCellsTMENS.ct [,c("arch1","arch2","arch3","arch4")]))
  markerCells <- MarkersCellsTMENS.ct[[marker]] #pull (MarkersCellsTMENS.ct,marker )
  distr <- quantile(markerCells,probs = c(0,0.1,0.5,0.8,0.99))
  if(distr[5] > min(markerCells)){
    markerCells[markerCells>distr[5]] <- distr[5]
  }
  print(distr)
  # simplex <- cbind(seq(0:1,0.001),seq(0:1,0.001),seq(0:1,0.001))
  fig1 <- plot_ly(x = Coords.3D[,1],
                  y = Coords.3D[,2],
                  z = Coords.3D[,3],
                  type="scatter3d",
                  color = ~markerCells,
                  alpha = 0.5,
                  mode="markers",
                  name= "cells",
                  marker = list(size= 3),
                  showlegend=TRUE)%>%
    colorbar(title = paste(marker, "in",cellType))%>%
    add_trace(x = X[,1],
              y = X[,2],
              z = X[,3],
              type = "scatter3d",mode = "markers+text",
              marker = list(size = 8),
              color = "red",
              text = ~rownames(X),#c("arch. 1", "arch. 2","arch. 3","arch. 4"),#c("arch. 2", "arch. 1","arch. 3","arch. 4"),
              textposition = c('top right','top right','bottom left','bottom left'),
              textfont = list(color = '#000000', size = 16),
              name = "archetypes",
              inherit=FALSE,
              showlegend = TRUE)%>%
    layout(scene = list(xaxis = list(title = "x"), 
                        yaxis = list(title = "y"), 
                        zaxis = list(title = "z") ))
  #fig1
  return(fig1)
}

## Get df with alpha values of TMENs and marker intensity value of cells of 
# a given cell type
get_tmens_coord_marker <- function(MarkersCellsTMENS,cellType, ma){
  if(ma =="H3K9ac/H3K27me3"){
    MarkersCellsTMENS.ct <- MarkersCellsTMENS%>%
      filter(cell_type == cellType & marker == "H3K9ac_H3K27me3")
  }else{
    MarkersCellsTMENS.ct <- MarkersCellsTMENS%>%
      filter(cell_type == cellType & marker == ma)
  }
  ###### TO CLEAN #######
  #%>%
  #pivot_wider(id_cols = c("site_id","patient_id","arch1","arch2","arch3","arch4"),
  #            names_from = "marker",
  #            values_from = "value")
  # machine precision is 2.220446e-16: least positive number x such as 1+ x!  = 1
  # 0 is < machine precision
  # each time we find a nb that is null, assign this value to avoid Nan, etc...
  # MarkersCellsTMENS.ct <- MarkersCellsTMENS.ct  %>%mutate(arch1=ifelse(arch1<.Machine$double.eps,eps,arch1))%>%
  #   mutate(arch2=ifelse(arch2 <.Machine$double.eps,eps,arch2))%>%
  #   mutate(arch3=ifelse(arch3 <.Machine$double.eps,eps,arch3))%>%
  #   mutate(arch4=ifelse(arch4 <.Machine$double.eps,eps,arch4))
  #return(MarkersCellsTMENS.ct%>%dplyr::select(-c(patient_id,site_id,cell_type_site,cell_type,marker)))#c(arch1,arch2,arch3,arch4,value)
  return(MarkersCellsTMENS.ct)
  
}

##### CORRECTION OF P VALUES ON SPEARMAN CORR TEST: FALSE DISCOVERY RATE
# Correleations between core/interface TMENs and marker expressed in cell type
correct_cor_fdr <- function(corMatrix,corCMtmens.pval,qThresh,corThresh){
  ### CORRECTIONS P VALUE FDR
  #dim(corCMtmens.pval)
  corCMtmens.pval.mat = corCMtmens.pval %>% column_to_rownames(var='names')
  fdrOut = 
    fdrtool(corCMtmens.pval.mat %>% as.matrix %>% as.numeric(), 
            statistic = "pvalue")
  corCMtmens.qval = 
    fdrOut$lfdr %>% 
    matrix(nrow=nrow(corCMtmens.pval.mat), ncol=ncol(corCMtmens.pval.mat),
           dimnames = list(rownames(corCMtmens.pval.mat), colnames(corCMtmens.pval.mat)))
  # dim(corCMtmens.qval)
  # dim(corCMtmens.pval.mat)
  # plot(corCMtmens.qval %>% as.matrix %>% as.numeric(), corCMtmens.pval.mat %>% as.matrix %>% as.numeric(), xlim=c(0,.1), ylim=c(0,.1))
  ## Qval and correlations thresholds
  # qThreshold <- 1/100
  # corThreshold <- .3
  # corThreshold2 <- .2
  filterMat <- corCMtmens.qval<qThresh & corMatrix[rownames(corCMtmens.qval),]>corThresh#corCMtmens.qval<qThresh & abs(corMatrix[rownames(corCMtmens.qval),])>corThresh
  return(filterMat)
}


##### COMPUTE SPEARMAN CORRELATIONS BETWEEN CELL PHENOTYPES AND TMENS 
# and filter out non-significant (FDR correction)/correlations under thresh
correlations_tmens_CM <- function(MarkersCellsTMENs,cellTypes, markers, qThresh=1/100, corThresh=0.3){
  corMatrix.pval <-
    matrix(nrow = length(cellTypes)*length(markers), ncol=4+6,
           dimnames = 
             list(c(),
                  c('a1','a2','a3',"a4","a1a2","a1a3","a1a4","a2a3","a2a4","a3a4")))#, ,"a1a2a3","a1a2a4","a1a3a4","a2a3a4","a1a2a3a4")))
  corMatrix <- corMatrix.pval
  cmNames <- c()
  id = 1
  ct = cellTypes[1]
  for(ct in cellTypes){
    #cat(sprintf("%s\n", ct))
    for(m in markers){
      data <- get_tmens_coord_marker(MarkersCellsTMENs,ct,m) 
      X <- data %>%dplyr::select(c(value))
      A1 <- data %>%dplyr::select(c(arch1,arch2,arch3,arch4))
      A2 <- data %>%dplyr::select(c(a1a2,a1a3,a1a4,a2a3,a2a4,a3a4))
      #A3 <- data %>%dplyr::select(c(a1a2a3,a1a2a4,a1a3a4,a2a3a4))
      #print(A1)
      cmNames[id] <- paste0(ct,";",m)
      corMatrix.pval[id,'a1'] <- cor.test(X%>%pull(value),A1%>%pull(arch1),method="spearman",exact=FALSE)$p.value
      corMatrix.pval[id,'a2'] <- cor.test(X%>%pull(value),A1%>%pull(arch2),method="spearman",exact=FALSE)$p.value
      corMatrix.pval[id,'a3'] <- cor.test(X%>%pull(value),A1%>%pull(arch3),method="spearman",exact=FALSE)$p.value
      corMatrix.pval[id,"a4"] <- cor.test(X%>%pull(value),A1%>%pull(arch4),method="spearman",exact=FALSE)$p.value
      
      corMatrix[id,'a1'] <- cor.test(X%>%pull(value),A1%>%pull(arch1),method="spearman",exact=FALSE)$estimate
      corMatrix[id,'a2'] <- cor.test(X%>%pull(value),A1%>%pull(arch2),method="spearman",exact=FALSE)$estimate
      corMatrix[id,'a3'] <- cor.test(X%>%pull(value),A1%>%pull(arch3),method="spearman",exact=FALSE)$estimate
      corMatrix[id,"a4"] <- cor.test(X%>%pull(value),A1%>%pull(arch4),method="spearman",exact=FALSE)$estimate
      #### A2 correlations ####
      corMatrix.pval[id,"a1a2"] <- cor.test(X%>%pull(value),A2%>%pull(a1a2),method="spearman",exact=FALSE)$p.value
      corMatrix.pval[id,"a1a3"] <- cor.test(X%>%pull(value),A2%>%pull(a1a3),method="spearman",exact=FALSE)$p.value
      corMatrix.pval[id,"a1a4"] <- cor.test(X%>%pull(value),A2%>%pull(a1a4),method="spearman",exact=FALSE)$p.value
      corMatrix.pval[id,"a2a3"] <- cor.test(X%>%pull(value),A2%>%pull(a2a3),method="spearman",exact=FALSE)$p.value
      corMatrix.pval[id,"a2a4"] <- cor.test(X%>%pull(value),A2%>%pull(a2a4),method="spearman",exact=FALSE)$p.value
      corMatrix.pval[id,"a3a4"] <- cor.test(X%>%pull(value),A2%>%pull(a3a4),method="spearman",exact=FALSE)$p.value
      
      corMatrix[id,"a1a2"] <- cor.test(X%>%pull(value),A2%>%pull(a1a2),method="spearman",exact=FALSE)$estimate
      corMatrix[id,"a1a3"] <- cor.test(X%>%pull(value),A2%>%pull(a1a3),method="spearman",exact=FALSE)$estimate
      corMatrix[id,"a1a4"] <- cor.test(X%>%pull(value),A2%>%pull(a1a4),method="spearman",exact=FALSE)$estimate
      corMatrix[id,"a2a3"] <- cor.test(X%>%pull(value),A2%>%pull(a2a3),method="spearman",exact=FALSE)$estimate
      corMatrix[id,"a2a4"] <- cor.test(X%>%pull(value),A2%>%pull(a2a4),method="spearman",exact=FALSE)$estimate
      corMatrix[id,"a3a4"] <- cor.test(X%>%pull(value),A2%>%pull(a3a4),method="spearman",exact=FALSE)$estimate
      #### A3 correlations ####
      # corMatrix.pval[id,"a1a2a3"] <- cor.test(X%>%pull(value),A3%>%pull(a1a2a3),method="spearman",exact=FALSE)$p.value
      # corMatrix.pval[id,"a1a2a4"] <- cor.test(X%>%pull(value),A3%>%pull(a1a2a4),method="spearman",exact=FALSE)$p.value
      # corMatrix.pval[id,"a1a3a4"] <- cor.test(X%>%pull(value),A3%>%pull(a1a3a4),method="spearman",exact=FALSE)$p.value
      # corMatrix.pval[id,"a2a3a4"] <- cor.test(X%>%pull(value),A3%>%pull(a2a3a4),method="spearman",exact=FALSE)$p.value
      # 
      # corMatrix[id,"a1a2a3"] <- cor.test(X%>%pull(value),A3%>%pull(a1a2a3),method="spearman",exact=FALSE)$estimate
      # corMatrix[id,"a1a2a4"] <- cor.test(X%>%pull(value),A3%>%pull(a1a2a4),method="spearman",exact=FALSE)$estimate
      # corMatrix[id,"a1a3a4"] <- cor.test(X%>%pull(value),A3%>%pull(a1a3a4),method="spearman",exact=FALSE)$estimate
      # corMatrix[id,"a2a3a4"] <- cor.test(X%>%pull(value),A3%>%pull(a2a3a4),method="spearman",exact=FALSE)$estimate
      # 
      # corMatrix[id,"a1a2a3a4"] <-
      #   cor.test(X%>%pull(value),data%>%pull(a1a2a3a4),method="spearman",exact=FALSE)$estimate
      # corMatrix.pval[id,"a1a2a3a4"] <-
      #   cor.test(X%>%pull(value),data%>%pull(a1a2a3a4),method="spearman",exact=FALSE)$p.value
      
      id = id+1
    }
  }
  rownames(corMatrix) <- cmNames
  rownames(corMatrix.pval) <- cmNames
  
  corCMtmens.pval <- corMatrix.pval%>%as_tibble(rownames=NA)%>%
    rownames_to_column(var="names")%>%
    #separate(names,into=c("cell_type","marker"),sep=";")%>%
    drop_na()
  
  corCMtmens <- corMatrix%>%as_tibble(rownames=NA)%>%
    rownames_to_column(var="names")%>%
    #separate(names,into=c("cell_type","marker"),sep=";")%>%
    drop_na()
  write_csv(corCMtmens,"./outputs/rawCorMatrix.csv")
  ### CORRECTIONS P VALUE FDR
  # library(fdrtool)
  # dim(corCMtmens.pval)
  # corCMtmens.pval.mat = corCMtmens.pval %>% column_to_rownames(var='names')
  # fdrOut = 
  #   fdrtool(corCMtmens.pval.mat %>% as.matrix %>% as.numeric(), 
  #           statistic = "pvalue")
  # corCMtmens.qval = 
  #   fdrOut$lfdr %>% 
  #   matrix(nrow=nrow(corCMtmens.pval.mat), ncol=ncol(corCMtmens.pval.mat),
  #          dimnames = list(rownames(corCMtmens.pval.mat), colnames(corCMtmens.pval.mat)))
  # # dim(corCMtmens.qval)
  # # dim(corCMtmens.pval.mat)
  # # plot(corCMtmens.qval %>% as.matrix %>% as.numeric(), corCMtmens.pval.mat %>% as.matrix %>% as.numeric(), xlim=c(0,.1), ylim=c(0,.1))
  # ## Qval and correlations thresholds
  # qThreshold <- 1/100
  # corThreshold <- .3
  # corThreshold2 <- .2
  #filterMat <- corCMtmens.qval<qThreshold & abs(corMatrix[rownames(corCMtmens.qval),])>corThreshold
  #saveRDS(corMatrix,"./corMatrix_raw.rds")
  filterMat <- correct_cor_fdr(corMatrix,corCMtmens.pval,qThresh,corThresh)
  mean(apply(filterMat, 1, sum) > 0)
  ## 90% of the marker are associated to at least one TMENs or interface!
  
  cMf <- corMatrix[names(which(apply(filterMat, 1, sum) > 0)),]
  #write_csv(cMf%>%as_tibble(rownames=NA)%>%rownames_to_column(var="cell_marker"), "./CM_corr_tmens_thresh0.3.csv")
  return(cMf)
}


#### Enriched cell types in each site 
# choose the 5% sites closest to the archetype
# compute the SD for each cell type
# compare the average density  d1 of a cell type from the top 5% closest sites with the av denisty of cell type from the rest of the sites d2
# if d1 > d2+ 1* SD ==> the cell type is representative  of the archetype
# do the same for the other archetypes

get_cell_type_enriched_arch <- function(cellAb.df=archsSitesCellAb,archetype="",thresh=0.99){
  sites_a <- cellAb.df%>%filter(.data[[archetype]] >quantile(sort(pull(cellAb.df,.data[[archetype]]),decreasing=TRUE),thresh))
  #print(sites_a)
  summaryA <- sites_a %>%group_by(cell_type)%>%summarise(meanCT=mean(cell_density),sdCT = sd(cell_density))%>%column_to_rownames(var="cell_type")
  #print(summaryA)
  sites_notA <- anti_join(cellAb.df ,sites_a,by = c("site_id","patient_id"))
  #print(sites_notA)
  summaryNA <- sites_notA%>%group_by(cell_type)%>%summarise(meanCT=mean(cell_density),sdCT = sd(cell_density))%>%column_to_rownames(var="cell_type")
  #print(summaryNA)
  #Arch.CT <-  data.frame(region = c(), cell_type)#c()
  cts <- c()
  for(c in CELLTYPES){
    if(as.numeric(summaryA[c,"meanCT"]) > (as.numeric(summaryNA[c,"meanCT"]) + as.numeric(summaryNA[c,"sdCT"]))){
      #print(c)
      #Arch.CT <-append(Arch.CT,c)
      cts <- append(cts,c)
    }
  }
  Archs.CT <- data.frame(niche = rep(archetype,length(cts)),cell_type=cts)
  return(Archs.CT)
}


### GET CELL TYPES FROM ALL ARCHEYTPES
# get cell types enriched in each archetype/TMEN
get_CT_enriched_all_archs<- function(archsSites =archetypesSitesKeren, cellAbundSites = cellAbSites,thresh=0.99){
  archs <- colnames(archetypesSitesKeren)[grepl("arch", colnames(archetypesSitesKeren))]# vector of archetypes
  # list of cell types enriched in each archetype
  #archetypes.ct <- list() 
  #names(archetypes.ct) <- archs
  archetypes_cts.all <- data.frame(niche=NULL, cell_type=NULL)
  archsSitesCellAb <- left_join(archetypesSitesKeren, cellAbSites,by=c("site_id","patient_id"))%>%pivot_longer(cols=append(CELLTYPES,"Unidentified"),names_to="cell_type",values_to= "cell_density")%>%ungroup()
  for (a in archs){
    res <- get_cell_type_enriched_arch(cellAb.df=archsSitesCellAb,archetype=a)
    archetypes_cts.all<- rbind(archetypes_cts.all,res)#archetypes.ct[[a]] <-res
  }
  archetypes_cts.all$niche <-  set_tmens_names(archetypes_cts.all$niche)
  return(archetypes_cts.all)#(archetypes.ct)
}



# correlations_tmens_CM <- function(MarkersCellsTMENs,cellTypes=CELLTYPES, markers=FuncMarkers, CMtmensPos, qThresh=1/100, corThresh=0.3){
#   corMatrix.pval <-
#     matrix(nrow = length(CELLTYPES)*length(FuncMarkers), ncol=4+6+4+1,
#            dimnames = 
#              list(c(),
#                   c('a1','a2','a3',"a4","a1a2","a1a3","a1a4","a2a3","a2a4", 
#                     "a3a4","a1a2a3","a1a2a4","a1a3a4","a2a3a4","a1a2a3a4"))
#     )
#   corMatrix <- corMatrix.pval
#   cmNames <- c()
#   id = 1
#   ct = cellTypes[1]
#   for(ct in cellTypes){
#     #cat(sprintf("%s\n", ct))
#     for(m in markers){
#       data <- get_tmens_coord_marker(MarkersCellsTMENs,ct,m) 
#       X <- data %>%dplyr::select(c(value))
#       A1 <- data %>%dplyr::select(c(arch1,arch2,arch3,arch4))
#       A2 <- data %>%dplyr::select(c(a1a2,a1a3,a1a4,a2a3,a2a4,a3a4))
#       A3 <- data %>%dplyr::select(c(a1a2a3,a1a2a4,a1a3a4,a2a3a4))
#       #print(A1)
#       cmNames[id] <- paste0(ct,";",m)
#       corMatrix.pval[id,'a1'] <- cor.test(X%>%pull(value),A1%>%pull(arch1),method="spearman",exact=FALSE)$p.value
#       corMatrix.pval[id,'a2'] <- cor.test(X%>%pull(value),A1%>%pull(arch2),method="spearman",exact=FALSE)$p.value
#       corMatrix.pval[id,'a3'] <- cor.test(X%>%pull(value),A1%>%pull(arch3),method="spearman",exact=FALSE)$p.value
#       corMatrix.pval[id,"a4"] <- cor.test(X%>%pull(value),A1%>%pull(arch4),method="spearman",exact=FALSE)$p.value
#       
#       corMatrix[id,'a1'] <- cor.test(X%>%pull(value),A1%>%pull(arch1),method="spearman",exact=FALSE)$estimate
#       corMatrix[id,'a2'] <- cor.test(X%>%pull(value),A1%>%pull(arch2),method="spearman",exact=FALSE)$estimate
#       corMatrix[id,'a3'] <- cor.test(X%>%pull(value),A1%>%pull(arch3),method="spearman",exact=FALSE)$estimate
#       corMatrix[id,"a4"] <- cor.test(X%>%pull(value),A1%>%pull(arch4),method="spearman",exact=FALSE)$estimate
#       #### A2 correlations ####
#       corMatrix.pval[id,"a1a2"] <- cor.test(X%>%pull(value),A2%>%pull(a1a2),method="spearman",exact=FALSE)$p.value
#       corMatrix.pval[id,"a1a3"] <- cor.test(X%>%pull(value),A2%>%pull(a1a3),method="spearman",exact=FALSE)$p.value
#       corMatrix.pval[id,"a1a4"] <- cor.test(X%>%pull(value),A2%>%pull(a1a4),method="spearman",exact=FALSE)$p.value
#       corMatrix.pval[id,"a2a3"] <- cor.test(X%>%pull(value),A2%>%pull(a2a3),method="spearman",exact=FALSE)$p.value
#       corMatrix.pval[id,"a2a4"] <- cor.test(X%>%pull(value),A2%>%pull(a2a4),method="spearman",exact=FALSE)$p.value
#       corMatrix.pval[id,"a3a4"] <- cor.test(X%>%pull(value),A2%>%pull(a3a4),method="spearman",exact=FALSE)$p.value
#       
#       corMatrix[id,"a1a2"] <- cor.test(X%>%pull(value),A2%>%pull(a1a2),method="spearman",exact=FALSE)$estimate
#       corMatrix[id,"a1a3"] <- cor.test(X%>%pull(value),A2%>%pull(a1a3),method="spearman",exact=FALSE)$estimate
#       corMatrix[id,"a1a4"] <- cor.test(X%>%pull(value),A2%>%pull(a1a4),method="spearman",exact=FALSE)$estimate
#       corMatrix[id,"a2a3"] <- cor.test(X%>%pull(value),A2%>%pull(a2a3),method="spearman",exact=FALSE)$estimate
#       corMatrix[id,"a2a4"] <- cor.test(X%>%pull(value),A2%>%pull(a2a4),method="spearman",exact=FALSE)$estimate
#       corMatrix[id,"a3a4"] <- cor.test(X%>%pull(value),A2%>%pull(a3a4),method="spearman",exact=FALSE)$estimate
#       #### A3 correlations ####
#       corMatrix.pval[id,"a1a2a3"] <- cor.test(X%>%pull(value),A3%>%pull(a1a2a3),method="spearman",exact=FALSE)$p.value
#       corMatrix.pval[id,"a1a2a4"] <- cor.test(X%>%pull(value),A3%>%pull(a1a2a4),method="spearman",exact=FALSE)$p.value
#       corMatrix.pval[id,"a1a3a4"] <- cor.test(X%>%pull(value),A3%>%pull(a1a3a4),method="spearman",exact=FALSE)$p.value
#       corMatrix.pval[id,"a2a3a4"] <- cor.test(X%>%pull(value),A3%>%pull(a2a3a4),method="spearman",exact=FALSE)$p.value
#       
#       corMatrix[id,"a1a2a3"] <- cor.test(X%>%pull(value),A3%>%pull(a1a2a3),method="spearman",exact=FALSE)$estimate
#       corMatrix[id,"a1a2a4"] <- cor.test(X%>%pull(value),A3%>%pull(a1a2a4),method="spearman",exact=FALSE)$estimate
#       corMatrix[id,"a1a3a4"] <- cor.test(X%>%pull(value),A3%>%pull(a1a3a4),method="spearman",exact=FALSE)$estimate
#       corMatrix[id,"a2a3a4"] <- cor.test(X%>%pull(value),A3%>%pull(a2a3a4),method="spearman",exact=FALSE)$estimate
#       
#       corMatrix[id,"a1a2a3a4"] <-
#         cor.test(X%>%pull(value),data%>%pull(a1a2a3a4),method="spearman",exact=FALSE)$estimate
#       corMatrix.pval[id,"a1a2a3a4"] <-
#         cor.test(X%>%pull(value),data%>%pull(a1a2a3a4),method="spearman",exact=FALSE)$p.value
#       
#       id = id+1
#     }
#   }
#   rownames(corMatrix) <- cmNames
#   rownames(corMatrix.pval) <- cmNames
#   
#   corCMtmens.pval <-corMatrix.pval%>%as_tibble(rownames=NA)%>%
#     rownames_to_column(var="names")%>%
#     #separate(names,into=c("cell_type","marker"),sep=";")%>%
#     drop_na()
#   
#   corCMtmens <-corMatrix%>%as_tibble(rownames=NA)%>%
#     rownames_to_column(var="names")%>%
#     #separate(names,into=c("cell_type","marker"),sep=";")%>%
#     drop_na()
#   
#   ### CORRECTIONS P VALUE FDR
#   library(fdrtool)
#   dim(corCMtmens.pval)
#   corCMtmens.pval.mat = corCMtmens.pval %>% column_to_rownames(var='names')
#   fdrOut = 
#     fdrtool(corCMtmens.pval.mat %>% as.matrix %>% as.numeric(), 
#             statistic = "pvalue")
#   corCMtmens.qval = 
#     fdrOut$lfdr %>% 
#     matrix(nrow=nrow(corCMtmens.pval.mat), ncol=ncol(corCMtmens.pval.mat),
#            dimnames = list(rownames(corCMtmens.pval.mat), colnames(corCMtmens.pval.mat)))
#   # dim(corCMtmens.qval)
#   # dim(corCMtmens.pval.mat)
#   # plot(corCMtmens.qval %>% as.matrix %>% as.numeric(), corCMtmens.pval.mat %>% as.matrix %>% as.numeric(), xlim=c(0,.1), ylim=c(0,.1))
#   
#   
#   ## Qval and correlations thresholds
#   qThreshold <- 1/100
#   corThreshold <- .3
#   corThreshold2 <- .2
#   filterMat <- corCMtmens.qval<qThresh & abs(corMatrix[rownames(corCMtmens.qval),])>corThresh
#   
#   mean(apply(filterMat, 1, sum) > 0)
#   ## 90% of the marker are associated to at least one TMENs or interface!
#   
#   cMf <- corMatrix[names(which(apply(filterMat, 1, sum) > 0)),]
#   #write_csv(cMf%>%as_tibble(rownames=NA)%>%rownames_to_column(var="cell_marker"), "./CM_corr_tmens_thresh0.3.csv")
#   
#   return(cMf)
# }


