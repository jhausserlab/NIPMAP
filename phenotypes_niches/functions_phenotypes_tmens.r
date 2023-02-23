library(tidyverse)
library(fdrtool)
#library(geometry)
library(cluster)
library(broom)
library(pROC)
library(ggpubr)
library(pheatmap)


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

## Get df with niches weights and marker intensity value of cells of 
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

## Scale the values of intensity markers for each cluster
ZscoreScale <- function(M, cutOff=3) {
  M0 = apply(M, 2, function(x) {(x-mean(x)) / sd(x)}) %>% t
  M0[M0>cutOff] = cutOff
  M0[M0<(-cutOff)] = -cutOff
  return(M0)
}

## Hierarchical clustering on cells from the same type
# ward's method (variance minimization, scaled with cut-off value of markers intensity =3)
# @param dfCellsMarkers: df with cell_type, cell_id columns and markers columns and cells as rows
# @param markers: string vector of functional markers used for clustering
# @param cellTyoes: string vector of cell types to be clustered 
# @param nClusts: number of clusters for hierarchical clustering(same as number of core/interfaces used in NIPMAP for comparing both methods)
hclust_cellTypes <- function(dfCellsMarkers = cell_markers,cellTypes,markers=FuncMarkers,nClusts=10){
  # select cell types from cellTypes variable in df 
  CellsFunctMarkers.df <- dfCellsMarkers[which(dfCellsMarkers$cell_type %in% cellTypes),append(append(append(all_of(markers),"cell_type"),"cell_id"),"patient_id")]# "Keratin-positive tumor"
  
  rownames(CellsFunctMarkers.df) <- make.names(paste0(CellsFunctMarkers.df$cell_type,"_",CellsFunctMarkers.df$cell_id,"_",CellsFunctMarkers.df$patient_id),unique=TRUE)
  # create emarker expression df for each cell type
  cells.markers <- group_split(CellsFunctMarkers.df%>%rownames_to_column(var="id")%>%group_by(cell_type))
  #print(cells.markers)
  list_heatmaps <- lapply(cells.markers, FUN = function(x) {
    ct <- unique(x$cell_type)
    x <- x%>%column_to_rownames(var="id")
    
    #Z-score normalization
    cm.mat <- x%>%dplyr::select(-c(cell_type,cell_id,patient_id))%>%ZscoreScale()%>%t()
    # hierarchical clustering
    CM.dist <- dist(cm.mat,method="euclidian")# computing euclidean distances in markers space 
    Cells.hclust <- hclust(CM.dist,method = "ward.D")
    # Cluster ids assigned to each cell: using the same nb of clusts 
    # as nb of niches/interfaces in NIPMAP
    sub_grp <- data.frame(cell_id=rownames(x),cluster=cutree(Cells.hclust, k =nClusts),row.names=NULL)#(cell_id=seq(1,nrow(x),1),cluster=cutree(Cells.hclust, k = 10),row.names=NULL)
    #write clusters ids in csv file for each cell
    write_csv(sub_grp%>%separate(cell_id,into=c("cell_type","label","SampleID"),sep="_"),paste0("./outputs/",ct,"_hclustIDs.csv"))
    #computing mean value of Z-scores for each cluster for heatmap
    clustMarkers.zscore <- cm.mat %>%as_tibble()%>%mutate(cell_id=rownames(x))%>% #as_tibble(rownames=NA)%>%rownames_to_column(var="cell_id")%>%
      left_join(.,sub_grp, by="cell_id")%>%
      pivot_longer(cols = FuncMarkers, names_to="marker",values_to = "value")%>%
      group_by(cluster,marker)%>%
      mutate(meanMClust = mean(value))%>%
      distinct(cluster,marker,meanMClust)%>%
      arrange(marker)%>%
      pivot_wider( names_from='marker',values_from = "meanMClust")%>%column_to_rownames(var="cluster")
    pdf(paste0("./figs/",ct,"cells_clusters.pdf"),width=6.5,height=5)
    pheatmap(as.matrix(clustMarkers.zscore%>%t()),cluster_cols = FALSE,cluster_rows = FALSE,main=ct,cell_width=1.5,cell_height = 2)
    dev.off()
  })
  return(list_heatmaps)
}

## Preliminary analysis of local neighborhood and cells phenotypes
# Scatter plot of expression of marker in a cell type across niche weight 
# of cells (including also interfaces with niches)
# including correlation stats
# CellsNiches.df: df of expression of markers of all cells & their niches weights
# cellType: str cell type
# Marker: str marker expression in cellType
# niche: (if str) name of niche or (if vector of str) interface of niches
plot_cellMarker_niche <- function(CellsNiches.df=MarkersCellsTMENS2,cellType, Marker,niche){
  CellsNiches.df1 <- CellsNiches.df%>%dplyr::rename(TLS=arch1)%>%
    dplyr::rename(inflammatory = arch2)%>%dplyr::rename(fibrotic=arch4)%>%
    dplyr::rename(cancer=arch3)
  if (length(niche)==1){
    plot1 <- ggplot(data=CellsNiches.df1%>%filter(cell_type==cellType & marker==Marker),aes_string(x=niche,y="value"))+
      geom_point(alpha=.1)+
      stat_cor(method="spearman",cor.coef.name="rho",label.x=0.1)+
      stat_cor(method="pearson",cor.coef.name="R",label.y=4.0,label.x=0.1)+
      xlab(paste0(niche,collapse=" x "))+ ylab(paste0(Marker," expression"))+
      ggtitle(paste0(Marker,"⁺ ",cellType," VS ",niche))
  }
  else{
    #print(head(CellsNiches.df1%>%filter(cell_type==cellType & marker==Marker)%>%mutate(product=eval(rlang::parse_expr((paste0(paste0(niche,collapse="*"))))))))
    plot1 <- ggplot(data=CellsNiches.df1%>%filter(cell_type==cellType & marker==Marker)%>%mutate(product=eval(rlang::parse_expr((paste0(paste0(niche,collapse="*")))))),aes(x=product,y=value))+
      geom_point(alpha=.1)+
      stat_cor(method="spearman",cor.coef.name="rho",label.x=0.1)+
      stat_cor(method="pearson",cor.coef.name="R",label.y=4.0,label.x=0.1)+
      xlab(paste0(niche,collapse=" x "))+ ylab(paste0(Marker," expression"))+
      ggtitle(paste0(Marker,"⁺ ",cellType," VS ",paste0(niche,collapse=" x ")))
    
  }
  return(plot1)
}


## Linear regression of niches weight over marker expression of cell type
# markersCells.niches: 
# markers: 
# cell_types: 
# coreIntfNiches: string vector of archetypes names 
lm_niches_phen <- function(markersCells.niches=MarkersCellsTMENS2,markers=FuncMarkers, cell_types=CellT,coreIntfNiches=coreIntf2){
  MarkersCellsTMENs.byNiche <- markersCells.niches%>%
    filter(marker %in% markers)%>%
    pivot_wider(names_from="marker",values_from="value")%>%
    pivot_longer(cols = coreIntfNiches,names_to="niches",values_to="alpha")%>%
    filter(cell_type %in% cell_types )%>%
    filter(niches %in%  coreIntfNiches)%>%
    mutate(id = paste(patient_id,site_id,cell_type,niches,sep="_"))%>%
    column_to_rownames(var="id")
  
  ### Split table of cells niches into sub df by cell type and niche(core + interfaces)
  NichesCellPhen <- group_split(MarkersCellsTMENs.byNiche%>%group_by(niches,cell_type))
  
  ## Fit linear regression
  lmform <- paste0("`", markers, "`", collapse = " + ")# make linear regression formula
  formula1 <- as.formula(paste0("alpha~",lmform))
  lmfit.stats <- lapply(NichesCellPhen ,function(x){
    #nichesNames<-append(nichesNames,unique(pull(x,niches)))
    #cellT<-append(cellT,unique(pull(x,cell_type)))
    x <- x%>%mutate(id = paste(patient_id,site_id,cell_type,niches,sep="_"))%>%
      column_to_rownames(var="id")
    #print(x$count)
    #print(unique(pull(x,niches)))
    LinReg <- lm(formula1,data=x)
    tidy(LinReg)%>%mutate(term=paste(unique(pull(x,cell_type)),term,sep=";"))%>%mutate(niche=rep(unique(x$niches)))
  })
  ## Remove intercept estimate from linear regression
  dfStats <- lmfit.stats%>%reduce(full_join, by = c("term","niche","estimate","std.error","statistic","p.value"))
  dfStats2 <- dfStats[!grepl("(Intercept)",dfStats$term),]
  #head(dfStats2)
  #dfStats%>%pull(term)%>%unique
  return(dfStats2)
}

## FDR correction of estimates p values from linear regression
fdr_lm_niches <- function(lmStats,qvalTresh=0.01){
  mat.tstat <- as.matrix(lmStats%>%dplyr::select(-c(estimate,std.error,p.value))%>%pivot_wider(names_from = "niche",values_from="statistic")%>%column_to_rownames(var="term"))
  
  mat.pval <- as.matrix(lmStats%>%dplyr::select(-c(estimate,std.error,statistic))%>%pivot_wider(names_from = "niche",values_from="p.value")%>%column_to_rownames(var="term"))
  
  ## Check NA values are from the same rows from both p value and t stat matrices
  rownames(mat.pval[rowSums(is.na(mat.pval)) > 0,]) ==rownames(mat.tstat[rowSums(is.na(mat.tstat))> 0,])
  
  mat.tstat2 <- mat.tstat[rowSums(is.na(mat.tstat))==0,]#as_tibble(mat.tstat,rownames=NA)%>%drop_na()%>%as.matrix(rownames=TRUE)
  
  mat.pval2 <-  mat.pval[rowSums(is.na(mat.pval))==0,]#as_tibble(mat.pval,rownames=NA)%>%drop_na()%>%as.matrix(rownames=TRUE)
  ## FDR correction for p-value
  
  fdrCorr <- fdrtool(mat.pval2%>% as.numeric(),
                     statistic = "pvalue")
  CMtmens.qval<- 
    fdrCorr$qval %>% 
    matrix(nrow=nrow(mat.pval2), ncol=ncol(mat.pval2),
           dimnames = list(rownames(mat.pval2), colnames(mat.pval2)))
  ## Set to 0 non-significant estimates
  mat.tstat2[CMtmens.qval>qvalTresh]=0
  ## Check estimates that are non-significant for none of the niches
  mat.tstat2[rowSums(mat.tstat2==0)==ncol(mat.tstat2),]
  ## Filter matrix of t stats ==> remove rows that are all non-signif, keep with at least 1 non signif association with a tmen
  filtered.tstat <- mat.tstat2[rowSums(mat.tstat2==0)<ncol(mat.tstat2),]#CMtmens.qval<0.01 & 
  return(filtered.tstat)
}


##Plot ROC curves of niche classification for cell types
#clustFiles: name of .csv hclusts files
# CT: cell types to plot
# markers: functional markers= FuncMarkers
# clustIDs: list of number of clusters(obtained from hclust) for each cell type(names of elements: cell type)
# MarkersCells.niches: markers epxression of all cells + their niches weights
# nicheNames: str vector of archetypes names
# cellsNiches: cells niches weights
plot_roc_niches_cells <- function(CT,clustFiles,clustIDs, nichesNames,markers,cellsNiches,MarkersCells.niches){
  # lmform<- paste0("`", markers, "`", collapse = " + ")# linear regression of niche weight over functional markers expressions
  # formula1 <- as.formula(paste0("alpha~",lmform))
  lmform <- paste0("`", markers, "`", collapse = " + ")# make linear regression formula
  formula1 <- as.formula(paste0("alpha~",lmform))
  #print(lmform)
  #print(formula1)
  MarkersCellsTMENs.byNiche <- MarkersCells.niches%>%
    filter(marker %in% markers)%>%
    pivot_wider(names_from="marker",values_from="value")%>%
    pivot_longer(cols = coreIntfNiches,names_to="niches",values_to="alpha")%>%
    filter(cell_type %in% CT)%>%
    filter(niches %in%  coreIntfNiches)%>%
    mutate(id = paste(patient_id,site_id,cell_type,niches,sep="_"))%>%
    column_to_rownames(var="id")
  
  ROCplots <- lapply(CT, function(x){
    ct <- x
    for (n in nichesNames){
      ar <- n
      print(ar)
      print(ct)
      #print(ar)
      S1 <- c()
      S2 <- c()
      ## GET hclusters of cell type
      hclusts.cells <- read.csv(paste0("./outputs/",ct,fileExt))%>%mutate(value=1)%>%pivot_wider(id_cols=c("cell_type","label","SampleID"),names_from="cluster",values_from="value",values_fill=0)
      
      ### GET observed nniche weight and set thresh to niche positivity to .5
      hclustsNiches.cells <- left_join(hclusts.cells%>%dplyr::rename(patient_id = SampleID)%>%dplyr::rename(site_id = label),cellsNiches%>%dplyr::rename(cell_type = cell_type_site),by = c("patient_id","cell_type","site_id"))%>%
        mutate(AOI = get(ar))%>%
        drop_na(AOI)%>%
        mutate(archPos = ifelse(AOI>=.5,1,0))## We define cells that belong or not to archetype if their niche wieght is >.5
      ### Get predicted proportions of inflammatory niche of DCs 
      
      Predlm <- predict.lm(lm(formula1, data =MarkersCellsTMENs.byNiche%>%
                                filter(cell_type==ct & niches==ar)))#mutate(patient_id=str_split(id,"_")[[1]])
      Pred.df <- data.frame(id=names(Predlm), value = Predlm,row.names = NULL)%>%mutate(niche_pred = ifelse(value>=.5,1,0))%>%
        separate(col=id,into= c("SampleID","label","cell_type","arch"),sep="_",convert=TRUE)
      
      NichesPredClust <- left_join(hclustsNiches.cells, Pred.df%>%dplyr::rename(patient_id=SampleID)%>%dplyr::rename(site_id=label),by=c("patient_id","cell_type","site_id"))
      ## Compute sensitivity, specificity of each cluster
      #clustCells <- pull(SensSpe, clustersCells)
      for(i in as.vector(clustIDs[[ct]])){
        #print("ok")
        #print(ct)
        #print(i)
        clustCT <- hclustsNiches.cells%>%mutate(cell_type==ct)
        #print(clustCT[,c(i,"archPos")])
        counts <-table(clustCT[,c(i,"archPos")])
        #print(counts)
        sens <- counts['1','1']/sum(counts[,'1'])
        #print(sens)
        S1<- append(S1, sens) # sensitivity
        spe <- counts['0','0']/sum(counts[,'0'])
        S2 <- append(S2, spe) # specificity
      }
      print(auc(NichesPredClust%>%roc(response=archPos,predictor=value)))
      pdf(paste0("./figs/",ct,"_",ar,"ROC.pdf"),height=4,width=4);
      plot(NichesPredClust%>%roc(response=archPos,predictor=value),main=paste(ct,"in",ar));
      points(S2,S1,pch=19, col="red");#(spe,sens,pch=19,col="red")
      text(S2+.01,S1+.01,col="red",labels=as.vector(clustIDs[[ct]]))
      #legend(0.2,.2,legend=c("NIPMAP","hclust"),fill=c("black","red"));
      dev.off()
    }
  })
}
##### CORRECTION OF P VALUES ON SPEARMAN CORR TEST: FALSE DISCOVERY RATE
# Correleations between core/interface TMENs and marker expressed in cell type
correct_cor_fdr <- function(corMatrix,corCMtmens.pval,qThresh,corThresh){
  ### CORRECTIONS P VALUE FDR
  #dim(corCMtmens.pval)
  print("ok")
  #write_csv(as.data.frame(corMatrix,row.names=rownames(corMatrix)),"./NiPhcorrelations.csv")
  saveRDS(as.data.frame(corMatrix,row.names=rownames(corMatrix)),"./NiPhcorrs.rds")
  corCMtmens.pval.mat = corCMtmens.pval %>% column_to_rownames(var='names')
  fdrOut = 
    fdrtool(corCMtmens.pval.mat %>% as.matrix %>% as.numeric(), 
            statistic = "pvalue")
  #write_csv(as.data.frame(corCMtmens.pval.mat,row.names=rownames(corCMtmens.pval.mat)),"./NiPhpvalue.csv")
  saveRDS(as.data.frame(corCMtmens.pval.mat,row.names=rownames(corCMtmens.pval.mat)),"./NiPhpvalue.rds")
  corCMtmens.qval = 
    fdrOut$qval %>% 
    matrix(nrow=nrow(corCMtmens.pval.mat), ncol=ncol(corCMtmens.pval.mat),
           dimnames = list(rownames(corCMtmens.pval.mat), colnames(corCMtmens.pval.mat)))
  # dim(corCMtmens.qval)
  # dim(corCMtmens.pval.mat)
  # plot(corCMtmens.qval %>% as.matrix %>% as.numeric(), corCMtmens.pval.mat %>% as.matrix %>% as.numeric(), xlim=c(0,.1), ylim=c(0,.1))
  ## Qval and correlations thresholds
  # qThreshold <- 1/100
  # corThreshold <- .3
  # corThreshold2 <- .2
  #write_csv(as.data.frame(corCMtmens.qval,row.names=rownames(corCMtmens.qval)),"./NiPhqvalue.csv")
  saveRDS(as.data.frame(corCMtmens.qval,row.names=rownames(corCMtmens.qval)),"./NiPhqvalue.rds")
  filterMat <- corCMtmens.qval<qThresh & corMatrix[rownames(corCMtmens.qval),]>corThresh#corCMtmens.qval<qThresh & abs(corMatrix[rownames(corCMtmens.qval),])>corThresh
  return(filterMat)
}


##### COMPUTE SPEARMAN CORRELATIONS BETWEEN CELL PHENOTYPES AND TMENS 
# and filter out non-significant (FDR correction)/correlations under thresh
#FIXME corCMtmens is useless, why do we keep it ?
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
  #write_csv(corCMtmens,"./outputs/rawCorMatrix.csv")
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

#### Get interface and niches names from nb of niches and interface order
getInterNiches <- function(nIntf,nbNiches){
  interfaces <- combn(paste0("a",as.vector(seq(1,nbNiches,1))),nIntf)
  coreIntf <- apply(interfaces,2,function(x) paste0(x,collapse=""))
  return(coreIntf)
}

#### Compute correlations between niches (and interfaces) and cell phenotypes
# markersCells.niches : df with 
# Markers: string vector of functional markers that identifies cell phenotypes
# corrMeth: either "spearman" or "pearson", method to compute correlations
# coreIntf: str vector of niches and interfaces to assess correlations (archetypes names from markersCells.niches)
# qThresh: double, cut-off of q value for FDR correction of p values in correlation test
# corThresh: double, cut-off to select high correlation
# correlation_niches_CM <- function(markersCells.niches,Markers,corrMeth="spearman",coreIntf,qThresh=1/100,corThresh=0.3,nbNiches=NBNICHES,intf=NINTERFACE){
#   rareCells<- markersCells.niches%>%group_by(cell_type)%>%summarise(total=n())%>%filter(total<=3*length(Markers))%>%pull(cell_type)
#   if (length(rareCells)>=1){
#     print("These cell types are rare (count<10),")
#     print(rareCells)
#   }
#   #print("They will be removed from the correlation analysis")
#   CM.gps <- markersCells.niches%>%filter(!cell_type %in% rareCells)%>% #removing rare cell types
#     filter(marker %in% Markers)%>%
#     mutate(id=paste(cell_type,marker,sep=";"))%>%
#     group_by(cell_type, marker)%>%split(f= as.factor(.$id))
#   #print(CM.gps)
#   res <-lapply(CM.gps,function(x){ #a%>%group_by(cell_type, marker)%>%group_split()%>%setNames(id)
#     corrs <- c()
#     corrpvals <-c()
#     #id <- paste0(unique(pull(x,cell_type)),unique(pull(x,marker)),sep=";")
#     #print(id)
#     #cellPhen.names <- append(cellPhen.names,id)
#     for(n in coreIntf){
#       #print(n)
#       #print(length(pull(x,value)))
#       if(any(is.na(pull(x,value))|sd(pull(x,value))==0 |length(pull(x,value))<3)){#|sd(pull(x,value))==0){
#         #print("ok")
#         corrValue <- NA
#         corrp <- NA
#       }
#       # else{
#       #   corrValue <- cor.test(pull(x,value),pull(x,get(n)),method = corrMeth,exact=FALSE)$estimate
#       #   corrp <- cor.test(pull(x,value),pull(x,get(n)),method = corrMeth,exact=FALSE)$p.value 
#       # }
#       # corrs <-append(corrs,corrValue)
#       # corrpvals <- append(corrpvals,corrp)
#       else{
#         # Checking if a cell phenotype has at least one cell that is present in the niche (niche weigth>1/2)
#         #or in the interface(niche weight>1/8), if not, set correlation/p-value to 0 (no correlation)
#         if ((n %in%getInterNiches(intf,nbNiches) &(length(x%>%filter(get(n)>(1/8))%>%pull(get(n)))==0))|((n %in%paste0("a",as.vector(seq(1,nbNiches,1))) & (length(x%>%filter(get(n)>0.5)%>%pull(get(n)))==0)))){
#           #print(unique(x$id))
#           #print(n)
#           corrValue <- 0
#           corrp <- 0.5
#         }
#         else{
#           corrValue <- cor.test(pull(x,value),pull(x,get(n)),method = corrMeth,exact=FALSE)$estimate
#           corrp <- cor.test(pull(x,value),pull(x,get(n)),method = corrMeth,exact=FALSE)$p.value 
#         }
#       }
#       corrs <-append(corrs,corrValue)
#       corrpvals <- append(corrpvals,corrp)
#     }
#     names(corrs) <- coreIntf2
#     names(corrpvals) <- coreIntf2
#     list("corr_value"=corrs,"p_value"= corrpvals)
#   })
#   #print(names(res))
#   CorrMat <- do.call(rbind,lapply(res,'[[',"corr_value"))#do.call(rbind,res["corr_value"])
#   CorrpVal.mat <- do.call(rbind,lapply(res,'[[',"p_value"))#do.call(rbind,res[[2]])
#   #saveRDS(CorrMat,"./outputs/corMatrix_raw.rds")
#   corCMniches.pval <- CorrpVal.mat%>%as_tibble(rownames=NA)%>%
#     rownames_to_column(var="names")%>%
#     drop_na()
#   #write_csv(as_tibble(CorrMat,rownames=NA), "./rawCorrMatrix.csv")
#   saveRDS(CorrMat,"./corMatrix_raw.rds")
#   filtMat <- correct_cor_fdr(CorrMat,corCMniches.pval,qThresh,corThresh)
#   #mean(apply(filtMat, 1, sum) > 0)
#   cM.filt <- CorrMat[names(which(apply(filtMat, 1, sum) > 0)),]
#   return(cM.filt)
# }
# 
########-----CORRELATION+ FDR NEW VERSION ------###########
correlation_niches_CM <- function(markersCells.niches,Markers,corrMeth="spearman",coreIntf,qThresh=1/100,corThresh=0.3,nbNiches){
  rareCells<- markersCells.niches%>%group_by(cell_type)%>%summarise(total=n())%>%filter(total<=3*length(Markers))%>%pull(cell_type)
  # Remove rare cell types, not enough cells to compute correlations
  if (length(rareCells)>=1){
    print("These cell types are rare (count<10),")
    print(rareCells)
  }
  #Groups of all possible cell phenotypes-niches combinations
  CM.gps <- markersCells.niches%>%filter(!cell_type %in% rareCells)%>% #removing rare cell types
    filter(marker %in% Markers)%>%
    pivot_longer(cols=coreIntf,names_to = "niche",values_to = "weight")%>%
    mutate(id=paste(cell_type,marker,niche,sep=";"))%>%
    mutate(threshNiche= ifelse(niche%in%paste0("a",as.vector(seq(1,nbNiches,1))),1/2,1/8))%>% # determine cutoff for cells presence in niche (1/2) or interface(1/8)
    group_by(cell_type, marker,niche)%>%split(f= as.factor(.$id))
  #print(CM.gps)
  #print("ok")
  
  #Filter out cells phenotyes-niches that don't have at least one cell present in niche or interface
  CMgps.filt <- CM.gps[lapply(CM.gps,function(x) length(x%>%filter(weight>threshNiche)%>%pull(weight)))>1]
  #print(length(CM.gps))
  #print(length(CMgps.filt))
  
  #For each cell phenotypes-niches combination that match the criteria, compute correlation test
  res <-lapply(CMgps.filt,function(x){ #a%>%group_by(cell_type, marker)%>%group_split()%>%setNames(id)
    #corrs <- c()
    #corrpvals <-c()
    #id <- paste0(unique(pull(x,cell_type)),unique(pull(x,marker)),sep=";")
    #print(id)
    #cellPhen.names <- append(cellPhen.names,id)
    #Compute correlations and p-value test
    corrValue <- cor.test(pull(x,value),pull(x,weight),method = corrMeth,exact=FALSE)$estimate
    corrp <- cor.test(pull(x,value),pull(x,weight),method = corrMeth,exact=FALSE)$p.value 
    dfres <- data.frame(rho = corrValue,pvalue = corrp)
    rownames(dfres) <- unique(x$id)
    
    dfres
    
  })
  #print(names(res))
  #print(names(res))
  # Merge list of dataframes
  CM2 <- do.call(rbind,res)
  corrOut <- bind_rows(CM2, .id = "id")
  corMatrix <- as.matrix(corrOut%>%dplyr::select(rho))
  
  ##########--------------FDR CORRECTION--------------##########
  corCMtmens.pval.mat = corrOut%>%
    #dplyr::select(c(pvalue))%>%
    drop_na()#%>%
  #column_to_rownames(var="id")#corCMtmens.pval %>% column_to_rownames(var='names')
  fdrOut =fdrtool(as.matrix(corCMtmens.pval.mat%>%dplyr::select((pvalue)))%>%as.numeric,statistic = "pvalue")
  #write_csv(as.data.frame(corCMtmens.pval.mat,row.names=rownames(corCMtmens.pval.mat)),"./NiPhpvalue.csv")
  #saveRDS(as.data.frame(corCMtmens.pval.mat,row.names=rownames(corCMtmens.pval.mat)),"./NiPhpvalue.rds")
  corCMtmens.qval <-fdrOut$qval %>%
    matrix(nrow=nrow(corCMtmens.pval.mat), ncol=1,
           dimnames = list(rownames(corCMtmens.pval.mat), "qvalue"))
  #write_csv(as.data.frame(corCMtmens.qval,row.names=rownames(corCMtmens.qval)),"./NiPhqvalue.csv")
  #saveRDS(as.data.frame(corCMtmens.qval,row.names=rownames(corCMtmens.qval)),"./NiPhqvalue.rds")
  #hist(corCMtmens.qval)
  # Build q-value matrix for all niches/interfaces
  corQval <- corCMtmens.qval%>%as_tibble(rownames=NA)%>%
    rownames_to_column(var="id")%>%
    separate(col=id, into=c("ct","marker","niche"),sep=";")%>%
    mutate(id= paste0(ct,";",marker),.keep="unused")%>%#column_to_rownames(var="id")%>%
    pivot_wider(id_cols=id, names_from= niche,values_from=qvalue)%>%
    column_to_rownames(var="id")%>%
    replace(is.na(.),1) # Set q value to 1 (to bel ater filtered out) to NA values, these associations didn't match aforementioned criteria
  print(head(corQval))
  
  #Recover all possible niches and interfaces with the cell phenotypes
  if(ncol(corQval)!=length(coreIntf)){
    corQval[,coreIntf[!(coreIntf %in% colnames(corQval))]]<-0
  }
  
  ##Build correlation matrix for all niches/interfaces
  colnames(corMatrix)<-"correlation"
  corMatrix2 <- corMatrix%>%as_tibble(rownames=NA)%>%
    rownames_to_column(var="id")%>%
    separate(col=id, into=c("ct","marker","niche"),sep=";")%>%
    mutate(id= paste0(ct,";",marker),.keep="unused")%>%#column_to_rownames(var="id")%>%
    pivot_wider(id_cols=id, names_from= niche,values_from=correlation)%>%
    column_to_rownames(var="id")%>%
    replace(is.na(.),0)# replace NA correlations by 0 (these correlations didn't match aforementioned criteria)
  #Recover all possible niches and interfaces with the cell phenotypes
  if(ncol(corMatrix2 )!=length(coreIntf)){
    corMatrix2 [,coreIntf[!(coreIntf %in% colnames(corMatrix2 ))]]<-0
  }
  
  #print(dim(corMatrix2 ))
  ## Filter correlations by qvalue and correlation threshold
  filterMat <- corQval<qThresh & corMatrix2[rownames(corQval),]>corThresh#corCMtmens.qval<qThresh & corMatrix2[rownames(corCMtmens.qval),]>corThresh
  # Keep cell phenotypes that are associated with at least one niche/interface, association being determined by correlation value >0.3 , q value<1/100
  cM.filt <- corMatrix2[names(which(apply(filterMat, 1, sum) > 0)),]%>%as.matrix
  #cMf2 <- cMf%>%replace(is.na(.),0)
  
  return(cM.filt)#(corQval)#(cM.filt)
}


#########---- HEATMAP ORGANIZED BY CELL TYPES ----###########
# nichsIntf: string vector of archetypes and interfacs as named in correlation matrix
plot_heatmap_CT <- function(CM.mat,nichesIntf,figPath="./figs/cM_byCells3.pdf"){
  cts <- unique(pull(as_tibble(CM.mat,rownames=NA)%>%rownames_to_column(var="names")%>%separate(names,into=c("cell_type","marker"),sep=";"),cell_type))
  #print(cts)
  CM_TMENs_ct <- as_tibble(CM.mat,rownames=NA)%>%
    rownames_to_column(var="names")%>%
    separate(names,into=c("cell_type","marker"),sep=";")%>%#mutate(marker=paste0(marker,"+"))%>%
    pivot_longer(cols=all_of(nichesIntf),names_to = "region",values_to="corr_val")%>%
    mutate(idx =match(cell_type, cts))%>%
    group_by(cell_type)%>%mutate(corr_val = corr_val+ (idx-1)*10)%>%pivot_wider(names_from="region", values_from="corr_val")%>%mutate(names = paste(cell_type,marker,sep=";"))%>%column_to_rownames(var="names")#%>%dplyr::select(-c(cell_type,idx))
  #print(head(CM_TMENs_ct))
  dists <- dist(as.matrix(CM_TMENs_ct%>%dplyr::select(-c(cell_type,marker,idx))),method="euclidean")
  hclustCells <- hclust(dists,method="ward.D") 
  #plot(as.dendrogram(hclustCells))
  cellTypes <- data.frame(cell_type =pull(CM_TMENs_ct ,cell_type))
  Markers<- pull(CM_TMENs_ct ,marker)
  ## Annotations colors = length of cell types
  colorCount = length(cts)#length(CELLTYPES)
  getPalette = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
  getPalette(colorCount)
  # List with colors for each annotation.
  CTcolors <- list(cell_type = getPalette(colorCount))
  names(CTcolors$cell_type) <- cts#CELLTYPES
  rownames(cellTypes) <- rownames(CM_TMENs_ct)
  #rownames(Markers) <-rownames(CM_TMENs_ct)
  figWidth <- nrow(CM.mat) * 30/267
  pdf(figPath,height=8, width=ifelse(figWidth<6, 6, figWidth+2))
  pheatmap(t(CM.mat),cluster_cols = hclustCells,cluster_rows = FALSE,annotation_col=cellTypes,annotation_colors= CTcolors,labels_col = Markers)#
  dev.off()
}

#########---- HEATMAP ORGANIZED BY MARKERS ----###########
plot_heatmap_markers <- function(CM.mat,nichesIntf,figPath="./figs/cM_byMarkers2.pdf"){
  markersCorr <- unique(pull(as_tibble(CM.mat,rownames=NA)%>%rownames_to_column(var="names")%>%separate(names,into=c("cell_type","marker"),sep=";"),marker))
  CM_TMENs_ph <- as_tibble(CM.mat,rownames=NA)%>%
    rownames_to_column(var="names")%>%
    separate(names,into=c("cell_type","marker"),sep=";")%>%#mutate(marker=paste0(marker,"+"))%>%
    pivot_longer(cols=all_of(nichesIntf),names_to = "region",values_to="corr_val")%>%
    mutate(idx =match(marker, markersCorr))%>%arrange(marker)%>%
    group_by(marker)%>%mutate(corr_val = corr_val+ (idx-1)*10)%>%pivot_wider(names_from="region", values_from="corr_val")%>%mutate(names = paste(marker,cell_type,sep=";"))%>%column_to_rownames(var="names")#%>%dplyr::select(-c(marker,idx))
  
  Markers2 <- data.frame(marker= pull(CM_TMENs_ph,marker))
  rownames(Markers2) <- rownames(CM_TMENs_ph)
  cellTypes2 <- pull(CM_TMENs_ph,cell_type) #data.frame(cell_type =pull(CM_TMENs_ph ,cell_type))
  
  colorCount = length(unique(pull(CM_TMENs_ph,marker)))
  getPalette = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
  getPalette(colorCount)
  CTcolors <- list(marker = getPalette(colorCount))
  names(CTcolors$marker) <- unique(pull(CM_TMENs_ph,marker))
  #rownames(cellTypes2) <- rownames(CM_TMENs_ph)
  
  dists2 <- dist(as.matrix(CM_TMENs_ph%>%dplyr::select(-c(marker,cell_type,idx))),method="euclidean")
  hclustPhen <- hclust(dists2,method="ward.D") 
  #plot(as.dendrogram(hclustPhen))
  cMf.copy <- CM.mat
  rownames(cMf.copy) <- unname(unlist(sapply(rownames(cMf.copy),function(x){paste0(strsplit(x,";")[[1]][2],";",strsplit(x,";")[[1]][1])})))#sub('(\\w|\\.|\\_|\\ / |+);(\\w|.|_| / |+)', '\\2;\\1', rownames(cMf))#sub('([^;]);([^;])', '\\2;\\1^', rownames(cMf))
  # 
  cMf_ordered <- cMf.copy%>%as_tibble(rownames=NA)%>%rownames_to_column(var="names")%>%separate(col="names",into=c("marker","cell_type"),sep=";")%>%group_by(marker)%>%arrange(marker)%>%mutate(names=paste0(marker,";",cell_type))%>%column_to_rownames(var="names")%>%dplyr::select(-c(cell_type,marker))
  figWidth <- nrow(CM.mat) * 30/267
  pdf(figPath,height=8, width=ifelse(figWidth<6, 6, figWidth+2))
  pheatmap(as.matrix(cMf_ordered%>%t),cluster_cols = hclustPhen,cluster_rows = FALSE,annotation_col = Markers2,labels_col = cellTypes2,annotation_colors=CTcolors)
  dev.off()
}



#### Enriched cell types in each site 
# choose the 1% sites closest to the archetype
# compute the SD for each cell type
# compare the average density  d1 of a cell type from the top 1% closest sites with the av denisty of cell type from the rest of the sites d2
# if d1 > d2+ k* SD ==> the cell type is representative  of the archetype
# do the same for the other archetypes
# @param cellAb.df: df of cell abundance of sites and their archetypes weights
#                   with columns of cell_type and cell_density and
#                   archetypes as columns and sites as rows(with their ids)
# @param archetype: str of nae of archetype
# @param k: int, value to multiply standard deviation of cell type abundance
# @param thresh: double, threshold value to select top sites(between 0 and 1)
# @return Archs.CT: df of niches and their enriched cell types
get_cell_type_enriched_arch <- function(cellAb.df,archetype,k=1,thresh=0.99){
  
  # sitesNiche <- cellAb.df %>%filter(.data[[archetype]]>= cutOff)
  # print(dim(sitesNiche))
  # sites_a <- sitesNiche
  sites_a <- cellAb.df%>%filter(.data[[archetype]] > quantile(sort(pull(cellAb.df,.data[[archetype]]),decreasing=TRUE),thresh))
  #print(sites_a)
  summaryA <- sites_a %>%
    group_by(cell_type)%>%
    summarise(meanCT=mean(cell_density),sdCT = sd(cell_density))%>%
    arrange(desc(meanCT))#%>% ## arrange by decreasing order of cell density ==> most abundant cell types of this niche

  summaryA2 <- summaryA%>%column_to_rownames(var="cell_type")
  cts <- c()
  for(c in summaryA%>%pull(cell_type)){ #(c in CELLTYPES)
    if(as.numeric(summaryA2[c,"meanCT"]) > k *as.numeric(summaryA2[c,"sdCT"])) {#if(as.numeric(summaryA[c,"meanCT"]) > (as.numeric(summaryNA[c,"meanCT"]) + 1*as.numeric(summaryNA[c,"sdCT"]))){
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
# get_CT_enriched_all_archs <- function(archsSites, cellAbundSites,thresh=0.99){
#   archs <- colnames(archsSites)[grepl("a", colnames(archsSites))]# vector of archetypes
#   # list of cell types enriched in each archetype
#   #archetypes.ct <- list() 
#   #names(archetypes.ct) <- archs
#   archetypes_cts.all <- data.frame(niche=NULL, cell_type=NULL)
#   archsSitesCellAb <- left_join(archsSites, cellAbundSites,by=c("site_id","patient_id"))%>%pivot_longer(cols=append(CELLTYPES,"Unidentified"),names_to="cell_type",values_to= "cell_density")%>%ungroup()
#   for (a in archs){
#     res <- get_cell_type_enriched_arch(cellAb.df=archsSitesCellAb,archetype=a)
#     archetypes_cts.all<- rbind(archetypes_cts.all,res)#archetypes.ct[[a]] <-res
#   }
#   archetypes_cts.all$niche <-  set_tmens_names(archetypes_cts.all$niche)
#   return(archetypes_cts.all)#(archetypes.ct)
# }

# @param archsSitesCellAb: df with niches weights (from a1 to a1 in columns)+ column of cell_type and cell_density
# @param archNames: named vector where names are archetypes indexes (a1 to an) and values are the explicit names of niches
# @param thresh  thresh to select the closest sites to a niche
get_CT_enriched_all_archs <- function(archsSitesCellAb,archNames,thresh=0.99){
  # list of cell types enriched in each archetype
  #archetypes.ct <- list() 
  #names(archetypes.ct) <- archs
  archs <- names(archNames)
  archetypes_cts.all <- data.frame(niche=NULL, cell_type=NULL)
  for (a in archs){
    res <- get_cell_type_enriched_arch(cellAb.df=archsSitesCellAb,archetype=a)
    archetypes_cts.all<- rbind(archetypes_cts.all,res)#archetypes.ct[[a]] <-res
  }
  archetypes_cts.all$niche <-  str_replace_all(archetypes_cts.all$niche,archNames)#set_tmens_names(archetypes_cts.all$niche)
  return(archetypes_cts.all)#(archetypes.ct)
}


#FIXME unitary tests + exception/error handling
# Plots cell type profile of niches found by Archetype Analysis
# @param sitesCA: df of sites(rows) cell types(columns) abundance
# @param Arch: ArchetypalAnalysis object after AA on PCA onf sites cell abundance
# @param pcaObj: PCA object of sites cell abundance
# @return barplot of archetypes cell type composition ==> niches
NichesCellProfile <- function(sitesCA = sites,Arch = AA,pcaObj = pca_obj,nbNiches=NBNICHES){
  NichesCellProf <- get_niches_cell_abund(sitesCellAb = sitesCA,pcaSites = pcaObj,ArchObj = Arch,nComp = as.integer(nbNiches-1))
  colnames(NichesCellProf) <- CELLTYPES
  rownames(NichesCellProf) <- paste0("arch",seq(1,nbNiches))
  
  NichesCellProp <- NichesCellProf%>%as_tibble(rownames = NA)%>%
    rownames_to_column(var="archetype")%>%
    pivot_longer(cols=all_of(CELLTYPES),names_to="cell_type",values_to = "cell_density")
  NichesCellProp[NichesCellProp<0] <-0
  barplot1 <- ggplot(data = NichesCellProp, aes(x = cell_type, y = cell_density,fill = archetype)) +
    geom_bar(stat = "identity",position = position_dodge(),width = 0.6) + 
    theme(axis.text.x = element_text(angle = 90, vjust = .2))+
    scale_fill_manual(values = COLARCHS)+
    xlab ("") + ylab("cell density")
  #barplot1
  ggsave("./barplotNiches.pdf",barplot1,height=3,width=4)
  return(barplot1)
}


#FIXME unitary tests + exception/error handling
# Creates table of cell phenotypes/cell types associated with niches and interfaces between two niches
# @param CM: matrix of correlation between cell phenotypes(rows) and niches/interfaces(columns)
# @param NichesCT: dataframe of cell types enriched in each niche columns: niche, cell_type
# @param NichesNames: named vector of archetypes(names) and niches names(value)
# @param  nichesCA.sorted: dataframe of cell abundance of niches, cell types sorted by desc order
# @param pathFigs: str of path where to save table figure (pdf)
# @return tabCellPhen3: dataframe of cell phenotypes of each niche/interface
#test str_detect(names(c("a1"="TLS","a3"="cancer","a2"="inflammatory","a!"=0)),"^a[:digit:]")
TableNichesPhenotypes <- function(CM,NichesCT,Niches.names,nichesCA.sorted,pathFigs){
  intfNames <-apply(combn(Niches.names,2),2,function(x) {paste0(x,collapse = " x ")})
  names(intfNames) <-apply(combn(names(Niches.names),2),2,function(x) {paste0(x,collapse = "")})
  
  tabTMENs2 <- as_tibble(CM,rownames=NA)%>%
    rownames_to_column(var="names")%>%
    separate(names,into=c("cell_type","marker"),sep=";")%>%
    mutate(marker=paste0(marker,"+"))%>%
    pivot_longer(cols=all_of(coreIntf2),names_to = "niche",values_to="corr_val")#%>%
  #mutate(region =)%>%#(region = set_tmens_names(region))%>%
  tabTMENs2$niche <- str_replace_all(pull(tabTMENs2,niche),intfNames)
  tabTMENs2$niche <- str_replace_all(pull(tabTMENs2,niche),Niches.names)
  
  tabTMENs2 <- tabTMENs2%>% filter(corr_val >=0.3)%>%
    #filter(!(grepl("cancer",niche,fixed=TRUE) & (marker%in% c("Keratin6+","Beta catenin+"))))%>%
    left_join(.,nichesCA.sorted, by=c("cell_type","niche"))%>%group_by(niche)%>%
    arrange(desc(cell_density))%>%dplyr::select(-c(cell_density))#### Select correlations that are >0.35 + remove correlations that arise from bleed over of structural proteins in cancer and interfaces with cancer.
  
  ######### Compare correlations from core VS from interfaces
  # group by cell type, marker
  #find rows with a common determinant: a niche
  # split region into two (for 2- interfaces)
  # find intersect between two subregions ==> common TMEN
  # OR find the niche that has the highest occurence in a subregion (core and interfacees) across same CM VStmens associations
  # OR the same but for the second subregion
  # IF there is an association between CM and one region ==> common region is the latter
  # group again by cell type, marker, and by this common region
  # select the maximum correlation value rho
  # Get correlations among core/interfaces for a given cell type/marker that are the highest
  # Ex: rho (Beta-cat+ Tregs; cancer) < rho( Beta-cat+ Tregs; cancer x fibrotic) 
  # ==> select rho( Beta-cat+ Tregs, cancer x fibrotic)
  #namesIntf <- paste0(rep("type",InftOrder),as.character(as.vector(seq(1,InftOrder))))
  
  tabCellPhen <- tabTMENs2%>%
    separate(col=niche, into = c("type1","type2"),sep="x")%>%
    group_by(cell_type,marker)%>%
    mutate(len= length(type1))%>%
    mutate(inter_region=ifelse(len==1,type1,
                               ifelse(head(sort(table(type1),decreasing=TRUE),1)>1, names(head(sort(table(type1),decreasing=TRUE),1)),
                                      ifelse(head(sort(table(type2),decreasing=TRUE),1)>1,  names(head(sort(table(type2),decreasing=TRUE),1)),
                                             ifelse(length(intersect(type1,type2))==0,paste(type1,type2,sep=" x "),intersect(type1,type2) )))))%>%
    group_by(cell_type, marker,inter_region)%>%filter(corr_val==max(corr_val))
  
  tabCellPhen2 <- tabCellPhen %>%
    mutate(niche = paste(type1,type2,sep=" x "))%>%
    #full_join(.,archs.CT,by=("region"))%>%
    mutate(niche=str_replace_all(niche," x NA",""))%>%
    filter(!(cell_type %in% c("DC / Mono", "CD3-T", "Mono / Neu", "Other immune"))) %>%
    group_by(niche)%>%mutate(cells = paste(unique(cell_type),collapse='\n'))%>%
    group_by(niche,cell_type)%>% 
    mutate(cell_phenotype= paste0(paste(marker, collapse = " "),cell_type)) %>%
    ungroup()%>%
    group_by(niche)%>%
    mutate(cell_phenotype=paste(unique(cell_phenotype),collapse="\n")) %>%
    distinct(niche, .keep_all = TRUE)%>%
    arrange(niche)%>%dplyr::select(niche, cells,cell_phenotype)
  
  tabCellPhen3 <- tabCellPhen2%>%
    mutate(cellPh =ifelse(niche %in% pull(NichesCT,niche),paste(NichesCT[NichesCT$niche==niche ,"cell_type"],collapse="\n"),""))%>%
    rowwise()%>%
    mutate(ct = ifelse(grepl("\n",cellPh,fixed=TRUE)==FALSE,paste(setdiff(as.vector(cellPh),str_split(cells,"\n")[[1]]),collapse="\n"),ifelse(grepl("\n",cells,fixed=TRUE)==FALSE,paste(setdiff(str_split(cellPh,"\n")[[1]], as.vector(cells)),collapse="\n"),paste(setdiff(str_split(cellPh,"\n")[[1]],str_split(cells,"\n")[[1]]),collapse="\n"))))%>%
    mutate(cell_phenotype = paste(cell_phenotype,ct,sep="\n"),.keep="unused")%>%
    dplyr::select(-c(cells, cellPh))%>%#mutate(sign  = "+")
    #distinct(region, .keep_all = TRUE)%>%
    arrange(niche)
  
  ### Display and save table in pdf
  th1 <- ttheme_default()
  g1 <- tableGrob(tabCellPhen3,rows=NULL,theme = th1)#(tabCellPhen2%>%mutate(region=str_replace_all(region," x NA","")),rows=NULL,theme = th1)
  #gpar.corefill = gpar(fill = "white", col = "grey95"))
  grid.newpage()
  grid.draw(g1)
  ggsave(paste0(pathFigs,"/tabNichePhen.pdf"),g1,width=7,height=10)
  
  return(tabCellPhen3)
  
}









