gc()
rm(list=ls())
#.libPaths("/scratch/anissa.el/R_old/x86_64-redhat-linux-gnu-library/4.0")
.libPaths("/home/common/R")
library(rjson)
library(tidyverse)
library(purrr)

### SET WORKING DIRECTORY
dirName <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dirName)
source("./phenotypes_niches/functions_phenotypes_tmens.r")
#Master\\\ File
jsonparams <- fromJSON(file="./params.json")
CELLTYPES <-jsonparams$cellTypes
ImageIDs <- jsonparams$ImageID
NSITES <- jsonparams$nbsites
RADIUS <- jsonparams$radiusSize
NBNICHES <- jsonparams$nbniches
METHOD <-jsonparams$countMeth
W <-jsonparams$xsize
H <-jsonparams$ysize
ROOT_DATA_PATH <- jsonparams$rootDataPath
ROOT_OUTPUT_PATH <-jsonparams$rootOutPath
COLNICHES <- jsonparams$colNiches
pathFigs <- jsonparams$pathFigs


#file1 = "./pca_sites_2niches.json" # pca object on sites elements
file2 = "./AA_sites.json" # archetype Analysis object based on sites cell abundance
#file3 = "./ca_sites_2niches.json" # cell abundance of randomly generated sites
file4 = "./cells_niches.json" # sites centered on cells and niches weights

#######---- Open .json files ----#######
#json_data <- fromJSON(file=file1)
json_data2 <- fromJSON(file=file2)
#json_data3 <- fromJSON(file=file3)
json_data4 <- fromJSON(file=file4)

##### LOAD OUTPUT OBJECTS
## Cell abundance in sites
# sitesCellAb <- as_tibble(lapply(json_data3$cellAbSites,unlist))
# write_csv(sitesCellAb%>%dplyr::select(-c(index, patient_id,site_id)),"sitesCA.csv")

niches <- paste0("a",as.vector(seq(1,NBNICHES,1)))
names(COLNICHES) <- niches
colNiches.hex <-unlist(lapply(COLNICHES, function(x){rgb(x[1],x[2],x[3],maxColorValue = 255)}))

# ## Archetypes coordinates in reduced PC space
# Archs_3D <- do.call(cbind,lapply(json_data2$archs_coord,unlist))
# ## Projection of sites cell abundance in reduced PC space
# pca3D <- matrix(unlist(json_data$PC_proj),nrow=17)[1:3,] #length(CELLTYPES)
# plotly::plot_ly(x=pca3D[1,],
#                 y=pca3D[2,],
#                 z=pca3D[3,],
#                 type = "scatter3d", mode = "markers",
#                 marker = list(symbol = "triangle",size = 4),
#                 name="sites",
#                 mode = "text")%>%
#   add_trace(x = Archs_3D[,1],
#             y = Archs_3D[,2],
#             z =Archs_3D[,3],
#             type = "scatter3d",
#             mode = "markers+text",
#             text = niches,
#             textposition = c('top right','bottom right','top left','top right'),
#             textfont = list(color = '#000000', size = 16),
#             showlegend = TRUE,
#             name = "niches",
#             marker = list(color=~colNiches.hex,symbol = "star-diamond",size = 12),
#             inherit = FALSE)%>%
#   layout(scene = list(xaxis = list(title = "PC1"),
#                       yaxis = list(title = "PC2"),
#                       zaxis = list(title = "PC3")))


# arche_list <- c()
# for (i in seq(1,NBNICHES)) {
#   arche_list <- arche_list %>% append(paste("arch",i) %>% str_replace_all(" ",""))
# }
# 
# #######---- The pc composition of all niches ----#######
# 
# niche_pccompo <- read.csv(paste(RESULT_DATA_PATH,"niche_pcabundance.csv") %>% str_replace_all(" ","")) %>% 
#   dplyr::select(-c(X)) %>% 
#   cbind(arche_list = arche_list) %>% 
#   column_to_rownames(var = "arche_list")
# niche_alfas <- matrix(unlist(file2$alfas), ncol=length(arche_list)) %>% as.data.frame()
# colnames(niche_alfas) <- arche_list
# 
# plist <- Faceprojection(arche_list=, niche_pccompo, niche_alfas)
# plot <- do.call("grid.arrange", c(plist, ncol = 4))
# ggsave("./face_projection_niches.pdf", plot, height=10,width=12)

######--- NICHE IDENTIFICATION 

NichesCellProf <- do.call(cbind,lapply(json_data2$nichesCA,unlist))
rownames(NichesCellProf) <- CELLTYPES
colnames(NichesCellProf) <- niches
NichesCellProp <- NichesCellProf%>%t%>%as_tibble(rownames = NA)%>%
  rownames_to_column(var="archetype")%>%
  pivot_longer(cols=all_of(CELLTYPES),names_to="cell_type",values_to = "cell_density")
NichesCellProp[NichesCellProp<0] <-0
barplot1 <- ggplot(data = NichesCellProp, aes(x = cell_type, y = cell_density,fill = archetype)) +
  geom_bar(stat = "identity",position = position_dodge(),width = 0.6) +
  scale_fill_manual(values = colNiches.hex)+
  theme(axis.text.x = element_text(angle = 90, vjust = .2))#+
  #xlab ("") + ylab("cell density")
ggsave("./barplot_niches.pdf",barplot1,height=3,width=4)

##########--------- NICHE-PHENOTYPE MAPPING --------##########
## Niches weights(proportions) of all cells from all images
niches<- paste0("a",as.vector(seq(1,NBNICHES,1)))
NINTERFACES <- 2
MARKERS <- c("dsDNA","Vimentin","SMA","B7H3","FoxP3","Lag3","CD4","CD16","CD56",
             "OX40","PD1","CD31","PD-L1","EGFR","Ki67","CD209","CD11c","CD138",
             "CD163","CD68","CSF-1R","CD8","CD3","IDO","Keratin17","CD63","CD45RO",
             "CD20","p53","Beta catenin","HLA-DR","CD11b","CD45","H3K9ac","Pan-Keratin",
             "H3K27me3","H3K9ac_H3K27me3","phospho-S6","MPO","Keratin6","HLA_Class_1")

markers <- read.csv("./phenotypes_niches/data/proteins_by_frame.csv")%>%filter(Purpose=="Functional")%>%pull(Biomarker)
cellsNiches <- as_tibble(lapply(json_data4$cells_niches,unlist))%>%
  rename_at(vars(matches("[0-9]")),~niches)%>%
  mutate(cell_id=as.numeric(cell_id))
#Make combinations of niches interfaces of order nIntf
getInterNiches <- function(nIntf,nbNiches){
  interfaces <- combn(paste0("a",as.vector(seq(1,nbNiches,1))),nIntf)
  coreIntf <- apply(interfaces,2,function(x) paste0(x,collapse=""))
  return(coreIntf)
}
coreIntf2 <- append(paste0("a",as.vector(seq(1,NBNICHES,1))),getInterNiches(NINTERFACES,NBNICHES))

# Get markers expression and niche weigths of cells
# df with following columns: cell_type, SampleID, cell_id, a1....an & interfaces, marker,value
cellsPhen.niches <- read.csv("./TMENS_analysis/data/cellData.csv",check.names=FALSE,header = TRUE, sep =',')%>%
  dplyr::select(-c(cellSize,C,Na,Si,P,Ca,Fe,immuneCluster,Ta,Au))%>%
  mutate(immuneGroup = recode(immuneGroup,`0`= 'None',`1`='Tregs', `2`='CD4-T',
                              `3`='CD8-T', `4`='CD3-T', `5`='NK',
                              `6`='B', `7`='Neutrophils', `8`='Macrophages', `9`='DC',
                              `10`='DC / Mono', `11`='Mono / Neu', `12`='Other immune')) %>%
  mutate(Group = recode(Group,`1`='Unidentified', `2`='Immune',
                        `3`='Endothelial', `4`='Mesenchymal-like',
                        `5` = 'Tumor',
                        `6` = 'Keratin-positive tumor'))%>%
  mutate(cell_type = ifelse(Group == 'Immune', cell_type<- immuneGroup,cell_type <- Group))%>%
  dplyr::select(-c(tumorYN,tumorCluster,Group,immuneGroup))%>%filter(cell_type!="Unidentified")%>%
  #dplyr::rename(patient_id = SampleID)%>%
  dplyr::rename(cell_id = cellLabelInImage)%>%
  mutate(H3K9ac_H3K27me3 = H3K9ac/H3K27me3)%>%
  left_join(cellsNiches%>%filter(cell_type!="Unidentified"),.,by=c("SampleID","cell_id","cell_type"))%>%
  filter(!(is.na(a1)| is.na(a2)| is.na(a3) | is.na(a4)))%>%
  mutate(a1a2 = a1*a2)%>%
  mutate(a1a3 = a1*a3)%>%
  mutate(a1a4 = a1*a4)%>%
  mutate(a2a3 = a2*a3)%>%
  mutate(a2a4 = a2*a4)%>%
  mutate(a3a4 = a3*a4)%>%
  mutate(a1a2a3 = a1*a2*a3)%>%
  mutate(a1a2a4 = a1*a2*a4)%>%
  mutate(a1a3a4 = a1*a3*a4)%>%
  mutate(a2a3a4 = a2*a3*a4)%>%
  mutate(a1a2a3a4 = a1*a2*a3*a4)%>%
  pivot_longer(cols=all_of(MARKERS),
               names_to="marker",values_to="value")

source("./phenotypes_niches/functions_phenotypes_tmens.r")
CM <- correlation_niches_CM(markersCells.niches=cellsPhen.niches,Markers=markers,corrMeth="spearman",coreIntf2,1/100,0.3,nbNiches=NBNICHES)

#Filter the cell phenotypes associated with cancer niche, contaminated by Keratin6 and Beta-catenin
Krt.filt <- names(which(CM[rownames(CM)[which(grepl("Keratin6",rownames(CM),fixed=TRUE))],"a3"]>0.3))
BCat.filt <- names(which(CM[rownames(CM)[which(grepl("Beta catenin",rownames(CM),fixed=TRUE))],"a3"]>0.2))

cMf <- CM[!rownames(CM)%in%BCat.filt& !rownames(CM)%in% Krt.filt,]

############# PLOT HEATMAPS AND TABLE OF SPATIAL CELL PHENOTYPES 
plot_heatmap_CT(CM.mat=cMf,coreIntf2,paste0(pathFigs,"/CMbyCells2.pdf"))
plot_heatmap_markers(CM.mat=cMf,coreIntf2,paste0(pathFigs,"/CMbyMarkers.pdf"))


##TABLE OF NICHE-ASSOCIATED CELL PHENOTYPES
#Archetype weights of randomly generated sites over TNBC patients images
archetypes_sites <- as.data.frame(do.call(cbind,lapply(json_data2$alfas,unlist)))

NichesNames <- c("a1"="TLS","a2"="inflammatory","a3"="cancer","a4"="fibrotic")
nichesCellAb <- NichesCellProf%>%t%>%as_tibble(rownames=NA)%>%rownames_to_column(var="archetype")%>%
  mutate(archetype=str_replace_all(archetype,NichesNames))%>%#str_replace_all(archetype,NichesNames))%>%
  column_to_rownames(var="archetype")%>%t()%>%
  as_tibble(rownames=NA)%>%rownames_to_column(var="cell_type")
# #Sort cell types per niche by decreasing cell abundance
nichesCA.sort <- nichesCellAb%>%
   pivot_longer(cols = as.vector(NichesNames),names_to="niche",values_to="cell_density")%>%
   group_by(niche)%>%arrange(desc(cell_density))

colnames(archetypes_sites) <- niches
archsSitesCellAb <- cbind(sitesCellAb,archetypes_sites)%>%
  pivot_longer(cols=append(CELLTYPES,"Unidentified"),names_to="cell_type",values_to= "cell_density")

### Get cell types enriched in each niche (take top 1% sites closest to niche)
archs.CT <- get_CT_enriched_all_archs(archsSitesCellAb,NichesNames)%>%#(archetypes_sites%>%t%>%as_tibble(rownames=NA),cellAbSites,thresh = 0.99)%>%
  group_by(niche)%>%
  filter(!(cell_type %in% c("DC / Mono","CD3-T", "Mono / Neu", "Other immune","Unidentified")))%>%
  mutate(cell_type = paste(unique(cell_type),collapse="\n"))%>%
  distinct(niche,cell_type)


source("./phenotypes_niches/functions_phenotypes_tmens.r")
## Get table of niches/interfaces-associated cell phenotypes
TableNichesPhenotypes(CM = cMf,NichesCT = archs.CT,Niches.names = NichesNames,nichesCA.sorted = nichesCA.sort,pathFigs = pathFigs)



