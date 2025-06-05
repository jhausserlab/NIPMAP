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


file1 = "./pca_sites.json" # pca object on sites elements
file2 = "./AA_sites.json" # archetype Analysis object based on sites cell abundance
file3 = "./ca_sites.json" # cell abundance of randomly generated sites
# file4 = "./cells_niches.json" # sites centered on cells and niches weights

#######---- Open .json files ----#######
json_data <- fromJSON(file=file1)
json_data2 <- fromJSON(file=file2)
json_data3 <- fromJSON(file=file3)
#json_data4 <- fromJSON(file=file4)


##### LOAD OUTPUT OBJECTS
## Cell abundance in sites

sitesCellAb <- as_tibble(lapply(json_data3$cellAbSites,unlist))
write_csv(sitesCellAb%>%dplyr::select(-c(index, patient_id,site_id)),"sitesCA.csv")

niches <- paste0("a",as.vector(seq(1,NBNICHES,1)))
#niches <- c('cancer', 'inflammatory', 'B follicle', 'Other', 'Low density')
names(COLNICHES) <- niches

colNiches.hex <-unlist(lapply(COLNICHES, function(x){rgb(x[1],x[2],x[3],maxColorValue = 255)}))
## Archetypes coordinates in reduced PC space
Archs_3D <- do.call(cbind,lapply(json_data2$archs_coord,unlist))


## Projection of sites cell abundance in reduced PC space
pca3D <- matrix(unlist(json_data$PC_proj),nrow=length(CELLTYPES))[1:3,]
plotly::plot_ly(x=pca3D[1,],
                y=pca3D[2,],
                z=pca3D[3,],
                type = "scatter3d", mode = "markers",
                marker = list(symbol = "triangle",size = 1),
                name="sites",
                mode = "text")%>%
  plotly::add_trace(x = Archs_3D[,1],
                    y = Archs_3D[,2],
                    z =Archs_3D[,3],
                    type = "scatter3d",
                    mode = "markers+text",
                    text = niches,
                    # textposition = c('top right','bottom right','top left','top right'),
                    textfont = list(color = '#000000', size = 16),
                    showlegend = TRUE,
                    name = "niches",
                    marker = list(color=~colNiches.hex,symbol = "star-diamond",size = 12),
                    inherit = FALSE)%>%
  plotly::layout(scene = list(xaxis = list(title = "PC1"),
                              yaxis = list(title = "PC2"),
                              zaxis = list(title = "PC3")))

######--- NICHE IDENTIFICATION 

NichesCellProf <- do.call(cbind,lapply(json_data2$nichesCA,unlist))
rownames(NichesCellProf) <- CELLTYPES
colnames(NichesCellProf) <- niches
NichesCellProp <- NichesCellProf%>%t%>%as_tibble(rownames = NA)%>%
  rownames_to_column(var="archetype")%>%
  pivot_longer(cols=all_of(CELLTYPES),names_to="cell_type",values_to = "cell_density")
NichesCellProp[NichesCellProp<0] <-0
# write_csv(NichesCellProp,"./bad_segmentation/nicheComposition_12.csv")

barplot1 <- ggplot(data = NichesCellProp, aes(x = cell_type, y = cell_density,fill = archetype)) +
  geom_bar(stat = "identity",position = position_dodge(),width = 0.6) +
  scale_fill_manual(values = colNiches.hex)+
  theme(axis.text.x = element_text(angle = 90, vjust = .2))#+
#xlab ("") + ylab("cell density")
ggsave("./barplotNiches.pdf",barplot1,height=3,width=4)


# save(Archs_3D, file = "./bad_segmentation/with_others/Archs_3D_4.RData")
# save(pca3D, file = "./bad_segmentation/with_others/pca3D_4.RData")



######--- Define the NICHE names
niches <- c('inflammatory', 'cancer', 'B follicle', 'Other', 'Low density')
names(COLNICHES) <- niches

colNiches.hex <-unlist(lapply(COLNICHES, function(x){rgb(x[1],x[2],x[3],maxColorValue = 255)}))
## Archetypes coordinates in reduced PC space
Archs_3D <- do.call(cbind,lapply(json_data2$archs_coord,unlist))


## Projection of sites cell abundance in reduced PC space
pca3D <- matrix(unlist(json_data$PC_proj),nrow=length(CELLTYPES))[1:3,]
plotly::plot_ly(x=pca3D[1,],
                y=pca3D[2,],
                z=pca3D[3,],
                type = "scatter3d", mode = "markers",
                marker = list(symbol = "triangle",size = 1),
                name="sites",
                mode = "text")%>%
  plotly::add_trace(x = Archs_3D[,1],
                    y = Archs_3D[,2],
                    z =Archs_3D[,3],
                    type = "scatter3d",
                    mode = "markers+text",
                    text = niches,
                    # textposition = c('top right','bottom right','top left','top right'),
                    textfont = list(color = '#000000', size = 16),
                    showlegend = TRUE,
                    name = "niches",
                    marker = list(color=~colNiches.hex,symbol = "star-diamond",size = 12),
                    inherit = FALSE)%>%
  plotly::layout(scene = list(xaxis = list(title = "PC1"),
                              yaxis = list(title = "PC2"),
                              zaxis = list(title = "PC3")))


NichesCellProf <- do.call(cbind,lapply(json_data2$nichesCA,unlist))
rownames(NichesCellProf) <- CELLTYPES
colnames(NichesCellProf) <- niches
NichesCellProp <- NichesCellProf%>%t%>%as_tibble(rownames = NA)%>%
  rownames_to_column(var="archetype")%>%
  pivot_longer(cols=all_of(CELLTYPES),names_to="cell_type",values_to = "cell_density")
NichesCellProp[NichesCellProp<0] <-0
# write_csv(NichesCellProp,"./bad_segmentation/nicheComposition_12.csv")

barplot1 <- ggplot(data = NichesCellProp, aes(x = cell_type, y = cell_density,fill = archetype)) +
  geom_bar(stat = "identity",position = position_dodge(),width = 0.6) +
  scale_fill_manual(values = colNiches.hex)+
  theme(axis.text.x = element_text(angle = 90, vjust = .2))#+
#xlab ("") + ylab("cell density")
ggsave("./barplotNiches.pdf",barplot1,height=3,width=4)

