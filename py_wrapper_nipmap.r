.libPaths("/scratch/anissa.el/R_old/x86_64-redhat-linux-gnu-library/4.0")
library(rjson)
library(tidyverse)
library(purrr)

### SET WORKING DIRECTORY
dirName <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dirName)
#Master\\\ File

#TODO add input: named vectors of colors codes(html) for visualizing niches
#TODO add source code: functions_phenotypes_tmens.r ==> functions_nipmap.r 

#TODO plot correlations heatmaps organized by cell types & by markers
#TODO plot table of cell phenotypes for each niche/interface

### CONTROL PANEL: SET PARAMETERS AND CONSTANTS
CELLTYPES = c('CD8-T', 'Other\\\ immune', 'DC\\\ /\\\ Mono', 'CD3-T', 'B', 'NK', 'Keratin-positive\\\ tumor', 'Tumor', 
              'CD4-T', 'Mesenchymal-like', 'Macrophages', 'Endothelial', 'Tregs', 'Unidentified', 'DC', 'Mono\\\ /\\\ Neu', 
              'Neutrophils')
ImageIDs <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
              20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32, 33, 34, 35, 36, 37,
              38, 39, 40, 41)
NSITES=100 # number of sites generated per image
RADIUS= 25# radius of site in micrometer²
NBNICHES = 4 # number of niches to find (in PC space of NBNICHES-1 dimensions)
METHOD ="gaussian"
#Method = "gaussian"
ROOT_DATA_PATH="./TMENS_analysis/data/cell_positions_data" 
#rootDataPath = 
ROOT_OUTPUT_PATH="./TMENS_analysis/output"

CellTypes=paste(CELLTYPES,collapse=",")# 
imageID=paste(ImageIDs,collapse = ",")
nsites = as.character(NSITES)
Radius = as.character(RADIUS)
nbNiches =  as.character(NBNICHES)


COLARCHS=c()#vector of HTML color codes of niches for plots TODO add sys.argv[9]
#system("source /scratch/anissa.el/miniconda3/etc/profile.d/conda.sh")
condaPath = "/scratch/anissa.el/miniconda3/envs/building-blocks"
pythonPath = "scratch/anissa.el/miniconda3/bin/python3"

system(paste("source /scratch/anissa.el/miniconda3/bin/activate", condaPath))
system(paste(pythonPath,"./main_nipmap.py",CellTypes," ",imageID," ",nsites," ",Radius," ",nbNiches," ",METHOD," ",ROOT_DATA_PATH," ",ROOT_OUTPUT_PATH))



file1 = "./pca_sites.json" # pca object on sites elements
file2 = "./AA_sites.json" # archetype Analysis object based on sites cell abundance
file3 = "./ca_sites.json" # cell abundance of randomly generated sites
file4 = "./cells_niches.json" # sites centered on cells and niches weights

#######---- Open .json files ----#######
json_data <- fromJSON(file=file1)
json_data2 <- fromJSON(file=file2)
json_data3 <- fromJSON(file=file3)
json_data4 <- fromJSON(file=file4)

##### LOAD OUTPUT OBJECTS
## Cell abundance in sites
sitesCellAb <- as_tibble(lapply(json_data3$cellAbSites,unlist))
write_csv(sitesCellAb%>%dplyr::select(-c(index, patient_id,site_id)),"sitesCA.csv")
## Archetypes coordinates in reduced PC space
Archs_3D <- do.call(cbind,lapply(json_data2$archs_coord,unlist))
## Projection of sites cell abundance in reduced PC space
pca3D <- matrix(unlist(json_data$PC_proj),nrow=17)[1:3,]
plotly::plot_ly(x=pca3D[1,],
                y=pca3D[2,],
                z=pca3D[3,],
                type="scatter3d",
                mode="marker")

## Niches weights(proportions) of all cells from all images 
cellsNiches <- as_tibble(lapply(json_data4$cells_niches,unlist))%>%mutate(site_id=as.numeric(site_id))





