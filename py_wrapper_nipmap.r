# install.packages("reticulate")
# remotes::install_github("reticulate")
.libPaths("/scratch/anissa.el/R_old/x86_64-redhat-linux-gnu-library/4.0")
install.packages("rjson")
library(rjson)
library(tidyverse)
library(purrr)
#install.packages("reticulate", dependencies = TRUE, INSTALL_opts = '--no-lock')
# library(reticulate)
# update.packages(instlib = "local")
# reticulate::py_config()
#reticulate::py_install("qpsolvers")
#conda_list()

# Sys.setenv(RETICULATE_PYTHON = "/scratch/anissa.el/miniconda3/envs/building-blocks/bin/python3")
# reticulate::use_python("/scratch/anissa.el/miniconda3/bin/python")
# reticulate::py_config()
# reticulate::use_condaenv(condaenv="building-blocks",
#                          conda = "/scratch/anissa.el/miniconda3/bin/conda",required=TRUE)
#reticulate::use_condaenv(condaenv="/scratch/anissa.el/miniconda3/envs/building-blocks",
#                         conda = "/scratch/anissa.el/miniconda3/bin/conda",required=TRUE)

### SET WORKING DIRECTORY
dirName <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dirName)
#Master\\\ File
### CONTROL PANEL: SET PARAMETERS AND CONSTANTS
CELLTYPES=paste(c('CD8-T', 'Other\\\ immune', 'DC\\\ /\\\ Mono', 'CD3-T', 'B', 'NK', 'Keratin-positive\\\ tumor', 'Tumor', 
            'CD4-T', 'Mesenchymal-like', 'Macrophages', 'Endothelial', 'Tregs', 'Unidentified', 'DC', 'Mono\\\ /\\\ Neu', 
            'Neutrophils'),collapse=",")
ImageIDs=paste(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
           20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32, 33, 34, 35, 36, 37,
           38, 39, 40, 41),collapse = ",")
NSITES=as.character(100)
RADIUS=as.character(25)
NBNICHES = as.character(4)
METHOD ="gaussian"
ROOT_DATA_PATH="./TMENS_analysis/data/cell_positions_data" 
ROOT_OUTPUT_PATH="./TMENS_analysis/output"
#system("source /scratch/anissa.el/miniconda3/etc/profile.d/conda.sh")
system("source /scratch/anissa.el/miniconda3/bin/activate /scratch/anissa.el/miniconda3/envs/building-blocks")
system(paste("/scratch/anissa.el/miniconda3/envs/building-blocks/bin/python3 ./main_nipmap.py",CELLTYPES," ",ImageIDs," ",NSITES," ",RADIUS," ",NBNICHES," ",METHOD," ",ROOT_DATA_PATH," ",ROOT_OUTPUT_PATH))

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
cellsNiches <- as_tibble(lapply(json_data4$cells_niches,unlist))%>%
  mutate(site_id=as.numeric(site_id))


### IMPORT PYTHON LIBRARIES

# conda_install("building-blocks", "qpsolvers")
# reticulate::py_install("qpsolvers")
# #### SOURCE PYTHON SCRIPTS 
# reticulate::source_python("./TMENS_analysis/src/utils/archetypes.py")
# reticulate::source_python("./TMENS_analysis/src/CellAbundance.py")
# reticulate::source_python("./TMENS_analysis/src/utils/equations.py")
# #qpSolve <- import("qpsolvers")
# reticulate::py_run_string("from qsolvers import solve_qp")
# reticulate::py_module_available("qpsolvers")
# reticulate::py_module_available("sklearn.decomposition")
# 
# ######--- GENERATE SITES WITH CELL ABUNDANCE FROM SPATIAL OMICS DATA ---######
# CellAb_list = generate_abundance_matrix(CELLTYPES, as.integer(ImageIDs), as.integer(NSITES),
#                                         as.integer(RADIUS),
#                                         method=METHOD, snr=3,
#                                         center_sites_cells=FALSE,
#                                         root=ROOT_DATA_PATH)
# 
# #list(sites, patients_ids,sites_ids) = join_abundance_matrices(CellAb_list )
# JointAbundance.mat = join_abundance_matrices(CellAb_list )
# #rm(resultJointAbundance.df)
# sites <- JointAbundance.mat[[1]]
# colnames(sites)<- CELLTYPES
# ##########----- PCA + ARCHETYPE ANALYSIS ON SITES CELL ABUND ----##########
# py_run_string("pca_3d = PCA()")
# py_run_string("pc3d = pca_3d.fit_transform(sites)")
# AA <- ArchetypalAnalysis(n_archetypes = as.integer(NBNICHES), 
#                    tolerance = 0.001, 
#                    max_iter = 200, 
#                    random_state = 0, 
#                    C = 0.0001, 
#                    initialize = "random",
#                    redundancy_try = 30)
# AA.fit_transform()

##########----- GENERATE SITES CENTERED ON CELLS AND THEIR NICHE WEIGHTS ----##########

##########----- SAVE OUTPUTS IN CSV ----##########


