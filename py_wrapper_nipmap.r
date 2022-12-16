 
#install.packages("reticulate")
# remotes::install_github("reticulate")
.libPaths("/scratch/anissa.el/R_old/x86_64-redhat-linux-gnu-library/4.0")
require(devtools)
#install_version("reticulate", version = "1.22", repos = "http://cran.us.r-project.org")
#install.packages("rjson")
library(rjson)
library(tidyverse)
library(purrr)
#install.packages("reticulate", dependencies = TRUE, INSTALL_opts = '--no-lock')
library(reticulate)
py_config()
#conda_# update.packages(instlib = "local")
# reticulate::py_config()
#reticulate::py_install("qpsolvers")
#conda_list()

Sys.setenv(RETICULATE_PYTHON = "/scratch/anissa.el/miniconda3/envs/building-blocks/bin/python3")
 #reticulate::use_python("/scratch/anissa.el/miniconda3/bin/python")
reticulate::use_python("/scratch/anissa.el/miniconda3/envs/building-blocks/bin/python3")
# reticulate::py_config()
reticulate::use_condaenv(condaenv="building-blocks",
                         conda="/scratch/anissa.el/miniconda3/bin/conda",required=TRUE)
#
#reticulate::use_condaenv(condaenv="/scratch/anissa.el/miniconda3/envs/building-blocks",
#                         conda = "/scratch/anissa.el/miniconda3/bin/conda",required=TRUE)

### SET WORKING DIRECTORY
dirName <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dirName)
#Master\\\ File

#TODO add input: named vectors of colors codes(html) for visualizing niches
#TODO add source code: functions_phenotypes_tmens.r ==> functions_nipmap.r 

#TODO plot correlations heatmaps organized by cell types & by markers
#TODO plot table of cell phenotypes for each niche/interface

### CONTROL PANEL: SET PARAMETERS AND CONSTANTS
# CELLTYPES=paste(c('CD8-T', 'Other\\\ immune', 'DC\\\ /\\\ Mono', 'CD3-T', 'B', 'NK', 'Keratin-positive\\\ tumor', 'Tumor', 
#             'CD4-T', 'Mesenchymal-like', 'Macrophages', 'Endothelial', 'Tregs', 'Unidentified', 'DC', 'Mono\\\ /\\\ Neu', 
#             'Neutrophils'),collapse=",")
# ImageIDs=paste(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
#            20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32, 33, 34, 35, 36, 37,
#            38, 39, 40, 41),collapse = ",")
# NSITES=as.character(100) # number of sites generated per image
# RADIUS=as.character(25) # radius of site in micrometer²
# NBNICHES = as.character(4) # number of niches to find (in PC space of NBNICHES-1 dimensions)
# METHOD ="gaussian"
# ROOT_DATA_PATH="./TMENS_analysis/data/cell_positions_data" 
# ROOT_OUTPUT_PATH="./TMENS_analysis/output"
# COLARCHS=c()#vector of HTML color codes of niches for plots TODO add sys.argv[9] 
# 
# #system("source /scratch/anissa.el/miniconda3/etc/profile.d/conda.sh")
# system("source /scratch/anissa.el/miniconda3/bin/activate /scratch/anissa.el/miniconda3/envs/building-blocks")
# system(paste("/scratch/anissa.el/miniconda3/envs/building-blocks/bin/python3 ./main_nipmap.py",CELLTYPES," ",ImageIDs," ",NSITES," ",RADIUS," ",NBNICHES," ",METHOD," ",ROOT_DATA_PATH," ",ROOT_OUTPUT_PATH))
# 
# 
# 
# file1 = "./pca_sites.json" # pca object on sites elements 
# file2 = "./AA_sites.json" # archetype Analysis object based on sites cell abundance
# file3 = "./ca_sites.json" # cell abundance of randomly generated sites
# file4 = "./cells_niches.json" # sites centered on cells and niches weights
# 
# #######---- Open .json files ----#######
# json_data <- fromJSON(file=file1)
# json_data2 <- fromJSON(file=file2)
# json_data3 <- fromJSON(file=file3)
# json_data4 <- fromJSON(file=file4)
# 
# ##### LOAD OUTPUT OBJECTS
# ## Cell abundance in sites
# sitesCellAb <- as_tibble(lapply(json_data3$cellAbSites,unlist))
# write_csv(sitesCellAb%>%dplyr::select(-c(index, patient_id,site_id)),"sitesCA.csv")
# ## Archetypes coordinates in reduced PC space 
# Archs_3D <- do.call(cbind,lapply(json_data2$archs_coord,unlist))
# ## Projection of sites cell abundance in reduced PC space
# pca3D <- matrix(unlist(json_data$PC_proj),nrow=17)[1:3,]
# plotly::plot_ly(x=pca3D[1,],
#                 y=pca3D[2,],
#                 z=pca3D[3,],
#                 type="scatter3d",
#                 mode="marker")


########--------EXPORT NICHE SEGMENTED IMAGES IN SVG FIGURES--------########
reticulate::source_python("/scratch/anissa.el/macro_micro_niches/macro_micro_niches2022/TMENS_analysis/src/utils/visualization.py")
reticulate::source_python("/scratch/anissa.el/macro_micro_niches/macro_micro_niches2022/TMENS_analysis/src/CellAbundance.py")
reticulate::source_python("/scratch/anissa.el/macro_micro_niches/macro_micro_niches2022/TMENS_analysis/src/utils/archetypes.py")

reticulate::py_run_string("import os")
reticulate::py_run_string("import pandas as pd")
reticulate::py_run_string("import numpy as np")
reticulate::py_run_string("from shutil import make_archive ")
reticulate::py_run_string("from pandas import read_csv")
reticulate::py_run_string("myDir = os.getcwd()")
reticulate::py_run_string("module_path = myDir + '/TMENS_analysis/'")
reticulate::py_run_string("if module_path not in sys.path:
                              sys.path.append(module_path)")
#reticulate::py_run_string("from src.utils.visualization import plot_cells_positions")
reticulate::py_run_string("from src.utils.archetypes import ArchetypalAnalysis")
reticulate::py_run_string("from src.CellAbundance import CellAbundance, join_abundance_matrices, generate_abundance_matrix")
reticulate::py_run_string("from src.utils.equations import compute_cells_niches_weights,get_niches_cell_abund")
reticulate::py_run_string("from src.utils.visualization import plot_cells_positions")
reticulate::py_run_string("from sklearn.decomposition import PCA")

ImageIDs = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
           20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32, 33, 34, 35, 36, 37,
           38, 39, 40, 41)
CELLTYPES = c('CD8-T', 'Other immune', 'DC / Mono', 'CD3-T', 'B', 'NK', 'Keratin-positive tumor', 'Tumor', 
              'CD4-T', 'Mesenchymal-like', 'Macrophages', 'Endothelial', 'Tregs', 'Unidentified', 'DC', 'Mono / Neu', 
              'Neutrophils')#c('CD8-T', 'Other\\\ immune', 'DC\\\ /\\\ Mono', 'CD3-T', 'B', 'NK', 'Keratin-positive\\\ tumor', 'Tumor', 
              #'CD4-T', 'Mesenchymal-like', 'Macrophages', 'Endothelial', 'Tregs', 'Unidentified', 'DC', 'Mono\\\ /\\\ Neu', 
              #'Neutrophils'
              #'
METHOD ="gaussian"
path_figs = "/figs_niches"
#reticulate::py_run_string("path_toFigs = os.path.join(myDir,path_figs)")
path_toFigs <- paste0(myDir,path_figs)
RADIUS = 25
NSITES = 100
NBNICHES <- as.integer(4)
ROOT_DATA_PATH="./TMENS_analysis/data/cell_positions_data" 
ROOT_OUTPUT_PATH="./TMENS_analysis/output"
# reticulate::py_run_string("CellAb_list = generate_abundance_matrix(CELLTYPES ,ImageIDs,NSITES,RADIUS,method=METHOD, snr=3,center_sites_cells=False,root=ROOT_DATA_PATH)
# sites, patients_ids,sites_ids, _ = join_abundance_matrices(CellAb_list)
# CellAb_df = pd.DataFrame()
# print('Generating sites with cell abundance...')
# for ca in CellAb_list:
#   abundance_df = pd.DataFrame(ca.abundance_matrix,columns = CELLTYPES)
#   abundance_df['site_id'] = np.arange(len(abundance_df))
#   abundance_df['patient_id'] = ca.patient_id
#   CellAb_df = CellAb_df.append(abundance_df)
# CellAb_df = CellAb_df.reset_index()


# #####----- PCA ON SITES ABUNDANCE ----#####
# print('Dimension reduction of sites cell abundances...')
# pca_obj = PCA()
# pc_proj = pca_obj.fit_transform(sites)
# #####----- ARCHETYPE ANALYSIS ----#####
# print('Finding niches...')
# AA = ArchetypalAnalysis(n_archetypes = NBNICHES, 
#                      tolerance = 0.001, 
#                      max_iter = 200, 
#                      random_state = 0, 
#                      C = 0.0001, 
#                      initialize = 'random',
#                      redundancy_try = 30)
# AA.fit_transform(pc_proj[:,:NBNICHES-1])
# print(str(NBNICHES)+' found !')
py_run_string("path_figs = 'figs_niches'")
#reticulate::py_run_string("path_toFigs = os.path.join(myDir,path_figs)")
py_run_string("path_toFigs = os.path.join(myDir,path_figs)")
py_run_string("if os.path.isdir(path_toFigs)==False:
   os.mkdir(path_toFigs)")

########------- GENERATE SITES WITH CELL ABUNDANCE FROM SPATIAL OMICS DATA -------######
#CellAb_list = generate_abundance_matrix(CELLTYPES, ImageIDs, NSITES,RADIUS,method=METHOD, snr=3,center_sites_cells=False,root=ROOT_DATA_PATH)
#sites, patients_ids,sites_ids, _ = join_abundance_matrices(CellAb_list)

print('Generating sites with cell abundance...')
CellAb_list <- generate_abundance_matrix(CELLTYPES, as.integer(ImageIDs), as.integer(NSITES),as.integer(RADIUS),method=METHOD, snr=as.integer(3),center_sites_cells=FALSE,root=ROOT_DATA_PATH)
outCellAb <- join_abundance_matrices(CellAb_list)
sites <- outCellAb[[1]]
colnames(sites)<-CELLTYPES
patients_ids <- outCellAb[[2]]
sites_ids <- outCellAb[[3]]
CellAb_df <- data.frame()#pd.DataFrame()


CellAb_df <-lapply(CellAb_list,function(x){
  abundance_df = data.frame(x['abundance_matrix']) #pd.DataFrame(ca.abundance_matrix,columns = CELLTYPES)
  colnames(abundance_df)= CELLTYPES
  abundance_df[,"site_id"] <- as.vector(seq(1,nrow(abundance_df)))#abundance_df['site_id'] = np.arange(len(abundance_df))
  abundance_df["patient_id"] <- x["patient_id"] #abundance_df['patient_id'] = ca.patient_id
  CellAb_df<- rbind(CellAb_df,abundance_df) #= CellAb_df.append(abundance_df)
})
CellAb.df <- do.call(rbind,CellAb_df)
#CellAb_df = CellAb_df.reset_index()

##########----- PCA + ARCHETYPE ANALYSIS ON SITES CELL ABUND ----##########

print('Dimension reduction PCA on sites cell abundance...')
pca_obj = PCA()
pc_proj = pca_obj$fit_transform(sites)

#pca3D <- matrix(pc_proj ,nrow=length(CELLTYPES))[1:NBNICHES-1,]
plotly::plot_ly(x=pc_proj[,1],#pca3D[1,],
                y=pc_proj[,2],#pca3D[2,],
                z=pc_proj[,3],#pca3D[3,],
                type="scatter3d",
                mode="marker")

######----- ARCHETYPE ANALYSIS ----######
print('Finding niches...')
np <- reticulate::import("numpy", convert = FALSE)
pd <- reticulate::import("pandas", convert = FALSE)
sh <- reticulate::import("shutil",convert=FALSE)
NBNICHES <- as.integer(NBNICHES)
#reticulate::py_run_string("AA = ArchetypalAnalysis(n_archetypes = NBNICHES,tolerance = 0.001,max_iter = 200,random_state = 0,C = 0.0001,initialize = 'random',redundancy_try = 30")
AA = ArchetypalAnalysis(n_archetypes = NBNICHES,
                        tolerance = 0.001,max_iter = as.integer(200),random_state = as.integer(0),
                        C = 0.0001,initialize = 'random',redundancy_try = as.integer(30))


pc_data<- np$array(pc_proj[,1:NBNICHES-1],dtype=np$float64)#np$array(unname(as.list(as.data.frame(pc_proj[,1:NBNICHES-1]))),dtype =np$float64)
#py_run_string("AA.fit(pc_data)")
####### Là où le bât blesse,... Il ne blesse plus niahahahaha

AA$fit_transform(pc_data)#as.matrix(pc_proj[,1:NBNICHES]) #np$array(as.numeric(pc_proj[,1:NBNICHES-1]),dtype=np$float64)
print(paste0(as.character(AA$n_archetypes),' niches found !'))
 # if os.path.isdir(path_toFigs)==False:
 #   os.mkdir(path_toFigs)

###----PLOT CELL ENRICHMENT
## Bar plot of cells proportions for each archetype
NichesCellProf <- get_niches_cell_abund(sitesCellAb=sites,pcaSites=pca_obj,ArchObj=AA,nComp=as.integer(NBNICHES-1))
colnames(NichesCellProf) <- CELLTYPES
rownames(NichesCellProf) <- paste0("arch",seq(1,NBNICHES))
NichesCellProp <- NichesCellProf%>%as_tibble(rownames=NA)%>%
  rownames_to_column(var="archetype")%>%
  pivot_longer(cols=all_of(CELLTYPES),names_to="cell_type",values_to = "cell_density")

##TODO marke fucntion barplot of cellular profile of archetypes prior to defining niches
barplot1 <- ggplot(data = NichesCellProp, aes(x = cell_type, y = cell_density,fill = archetype)) +
  geom_bar(stat = "identity",position = position_dodge(),width = 0.6) + 
  #theme(axis.text.x = element_text(angle = 35,hjust = 1.2,vjust = 0.9))+
  theme(axis.text.x = element_text(angle = 90, vjust = .2))+
  xlab ("") + ylab("cell density")
barplot1
ggsave("./testBarplot1.pdf",barplot1,height=3,width=4)


##########----- NICHE SEGMENTATION OF IMAGES ----##########
#For all images
print("Segmenting images into niches...")
COLARCHS = c(c(255, 0, 223),c(255,0,0),c(70,203,236),c(0,0,0))# RGB code
for (i in c(1,2)){ #ImageIDs
  print(i)
  GRANULARITY = as.integer(5)
  cell_data = pd$read_csv(paste0(ROOT_DATA_PATH,"/patient",as.character(i),"_cell_positions.csv"))
  print(cell_data)
  #FIXME add legends to plot in python
  
  fig <- plot_cells_positions(cell_data, CELLTYPES, segment_image=TRUE, counting_type=METHOD,
                          color_vector=COLARCHS,segmentation_type='colors', granularity=GRANULARITY, radius=as.integer(RADIUS),
                          pca_obj=pca_obj, AA_obj=AA , to_plot = 'None',
                          path_fig= paste0(path_toFigs,"/nichesSeg_patient",as.character(i),".svg"))
}
make_archive("figs_niches","zip", path_toFigs)  



##########----- GENERATE SITES CENTERED ON CELLS AND THEIR NICHE WEIGHTS ----##########
print("Computing cells' niche weights, the operation might take some time...")
CellAbCC_list = generate_abundance_matrix(CELLTYPES, as.integer(ImageIDs), as.integer(NSITES),as.integer(RADIUS),method=METHOD, snr=3,center_sites_cells=TRUE,root=ROOT_DATA_PATH)
outCellAbCC <- join_abundance_matrices(CellAbCC_list,center_sites_cells=TRUE)

sitesCC <- outCellAbCC[[1]]
colnames(sitesCC) <- CELLTYPES
patients_ids2 <- outCellAbCC[[2]]
sites_ids2 <- outCellAbCC[[3]]
#sitesCC, patients_ids2,sites_ids2, _ =
#CellAbCC_df = pd.DataFrame()
CellNeighb.df <- data.frame()
CellAbCC_df<- lapply(CellAbCC_list,function(x){
  df_ca <- x$get_site_cell_id_df()
  df_ca[,'patient_id'] <- as.integer(x$patient_id)
  CellNeighb.df <- rbind(CellNeighb.df,df_ca)
})
CellAbCC.df2 <- do.call(rbind,CellNeighb.df)

# for ca in CellAbCC_list:
#   df_ca = ca.get_site_cell_id_df()
# df_ca['patient_id'] = int(ca.patient_id)
# CellAbCC_df = CellAbCC_df.append(df_ca)
# CellAbCC_df = CellAbCC_df.reset_index(drop = True)


sites_alfa <- compute_cells_niches_weights(niches=NichesCellProf,cellsSites=sitesCC,nbNiches=as.integer(NBNICHES))

sites_archs <- data.frame(sites_alfa) #pd.DataFrame(sites_alfa)
sites_archs[,'SampleID'] = patients_ids2
sites_archs[,"site_id"] = sites_ids2[,"site_id"]
sites_archs[,"cell_type_site"] = sites_ids2[,"cell_type_site"]
sites_archs[,"TOT_cell_dens"]= rowSums(sitesCC) #sitesCC.sum(axis=1)

## Niches weights(proportions) of all cells from all images 

#cellsNiches <- as_tibble(lapply(json_data4$cells_niches,unlist))%>%mutate(site_id=as.numeric(site_id))
#TODO add constant: functional markers
MARKERS = c() #phenotypic markers

#TODO add function to select niches / interfaces to map to cell phenotypes
CM <- correlation_niches_CM(markersCells.niches=MarkersCellsTMENS2,Markers=FuncMarkers,corrMeth="spearman",coreIntf=coreIntf2)


#TODO add  input: named vectors for naming the archetypes





####################--------------GARBAGE--------------####################

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

##########----- SAVE OUTPUTS IN CSV ----##########


