#reticulate::use_condaenv(condaenv = "building-blocks",
#                         conda = "/scratch/anissa.el/miniconda3/bin/conda", required = TRUE)
######### NIPMAP SCRIPT
.libPaths("/scratch/anissa.el/R_old/x86_64-redhat-linux-gnu-library/4.0") ## add path to R libraries
require(devtools)
library(rjson)
library(tidyverse)
library(purrr)
library(reticulate) ## version 1.22 or  1.21 of the package !!
py_config()


### SET PYTHON ENVIRONMENT TO EITHER MINICONDA OR LOCAL ENVIRONMENT WHERE ALL PACKAGES ARE AVAILABLE
Sys.setenv(RETICULATE_PYTHON = "/scratch/anissa.el/miniconda3/envs/building-blocks/bin/python3")
#reticulate::use_python("/scratch/anissa.el/miniconda3/bin/python")
reticulate::use_python("/scratch/anissa.el/miniconda3/envs/building-blocks/bin/python3")
reticulate::use_condaenv(condaenv = "building-blocks",
                         conda = "/scratch/anissa.el/miniconda3/bin/conda",required=TRUE)

### SET WORKING DIRECTORY
dirName <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dirName)
source("./phenotypes_niches/functions_phenotypes_tmens.r")
SRCPYPATH = "./TMENS_analysis/src/"
py_run_string("import sys")
py_run_string("sys.path.insert(0, './TMENS_analysis/src')")
########-------- SOURCE PYTHON SCRIPTS FROM /TMENS_analysis/src -------########
#py_run_string("import sys")
#py_run_string("sys.path.insert(0, './TMENS_analysis/src')")
# eq <- reticulate::import_from_path("utils.equations",path=SRCPYPATH,convert=FALSE)#reticulate::source_python("/scratch/anissa.el/macro_micro_niches/macro_micro_niches2022/TMENS_analysis/src/utils/equations.py")
# ar <- reticulate::import_from_path("utils.archetypes",path=SRCPYPATH,convert=FALSE)#reticulate::source_python("/scratch/anissa.el/macro_micro_niches/macro_micro_niches2022/TMENS_analysis/src/utils/archetypes.py")
# reticulate::source_python("/scratch/anissa.el/macro_micro_niches/macro_micro_niches2022/TMENS_analysis/src/CellAbundance.py")
# ca <- reticulate::import_from_path("CellAbundance",path="./TMENS_analysis/src/",convert=FALSE)
# viz <- reticulate::import_from_path("utils.visualization",path=SRCPYPATH,convert=FALSE)
# #reticulate::source_python("/scratch/anissa.el/macro_micro_niches/macro_micro_niches2022/TMENS_analysis/src/utils/visualization.py")
# np <- reticulate::import("numpy", convert = FALSE)
# pd <- reticulate::import("pandas", convert = FALSE)
# sh <- reticulate::import("shutil",convert=FALSE)
# skl <- reticulate::import("sklearn.decomposition",convert=FALSE)

#reticulate::py_run_string("from src.utils.visualization import plot_cells_positions,colors")

############# CONTROL PANEL - SET PARAMETERS #############
ImageIDs = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
             20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32, 33, 34, 35, 36, 37,
             38, 39, 40, 41)
CELLTYPES = c('CD8-T', 'Other immune', 'DC / Mono', 'CD3-T', 'B', 'NK', 'Keratin-positive tumor', 'Tumor', 
               'CD4-T', 'Mesenchymal-like', 'Macrophages', 'Endothelial', 'Tregs', 'Unidentified', 'DC', 'Mono / Neu', 
               'Neutrophils')

METHOD ="gaussian" #only valid method to use for counting cells (gaussian kernel density)
path_figs = "/figs_niches"
path_toFigs <- paste0(dirName,path_figs)
RADIUS = 25
NSITES = 100
NBNICHES = 4
XSIZE = 800
YSIZE = 800
ROOT_DATA_PATH = "./TMENS_analysis/data/cell_positions_data"
ROOT_OUTPUT_PATH = "."#"./TMENS_analysis/output"
COLARCHS=c(rgb(255/255, 0/255, 223/255),rgb(255/255,0/255,0/255),rgb(70/255,203/255,236/255),rgb(0,0,0))# RGB code of colors for each archetype


##############------ISS data ------################
CELLTYPES = c('secretory', 'epithelial', 'proliferating smooth muscle', 'lymphatic endothelial', 'arterial', 
              'mesenchymal', 'ciliated epithelial', 'neuronal', 'pericyte', 'ASM', 
              'fibroblast', 'megakaryocyte', 'proximal progenitor', 'endothelial', 'neuroendocrine', 
              'chondroblast', 'Dendritic', 'macrophage', 'myeloid', 
              'mesothelial', 'monocyte', 'capillary', 'ILC3', 'B-cell', 
              'erythrocyte', 'myeloid progenitor', 'mast/basophil', 'venous', 
              'neutrophil', 'ILC2', 'lymphoid progenitor', 'natural killer')#c('secretory', 'epithelial', 'proliferating smooth muscle', 'lymphatic endothelial', 'arterial', 
#               'mesenchymal', 'ciliated epithelial', 'neuronal', 'pericyte', 'ASM', 
#               'fibroblast', 'megakaryocyte', 'proximal progenitor', 'endothelial', 'neuroendocrine', 
#               'chondroblast', 'airway smooth muscle', 'Dendritic', 'macrophage', 'myeloid', 
#               'mesothelial', 'monocyte', 'capillary', 'ILC3', 'B-cell', 
#               'erythrocyte', 'myeloid progenitor', 'dendritic', 'mast/basophil', 'venous', 
#               'neutrophil', 'ILC2', 'lymphoid progenitor', 'natural killer')

ImageIDs = c(1)
METHOD = "gaussian"
RADIUS = 25
NSITES = 20000 
NBNICHES = 5
XSIZE = 6500
YSIZE = 6500

COLARCHS=c(rgb(1, 0, 223/255),
           rgb(0, 0, 1),
           rgb(70/255, 203/255, 236/255),
           rgb(1, 0, 0),
           rgb(0, 0, 0),
           rgb(0, 153/255, 51/255),
           rgb(1, 1, 0))
ROOT_DATA_PATH= "./ISS_analysis/data"
ROOT_OUTPUT_PATH="./outputsISS/wrapper"
MARKERS <- read_rds("./ISS_analysis/data/ISS_marker_list.rds")
###############--------- RADIUS ANALYSIS ---------###############

# radiusAnalysis <-function(saveFig=TRUE, pathFig=paste0(path_toFigs,"./plot_radius.svg")){
#   nbs <- seq(log(5),log(190),length.out=10)#np.linspace(np.log(5), np.log(190), num=10)
#   radiusVec <- round(e^nbs)# range of radius values
#   #Generate sites cell abundance of each radius value + PCA
#   expl_var_ratio_gauss <- lapply(radiusVec,function(x){
#     gaussian_count_list = ca$generate_abundance_matrix(CELLTYPES, r_to_py(list(as.integer(ImageIDs))), as.integer(NSITES), 
#                                                     as.integer(x), method=METHOD, snr=as.integer(3), 
#                                                     root=ROOT_DATA_PATH,
#                                                     image_x_size=as.integer(XSIZE),
#                                                     image_x_size=as.integer(YSIZE))
#     sites <- ca$join_abundance_matrices(gaussian_count_list)[[1]]
#     pca = PCA()
#     pc = pca$fit_transform(sites)
#     cumsum(pca$explained_variance_ratio_)
#   })
#   names(expl_var_ratio_gauss) <- radiusVec
#   #plot variance explained by each PC for each radius value
#   plotRadius <- viz$radius_pc_all_variance(reticulate::dict(expl_var_ratio_gauss),radius_lim=as.integer(RADIUS),nPC_lim=as.integer(NBNICHES-1),cells_number = as.integer(length(CELLTYPES)+1),save_fig = saveFig, path_fig = pathFig)
#   return(plotRadius)
# }
# 
# radiusAnalysis(saveFig=TRUE, pathFig="./plot_radius.svg")


########------- GENERATE SITES WITH CELL ABUNDANCE FROM SPATIAL OMICS DATA -------######
#TODO function CellAbundanceSites
print('Generating sites with cell abundance...')
CellAb_list <- ca$generate_abundance_matrix(CELLTYPES, r_to_py(list(as.integer(ImageIDs))), as.integer(NSITES),
                                            as.integer(RADIUS),method=METHOD, snr=as.integer(3),image_x_size=as.integer(XSIZE),
                                         image_y_size=as.integer(YSIZE),center_sites_cells=FALSE,root=ROOT_DATA_PATH)

outCellAb <- ca$join_abundance_matrices(CellAb_list)
sites <- outCellAb[0]#outCellAb[[1]])
colnames(sites) <- CELLTYPES
patients_ids <- outCellAb[1] #outCellAb[[2]]
sites_ids <- outCellAb[2] #outCellAb[[3]]
CellAb_df <- data.frame()

## Get sites cell abundance data frame
CellAb_df <- lapply(CellAb_list,function(x){
  abundance_df <- data.frame(x['abundance_matrix']) #pd.DataFrame(ca.abundance_matrix,columns = CELLTYPES)
  colnames(abundance_df) <- CELLTYPES
  abundance_df[,"site_id"] <- as.vector(seq(1,nrow(abundance_df)))#abundance_df['site_id'] = np.arange(len(abundance_df))
  abundance_df["patient_id"] <- x["patient_id"]
  CellAb_df <- rbind(CellAb_df,abundance_df)
})
CellAb.df <- do.call(rbind,CellAb_df)

######Python script
py_run_string("import os")
py_run_string("import sys")
py_run_string("myDir = os.getcwd()")
#print(sys.path)
py_run_string("module_path = myDir + '/TMENS_analysis/'")
py_run_string("if module_path not in sys.path:
                  sys.path.append(module_path)")
#print(sys.path)
py_run_string("import pandas as pd")
py_run_string("import numpy as np")
py_run_string("from src.utils.archetypes import ArchetypalAnalysis")
py_run_string("from src.CellAbundance import CellAbundance, join_abundance_matrices, generate_abundance_matrix")
py_run_string("from src.utils.equations import compute_cells_niches_weights,get_niches_cell_abund")
py_run_string("from src.utils.visualization import plot_cells_positions")
py_run_string("from sklearn.decomposition import PCA")
py_run_string("import json")
py_run_string("import shutil")

if(length(ImageIDs)==1){
  ImageIDsPY <- r_to_py(list(as.integer(ImageIDs)))
}else{
  ImageIDsPY <- as.integer(ImageIDs)
}


nsites <- as.integer(NSITES)
radius <- as.integer(RADIUS)
xsize <- as.integer(XSIZE)
ysize <- as.integer(YSIZE)
py_run_string("CellAb_list=generate_abundance_matrix(r.CELLTYPES, r.ImageIDsPY,r.nsites ,
                                           r.radius ,method=r.METHOD, snr=3,image_x_size=r.xsize,
                                         image_y_size=r.ysize,center_sites_cells=False,root=r.ROOT_DATA_PATH)")

py_run_string("sites, patients_ids,sites_ids, _ = join_abundance_matrices(CellAb_list)")


##########----- PCA + ARCHETYPE ANALYSIS ON SITES CELL ABUND ----##########

# print('Dimension reduction PCA on sites cell abundance...')
# pca_obj = skl$PCA()   #PCA()
# pc_proj = pca_obj$fit_transform(sites)

# plotly::plot_ly(x = pc_proj[,1],
#                 y = pc_proj[,2],
#                 z = pc_proj[,3],
#                 type = "scatter3d",
#                 mode = "marker")

##### Python script of PCA
py_run_string("pca_obj=PCA()")
py_run_string("pc_proj=pca_obj.fit_transform(sites)")

######----- ARCHETYPE ANALYSIS -----######
print('Finding niches...')

# AA = ar$ArchetypalAnalysis(n_archetypes = as.integer(NBNICHES),
#                         tolerance = 0.001,max_iter = as.integer(200),random_state = as.integer(0),
#                         C = 0.0001,initialize = 'random',redundancy_try = as.integer(30))
# 
# pc_data<- np$array(pc_proj[,1:NBNICHES-1],dtype=np$float64)#np$array(unname(as.list(as.data.frame(pc_proj[,1:NBNICHES-1]))),dtype =np$float64)
# #py_run_string("AA.fit(pc_data)")
# NBNICHES=as.integer(NBNICHES)
# py_run_string("r.AA.fit_transform(r.pc_proj[:,:r.NBNICHES-1])",convert=FALSE)
# AA$fit_transform(r_to_py(pc_proj[,1:NBNICHES-1]))     #pc_data #as.matrix(pc_proj[,1:NBNICHES]) #np$array(as.numeric(pc_proj[,1:NBNICHES-1]),dtype=np$float64)
# print(paste0(as.character(AA$n_archetypes),' niches found !'))
# 
# # Display barplot of niches cell abundance
# NichesCellProfile(sitesCA = sites, Arch = AA, pcaObj = pca_obj)


##### Python script of Archetype ANalysis
nbniches <- as.integer(NBNICHES)
py_run_string("AA=ArchetypalAnalysis(n_archetypes = r.nbniches,
                        tolerance = 0.001,max_iter = 200,random_state = 0,
                        C = 0.0001,initialize = 'random',redundancy_try = 30) ")
py_run_string("AA.fit_transform(pc_proj[:,:r.nbniches-1])")

##########----- NICHE SEGMENTATION OF IMAGES ----##########
#For all images
print("Segmenting images into niches...")
#FIXME colors for niches segmentations, creat one unique variable
# COLARCHS2 = list(list(255, 0, 223),
# list(0, 0, 255),
# list(255, 0, 0),
# list(70, 203, 236),
# list(0, 0, 0))
# colors = c('#629563', '#044E75', '#CA8F04', '#645D0D','#43BC52', 
#           '#B25E89', '#2E3790', '#F118BE','#50974E', '#3273D6', 
#           '#0AF24B', '#A3F159','#933835', '#CEB134', '#226BCF', 
#           '#856218','#831CCB','#EEE8AA', '#FFA500', '#FFFF00', 
#           '#7CFC00', '#9ACD32', '#32CD32', '#20B2AA','#00FFFF',
#           '#191970', '#8A2BE2','#DDA0DD','#FF1493','#8B4513',
#           '#B0C4DE','#CD853F','#9400D3','#0000FF')
# for (i in ImageIDs){  #set here the Images IDs to plot
#   #print(i)
#   GRANULARITY = as.integer(5)
#   cell_data = pd$read_csv(paste0(ROOT_DATA_PATH,"/patient",as.character(i),"_cell_positions.csv"))
#   #print(cell_data)
#   #FIXME add legends to plot in python
#   fig <- viz$plot_cells_positions(cell_data, CELLTYPES, segment_image=TRUE, counting_type=METHOD,h=as.integer(1000),w=as.integer(1000),
#                               color_vector = r_to_py(COLARCHS2),segmentation_type='smooth', granularity=GRANULARITY, radius=as.integer(RADIUS),
#                              pca_obj=pca_obj, AA_obj=AA , to_plot = 'None',
#                               path_fig= paste0(path_toFigs,"/nichesSeg_patient",as.character(i),".svg"))
#   fig2 <- viz$plot_cells_positions(cell_data, CELLTYPES, segment_image=FALSE, counting_type=METHOD ,color_vector=colors,to_plot = 'all',
#                               path_fig= paste0(path_toFigs,"/cells_patient",as.character(i),".svg"))
# }
# 



##########----- GENERATE SITES CENTERED ON CELLS AND THEIR NICHE WEIGHTS ----##########
#TODO make functuon to compute niches weights of cells
print("Computing cells' niche weights, the operation might take some time...")
# CellAbCC_list = ca$generate_abundance_matrix(CELLTYPES, r_to_py(list(as.integer(ImageIDs))), 
#                                              as.integer(NSITES),as.integer(RADIUS),method=METHOD, 
#                                              snr = as.integer(3),image_x_size = as.integer(XSIZE),
#                                              image_y_size = as.integer(YSIZE),center_sites_cells = TRUE,root = ROOT_DATA_PATH)
# outCellAbCC <- ca$join_abundance_matrices(CellAbCC_list)#,center_sites_cells=TRUE
# outCellAbCC2 <- ca$join_abundance_matrices(CellAbCC_list,center_sites_cells=TRUE)
# sitesCC <- outCellAbCC[[1]]
# sitesCC <- outCellAbCC[0]
# colnames(sitesCC) <- CELLTYPES
# patients_ids2 <- outCellAbCC[[2]]
# patients_ids2<- outCellAbCC[1]
# sites_ids2 <- outCellAbCC[[3]]
# sites_ids2 <- outCellAbCC[2]
# 
# sites_ids2cc <- outCellAbCC2[[3]]


# CellAbCC = r_to_py(ca$generate_abundance_matrix(CELLTYPES, r_to_py(list(as.integer(ImageIDs))), 
#                                                 as.integer(NSITES),as.integer(RADIUS),method=METHOD,
#                                                 snr = as.integer(3),image_x_size = as.integer(XSIZE),
#                                                 image_y_size = as.integer(YSIZE),center_sites_cells = TRUE,
#                                                 root = ROOT_DATA_PATH))
# reticulate::py_run_string("sites_cc, patientsIDs, sitesIDs, _ = join_abundance_matrices(CellAbCC)")
# resCC <- r_to_py(ca$join_abundance_matrices(CellAbCC))
# 
# resNichesCC<-fromJSON(file="/scratch/ziqi.kang/NIPMAP_updated/macro_micro_niches/ISS_analysis/data/niche5_radius25_patientid[1]/cells_niches.json")
# resNiches <- fromJSON(file="./ISS_analysis/data/niche5_radius25_patientid[1]/cells_niches.json")
# cells_alfas <- do.call(rbind, resNichesCC[[1]])%>%t() %>% 
#   as.data.frame() %>% 
#   `colnames<-`(c("a1","a2","a3","a4","a5", "patient_id", "site_id",
#                  "cell_type_site", "TOT_cell_dens"))
# 
# cells_alfas2<- do.call(rbind, resNiches[[1]])%>%t() %>% 
#   as.data.frame() %>% 
#   `colnames<-`(c("a1","a2","a3","a4","a5", "patient_id", "site_id",
#                  "cell_type_site", "TOT_cell_dens"))

#sites_ids2cc <- outCellAbCC2[[3]]
# CellAbCC_list = generate_abundance_matrix(CELL_TYPES, patient_id, N_SITE, RADIUS, 
#                                           method=METHOD, snr=3, image_x_size = x_size, image_y_size = y_size, 
#                                           center_sites_cells=True, root=ROOT_DATA_PATH)
# sitesCC, patients_ids2, sites_ids2, _ = join_abundance_matrices(CellAbCC_list)
# NichesProf = get_niches_cell_abund(sitesCellAb=sites, pcaSites=pca_obj, ArchObj=AA, nComp=NBNICHES-1)
# sites_alfa = compute_cells_niches_weights(niches=NichesProf,cellsSites=sitesCC,nbNiches=NBNICHES)
# 
# sites_archs = pd.DataFrame(sites_alfa)
# sites_archs['patient_id'] = patients_ids2
# sites_archs["site_id"] = sites_ids2[:,0]
# sites_archs["cell_type_site"] = sites_ids2[:,1]
# sites_archs["TOT_cell_dens"]= sitesCC.sum(axis=1)
# 
# CellAbCC_df = pd.DataFrame(columns=[])
# for ca in CellAbCC_list:
#   df_ca = ca.get_site_cell_id_df()
# df_ca['patient_id'] = np.repeat(int(ca.patient_id), df_ca.shape[0], axis=0)
# # CellAbCC_df = pd.concat([CellAbCC_df, df_ca], axis=1)
# CellAbCC_df = CellAbCC_df.append(df_ca)
# CellAbCC_df = CellAbCC_df.reset_index(drop = True)

######## USELESS ==> to clean ??
# CellNeighb.df <- data.frame()
# CellAbCC_df<- lapply(CellAbCC_list,function(x){
#   df_ca <- x$get_site_cell_id_df()
#   df_ca[,'patient_id'] <- as.integer(x$patient_id)
#   CellNeighb.df <- rbind(CellNeighb.df,df_ca)
# })
# CellAbCC.df2 <- do.call(rbind,CellAbCC_df)

#Get cell abundance of niches 
# NichesCellProf <- eq$get_niches_cell_abund(sitesCellAb = sites,pcaSites = pca_obj,ArchObj = AA,nComp = as.integer(NBNICHES-1))
# # colnames(NichesCellProf) <- CELLTYPES
# # rownames(NichesCellProf) <- paste0("a",seq(1,NBNICHES))
# 
# # Compute niches weights of cells (excluding borders)
# sites_alfa <- eq$compute_cells_niches_weights(niches=NichesCellProf,cellsSites=sitesCC,nbNiches=as.integer(NBNICHES))
# sitesAlfa <- py_to_r(sites_alfa)
# colnames(sitesAlfa) <- paste0("a",seq(1,NBNICHES))
# sitesIDs2 <- py_to_r(pd$DataFrame(sites_ids2,columns=c("cell_id","cell_type")))
# sites_archs <- data.frame(sitesAlfa) #pd.DataFrame(sites_alfa)
# sites_archs[,'SampleID'] = py_to_r(patients_ids2)
# sites_archs[,"cell_id"] = as.numeric(sitesIDs2$cell_id)
# sites_archs[,"cell_type"] = sitesIDs2$cell_type

# colnames(sites_alfa) <- paste0("a",seq(1,NBNICHES))
# sites_archs <- data.frame(sites_alfa) #pd.DataFrame(sites_alfa)
# sites_archs[,'SampleID'] = patients_ids2
# sites_archs[,"cell_id"] = sites_ids2[,"site_id"]
# sites_archs[,"cell_type"] = sites_ids2[,"cell_type_site"]
#write_csv(sites_archs,"./cells_niches_weights_ISS.csv")
#sites_archs[,"TOT_cell_dens"]= rowSums(sitesCC) #sitesCC.sum(axis=1)


#### Python script for cells niche weight computing
py_run_string("print(type(r.nsites))")
py_run_string("print(type(r.CELLTYPES))")
py_run_string("print(type(r.ImageIDsPY))")
py_run_string("print(type(r.radius))")
py_run_string("print(type(r.ysize))")
py_run_string("print(type(r.METHOD))")
py_run_string("print(type(AA))")
py_run_string("print(type(r.ysize))")
py_run_string("CellAbCC = generate_abundance_matrix(r.CELLTYPES, r.ImageIDsPY, 
                                                 int(r.nsites),int(r.radius),method=r.METHOD,
                                                 snr = 3,image_x_size =int(r.xsize),
                                                 image_y_size =int(r.ysize),center_sites_cells =True,
                                                 root = r.ROOT_DATA_PATH)")
py_run_string("sites_cc, patientsIDcc, sitesIDcc, _ = join_abundance_matrices(CellAbCC)")
py_run_string("niches_cell_prof = get_niches_cell_abund(sitesCellAb = sites,pcaSites = pca_obj,ArchObj = AA,nComp = (r.nbniches-1))")
py_run_string("sites_alfa = compute_cells_niches_weights(niches=niches_cell_prof,cellsSites=sites_cc,nbNiches=int(r.nbniches))")
py_run_string("sites_archs = pd.DataFrame(sites_alfa)")
py_run_string("sites_archs['SampleID'] = patientsIDcc")
py_run_string("sites_archs['cell_id'] = sitesIDcc[:,0]")
py_run_string("sites_archs['cell_type'] = sitesIDcc[:,1]")
py_run_string("sites_archs.to_csv('./sites_cells_niches.csv',index=False)")


# from src.utils.equations import compute_cells_niches_weights,get_niches_cell_abund
# 
# n_comp = 3
# sites2, patients_ids2,sites_ids2,gradients2 = join_abundance_matrices(abundance_matrix_center_cells)
# 
# #tmens = np.dot(AA_3D.archetypes.T, pca.components_[:n_comp,:])+np.mean(sites, axis=0)
# 
# tmens = get_niches_cell_abund(sitesCellAb = sites,pcaSites = pca,ArchObj=AA_3D,nComp=n_comp)
# #print(tmens)
# NB_TMENS = tmens.shape[0]

#sites_alfas = compute_cells_niches_weights(nbNiches=NB_TMENS,niches=tmens,cellsSites =sites2)
###############--------SOLVE HERE--------- #################
# sites_archs<- as.data.frame(do.call(cbind,lapply(py$sites_archs,data.frame)))#as.data.frame(do.call(cbind,py$sites_archs)) #do.call(cbind, lapply(py$sites_archs,as.data.frame))#as_tibble(do.call(rbind,py$sites_archs)%>%t,rownames=NA)%>%mutate
# names(sites_archs) <- names(py$sites_archs)
niches<- paste0("a",as.vector(seq(1,NBNICHES,1)))
#py_run_file("main_nipmap.py")
# sites_archs<- sites_archs%>%rename_at(vars(matches("[0-9]")),~niches)%>%mutate(cell_id=as.numeric(cell_id))
sites_archs <- read_csv("./sites_cells_niches.csv")%>%rename_at(vars(matches("[0-9]")),~niches)
################---------NICHE-PHENOTYPE MAPPING ------------###################
NINTERFACES = 2

# Phenotypic markers
# ALL MARKERS
 MARKERS <- c("dsDNA","Vimentin","SMA","B7H3","FoxP3","Lag3","CD4","CD16","CD56",
              "OX40","PD1","CD31","PD-L1","EGFR","Ki67","CD209","CD11c","CD138",
              "CD163","CD68","CSF-1R","CD8","CD3","IDO","Keratin17","CD63","CD45RO",
              "CD20","p53","Beta catenin","HLA-DR","CD11b","CD45","H3K9ac","Pan-Keratin",
              "H3K27me3","H3K9ac_H3K27me3","phospho-S6","MPO","Keratin6","HLA_Class_1")

markers <- read.csv("./phenotypes_niches/data/proteins_by_frame.csv")%>%filter(Purpose=="Functional")%>%pull(Biomarker)

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
  left_join(sites_archs%>%filter(cell_type!="Unidentified"),.,by=c("SampleID","cell_id","cell_type"))%>%
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
  pivot_longer(cols=MARKERS,
               names_to="marker",values_to="value")
source("./phenotypes_niches/functions_phenotypes_tmens.r")

# scale2 <- function(x, na.rm = TRUE,cutOff=3){
#   res =(x - mean(x, na.rm = na.rm)) / sd(x, na.rm)
#   if(res>cutOff){
#     res = cutOff}
#   if (res< (-cutOff)){
#     res= -cutOff
#   } 
#   return(res)
# }
# 
cell_markers <- read.csv(paste(ROOT_DATA_PATH, "/cellData_",ImageIDs,".csv") %>% str_replace_all(" ", ""))

cell_markers$cell_type <- replace(cell_markers$cell_type, cell_markers$cell_type == "Erythroblast/RBC", "Erythroblast.RBC")

markers <- cell_markers %>%
  dplyr::select(-c(SampleID, cellLabelInImage, X, Y, cell_type, complex_celltype)) %>%
  ZscoreScale() %>%
  t() %>%
  as.data.frame()

cell_markers <- cell_markers %>%
  dplyr::select(c(SampleID, cellLabelInImage, X, Y, cell_type, complex_celltype)) %>%
  rename(cell_id=cellLabelInImage)%>%
  #mutate_at(MARKERS, scale2) # mutate_at(MARKERS, ZscoreScale)#
  cbind(markers)
#cell_markers[,MARKERS]
#


cellsPhen.niches <- left_join(sites_archs,cell_markers,by=c("SampleID","cell_id","cell_type"))%>%
  filter(!(is.na(a1)| is.na(a2)| is.na(a3) | is.na(a4) |is.na(a5)))%>%
  mutate(a1a2 = a1*a2)%>%
  mutate(a1a3 = a1*a3)%>%
  mutate(a1a4 = a1*a4)%>%
  mutate(a2a3 = a2*a3)%>%
  mutate(a2a4 = a2*a4)%>%
  mutate(a3a4 = a3*a4)%>%
  mutate(a1a5 = a1*a5)%>%
  mutate(a2a5 = a2*a5)%>%
  mutate(a3a5 = a3*a5)%>%
  mutate(a4a5=a4*a5)%>%
  #mutate(a1a2a3 = a1*a2*a3)%>%
  #mutate(a1a2a4 = a1*a2*a4)%>%
  #mutate(a1a3a4 = a1*a3*a4)%>%
  #mutate(a2a3a4 = a2*a3*a4)%>%
  #mutate(a1a2a3a4 = a1*a2*a3*a4)%>%
  pivot_longer(cols=all_of(MARKERS),
               names_to="marker",values_to="value")
# library(rjson)
# file4 = "./ISS_analysis/data/niche5_radius25_patientid[1]/cells_niches.json"# sites centered on cells and niches weights
# cells_niches <- fromJSON(file=file4)
# tmens_per_cell <- cells_niches$cells_niches
# 
# tmens_per_cell <- do.call(cbind, tmens_per_cell) %>% 
#   as.data.frame() %>% 
#   `colnames<-`(c("a1","a2","a3","a4","a5", "patient_id", "site_id",
#                  "cell_type_site", "TOT_cell_dens"))
# tmens_per_cell$patient_id <- as.integer(tmens_per_cell$patient_id)
# tmens_per_cell$site_id <- as.integer(tmens_per_cell$site_id)
# for (i in c("a1","a2","a3","a4","a5")) {
#   tmens_per_cell[,i] <- as.double(tmens_per_cell[,i])
# }
# tmens_per_cell$cell_type_site <- as.character(tmens_per_cell$cell_type_site)
# 
# cell_markers <- read.csv(paste(ROOT_DATA_PATH, "/cellData_",ImageIDs,".csv") %>% str_replace_all(" ", ""))
# cell_markers$cell_type <- replace(cell_markers$cell_type, cell_markers$cell_type == "Erythroblast/RBC", "Erythroblast.RBC")
# markers <- cell_markers %>% 
#   dplyr::select(-c(SampleID, cellLabelInImage, X, Y, cell_type, complex_celltype)) %>% 
#   ZscoreScale() %>% 
#   t() %>% 
#   as.data.frame()
# cell_markers <- cell_markers %>% 
#   dplyr::select(c(SampleID, cellLabelInImage, X, Y, cell_type, complex_celltype)) %>% 
#   rename(cell_id=cellLabelInImage)%>%
#   #mutate_at(MARKERS, scale2) # mutate_at(MARKERS, ZscoreScale)#
#   cbind(markers)
# # tmens_per_cell <- read_csv("/scratch/anissa.el/macro_micro_niches/macro_micro_niches2022/cells_niches_weights_ISS.csv")
# cellsPhen.niches <- left_join(tmens_per_cell %>% rename(SampleID=patient_id, cell_id=site_id, cell_type=cell_type_site),
#                               cell_markers,
#                               by=c("SampleID","cell_id","cell_type"))%>%
#   filter(!(is.na(a1)| is.na(a2)| is.na(a3) | is.na(a4) |is.na(a5)))%>%
#   mutate(a1a2 = a1*a2)%>%
#   mutate(a1a3 = a1*a3)%>%
#   mutate(a1a4 = a1*a4)%>%
#   mutate(a2a3 = a2*a3)%>%
#   mutate(a2a4 = a2*a4)%>%
#   mutate(a3a4 = a3*a4)%>%
#   mutate(a1a5 = a1*a5)%>%
#   mutate(a2a5 = a2*a5)%>%
#   mutate(a3a5 = a3*a5)%>%
#   mutate(a4a5 = a4*a5)%>%
#   #mutate(a1a2a3 = a1*a2*a3)%>%
#   #mutate(a1a2a4 = a1*a2*a4)%>%
#   #mutate(a1a3a4 = a1*a3*a4)%>%
#   #mutate(a2a3a4 = a2*a3*a4)%>%
#   #mutate(a1a2a3a4 = a1*a2*a3*a4)%>%
#   pivot_longer(cols=all_of(MARKERS),
#                names_to="marker",values_to="value")

source("./phenotypes_niches/functions_phenotypes_tmens.r")

CM <- correlation_niches_CM(markersCells.niches=cellsPhen.niches,Markers=MARKERS,corrMeth="spearman",coreIntf2,1/100,0.3)

CM <- correlation_niches_CM(markersCells.niches=cellsPhen.niches,Markers=markers,corrMeth="spearman",coreIntf2,1/100,0.3)

##Filtering step: remove markers whose association with cancer niche is probably due to spill over or cell segmentation error
# All Keratin6/Beta-catenin + phenotypes associated with cancer(and inteface with cancer) niche
Krt.filt <- names(which(CM[rownames(CM)[which(grepl("Keratin6",rownames(CM),fixed=TRUE))],"a3"]>0.3))
#Krt.filt <- names(which(CM[rownames(CM)[which(grepl("Keratin6",rownames(CM),fixed=TRUE))],"a3"]>0.3))

BCat.filt <- names(which(CM[rownames(CM)[which(grepl("Beta catenin",rownames(CM),fixed=TRUE))],"a3"]>0.2))

cMf <- CM[!rownames(CM)%in%BCat.filt& !rownames(CM)%in% Krt.filt,] # remove Keratin and Beta catenin + cell phenotypes associated with cancer niche
cMf <- CM
# Plot heatmaps of niches-phenotypes associations
plot_heatmap_CT(CM.mat=cMf,coreIntf2,"./outputsISS/wrapper/CMbyCells.pdf")

plot_heatmap_CT(CM.mat=cMf,coreIntf2,"./CMbyCells2.pdf")

plot_heatmap_markers(CM.mat=cMf,coreIntf2,paste0(path_toFigs,"/CMbyMarkers.pdf"))


#Named vector of archetypes: niches names
NichesNames <- c("a1"="TLS","a2"="inflammatory","a3"="cancer","a4"="fibrotic")

###########....... TABLE OF NICHE-ASSOCIATED CELL PHENOTYPES .......###########

#Archetype weights of randomly generated sites over TNBC patients images
archetypes_sites <-  AA$alfa
rownames(archetypes_sites) <- names(NichesNames)
#Cell abundance of sites
cellAbSites <- CellAb.df

nichesCellAb <- NichesCellProf%>%as_tibble(rownames=NA)%>%rownames_to_column(var="archetype")%>%
  mutate(archetype=str_replace_all(archetype,NichesNames))%>%
  column_to_rownames(var="archetype")%>%t()%>%
  as_tibble(rownames=NA)%>%rownames_to_column(var="cell_type")

#Sort cell types per niche by decreasing cell abundance
nichesCA.sort <- nichesCellAb%>%
  pivot_longer(cols = as.vector(NichesNames),names_to="niche",values_to="cell_density")%>%
  group_by(niche)%>%arrange(desc(cell_density))


archsSitesCellAb <- cbind(cellAbSites,archetypes_sites%>%t)%>% 
  pivot_longer(cols=append(CELLTYPES,"Unidentified"),names_to="cell_type",values_to= "cell_density")%>%ungroup()

### Get cell types enriched in each niche (take top 1% sites closest to niche)
archs.CT <- get_CT_enriched_all_archs(archsSitesCellAb,NichesNames)%>%#(archetypes_sites%>%t%>%as_tibble(rownames=NA),cellAbSites,thresh = 0.99)%>%
  group_by(niche)%>%
  filter(!(cell_type %in% c("DC / Mono","CD3-T", "Mono / Neu", "Other immune","Unidentified")))%>%mutate(cell_type = paste(unique(cell_type),collapse="\n"))%>%distinct(niche,cell_type)


## Get table of niches/interfaces-associated cell phenotypes
TableNichesPhenotypes(cMf,archs.CT,NichesNames,path_toFigs) #tabNichePhen.pdf saved


