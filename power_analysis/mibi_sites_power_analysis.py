import os
import sys
myDir = os.getcwd()
#print(sys.path)
module_path = "../TMENS_analysis/"#myDir + "/TMENS_analysis/"
if module_path not in sys.path:
    sys.path.append(module_path)
#print(sys.path)
import pandas as pd
import numpy as np
from src.utils.archetypes import ArchetypalAnalysis
from src.CellAbundance import CellAbundance, join_abundance_matrices, generate_abundance_matrix
from src.utils.equations import compute_cells_niches_weights,get_niches_cell_abund
from src.utils.visualization import plot_cells_positions,radius_pc_all_variance,plot_cells_positions,archetypes_bar_plot
from sklearn.decomposition import PCA
from sklearn.metrics import mean_squared_error
import json
import shutil
import math
from shutil import make_archive
from multiprocessing import Pool
from functools import partial
import functools
import json
import re
from itertools import permutations

#### PARAMETERS
CELLTYPES = ['CD8-T', 'Other immune', 'DC / Mono', 'CD3-T', 'B', 'NK', 'Keratin-positive tumor', 'Tumor','CD4-T', 'Mesenchymal-like', 'Macrophages', 'Endothelial', 'Tregs', 'Unidentified', 'DC', 'Mono / Neu','Neutrophils']
ImageIDs = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32, 33, 34, 35, 36, 37,38, 39, 40, 41]
#NSITES = 100
RADIUS = 25
NBNICHES = 4
METHOD ="gaussian"
XSIZE= 800
YSIZE = 800
ROOT_DATA_PATH = "./TMENS_analysis/data/cell_positions_data"
ROOT_OUTPUT_PATH = "./TMENS_analysis/output"
COLARCHS = np.array([[255, 0, 223],[255,0,0],[70,203,236],[0,0,0]])
NBITER = 100


def get_niches(sites):
  # NSITES = int((SamplingInt/100)*Atot/As)
  # #print(NSITES)
  # ### Generate sites cell abundance
  # CellAb_list = generate_abundance_matrix(CELLTYPES, ImageIDs, NSITES,RADIUS,method=METHOD,random_seed=j,snr=3,center_sites_cells=False,root=ROOT_DATA_PATH)#generate_abundance_matrix(CELLTYPES, ImageIDs, NSITES,RADIUS,method=METHOD,random_seed=1022,snr=3,center_sites_cells=False,root=ROOT_DATA_PATH)
  # sites, patients_ids,sites_ids, _ = join_abundance_matrices(CellAb_list)
  # CellAb_df = pd.DataFrame()
  # print("Generating sites with cell abundance...")
  # for ca in CellAb_list:
  #     abundance_df = pd.DataFrame(ca.abundance_matrix,columns = CELLTYPES)
  #     abundance_df['site_id'] = np.arange(len(abundance_df))
  #     abundance_df['patient_id'] = ca.patient_id
  #     CellAb_df = CellAb_df.append(abundance_df)
  # CellAb_df = CellAb_df.reset_index()
  #print(sites)
  #####----- PCA ON SITES ABUNDANCE ----#####
  print("Dimension reduction of sites cell abundances...")
  pca_obj = PCA()
  pc_proj = pca_obj.fit_transform(sites)
  #####----- ARCHETYPE ANALYSIS ----#####
  print("Finding niches...")
  AA = ArchetypalAnalysis(n_archetypes = NBNICHES,
                      tolerance = 0.0001, #0.001,
                      max_iter = 500,#10000,#2000,#200,
                      random_state = 0,
                      C = 0.0001,
                      initialize = "furthest_sum",#'random',
                      redundancy_try = 30)
  AA.fit_transform(pc_proj[:,:NBNICHES-1])
  print(str(NBNICHES)+" niches found!")
  
  niches_cell_profile = get_niches_cell_abund(sites.to_numpy(),pca_obj,AA,NBNICHES-1)
  
  #print(niches_cell_profile)
  niches_cellAb = pd.DataFrame(niches_cell_profile,columns=CELLTYPES)
  #print(niches_cellAb.head)
#archetype_colors = [[255/255, 0., 223/255],[255/255,0.,0.],[70/255,203/255,236/255],[0.,0.,0.]]
#archetypes_bar_plot(niches_cellAb, CELLTYPES,None, y_axis = 'density', radius = RADIUS,path_fig=None)

  return niches_cell_profile

def compute_niche_error(inDir,nbIter,nichesStd):
  sitesCA = pd.read_csv("./"+ inDir +"/sitesCellAb_MIBI_sim"+str(nbIter)+".csv")
  niches = get_niches(sitesCA)
  pd.DataFrame(niches).to_csv("./"+ inDir +"/nichesCellAb_MIBI_sim"+str(nbIter)+".csv")
  errors = []
  for i in range(NBNICHES): # compare cell abundance to reference TLS niche ==> compute RMSE
    err = np.sqrt(np.sum((niches[i,:] - nichesStd[0,:])**2))
    errors.append(err)
    #print(err)
  print("min error: ",str(min(errors)))
  return (1/len(CELLTYPES))*min(errors)

def niche_errors(files):
  nichePrev = float(files.split("_")[0].replace("NEWnichePrev",""))
  nbImg = int(files.split("_")[1].replace("img",""))
  #nbIter = NBITER
  fileJson = open('../AA_sites.json')
  dataNiches = json.load(fileJson)
  #print(dataNiches)
  nichesStd = np.array(dataNiches["nichesCA"])
  errors = []
  for i in range(1,NBITER+1):
    err = compute_niche_error(files,i,nichesStd)
    errors.append(err)
  pd.DataFrame({"simulation_nb":[i for i in range(1,NBITER+1)],"error":errors}).to_csv("./"+files+"/errors_nichePrev"+str(nichePrev)+"_"+str(nbImg)+"img.csv")
  return errors
  



if __name__=="__main__":
  #sitesCA = pd.read_csv("./nichePrev0.0024_5img/sitesCellAb_MIBI_sim5.csv")
  #niches = get_niches(sitesCA)
  #print(niches[0,:])
  # Opening JSON file
  # fileJson = open('../AA_sites.json')
  # 
  # #dataNiches = 
  # dataNiches = json.load(fileJson)
  # #print(dataNiches)
  # nichesStd = np.array(dataNiches["nichesCA"])
  # print(nichesStd[0,:])
  # print(nichesStd.shape)
  # error = compute_niche_error(0.0024,5,2,nichesStd)
  # print(error)
  root='./'
  foldersList = [ item for item in os.listdir(root) if os.path.isdir(os.path.join(root, item)) ]
  print(len(foldersList))
  foldersList = [ item for item in os.listdir(root) if (os.path.isdir(os.path.join(root, item)) and re.match("NEW",item)) ]
  print(len([i for i in foldersList if "NEW" in i])==len(foldersList))
  print(len(foldersList))
  print(foldersList)
  #foldersList = ["NEWnichePrev0.008_32img"]
  # res= niche_errors(files=foldersList[1])
  # print(res)
  
  ##Multiprocessing: compute errors niche TLS over all directories/simulations
  #
  # pool = Pool()
  # args = foldersList #range(1,NITER+1)
  # results = pool.map(niche_errors, args)
  # print(results)
  # pool.close()
  # pool.join()
  
  pool = Pool()
  args = foldersList #range(1,NITER+1)
  results = pool.map(niche_errors, args)
  print(results)
  pool.close()
  pool.join()
  
  
  
