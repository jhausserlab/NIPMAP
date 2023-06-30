import os
import sys
myDir = os.getcwd()
#print(sys.path)
module_path = myDir + "/"  #"/TMENS_analysis/"
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
ROOT_DATA_PATH = "./data/cell_positions_data"
ROOT_OUTPUT_PATH = "./output"
COLARCHS = np.array([[255, 0, 223],[255,0,0],[70,203,236],[0,0,0]]) #TODO  sys.argv[9]
## create figs directory for niches
path_figs = "figs_niches"
path_toFigs = os.path.join(myDir,path_figs)

# Sites and image areas to compute the number of sites according to sampling intensity (in%)
As = math.pi * RADIUS**2
Atot = XSIZE * YSIZE

sampling_int = lambda x: x*As/Atot
nSamp = [1,3,10,30,100,300] #sampling intensity in %
NITER= 100 

# Identification of niches
# Generate sampling sites, perform PCA and Archetype Analysis
def get_niches(SamplingInt,j):
  NSITES = int((SamplingInt/100)*Atot/As)
  #print(NSITES)
  ### Generate sites cell abundance
  CellAb_list = generate_abundance_matrix(CELLTYPES, ImageIDs, NSITES,RADIUS,method=METHOD,random_seed=j,snr=3,center_sites_cells=False,root=ROOT_DATA_PATH)#generate_abundance_matrix(CELLTYPES, ImageIDs, NSITES,RADIUS,method=METHOD,random_seed=1022,snr=3,center_sites_cells=False,root=ROOT_DATA_PATH)
  sites, patients_ids,sites_ids, _ = join_abundance_matrices(CellAb_list)
  CellAb_df = pd.DataFrame()
  print("Generating sites with cell abundance...")
  for ca in CellAb_list:
      abundance_df = pd.DataFrame(ca.abundance_matrix,columns = CELLTYPES)
      abundance_df['site_id'] = np.arange(len(abundance_df))
      abundance_df['patient_id'] = ca.patient_id
      CellAb_df = CellAb_df.append(abundance_df)
  CellAb_df = CellAb_df.reset_index()
  
  #####----- PCA ON SITES ABUNDANCE ----#####
  print("Dimension reduction of sites cell abundances...")
  pca_obj = PCA()
  pc_proj = pca_obj.fit_transform(sites)
  #####----- ARCHETYPE ANALYSIS ----#####
  print("Finding niches...")
  AA = ArchetypalAnalysis(n_archetypes = NBNICHES,
                      tolerance = 0.001,
                      max_iter = 200,
                      random_state = 0,
                      C = 0.0001,
                      initialize = 'random',
                      redundancy_try = 30)
  AA.fit_transform(pc_proj[:,:NBNICHES-1])
  print(str(NBNICHES)+" niches found!")
  #get niches cell abundance
  niches_cell_profile = get_niches_cell_abund(sites,pca_obj,AA,NBNICHES-1)
  #print(niches_cell_profile)
  niches_cellAb = pd.DataFrame(niches_cell_profile,columns=CELLTYPES)
  #print(niches_cellAb.head)

  return niches_cell_profile#niches_cellAb

## Compute RMSE
# reference is niches cell abundance from 1000% sampling of images
def compute_rmse(C1000, Csamp,nbIter= 10):
  err = np.sum(np.sqrt(np.sum(C1000- Csamp)**2))
  #rms = mean_squared_error(C1000, Csamp, squared=False)
  return  err/nIter #rms #err/nIter

# Compute niches cell composition of niches from 1000% sampling intensity of images
def get_std_niches(NITER):
  StdNiches = []
  for i in range(NITER):
    Cstd = get_niches(SamplingInt=1000)
    StdNiches.append(Cstd)
  Cstd = np.concatenate(StdNiches)
  pd.DataFrame(Cstd).to_csv("./niches_std1000.csv")
  return Cstd
    

def sampling_rmse(nSamp,Cstd,j):
  
  Csampl = get_niches(SamplingInt=nSamp,j=j)
  print(Cstd.shape)
  print(Csampl.shape)
  # Permutations of niches indices to compare the same niches
  perms = list(permutations(np.arange(0,NBNICHES)))
  rse_perms = []
  for i in range(len(perms)):
    print(Cstd.shape)
    rse = np.sum(Cstd.take(perms[i], axis = 0) - Csampl)**2
    rse_perms.append(np.sqrt(np.sum(rse)))
  #pd.DataFrame({"simulation_nb": j, "RSE": errors}).to_csv("./rse_sampling_mibi_"+str(nSamp)+".csv")
  return min(rse_perms) #np.sum(rse) #np.sum(np.sqrt(np.sum(Cstd[(NBNICHES*j-NBNICHES):(NBNICHES*j),:]- Csampl)**2)) #errors # None

''''
print(Cstd)
#perms = list(permutations([0,1,2,3]))
# Get permutations of niches indices
perms = list(permutations(np.arange(0,NBNICHES)))
print(perms)
print(len(perms))
print(Cstd.take(perms[1], axis = 0))
'''

if __name__=="__main__":
  #Cstd = get_std_niches(NITER)
  Cstd = pd.read_csv(ROOT_OUTPUT_PATH+"/niches_std1000.csv").to_numpy()[0:NBNICHES,1:]

  for i in nSamp:
    pool = Pool() # multiprocessing, 1 thread for one iteration
    args = range(1,NITER+1)
    #results = pool.map(sampling_rmse(nSamp), args)
    #lock = multiprocessing.Lock()
    func = functools.partial(sampling_rmse, i,Cstd)
    errors = []
    results = pool.map(func, args)
    errors.append(results)
    pd.DataFrame({"simulation_nb": args, "RSE": results}).to_csv(ROOT_OUTPUT_PATH+"/rse_sampling_mibi_"+str(i)+".csv")
    #print(results)
    pool.close()
    pool.join()

# nSamp=30
# j = 1
# Csampl = get_niches(SamplingInt=nSamp,j=j)
# pd.DataFrame(Csampl).to_csv("./niches_sampling"+str(nSamp)+".csv")

# Cstd = get_niches(SamplingInt=10)
# Csampl = get_niches(SamplingInt=1)
# rmse = compute_rmse(C1000=Cstd, Csamp=Csampl,nbIter= 10)






