import os
import sys
myDir = os.getcwd()
#print(sys.path)
module_path = myDir + "/TMENS_analysis/"
if module_path not in sys.path:
    sys.path.append(module_path)
#print(sys.path)
import pandas as pd
import numpy as np
from src.utils.archetypes import ArchetypalAnalysis
from src.CellAbundance import CellAbundance, join_abundance_matrices, generate_abundance_matrix
from src.utils.equations import compute_cells_niches_weights,get_niches_cell_abund
from src.utils.visualization import plot_cells_positions
from sklearn.decomposition import PCA
import json
import shutil
from shutil import make_archive



#### ARGUMENTS
CELLTYPES = list(sys.argv[1].split(",")) #['CD8-T', 'Other immune', 'DC / Mono', 'CD3-T', 'B', 'NK', 'Keratin-positive tumor', 'Tumor','CD4-T', 'Mesenchymal-like', 'Macrophages', 'Endothelial', 'Tregs', 'Unidentified', 'DC', 'Mono / Neu','Neutrophils']
ImageIDs = [int(i) for i in sys.argv[2].split(",")] #[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32, 33, 34, 35, 36, 37,38, 39, 40, 41]
NSITES = int(sys.argv[3]) #100
RADIUS = int(sys.argv[4]) #25
NBNICHES = int(sys.argv[5]) #4
METHOD = sys.argv[6] #"gaussian"
ROOT_DATA_PATH = sys.argv[7]  # "./TMENS_analysis/data/cell_positions_data" 
ROOT_OUTPUT_PATH = sys.argv[8] #"./TMENS_analysis/output"
COLARCHS = np.array([[255, 0, 223],[255,0,0],[70,203,236],[0,0,0]]) #TODO  sys.argv[9]

## create figs directory for niches
path_figs = "figs_niches"
path_toFigs = os.path.join(myDir,path_figs)


if __name__ == "__main__":
  #####----- GENERATE SITES AND COMPUTE CELL ABUNDANCE ----#####
  CellAb_list = generate_abundance_matrix(CELLTYPES, ImageIDs, NSITES,RADIUS,method=METHOD, snr=3,center_sites_cells=False,root=ROOT_DATA_PATH)
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
  print(str(NBNICHES)+" found !")
  
  if os.path.isdir(path_toFigs)==False:
    os.mkdir(path_toFigs)
  
  ##########----- SAVE PLOTS IN ZIP DIRECTORY ----##########
  # print("Segmenting images into niches...")
  # for i in ImageIDs:
  #   GRANULARITY = 5
  #   cell_data = pd.read_csv(ROOT_DATA_PATH+"/patient{}_cell_positions.csv".format(i))
  #   fig= plot_cells_positions(cell_data, CELLTYPES, segment_image=True, counting_type=METHOD,
  #                          color_vector=COLARCHS,segmentation_type='colors', granularity=GRANULARITY, radius=RADIUS,
  #                          pca_obj=pc_proj[:,:3], AA_obj=AA, to_plot = 'None',
  #                          path_fig= path_toFigs+"nichesSeg_patient{}.svg".format(i))
  #                          
  # shutil.make_archive("figs_niches","zip", path_toFigs)  
  
  #####----- GENERATE SITES CENTERED ON CELLS AND THEIR NICHE WEIGHTS ----#####
  print("Computing cells' niche weights, the operation might take some time...")
  CellAbCC_list = generate_abundance_matrix(CELLTYPES, ImageIDs, NSITES,RADIUS,method=METHOD, snr=3,center_sites_cells=True,root=ROOT_DATA_PATH)
  sitesCC, patients_ids2,sites_ids2, _ = join_abundance_matrices(CellAbCC_list)
  CellAbCC_df = pd.DataFrame()
  for ca in CellAbCC_list:
    df_ca = ca.get_site_cell_id_df()
    df_ca['patient_id'] = int(ca.patient_id)
    CellAbCC_df = CellAbCC_df.append(df_ca)
  CellAbCC_df = CellAbCC_df.reset_index(drop = True)
  
  NichesProf = get_niches_cell_abund(sitesCellAb=sites,pcaSites=pca_obj,ArchObj=AA,nComp=NBNICHES-1)
  sites_alfa = compute_cells_niches_weights(niches=NichesProf,cellsSites=sitesCC,nbNiches=NBNICHES)
  sites_archs = pd.DataFrame(sites_alfa)
  sites_archs['patient_id'] = patients_ids2
  sites_archs["site_id"] = sites_ids2[:,0]
  sites_archs["cell_type_site"] = sites_ids2[:,1]
  sites_archs["TOT_cell_dens"]= sitesCC.sum(axis=1)
  
  # print(CellAb_df[0:9,: ])
  # print(CellAbCC_df[0:9,: ])
  # print(sites_archs[0:9,: ])
  ##TODO add plot_cells_positions for all images
  ##TODO add plot_cells_TI_border for all images and only for tumor-immune borders (+ one parameter)
  ##TODO directory with figures of segmented images and zip it  
  ##TODO 

                           
 
  ##########----- SAVE OUTPUTS IN CSV ----##########
  ## SAVE PCA AND Archetype Analysis OBJECTS
  print("Saving outputs in json files...")
  dict_pca = {"PC_proj":pc_proj.tolist(),"components":pca_obj.components_.tolist(),"expl_variance":pca_obj.explained_variance_.tolist(),
    "expl_var_ratio":pca_obj.explained_variance_ratio_.tolist(),"mean":pca_obj.mean_.tolist()}
    
  dict_AA = {"archs_coord": AA.archetypes.tolist(), "alfas": AA.alfa.tolist()}
    
  dict_caSites = {"cellAbSites": CellAb_df.to_dict()}
  dict_caSitesCC = {"cells_niches": sites_archs.to_dict(),"cellAb_sitesCC": CellAbCC_df.to_dict()}
    
  # Serializing json objects
  PCA_json = json.dumps(dict_pca, indent=4)
  AA_json = json.dumps(dict_AA, indent=4)
  caSites_json = json.dumps(dict_caSites,indent=4)
  cellsNiches_json = json.dumps(dict_caSitesCC,indent=4)
     
  # Writing to .json files
  with open("./pca_sites.json", "w") as outfile:
    outfile.write(PCA_json)
    
  with open("./AA_sites.json","w") as outfile2:
    outfile2.write(AA_json)
      
  with open("./ca_sites.json","w") as outfile3:
    outfile3.write(caSites_json)
      
  with open("./cells_niches.json","w") as outfile4:
    outfile4.write(cellsNiches_json)

