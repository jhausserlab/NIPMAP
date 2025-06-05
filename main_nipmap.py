import os
import sys
import argparse
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
from src.utils.visualization import plot_cells_positions,radius_pc_all_variance,plot_cells_positions
from sklearn.decomposition import PCA
import json
import shutil
from shutil import make_archive



#### ARGUMENTS
# CELLTYPES = list(sys.argv[1].split(",")) #['CD8-T', 'Other immune', 'DC / Mono', 'CD3-T', 'B', 'NK', 'Keratin-positive tumor', 'Tumor','CD4-T', 'Mesenchymal-like', 'Macrophages', 'Endothelial', 'Tregs', 'Unidentified', 'DC', 'Mono / Neu','Neutrophils']
# ImageIDs = [int(i) for i in sys.argv[2].split(",")] #[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32, 33, 34, 35, 36, 37,38, 39, 40, 41]
# NSITES = int(sys.argv[3]) #100
# RADIUS = int(sys.argv[4]) #25
# NBNICHES = int(sys.argv[5]) #4
# METHOD = int(sys.argv[6]) #"gaussian"
# XSIZE= int(sys.argv[7])
# YSIZE = int(sys.argv[8])
# ROOT_DATA_PATH = sys.argv[9]  #"./TMENS_analysis/data/cell_positions_data"
# ROOT_OUTPUT_PATH = sys.argv[10] #"./TMENS_analysis/output"
# COLARCHS = np.array([[255, 0, 223],[255,0,0],[70,203,236],[0,0,0]]) #TODO  sys.argv[9]
CELLTYPES = ['Other', 'CD4-T', 'Treg', 'Plasma cell', 'B cell', 'pDC', 'CD8-T',
       'Keratin-positive tumor', 'mature DC', 'cDC2', 'cDC1',
       'mature cDC1', 'Folicular DC']
ImageIDs = [1, 2]
NSITES = 1000
RADIUS = 50
NBNICHES = 5 #7 #3 #2 #7 #6#7 #4 #5
METHOD ="gaussian"
XSIZE= 18000
YSIZE = 16000
ROOT_DATA_PATH = "./TMENS_analysis/data/cell_positions_data_mIF"
ROOT_OUTPUT_PATH = "./TMENS_analysis/output"


# COLARCHS = np.array([[255, 0, 223],[255,0,0], [190, 230, 90],[0,0,0]])
COLARCHS = np.array([[255, 0, 223],[255,0,0],[70,203,236], [190, 230, 90],[0,0,0]])
# COLARCHS = np.array([[255, 0, 223],[255,0,0],[70,203,236],[0, 0, 255],[0,0,0]])
# COLARCHS = np.array([[255, 0, 223],[255,0,0],[70,203,236],[0, 0, 255],[255, 165, 0], [0, 255, 0],[0,0,0]])

## create figs directory for niches
path_figs = "figs_niches"
path_toFigs = os.path.join(myDir,path_figs)


def parse_arguments():
    parser = argparse.ArgumentParser(description="NIche Phenotype MAPping (NIPMAP) analysis from spatial multiplex data: Multiplex Ion Imaging on 41 Triple Negative Breast tumors from Keren et al, Cell(2018) and In Situ Sequencing data on human lung development from Sountoulidis et al, 2022")
    parser.add_argument("-mask", "--mask", action="store_true", help="Enable mask")

    return parser.parse_args()


if __name__ == "__main__":
  args = parse_arguments()
  mask = args.mask
  if mask:
    print(f"Using mask from .geojson files in {ROOT_DATA_PATH}")
  '''
  #####----- RADIUS ANALYSIS -----####
  rad = np.linspace(np.log(5), np.log(127), num=10)# get interval of radiuses
  radius = np.rint(np.power(np.e, rad).astype(float))
  METHOD = 'gaussian'

  expl_var_ratio_gauss = {}
  # Get explained variace of each PC on sites cell abundance interating on radius size of sites
  for r in radius:
    print(r)
    gaussian_count_list = generate_abundance_matrix(CELLTYPES, ImageIDs, NSITES, r, center_sites_cells=False,method=METHOD, snr=3, root=ROOT_DATA_PATH, mask = mask)
    sites, patient_ids,s_ids, _ = join_abundance_matrices(gaussian_count_list)
    pca = PCA()
    pc = pca.fit_transform(sites)
    expl_var_ratio_gauss[r] = np.cumsum(pca.explained_variance_ratio_)
  radius_pc_all_variance(expl_var_ratio_gauss,radius_lim=50,nPC_lim=4,cells_number=len(CELLTYPES)+1,save_fig=True, path_fig="./plot_rad_var_gauss.svg", pca_limit=6)
  '''
  print(RADIUS)

  #####----- GENERATE SITES AND COMPUTE CELL ABUNDANCE ----#####
  CellAb_list = generate_abundance_matrix(CELLTYPES, ImageIDs, NSITES,RADIUS,method=METHOD,random_seed=1022,snr=3,center_sites_cells=False,root=ROOT_DATA_PATH, mask = mask)
  sites, patients_ids,sites_ids, _ = join_abundance_matrices(CellAb_list)
  CellAb_df = pd.DataFrame()
  print("Generating sites with cell abundance...")
  for ca in CellAb_list:
    abundance_df = pd.DataFrame(ca.abundance_matrix,columns = CELLTYPES)
    abundance_df['site_id'] = np.arange(len(abundance_df))
    abundance_df['patient_id'] = ca.patient_id
    CellAb_df = pd.concat([CellAb_df, abundance_df], ignore_index=True)
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
    
  niches_cell_profile = get_niches_cell_abund(sites,pca_obj,AA,NBNICHES-1)
  #pd.DataFrame(niches_cell_profile)
  
  if os.path.isdir(path_toFigs)==False:
    os.mkdir(path_toFigs)
  ImageId = [3]
  # ##########----- SAVE PLOTS IN ZIP DIRECTORY ----##########
  '''
  print("Segmenting images into niches...")
  for i in ImageId:
    GRANULARITY = 100 
    cell_data = pd.read_csv(ROOT_DATA_PATH+"/patient{}_cell_positions.csv".format(i))
    fig= plot_cells_positions(cell_data,cell_types=CELLTYPES, segment_image=True, counting_type=METHOD,
                           color_vector=COLARCHS,segmentation_type='colors', granularity=GRANULARITY, radius=RADIUS,
                           h=YSIZE,w=XSIZE,
                           pca_obj=pca_obj, AA_obj=AA, to_plot = 'None',
                           path_fig= path_toFigs+"/"+str(NBNICHES)+"nichesSeg_patient{}.svg".format(i))
  
  # shutil.make_archive("/figs_niches","zip", path_toFigs)
  '''
  #####----- GENERATE SITES CENTERED ON CELLS AND THEIR NICHE WEIGHTS ----#####
  print("Computing cells' niche weights, the operation might take some time...")
  # Call the function (without modifying it)
  CellAbCC_list = generate_abundance_matrix(
      CELLTYPES, ImageIDs, NSITES, RADIUS, method=METHOD, snr=3,
      center_sites_cells=True, border=True, root=ROOT_DATA_PATH, mask=mask
  )
  sitesCC, patients_ids2,sites_ids2, _ = join_abundance_matrices(CellAbCC_list)
  NichesProf = get_niches_cell_abund(sitesCellAb=sites,pcaSites=pca_obj,ArchObj=AA,nComp=NBNICHES-1)
  sites_alfa = compute_cells_niches_weights(niches=NichesProf,cellsSites=sitesCC,nbNiches=NBNICHES)

  sites_archs = pd.DataFrame(sites_alfa)
  sites_archs['SampleID'] = patients_ids2
  sites_archs["cell_id"] = sites_ids2[:,0]
  sites_archs["cell_type"] = sites_ids2[:,1]
  sites_archs["TOT_cell_dens"]= sitesCC.sum(axis=1)
  sites_archs.to_csv("./sites_cells_archs.csv",index=False)
  
  #TODO add plot_cells_TI_border for all images and only for tumor-immune borders (+ one parameter)
  #TODO directory with figures of segmented images and zip it  
  
  ##########----- SAVE OUTPUTS IN CSV ----##########
  ## SAVE PCA AND Archetype Analysis OBJECTS
  print("Saving outputs in json files...")
  dict_pca = {"PC_proj":pc_proj.tolist(),"components":pca_obj.components_.tolist(),"expl_variance":pca_obj.explained_variance_.tolist(),
              "expl_var_ratio":pca_obj.explained_variance_ratio_.tolist(),"mean":pca_obj.mean_.tolist()}

  dict_AA = {"archs_coord": AA.archetypes.tolist(), "alfas": AA.alfa.tolist(),"nichesCA":niches_cell_profile.tolist()}

  dict_caSites = {"cellAbSites": CellAb_df.to_dict()}
  dict_caSitesCC = {"cells_niches": sites_archs.to_dict()}
  dict_params = {"cellTypes":CELLTYPES,"ImageID":ImageIDs,"nbsites":[NSITES],"radiusSize":[RADIUS],"nbniches":[NBNICHES],"countMeth":[METHOD],"xsize":[XSIZE],"ysize":[YSIZE],"rootDataPath":[ROOT_DATA_PATH],"rootOutPath":[ROOT_OUTPUT_PATH],"colNiches":COLARCHS.tolist(),"pathFigs":[path_figs]}

  # Serializing json objects
  PCA_json = json.dumps(dict_pca, indent=4)
  AA_json = json.dumps(dict_AA, indent=4)
  caSites_json = json.dumps(dict_caSites,indent=4)
  cellsNiches_json = json.dumps(dict_caSitesCC,indent=4)
  params_json = json.dumps(dict_params,indent=4)

  # Writing to .json files
  with open("./pca_sites.json", "w") as outfile:
      outfile.write(PCA_json)

  with open("./AA_sites.json","w") as outfile2:
      outfile2.write(AA_json)

  with open("./ca_sites.json","w") as outfile3:
      outfile3.write(caSites_json)

  with open("./cells_niches.json","w") as outfile4:
      outfile4.write(cellsNiches_json)

  with open("./params.json","w") as outfile5:
      outfile5.write(params_json)
