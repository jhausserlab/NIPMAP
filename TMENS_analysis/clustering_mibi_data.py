import os
import sys
module_path = os.path.abspath(os.path.join('../..'))
if module_path not in sys.path:
    sys.path.append(module_path)
os.getcwd()
import numpy as np
import pandas as pd
import math
import matplotlib
from sklearn.decomposition import PCA
from sklearn.metrics import r2_score
from scipy.spatial.distance import cdist,pdist
import matplotlib.pyplot as plt
#sys.path.insert(1,'/scratch/anissa.el/macro_micro_niches/macro_micro_niches2022/TMENS_analysis/')
from src.CellAbundance import CellAbundance, join_abundance_matrices, generate_abundance_matrix
from src.utils.archetypes import ArchetypalAnalysis
from src.utils.visualization import plot_scatter_pca, plot_3Dscatter_pca, archetypes_bar_plot, archetype_simple_plot
from src.utils.equations import arch2color, alfa2rgb, scale, color_mapper
from sklearn.cluster import MiniBatchKMeans
from mpl_toolkits import mplot3d
import random


## CREATE SITES AND cells abundance in each site (gaussian kernel)
#random.seed(10)
CELL_TYPES = ['CD8-T', 'Other immune', 'DC / Mono', 'CD3-T', 'B', 'NK','Keratin-positive tumor', 'Tumor','CD4-T', 'Mesenchymal-like', 'Macrophages', 'Endothelial', 'Tregs', 'Unidentified', 'DC', 'Mono / Neu', 
              'Neutrophils']
patient_ids = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41]
N_SITE = 100
RADIUS = 25
ROOT_DATA_PATH = "./data/cell_positions_data"
METHOD = "gaussian"
SEED = 3
NB_COMP = 17

## GENEREATE SITES CELL ABUNDANCE +PCA
# generate sites cell abundance from images  cell data
# @param patient_id {list} of Image/patient IDs
# @return sites {numpy array} matrix of sites (rows) cell types abundance (columns) 
def generate_sites_ca(patient_ids):
  cell_ab_list = generate_abundance_matrix(CELL_TYPES, patient_ids, N_SITE, 
  RADIUS, method = METHOD, snr=3, center_sites_cells=False,root = ROOT_DATA_PATH)
  sites, patients_ids,sites_ids, _ = join_abundance_matrices(cell_ab_list) 
  
  return sites



##### Kmeans clustering VS TMENs: variance explained
### Number of PCs
## Compute SS cell abund in sites - cell abundfrom n PCs x archetypes prop in each PC
## Compute EV = SS(Model)/TSS
## Plot it 
# @param sites {numpy array} 
# @param n_archs {list} of numbers of archetypes for which we compute expl. var.
# @param pca_pts {numpy array} sites projected onto PCA space
# @param pca_obj {PCA object} PCA object of sites cell abundance
# @returm
def plot_EV(sites, n_archs, pca_pts, pca_obj):
  #n_archs = list(range(2, NB_COMP))#NB_COMP-7-14
  #n_archs=[2]#,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
  #print(np.mean(sites,axis=0),np.mean(sites,axis=0).shape)
  TSS = ((sites - np.mean(sites,axis=0))**2).sum() # The total sum of squares
  #TSS = pca_sites.explained_variance_.sum()
  rows, cols = sites.shape
  print("Total variance",TSS)
  #print("Exp var PCA: ", pca_3d.explained_variance_ratio_)
  ## TMENs
  decenter_func = lambda x: x + np.mean(sites,axis=0)
  exp_var_tmens = []
  exp_var_kmeans=[]
  for n in n_archs:
       print(n)
       # Archetype analysis
       ArchAn = archetypes_analysis_sites(pc_points=pca_pts, nArchs=n)
       predictedCA = np.dot(pca_obj.components_[:, :(n-1)], np.dot(ArchAn.archetypes,ArchAn.alfa))
       #print(predictedCA)
       #print(np.mean(sites,axis=0).shape)
       predCA_cent = decenter_func(predictedCA.T) #predictedCA.T + np.mean(sites,axis=0)
       #print("predicted: ",predCA_cent)
       ESS = (((predCA_cent - np.mean(sites,axis=0))**2).sum()/TSS)*100
       print("EV: ", ESS )
       exp_var_tmens.append(ESS)
       
       ### K MEANS CLUSTERING
       KMeans = MiniBatchKMeans(n_clusters = n,random_state=SEED).fit(sites)
       centroids = KMeans.cluster_centers_
       D_k = cdist(sites, centroids, 'euclidean')
       dists = np.min(D_k,axis=1)
       tot_CSS = sum(dists**2)  # within-cluster SS
       bet_CSS = TSS - tot_CSS  # between clusters SS
       print("EV Kmeans: ",(bet_CSS/TSS)*100)
       exp_var_kmeans.append((bet_CSS/TSS)*100)
     
    ###Scatter plot of Explained variance
  plt.plot(n_archs,exp_var_tmens, '-o')
  plt.plot(n_archs, exp_var_kmeans, '-o',c='red')
  plt.ylabel("% variance explained")
  plt.xlabel("number of niches")
  plt.ylim(0, 100)
  plt.savefig("./output/figs/EV_tmensVSkmeans.svg",format="svg")
  plt.show()
  
  return None  #exp_var_tmens #(exp_var_tmens,exp_var_kmeans)


#### --- Kmeans clustering: 10 clusters (or 9 ) -- ####
#k = 10
#neighborhood_name = "neighborhood"+str(k)

# k-means clustering of sites cell abundance
#@param sites_toclust {numpy array} matrix of sites(rows) x cell types(columns)
#@return labelskm {list} list of cluster labels for each site
def kmeans_clust_sites(sites_toclust):
  k = NB_NEIGHBORHOODS
  k_centroids = {}
  km = MiniBatchKMeans(n_clusters = NB_NEIGHBORHOODS,random_state=SEED)
  labelskm = km.fit_predict(sites_toclust)# pc3d[:, :3]
  k_centroids[k] = km.cluster_centers_
  return labelskm



# computes archetypes from sites cell abundance on PC space
# @param pc_points {numpy array} sites(rows) coordinates of each PC (column)
# @param nArchs {int} number of archetypes= nb of PC +1
# @return AA_3D {Archetype Analysis object} AA object from sites cell abundance
def archetypes_analysis_sites(pc_points, nArchs):
  ### ARCHETYPE ANALYSIS
  AA_3D = ArchetypalAnalysis(n_archetypes = nArchs, 
                          tolerance = 0.001, 
                          max_iter = 1000,#200, 
                          random_state = 0, 
                          C = 0.0001, 
                          initialize = 'random',
                          redundancy_try = 30)
  AA_3D.fit_transform(pc_points[:, :nArchs-1])

  return AA_3D
#pd.DataFrame(AA_3D.alfa).to_csv("./outputs/sites_neighborhoods_alfa.csv")



#### Compare AIC niches VS kmeans
def plot_aic_niches_kmeans(sites,pca_pts,pca_obj,TSS):
  n_archs = list(range(2,len(CELL_TYPES)+1))
  aic_niches =[]
  aic_kmeans = []
  for n in n_archs:
    print(n)
    # Archetype analysis
    ArchAn = archetypes_analysis_sites(pc_points=pca_pts, nArchs=n)
    predictedCA = np.dot(pca_obj.components_[:, :(n-1)], np.dot(ArchAn.archetypes,ArchAn.alfa))

    predCA_cent = decenter_func(predictedCA.T) #predictedCA.T + np.mean(sites,axis=0)
    ESS = (((predCA_cent - np.mean(sites,axis=0))**2).sum()/TSS)*100
    print("EV: ", ESS )
    RSS = TSS - ((predCA_cent - np.mean(sites,axis=0))**2).sum()
    print("RSS: ",RSS)
    pniches = n*(n-1) +sites.shape[0]*(n-1) 
    #pniches= n*(n-1)*len(CELL_TYPES)
    freeParam_pca = list(range(len(CELL_TYPES)-(n-2),len(CELL_TYPES)+1))
    #print(freeParam_pca)
    #pniches = n*(sum(freeParam_pca)) ## PCA orthogonality ==> each component, we loose 1 dimension
    #pniches =n*((len(CELL_TYPES)*len(CELL_TYPES)+1)/2)
    pniches = n*(n-1) + len(CELL_TYPES) + (n-1)*(len(CELL_TYPES)-2)
    pniches = n*(n-1) + len(CELL_TYPES) + len(CELL_TYPES)*(len(CELL_TYPES)+1)/2
    pniches = n*(n-1) + len(CELL_TYPES) + sum(freeParam_pca)
    AICn = 2*pniches + sites.shape[0]*math.log(RSS)
    print("AIC niches: ", AICn)
    aic_niches.append(AICn)
    
    ### K MEANS CLUSTERING
    KMeans = MiniBatchKMeans(n_clusters = n,random_state=SEED).fit(sites)
    centroids = KMeans.cluster_centers_
    D_k = cdist(sites, centroids, 'euclidean')
    dists = np.min(D_k,axis=1)
    tot_CSS = sum(dists**2)  # within-cluster SS not explained by model
    bet_CSS = TSS - tot_CSS  # between clusters SS explained by k-means
    print("EV Kmeans: ",(bet_CSS/TSS)*100)
    print("RSS: ",tot_CSS)
    pclust = n*len(CELL_TYPES) + sites.shape[0]
    pclust = n*len(CELL_TYPES)
    AICk = 2*pclust +sites.shape[0]*math.log(tot_CSS)
    AICk = pclust + 0.5*tot_CSS
    print("AIC kmeans: ", AICk)
    aic_kmeans.append(AICk)
  print(aic_niches)
  print(aic_kmeans)
  plt.plot(n_archs,aic_niches, '-o')
  plt.plot(n_archs, aic_kmeans, '-o',c='red')
  plt.ylabel("AIC")
  plt.xlabel("number of niches")
  #plt.ylim(0, 100)
  plt.savefig("./output/figs/AIC_tmensVSkmeans4.png",format="png")
  plt.show()
  return None
       



##### PLOT ######

# Plots sites in PCA on of cell abundance colored by kmeans clusters
# previously found by sites cell abundance
# @param pc3d {numpy array} matrix of sites(rows) x PC(columns)
# @param km_labs {numpy array} matrix of sites(rows) x cluster labels
# @param pca_sites_EV {list} of explained variance (%) for each PC
# @return None
def plot_sites_kmeans(pc3d, km_labs,pca_sites_EV):#archetypes= AA_3D.archetypes, 
  fig = plt.figure()
  ax = plt.axes(projection='3d')
  ax = plt.axes(projection='3d')
  color_map = plt.get_cmap('tab10')
  
  scatter_plot=ax.scatter3D(pc3d[:,0],pc3d[:,1],pc3d[:,2],alpha=0.6,c=km_labs,cmap = color_map)#
  #av_clusters = ax.scatter3D(av_clust_proj[:,0],av_clust_proj[:,1],av_clust_proj[:,2],
  #color = "red", marker = "D",alpha=1)
  #archetypes3d = ax.scatter(AA_3D.archetypes[0, :],AA_3D.archetypes[1, :],AA_3D.archetypes[2, :],c="black",s=50,marker="d",alpha=1)
  ax.set_xlabel("PC1 ({} %)".format(round(pca_sites_EV[0]*100,1)))
  ax.set_ylabel("PC2 ({} %)".format(round(pca_sites_EV[1]*100,1)))
  ax.set_zlabel("PC3 ({} %)".format(round(pca_sites_EV[2]*100,1)))
  ax.set_xticklabels([])
  ax.set_yticklabels([])
  ax.set_zticklabels([])
  
  plt.colorbar(scatter_plot)
  plt.show()



################------- MAIN -------################

if __name__ == '__main__':
  NB_NEIGHBORHOODS = 10
  sites = generate_sites_ca(patient_ids=patient_ids)
  pca_obj = PCA(n_components=NB_COMP)
  pca_pts = pca_obj.fit_transform(sites)
  
  labelskm = kmeans_clust_sites(sites_toclust =sites)
  sites_cells = pd.DataFrame(sites,columns = CELL_TYPES)
  sites_cells['neighborhood_id'] = labelskm # CLUSTERS IDS for each site
  #archetypes_analysis_sites(pc3d = pc3d)

  plot_sites_kmeans(pc3d=pca_pts, km_labs =labelskm,pca_sites_EV = pca_obj.explained_variance_ratio_)
  
  n_archs = list(range(2,18))
  print(n_archs)
  print("TSS: ", ((sites - np.mean(sites,axis=0))**2).sum()) 
  #TSS = ((sites - np.mean(sites,axis=0))**2).sum() 
  #print("Variance explained * 4000 nb sites: ",pca_obj.explained_variance_.sum() *4000)
  print(pca_obj.explained_variance_ratio_)
  
  plot_EV(sites, n_archs, pca_pts, pca_obj) #(sites, n_archs, pc3d_sites, pca_sites)
  #plot_aic_niches_kmeans(sites,pca_pts,pca_obj,TSS)




