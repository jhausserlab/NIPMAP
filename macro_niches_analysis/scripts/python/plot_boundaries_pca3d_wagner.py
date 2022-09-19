import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.spatial import HalfspaceIntersection
import matplotlib.patches as mpatches
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import ConvexHull
import mpl_toolkits.mplot3d as a3
import matplotlib.colors
from sys import path
path.append("/srv/mfs/hausserlab/fabio/data_analysis/src/utils/")
from archetypes import ArchetypalAnalysis

filename1 = "/srv/mfs/hausserlab/anissa.el/ImmuneStates/wagner_analysis/PCA_Wagner_3D.csv"
filename2 = "/srv/mfs/hausserlab/anissa.el/ImmuneStates/wagner_analysis/PCA_Wagner_components.csv" 
filename3 = "/srv/mfs/hausserlab/anissa.el/ImmuneStates/wagner_analysis/PCA_Wagner_means.csv"
pcData = pd.read_csv(filename1) #"PCA_Wagner_3D.csv"
pc3dWagner = np.array(pcData, dtype = "float64")
pca_3d_components = pd.read_csv(filename2)# Scores of PCs
n_components = 3 #Number of kept PCs
n_of_std = 1.5 ## ?? Good question
X = pc3dWagner #Our obseervations projected on the 3PCs space
V = np.array(pca_3d_components, dtype = "float64")
pca_means = pd.read_csv(filename3)# List of means (one for each cell type) used for centering in PCA
#print(pca_means)
means = np.array([pca_means["V1"][i] for i in range(0,10)],dtype = "float64")
#print(V.shape)
print(means.shape)
## Compute standard deviation for the PCs that were not kept (from the 3rd to the last PC)
std_vx_last = np.matmul(V[:, n_components:], X.T[n_components:,:]).std(axis=1)
std_vx = np.matmul(V, X.T).std(axis=1)
### FIRST CONSTRAINT: Non-negativity of xi (percentages of cell types)
# Forst compute A and bb from the forst halfspace such as AX + b >=0
A = -V[:, :3]
b = -means-n_of_std*std_vx_last

ones = np.array([1,1,1,1,1,1,1,1,1,1])
std_vx_area = np.dot(ones, np.matmul(V[:, n_components:], X.T[n_components:,:])).std()
#### SECOND CONSRTRAINT: Sum of xi  for each observation must be equal to 100
A_bis = np.dot(ones.T, V[:, :n_components])#V[:, :3]#
b_bis =  np.dot(ones.T, means) - 100 - n_of_std*std_vx_area
print("A and b shapes: ",A.shape,b.shape)
print("A bis shape: ",A_bis.shape)
#print(np.c_[A,b].shape)
#print(np.append(A_bis, b_bis))
halfspaces = np.vstack((np.c_[A,b], np.append(A_bis, b_bis))).astype('double')#np.vstack((np.c_[A,b], np.append(A_bis, b_bis))).astype('double')
print(halfspaces.shape)
print("v shape: ",V[:, :n_components].shape )
print("Origin point:",np.zeros(10).reshape(10,1).T, np.zeros(10).reshape(10,1).T.shape)
zeros = np.array([0,0,0,0,0,0,0,0,0,0])
feasible_point = np.matmul(zeros.reshape(10,1).T- means,V[:, :n_components]).T #
print("feasible point: ",feasible_point[:,0])
hs = HalfspaceIntersection(halfspaces, feasible_point[:,0])

def plot_3Dscatter_pca(principal_components, evr, labels, cmap=plt.cm.nipy_spectral, archetypes=None, halfspaces=None, original_axis = None, cell_type=None):
    #tot_evr = np.sum(evr[0:3])
    #print("{:.2f}% Total Exp. Var.".format(tot_evr))
    fig = plt.figure(1, figsize=(8, 8))
    plt.clf()
    ax = Axes3D(fig, rect=(0, 0, 1, 1), elev=21, azim=-58)

    #minimum, maximum = _calculate_pca_max_min(principal_components[:, 0:3])
    #ax.set_xlim(minimum, maximum)
    #ax.set_ylim(minimum, maximum)
    #ax.set_zlim(minimum, maximum)
    ax.set_xlabel("PC1 ")#"PC1 ({:.2f}% Exp. Var.)".format(evr[0] * 100))
    ax.set_ylabel("PC2 ")#"PC2 ({:.2f}% Exp. Var.)".format(evr[1] * 100))
    ax.set_zlabel("PC3 ")#"PC3 ({:.2f}% Exp. Var.)".format(evr[2] * 100))

    alpha = 0.1 if halfspaces is not None else 1
    alpha = 0.05 if original_axis is not None else 1
    ax.scatter(principal_components[:, 0], principal_components[:, 1], principal_components[:, 2],
               c="blue", cmap=cmap, edgecolor='k', s=30,alpha = alpha)

    if archetypes is not None:
        ax.scatter(archetypes[0, :], archetypes[1, :], archetypes[2, :],
                   marker = "^", s = 100, color='orange')

    if original_axis is not None:
        origin = original_axis[len(original_axis)-1, :]
        for i in range(len(original_axis)-1):
            plt.plot([origin[0], original_axis[i, 0]], [origin[1], original_axis[i, 1]], [origin[2], original_axis[i, 2]], 'k-', color='b', alpha=0.5)
            ax.text(original_axis[i, 0], original_axis[i, 1], original_axis[i, 2], cell_type[i], color='red')

    if halfspaces is not None:
        x, A, b = halfspaces
        verts = x.intersections
        hull = ConvexHull(verts)
        faces = hull.simplices

        for i, s in enumerate(faces):
            sq = [
                [verts[s[0], 0], verts[s[0], 1], verts[s[0], 2]],
                [verts[s[1], 0], verts[s[1], 1], verts[s[1], 2]],
                [verts[s[2], 0], verts[s[2], 1], verts[s[2], 2]]
            ]

            f = a3.art3d.Poly3DCollection([sq])

            eps=1e-6
            if np.linalg.norm(np.dot(A, np.array(sq).T) + b) < eps:
                f.set_color("#ff0000")
                #f.set_alpha(0.1)
            else:
                f.set_color(matplotlib.colors.rgb2hex(np.random.rand(3)))
                #f.set_alpha(0.2)
            f.set_edgecolor([(0,0.5,0,0.8)])#k
            f.set_alpha(0.1)
            ax.add_collection3d(f)

    return plt
    
    
if __name__ == "__main__":
  
  plt = plot_3Dscatter_pca(pc3dWagner, evr=0.93, labels=np.zeros(12), archetypes=None, halfspaces=(hs,A_bis,b_bis))
  #plt.savefig("../../../output/grants_images/3d_gaussian_boudaries_plot.svg")
  plt.show()
