import numpy as np
from qpsolvers import solve_qp
from sklearn.decomposition import PCA

def tumor_misfit_error(original_sites, approx_sites):
    denom = np.sum((original_sites - np.mean(original_sites, axis=0))**2)
    print(denom)
    error = np.sum((np.subtract(approx_sites, original_sites))**2) / denom
    return error


def arch2color(arch):
    d = {0: 'red', 1: 'blue', 2: 'green'}
    return d[arch]


def alfa2rgb(v):
    rgb_vect = (v * 255).astype('uint8')
    return rgb_vect


def alfa2color(colors, v):
    # d = {0: [255, 0, 0], 1: [0, 255, 0], 2: [0, 0, 255], 3: [255, 255, 0]}
    d = {i: c for i, c in enumerate(colors)}
    return np.array(d[np.argmax(v)])


def color_mapper(color_vector, alfa):
    #print(color_vector)
    #print("ALFA:",alfa)
    v = np.array([alfa[i] * color_vector[:, i] for i in range(len(alfa))]).T
    rgb_code = np.sum(v, axis=1)
    #print(rgb_code)
    return rgb_code


def scale(x, min, max, new_min=0, new_max=1):
    return (new_max-new_min) / (max - min) * (x - max) + max


def get_niches_cell_abund(sitesCellAb,pcaSites,ArchObj,nComp):
    niches = np.dot(ArchObj.archetypes.T, pcaSites.components_[:nComp,:])+np.mean(sitesCellAb, axis=0)
    return niches

# compute niches weights for a dataset of sites centered on cells
# Using quadratic programming, find alpha such that:
# X = alpha x B
# with X being the  cell abundance of sites centered on cells,
# and B being the cell abundance of each niche
def compute_cells_niches_weights(niches,cellsSites,nbNiches=4):
    sites_alfas = np.zeros(shape = (cellsSites.shape[0],nbNiches))
    for i in np.arange(cellsSites.shape[0]):
        site = cellsSites[i,:].T
        #print(sites_ids2[i])
        P = np.dot(niches, niches.T)
        q = -np.dot(niches,site)
        A = np.ones(nbNiches)
        b = 1
        G = np.diag(-np.ones(nbNiches))
        h = np.zeros(nbNiches)
        #print(A)
        res = solve_qp(P, q, G, h, A, b)
        sites_alfas[i] = res
        #print("QP solution: x = {}".format(res))
    return sites_alfas

def PCA_sites(X):
    pca_obj = PCA()
    pcs = pca_obj.fit_transform(X)
    return pcs
