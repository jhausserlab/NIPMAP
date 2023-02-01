import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import KDTree
import matplotlib.pyplot as plt


def generate_abundance_matrix(cell_types, patient_ids, n_site, radius, method, snr=1, root="../../output/",image_x_size=800, image_y_size=800,
                              random_seed=False,center_sites_cells=False,border=False):
    return [CellAbundance(p, n_site, radius, cell_types, method=method, snr=snr, root=root,image_x_size=image_x_size, image_y_size=image_y_size, random_seed=random_seed,center_sites_cells=center_sites_cells,border=border) for
            p in patient_ids]


def join_abundance_matrices(cell_abundances_list,center_sites_cells=False):
    sites = []
    patients_ids = []
    gradients = []
    cell_sites_ids = []
    for ca in cell_abundances_list:
        sites.extend(ca.abundance_matrix)
        patients_ids.extend([ca.patient_id] * len(ca.abundance_matrix))
        cell_sites_ids.extend(ca.sites_cell_ids)
        gradients.extend(ca.gradient)
    if center_sites_cells ==False:
      return np.array(sites), np.array(patients_ids),np.array(cell_sites_ids), np.array(gradients)
    else:
      return np.array(sites), np.array(patients_ids),pd.DataFrame(cell_sites_ids,columns=["site_id","cell_type_site"]), np.array(gradients)


class CellAbundance:
    def __init__(self, patient_id, n_sites,
                 sites_radius, cell_types, root='',
                 method='abs', norm_method=None, snr=1, random_seed=False, image_x_size=800, image_y_size=800,center_sites_cells=False,border=False):
        if random_seed:
            np.random.seed(random_seed)

        self.patient_id = patient_id
        self.n_sites = n_sites
        self.sites_radius = sites_radius
        self.root = root
        self.cell_positions_df = self._load_cell_position_from_csv()
        self.method = method
        self.norm_method = norm_method
        self.snr = snr
        self.k = self.snr * self.snr
        self.pca = None
        self.image_x_size = image_x_size
        self.image_y_size = image_y_size
        self.center_sites_cells = center_sites_cells
        self.border = border
        # for gaussian weigthing we could select a bigger radius
        self.radius_coeff = 2
        self.sites_cell_ids = self.calculate_site_groups().keys()
        self.cell_types = np.array(cell_types)
        self.abundance_matrix, self.gradient = self.calculate_abundace_matrix()

    @staticmethod
    def is_in_cirle(x, y, center_x, center_y, rad):
        return (x - center_x) * (x - center_x) + (y - center_y) * (y - center_y) <= rad * rad

    @staticmethod
    def calculate_frequencies(site, types):
        counts = CellAbundance.calculate_cells_count(site, types)
        freq = counts / np.sum(counts)
        assert freq.shape == types.shape
        return freq

    @staticmethod
    def calculate_cells_count(site, types, radius):
        counts = np.array([np.count_nonzero(site[:, 2] == t) for t in types]) #why we normalized? -> / (np.pi*radius**2)
        return counts

    @staticmethod
    def calculate_gaussian_count(site, types, c, radius):
        '''
        @param: site {numpy array} cells in a site 
        @param: types {list} cell types
        @param: c {numpy array} x,y positions of center of site
        @param: radius {double} bandwidth for density estimation
        @return: density for each cell type within site
        '''
        #score = CellAbundance.get_gaussian_weights(site[:, 0], site[:, 1], site[:, 4], site[:, 5], radius)
        #c = site[0, 4:6].astype(float)
        # x,y coordinates of cells that are in the input site
        ps = site[:, :2].astype(float)
        #x =-((ps-c)**2).sum(1)/(2*radius**2)
        #print(x.dtype)
        # c stands for the x,y coordinates of the center of the site 
        score = np.exp(-((ps-c)**2).sum(1)/(2*radius**2)) / (2*np.pi*radius**2)
        gradientScore = ((ps - c) / (radius ** 2)) * score[:, np.newaxis]
        # print(gradientScore.shape) 9, 2

        #print(gradientScore)

        counts = np.array([np.sum(score[np.where(site[:, 2] == t)[0]]) for t in types])
        #print(counts)
        # gradient calculation
        gradient = np.array([np.sum(gradientScore[np.where(site[:, 2] == t)[0], :], axis=0) for t in types])
        #print(gradient)
        return counts, np.linalg.norm(gradient, axis=1)

    #@staticmethod
    #def get_gaussian_weights(x, y, x_centers, y_centers, radius):
    #    d = np.sqrt(((x-x_centers)**2 + (y-y_centers)**2).astype(float))
    #    gaussian_weigths = np.exp(-(d**2)/(2*radius**2))/(2*np.pi*radius**2)
    #    return gaussian_weigths


    def _load_cell_position_from_csv(self):
        '''
        creates a df from csv file.
        Contains cells data in an image (patient):
        x,y coordinates,label and cell type
        '''
        df = pd.read_csv("{}/patient{}_cell_positions.csv".format(self.root, self.patient_id))

        # Load the image size (for Nature paper)
        if 'x_size' in df.columns and 'y_size' in df.columns:
            self.image_x_size = df['x_size'][0]
            self.image_y_size = df['y_size'][0]
        return df

    def calculate_site_groups(self):
        """
        generates sites given a number 
        :@return: sites list
        """
        # TODO: to improve, not efficient at all! multiprocessing?
        x = self.cell_positions_df['x'].to_numpy().astype(float)
        y = self.cell_positions_df['y'].to_numpy().astype(float)
        t = self.cell_positions_df['cell_type'].to_numpy()
        cell_id = self.cell_positions_df['label']
        #print(cell_id.tolist())
        radius = self.sites_radius
        if self.method == "gaussian":
            radius = self.radius_coeff*radius

        if radius > min(self.image_x_size, self.image_y_size) - radius: # probably twice radius
            raise ValueError("radius too big!")
        # Center sites on cells or generate randomly(uniform) 100 sites per image    
        if self.center_sites_cells == False:
            x_centers = np.random.uniform(low=radius, high=self.image_x_size - radius,size=self.n_sites)
            y_centers = np.random.uniform(low=radius, high=self.image_y_size - radius,size=self.n_sites)
            centers = np.stack((x_centers, y_centers), axis=-1)
            points = np.stack((x, y), axis=-1)
            tree = KDTree(points, leaf_size=5)
            idx_b = tree.query_radius(centers, r=radius, count_only=False)
            sites = {c_idx: (np.stack((x[idx_b[c_idx]], y[idx_b[c_idx]], t[idx_b[c_idx]], cell_id[idx_b[c_idx]],
                                      [x_centers[c_idx]]*len(idx_b[c_idx]), [y_centers[c_idx]]*len(idx_b[c_idx])), axis=-1),
                             centers[c_idx])
                     for c_idx in range(x_centers.shape[0])}
        else:
            #print(self.border)
            
            if self.border==False:
                # Filter cells that are in the border of the image
                threshX = self.image_x_size - radius
                threshY = self.image_y_size - radius
                cells_positions = self.cell_positions_df.loc[((self.cell_positions_df['x']<= threshX) & (self.cell_positions_df['x']>=radius)) & ((self.cell_positions_df['y']<= threshY) & (self.cell_positions_df['y']>=radius))]
            #print(self.cell_positions_df.shape, cells_positions.shape)
            #x = self.cell_positions_df['x'].to_numpy().astype(float)
        #y = self.cell_positions_df['y'].to_numpy().astype(float)
            else:
                cells_positions = self.cell_positions_df 
#                cells_positions = self.cell_positions_df.loc[((self.cell_positions_df['x']<= self.image_x_size) & #(self.cell_positions_df['x']>=0)) & ((self.cell_positions_df['y']<= self.image_y_size) & (self.cell_positions_df['y']>=0))]
                
            x_centers = cells_positions['x'].to_numpy().astype(float)
            y_centers = cells_positions['y'].to_numpy().astype(float)
            
            t_center = cells_positions['cell_type'].to_numpy()
            cell_id_center = cells_positions['label']
            centers = np.stack((x_centers, y_centers), axis=-1)
            #print(x_centers.shape[0])
            points = np.stack((x, y), axis=-1)
            tree = KDTree(points, leaf_size=5)
            idx_b = tree.query_radius(centers, r=radius, count_only=False)
            sites = {(cell_id_center.tolist()[c_idx],t_center.tolist()[c_idx]): (np.stack((x[idx_b[c_idx]], y[idx_b[c_idx]], t[idx_b[c_idx]], cell_id[idx_b[c_idx]],
                                      [x_centers[c_idx]]*len(idx_b[c_idx]), [y_centers[c_idx]]*len(idx_b[c_idx])), axis=-1),
                             centers[c_idx])
                     for c_idx in range(x_centers.shape[0])}

        return sites

    def get_site_cell_map_dataframe(self):
        sites = self.calculate_site_groups()
        patient_cell_df = pd.DataFrame()
        for site_idx, (site,center) in sites.items():
            patient_cell_df = patient_cell_df.append(pd.DataFrame(
                {'site_id': site_idx,'site_x_centers': site[:, 4],  'site_y_centers': site[:, 5], 'x_cell': site[:, 0],
                 'y_cell': site[:, 1], 'cell_type': site[:, 2], 'cell_id': site[:, 3]}))
        return patient_cell_df
    
    def get_site_cell_id_df(self):
        sites = self.calculate_site_groups()
        patient_cell_df = pd.DataFrame()
        for c in sites.keys():
            #print(c)
            #print(pd.DataFrame({'cell_id': [c[0]], 'cell_type': [c[1]]}))
            patient_cell_df = patient_cell_df.append(pd.DataFrame({'cell_id': [c[0]], 'cell_type': [c[1]]}))
        return patient_cell_df

    def calculate_abundace_matrix(self):
        """
        calculates the sites abundance matrix, the counting is based on the argument passed in the constructor of the class
        :@return: {numpy array}(n_sites, #cell_types) matrix containing the abundance of cells for each site. In case of gaussian counting
        method selected, the gradient matrix is also returned.
        """
        sites = self.calculate_site_groups()
        abundance_matrix = []
        cnt_zeroes = 0
        gradients = []
        for site_idx, (site, center) in sites.items():
            #if len(site) != 0:
            if self.method == 'norm':
                x = self.calculate_frequencies(site, self.cell_types)
            elif self.method == 'abs':
                x = self.calculate_cells_count(site, self.cell_types, self.sites_radius)
            elif self.method == 'abs_log':
                x = np.log(self.calculate_cells_count(site, self.cell_types) + self.k)
            elif self.method == 'gaussian':
                x, gradient = self.calculate_gaussian_count(site, self.cell_types, center, self.sites_radius)
                gradients.append(gradient)
            else:
                ValueError("Wrong Method! method should be norm or abs or abs_log.")
            abundance_matrix.append(x)
            #else:
            #    cnt_zeroes += 1
        abundance_matrix = np.array(abundance_matrix)
        #if self.norm_method is not None:
        #    if self.norm_method == 'avg':
        #        print(abundance_matrix.mean())
        #        abundance_matrix = abundance_matrix / abundance_matrix.mean(0)[None, :]

        if abundance_matrix.shape[0] == 0:
            assert abundance_matrix.shape[0] == len(sites) - cnt_zeroes
        else:
            assert abundance_matrix.shape == (len(sites) - cnt_zeroes, len(self.cell_types))

        return abundance_matrix, gradients

    def perform_PCA(self, scale=True):
        self.pca = PCA()
        X = self.abundance_matrix
        if scale:
            X = StandardScaler().fit_transform(X)
        principal_components = self.pca.fit_transform(X)
        return principal_components
