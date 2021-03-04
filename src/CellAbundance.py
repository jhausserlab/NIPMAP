import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


def generate_abundance_matrix(cell_types, patient_ids, n_site, radius, method, snr=1, root="../../output/"):
    return [CellAbundance(p, n_site, radius, cell_types, method=method,snr=snr, root=root) for p in patient_ids]


def join_abundance_matrices(cell_abundances_list):
    sites = []
    patients_ids = []
    for ca in cell_abundances_list:
        sites.extend(ca.abundance_matrix)
        patients_ids.extend([ca.patient_id]*len(ca.abundance_matrix))

    return np.array(sites), np.array(patients_ids)


class CellAbundance:
    def __init__(self, patient_id, n_sites,
                 sites_radius, cell_types, root='',
                 method='abs', snr=1, random_seed=False, image_x_size=800, image_y_size=800):
        if random_seed:
            np.random.seed(random_seed)

        self.patient_id = patient_id
        self.n_sites = n_sites
        self.sites_radius = sites_radius
        self.root = root
        self.cell_positions_df = self._load_cell_position_from_csv()
        self.method = method
        self.snr = snr
        self.k = self.snr * self.snr
        self.pca = None
        self.image_x_size = image_x_size
        self.image_y_size = image_y_size

        self.cell_types = np.array(cell_types)
        self.abundance_matrix = self.calculate_abundace_matrix()

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
    def calculate_cells_count(site, types):
        counts = np.array([np.count_nonzero(site[:, 2] == t) for t in types])
        return counts

    def _load_cell_position_from_csv(self):
        df = pd.read_csv("{}/patient{}_cell_positions.csv".format(self.root, self.patient_id))
        return df

    def calculate_site_groups(self):
        # TODO: to improve, not efficient at all! multiprocessing?
        x = self.cell_positions_df['x'].to_numpy()
        y = self.cell_positions_df['y'].to_numpy()
        t = self.cell_positions_df['cell_type'].to_numpy()
        if self.sites_radius > min(self.image_x_size, self.image_y_size)-self.sites_radius:
            raise ValueError("radius too big!")
        x_centers = np.random.uniform(low=self.sites_radius, high=self.image_x_size-self.sites_radius, size=self.n_sites)
        y_centers = np.random.uniform(low=self.sites_radius, high=self.image_y_size-self.sites_radius, size=self.n_sites)
        sites = {}
        for c_idx in range(x_centers.shape[0]):
            idx = np.where(self.is_in_cirle(x, y, x_centers[c_idx], y_centers[c_idx], self.sites_radius))
            sites[c_idx] = np.array([(x[i], y[i], t[i]) for i in idx])

        return sites

    def calculate_abundace_matrix(self):
        sites = self.calculate_site_groups()
        abundance_matrix = []
        cnt_zeroes = 0
        for site_idx, site in sites.items():
            if len(site) != 0:
                if self.method == 'norm':
                    x = self.calculate_frequencies(site, self.cell_types)
                elif self.method == 'abs':
                    x = self.calculate_cells_count(site, self.cell_types)
                elif self.method == 'abs_log':
                    x = np.log(self.calculate_cells_count(site, self.cell_types) + self.k)
                else:
                    ValueError("Wrong Method! method should be norm or abs or abs_log.")
                abundance_matrix.append(x)
            else:
                cnt_zeroes += 1

        abundance_matrix = np.array(abundance_matrix)

        assert abundance_matrix.shape == (len(sites)-cnt_zeroes, len(self.cell_types))

        return abundance_matrix

    def perform_PCA(self, scale=True):
        self.pca = PCA()
        X = self.abundance_matrix
        if scale:
            X = StandardScaler().fit_transform(X)
        principal_components = self.pca.fit_transform(X)
        return principal_components

