import numpy as np


class mPCA():
    def __init__(self, centered=True):
        self.centered = centered
        self.eig_vect = None
        self.centered_data = None
        self.columns_mean = None

    def fit_PCA(self, data):
        self.centered_data = data
        if self.centered:
            self.columns_mean = np.mean(data.T, axis=1)
            self.centered_data = data-self.columns_mean
        cov_matrix = np.cov(self.centered_data.T)
        eig_val, eig_vect = np.linalg.eig(cov_matrix)
        perm = np.argsort(-eig_val)
        self.eig_vect = eig_vect[:, perm]

    def fit_transform(self, data):
        self.fit_PCA(data)
        res = self.eig_vect.T.dot(self.centered_data.T)
        return res.T

    def project_data(self, x):
        if self.eig_vect is None:
            raise ValueError("You must fit pca before")

        c_x = x
        if self.centered:
            c_x = x-self.columns_mean
        res = self.eig_vect.T.dot(c_x.T)
        return res.T
