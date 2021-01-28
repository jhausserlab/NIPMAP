import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from src.Cell import Cell

ROOT_DATA_PATH = "../../david_data"


def _get_centroid(x_indices, y_indices, N):
    return np.sum(x_indices)/N, np.sum(y_indices)/N


def _get_boundingbox_center(x_indices, y_indices):
    return (np.min(x_indices) + np.max(x_indices))/2, \
           (np.min(y_indices) + np.max(y_indices))/2


class CellsImage:
    def __init__(self, patientID, micrometer=True):
        self.patientID = patientID
        self.labels_matrix = self._load_image()
        self.cells = None
        self.micrometer = micrometer

    def _load_image(self):
        return plt.imread("{}/p{}_labeledcellData.tiff".format(ROOT_DATA_PATH, self.patientID))

    def _get_single_cell(self, cell_label, cell_type, center='bb_center'):
        if center not in ['bb_center', 'centroid']:
            raise ValueError('Wrong center. Center parameters must be one of '
                             'the following values: bb_center, centroid.')

        bin_matrix = self.labels_matrix == cell_label
        indices = np.array(np.where(bin_matrix)).astype(float)

        if self.micrometer:
            indices[0] = indices[0] * 0.39
            indices[1] = indices[1] * 0.39
        if center == 'bb_center':
            x, y = _get_boundingbox_center(indices[0], indices[1])
        else:
            N_pixels = np.count_nonzero(bin_matrix)
            x, y = _get_centroid(indices[0], indices[1], N_pixels)

        # TODO: human readable cell type
        return Cell(x, y, cell_label, cell_type)

    def get_cells(self, cells_idx, cell_types):
        if not len(cells_idx) == len(cell_types):
            raise ValueError("cells_idx and cell type lists must have the same length")

        if self.cells is None:
            self.cells = [self._get_single_cell(cell_id, ct)
                          for cell_id, ct in zip(cells_idx, cell_types)]

        return self.cells

    def get_colored_matrix(self, cells_idx, cell_types, id_map):
        if not len(cells_idx) == len(cell_types):
            raise ValueError("cells_idx and cell type lists must have the same length")

        cell_type_map = np.copy(self.labels_matrix)
        for cell_id, ct in zip(cells_idx, cell_types):
            cell_type_map[cell_type_map == cell_id] = id_map[ct]

        return cell_type_map

    def to_csv(self, filename):
        if self.cells is not None:
            df = pd.DataFrame([c.__dict__ for c in self.cells])
            df.to_csv(filename, header=True, index=False, sep=',')
        else:
            raise ValueError("Variable cells is None, try to call get_cells method before.")






