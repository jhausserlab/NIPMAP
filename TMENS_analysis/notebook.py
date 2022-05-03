import pandas as pd
from skimage import io
import matplotlib.pyplot as plt
import numpy as np

DATA_PATH = "../david_data/"

# %%
# SampleID = patientID
# cellLabelInImage = Cell label
# dsDNA = how much dsdna is present
# tumorYN = tumor yes/no
# group =
data = pd.read_csv(DATA_PATH + "cellData.csv")
columns = data.columns

# %%
patients_idx = data['SampleID'].unique()

# %%

p10 = io.imread(DATA_PATH + "configP12General.png")
io.imshow(p10)
plt.show()

# %%

evalvec = np.load(DATA_PATH + "types_inf_Grained_4.npy")

# %%
Patient_ID = 16

sample1 = data[data['SampleID'] == Patient_ID]

p_image = io.imread(DATA_PATH + "p{}_labeledcellData.tiff".format(Patient_ID))
immune_group_image = np.asarray(p_image)

# Building immuneGroup image
_map = {int(row['cellLabelInImage']): int(row['immuneGroup']) for _, row in sample1.iterrows()}
for k, v in _map.items():
    immune_group_image[p_image == k] = v

# %%
import matplotlib.image as mpimg
plt.imshow(immune_group_image, cmap=plt.get_cmap('Paired'))
plt.show()

# %%
from src.utils.visualization import plot_cell_positions

patientID = 35
cell_data = pd.read_csv("../output/cell_positions_data/patient{}_cell_positions.csv".format(patientID))
plot_cell_positions(cell_data, "immune_group")
