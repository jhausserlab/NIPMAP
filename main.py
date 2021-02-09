from src.Cell import group_map
import pandas as pd
import numpy as np
from multiprocessing import Pool
from src.CellsImage import CellsImage


def generate_xy_position_csv(df, patientID):
    patient_data = df[data['SampleID'] == patientID]
    image = CellsImage(patientID)
    print("Getting cells from image for patient: {}".format(patientID))
    group = patient_data['Group']
    immune_group = patient_data['immuneGroup']
    cell_type = [group_map(g, ig) for g, ig in zip(group, immune_group)]
    #print(cell_type)
    _ = image.get_cells(list(patient_data['cellLabelInImage']), cell_type)
    image.to_csv("{}/patient{}_cell_positions.csv".format(OUTPUT_PATH, patientID))

    return 1


if __name__ == "__main__":
    DATA_PATH = "../david_data"
    OUTPUT_PATH = "../output/cell_positions_data"
    data = pd.read_csv("{}/cellData.csv".format(DATA_PATH))
    patients_idx = data['SampleID'].unique()
    to_remove = np.array([42, 43, 44])
    patients_idx = np.setdiff1d(patients_idx, to_remove)
    with Pool(processes=20) as pool:
        results = [pool.apply_async(generate_xy_position_csv, args=(data, patientID)) for patientID in patients_idx]
        results = [p.get() for p in results]

    print("Finished!")

