# README

Analysis of Tumor MicroEnvironmental Niches (TMENs) from images of samples of Breast Cancer (BC).
Data used from [Keren et al, Cell(2018)](10.1016/j.cell.2018.08.039)

## Prerequisites

Jupyter notebook installed

Python 3.8 installed
List of python libraries to install 

```bash
pip install matplotlib
pip install scipy
pip install pandas
pip install numpy
pip install scikit-learn
pip install seaborn
pip install qpsolvers
```

## FILES

```
└───/main_nipmap.py: calls py functions for NIPMAP pipeline + saves outputs into .json files
└───/py_wrapper_nipmap.r: get NIPMAP outputs from .json files and loads R objects for further analyses
|
/macro_niches_analysis
|
└───/data
|
└───/figs
|
└───/outputs
|
└───/scripts
/phenotype_niches
|
/TMENS_analysis
|   README.md
|
└───/data
|   |   cellData.csv: normalized data (markers intensity) from CyTOF expfrom Keren et al.,Cell(2018)
|   └───/tnbc_nature_cancer_dataset: contains csv files of cell positions and label for each patient (patientID_cell_position.csv)
|   |
|   └───/cell_positions_data: .csv files for that contain for each patientx,y positions and label of cells from Keren et al.,Cell(2018)
|
└───/notebooks
|   |
|   └───/keren_tmens_analysis
|   |   |   archetype_analysis.ipynb
|   |   |   boundaries_analysis.ipynb
|   |   |   cell_positions_visualization.ipynb
|   |   |   csv_reading.ipynb
|   |   |   distribution_cell_types.ipynb
|   |   |   low_abundance_archetype_analysis.ipynb
|   |   |   misfit_tumor_errors.ipynb
|   |   |   NNMF.ipynb 
|   |   |   normalization.ipynb
|   |   |   PCA_all_tumors_sites_together.ipynb
|   |   |   PCA_log_abundance_analysis.ipynb
|   |   |   PCA_tumor_site_analysis.ipynb
|   |   |   performance_testing.ipynb
|   |   |   site_radius_analysis.ipynb
|   |   |   sites_patient_mapper.ipynb
|   |   |   TME_architecture_analysis
|		|
|   └───/nature_cancer_building_bloks_analysis
|       |   adapt_csv.ipynb
|       |   cells_positions.ipynb
|       |   gaussian_nature_archetype_analysis.ipynb
|       |   gaussian_nature_cells_positions.ipynb
|       |   nature_all_tumors_analysis.ipynb
|       |   nature_cells_positions.ipynb
|       |   nature_PCA_site_analysis.ipynb 
|       |   nature_radius_analysis.ipynb
|       |   nature_sites_patient_mapper.ipynb
|       |   nature_tme_arch_analysis.ipynb 
|
└───/output
|   |
|   └───/csv_files_nature
| 
└───/src
    |	CellAbundance.py
    |   Cell.py
    |   CellsImage.py
    |   mPCA.py
    └───/utils
        |   archetypes.py
        |   equations.py
        |   visualization.py
```


## License


## Contact
Fabio Lipreri - 
Anissa El Marrahi - <anissa.el@scilifelab.se>
Jean Hausser - <jean.hausser@scilifelab.se>

## Acknowledgments
Fabio Lipreri
David Alber
