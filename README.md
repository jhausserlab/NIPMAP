# README

NIche Phenotype MAPping (NIPMAP) analysis from spatial multiplex data: Multiplex Ion Imaging on 41 Triple Negative Breast tumors from [Keren et al, Cell(2018)](10.1016/j.cell.2018.08.039) and In Situ Sequencing data on human lung development from [Sountoulidis et al, 2022](https://doi.org/10.1101/2022.01.11.475631)

## Prerequisites

* Jupyter notebook installed

* Python 3.7 installed
    List of python libraries to install 

    ```bash
    pip install matplotlib
    pip install scipy
    pip install pandas
    pip install numpy
    pip install scikit-learn
    pip install seaborn
    pip install matplotlib
    pip install qpsolvers
    ```
* R 4.1.3 with RStudio installed
    List of R libraries to install 
    ```
    pckgs <- c("tidyverse","ggplot2","ade4","factoextra","plotly","igraph","reshape2","ggrepel","viridis","fdrtool","pheatmap","cluster","broom","pROC","ggpubr","devtools","ggridges")
    install.packages(pckgs)
    ```
    **//!\\ CAUTION** reticulate library version 1.22 should be downloaded, not the latest !! 

## Quick start
The input is a CSV file where rows represent cells and columns represent cell_type, x position, y position, tissue_sample_of_origin, XXX. See XXX.csv for an example.

Geometrically, the samples should have the shape of a rectangle to generate sites and compute cell abundance. Cell types must have been already identified prior to NIPMAP analysis: NIPMAP considers cell types as a given. The list of markers should include phenotypic markers, that is markers (proteins or other) that were not used to identify cell types.

Run script called nipmap.r with your own data. Witin this script, adjust parameters (number of niches, size of sampling sites, number of sites, ...). 
Running the analysis will produce the following outputs:
* Generation of sites and cellular abundances within them (+ radius size selection)
* PCA and Archetype Analysis
* Niche identification
* Niche segmentation of images
* Niche-phenotype mappping

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
└───/data
|
└───/figs
|
└───/outputs
|
/TMENS_analysis
|   README.md
|
└───/data
|   |   cellData.csv: normalized data (markers intensity) from CyTOF expfrom Keren et al.,Cell(2018)
|   └───/tnbc_nature_cancer_dataset: contains csv files of cell positions and label for each patient (patientID_cell_position.csv)
|   |
|   └───/cell_positions_data: .csv files for that contain for each patientx,y positions and label of cells from Keren et al.,Cell(2018)
|   |
|   └───/keren_tmens_analysis
|   |   |   archetype_analysis.ipynb
|   |   |   boundaries_analysis.ipynb
|   |   |   cell_positions_visualization.ipynb
|   |   |   site_radius_analysis.ipynb
|   |   |   sites_patient_mapper.ipynb
|   |   |   TME_architecture_analysis
|
└───/output
|   |
|   └───/csv_files_nature
| 
└───/src
    |	CellAbundance.py
    └───/utils
        |   archetypes.py
        |   equations.py
        |   visualization.py
```

## Reproducing analysis from El Marrahi et al.
* **Niches identification**: Open /TMENS_analysis/notebooks/keren_building_blocks_analysis/archetype_analysis.ipynb and excute "3D Archetypes- All tumors- 4 archetypes", "Visualization and interpretation of archetypes" and "Export data from PCA and Archetype analysis" sections. Open in /TMENS_analysis/keren_building_blocks_analysis sites_patient_mapper.ipynb and execute each chunk and name .csv files
* **Images segmentation into niches**: Open /TMENS_analysis/notebooks/cell_positions_visualization.ipynb
* **Macroscopic analysis of niches from CyTOF data (Wagner et al,2019)**: Open and excute the following R scripts from /macro_niches_analysis folder: 1. Processing of CyTOF data: scBC_analysis.Rmd, 2.Macro-microscopic cell composition of tumors mapping:  scBC_newCells.Rmd, 3. Linear regression of macroscopic cellular abundance over niches: lm_TMENS.Rmd Figures are found in /figs folder from /macro_niches_analysis
* **Niche-phenotype mapping**: Open and excute this file from /phenotypes_niches: niche_phenotype_mapping.Rmd. Figures are found in /phenotypes_niches/figs

## License

## Contact
Fabio Lipreri - <>
Anissa El Marrahi - <anissa.el@scilifelab.se>
Ziqi Kang - <>
Jean Hausser - <jean.hausser@scilifelab.se>

## Acknowledgments
Fabio Lipreri
David Alber
