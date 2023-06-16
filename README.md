# README

NIche Phenotype MAPping (NIPMAP) analysis from spatial multiplex data: Multiplex Ion Imaging on 41 Triple Negative Breast tumors from [Keren et al, Cell(2018)](10.1016/j.cell.2018.08.039) and In Situ Sequencing data on human lung development from [Sountoulidis et al, 2022](https://doi.org/10.1101/2022.01.11.475631)

## Prerequisites

* Jupyter notebook installed

* Python 3.7.13 installed
    List of python libraries to install 

    ```bash
    pip install matplotlib
    pip install scipy
    pip install pandas==1.3.5
    pip install numpy==1.8.1
    pip install scikit-learn
    pip install seaborn
    pip install matplotlib
    pip install qpsolvers==1.9.0
    pip install mpl_toolkits
    ```
* R 4.1.3 with RStudio installed
    List of R libraries to install 
    ```
    pkgs <- c("tidyverse","ggplot2","ade4","factoextra","plotly","igraph","reshape2","ggrepel","viridis","fdrtool","pheatmap","cluster","broom","pROC","ggpubr","devtools","ggridges")
    install.packages(pkgs)
    ```
Installation will take around 20 min in a environment equipped for standard data science.

## Quick start
NIPMAP is a multiplex histology data analysis tool to unravel tissue architecture. It is required prior to starting the analysis to have one .csv file for each Sample image(one image = one patient or sample) (named patient\<Patient ID or number\>_cell_positions.csv) with cells as rows and their data in these columns: (x,y) positions, cell ID inthe image and its cell type respectively named  *x,y,label* and *cell_type*. 
To assess niche-phenotype associations, the input file should be a .csv data table named *cellData.csv* reporting for each cell (rows) their label in the image,the sample of origin and their marker intensity value reported in columns respectively named *cellLabelInImage*, *SampleID*, *marker1*, *marker2*, etc...

1. Open main_nipmap.py and set the parameters. Execute it and enter the size (in micrometers) of the radius of sites, after radius analysis. This script will generate these outputs as json files:
* Generation of sites and cellular abundances within them (+ radius size selection)
* PCA and Archetype Analysis
* Niche identification
* Niche segmentation of images
2. Open py_wrapper_nipmap.r script,enter the parameters and execute it. This script produces these outputs including figures: 
* Niche-phenotype mappping

Run time is around 20 min for one sample and around 60 min for a 40-patient sample MIBI dataset.

Note: in this version of NIPMAP, we consider that one sample = one image
Note #2: NIPMAP doesn't aim to correct cell segmentation error or cell type mis-assignments, this needs to be addressed prior to niche-phenotype mapping. 

## FILES

```
└───/main_nipmap.py: calls py functions for NIPMAP pipeline + saves outputs into .json files
└───/py_wrapper_nipmap.r: get NIPMAP outputs from .json files and loads R objects for further analyses
└───/README.md
|
/ISS_analysis
|
└───/data
|
└───/notebooks
|
└───/output
|
└───/scripts
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
|
└───/clustering_MIBI_data.py: clustering (k-means) of sites cell abundance compared with niche-finding in NIPMAP
└───/building_blocks_mapping.Rmd: Marker expression analysis in niches
└───/data
|   |   cellData.csv: normalized data (markers intensity) from CyTOF expfrom Keren et al.,Cell(2018)
|   |
|   └───/cell_positions_data: .csv files for that contain for each patientx,y positions and label of cells from Keren et al.,Cell(2018)
|   |
|   └───/keren_building_blocks_analysis
|   |   |   archetype_analysis.ipynb
|   |   |   boundaries_analysis.ipynb
|   |   |   cell_positions_visualization.ipynb
|   |   |   site_radius_analysis.ipynb
|   |   |   sites_patient_mapper.ipynb
|
└───/output
| 
└───/src
    |	CellAbundance.py
    └───/utils
        |   archetypes.py
        |   equations.py
        |   visualization.py
```

## Reproducing analysis from El Marrahi et al.
* **Niches identification**: Open TMENS_analysis/notebooks/keren_building_blocks_analysis/sites_radius_analysis.ipynb and execute the chunks to select the radius of sites that allows to explain the most the covariance structure of the dataset. Open TMENS_analysis/notebooks/keren_building_blocks_analysis/archetype_analysis.ipynb and excute "3D Archetypes- All tumors- 4 archetypes", "Visualization and interpretation of archetypes" and "Export data from PCA and Archetype analysis" sections. Open in TMENS_analysis/keren_building_blocks_analysis sites_patient_mapper.ipynb and execute each chunk and name .csv files
* **Images segmentation into niches**: Open TMENS_analysis/notebooks/cell_positions_visualization.ipynb
* **Niche-phenotype mapping**: Open and excute this file from phenotypes_niches: niche_phenotype_mapping.Rmd. Figures are found in /phenotypes_niches/figs
* **Comparing NIPMAP and clustering approaches in niche identification**: run TMENS_analysis/clustering_MIBI_data.py script
* **Macroscopic analysis of niches from CyTOF data (Wagner et al,2019)**: Open and excute the following R scripts from /macro_niches_analysis folder: 1. Processing of CyTOF data: scBC_analysis.Rmd, 2.Macro-microscopic cell composition of tumors mapping:  scBC_newCells.Rmd, 3. Linear regression of macroscopic cellular abundance over niches: lm_TMENS.Rmd Figures are found in /figs folder from /macro_niches_analysis
* **NIPMAP on ISS dataset from Sountoulidis et al**: Open ISS_analysis/notebooks/archetype_analysis.ipynb and excute "Radius Analysis", "Archetype Analysis", "Visualization" and "Save data". The "Save data" part generates all the files needed for downstream analyses. Then, open ISS_analysis/notebooks/HybISS_niche_explore.Rmd and excute the whole Rmarkdown file to generate all the figures included in the section "NIPMAP identifies the cellular and phenotypic architecture of developing lung from in situ RNA sequencing" of the paper. Figures generated for this analysis are found in ISS_analysis/output.

## License

GNU General Public License version 3

## Contact
Anissa El Marrahi - <anissa.el@scilifelab.se>   <anissel12@gmail.com>
Ziqi Kang - <ziqi.kang@scilifelab.se>
Jean Hausser - <jean.hausser@scilifelab.se>

## Acknowledgments
Fabio Lipreri
David Alber
