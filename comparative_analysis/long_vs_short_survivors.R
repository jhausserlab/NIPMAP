gc()
rm(list=ls())
#.libPaths("/scratch/anissa.el/R_old/x86_64-redhat-linux-gnu-library/4.0")
.libPaths("/home/common/R")
library(rjson)
library(tidyverse)
library(fdrtool)
# library(purrr)
# library(plotly)
# library(tidyr)
# library(dplyr)

### SET WORKING DIRECTORY
dirName <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
setwd(dirName)
# Get a list of all R files in the specified directory
r_files <- list.files(path = "./comparative_analysis/src/", pattern = "\\.R$", full.names = TRUE)
# Source each R file using a for loop
for (file in r_files) {
  source(file)
}



jsonparams <- fromJSON(file="./params.json")
CELLTYPES <-jsonparams$cellTypes
ImageIDs <- jsonparams$ImageID
NSITES <- jsonparams$nbsites
RADIUS <- jsonparams$radiusSize
NBNICHES <- jsonparams$nbniches
METHOD <-jsonparams$countMeth
W <-jsonparams$xsize
H <-jsonparams$ysize
ROOT_DATA_PATH <- jsonparams$rootDataPath
ROOT_OUTPUT_PATH <-jsonparams$rootOutPath
COLNICHES <- jsonparams$colNiches
pathFigs <- jsonparams$pathFigs


file1 = "./pca_sites.json" # pca object on sites elements
file2 = "./AA_sites.json" # archetype Analysis object based on sites cell abundance
file3 = "./ca_sites.json" # cell abundance of randomly generated sites
file4 = "./cells_niches.json" # sites centered on cells and niches weights

#######---- Open .json files ----#######
json_data <- fromJSON(file=file1)
json_data2 <- fromJSON(file=file2)
json_data3 <- fromJSON(file=file3)
json_data4 <- fromJSON(file=file4)



#######---- Comparative analysis Short vs long survivors ----#######


### SET VARIABLES
## Niche indentification
# Define a niche for each weight from barplotNiches.pdf
custom_nichesLabels <- c("TLS", 'inflammatory', 'cancer', 'necrotic')
# Shorten interfaces name for barplot visibility
short_interfaces_names <- c('TLS.inflam', 'TLS.cancer', 'TLS.necr', 'inflam.cancer', 'inflam.necr', 'cancer.necr')
# markers that we don't want to use for the comparative analysis on MFI values (lineage + B7H3, OX40, CD163, `CSF-1R`)
Unwanted_markers <- c("CD11b", "CD11c", "CD16", "CD20", "CD209", "CD3", "CD31", "CD4", "CD45", "CD56",
                     "CD68", "CD8", "dsDNA", "EGFR", "MPO", "Pan-Keratin", "SMA", "Vimentin",
                     "B7H3", "OX40", "CD163", "CSF-1R")

## Long vs short survivors
# Get the patients indices that have over 4000 days of survivals into a list
long_survivors4000 <- c(12, 14, 18, 20, 25, 26)
# Define treshold to associate each cell to a niche / interface
treshold_niches <- 0.5
treshold_interfaces <- 0.125



### VISUALIZE SIMPLEX LONG VS SHORT SURVIVORS ###

## LOAD OUTPUT OBJECTS
# Cell abundance in sites
sitesCellAb <- as_tibble(lapply(json_data3$cellAbSites,unlist))
write_csv(sitesCellAb%>%dplyr::select(-c(index, patient_id,site_id)),"sitesCA.csv")
niches <- paste0("a",as.vector(seq(1,NBNICHES,1)))
names(COLNICHES) <- niches
colNiches.hex <-unlist(lapply(COLNICHES, function(x){rgb(x[1],x[2],x[3],maxColorValue = 255)}))

## Archetypes coordinates in reduced PC space
Archs_3D <- do.call(cbind,lapply(json_data2$archs_coord,unlist))
## Projection of sites cell abundance in reduced PC space
pca3D <- matrix(unlist(json_data$PC_proj),nrow=17)[1:3,] #length(CELLTYPES)

## Create a color map assigned to long vs short survivors for the simplex
color <- create_color_map(long_survivors4000, pca3D)
## Plot simplex long vs short survivors with vertices
plot_simplex(pca3D, Archs_3D, color, custom_nichesLabels, colNiches.hex)



### NICHE ABUNDANCE ###

## Weights for each sites
alfas = json_data2$alfas
## Compute mean alpha (niches and interfaces) for each patient (long and short survivors)
weights_N.I <- Niche_Interfaces_ShortVSLong_survivors(alfas, custom_nichesLabels, short_interfaces_names, NSITES, long_survivors4000, ImageIDs)
long_survivors <- weights_N.I$long_survivors
short_survivors <- weights_N.I$short_survivors

## NICHES barplot
# Select niche and 'patientIDs' columns
long_survivors_niches <- long_survivors[, c(1:NBNICHES, max(ncol(long_survivors))-1)]
short_survivors_niches <- short_survivors[, c(1:NBNICHES, max(ncol(short_survivors))-1)]
# Plot barplot niche abundance
barplot_N.I_abundance(long_survivors_niches, short_survivors_niches, "niche")

## INTERFACES barplot
# Select interfaces and 'patientIDs' columns
long_survivors_interfaces <- long_survivors[, -c(1:NBNICHES, max(ncol(long_survivors)))]
short_survivors_interfaces <- short_survivors[, -c(1:NBNICHES, max(ncol(short_survivors)))]
# Plot barplot interface abundance
barplot_N.I_abundance(long_survivors_interfaces, short_survivors_interfaces, "interface")



### NICHE-PHENOTYPE MAPPING LONG VS SHORT SURVIVORS

# Niches weights(proportions) of all cells from all images
niches<- paste0("a",as.vector(seq(1,NBNICHES,1)))
cellsNichesInterfaces <- as_tibble(lapply(json_data4$cells_niches,unlist))%>%
  rename_at(vars(matches("[0-9]")),~niches)%>%
  mutate(cell_id=as.numeric(cell_id))%>%
  select(-TOT_cell_dens)
# Rename niches and add interfaces weight columns
cellsNichesInterfaces <- createInterfaces(cellsNichesInterfaces, custom_nichesLabels, short_interfaces_names)
# Associate each cell to a niche or interfaces
cellsNichesInterfaces <- associateCellsToNichesInterfaces(cellsNichesInterfaces, custom_nichesLabels, treshold_niches, short_interfaces_names, treshold_interfaces)
# Associate cells with niches AND MFI of functional markers
cells.NichesInterface.Phen <- associateCellsToFunctionalMarkers(cellsNichesInterfaces, Unwanted_markers)

## Continuous analyses
Niches_Interfaces <- unique(cells.NichesInterface.Phen$niche)
cell_types <- unique(cells.NichesInterface.Phen$cell_type)
Functionnal_markers <- read.csv("./phenotypes_niches/data/proteins_by_frame.csv")%>%filter(Purpose=="Functional")%>%pull(Biomarker)
# Define min MFI value  
min_MFI <- compute_minMFI(cells.NichesInterface.Phen, long_survivors4000, Niches_Interfaces, cell_types, Functionnal_markers)
# Create table with p values comparing long and short survivors MFI values of functionnal markers for niche and cell type
log_ratio_LS.SS <- compute_logRatio_and_pvalues(cells.NichesInterface.Phen, long_survivors4000, min_MFI, Niches_Interfaces, cell_types, Functionnal_markers)
# Remove NaN p values (same MFI values for all cells in both long and short survivors)
nan_pvalues <- which(is.nan(log_ratio_LS.SS$pvalue))
log_ratio_LS.SS_ <- log_ratio_LS.SS[complete.cases(log_ratio_LS.SS$pvalue), ]

# Number of niches before and after filtering niches that contain at least 100 cells in long or short survivors
niches_uniques_for_MFI_L_and_S <- unique(log_ratio_LS.SS$niche)
plot_nicheCount(Niches_Interfaces, niches_uniques_for_MFI_L_and_S)

# Correction for multiple testing
# Compute Q values and remove combination for which it is greater than 0.1
FDR <- fdrtool(log_ratio_LS.SS_$pvalue, statistic="pvalue")
qvalues <- FDR$qval
log_ratio_LS.SS_$qvalue <- qvalues
log_ratio_LS.SS_qvalTRESH <- log_ratio_LS.SS_[log_ratio_LS.SS_$qvalue <= 0.1, ]

# If there is still too much significant combination 
# Select log ratio that are greater or lower than logRatio treshold
tresh_logRatio <- log10(1.3)
log_ratio_LS.SS_qvalTRESH_logratioTRESH <- log_ratio_LS.SS_qvalTRESH[log_ratio_LS.SS_qvalTRESH$log_ratioLS > tresh_logRatio 
                                                                     | log_ratio_LS.SS_qvalTRESH$log_ratioLS < -tresh_logRatio, ]

# Heatmap du log ratio long vs short for significant q values
heatmap_logRatio_LvsS_significantQvalues(log_ratio_LS.SS_qvalTRESH_logratioTRESH, log_ratio_LS.SS_qvalTRESH, niches_uniques_for_MFI_L_and_S)

# plot distribution of MFI values for a specific combination marker-celltype-niche
# You can find the Comb_ID corresponding to a heatmap value / combination marker-celltype-niche in table_heatmap_complete$Combination_ID
Comb_ID <- 207
plot_MFI_distribution(Comb_ID, cells.NichesInterface.Phen, long_survivors4000, min_MFI, Niches_Interfaces, cell_types, Functionnal_markers)
