gc()
rm(list=ls())
#.libPaths("/scratch/anissa.el/R_old/x86_64-redhat-linux-gnu-library/4.0")
.libPaths("/home/common/R")
library(rjson)
library(tidyverse)
library(fdrtool)
library(readxl)
library(survival)
library(data.table)
library(dbscan)
library(progress)
library(ggplot2)
library(dplyr)
library(tidyr)
library(survival)
library(data.table)
library(progress)
library(dbscan)
# library(purrr)
# library(plotly)

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
#!!!!! no space in the names
custom_nichesLabels <- c('Inflammatory', 'Cancer', 'Bfollicle', 'Other', 'Lowdensity')

# Generate all pairwise combinations and combine with '.'
interfaces_names <- combn(custom_nichesLabels, 2, FUN = function(x) paste(x, collapse = "."))
triple_interfaces_names <- combn(custom_nichesLabels, 3, FUN = function(x) paste(x, collapse = "."))

# Define functinal markers
Functional_markers_single <- c("Ki67", "PD1", "IgA/IgG")
Functional_markers <- unlist(lapply(1:length(Functional_markers_single), function(i) combn(Functional_markers_single, i, simplify = FALSE)), recursive = FALSE)

## Long vs short survivors
# Get the patients that do not relapse
long_survivors4000 <- c(24)
# Get the patient PFS
survival_data <- data.frame(
  patient_id = c(9, 24),
  PFS_months = c(4.3035, 32.7858),
  event_status = c(1, 0)
)

# Define treshold to associate each cell to a niche / interface
treshold_niches <- 0.5
treshold_interfaces <- 0.125
treshold_triple_interfaces <- 0.0185



##### VISUALIZE SIMPLEX LONG VS SHORT SURVIVORS #####

## LOAD OUTPUT OBJECTS
# Cell abundance in sites
# sitesCellAb <- as_tibble(lapply(json_data3$cellAbSites,unlist))
# write_csv(sitesCellAb%>%dplyr::select(-c(index, patient_id,site_id)),"sitesCA.csv")
niches <- paste0("a",as.vector(seq(1,NBNICHES,1)))
names(COLNICHES) <- niches
colNiches.hex <-unlist(lapply(COLNICHES, function(x){rgb(x[1],x[2],x[3],maxColorValue = 255)}))

## Archetypes coordinates in reduced PC space
Archs_3D <- do.call(cbind,lapply(json_data2$archs_coord,unlist))
## Projection of sites cell abundance in reduced PC space
pca3D <- matrix(unlist(json_data$PC_proj),nrow=13)[1:4,] #length(CELLTYPES)

## Create a color map assigned to long vs short survivors for the simplex
color <- create_color_map(long_survivors4000, pca3D, NSITES, ImageIDs)
## Plot simplex long vs short survivors with vertices
plot_simplex(pca3D, Archs_3D, color, custom_nichesLabels, colNiches.hex)



##### NICHE ABUNDANCE #####
## Weights for each sites
alfas = json_data2$alfas
## Compute mean alpha (niches and interfaces) for each patient (long and short survivors)
weights_N.I <- Niche_Interfaces_ShortVSLong_survivors(alfas, custom_nichesLabels, interfaces_names, NSITES, long_survivors4000, ImageIDs)
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
colnames(cellsNichesInterfaces)[1:NBNICHES] <- custom_nichesLabels


# # Set the 'LowDensity' column to 0
# cellsNichesInterfaces$LowDensity <- 0
# ## Normalize the other columns ('Inflammatory', 'Cancer', 'Bfollicle', 'Other')
# # Divide each value by the row sum to normalize at 1
# cellsNichesInterfaces[c('Inflammatory', 'Cancer', 'Bfollicle', 'Other')] <- 
#   cellsNichesInterfaces[c('Inflammatory', 'Cancer', 'Bfollicle', 'Other')] / 
#   rowSums(cellsNichesInterfaces[c('Inflammatory', 'Cancer', 'Bfollicle', 'Other')])

# Rename niches and add interfaces weight columns
#cellsNichesInterfaces <- createInterfaces(cellsNichesInterfaces, custom_nichesLabels, interfaces_names, triple_interfaces_names)
cellsNichesInterfaces <- createInterfaces(cellsNichesInterfaces, custom_nichesLabels)
# Associate each cell to a niche or interfaces
cellsNichesInterfaces <- associateCellsToNichesInterfaces(cellsNichesInterfaces, custom_nichesLabels, treshold_niches, interfaces_names, treshold_interfaces, triple_interfaces_names, treshold_triple_interfaces)
# Proportion of niches and interfaces
table(cellsNichesInterfaces$niche)
# Associate cells with niches AND Markers / full phenotype
directory_path <- "./TMENS_analysis/data/cell_positions_data_mIF/"
cells.NichesInterface.Phen <- associateCellsToFunctionalMarkersBI3(cellsNichesInterfaces, directory_path)
colSums(is.na(cells.NichesInterface.Phen))
unique(cells.NichesInterface.Phen$niche)
table(cells.NichesInterface.Phen$niche)



##### VISUALIZATION OF THE CELLS MAPPED BY NICHE #####

## Create colormap for interfaces and triple interfaces
names(COLNICHES) <- custom_nichesLabels
# Function to average RGB and convert to hex
average_rgb_to_hex <- function(colors) {
  rgb_matrix <- do.call(rbind, colors)
  avg_rgb <- colMeans(rgb_matrix)
  rgb(avg_rgb[1], avg_rgb[2], avg_rgb[3], maxColorValue = 255)
}
# Original colors in hex
base_colors <- sapply(COLNICHES, function(rgb_val) {
  rgb(rgb_val[1], rgb_val[2], rgb_val[3], maxColorValue = 255)
})
# Generate pairwise interface names and colors
interfaces_names <- combn(custom_nichesLabels, 2, FUN = function(x) paste(x, collapse = "."))
interface_colors <- sapply(interfaces_names, function(name) {
  labels <- strsplit(name, "\\.")[[1]]
  average_rgb_to_hex(COLNICHES[labels])
}, USE.NAMES = TRUE)
# Generate triple interface names and colors
triple_interfaces_names <- combn(custom_nichesLabels, 3, FUN = function(x) paste(x, collapse = "."))
triple_interface_colors <- sapply(triple_interfaces_names, function(name) {
  labels <- strsplit(name, "\\.")[[1]]
  average_rgb_to_hex(COLNICHES[labels])
}, USE.NAMES = TRUE)
# Add custom 'mixed' color
mixed_color <- c(mixed = rgb(150, 150, 150, maxColorValue = 255))
# Combine all colors
color_vector <- c(base_colors, interface_colors, triple_interface_colors, mixed_color)
# Ensure niche is a factor for correct legend mapping
cells.NichesInterface.Phen$niche <- factor(cells.NichesInterface.Phen$niche, levels = names(color_vector))
unique(cells.NichesInterface.Phen$niche)

# Create the plot with properly assigned colors
p <- ggplot(cells.NichesInterface.Phen, aes(x = x, y = y, color = niche)) +
  geom_point(size = 0.1) +                         
  scale_y_reverse() +                              
  labs(title = "Cell Locations by Niche Across All Samples",
       x = NULL,
       y = NULL,
       color = "Niche") +                           
  scale_color_manual(values = color_vector) +  # Correct color assignment
  guides(color = guide_legend(override.aes = list(shape = 15, size = 7))) +                                            
  facet_wrap(~SampleID, scales = "free") +  # Adaptive scales for each sample
  theme_minimal(base_size = 14) + 
  theme(
    panel.background = element_rect(fill = "white", color = "white"),  # Set white background
    plot.background = element_rect(fill = "white", color = "white"),
    panel.grid.major = element_line(color = "gray90"),  # Light grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.text = element_blank(),        # Remove axis text (tick labels)
    axis.ticks = element_blank(),       # Remove axis ticks
    axis.title = element_blank()        # Remove axis titles
  )
dir.create("comparative_analysis/figs", recursive = TRUE, showWarnings = FALSE)
ggsave("comparative_analysis/figs/Cell_Locations_Niche.png", plot = p, width = 20, height = 10, dpi = 300, bg = "white")


#######---- CONTINUOUS ANALYSES ----#######
##### LINEAGE MARKERS #####
Niches_Interfaces <- unique(cells.NichesInterface.Phen$niche)
cell_types <- unique(cells.NichesInterface.Phen$cell_type)
unique_sample_ids <- unique(cells.NichesInterface.Phen$SampleID)

result <- compute_count_logRatio_and_pvaluesBI3_allcells(cells.NichesInterface.Phen, long_survivors4000, Niches_Interfaces, cell_types, Functional_markers, unique_sample_ids, survival_data)
log_ratio_LS.SS_ <- result$log_ratio_LS.SS
log_ratio_LS.SS <- log_ratio_LS.SS_[ 
  log_ratio_LS.SS_$pvalue <= 1, 
]
survival_data_ratios <- result$survival_data

## If no pseudo count
# Find the max and min values for the hr column (ignoring Inf and -Inf)
max_val <- max(log_ratio_LS.SS$hr[is.finite(log_ratio_LS.SS$hr)], na.rm = TRUE)
min_val <- min(log_ratio_LS.SS$hr[is.finite(log_ratio_LS.SS$hr)], na.rm = TRUE)
# Replace Inf and -Inf with the respective max and min values
log_ratio_LS.SS$hr[log_ratio_LS.SS$hr == Inf] <- max_val
log_ratio_LS.SS$hr[log_ratio_LS.SS$hr == -Inf] <- min_val
# Clamp hr values between 0.5 and 2
log_ratio_LS.SS$hr <- pmax(pmin(log_ratio_LS.SS$hr, 2), 0.5)


##### BUBBLE MAP OF PROPORTIONS
# For the heatmap: Get niche as row, CT+marker as column and hr as values
table_heatmap <- log_ratio_LS.SS %>%
  mutate(cell_type_marker = cell_type) %>%
  select(niche, hr, cell_type_marker, nb_cells, Combination_ID, pvalue)

# Create a reference table with all possible combination of CT+marker with niches
# In order to have niche combination for every CT+marker (also the ones that are not significant -> hr value set to 0)
# Number of niches before and after filtering niches
niches_uniques_for_MFI_L_and_S <- unique(log_ratio_LS.SS$niche)
reference_table <- expand.grid(
  cell_type_marker = unique(table_heatmap$cell_type_marker),
  niche = niches_uniques_for_MFI_L_and_S
)
# Join with table_heatmap (that contains hr value, if the combination is not in table_heatmap -> hr set to 0)
table_heatmap_complete <- reference_table %>%
  left_join(table_heatmap, by = c("cell_type_marker", "niche")) %>%
  mutate(hr = coalesce(hr, 0))  # Replace missing values with 0
table_heatmap_complete <- table_heatmap_complete %>%
  rename(log_ratio_R_NR = hr)
# Ensure pvalue is numeric
table_heatmap_complete$pvalue <- as.numeric(table_heatmap_complete$pvalue)

bubble_map_plot <- table_heatmap_complete %>%
  ggplot(aes(x = cell_type_marker, y = niche)) +
  
  # First layer: all points, no outline
  geom_point(
    aes(color = log_ratio_R_NR, size = 1 / pvalue),
    alpha = 0.7
  ) +
  
  # Second layer: only p < 0.05 points, with black outline
  geom_point(
    data = subset(table_heatmap_complete, pvalue < 0.05),
    aes(color = log_ratio_R_NR, size = 1 / pvalue),
    shape = 21,
    stroke = 0.3,
    fill = NA,
    color = "black"
  ) +
  
  # Third layer: add text labels for median nb_cells when p < 0.05
  geom_text(
    #data = subset(table_heatmap_complete, pvalue < 0.05),
    data = table_heatmap_complete,
    aes(label = nb_cells),
    color = "black",
    size = 2
  ) +
  
  # Color scale for log hazard ratio
  scale_color_gradientn(
    colours = c("#237c04", "white", "red"),
    values = scales::rescale(c(0.5, 1, 2)),
    limits = c(0.5, 2),
    name = "HR"
  ) +
  
  # Size scale for inverse p-value
  scale_size_continuous(
    name = "p-value", 
    breaks = c(1/0.01, 1/0.05, 1/0.1, 1/0.5),
    labels = c("0.01", "0.05", "0.1", "0.5"),
    range = c(2, 10)
  ) +
  
  labs(
    title = "Bubble Map of log ratio vs p-value",
    x = "Cell Type Marker",
    y = "Niche"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Print bubble map
print(bubble_map_plot)
ggsave(
  filename = "comparative_analysis/figs/Bubblemap_ratio_lineage_hazardratio.png",
  plot     = bubble_map_plot,
  dpi      = 300,
  width    = 7,    # adjust width in inches
  height   = 5,    # adjust height in inches
  units    = "in"
)

##### BUBBLE MAP OF SPECIFIC COMBINATION
NI <- "Bfollicle"
CT <- "CD4-T"

NI_table <- cells.NichesInterface.Phen[cells.NichesInterface.Phen$niche == NI, ]
NI_CT_table <- NI_table[NI_table$cell_type == CT, ]
NI_CT_FM_table <- NI_CT_table

## Get number of cell for the specified CT with corresponding FM and in the corresponding NI
NI_CT_FM_table_num_cell <- NI_CT_FM_table %>%
  group_by(SampleID) %>%
  summarise(num_cells = n()) 
# Check if any of the unique SampleIDs are missing in NI_CT_FM_table_num_cell$SampleID
missing_sample_ids <- setdiff(unique_sample_ids, NI_CT_FM_table_num_cell$SampleID)
# If there are missing SampleIDs, create a new data frame with them and num_cells set to 0
if (length(missing_sample_ids) > 0) {
  new_rows <- data.frame(SampleID = missing_sample_ids, num_cells = 0)
  NI_CT_FM_table_num_cell <- rbind(NI_CT_FM_table_num_cell, new_rows)
}

## Get number of total cells in the specified niches
NI_table_num_cell <- NI_table %>%
  group_by(SampleID) %>%
  summarise(num_cells = n()) 
# Check if any of the unique SampleIDs are missing in NI_CT_FM_table_num_cell$SampleID
missing_sample_ids <- setdiff(unique_sample_ids, NI_table_num_cell$SampleID)
# If there are missing SampleIDs, create a new data frame with them and num_cells set to 0
if (length(missing_sample_ids) > 0) {
  new_rows <- data.frame(SampleID = missing_sample_ids, num_cells = 0)
  NI_table_num_cell <- rbind(NI_table_num_cell, new_rows)
}



count_table <- NI_table_num_cell %>%
  rename(num_cells_NI = num_cells) %>% 
  left_join(NI_CT_FM_table_num_cell %>% rename(num_cells_NI_CT_FM = num_cells), by = "SampleID")
# Remove patients where number of cells in the niche is less than 100
# count_table <- count_table %>%
#   filter(num_cells_NI >= 100)

count_table$ratio_num_cells <- count_table$num_cells_NI_CT_FM / count_table$num_cells_NI

LS_count_table <- subset(count_table, SampleID %in% long_survivors4000)
SS_count_table <- subset(count_table, !(SampleID %in% long_survivors4000))


LS_filtered <- na.omit(LS_count_table$ratio_num_cells)
SS_filtered <- na.omit(SS_count_table$ratio_num_cells)
length(LS_filtered)
length(SS_filtered)

# Only run Wilcoxon test if both groups have at least two observations
if (length(LS_filtered) >= 1 && length(SS_filtered) >= 1) {
  p_value <- wilcox.test(LS_filtered, SS_filtered, exact = FALSE)$p.value
} else {
  p_value <- "No niche for one group of patients"
}
LS_count_table <- subset(count_table, SampleID %in% long_survivors4000)
SS_count_table <- subset(count_table, !(SampleID %in% long_survivors4000))
LS_ratio <- LS_count_table$ratio_num_cells
SS_ratio <- SS_count_table$ratio_num_cells
# Create a combined data frame
data_combined <- data.frame(
  Ratio = c(LS_ratio, SS_ratio),
  Group = c(rep("Non-relapsed", length(LS_ratio)), rep("Relapsed", length(SS_ratio)))
)

proportion_violin <- ggplot(data_combined, aes(x = Group, y = Ratio, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.5) + 
  geom_jitter(color = "black", alpha = 0.5, width = 0.2) +
  labs(title = "Comparison of LS and SS Ratio Distributions",
       x = "Group",
       y = "Proportion of cells") +
  scale_fill_manual(values = c("#237c04", "red")) +
  theme_minimal()
print(proportion_violin)
ggsave(
  filename = sprintf("comparative_analysis/figs/Proportion_%s_in_%s.png", CT, NI),
  plot     = proportion_violin,
  dpi      = 300,
  width    = 4,    # adjust width in inches
  height   = 4,    # adjust height in inches
  units    = "in"
)



##### FUNCTIONAL MARKERS #####
result <- compute_count_logRatio_and_pvaluesBI3_final(cells.NichesInterface.Phen, long_survivors4000, Niches_Interfaces, cell_types, Functional_markers, unique_sample_ids, survival_data)
log_ratio_LS.SS_ <- result$log_ratio_LS.SS
log_ratio_LS.SS <- log_ratio_LS.SS_[ 
  log_ratio_LS.SS_$pvalue <= 1, 
]
survival_data_ratios <- merge(survival_data_ratios, result$survival_data, by = c("patient_id", "PFS_months", "event_status"), all.x = TRUE)
log_ratio_LS.SS <- na.omit(log_ratio_LS.SS)

## If no pseudo count
# Find the max and min values for the hr column (ignoring Inf and -Inf)
max_val <- max(log_ratio_LS.SS$hr[is.finite(log_ratio_LS.SS$hr)], na.rm = TRUE)
min_val <- min(log_ratio_LS.SS$hr[is.finite(log_ratio_LS.SS$hr)], na.rm = TRUE)
# Replace Inf and -Inf with the respective max and min values
log_ratio_LS.SS$hr[log_ratio_LS.SS$hr == Inf] <- max_val
log_ratio_LS.SS$hr[log_ratio_LS.SS$hr == -Inf] <- min_val
# Clamp hr values between 0.5 and 2
log_ratio_LS.SS$hr <- pmax(pmin(log_ratio_LS.SS$hr, 2), 0.5)

## Remove Tumor cells with IgG/IgA and Other with IgGIgA
log_ratio_LS.SS <- log_ratio_LS.SS %>%
  filter(
    !(cell_type %in% c("B cell", "Keratin-positive tumor", "Other") & grepl("IgA/IgG", marker))
  )


##### BUBBLE MAP OF PROPORTIONS
# For the heatmap: Get niche as row, CT+marker as column and hr as values
table_heatmap <- log_ratio_LS.SS %>%
  mutate(cell_type_marker = paste(cell_type, marker, sep = "_")) %>%
  select(niche, hr, cell_type_marker, nb_cells, Combination_ID, pvalue)

# Create a reference table with all possible combination of CT+marker with niches
# In order to have niche combination for every CT+marker (also the ones that are not significant -> logratio value set to 0)
# Number of niches before and after filtering niches
niches_uniques_for_MFI_L_and_S <- unique(log_ratio_LS.SS$niche)
reference_table <- expand.grid(
  cell_type_marker = unique(table_heatmap$cell_type_marker),
  niche = niches_uniques_for_MFI_L_and_S
)
# Join with table_heatmap (that contains logratio value, if the combination is not in table_heatmap -> logratio set to 0)
table_heatmap_complete <- reference_table %>%
  left_join(table_heatmap, by = c("cell_type_marker", "niche")) %>%
  mutate(hr = coalesce(hr, 0))  # Replace missing values with 0
table_heatmap_complete <- table_heatmap_complete %>%
  rename(log_ratio_R_NR = hr)
# Ensure pvalue is numeric
table_heatmap_complete$pvalue <- as.numeric(table_heatmap_complete$pvalue)


bubble_map_plot <- table_heatmap_complete %>%
  ggplot(aes(x = cell_type_marker, y = niche)) +
  
  # First layer: all points, no outline
  geom_point(
    aes(color = log_ratio_R_NR, size = 1 / pvalue),
    alpha = 0.7
  ) +
  
  # Second layer: only p < 0.05 points, with black outline
  geom_point(
    data = subset(table_heatmap_complete, pvalue < 0.05),
    aes(color = log_ratio_R_NR, size = 1 / pvalue),
    shape = 21,
    stroke = 0.3,
    fill = NA,
    color = "black"
  ) +
  
  # Third layer: add text labels for median nb_cells when p < 0.05
  geom_text(
    #data = subset(table_heatmap_complete, pvalue < 0.05),
    data = table_heatmap_complete,
    aes(label = nb_cells),
    color = "black",
    size = 2
  ) +
  
  # Color scale for hazard ratio
  scale_color_gradientn(
    colours = c("#237c04", "white", "red"),
    values = scales::rescale(c(0.5, 1, 2)),
    limits = c(0.5, 2),
    name = "HR"
  ) +
  
  # Size scale for inverse p-value
  scale_size_continuous(
    name = "p-value", 
    breaks = c(1/0.01, 1/0.05, 1/0.1, 1/0.5),
    labels = c("0.01", "0.05", "0.1", "0.5"),
    range = c(2, 10)
  ) +
  
  labs(
    title = "Bubble Map of log ratio vs p-value",
    x = "Cell Type Marker",
    y = "Niche"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Print bubble map
print(bubble_map_plot)
ggsave(
  filename = "comparative_analysis/figs/Bubblemap_ratio_functional_hazardratio.png",
  plot     = bubble_map_plot,
  dpi      = 300,
  width    = 7,   
  height   = 5,    
  units    = "in",
  limitsize = FALSE
)


##### BUBBLE MAP OF SPECIFIC COMBINATION
NI <- "Bfollicle"
NI_table <- cells.NichesInterface.Phen[cells.NichesInterface.Phen$niche == NI, ]
CT <- "CD4-T"
NI_CT_table <- NI_table[NI_table$cell_type == CT, ]
FM <- c("PD1")

# Remove cases where the cell type is Plasma cells and functional marker "IgA/IgG"
if (CT == "Plasma cell") {
  # Case 1: If the list FM contains only "IgA/IgG", skip to next iteration
  if (length(FM) == 1 && FM == "IgA/IgG") {
    next 
  }
  
  # Case 2: If the list FM contains "IgA/IgG", remove it from the list
  if ("IgA/IgG" %in% FM) {
    FM <- FM[FM != "IgA/IgG"] 
  }
}
match_matrix <- do.call(cbind, lapply(FM, function(fm) str_detect(NI_CT_table$Phenotype, fm)))
match_matrix <- as.matrix(match_matrix)
match_matrix[is.na(match_matrix)] <- FALSE
NI_CT_FM_table <- NI_CT_table[rowSums(match_matrix) == length(FM), ]
## Get number of cell for the specified CT with corresponding FM and in the corresponding NI
NI_CT_FM_table_num_cell <- NI_CT_FM_table %>%
  group_by(SampleID) %>%
  summarise(num_cells = n()) 
# Check if any of the unique SampleIDs are missing in NI_CT_FM_table_num_cell$SampleID
missing_sample_ids <- setdiff(unique_sample_ids, NI_CT_FM_table_num_cell$SampleID)
# If there are missing SampleIDs, create a new data frame with them and num_cells set to 0
if (length(missing_sample_ids) > 0) {
  new_rows <- data.frame(SampleID = missing_sample_ids, num_cells = 0)
  NI_CT_FM_table_num_cell <- rbind(NI_CT_FM_table_num_cell, new_rows)
}
## Get number of total cells in the specified niches
NI_table_num_cell <- NI_table %>%
  group_by(SampleID) %>%
  summarise(num_cells = n()) 
# Check if any of the unique SampleIDs are missing in NI_CT_FM_table_num_cell$SampleID
missing_sample_ids <- setdiff(unique_sample_ids, NI_table_num_cell$SampleID)

# Remove patients where there is no cells in the niche
if (length(missing_sample_ids) > 0) {
  NI_CT_FM_table_num_cell <- NI_CT_FM_table_num_cell %>%
    filter(!SampleID %in% missing_sample_ids)
}

count_table <- NI_table_num_cell %>%
  rename(num_cells_NI = num_cells) %>% 
  left_join(NI_CT_FM_table_num_cell %>% rename(num_cells_NI_CT_FM = num_cells), by = "SampleID")
# Remove patients where number of cells in the niche is less than 100
# count_table <- count_table %>%
#   filter(num_cells_NI >= 100)

count_table$ratio_num_cells <- count_table$num_cells_NI_CT_FM / count_table$num_cells_NI

LS_count_table <- subset(count_table, SampleID %in% long_survivors4000)
SS_count_table <- subset(count_table, !(SampleID %in% long_survivors4000))

LS_filtered <- na.omit(LS_count_table$ratio_num_cells)
SS_filtered <- na.omit(SS_count_table$ratio_num_cells)
length(LS_filtered)
length(SS_filtered)

# Only run Wilcoxon test if both groups have at least two observations
if (length(LS_filtered) >= 1 && length(SS_filtered) >= 1) {
  p_value <- wilcox.test(LS_filtered, SS_filtered, exact = FALSE)$p.value
} else {
  p_value <- "No niche for one group of patients"
}
print(p_value)

LS_ratio <- LS_count_table$ratio_num_cells
SS_ratio <- SS_count_table$ratio_num_cells
data_combined <- data.frame(
  Ratio = c(LS_ratio, SS_ratio),
  Group = c(rep("Non-relapsed", length(LS_ratio)), rep("Relapsed", length(SS_ratio)))
)

proportion_violin <- ggplot(data_combined, aes(x = Group, y = Ratio, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.5) + 
  geom_jitter(color = "black", alpha = 0.5, width = 0.2) +
  labs(title = "Comparison of LS and SS Ratio Distributions",
       x = "Group",
       y = "Proportion of cells") +
  scale_fill_manual(values = c("#237c04", "red")) +
  theme_minimal()
print(proportion_violin)
FM_clean <- gsub(",", "_", gsub("/", "", FM))
ggsave(
  filename = sprintf("comparative_analysis/figs/Proportion_%s_%s_in_%s.png", CT, FM_clean, NI),
  plot     = proportion_violin,
  dpi      = 300,
  width    = 4,    # adjust width in inches
  height   = 4,    # adjust height in inches
  units    = "in"
)




##### INTERACTIONS #####

combinations <- combn(CELLTYPES, 2)
possible_interactions <- as.data.frame(t(combinations))
colnames(possible_interactions) <- c("Reference", "Target")
# 1) Define the two axes of interest:
refs <- c("CD4-T","CD8-T","Treg",
          "B cell","Plasma cell",
          "Keratin-positive tumor")
tgts <- c("CD4-T","CD8-T","Treg",
          "cDC1","cDC2","mature DC","pDC",
          "B cell","Folicular DC")
# 2) Build the directed interaction list:
possible_interactions <- expand.grid(
  Reference = refs,
  Target    = tgts,
  stringsAsFactors = FALSE
)

result <- compute_interactions_and_pval(cells.NichesInterface.Phen, Niches_Interfaces, possible_interactions, survival_data)
cox_model <- result[[1]]
survival_data_ratios <- merge(survival_data_ratios, result[[2]], by = c("patient_id", "PFS_months", "event_status"), all.x = TRUE)

log_ratio_LS.SS_ <- cox_model
log_ratio_LS.SS <- log_ratio_LS.SS_[ 
  log_ratio_LS.SS_$p.value <= 1, 
]
log_ratio_LS.SS <- na.omit(log_ratio_LS.SS)
  
## If no pseudo count
# Find the max and min values for the hr column (ignoring Inf and -Inf)
max_val <- max(log_ratio_LS.SS$hr[is.finite(log_ratio_LS.SS$hr)], na.rm = TRUE)
min_val <- min(log_ratio_LS.SS$hr[is.finite(log_ratio_LS.SS$hr)], na.rm = TRUE)
# Replace Inf and -Inf with the respective max and min values
log_ratio_LS.SS$hr[log_ratio_LS.SS$hr == Inf] <- max_val
log_ratio_LS.SS$hr[log_ratio_LS.SS$hr == -Inf] <- min_val
# Clamp hr values between 0.5 and 2
log_ratio_LS.SS$hr <- pmax(pmin(log_ratio_LS.SS$hr, 2), 0.5)


log_ratio_LS.SS <- log_ratio_LS.SS %>%
  rename(interaction_NI = interaction) %>%  # rename column
  mutate(
    # Remove "Interactions_" prefix first
    temp = str_remove(interaction_NI, "^Interactions_"),
    # Now split into interaction and niche
    interaction = str_extract(temp, ".*(?=_[^_]+$)"),  # everything up to the last "_"
    niche = str_extract(temp, "[^_]+$")                # everything after the last "_"
  ) %>%
  select(-temp)  


##### BUBBLE MAP OF INTERACTIONS
# For the heatmap: Get niche as row, CT+marker as column and log_ratioLS as values
table_heatmap <- log_ratio_LS.SS %>%
  mutate(cell_type_marker = interaction) %>%
  mutate(pvalue = p.value) %>%
  select(hr, cell_type_marker, pvalue, niche, median_nb_interactions)
table_heatmap_complete <- table_heatmap %>%
  rename(log_ratio_R_NR = hr)
# Ensure pvalue is numeric
table_heatmap_complete$pvalue <- as.numeric(table_heatmap_complete$pvalue)

## BUBBLE MAP OF T/B cells INTERACTIONS
filtered_table_heatmap <- table_heatmap_complete %>%
  filter(str_extract(cell_type_marker, "^[^_]+") %in% c("CD4-T", "CD8-T", "Treg", "B cell"))

bubble_map_plot <- filtered_table_heatmap %>%
  ggplot(aes(x = cell_type_marker, y = niche)) +
  
  # First layer: all points, no outline
  geom_point(
    aes(color = log_ratio_R_NR, size = 1 / pvalue),
    alpha = 0.7
  ) +
  
  # Second layer: only p < 0.05 points, with black outline
  geom_point(
    data = subset(filtered_table_heatmap, pvalue < 0.05),
    aes(color = log_ratio_R_NR, size = 1 / pvalue),
    shape = 21,
    stroke = 0.3,
    fill = NA,
    color = "black"
  ) +
  
  # Third layer: add text labels for median nb_cells when p < 0.05
  geom_text(
    data = subset(filtered_table_heatmap, pvalue < 0.05),
    #data = filtered_table_heatmap,
    aes(label = median_nb_interactions),
    color = "black",
    size = 2
  ) +
  
  # Color scale for hazard ratio
  scale_color_gradientn(
    colours = c("#237c04", "white", "red"),
    values = scales::rescale(c(0.5, 1, 2)),
    limits = c(0.5, 2),
    name = "HR"
  ) +
  
  # Size scale for inverse p-value
  scale_size_continuous(
    name = "p-value", 
    breaks = c(1/0.01, 1/0.05, 1/0.1, 1/0.5),
    labels = c("0.01", "0.05", "0.1", "0.5"),
    range = c(2, 10)
  ) +
  
  labs(
    title = "Bubble Map of log ratio vs p-value",
    x = "Cell Type Marker",
    y = "Niche"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Print bubble map
print(bubble_map_plot)
ggsave(
  filename = "comparative_analysis/figs/Bubblemap_TB_interactions_hazardratio.png",
  plot     = bubble_map_plot,
  dpi      = 300,
  width    = 8,   
  height   = 5,    
  units    = "in",
  limitsize = FALSE
)


## BUBBLE MAP OF TUMOR cells INTERACTIONS
filtered_table_heatmap <- table_heatmap_complete %>%
  filter(str_extract(cell_type_marker, "^[^_]+") %in% c("Keratin-positive tumor"))

bubble_map_plot <- filtered_table_heatmap %>%
  ggplot(aes(x = cell_type_marker, y = niche)) +
  
  # First layer: all points, no outline
  geom_point(
    aes(color = log_ratio_R_NR, size = 1 / pvalue),
    alpha = 0.7
  ) +
  
  # Second layer: only p < 0.05 points, with black outline
  geom_point(
    data = subset(filtered_table_heatmap, pvalue < 0.05),
    aes(color = log_ratio_R_NR, size = 1 / pvalue),
    shape = 21,
    stroke = 0.3,
    fill = NA,
    color = "black"
  ) +
  
  # Third layer: add text labels for median nb_cells when p < 0.05
  geom_text(
    data = subset(filtered_table_heatmap, pvalue < 0.05),
    #data = filtered_table_heatmap,
    aes(label = round(median_nb_interactions, 1)),
    color = "black",
    size = 2
  ) +
  
  # Color scale for hazard ratio
  scale_color_gradientn(
    colours = c("#237c04", "white", "red"),
    values = scales::rescale(c(0.5, 1, 2)),
    limits = c(0.5, 2),
    name = "HR"
  ) +
  
  # Size scale for inverse p-value
  scale_size_continuous(
    name = "p-value", 
    breaks = c(1/0.01, 1/0.05, 1/0.1, 1/0.5),
    labels = c("0.01", "0.05", "0.1", "0.5"),
    range = c(2, 10)
  ) +
  
  labs(
    title = "Bubble Map of log ratio vs p-value",
    x = "Cell Type Marker",
    y = "Niche"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Print bubble map
print(bubble_map_plot)
ggsave(
  filename = "comparative_analysis/figs/Bubblemap_tumorcells_interactions_hazardratio.png",
  plot     = bubble_map_plot,
  dpi      = 300,
  width    = 6,   
  height   = 5,    
  units    = "in",
  limitsize = FALSE
)



##### BUBBLE MAP OF SPECIFIC INTERACTIONS
NI = 'Inflammatory'
In = 'Keratin-positive tumor_CD4-T'

interactions_niche <- data.frame(
  patient_id = survival_data_ratios$patient_id,
  Interactions = survival_data_ratios[[sprintf("Interactions_%s_%s", In, NI)]]
)

LS_count_table <- subset(interactions_niche, patient_id %in% long_survivors4000)
SS_count_table <- subset(interactions_niche, !(patient_id %in% long_survivors4000))

LS_filtered <- na.omit(LS_count_table$Interactions)
SS_filtered <- na.omit(SS_count_table$Interactions)

length(LS_filtered)
length(SS_filtered)

# Only run Wilcoxon test if both groups have at least two observations
if (length(LS_filtered) >= 1 && length(SS_filtered) >= 1) {
  p_value <- wilcox.test(LS_filtered, SS_filtered, exact = FALSE)$p.value
} else {
  p_value <- "No niche for one group of patients"
}

print(p_value)

LS_ratio <- LS_count_table$Interactions
SS_ratio <- SS_count_table$Interactions

data_combined <- data.frame(
  Ratio = c(LS_ratio, SS_ratio),
  Group = c(rep("Non-relapsed", length(LS_ratio)), rep("Relapsed", length(SS_ratio)))
)

proportion_violin <- ggplot(data_combined, aes(x = Group, y = Ratio, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.5) + 
  geom_jitter(color = "black", alpha = 0.5, width = 0.2) + 
  labs(title = "Comparison of LS and SS Ratio Distributions",
       x = "Group",
       y = "Normalized Interactions Score") +
  scale_fill_manual(values = c("#237c04", "red")) +
  theme_minimal()
print(proportion_violin)
ggsave(
  filename = sprintf("comparative_analysis/figs/Interactions_%s_in_%s.png", In, NI),
  plot     = proportion_violin,
  dpi      = 300,
  width    = 4,    # adjust width in inches
  height   = 4,    # adjust height in inches
  units    = "in"
)
