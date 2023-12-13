gc()
rm(list=ls())
#.libPaths("/scratch/anissa.el/R_old/x86_64-redhat-linux-gnu-library/4.0")
.libPaths("/home/common/R")
library(rjson)
library(tidyverse)
library(purrr)
library(plotly)
library(tidyr)
library(ggpubr)
library(dplyr)
library(fdrtool)

### SET WORKING DIRECTORY
dirName <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
setwd(dirName)
source("./phenotypes_niches/functions_phenotypes_tmens.r")

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
custom_nichesLabels <- list("TLS", 'inflammatory', 'cancer', 'necrotic')
# Shorten interfaces name for barplot visibility
short_interfaces_names <- c('TLS.inflam', 'TLS.cancer', 'TLS.necr', 'inflam.cancer', 'inflam.necr', 'cancer.necr')

## Long vs short survivors
# Get the patients indices that have over 4000 days of survivals into a list
long_survivors4000 <- c(12, 14, 18, 20, 25, 26)
# Define treshold to associate each cell to a niche / interface
treshold_niches <- 0.5
treshold_interfaces <- 0.125



### VISUALIZE SIMPLEX LONG VS SHORT SURVIVORS

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
site_groups <- rep(1, length(pca3D[1,]))
value_to_set <- 2
for (i in seq_along(long_survivors4000)) {
  start_index <- 100 * (long_survivors4000[i] - 1) + 1
  end_index <- 100 * long_survivors4000[i]
  site_groups[start_index:end_index] <- value_to_set
}
map_site_group_to_color <- function(site_group) {
  if (site_group == 1) {
    return("rgba(200, 0, 0, 0.2)")
  } else if (site_group == 2) {
    return("green")
  }
}
color <- sapply(site_groups, map_site_group_to_color)

## Plot simplex
# Define simplex vertices coordinates
simplex_x <- c(Archs_3D[1, 1], Archs_3D[2, 1], Archs_3D[3, 1], Archs_3D[4, 1], Archs_3D[1, 1])
simplex_y <- c(Archs_3D[1, 2], Archs_3D[2, 2], Archs_3D[3, 2], Archs_3D[4, 2], Archs_3D[1, 2])
simplex_z <- c(Archs_3D[1, 3], Archs_3D[2, 3], Archs_3D[3, 3], Archs_3D[4, 3], Archs_3D[1, 3])
simplex_x <- c(simplex_x, Archs_3D[c(1, 3), 1], Archs_3D[c(2, 4), 1])
simplex_y <- c(simplex_y, Archs_3D[c(1, 3), 2], Archs_3D[c(2, 4), 2])
simplex_z <- c(simplex_z, Archs_3D[c(1, 3), 3], Archs_3D[c(2, 4), 3])

plotly::plot_ly(x = pca3D[1, ],
                y = pca3D[2, ],
                z = pca3D[3, ],
                type = "scatter3d",
                mode = "markers",
                marker = list(symbol = "triangle", size = 4, color = color),
                name = "sites",
                mode = "text") %>%
  add_trace(x = Archs_3D[, 1],
            y = Archs_3D[, 2],
            z = Archs_3D[, 3],
            type = "scatter3d",
            mode = "markers+text",
            text = custom_nichesLabels,
            textposition = c('top right', 'bottom right', 'top left', 'top right'),
            textfont = list(color = '#000000', size = 16),
            showlegend = TRUE,
            name = "niches",
            marker = list(color = ~colNiches.hex, symbol = "star-diamond", size = 12),
            inherit = FALSE) %>%
  add_trace(x = simplex_x,  # Use the defined simplex coordinates
            y = simplex_y,  # Use the defined simplex coordinates
            z = simplex_z,  # Use the defined simplex coordinates
            type = "scatter3d",
            mode = "markers+lines",  # Draw lines
            name = "Simplex",
            line = list(color = "blue", width = 2)) %>%
  layout(scene = list(xaxis = list(title = "PC1"),
                      yaxis = list(title = "PC2"),
                      zaxis = list(title = "PC3"))
  )



### NICHE ABUNDANCE
## Compute mean alpha (niches and interfaces) for each patient
alfas = json_data2$alfas
# Define the column names
col_names <- unlist(custom_nichesLabels)
# Convert the list to a data frame with specified column names
alfas_NI <- data.frame(
  `colnames<-`(do.call(cbind, alfas), col_names)
)

# Get all combinations of column names
column_combinations <- combn(col_names, 2, simplify = FALSE)
# Iterate through the column combinations
for (cols in column_combinations) {
  col1 <- cols[1]
  col2 <- cols[2]
  # Create a new column name
  new_col_name <- paste(col1, col2, sep = ".")
  # Create the new column by multiplying the two existing columns
  alfas_NI[new_col_name] <- alfas_NI[[col1]] * alfas_NI[[col2]]
}

# Rename interfaces with shorten names
names(alfas_NI)[(ncol(alfas_NI) - (choose(NBNICHES, 2)-1)):(ncol(alfas_NI))] <- short_interfaces_names
# Number of rows in each chunk
chunk_size <- NSITES
# Create a data frame to store the means or medians
alfas_NI_patients <- data.frame(matrix(NA, ncol = ncol(alfas_NI), nrow = 0))
# Iterate through the original data frame in chunks
for (i in seq(1, nrow(alfas_NI), by = chunk_size)) {
  chunk <- alfas_NI[i:min(i + chunk_size - 1, nrow(alfas_NI)), ]
  # Calculate the mean for each column in the chunk
  chunk_mean <- colMeans(chunk, na.rm = TRUE)
  # Alternatively, you can calculate the median using colMedians from the 'matrixStats' package
  # chunk_median <- matrixStats::colMedians(chunk, na.rm = TRUE)
  # Append the chunk_mean or chunk_median to the alfas_NI_patients data frame
  alfas_NI_patients <- rbind(alfas_NI_patients, chunk_mean)
}
colnames(alfas_NI_patients) <- colnames(alfas_NI)

# Add the patientIDs as a new column to alfas_NI_patients data frame
alfas_NI_patients$patientIDs <- ImageIDs
# Use the %in% operator to check if each value of "patientIDs" is in the list
alfas_NI_patients$surv <- as.integer(alfas_NI_patients$patientIDs %in% long_survivors4000)
# Create df with weight of niches for each long survivors
long_survivors <- alfas_NI_patients[alfas_NI_patients$patientIDs %in% long_survivors4000, ]
long_survivors <- long_survivors[order(long_survivors$patientIDs), ]
row.names(long_survivors) <- NULL
# Create df with weight of niches for each short survivors
short_survivors <- alfas_NI_patients[!(alfas_NI_patients$patientIDs %in% long_survivors4000), ]
short_survivors <- short_survivors[order(short_survivors$patientIDs), ]
row.names(short_survivors) <- NULL



### NICHES barplot
# Select niches and the 'patientIDs' columns
long_survivors_niches <- long_survivors[, c(1:NBNICHES, max(ncol(long_survivors))-1)]
short_survivors_niches <- short_survivors[, c(1:NBNICHES, max(ncol(short_survivors))-1)]
# Exclude the "patientIDs" column for both DataFrames
long_df <- long_survivors_niches %>% select(-patientIDs)
short_df <- short_survivors_niches %>% select(-patientIDs)

combined_df <- rbind(
  cbind(surv = "short_survivor", short_df),
  cbind(surv = "long_survivor", long_df)
)
barplot_table <- combined_df %>%
  pivot_longer(cols = -surv, names_to = "niche", values_to = "alfa")

x_axis_labels <- list()
# Split the data by 'niche' and perform Wilcoxon test for each niche
niche_groups <- split(barplot_table, barplot_table$niche)
# Get the unique niches in the order they appear
unique_niches <- unique(barplot_table$niche)
# Reorder the niche_groups list to match the original order
niche_groups <- niche_groups[unique_niches]

for (niche_name in names(niche_groups)) {
  niche_data <- niche_groups[[niche_name]]
  wilcox_result <- wilcox.test(niche_data$alfa ~ niche_data$surv)
  
  # Extract p-value and create x-axis label with stars
  p_value <- wilcox_result$p.value
  star_coeff <- floor(-log10(p_value))
  star <- strrep("*", star_coeff)
  if (star!=""){
    cat(rep("-", 30), "\n")
    cat(niche_name, star, "\n")
    cat(rep("-", 30), "\n")
  }
  x_axis_labels <- c(x_axis_labels, paste(niche_name, star))
}
ggbarplot(barplot_table, x = "niche", y = "alfa", 
          add = c("mean_se", "point"),
          add.params = list(color = "black", size = 0.5),
          fill = "surv", color = "surv",
          palette = c("lightgreen", "red"),
          position = position_dodge(0.8)) +
  scale_x_discrete(labels = x_axis_labels) +
  xlab(NULL) +  
  ylab("mean(alfas) per patient")



### INTERFACES barplot
# Select interfaces and 'patientIDs' columns
long_survivors_interfaces <- long_survivors[, -c(1:NBNICHES, max(ncol(long_survivors)))]
short_survivors_interfaces <- short_survivors[, -c(1:NBNICHES, max(ncol(short_survivors)))]
# Exclude the "patientIDs" column for both DataFrames
long_df <- long_survivors_interfaces %>% select(-patientIDs)
short_df <- short_survivors_interfaces %>% select(-patientIDs)
# Create a data frame for plotting
combined_df <- rbind(
  cbind(surv = "short_survivor", short_df),
  cbind(surv = "long_survivor", long_df)
)
barplot_table <- combined_df %>%
  pivot_longer(cols = -surv, names_to = "interface", values_to = "alfa")

x_axis_labels <- list()
# Split the data by 'interface' and perform Wilcoxon test for each interface
interface_groups <- split(barplot_table, barplot_table$interface)
# Get the unique interfaces in the order they appear
unique_interfaces <- unique(barplot_table$interface)
# Reorder the interface_groups list to match the original order
interface_groups <- interface_groups[unique_interfaces]

for (interface_name in names(interface_groups)) {
  interface_data <- interface_groups[[interface_name]]
  wilcox_result <- wilcox.test(interface_data$alfa ~ interface_data$surv)
  
  # Extract p-value and create x-axis label with stars
  p_value <- wilcox_result$p.value
  star_coeff <- floor(-log10(p_value))
  star <- strrep("*", star_coeff)
  if (star!=""){
    cat(rep("-", 30), "\n")
    cat(interface_name, star, "\n")
    cat(rep("-", 30), "\n")
  }
  x_axis_labels <- c(x_axis_labels, paste(interface_name, star))
}

ggbarplot(barplot_table, x = "interface", y = "alfa", 
          add = c("mean_se", "point"),
          add.params = list(color = "black", size = 0.5),
          fill = "surv", color = "surv",
          palette = c("lightgreen", "red"),
          position = position_dodge(0.8)) +
  scale_x_discrete(labels = x_axis_labels) +
  xlab(NULL) +  
  ylab("mean(alfas) per patient")




### NICHE-PHENOTYPE MAPPING LONG VS SHORT SURVIVORS

## Niches weights(proportions) of all cells from all images
niches<- paste0("a",as.vector(seq(1,NBNICHES,1)))
cellsNiches <- as_tibble(lapply(json_data4$cells_niches,unlist))%>%
  rename_at(vars(matches("[0-9]")),~niches)%>%
  mutate(cell_id=as.numeric(cell_id))
## Create a new table that will contain niches and interfaces weights
cellsNichesInterfaces <- cellsNiches %>%
  select(-TOT_cell_dens)
# Rename niches
colnames(cellsNichesInterfaces)[1:NBNICHES] <- col_names
# Get all combinations of column names
column_combinations <- combn(col_names, 2, simplify = FALSE)
# Iterate through the column combinations
for (cols in column_combinations) {
  col1 <- cols[1]
  col2 <- cols[2]
  # Create a new column name
  new_col_name <- paste(col1, col2, sep = ".")
  # Create the new column by multiplying the two existing columns
  cellsNichesInterfaces[new_col_name] <- cellsNichesInterfaces[[col1]] * cellsNichesInterfaces[[col2]]
}
# Rename interfaces with shorten names
names(cellsNichesInterfaces)[(ncol(cellsNichesInterfaces) - (choose(NBNICHES, 2)-1)):(ncol(cellsNichesInterfaces))] <- short_interfaces_names

## Associate each cell to a niche or interfaces
# Initialize a variable to store the result
cellsNichesInterfaces$niche <- "undefined"
# Loop through each niche and check the condition
for (niche in col_names) {
  condition <- cellsNichesInterfaces[[niche]] > treshold_niches
  cellsNichesInterfaces$niche[condition] <- niche
}
# Repeat the process for interfaces
for (interface in short_interfaces_names) {
  condition <- cellsNichesInterfaces[[interface]] > treshold_interfaces & cellsNichesInterfaces$niche == "undefined"
  cellsNichesInterfaces$niche[condition] <- interface
}

# Proportion of niches and interfaces
table(cellsNichesInterfaces$niche)



## Define markers to remove
Lineage_markers <- read.csv("./phenotypes_niches/data/proteins_by_frame.csv") %>%filter(Purpose=="Lineage") %>%pull(Biomarker)
# c(B7H3, OX40, CD163, `CSF-1R`)

## Associate cells with niches AND MFI of functional markers
cells.NichesInterface.Phen <- read.csv("./TMENS_analysis/data/cellData.csv",check.names=FALSE,header = TRUE, sep =',')%>%
  dplyr::select(-c(cellSize,Background,C,Na,Si,P,Ca,Fe,immuneCluster,Ta,Au))%>%
  mutate(immuneGroup = recode(immuneGroup,`0`= 'None',`1`='Tregs', `2`='CD4-T',
                              `3`='CD8-T', `4`='CD3-T', `5`='NK',
                              `6`='B', `7`='Neutrophils', `8`='Macrophages', `9`='DC',
                              `10`='DC / Mono', `11`='Mono / Neu', `12`='Other immune')) %>%
  mutate(Group = recode(Group,`1`='Unidentified', `2`='Immune',
                        `3`='Endothelial', `4`='Mesenchymal-like',
                        `5` = 'Tumor',
                        `6` = 'Keratin-positive tumor'))%>%
  mutate(cell_type = ifelse(Group == 'Immune', cell_type<- immuneGroup,cell_type <- Group))%>%
  dplyr::select(-c(tumorYN,tumorCluster,Group,immuneGroup))%>%filter(cell_type!="Unidentified")%>%
  #dplyr::rename(patient_id = SampleID)%>%
  dplyr::rename(cell_id = cellLabelInImage)%>%
  left_join(cellsNichesInterfaces%>%filter(cell_type!="Unidentified"),.,by=c("SampleID","cell_id","cell_type"))%>%
  select(-one_of(Lineage_markers), -c(B7H3, OX40, CD163, `CSF-1R`))



## Continuous analyses

# Define min MFI value  
Niches_Interfaces <- unique(cells.NichesInterface.Phen$niche)
cell_types <- unique(cells.NichesInterface.Phen$cell_type)
Functionnal_markers <- read.csv("./phenotypes_niches/data/proteins_by_frame.csv")%>%filter(Purpose=="Functional")%>%pull(Biomarker)

min_MFI <- Inf
for (NI in Niches_Interfaces) {
  # Filter lines for current niche / interface
  NI_table <- cells.NichesInterface.Phen[cells.NichesInterface.Phen$niche == NI, ]
  for (CT in cell_types) {
    # Filter lines for current niche / interface and cell type
    NI_CT_table <- NI_table[NI_table$cell_type == CT, ]
    for (FM in Functionnal_markers) {
      NI_CT_FM_table <- NI_CT_table[, c('SampleID','cell_id','cell_type','niche',FM)]
      
      LS_NI_CT_FM_table <- subset(NI_CT_FM_table, SampleID %in% long_survivors4000)
      SS_NI_CT_FM_table <- subset(NI_CT_FM_table, !(SampleID %in% long_survivors4000))
      # Filter combination that contain at least 100 cells in long or short survivors
      if (nrow(LS_NI_CT_FM_table) > 100 & nrow(SS_NI_CT_FM_table) > 100) {
        MFI_LS <- unname(unlist(LS_NI_CT_FM_table[, FM]))
        MFI_SS <- unname(unlist(SS_NI_CT_FM_table[, FM]))
        min_val <- min(c(MFI_SS, MFI_LS))
        
        if (min_val < min_MFI) {
          min_MFI <- min_val
        }
        
      }
    }
  }
}


## Create table with p values comparing long and short survivors MFI values of functionnal markers
# For niche and cell type
pvalues <- numeric(0)
log_ratio_LS.SS <- data.frame(
  niche = character(),
  cell_type = character(),
  marker = character(),
  MFI_LS = numeric(),
  MFI_SS = numeric(),
  log_ratioLS = numeric(),
  pvalue = numeric()
)
for (NI in Niches_Interfaces) {
  # Filter lines for current niche / interface
  NI_table <- cells.NichesInterface.Phen[cells.NichesInterface.Phen$niche == NI, ]
  for (CT in cell_types) {
    # Filter lines for current niche / interface and cell type
    NI_CT_table <- NI_table[NI_table$cell_type == CT, ]
    for (FM in Functionnal_markers) {
      NI_CT_FM_table <- NI_CT_table[, c('SampleID','cell_id','cell_type','niche',FM)]
      
      LS_NI_CT_FM_table <- subset(NI_CT_FM_table, SampleID %in% long_survivors4000)
      SS_NI_CT_FM_table <- subset(NI_CT_FM_table, !(SampleID %in% long_survivors4000))
      # Filter combination that contain at least 100 cells in long or short survivors
      if (nrow(LS_NI_CT_FM_table) > 100 & nrow(SS_NI_CT_FM_table) > 100) {
        # Make all MFI values positive
        if (min_MFI< 0) {
          MFI_LS <- unname(unlist(LS_NI_CT_FM_table[, FM])) + (-min_MFI)
          MFI_SS <- unname(unlist(SS_NI_CT_FM_table[, FM])) + (-min_MFI)
        }
        else {
          MFI_LS <- unname(unlist(LS_NI_CT_FM_table[, FM]))
          MFI_SS <- unname(unlist(SS_NI_CT_FM_table[, FM]))
        }
        
        MWtest <- wilcox.test(MFI_LS, MFI_SS)
        pvalues <- c(pvalues, MWtest$p.value)
        
        mean_LS_NI_CT_FM <- mean(MFI_LS)
        mean_SS_NI_CT_FM <- mean(MFI_SS)
        ratioLS <- mean_LS_NI_CT_FM/mean_SS_NI_CT_FM
        log_ratio_LS.SS <- rbind(log_ratio_LS.SS, 
                                 data.frame(niche = NI, cell_type = CT, marker = FM, 
                                            MFI_LS = length(MFI_LS), MFI_SS = length(MFI_SS), log_ratioLS = log10(ratioLS), pvalue = MWtest$p.value))
      }
    }
  }
}
log_ratio_LS.SS$Combination_ID <- seq_len(nrow(log_ratio_LS.SS))


## Number of niches before and after filtering niches that contain at least 100 cells in long or short survivors
niches_uniques_for_MFI_L_and_S <- unique(log_ratio_LS.SS$niche)
# Compute length of both lists
taille_Niches_Interfaces <- length(Niches_Interfaces)
taille_niches_uniques <- length(niches_uniques_for_MFI_L_and_S)
data <- data.frame(
  Liste = c("All niches", "Niches long&short surv"),
  Taille = c(taille_Niches_Interfaces, taille_niches_uniques)
)
# Create barplot with ggplot2
ggplot(data, aes(x = Liste, y = Taille, fill = Liste)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +  
  geom_text(aes(label = Taille), position = position_dodge(width = 0.7), vjust = -0.5) +  #
  labs(title = "Nb of niches with at least 100cells in both long & short survivors",
       x = "",
       y = "Nb niches") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(),  
    axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0)), 
    plot.title = element_text(hjust = 0.25), 
  ) +
  coord_cartesian(ylim = c(0, 15))


# Remove NaN p values (same MFI values for all cells in both long and short survivors)
nan_pvalues <- which(is.nan(log_ratio_LS.SS$pvalue))
log_ratio_LS.SS_ <- log_ratio_LS.SS[complete.cases(log_ratio_LS.SS$pvalue), ]

## Correction for multiple testing
# Compute Q values and remove combination for which it is greater than 0.1
FDR <- fdrtool(log_ratio_LS.SS_$pvalue, statistic="pvalue")
qvalues <- FDR$qval
log_ratio_LS.SS_$qvalue <- qvalues
log_ratio_LS.SS_qvalTRESH <- log_ratio_LS.SS_[log_ratio_LS.SS_$qvalue <= 0.1, ]


## If there is still too much significant combination 
# Select log ratio that are greater or lower than logRatio treshold
tresh_logRatio <- log10(1.3)
log_ratio_LS.SS_qvalTRESH_logratioTRESH <- log_ratio_LS.SS_qvalTRESH[log_ratio_LS.SS_qvalTRESH$log_ratioLS > tresh_logRatio 
                                                                     | log_ratio_LS.SS_qvalTRESH$log_ratioLS < -tresh_logRatio, ]


## Heatmap du log ratio long vs short for significant q values
# For combination CT + marker that have logRatio greater or lower than the logRatio treshold, add the niches combination values that are significant
table_heatmap <- log_ratio_LS.SS_qvalTRESH_logratioTRESH %>%
  distinct(cell_type, marker) %>%
  left_join(log_ratio_LS.SS_qvalTRESH, by = c("cell_type", "marker"))
# For the heatmap: Get niche as row, CT+marker as column and log_ratioLS as values
table_heatmap <- table_heatmap %>%
  mutate(cell_type_marker = paste(cell_type, marker, sep = "_")) %>%
  select(niche, log_ratioLS, cell_type_marker, Combination_ID)
# Create a reference table with all possible combination of CT+marker with niches
# In order to have niche combination for every CT+marker (also the ones that are not significant -> logratio value set to 0)
reference_table <- expand.grid(
  cell_type_marker = unique(table_heatmap$cell_type_marker),
  niche = niches_uniques_for_MFI_L_and_S
)
# Join with table_heatmap (that contain logratio value, if the combination is not in table_heatmap -> logratio set to 0)
table_heatmap_complete <- reference_table %>%
  left_join(table_heatmap, by = c("cell_type_marker", "niche")) %>%
  mutate(log_ratioLS = coalesce(log_ratioLS, 0))  # Remplacer les valeurs manquantes par 0
# Create heatmap
heatmap_plot <- table_heatmap_complete %>%
  ggplot(aes(x = cell_type_marker, y = niche, fill = log_ratioLS)) +
  geom_tile() +
  scale_fill_gradient2(low = "#66A3FF", mid = "white", high = "#FF6666", midpoint = 0) +
  labs(title = "Heatmap of log_ratioLS",
       x = "Cell Type Marker",
       y = "Niche") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  geom_hline(yintercept = seq(0.5, length(unique(table_heatmap_complete$niche)) + 0.5), color = "#777777", linetype = "solid", linewidth = 1) +
  geom_vline(xintercept = seq(0.5, length(unique(table_heatmap_complete$cell_type_marker)) + 0.5), color = "#777777", linetype = "solid", linewidth = 1)
# Plot heatmap
print(heatmap_plot)


## plot distribution of MFI values for a specific combination marker-celltype-niche
# You can find the Comb_ID corresponding to a heatmap value / combination marker-celltype-niche in table_heatmap_complete$Combination_ID
Comb_ID <- 207
i = 1
for (NI in Niches_Interfaces) {
  # Filter lines for current niche / interface
  NI_table <- cells.NichesInterface.Phen[cells.NichesInterface.Phen$niche == NI, ]
  for (CT in cell_types) {
    # Filter lines for current niche / interface and cell type
    NI_CT_table <- NI_table[NI_table$cell_type == CT, ]
    for (FM in Functionnal_markers) {
      NI_CT_FM_table <- NI_CT_table[, c('SampleID','cell_id','cell_type','niche',FM)]
      
      LS_NI_CT_FM_table <- subset(NI_CT_FM_table, SampleID %in% long_survivors4000)
      SS_NI_CT_FM_table <- subset(NI_CT_FM_table, !(SampleID %in% long_survivors4000))
      # Filter combination that contain at least 100 cells in long or short survivors
      if (nrow(LS_NI_CT_FM_table) > 100 & nrow(SS_NI_CT_FM_table) > 100) {
        if (i == Comb_ID){
          print(i)
          print(NI)
          print(CT)
          print(FM)
          MFI_LS <- unname(unlist(LS_NI_CT_FM_table[, FM])) + (-min_MFI)
          MFI_SS <- unname(unlist(SS_NI_CT_FM_table[, FM])) + (-min_MFI)
        }
        i <- i + 1
      }
    }
  }
}
# Data frame for distribution long and short survivor for Comb_ID
table_MFI_LSandSS <- data.frame(
  Groupe = rep(c("MFI_LS", "MFI_SS"), times = c(length(MFI_LS), length(MFI_SS))),
  Valeur = c(MFI_LS, MFI_SS)
)
# Plot density graph
ggplot(table_MFI_LSandSS, aes(x = Valeur, fill = Groupe)) +
  geom_density(alpha = 0.5) +
  labs(title = "Distribution of MFI values",
       x = "MFI",
       y = "Density") +
  scale_fill_manual(values = c("green", "red"))

