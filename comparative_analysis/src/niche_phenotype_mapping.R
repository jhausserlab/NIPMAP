

createInterfaces <- function(cellsNichesInterfaces, col_names, short_interfaces_names) {
  
  # Rename niches
  colnames(cellsNichesInterfaces)[1:NBNICHES] <- col_names
  
  # Get all combinations of niches
  column_combinations <- combn(col_names, 2, simplify = FALSE)
  
  # Iterate through the niches combinations to create interfaces columns
  for (cols in column_combinations) {
    col1 <- cols[1]
    col2 <- cols[2]
    new_col_name <- paste(col1, col2, sep = ".")
    cellsNichesInterfaces[new_col_name] <- cellsNichesInterfaces[[col1]] * cellsNichesInterfaces[[col2]]
  }
  
  # Rename interfaces with shorten names
  names(cellsNichesInterfaces)[(ncol(cellsNichesInterfaces) - (choose(NBNICHES, 2)-1)):(ncol(cellsNichesInterfaces))] <- short_interfaces_names
  
  return(cellsNichesInterfaces)
}


associateCellsToNichesInterfaces <- function(cellsNichesInterfaces, col_names, treshold_niches, short_interfaces_names, treshold_interfaces) {
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
  
  return(cellsNichesInterfaces)
}

associateCellsToFunctionalMarkers <- function(cellsNichesInterfaces, Unwanted_markers) {
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
    select(-one_of(Unwanted_markers))
  
  return(cells.NichesInterface.Phen)
}

compute_minMFI <- function(cells.NichesInterface.Phen, long_survivors4000, Niches_Interfaces, cell_types, Functionnal_markers) {
  # Initialize min_MFI to positive infinity
  min_MFI <- Inf
  
  # Loop over Niches_Interfaces
  for (NI in Niches_Interfaces) {
    # Filter lines for the current niche/interface
    NI_table <- cells.NichesInterface.Phen[cells.NichesInterface.Phen$niche == NI, ]
    
    # Loop over cell_types
    for (CT in cell_types) {
      # Filter lines for the current niche/interface and cell type
      NI_CT_table <- NI_table[NI_table$cell_type == CT, ]
      
      # Loop over Functionnal_markers
      for (FM in Functionnal_markers) {
        # Extract relevant columns for the current combination
        NI_CT_FM_table <- NI_CT_table[, c('SampleID', 'cell_id', 'cell_type', 'niche', FM)]
        
        # Subset data for long survivors and short survivors
        LS_NI_CT_FM_table <- subset(NI_CT_FM_table, SampleID %in% long_survivors4000)
        SS_NI_CT_FM_table <- subset(NI_CT_FM_table, !(SampleID %in% long_survivors4000))
        
        # Filter combination that contains at least 100 cells in long or short survivors
        if (nrow(LS_NI_CT_FM_table) > 100 & nrow(SS_NI_CT_FM_table) > 100) {
          # Extract MFI values
          MFI_LS <- unname(unlist(LS_NI_CT_FM_table[, FM]))
          MFI_SS <- unname(unlist(SS_NI_CT_FM_table[, FM]))
          
          # Find the minimum MFI value
          min_val <- min(c(MFI_SS, MFI_LS))
          
          # Update min_MFI if the current min_val is smaller
          if (min_val < min_MFI) {
            min_MFI <- min_val
          }
        }
      }
    }
  }
  
  # Return the minimum MFI value
  return(min_MFI)
}



compute_logRatio_and_pvalues <- function(cells.NichesInterface.Phen, long_survivors4000, min_MFI, Niches_Interfaces, cell_types, Functionnal_markers) {
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
  return(log_ratio_LS.SS)
}


# Number of niches before and after filtering niches that contain at least 100 cells in long or short survivors
plot_nicheCount <- function(Niches_Interfaces, niches_uniques_for_MFI_L_and_S) {
  
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
    labs(title = "Nb of niches with at least 100 cells in both long & short survivors",
         x = "",
         y = "Nb niches") +
    theme_minimal() +
    theme(
      axis.text.y = element_text(),  
      axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0)), 
      plot.title = element_text(hjust = 0.25), 
    ) +
    coord_cartesian(ylim = c(0, 15))
}

# Heatmap du log ratio long vs short for significant q values
heatmap_logRatio_LvsS_significantQvalues <- function(log_ratio_LS.SS_qvalTRESH_logratioTRESH, log_ratio_LS.SS_qvalTRESH, niches_uniques_for_MFI_L_and_S) {
  # For combination CT + marker that have logRatio greater or lower than the logRatio threshold, add the niches combination values that are significant
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
  
  # Join with table_heatmap (that contains logratio value, if the combination is not in table_heatmap -> logratio set to 0)
  table_heatmap_complete <- reference_table %>%
    left_join(table_heatmap, by = c("cell_type_marker", "niche")) %>%
    mutate(log_ratioLS = coalesce(log_ratioLS, 0))  # Replace missing values with 0
  
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
  
  # Print heatmap
  print(heatmap_plot)
}


# Distribution of MFI values for a specific combination marker-celltype-niche
plot_MFI_distribution <- function(Comb_ID, cells.NichesInterface.Phen, long_survivors4000, min_MFI, Niches_Interfaces, cell_types, Functionnal_markers) {
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
}
  
