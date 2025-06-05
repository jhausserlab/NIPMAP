createInterfaces <- function(cellsNichesInterfaces, niche_cols = NULL) {
  
  # Automatically detect niches if not specified
  if (is.null(niche_cols)) {
    niche_cols <- c("Cancer", "Inflammatory", "Bfollicle", "Other", "Lowdensity")
  }
  
  NBNICHES <- length(niche_cols)
  
  # Get original column names to preserve
  original_cols <- colnames(cellsNichesInterfaces)
  
  # Rename niche columns to standard order if needed
  colnames(cellsNichesInterfaces)[match(niche_cols, colnames(cellsNichesInterfaces))] <- niche_cols
  
  # Create pairwise interfaces
  pairwise_combos <- combn(niche_cols, 2, simplify = FALSE)
  for (cols in pairwise_combos) {
    new_name <- paste(cols, collapse = ".")
    cellsNichesInterfaces[[new_name]] <- cellsNichesInterfaces[[cols[1]]] * cellsNichesInterfaces[[cols[2]]]
  }
  
  # Create triple interfaces
  triple_combos <- combn(niche_cols, 3, simplify = FALSE)
  for (cols in triple_combos) {
    new_name <- paste(cols, collapse = ".")
    cellsNichesInterfaces[[new_name]] <- cellsNichesInterfaces[[cols[1]]] * cellsNichesInterfaces[[cols[2]]] * cellsNichesInterfaces[[cols[3]]]
  }
  # Assign short names to the new columns
  interface_start <- length(original_cols) + 1
  interface_end <- ncol(cellsNichesInterfaces)
  
  return(cellsNichesInterfaces)
}




associateCellsToNichesInterfaces <- function(cellsNichesInterfaces, col_names, treshold_niches, short_interfaces_names, treshold_interfaces, triple_interfaces_names, treshold_triple_interfaces) {
  # Initialize a variable to store the result
  cellsNichesInterfaces$niche <- "mixed"
  
  # Loop through each niche and check the condition
  for (niche in col_names) {
    condition <- cellsNichesInterfaces[[niche]] > treshold_niches
    cellsNichesInterfaces$niche[condition] <- niche
  }
  
  # Repeat the process for interfaces
  for (interface in short_interfaces_names) {
    condition <- (cellsNichesInterfaces[[interface]] > treshold_interfaces) & (cellsNichesInterfaces$niche == "mixed")
    cellsNichesInterfaces$niche[condition] <- interface
  }
  
  # Repeat the process for triple interfaces
  for (triple_interface in triple_interfaces_names) {
    condition <- (cellsNichesInterfaces[[triple_interface]] > treshold_triple_interfaces) & (cellsNichesInterfaces$niche == "mixed")
    cellsNichesInterfaces$niche[condition] <- triple_interface
  }
  
  # Proportion of niches and interfaces
  table(cellsNichesInterfaces$niche)
  
  return(cellsNichesInterfaces)
}

associateCellsToFunctionalMarkersBI3 <- function(cellsNichesInterfaces, directory_path) {
  file_names <- list.files(directory_path, pattern = "\\.csv$", full.names = TRUE)
  # Initialize an empty list to store data frames
  data_list <- list()
  # Loop through each file
  for (file in file_names) {
    file_name <- basename(file)
    sample_id <- as.numeric(gsub("[^0-9]", "", file_name))  # Extract numeric SampleID
    data <- read.csv(file)  
    # Check if necessary columns exist
    if ("Phenotype" %in% colnames(data) & "label" %in% colnames(data)) {
      temp_data <- data[, c("Phenotype", "label", "x", "y")]
      temp_data$SampleID <- sample_id
      # Create unique key by merging SampleID and label
      temp_data$UniqueKey <- paste0(sample_id, "_", as.character(temp_data$label))
      data_list[[length(data_list) + 1]] <- temp_data
    }
  }
  # Combine all data frames into one
  final_table <- do.call(rbind, data_list)
  # Create a unique key in cellsNichesInterfaces
  cellsNichesInterfaces$UniqueKey <- paste0(cellsNichesInterfaces$SampleID, "_", format(as.numeric(cellsNichesInterfaces$cell_id), scientific = FALSE, trim = TRUE))
  # Perform a left join to add the Phenotype column to cellsNichesInterfaces
  cellsNichesInterfaces <- merge(cellsNichesInterfaces, final_table[, c("UniqueKey", "Phenotype", "x", "y")], 
                                 by = "UniqueKey", all.x = TRUE)
  
  return(cellsNichesInterfaces)
}





compute_interactions_and_pval <- function(cells.NichesInterface.Phen, Niches_Interfaces, possible_interactions, survival_data){
  survival_data_allinteractions <- survival_data
  for (NI in Niches_Interfaces) {
    print(NI)
    NI_table <- cells.NichesInterface.Phen[cells.NichesInterface.Phen$niche == NI, ]
    
    # convert to data.table
    DT  <- as.data.table(NI_table)
    PI  <- as.data.table(possible_interactions)
    
    # prepare output list
    out <- vector("list", length = 0)
    
    # loop over samples
    sample_ids <- sort(unique(DT$SampleID))
    pb <- progress_bar$new(
      total = length(sample_ids),
      format = "  Processing [:bar] :current/:total (:percent) - Patient :patient",
      clear = FALSE, width = 60
    )
    
    for (s in sample_ids) {
      pb$tick(tokens = list(patient = s))
      subDT <- DT[SampleID == s]
      coords <- as.matrix(subDT[, .(x, y)])
      
      # build a radius‐neighbor index for this sample
      fr <- frNN(coords, eps = 50)  # eps in same units as x,y
      
      # precompute cell‐type indices
      idx_by_type <- split(seq_len(nrow(subDT)), subDT$cell_type)
      
      # for each requested interaction
      for(i in seq_len(nrow(PI))) {
        ref_type <- PI$Reference[i]
        tgt_type <- PI$Target[i]
        
        ref_idx <- idx_by_type[[ref_type]] %||% integer(0)
        tgt_idx <- idx_by_type[[tgt_type]] %||% integer(0)
        
        if (length(ref_idx)==0 || length(tgt_idx)==0) {
          conn <- 0L
          norm <- 0
        } else {
          # All neighbor IDs of ref_idx
          neighbors_list <- fr$id[ref_idx]
          
          # Flatten all neighbor indices to one vector
          all_neighbors <- unlist(neighbors_list, use.names = FALSE)
          
          # Efficient counting: how many of them are target cells
          conn <- sum(all_neighbors %in% tgt_idx)
          
          norm <- conn / length(tgt_idx)
        }
        
        out[[length(out)+1]] <- list(
          SampleID         = s,
          Reference        = ref_type,
          Target           = tgt_type,
          interactions     = conn,
          norm_interactions = norm
        )
      }
    }
    
    # bind into one data.table
    result_DT <- rbindlist(out)
    
    
    norm_long <- result_DT %>%
      mutate(
        Reference_Target = paste0(
          "Interactions_",
          Reference, "_",
          Target, "_", NI
        )
      ) %>%
      select(SampleID, Reference_Target, norm_interactions)
    
    
    norm_wide <- norm_long %>%
      pivot_wider(
        names_from  = Reference_Target,
        values_from = norm_interactions
      ) %>%
      rename(patient_id = SampleID)
    
    survival_data_allinteractions <- merge(survival_data_allinteractions, norm_wide, by = "patient_id", all.x = TRUE)
  }
  ## Compute cox model on each interactions
  # 1) Which columns are your interactions?
  interaction_cols <- setdiff(
    names(survival_data_allinteractions),
    c("patient_id", "PFS_months", "event_status")
  )
  
  # 2) Fit one Cox model per interaction and extract HR + p-value
  cox_results <- map_dfr(interaction_cols, function(col) {
    # build formula with backticks
    f <- as.formula(paste0("Surv(PFS_months, event_status) ~ `", col, "`"))
    
    # fit safely
    fit <- tryCatch(coxph(f, data = survival_data_allinteractions),
                    error = function(e) NULL)
    # Compute median number of interactions for this column
    median_nb_interactions <- median(survival_data_allinteractions[[col]], na.rm = TRUE)
    
    if (is.null(fit)) {
      return(tibble(
        interaction = col,
        hr          = NA_real_,
        p.value     = NA_real_,
        median_nb_interactions = median_nb_interactions
      ))
    }
    
    # pull out the 1×5 numeric vector: coef, exp(coef), se, z, Pr(>|z|)
    stats <- summary(fit)$coefficients[1, ]
    
    tibble(
      interaction = col,
      hr          = stats["exp(coef)"],
      p.value     = stats["Pr(>|z|)"],
      median_nb_interactions = median_nb_interactions
    )
  })
  
  # 3) Sort by p-value
  cox_results %>% arrange(p.value)
  
  return(list(cox_results, survival_data_allinteractions))
}



compute_count_logRatio_and_pvaluesBI3_final <- function(cells.NichesInterface.Phen, long_survivors4000, Niches_Interfaces, cell_types, Functionnal_markers, unique_sample_ids, survival_data) {
  pvalues <- numeric(0)
  log_ratio_LS.SS <- data.frame(
    niche = character(),
    cell_type = character(),
    marker = character(),
    Density_R = numeric(),
    Density_NR = numeric(),
    log_ratio_density = numeric(),
    pvalue = numeric()
  )
  survival_data_allratios <- survival_data
  
  for (NI in Niches_Interfaces) {
    print(NI)
    # Filter lines for current niche / interface
    NI_table <- cells.NichesInterface.Phen[cells.NichesInterface.Phen$niche == NI, ]
    for (CT in cell_types) {
      # Filter lines for current niche / interface and cell type
      NI_CT_table <- NI_table[NI_table$cell_type == CT, ]
      for (FM in Functional_markers) {
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
          rename(num_cells_NI = num_cells) %>%  # Rename num_cells column
          left_join(NI_CT_FM_table_num_cell %>% rename(num_cells_NI_CT_FM = num_cells), by = "SampleID")
        # Remove patients where number of cells in the niche is less than 100
        # count_table <- count_table %>%
        #   filter(num_cells_NI >= 100)
        
        count_table$ratio_num_cells <- count_table$num_cells_NI_CT_FM / count_table$num_cells_NI
        
        
        ratio_data <- data.frame(
          patient_id = count_table$SampleID,
          ratio = count_table$ratio_num_cells  # Example ratio values
        )
        
        if (exists("survival_data") && nrow(ratio_data) >= 1) {
          survival_data_ <- merge(survival_data, ratio_data, by = "patient_id")
          FM_concatenated <- paste(FM, collapse = ",")
          FM_clean <- gsub(",", "_", gsub("/", "", FM_concatenated))
          ratio_col_name <- paste0("Ratio_", CT, "_", FM_clean, "_", NI)
          colnames(survival_data_)[colnames(survival_data_) == "ratio"] <- ratio_col_name
          
          # Find the column with "Ratio" in the name
          ratio_col <- grep("Ratio", names(survival_data_), value = TRUE)
          survival_data_allratios <- merge(survival_data_allratios, survival_data_, by = c("patient_id", "PFS_months", "event_status"), all.x = TRUE)
          
          # Build formula safely, handling special characters
          cox_formula <- as.formula(paste("Surv(PFS_months, event_status) ~", paste0("`", ratio_col, "`")))
          
          # Try fitting the Cox model, handle error gracefully
          p_value <- tryCatch({
            # Try to run the Cox model
            cox_model <- coxph(cox_formula, data = survival_data_)
            
            # Try to extract p-value
            p_value <- summary(cox_model)$coefficients[, "Pr(>|z|)"]
            
          }, error = function(e) {
            # If error happens, print message and the merged survival data
            message("❌ Cox model failed: ", e$message)
            message("🔍 Merged survival_data_ that caused the error:")
            print(survival_data_)
            
            # Return NA or NULL so your code continues
            return(NA)
          })
          
          
        } else {
          LS_count_table <- subset(count_table, SampleID %in% long_survivors4000)
          SS_count_table <- subset(count_table, !(SampleID %in% long_survivors4000))
          
          
          LS_filtered <- na.omit(LS_count_table$ratio_num_cells)
          SS_filtered <- na.omit(SS_count_table$ratio_num_cells)
          
          # Check again the length after removing NA values
          length(LS_filtered)
          length(SS_filtered)
          
          # Only run Wilcoxon test if both groups have at least two observations
          if (length(LS_filtered) >= 1 && length(SS_filtered) >= 1) {
            p_value <- wilcox.test(LS_filtered, SS_filtered, exact = FALSE)$p.value
          } else {
            p_value <- "No niche for one group of patients"
          }
          
          
          ## Add pseudo count for logratio
          # The pseudo count is computed based on the total number of cells in the niche
          # so if the niche has low number of cells the pseudo count don't impact too much compare to just adding 1
          count_table_for_logratio <- NI_table_num_cell %>%
            rename(num_cells_NI = num_cells) %>%  # Rename num_cells column
            left_join(NI_CT_FM_table_num_cell %>% rename(num_cells_NI_CT_FM = num_cells), by = "SampleID")
          # Remove patients where number of cells in the niche is less than 100
          # count_table_for_logratio <- count_table_for_logratio %>%
          #   filter(num_cells_NI >= 100)
          ## Add pseudo count max
          # max_nb_cells_in_niche <- max(count_table_for_logratio$num_cells_NI)
          # count_table_for_logratio$num_cells_NI_CT_FM <- count_table_for_logratio$num_cells_NI_CT_FM + (count_table_for_logratio$num_cells_NI / max_nb_cells_in_niche)
          # count_table_for_logratio$ratio_num_cells <- count_table_for_logratio$num_cells_NI_CT_FM / count_table_for_logratio$num_cells_NI
          
          ## Add pseudo count 1
          count_table_for_logratio$num_cells_NI_CT_FM <- count_table_for_logratio$num_cells_NI_CT_FM + 1
          count_table_for_logratio$ratio_num_cells <- count_table_for_logratio$num_cells_NI_CT_FM / count_table_for_logratio$num_cells_NI
          
          LS_count_table <- subset(count_table_for_logratio, SampleID %in% long_survivors4000)
          SS_count_table <- subset(count_table_for_logratio, !(SampleID %in% long_survivors4000))
          
          # Compute the mean density for each group
          LS_count <- median(LS_count_table$ratio_num_cells, na.rm = TRUE)
          SS_count <- median(SS_count_table$ratio_num_cells, na.rm = TRUE)
          # Calculate the log ratio (log of density for long survivors divided by short survivors)
          ratio_count <- LS_count / SS_count
        }
        
        

        
        # Create the data frame and append it to log_ratio_LS.SS
        # log_ratio_LS.SS <- rbind(log_ratio_LS.SS, 
        #                          data.frame(niche = NI, cell_type = CT, marker = FM_concatenated, 
        #                                     log_ratioLS = log10(ratio_count), 
        #                                     pvalue = p_value))
        log_ratio_LS.SS <- rbind(log_ratio_LS.SS, 
                                 data.frame(niche = NI, cell_type = CT, marker = FM_concatenated,
                                            nb_cells = as.integer(median(count_table$num_cells_NI_CT_FM[!is.nan(count_table$ratio_num_cells)], na.rm = TRUE)),
                                            hr = exp(coef(cox_model)), 
                                            pvalue = p_value))
      }
    }
  }
  log_ratio_LS.SS$Combination_ID <- seq_len(nrow(log_ratio_LS.SS))
  return(list(log_ratio_LS.SS = log_ratio_LS.SS, survival_data = survival_data_allratios))
}









compute_count_logRatio_and_pvaluesBI3_allcells <- function(cells.NichesInterface.Phen, long_survivors4000, Niches_Interfaces, cell_types, Functionnal_markers, unique_sample_ids, survival_data = NULL) {
  pvalues <- numeric(0)
  log_ratio_LS.SS <- data.frame(
    niche = character(),
    cell_type = character(),
    marker = character(),
    Density_R = numeric(),
    Density_NR = numeric(),
    log_ratio_density = numeric(),
    pvalue = numeric()
  )
  survival_data_allratios <- survival_data
  for (NI in Niches_Interfaces) {
    print(NI)
    # Filter lines for current niche / interface
    NI_table <- cells.NichesInterface.Phen[cells.NichesInterface.Phen$niche == NI, ]
    for (CT in cell_types) {
      # Filter lines for current niche / interface and cell type
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
        rename(num_cells_NI = num_cells) %>%  # Rename num_cells column
        left_join(NI_CT_FM_table_num_cell %>% rename(num_cells_NI_CT_FM = num_cells), by = "SampleID")
      # Remove patients where number of cells in the niche is less than 100
      # count_table <- count_table %>%
      #   filter(num_cells_NI >= 100)
      
      count_table$ratio_num_cells <- count_table$num_cells_NI_CT_FM / count_table$num_cells_NI
      
      
      ratio_data <- data.frame(
        patient_id = count_table$SampleID,
        ratio = count_table$ratio_num_cells  # Example ratio values
      )
      
      if (exists("survival_data") && nrow(ratio_data) >= 1) {
        survival_data_ <- merge(survival_data, ratio_data, by = "patient_id")
        ratio_col_name <- paste0("Ratio_", CT, "_", NI)
        colnames(survival_data_)[colnames(survival_data_) == "ratio"] <- ratio_col_name
        
        # Find the column with "Ratio" in the name
        ratio_col <- grep("Ratio", names(survival_data_), value = TRUE)
        survival_data_allratios <- merge(survival_data_allratios, survival_data_, by = c("patient_id", "PFS_months", "event_status"), all.x = TRUE)
        
        # Build formula safely, handling special characters
        cox_formula <- as.formula(paste("Surv(PFS_months, event_status) ~", paste0("`", ratio_col, "`")))
        
        # Try fitting the Cox model, handle error gracefully
        p_value <- tryCatch({
          # Try to run the Cox model
          cox_model <- coxph(cox_formula, data = survival_data_)
          
          # Try to extract p-value
          p_value <- summary(cox_model)$coefficients[, "Pr(>|z|)"]
          
        }, error = function(e) {
          # If error happens, print message and the merged survival data
          message("❌ Cox model failed: ", e$message)
          message("🔍 Merged survival_data_ that caused the error:")
          print(survival_data_)
          
          # Return NA or NULL so your code continues
          return(NA)
        })

        
      } else {
        LS_count_table <- subset(count_table, SampleID %in% long_survivors4000)
        SS_count_table <- subset(count_table, !(SampleID %in% long_survivors4000))
        
        
        LS_filtered <- na.omit(LS_count_table$ratio_num_cells)
        SS_filtered <- na.omit(SS_count_table$ratio_num_cells)
        
        # Check again the length after removing NA values
        length(LS_filtered)
        length(SS_filtered)
        
        # Only run Wilcoxon test if both groups have at least two observations
        if (length(LS_filtered) >= 1 && length(SS_filtered) >= 1) {
          p_value <- wilcox.test(LS_filtered, SS_filtered, exact = FALSE)$p.value
        } else {
          p_value <- "No niche for one group of patients"
        }
        
        
        ## Add pseudo count for logratio
        # The pseudo count is computed based on the total number of cells in the niche
        # so if the niche has low number of cells the pseudo count don't impact too much compare to just adding 1
        count_table_for_logratio <- NI_table_num_cell %>%
          rename(num_cells_NI = num_cells) %>%  # Rename num_cells column
          left_join(NI_CT_FM_table_num_cell %>% rename(num_cells_NI_CT_FM = num_cells), by = "SampleID")
        # Remove patients where number of cells in the niche is less than 100
        # count_table_for_logratio <- count_table_for_logratio %>%
        #   filter(num_cells_NI >= 100)
        ## Add pseudo count max
        # max_nb_cells_in_niche <- max(count_table_for_logratio$num_cells_NI)
        # count_table_for_logratio$num_cells_NI_CT_FM <- count_table_for_logratio$num_cells_NI_CT_FM + (count_table_for_logratio$num_cells_NI / max_nb_cells_in_niche)
        # count_table_for_logratio$ratio_num_cells <- count_table_for_logratio$num_cells_NI_CT_FM / count_table_for_logratio$num_cells_NI
        
        ## Add pseudo count 1
        count_table_for_logratio$num_cells_NI_CT_FM <- count_table_for_logratio$num_cells_NI_CT_FM + 1
        count_table_for_logratio$ratio_num_cells <- count_table_for_logratio$num_cells_NI_CT_FM / count_table_for_logratio$num_cells_NI
        
        LS_count_table <- subset(count_table_for_logratio, SampleID %in% long_survivors4000)
        SS_count_table <- subset(count_table_for_logratio, !(SampleID %in% long_survivors4000))
        
        # Compute the mean density for each group
        LS_count <- median(LS_count_table$ratio_num_cells, na.rm = TRUE)
        SS_count <- median(SS_count_table$ratio_num_cells, na.rm = TRUE)
        # Calculate the log ratio (log of density for long survivors divided by short survivors)
        ratio_count <- LS_count / SS_count
      }
      
      

      

      
      # Create the data frame and append it to log_ratio_LS.SS
      # log_ratio_LS.SS <- rbind(log_ratio_LS.SS, 
      #                          data.frame(niche = NI, cell_type = CT, 
      #                                     LS_nb_cells = as.integer(sum(LS_count_table$num_cells_NI_CT_FM)), 
      #                                     SS_nb_cells = as.integer(sum(SS_count_table$num_cells_NI_CT_FM)), 
      #                                     log_ratioLS = log10(ratio_count), 
      #                                     pvalue = p_value))
      log_ratio_LS.SS <- rbind(log_ratio_LS.SS, 
                               data.frame(niche = NI, cell_type = CT, 
                                          # LS_nb_cells = as.integer(sum(LS_count_table$num_cells_NI_CT_FM)), 
                                          # SS_nb_cells = as.integer(sum(SS_count_table$num_cells_NI_CT_FM)), 
                                          nb_cells = as.integer(median(count_table$num_cells_NI_CT_FM[!is.nan(count_table$ratio_num_cells)], na.rm = TRUE)), 
                                          hr = exp(coef(cox_model)), 
                                          pvalue = p_value))
    }
  }
  log_ratio_LS.SS$Combination_ID <- seq_len(nrow(log_ratio_LS.SS))
  return(list(log_ratio_LS.SS = log_ratio_LS.SS, survival_data = survival_data_allratios))
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
  
