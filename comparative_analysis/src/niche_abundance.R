library(tidyverse)
library(ggpubr)

Niche_Interfaces_ShortVSLong_survivors <- function(alfas, col_names, short_interfaces_names, NSITES, long_survivors4000, ImageIDs) {

  # Convert the list to a data frame with specified column names
  alfas_NI <- data.frame(`colnames<-`(do.call(cbind, alfas), col_names))
  
  # Get all combinations of column names
  column_combinations <- combn(col_names, 2, simplify = FALSE)
  
  # Create interfaces weights
  for (cols in column_combinations) {
    col1 <- cols[1]
    col2 <- cols[2]
    
    # Create a new column name
    new_col_name <- paste(col1, col2, sep = ".")
    
    # Create the new column by multiplying the two existing columns
    alfas_NI[new_col_name] <- alfas_NI[[col1]] * alfas_NI[[col2]]
  }
  
  # Rename interfaces with shortened names
  names(alfas_NI)[(ncol(alfas_NI) - (choose(length(col_names), 2) - 1)):(ncol(alfas_NI))] <- short_interfaces_names
  
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
  
  # Create df with the weight of niches for each long survivor
  long_survivors <- alfas_NI_patients[alfas_NI_patients$patientIDs %in% long_survivors4000, ]
  long_survivors <- long_survivors[order(long_survivors$patientIDs), ]
  row.names(long_survivors) <- NULL
  
  # Create df with the weight of niches for each short survivor
  short_survivors <- alfas_NI_patients[!(alfas_NI_patients$patientIDs %in% long_survivors4000), ]
  short_survivors <- short_survivors[order(short_survivors$patientIDs), ]
  row.names(short_survivors) <- NULL
  
  # Return the processed data
  return(list(long_survivors = long_survivors, short_survivors = short_survivors))
}



# Barplot niche / interfaces abundance
barplot_N.I_abundance <- function(long_df, short_df, group_variable) {
  
  # Exclude the "patientIDs" column for both DataFrames
  long_df <- long_df %>% select(-patientIDs)
  short_df <- short_df %>% select(-patientIDs)
  
  # Create a data frame for plotting
  combined_df <- rbind(
    cbind(surv = "short_survivor", short_df),
    cbind(surv = "long_survivor", long_df)
  )
  
  # Pivot long based on the provided group variable
  barplot_table <- combined_df %>%
    pivot_longer(cols = -surv, names_to = group_variable, values_to = "alfa")
  
  x_axis_labels <- list()
  
  # Split the data by the provided group variable and perform Wilcoxon test for each group
  group_variable_groups <- split(barplot_table, barplot_table[[group_variable]])
  
  # Get the unique values of the group variable in the order they appear
  unique_groups <- unique(barplot_table[[group_variable]])
  
  # Reorder the group_variable_groups list to match the original order
  group_variable_groups <- group_variable_groups[unique_groups]
  
  for (group_value in names(group_variable_groups)) {
    group_data <- group_variable_groups[[group_value]]
    wilcox_result <- wilcox.test(group_data$alfa ~ group_data$surv)
    
    # Extract p-value and create x-axis label with stars
    p_value <- wilcox_result$p.value
    star_coeff <- floor(-log10(p_value))
    star <- strrep("*", star_coeff)
    
    if (star != "") {
      cat(rep("-", 30), "\n")
      cat(group_value, star, "\n")
      cat(rep("-", 30), "\n")
    }
    
    x_axis_labels <- c(x_axis_labels, paste(group_value, star))
  }
  
  # Plotting
  ggbarplot(barplot_table, x = group_variable, y = "alfa", 
            add = c("mean_se", "point"),
            add.params = list(color = "black", size = 0.5),
            fill = "surv", color = "surv",
            palette = c("lightgreen", "red"),
            position = position_dodge(0.8)) +
    scale_x_discrete(labels = x_axis_labels) +
    xlab(NULL) +  
    ylab(paste("mean(alfas) per patient - Grouped by", group_variable))
}
