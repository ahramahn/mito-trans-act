rm(list =ls())
wkdir <- "/Users/ahramahn/OneDrive - University of Miami/1. Research/DMS-MaP-Seq/Yeast-trans_act/COX1/new_results/PET309"
setwd(wkdir) ; getwd() ; dir()

library(dplyr)

# List of file names
file_names <- c("PET309_COX1_5UTR", "PET309_COX1_CDS", "PET309_COX1_3UTR")

# Loop over the file names
for (file in file_names) {
  # Construct the file path
  file_path <- paste(file, "ratio.csv", sep="_")
  
  # Read the CSV file
  df <- read.csv(file_path, header = TRUE)
  
  # Apply the mutate operation
  df <- df %>% 
    mutate(Mutated = ifelse(Base %in% c("T", "G") | Mutated > 0.2, NA, Mutated))
  
  # Assign the mutated data frame to a variable with the same name as the file
  assign(file, df, envir = .GlobalEnv)
}

# Define the processing function
process_sample <- function(df) {
  
  # Combine reactivity data for nucleotide A and C
  A_C_reactivity <- c(df[df$Base == "A", "Mutated"], df[df$Base == "C", "Mutated"])
  A_C_reactivity <- as.numeric(A_C_reactivity)
  A_C_reactivity <- A_C_reactivity[!is.na(A_C_reactivity)]
  
  # Get the median of the top 10% of A and C combined
  top_10_percent_A_C <- sort(A_C_reactivity, decreasing = TRUE)[1:floor(0.1 * length(A_C_reactivity))]
  median_A_C <- median(top_10_percent_A_C)
  
  # Normalization function
  normalize_reactivity <- function(nucleotide, dms_react) {
    # Convert to numeric; if not possible, return NA
    reactivity <- suppressWarnings(as.numeric(dms_react))
    if (is.na(reactivity)) return(NA)
    
    if (nucleotide == "A" || nucleotide == "C") {
      return(reactivity / median_A_C)
    } else {
      return(NA) # return NA for other nucleotides
    }
  }
  
  # Apply the normalization function to the data and create a new column "normalized_DMS"
  df$normalized_DMS <- mapply(normalize_reactivity, df$Base, df$Mutated)
  
  # Cap the normalized_DMS values at 1
  df$normalized_DMS <- pmin(df$normalized_DMS, 1)
  
  return(df)
}

PET309_COX1_5UTR_processed <- process_sample(PET309_COX1_5UTR)
PET309_COX1_3UTR_processed <- process_sample(PET309_COX1_3UTR)
PET309_COX1_CDS_processed <- process_sample(PET309_COX1_CDS)



# List of processed data frame names
processed_names <- c("PET309_COX1_5UTR_processed", "PET309_COX1_CDS_processed", "PET309_COX1_3UTR_processed")

# Loop to add SHAPE_reactivity to each processed data frame
for (name in processed_names) {
  # Retrieve the processed data frame
  processed_df <- get(name)
  
  # Add the SHAPE_reactivity column
  processed_df$SHAPE_reactivity <- ifelse(processed_df$Base %in% c('G', 'T') | is.na(processed_df$normalized_DMS), -999, processed_df$normalized_DMS)
  
  # Reassign the modified data frame back to its variable
  assign(name, processed_df, envir = .GlobalEnv)
}


# FOR VARNA

PET309_COX1_5UTR_processed_varna <- process_sample(PET309_COX1_5UTR)
PET309_COX1_3UTR_processed_varna <- process_sample(PET309_COX1_3UTR)
PET309_COX1_CDS_processed_varna <- process_sample(PET309_COX1_CDS)



processed_names2 <- c("PET309_COX1_5UTR_processed_varna", "PET309_COX1_CDS_processed_varna", "PET309_COX1_3UTR_processed_varna")

for (name in processed_names2) {
  # Retrieve the processed data frame
  processed_df <- get(name)
  
  # Add the SHAPE_reactivity column
  processed_df$SHAPE_reactivity <- ifelse(processed_df$Base %in% c('G', 'T') | is.na(processed_df$normalized_DMS), 0, processed_df$normalized_DMS)
  
  # Reassign the modified data frame back to its variable
  assign(name, processed_df, envir = .GlobalEnv)
}

all_processed_names <- c("PET309_COX1_5UTR_processed", "PET309_COX1_CDS_processed", "PET309_COX1_3UTR_processed", 
"PET309_COX1_5UTR_processed_varna", "PET309_COX1_CDS_processed_varna", "PET309_COX1_3UTR_processed_varna")

# Loop to save each processed data frame
for (name in all_processed_names) {
  # Retrieve the processed data frame
  processed_df <- get(name)
  
  # Construct file name for saving
  file_name <- paste(name, "final.SHAPE", sep = "_")
  
  # Save the file
  write.table(processed_df[, c('Position', 'SHAPE_reactivity')],
              file_name, 
              quote = FALSE, 
              row.names = FALSE, 
              col.names = FALSE, 
              sep = "\t")
}
