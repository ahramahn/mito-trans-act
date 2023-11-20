rm(list =ls())
getwd()
wkdir <- "/Users/ahramahn/Dropbox/Mac/Documents/dreem_seq/Yeast-trans_act/COX1"
setwd(wkdir) ; getwd() ; dir()


library(tidyverse)

# Read the data
P1 <- read.csv("Corrected_P1-related_m-ratio_seqbar.csv", header =TRUE)
M1 <- read.csv("Corrected_MSS51-1_related_m-ratio_seqbar.csv", header = TRUE)
M2 <- read.csv("Corrected_MSS51-2_related_m-ratio_seqbar.csv", header = TRUE)
M3.1 <- read.csv("Corrected_M3-1_related_m-ratio_seqbar.csv", header =TRUE)
M3.2 <- read.csv("Corrected_M3-2_related_m-ratio_seqbar.csv", header =TRUE)
M4 <- read.csv("Corrected_M4_related_m-ratio_seqbar.csv", header = TRUE)
WT1 <- read.csv("Corrected_WT1-related_m-ratio_seqbar.csv", header = TRUE)
WT2 <- read.csv("Corrected_WT2_related_m-ratio_seqbar.csv", header = TRUE)


# Average WT1 and WT2
M1_M2_avg <- data.frame(
  Position = WT1$Position,
  Base = WT1$Base,
  Mutated = (M1$Mutated + M2$Mutated) / 2
)


M3.1_2_avg <- data.frame(
  Position = M3.1$Position,
  Base = M3.1$Base,
  Mutated = (M3.1$Mutated + M3.2$Mutated) / 2
)

# Average (M3.1, M3.2 combined) and M4
M3_M4_avg <- data.frame(
  Position = M3.1$Position,
  Base = M3.1$Base,
  Mutated = (M3.1_2_avg$Mutated + M4$Mutated) / 2
)

# Average WT1 and WT2
WT_avg <- data.frame(
  Position = WT1$Position,
  Base = WT1$Base,
  Mutated = (WT1$Mutated + WT2$Mutated) / 2
)

# Set values lower than 0.005 to zero for M3.1_2_avg, M3_M4_avg, and WT_avg
P1$Mutated[P1$Mutated < 0.005] <- 0
M1_M2_avg$Mutated[M1_M2_avg$Mutated < 0.005] <- 0
M3.1_2_avg$Mutated[M3.1_2_avg$Mutated < 0.005] <- 0
M3_M4_avg$Mutated[M3_M4_avg$Mutated < 0.005] <- 0
WT_avg$Mutated[WT_avg$Mutated < 0.005] <- 0

P1$Mutated[P1$Mutated > 0.3] <- NA
M1_M2_avg$Mutated[M1_M2_avg$Mutated > 0.3] <- NA
M3_M4_avg$Mutated[M3_M4_avg$Mutated > 0.3] <- NA
WT_avg$Mutated[WT_avg$Mutated > 0.3] <- NA

P1 <- P1[1:2065, ]
M1_M2_avg <- M1_M2_avg[1:2065, ]
M3_M4_avg <- M3_M4_avg[1:2065, ]
WT_avg <- WT_avg[1:2065, ]

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


# Apply the function to each data frame in the list
P1_process_sample <- process_sample(P1)
M1_M2_process_samples <- process_sample(M1_M2_avg)
WT_processed_samples <- process_sample(WT_avg)
M3_M4_processed_samples <- process_sample(M3_M4_avg)

WT_dot <- "...............((((((((((((((..........))))))))))))))...............(((((((.((((((((((.((((((.....)))))))))).)))))).))))))).((((((((((((((..(((((((........)))))))...))))))))))))))...................((((((((.((((......))))))))))))..........................................(((((((((((............(((((.....((((((((((((.........(((((....))))).........))))))))))))...))))).((((((((......))))))))............))))))))))).............................................................................................................(((((.....)))))..((((..........)))).((((((....)))))).....((((((((((((((.............((((((...(((.....((((((((((((..((((..((((((((((.....................((((.((......))))))..((((((.......)))))).........(......).))))))))))..))))..))))).)))))))....)))...))))))...((((((...((((......)))).))))))....))))))))))).)))...(((...........))).................(((...........))).(((((((..........))))))).............((((((((((.(((.((((.(((.((((..(((((((............((((((((.(((((((...............)))))))...(((((((.......)))).)))..))))))))............((((.(((.....))).)))).)))))))..)))).))).......((((..(((((((.....)))))))...)))).....((((.((((...................))))....)))).(((((.((((...........)))).))))).((((((((((.......)))))))))).............)))))))..))))))))))...........(((((.......(((((((.......................................((((((.........................))))))...(((((..(((((..((((((((....))))))..................))..)))))..))))).)))))))....))))).(((((((.....((((....)))).................)))))))........((((((((....)).)))))).....(((.(((((((((...((((.....(((((((((((.....))))))))))).))))...)))))))))..)))....(((((.(((((............))))))))))......................................................((((((((...............))))))))......((((((((.......))))))))..((((..(((((.((((((.........))))))...(((((...((((.(((((....((((((((.(((((((((((........................))))).)).))))))))))))....)))))))))..)))))..)))))..)))).(((....))).........................................((((((......................)))))).......((((((..................))))))..........."



MSS51_dot <- ".(((.((((((....((((((((((((((..........))))))))))))))...............(((((((.((((((((((.((((((.....)))))))))).)))))).))))))).((((((((((((((..(((((((........)))))))...))))))))))))))...................((((((((.((((......))))))))))))....................)))))).))).((((((((...(((((((((((............(((((.....(((((((((((((........(((((....)))))........)))))))))))))...))))).((((((((......))))))))............))))))))))).......)))))))).............................................................((((......))))...................(((.......(((((((((................(((((((....)))))))....((((((((((((((((...........((((((...(((.....((((((((((((..((((..((((((((((.....................((((.((......))))))..((((((.......)))))).........(......).))))))))))..))))..))))).)))))))....)))...))))))...((((((...((((......)))).)))))).)).))))))))))).)))...(((...........))).))).)))))).......)))............(((((((((..........)))))))))...........((((((((((.(((.((((.(((.((((..(((((((............((((((((.(((((((...............)))))))...(((((((.......))).))))..))))))))............((((.(((.....))).)))).)))))))..)))).))).......((((..(((((((.....)))))))...)))).....((((.((((...................))))....)))).(((((.((((...........)))).))))).((((((((((.......)))))))))).............)))))))..)))))))))).....((((..(((((.......(((((((.......................................((((((.........................))))))...(((((..(((((..((((((((....))))))..................))..)))))..))))).)))))))....))))).)))).((((((((.((....))..............((((.....)))).......((((((((....)).))))))........)))))))).....((((.....(((((((((((.....))))))))))).))))...(((((((...........(((((.(((((............))))))))))......................................................((((((((...............))))))))......((((((((.......))))))))............((((((((.((....)))))))))).(((((...((((.(((((....((((((((.(((((((((((........................))))).)).))))))))))))....)))))))))..))))).....(((.......)))..)))))))......................................((((((......................)))))).......((((((..................))))))..........."

MSS116_dot <- "((((.((((((....((((((((((((((..........))))))))))))))...............(((((((.((((((((((.((((((.....)))))))))).)))))).))))))).((((((((((((((..(((((((........)))))))...)))))))))))))).......))))))))))..((((((((.((((......))))))))))))..........................................(((((((((((......................(((((((((((((........(((((....)))))........)))))))))))))........(((((((((......)))))))))...........))))))))))).......................................................................................................................((((...(((((.............(((((((....)))))))....((((((((((((((.............((((((...(((.....((((((((((((..((((..((((((((((.....................((((.((......))))))..(((((.........))))).........(......).))))))))))..))))..))))).)))))))....)))...))))))...((((((...((((......)))).))))))....))))))))))).)))..........)))))...)))).................................(((((((..........)))))))..............(((((((((.(((.((((.(((.((((..((((((((...........((((((((.(((((((...............)))))))...(((((((.......))).))))..)))))))).............(((.(((.....))).))).))))))))..)))).))).......((((..(((((((.....)))))))...)))).....((((.((...........................)))))).(((((.((((...........)))).))))).((((((((((.......)))))))))).............)))))))..)))))))))......(((....((((.((..................................................((((((.........................))))))...................((((((....))))))..((((.......)))).((((((.....))).))).((((...........((((((......((((....))))..................))))))...........)))))).))))....))).....(((.(((((((((...((((.....(((((((((((.....))))))))))).))))...)))))))))..)))....(((((.(((((............)))))))))).((((((((.............................................((((((((...............))))))))......((((((((.......))))))))..((((..((((((.(((((.........)))))....(((((...((((.(((((....((((((((.(((((((((((........................))))).)).))))))))))))....)))))))))..))))).))))))..)))).(((....)))............)))))))).....................((((((......................)))))).......((((((..................))))))..........."

PET309_dot <- "...............((((((((((((((..........))))))))))))))...............(((((((.((((((((((.((((((.....)))))))))))).)))).))))))).((((((((((((((..(((((((........)))))))...))))))))))))))...................((((((((.((((......))))))))))))..........................................(((((((((((......................(((((((((((((........(((((....)))))........)))))))))))))........(((((((((......)))))))))...........))))))))))).........................................................................................................(((((((((.....))))).)))).......((......(((((((....)))))))....)).(((((((((((.............((((((...(((.....((((((((((((..((((..((((((((((.....................((((.((......))))))..(((((.........))))).........(......).))))))))))..))))..))))).)))))))....)))...))))))...((((((...((((......)))).))))))....))))))))))).(((...(((((...((((.......................)))))))))...))).((((((((..........))))))))............((((((((((.(((.((((.(((.((((..((((((((...........((((((((.(((((((...............)))))))..((((((((.......))).))))).)))))))).............(((.(((.....))).))).))))))))..)))).))).......((((..(((((((.....)))))))...))))....(((((.((((...................))))....))))).((((.((((...........)))).))))..((((((((((.......)))))))))).............)))))))..)))))))))).....((((..(((((.......(((((((.......................................((((((.........................))))))...(((((..(((((..((((((((....))))))..................))..)))))..))))).)))))))....))))).))))........((((....)))).............(((.....)))........((((((((....)).)))))).....(((.(((((((((...((((.....(((((((((((.....))))))))))).))))...)))))))))..)))....(((((.(((((............))))))))))......................................................((((((((...............))))))))......((((((((.......))))))))..((((..((((((((((((.........))))))...(((((...((((.(((((....((((((((.(((((((((((........................))))).)).))))))))))))....)))))))))..))))).))))))..)))).(((....))).........................................((((((......................)))))).......((((((..................))))))..........."


WT_dot_unlist <- unlist(strsplit(WT_dot, ""))

# Transform each character to "Negative" for unpaired bases and "Positive" for paired bases
WT_transformed <- sapply(WT_dot_unlist, function(x) {
  if (x == ".") {
    return ("FALSE")
  } else if (x == "(" || x == ")") {
    return ("TRUE")
  } else {
    return(NA)  # This handles any other unexpected characters
  }
})

# Append the transformed vector to the original dataset
WT_processed_samples$STRUCT <- WT_transformed



# MSS51 : the variable name is a bit different, so catious
M1_M2_dot_unlist <- unlist(strsplit(MSS51_dot, ""))

M1_M2_transformed <- sapply(M1_M2_dot_unlist, function(x) {
  if (x == ".") {
    return ("FALSE")
  } else if (x == "(" || x == ")") {
    return ("TRUE")
  } else {
    return(NA)  # This handles any other unexpected characters
  }
})

# Append the transformed vector to the original dataset
M1_M2_process_samples$STRUCT <- M1_M2_transformed



# MSS116
M3_M4_dot_unlist <- unlist(strsplit(MSS116_dot, ""))

# Transform each character to "Negative" for unpaired bases and "Positive" for paired bases
M3_M4_transformed <- sapply(M3_M4_dot_unlist, function(x) {
  if (x == ".") {
    return ("FALSE")
  } else if (x == "(" || x == ")") {
    return ("TRUE")
  } else {
    return(NA)  # This handles any other unexpected characters
  }
})

# Append the transformed vector to the original dataset
M3_M4_processed_samples$STRUCT <- M3_M4_transformed


P1_dot_unlist <- unlist(strsplit(PET309_dot, ""))

# Transform each character to "Negative" for unpaired bases and "Positive" for paired bases
P1_transformed <- sapply(P1_dot_unlist, function(x) {
  if (x == ".") {
    return ("FALSE")
  } else if (x == "(" || x == ")") {
    return ("TRUE")
  } else {
    return(NA)  # This handles any other unexpected characters
  }
})

# Append the transformed vector to the original dataset
P1_process_sample$STRUCT <- P1_transformed





# Classification function based on a threshold
classify_based_on_threshold <- function(df, threshold) {
  classification_col <- paste("classification_", threshold, sep="")
  df[[classification_col]] <- NA
  
  # FP condition
  indices <- which(df$STRUCT == "FALSE" & df[["normalized_DMS"]] < threshold)
  df[indices, classification_col] <- "FN"
  
  # TN condition
  indices <- which(df$STRUCT == "TRUE" & df[["normalized_DMS"]] >= threshold)
  df[indices, classification_col] <- "TN"
  
  # FN condition
  indices <- which(df$STRUCT == "FALSE" & df[["normalized_DMS"]] >= threshold)
  df[indices, classification_col] <- "FP"
  
  # TP condition
  indices <- which(df$STRUCT == "TRUE" & df[["normalized_DMS"]] < threshold)
  df[indices, classification_col] <- "TP"
  
  return(df)
}


# Thresholds
thresholds <- seq(0.05, 1, by=0.05)

# Apply function for each threshold for all samples
data_frames <- list(WT_processed_samples, M1_M2_process_samples, M3_M4_processed_samples, P1_process_sample)

# make a function
apply_thresholds <- function(data) {
  for (t in thresholds) {
    data <- classify_based_on_threshold(data, t)
  }
  return(data)
}

data_frames <- lapply(data_frames, apply_thresholds)
# extract data from data frame
WT_processed_samples <- data_frames[[1]]
M1_M2_process_samples <- data_frames[[2]]
M3_M4_processed_samples <- data_frames[[3]]
P1_process_sample <- data_frames[[4]]

process_dataset <- function(data) {
# set the window size
window_size <- 100
half_window <- window_size / 2

# Function to compute counts for a window and a threshold
get_counts <- function(window, classification_col) {
  TP <- sum(window[[classification_col]] == "TP", na.rm=TRUE)
  TN <- sum(window[[classification_col]] == "TN", na.rm=TRUE)
  FP <- sum(window[[classification_col]] == "FP", na.rm=TRUE)
  FN <- sum(window[[classification_col]] == "FN", na.rm=TRUE)
  return(list(TP=TP, TN=TN, FP=FP, FN=FN))
}


# Sliding window loop
for (i in 51:(nrow(data) - 49)) {
  window <- data[(i - 50):(i + 49), ] # Get the 80nt window
  
  for (j in thresholds) {
    classification_col <- paste("classification_", j, sep="")
    counts <- get_counts(window, classification_col)
    
    # Assign counts to the central nucleotide
    data[i, paste("TP_", j, sep="")] <- counts$TP
    data[i, paste("TN_", j, sep="")] <- counts$TN
    data[i, paste("FP_", j, sep="")] <- counts$FP
    data[i, paste("FN_", j, sep="")] <- counts$FN
  }
}


# Function to compute TPR and FPR
tpr <- function(TP, FN) {
  ifelse(is.na(TP) | is.na(FN) | TP + FN == 0, NA, TP / (TP + FN))
}

fpr <- function(FP, TN) {
  ifelse(is.na(FP) | is.na(TN) | FP + TN == 0, NA, FP / (FP + TN))
}

# Compute TPR and FPR for each threshold
for (t in thresholds) {
  TP_col <- paste("TP_", t, sep="")
  FP_col <- paste("FP_", t, sep="")
  FN_col <- paste("FN_", t, sep="")
  TN_col <- paste("TN_", t, sep="")
  
  TPR_col <- paste("TPR_", t, sep="")
  FPR_col <- paste("FPR_", t, sep="")
  
  data[, TPR_col] <- tpr(data[[TP_col]], data[[FN_col]])
  data[, FPR_col] <- fpr(data[[FP_col]], data[[TN_col]])
}

# Function to compute AUC for ROC curve
compute_auc_roc <- function(fpr_points, tpr_points) {
  auc = 0
  for(i in 2:length(fpr_points)) {
    auc = auc + (fpr_points[i] - fpr_points[i-1]) * (tpr_points[i] + tpr_points[i-1]) / 2
  }
  return(auc)
}

# Compute AUROC for each nucleotide and store in a new column
data$AUC_ROC_per_nucleotide <- apply(data, 1, function(row) {
  fpr_points <- as.numeric(row[grep("FPR_", names(row))])
  tpr_points <- as.numeric(row[grep("TPR_", names(row))])
  
  # Check if there are any missing values, and remove them
  valid_indices <- which(!is.na(fpr_points) & !is.na(tpr_points))
  fpr_points <- fpr_points[valid_indices]
  tpr_points <- tpr_points[valid_indices]
  
  # Sort by FPR
  sorted_indices <- order(fpr_points)
  fpr_points <- fpr_points[sorted_indices]
  tpr_points <- tpr_points[sorted_indices]
  
  return(compute_auc_roc(fpr_points, tpr_points))
})


list_columns <- sapply(data, is.list)
for (col in names(data)[list_columns]) {
  data[[col]] <- sapply(data[[col]], toString)
}

data$AUC_ROC_per_nucleotide <- sapply(data$AUC_ROC_per_nucleotide, function(x) if(length(x) == 0) NA else x)

return(data)
}


# Run a loop
data_frames_names <- c("WT_processed_samples", "M1_M2_process_samples", "M3_M4_processed_samples", "P1_process_sample")
data_frames_list <- list(WT_processed_samples, M1_M2_process_samples, M3_M4_processed_samples, P1_process_sample)

processed_data_frames <- lapply(data_frames_list, process_dataset)

# If you want to assign the processed datasets back to their original variable names:
for (i in 1:length(data_frames_names)) {
  assign(data_frames_names[i], processed_data_frames[[i]], envir = .GlobalEnv)
}

write.csv(WT_processed_samples, "WT_COX1_yeast_all.csv")
write.csv(M1_M2_process_samples, "Mss51_COX1_yeast_all.csv")
write.csv(M3_M4_processed_samples, "Mss116_COX1_yeast_all.csv")
write.csv(P1_process_sample, "Pet309_COX1_yeast_all.csv")

# make a figure

library(ggplot2)



create_plot <- function(plot_data_roc, plot_title_prefix, save_filename_prefix) {
  plot_data_roc$Position <- as.numeric(as.character(plot_data_roc$Position))
  plot_data_roc$AUC_ROC_per_nucleotide <- as.numeric(plot_data_roc$AUC_ROC_per_nucleotide)
  avg_auc_roc <- mean(plot_data_roc$AUC_ROC_per_nucleotide, na.rm = TRUE)
  p <- ggplot(plot_data_roc, aes(x = Position, y = AUC_ROC_per_nucleotide)) + 
    geom_line(color = "blue", na.rm = TRUE) +
    geom_hline(yintercept = avg_auc_roc, linetype = "dashed", color = "red") + 
    labs(title = paste0(plot_title_prefix, " AUROC per position"),
         y = "AUROC Value", x = "Position (nt)",
         subtitle = paste("Average AUROC:", round(avg_auc_roc, 2))) +
    theme_minimal() +
    scale_y_continuous(limits = c(0.5, 1.2)) +
    scale_x_continuous(breaks = seq(0, max(plot_data_roc$Position, na.rm = TRUE), by = 100)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(paste0(save_filename_prefix, "_AUROC_plot.png"), p, bg = "white")
}

datasets <- list(WT_processed_samples, M1_M2_process_samples, M3_M4_processed_samples, P1_process_sample)
titles <- c("WT", "M1_M2", "M3_M4", "P1")
filenames <- c("WT", "M1_M2", "M3_M4", "P1")

for (i in 1:length(datasets)) {
  create_plot(datasets[[i]], titles[i], filenames[i])
}



### the overall AUCROC values were too low
### that seems to be that the same vales across FPR makes a NA value so it is not applied, even though it should be. 