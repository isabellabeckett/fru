## ----setup, include=FALSE------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ------------------------------------------------------------------------------------------------------------
library(googlesheets4)
library(malecns)
library(dplyr)
library(nat.nblast)
library(fafbseg)
library(tidyverse)
library(coconatfly)
library(dendextend)
library(ggplot2)
library(arrow)
library(purrr)








## ------------------------------------------------------------------------------------------------------------
rds_dir <- "/net/flystore3/jdata/jdata5/JPeople/Bella/fru_clone_uploads/clone_sheets"
# Directory to save the dendrogram plots
conn_dir <- "/net/flystore3/jdata/jdata5/JPeople/Bella/fru_clone_uploads/conn_tables"

# Get all .rds filenames in the directory
rds_files <- list.files(path = rds_dir, pattern = "\\.rds$", full.names = TRUE)
rds_files <- rds_files[grepl("FW_.*", basename(rds_files))]

# Process each .rds file, read, and add a clone column
fw_combined_data <- rds_files %>%
  map_dfr(~ {
    clone <- tools::file_path_sans_ext(basename(.x))
    fru_sheets <- readRDS(.x)
    
    # Convert all columns to character to ensure consistency
    fru_sheets <- fru_sheets %>%
      mutate(across(everything(), as.character))
    
    # Add the clone column
    fru_sheets %>% mutate(clone = clone)
  })


FW.Y <- fw_combined_data %>% filter(Fru=="Y") %>% select(root_id)
  


## ------------------------------------------------------------------------------------------------------------
FW_all_hemilineages <- fw_combined_data %>% filter(!is.na(hemilineage)) %>% select(clone) %>% distinct()
FW_complete_hemilineages <- fw_combined_data %>% filter(clone %in% FW_all_hemilineages$clone)


## ------------------------------------------------------------------------------------------------------------
FW_clones_with_both_fru <- FW_complete_hemilineages %>%
  group_by(clone) %>%
  filter(any(Fru == "Y") & any(Fru == "N")) %>%
  pull(clone) %>%
  unique()

FW_complete_YandN <- FW_complete_hemilineages %>%
  filter(clone %in% FW_clones_with_both_fru)


## ------------------------------------------------------------------------------------------------------------
compute_totals_indv_pre <- function(df, direction, threshold = 2) {
  choose_segmentation('public')
  connection_table <- flywire_partner_summary(df, partners = direction, threshold = threshold)
  connection_table$query <- df
  connection_table <- connection_table %>%
  group_by(query) %>%
  mutate(total_weight = sum(weight)) %>%
  ungroup()
  connection_table$dimorphic_partner <- ifelse(connection_table$pre_id %in% FW.Y$root_id, "Y", "N")
  FW_dimorphic_summary <- connection_table %>%
  group_by(query, dimorphic_partner) %>%
  summarise(weight = sum(weight), .groups = 'drop') %>% 
    filter(dimorphic_partner == "Y")

# Step 2: Extract total_weight where dimorphic_partner is "Y"
dimorphic_total <- connection_table %>%
  filter(dimorphic_partner == "Y") %>%
  select(query, dimorphic_partner, total_weight) %>%
  distinct() # Ensure unique rows if there are duplicates

# Step 3: Join the summarized weight data with the total_weight
FWdm2 <- FW_dimorphic_summary %>%
  left_join(dimorphic_total, by = c("query", "dimorphic_partner"))

  
  return(FWdm2)
}

compute_totals_indv_post <- function(df, direction, threshold = 2) {
  choose_segmentation('public')
  connection_table <- flywire_partner_summary(df, partners = direction, threshold = threshold)
  connection_table$query <- df
  connection_table <- connection_table %>%
  group_by(query) %>%
  mutate(total_weight = sum(weight)) %>%
  ungroup()
  connection_table$dimorphic_partner <- ifelse(connection_table$post_id %in% FW.Y$root_id, "Y", "N")
  FW_dimorphic_summary <- connection_table %>%
  group_by(query, dimorphic_partner) %>%
  summarise(weight = sum(weight), .groups = 'drop') %>% 
    filter(dimorphic_partner == "Y")

dimorphic_total <- connection_table %>%
  filter(dimorphic_partner == "Y") %>%
  select(query, dimorphic_partner, total_weight) %>%
  distinct() 

# Step 3: Join the summarized weight data with the total_weight
FWdm2 <- FW_dimorphic_summary %>%
  left_join(dimorphic_total, by = c("query", "dimorphic_partner"))

  
  return(FWdm2)
}


## ------------------------------------------------------------------------------------------------------------
FW_complete_hemilineages$percent_dimorphic_input <- NA
FW_complete_hemilineages$percent_dimorphic_output <- NA

  # Filter data for the current clone
  FW_clone_Y_table <- FW_complete_hemilineages %>% filter(Fru == "Y") 
  FW_clone_N_table <- FW_complete_hemilineages %>% filter(Fru == "N") 
  FW_clone_Y_root_id <- FW_clone_Y_table$root_id
  FW_clone_N_root_id <- FW_clone_N_table$root_id
  FW_clone_Y <- as.character(FW_clone_Y_root_id)
  FW_clone_N <- as.character(FW_clone_N_root_id)


## Input to Fru neurons
s <- length(FW_clone_Y)
FW_clone_Y_in <- NULL
for(i in 1:s)
{
  if (i %% 10 == 0) {
    message(paste("Processing neuron", i))
  }
  
  FW_clone_Y_temp <- tryCatch({
    suppressWarnings({
      suppressMessages({
    compute_totals_indv_pre(FW_clone_Y[i], "i")
    })
  })
  }, error = function(e) {
    message(paste("Error for neuron", i, "- Skipping neuron"))
    return(NULL)  # Return NULL if there's an error
  })
  
  # Only bind the result if it's not NULL
  if (!is.null(FW_clone_Y_temp)) {
    FW_clone_Y_in <- rbind(FW_clone_Y_in, FW_clone_Y_temp)
  }
}  

## Input to non-Fru neurons
s <- length(FW_clone_N)
FW_clone_N_in <- NULL
for(i in 1:s)
{
  if (i %% 10 == 0) {
    message(paste("Processing neuron", i))
  }
  
  FW_clone_N_temp <- tryCatch({
    suppressWarnings({
      suppressMessages({
    compute_totals_indv_pre(FW_clone_N[i], "i")
    })
  })
  }, error = function(e) {
    message(paste("Error for neuron", i, "- Skipping neuron"))
    return(NULL)
  })
  if (!is.null(FW_clone_N_temp)) {
    FW_clone_N_in <- rbind(FW_clone_N_in, FW_clone_N_temp)
  }
}  

## Output to Fru neurons
s <- length(FW_clone_Y)
FW_clone_Y_out <- NULL
for(i in 1:s)
{
  if (i %% 10 == 0) {
    message(paste("Processing neuron", i))
  }
  
  FW_clone_Y_temp <- tryCatch({
    suppressMessages({
      suppressWarnings({
    compute_totals_indv_post(FW_clone_Y[i], "o")
    })
  })
  }, error = function(e) {
    message(paste("Error for neuron", i, "- Skipping neuron"))
    return(NULL)
  })
  if (!is.null(FW_clone_Y_temp)) {
    FW_clone_Y_out <- rbind(FW_clone_Y_out, FW_clone_Y_temp)
  }
}  

##Output to non-Fru neurons
s <- length(FW_clone_N)
FW_clone_N_out <- NULL
for(i in 1:s)
{
  if (i %% 10 == 0) {
    message(paste("Processing neuron", i))
  }
  
  FW_clone_N_temp <- tryCatch({
    suppressWarnings({
      suppressMessages({
    compute_totals_indv_post(FW_clone_N[i], "o")
    })
  })
  }, error = function(e) {
    message(paste("Error for neuron", i, "- Skipping neuron"))
    return(NULL)
  })
  if (!is.null(FW_clone_N_temp)) {
    FW_clone_N_out <- rbind(FW_clone_N_out, FW_clone_N_temp)
  }
}  
 
FW_complete_hemilineages2 <- FW_complete_hemilineages %>%
  filter(!is.na(Fru)) %>%
  distinct(root_id, .keep_all = T)

s <- unlist(FW_complete_hemilineages2[["root_id"]]) 
for(i in s) {
  # Filter the neuron row
  neuron <- FW_complete_hemilineages2 %>% filter(root_id == i)
  index <- which(FW_complete_hemilineages2$root_id == i)
  # Check the condition
  if (neuron$Fru == "Y") {
     if (i %in% FW_clone_Y_in[["query"]]){
    # Find the index of the current neuron in the data frame
    corresponding <- FW_clone_Y_in %>% filter(query == i)
    # Update the specific column
    FW_complete_hemilineages2[index, 6] <- 100*corresponding$weight / corresponding$total_weight
     } else {
      FW_complete_hemilineages2[index, 6] = 0
    }
  } else if (neuron$Fru == "N") {
    if(i %in% FW_clone_N_in[["query"]]){
    corresponding <- FW_clone_N_in %>% filter(query == i)
    FW_complete_hemilineages2[index, 6] <- 100*corresponding$weight / corresponding$total_weight
    } else {
      FW_complete_hemilineages2[index, 6] = 0
    }
  } else {
      FW_complete_hemilineages2[index, 6] = 0
  }
}

for(i in s) {
  # Filter the neuron row
  neuron <- FW_complete_hemilineages2 %>% filter(root_id == i)
  index <- which(FW_complete_hemilineages2$root_id == i)
  # Check the condition
  if (neuron$Fru == "Y") {
     if (i %in% FW_clone_Y_out[["query"]]){
    # Find the index of the current neuron in the data frame
    corresponding <- FW_clone_Y_out %>% filter(query == i)
    # Update the specific column
    FW_complete_hemilineages2[index, 7] <- 100*corresponding$weight / corresponding$total_weight
     } else {
      FW_complete_hemilineages2[index, 7] = 0
    }
  } else if (neuron$Fru == "N") {
    if(i %in% FW_clone_N_out[["query"]]){
    corresponding <- FW_clone_N_out %>% filter(query == i)
    FW_complete_hemilineages2[index, 7] <- 100*corresponding$weight / corresponding$total_weight
    } else {
      FW_complete_hemilineages2[index, 7] = 0
    }
  } else {
      FW_complete_hemilineages2[index, 7] = 0
  }
}




## ------------------------------------------------------------------------------------------------------------
# Ensure the directory exists (you can skip this if you know it exists)
if (!dir.exists(conn_dir)) {
  dir.create(conn_dir, recursive = TRUE)
}

# Save each table as an .rds file
saveRDS(FW_clone_Y_in, file = file.path(conn_dir, "FW_clone_Y_in.rds"))
saveRDS(FW_clone_N_in, file = file.path(conn_dir, "FW_clone_N_in.rds"))
saveRDS(FW_clone_Y_out, file = file.path(conn_dir, "FW_clone_Y_out.rds"))
saveRDS(FW_clone_N_out, file = file.path(conn_dir, "FW_clone_N_out.rds"))
saveRDS(FW_complete_hemilineages2, file = file.path(conn_dir, "FW_complete_hemilineages2.rds"))

