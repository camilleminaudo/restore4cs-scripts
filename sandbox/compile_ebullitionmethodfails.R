
# ---
# Authors: Miguel Cabrera
# Project: "RESTORE4Cs"
# date: "Dec 2023"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
# This script is used to retireve the incubations where ebullitionVSdiffussion approach for CH4 produces weird results(i.e. NAs or Cero flux in any of total estimate, ebullition or diffusion)

#Inputs: 
#per-season CSV files with computed fluxes
#fieldsheets 

rm(list = ls()) # clear workspace
cat("/014") # clear console


# ---- packages ----
library(tidyverse)
library(readxl)
library(lubridate)
library(zoo)
library(ggplot2)
library(grid)
library(egg)
require(dplyr)
require(purrr)
require(data.table)
require(tools)

repo_root <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
# files.sources = list.files(path = paste0(repo_root,"/functions"), full.names = T)
# for (f in files.sources){source(f)}


# ---- Directories and data loading ----
dropbox_root <- "C:/Users/Miguel/Dropbox/RESTORE4Cs - Fieldwork/Data" # You have to make sure this is pointing to the write folder on your local machine
results_path <- paste0(dropbox_root,"/GHG/Processed data/computed_flux/level_incubation/")

plots_path <- paste0(dropbox_root,"/GHG/Processed data/computed_flux/plots")

diagnostic_path<- paste0(dropbox_root, "/GHG/Working data/EbullitionVSdiffusion diagnosis/")

#Load csv files produced by raw2flux.R. script
setwd(results_path)

listf <- list.files(path = results_path, pattern = "ch4_fluxes", all.files = T, full.names = T, recursive = F)
table_results_all <- NULL

for (f in listf){
  table_results_all <- rbind(table_results_all, 
                             read.csv(file = f, header = T))}

#Add
table_results_all$sampling <- str_sub(table_results_all$UniqueID, start = 1, 2)
table_results_all$pilotsite <- str_sub(table_results_all$UniqueID, start = 4, 5)
table_results_all$subsite <- str_sub(table_results_all$UniqueID, start = 7, 8)
table_results_all$siteID <- str_sub(table_results_all$UniqueID, start = 4, 8)

Nabullition<- table_results_all %>% 
  filter(is.na(ebullition+diffusion+total_estimated))

ceroflux<- table_results_all %>% 
  filter(total_estimated==0)


write.csv(Nabullition, file = paste0(diagnostic_path, "ebullition_NAs.csv"),row.names = F)
write.csv(ceroflux, file = paste0(diagnostic_path, "cerototalflux.csv"), row.names = F)

Nabullition$UniqueID


#Proportion of incubations where method fails: ~52%
(dim(Nabullition)[1]+dim(ceroflux)[1])/dim(table_results_all)[1]


