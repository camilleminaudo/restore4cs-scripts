
# ---
# Authors: Miguel Cabrera
# Project: "RESTORE4Cs"
# date: "Dec 2023"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
# This script is used to retireve the incubations where ebullitionVSdiffussion approach for CH4 produces weird results(i.e. NAs in any of total estimate, ebullition or diffussion as well as Cero flux in total estimate). It only considers incubations where water_depth > 0

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
files.sources = list.files(path = paste0(repo_root,"/functions"), full.names = T)
for (f in files.sources){source(f)}


# ---- Directories and data loading ----
dropbox_root <- "C:/Users/Miguel/Dropbox/RESTORE4Cs - Fieldwork/Data" # You have to make sure this is pointing to the write folder on your local machine
results_path <- paste0(dropbox_root,"/GHG/Processed data/computed_flux/level_incubation/")
fieldsheetpath <- paste0(dropbox_root,"/GHG/Fieldsheets")
plots_path <- paste0(dropbox_root,"/GHG/Processed data/computed_flux/plots")

diagnostic_path<- paste0(dropbox_root, "/GHG/Working data/EbullitionVSdiffusion diagnosis/")

#Load csv files produced by raw2flux.R. script
# setwd(results_path)

listf <- list.files(path = results_path, pattern = "ch4_fluxes", all.files = T, full.names = T, recursive = F)
table_results_all <- NULL

for (f in listf){
  table_results_all <- rbind(table_results_all, 
                             read.csv(file = f, header = T))}


#Add fieldsheet metadata 
#Read all GHGfieldsheets (using function from Camille)
fieldsheets<- list.files(path= fieldsheetpath, recursive = T, pattern = "GHG.xlsx",full.names = T)
field<- read_GHG_fieldsheets(fieldsheets)

#Drop unwanted variables and rename key UniqueID:
field<- field %>% 
  select(-c(person_sampling,gas_analyzer,logger_floating_chamber,logger_transparent_chamber,start_time,initial_co2,initial_ch4,end_time,final_ch4,final_co2,chamber_height_cm, unix_start,unix_stop)) %>% 
  rename(UniqueID=uniqID)

table_results_meta<- table_results_all %>% 
  merge.data.frame(field, by="UniqueID",all.x = T)

#Filter for water depth > 0
table_water_results<- table_results_meta %>% filter(water_depth>0)

#NAs in flux:
Nabullition<- table_water_results %>% 
  filter(is.na(ebullition+diffusion+total_estimated)) 

#Cero in total flux:
ceroflux<- table_water_results %>% 
  filter(total_estimated==0)


write.csv(Nabullition, file = paste0(diagnostic_path, "ebullition_NAs.csv"),row.names = F)
write.csv(ceroflux, file = paste0(diagnostic_path, "cerototalflux.csv"), row.names = F)

#Proportion of incubations where method fails: ~3%
(dim(Nabullition)[1]+dim(ceroflux)[1])/dim(table_water_results)[1]*100





