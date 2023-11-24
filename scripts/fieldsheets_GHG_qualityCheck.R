
# ---
# Authors: Camille Minaudo
# Project: "RESTORE4Cs"
# date: "Oct 2023"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
# This script loads all GHG chamber measurement fieldsheets and runs a simple quality check:
# ==> makes sure that end time is always > start time
#


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

source(paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/get_unix_times.R"))
source(paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/read_GHG_fieldsheets.R"))

# ---- Directories ----
dropbox_root <- "C:/Users/Camille Minaudo/Dropbox/RESTORE4Cs - Fieldwork/Data"
fieldsheetpath <- paste0(dropbox_root,"/GHG/Fieldsheets")

# ---- List GHG chamber fieldsheets in Dropbox and read them ---
# list filenames
myfieldsheets_list <- list.files(fieldsheetpath, pattern = "Fieldsheet-GHG.xlsx", all.files = T, full.names = T, recursive = T)
# Read all fieldsheets and put them in a single dataframe
fieldsheet <- read_GHG_fieldsheets(myfieldsheets_list)


# --------- rows with stop time < start time
ind_erronous_times <- which(fieldsheet$unix_stop < fieldsheet$unix_start)

message("the following rows show erronous start/stop times")
as.data.frame(fieldsheet[ind_erronous_times,c("pilot_site","subsite","start_time","end_time")])


# --------- rows with possible error with CH4 units (ppb instead of ppm)
ind_suspicious_ch4 <- which(fieldsheet$final_ch4 > 2000)

message("the following rows show suspicious methane levels")
as.data.frame(fieldsheet[ind_suspicious_ch4,c("pilot_site","subsite","plot_id","final_ch4")])



# --------- rows with depth possibly reported in m instead of cm
fieldsheet_not_dry <- fieldsheet[fieldsheet$water_depth != 0,]
ind_suspicious_depth <- which(fieldsheet_not_dry$water_depth < 2)

message("the following rows show suspicious water depth ranges:")
as.data.frame(fieldsheet_not_dry[ind_suspicious_depth,c("pilot_site","subsite","plot_id","water_depth")])
message("which corresponds to the following subsites:")
unique(fieldsheet_not_dry$subsite[ind_suspicious_depth])



