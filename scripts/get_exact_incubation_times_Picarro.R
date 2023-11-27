# ---
# Authors: Camille Minaudo
# Project: "RESTORE4Cs"
# date: "Nov 2023"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
# This script extracts the exact datetime for each chamber measurements performed
# with the Picarro Gas Scouter. It loads the corresponding fieldworksheet and
# produce a new fieldwork sheet with corrected time stamps

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
library(GoFluxYourself)
require(dplyr)
require(purrr)

# ---- functions ----


source(paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/get_unix_times.R"))


read_map_incubations <- function(path2folder){

  my_maps_filenames <- list.files(data_folder, pattern = ".csv", all.files = T, full.names = T, recursive = F)
  isF <- T
  for(my_maps_filename in my_maps_filenames){

    map_incubations_temp <- read.csv(file = my_maps_filename,
                                     header = T)
    map_incubations_temp <- map_incubations_temp[map_incubations_temp$Species ==  "CH4",]

    if(isF){
      isF <- F
      map_incubations <- map_incubations_temp
    } else {
      map_incubations <- rbind(map_incubations, map_incubations_temp)
    }
  }

  map_incubations <- data.frame(subsite = basename(path2folder),
                                plot = map_incubations$Comment,
                                time_code = map_incubations$Time.Code,
                                start = map_incubations$start_fit,
                                stop = map_incubations$end_fit)

  return(map_incubations)
}




# ---- Directories ----
dropbox_root <- "C:/Users/Camille Minaudo/Dropbox/RESTORE4Cs - Fieldwork/Data"
datapath <- paste0(dropbox_root,"/GHG/RAW data/RAW Data Picarro")
fieldsheetpath <- paste0(dropbox_root,"/GHG/Fieldsheets")
path2processed <- paste0(dropbox_root,"/GHG/Processed data")

# ---- List folders with Picarro data in Dropbox ----
data_folders <- list.dirs(datapath, full.names = T, recursive = F)[-1]
subsites <- basename(data_folders)

k=0
for (subsite in subsites){
  message(paste0("processing data for ",subsite))
  k=k+1
  data_folder <- data_folders[k]

  # read all the csv files in data_folder, and group into a single one
  map_incubations <- read_map_incubations(path2folder = data_folder)

  # load corresponding fieldsheet
  pilotsite <- substr(subsite, 1, 5)
  f_name <- paste0(fieldsheetpath,"/",pilotsite,"/",subsite,"-Fieldsheet-GHG.xlsx")

  fieldsheet <- readxl::read_xlsx(f_name, col_names = T)# Read it first to get the headers
  my_headers <- names(fieldsheet)
  fieldsheet <- readxl::read_xlsx(f_name, col_names = F, range = "A3:V50", n_max = 100)# re-read and affect correctly the headers
  names(fieldsheet) <- my_headers
  fieldsheet <- fieldsheet[!is.na(fieldsheet$plot_id),]
  fieldsheet$date <- as.Date( fieldsheet$date, tryFormats = c("%d.%m.%Y", "%d/%m/%Y"))


  fieldsheet$unix_start <- get_unix_times(mydate = fieldsheet$date, mytime = fieldsheet$start_time)
  fieldsheet$unix_stop <- get_unix_times(mydate = fieldsheet$date, mytime = fieldsheet$end_time)


  # check the closest incubation in map_incubations for each row in fieldsheet.
  # if more than 3 minutes apart, we consider the row in map_incubations out of sampling

  corresponding_row <- NA*fieldsheet$plot_id

  for (i in seq_along(fieldsheet$plot_id)){
    ind <- which.min(abs(fieldsheet$unix_start[i] - map_incubations$start))
    if(abs(fieldsheet$unix_start[i] - map_incubations$start[ind])<3*60){corresponding_row[i] <- ind}
  }
  corresponding_row <- corresponding_row[!is.na(corresponding_row)]
  map_incubations <- map_incubations[corresponding_row, ]

  if(dim(map_incubations)[1] != dim(fieldsheet)[1]){
    fieldsheet <- fieldsheet[seq_along(corresponding_row),]
  }

  fieldsheet$start_time <- strftime(fieldsheet$start_time, format="%H:%M:%S", tz = 'utc')
  fieldsheet$end_time <- strftime(fieldsheet$end_time, format="%H:%M:%S", tz = 'utc')

  fieldsheet$unix_start_corrected <- map_incubations$start
  fieldsheet$unix_stop_corrected <- map_incubations$stop

  fieldsheet$timestamp_start <- as.POSIXct(fieldsheet$unix_start_corrected, tz = "UTC", origin = "1970-01-01")
  fieldsheet$timestamp_stop <- as.POSIXct(fieldsheet$unix_stop_corrected, tz = "UTC", origin = "1970-01-01")

  fieldsheet$start_time_corrected <- strftime(fieldsheet$timestamp_start, format="%H:%M:%S", tz = 'utc')
  fieldsheet$end_time_corrected <- strftime(fieldsheet$timestamp_stop, format="%H:%M:%S", tz = 'utc')

  setwd(paste0(path2processed,"/corrected_fieldsheets"))
  f_name <- paste0(subsite,"-Fieldsheet-GHG_corrected.csv")
  write.csv(file = f_name, x = fieldsheet, row.names = F)

}





