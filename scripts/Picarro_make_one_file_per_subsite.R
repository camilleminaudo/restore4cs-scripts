
# ---
# Authors: Camille Minaudo
# Project: "RESTORE4Cs"
# date: "Nov 2023"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
# This script loads all Picarro raw measurements (data level L0a) from a given folder
# and produces RData archives (data level L0b), stored with a dedicated file per
# campaign-subsite (e.g. S1-CA-P1)


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
require(pbapply)


source(paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/get_unix_times.R"))
source(paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/G2508_import.R"))
source(paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/import2RData.R"))
source(paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/read_GHG_fieldsheets.R"))



# ---- Directories ----
dropbox_root <- "C:/Users/Camille Minaudo/Dropbox/RESTORE4Cs - Fieldwork/Data"
datapath <- paste0(dropbox_root,"/GHG/RAW data/RAW Data Picarro")
fieldsheetpath <- paste0(dropbox_root,"/GHG/Fieldsheets")


# ---- List GHG chamber fieldsheets in Dropbox and read them ---
# list filenames
myfieldsheets_list <- list.files(fieldsheetpath, pattern = "Fieldsheet-GHG.xlsx", all.files = T, full.names = T, recursive = T)
# Read all fieldsheets and put them in a single dataframe
fieldsheet <- read_GHG_fieldsheets(myfieldsheets_list)

fieldsheet_picarro <- fieldsheet[fieldsheet$gas_analyzer=="Picarro",]


# ---- Import and store measurements to RData ----
data_folders <- list.dirs(datapath, full.names = T, recursive = F)
r <- grep(pattern = "RData",x=data_folders)
data_folders <- data_folders[-r]

message("Here is the list of data folders in here:")
print(data_folders)

map_analysers <- data.frame(instrument = c("G2508"),
                            date.format = c("ymd"))

mainDirRData = paste(datapath,"/RData",sep="")

for (data_folder in data_folders){
  setwd(data_folder)
  subsite = basename(data_folder)
  message(paste0("processing folder ",basename(data_folder)))
  import2RData(path = data_folder, instrument = map_analysers$instrument,
               date.format = map_analysers$date.format, timezone = 'UTC')

  # load all these R.Data into a single dataframe
  file_list <- list.files(path = paste(data_folder,"/RData",sep=""), full.names = T)
  isF <- T
  for(i in seq_along(file_list)){
    load(file_list[i])
    if(isF){
      isF <- F
      mydata_imp <- data.raw
    } else {
      mydata_imp <- rbind(mydata_imp, data.raw)
    }
    rm(data.raw)
  }

  # get read of possible duplicated data
  is_duplicate <- duplicated(mydata_imp$POSIX.time)
  mydata <- mydata_imp[!is_duplicate,]

  # creating a folder where to put this data
  dir.create(file.path(mainDirRData, subsite))
  setwd(file.path(mainDirRData, subsite))

  # save this dataframe as a new RData file
  save(mydata, file = paste0("data_",subsite,".RData"))
}







