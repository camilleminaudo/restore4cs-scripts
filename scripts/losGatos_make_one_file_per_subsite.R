
# ---
# Authors: Camille Minaudo
# Project: "RESTORE4Cs"
# date: "Nov 2023"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
# This script loads all Los Gatos raw measurements (data level L0a) from a given folder
# and produces RData archives (data level L0b), stored with a dedicated folder per
# by campaign-subsite


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
library(goFlux)
require(dplyr)
require(purrr)
require(pbapply)


repo_root <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
files.sources = list.files(path = paste0(repo_root,"/functions"), full.names = T)
for (f in files.sources){source(f)}



# ---- Directories ----
dropbox_root <- "C:/Users/Camille Minaudo/Dropbox/RESTORE4Cs - Fieldwork/Data"
datapath <- paste0(dropbox_root,"/GHG/RAW data/RAW Data Los Gatos")
fieldsheetpath <- paste0(dropbox_root,"/GHG/Fieldsheets")

# ---- List GHG chamber fieldsheets in Dropbox and read them ---
# list filenames
myfieldsheets_list <- list.files(fieldsheetpath, pattern = "Fieldsheet-GHG.xlsx", all.files = T, full.names = T, recursive = T)

#################################
sampling <- "S3"
#################################

i <- grep(pattern = sampling, x = myfieldsheets_list) # selecting the files corresponding to the selected sampling campaign
myfieldsheets_list <- myfieldsheets_list[i]

# Read all fieldsheets and put them in a single dataframe
fieldsheet <- read_GHG_fieldsheets(myfieldsheets_list)


fieldsheet_losgatos <- fieldsheet[fieldsheet$gas_analyzer=="Los Gatos",]


# ---- Import and store measurements to RData ----
data_folders <- list.dirs(datapath, full.names = T, recursive = F)
r <- grep(pattern = "RData",x=data_folders)
data_folders <- data_folders[-r]

message("Here is the list of data folders in here:")
print(data_folders)

map_analysers <- data.frame(instrument = c("LGR"),
                            date.format = c("mdy"))

mainDirRData = paste(datapath,"/RData",sep="")

for (data_folder in data_folders){
  setwd(data_folder)
  subsite = basename(data_folder)
  message(paste0("processing folder ",basename(data_folder)))
  import2RData(path = data_folder, instrument = map_analysers$instrument,
               date.format = map_analysers$date.format, timezone = 'UTC')

  # load all these R.Data into a single dataframe
  file_list <- list.files(path = paste(data_folder,"/RData",sep=""), full.names = T)
  if(length(file_list)>0){
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
    dir.create(file.path(mainDirRData,subsite))
    setwd(file.path(mainDirRData, subsite))

    # save this dataframe as a new RData file
    save(mydata, file = paste0("data_",subsite,".RData"))
  }
}







