
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



# ---- Directories ----
dropbox_root <- "C:/Users/Camille Minaudo/Dropbox/RESTORE4Cs - Fieldwork/Data"
datapath <- paste0(dropbox_root,"/GHG/RAW data/RAW Data Picarro")
fieldsheetpath <- paste0(dropbox_root,"/GHG/Fieldsheets")

# ---- List GHG chamber fieldsheets in Dropbox ----
myfieldsheets_list <- list.files(fieldsheetpath, pattern = "Fieldsheet-GHG.xlsx", all.files = T, full.names = T, recursive = T)


# ---- Read all fieldsheets and put them in a single dataframe ----

# Read the first one to get the headers
fieldsheet_temp <- readxl::read_xlsx(myfieldsheets_list[1],
                                     col_names = T)
my_headers <- names(fieldsheet_temp)

# Go through them all and keep the info
isF <- T
for (f in myfieldsheets_list){
  fieldsheet_temp <- readxl::read_xlsx(f, col_names = F, range = "A3:V30")
  names(fieldsheet_temp) <- my_headers
  fieldsheet_temp <- fieldsheet_temp[!is.na(fieldsheet_temp$plot_id),]
  fieldsheet_temp$date <- as.Date(fieldsheet_temp$date, tryFormats = c("%d.%m.%Y", "%d/%m/%Y"))
  fieldsheet_temp$subsiteID <- gsub(pattern = "-Fieldsheet-GHG.xlsx", replacement = "", x = basename(f))


  if(isF){
    isF <- F
    fieldsheet <- fieldsheet_temp
  } else {
    fieldsheet <- rbind(fieldsheet, fieldsheet_temp)
  }
}


fieldsheet <- fieldsheet[!is.na(fieldsheet$longitude),]

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
  mydata_imp_clean <- mydata_imp[!is_duplicate,]

  # save this dataframe as a new RData file
  setwd(mainDirRData)
  save(mydata_imp_clean, file = paste0(basename(data_folder),"_data_clean.RData"))
}







