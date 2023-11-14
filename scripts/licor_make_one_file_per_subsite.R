
# ---
# Authors: Camille Minaudo
# Project: "RESTORE4Cs"
# date: "Oct 2023"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
# This script loads all raw measurement files (data level L0a) from one of the gas
# analyzers used in the project, transform it into a unified harmonized csv
# file, and saves a separate file + a pdf plot out of each incubation performed.


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


source(paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/click.peak.R"))
source(paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/click.peak.loop.R"))
source(paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/flux.term.R"))
source(paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/read_Licor.R"))
source(paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/get_unix_times.R"))
source(paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/G2508_import.R"))
source(paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/import2RData.R"))



# ---- Directories ----
dropbox_root <- "C:/Users/Camille Minaudo/Dropbox/RESTORE4Cs - Fieldwork/Data"
datapath <- paste0(dropbox_root,"/GHG/RAW data")
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
  fieldsheet_temp$date <- as.Date( fieldsheet_temp$date, tryFormats = c("%d.%m.%y", "%d/%m/%y"))

  if(isF){
    isF <- F
    fieldsheet <- fieldsheet_temp
  } else {
    fieldsheet <- rbind(fieldsheet, fieldsheet_temp)
  }
}


fieldsheet <- fieldsheet[!is.na(fieldsheet$longitude),]




# ---- Import and store measurements to RData ----
data_folders <- list.dirs(datapath, full.names = T, recursive = F)[-2]
message("Here is the list of data folders in here:")
print(data_folders)


i_licor <- grep(pattern = "Licor", x = data_folders)
i_losgatos <- grep(pattern = "Los Gatos", x = data_folders)
i_picarro <- grep(pattern = "Picarro", x = data_folders)
map_analysers <- data.frame(instrument = c("LI-7810","LGR","G2508"),
                            i = c(i_licor, i_losgatos, i_picarro),
                            date.format = c("ymd","mdy","ymd"))
k=0
for (data_folder in data_folders){
  k=k+1
  setwd(data_folder)
  import2RData(path = data_folder, instrument = map_analysers$instrument[k],
               date.format = map_analysers$date.format[k], timezone = 'UTC')

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
  mydata_imp_clean <- mydata_imp %>%
    group_by(DATE,TIME) %>% # the complete group of interest
    mutate(duplicate = n()) %>% # count number in each group
    filter(duplicate == 1) %>% # select only unique records
    select(-duplicate) # remove group count column


  # save this dataframe as a new RData file
  setwd(data_folder)
  save(mydata_imp_clean, file = "data_all_clean.RData")

}





