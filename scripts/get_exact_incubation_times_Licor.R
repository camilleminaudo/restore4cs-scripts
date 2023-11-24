# ---
# Authors: Camille Minaudo
# Project: "RESTORE4Cs"
# date: "Nov 2023"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
# This script extracts the exact datetime for each chamber measurements performed
# with the LiCor. It loads the corresponding fieldworksheet and
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
source(paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/read_Licor.R"))



# ---- Directories ----
dropbox_root <- "C:/Users/Camille Minaudo/Dropbox/RESTORE4Cs - Fieldwork/Data"
datapathRAW <- paste0(dropbox_root,"/GHG/RAW data/RAW Data Licor-7810")
datapath <- paste0(dropbox_root,"/GHG/RAW data/RAW Data Licor-7810/RData")
fieldsheetpath <- paste0(dropbox_root,"/GHG/Fieldsheets")
path2processed <- paste0(dropbox_root,"/GHG/Processed data")

# ---- List folders with Licor data in Dropbox ----
data_folders <- list.dirs(datapath, full.names = T, recursive = F)
subsites <- basename(data_folders)
print(subsites)

# go through the RAW data
fs <- list.files(path = datapathRAW, pattern = c(".txt", ".data"), full.names = T, recursive = T)
r <- grep(pattern = ".RData",x=fs)
fs <- fs[-r]


isF = T
for (f in fs){
  message(paste0("extracting info from ",basename(f)))
  data <- read_Licor(file = f)
  labels <- unique(data$label)
  labels<- labels[-which(labels=="")]

  # get start and stop times for each label
  for (label in labels){
    start <- min(data$unixtime[data$label == label])
    stop <- max(data$unixtime[data$label == label])

    map_incubations_temp <- data.frame(date = first(data$date[data$label == label]),
                                       label = label,
                                       start = start,
                                       stop = stop,
                                       time_start = strftime(start, format="%H:%M:%S", tz = 'utc'),
                                       time_stop = strftime(stop, format="%H:%M:%S", tz = 'utc'),
                                       folder = basename(dirname(f)),
                                       file = basename(f))
    if(isF){
      isF = F
      map_incubations <- map_incubations_temp
    } else {
      map_incubations <- rbind(map_incubations,map_incubations_temp)
    }
  }
}


# uniqueID <- paste0(map_incubations$date, map_incubations$label)
map_incubations <- map_incubations[!duplicated(map_incubations$start),]
map_incubations <- map_incubations[order(map_incubations$start),]

# Save incubtion map to Dropbox
write.csv(x = map_incubations, file = paste0(datapathRAW,"/map_incubations.csv"), row.names = F)


# k=0
for (subsite in subsites){
  message(paste0("processing data for ",subsite))
  # k=k+1
  # data_folder <- data_folders[k]

  # load corresponding fieldsheet
  pilotsite <- substr(subsite, 1, 5)
  f_name <- paste0(fieldsheetpath,"/",pilotsite,"/",subsite,"-Fieldsheet-GHG.xlsx")
  message(f_name)

  fieldsheet <- readxl::read_xlsx(f_name, col_names = T)# Read it first to get the headers
  my_headers <- names(fieldsheet)
  fieldsheet <- readxl::read_xlsx(f_name, col_names = F, range = "A3:V50")# re-read and affect correctly the headers
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
    # corresponding_row[i] <- ind
    if(abs(fieldsheet$unix_start[i] - map_incubations$start[ind])<3*60){corresponding_row[i] <- ind}
  }
  isNA <- !is.na(corresponding_row)
  corresponding_row <- corresponding_row[!is.na(corresponding_row)]
  my_map_incubations <- map_incubations[corresponding_row, ]

  if(dim(my_map_incubations)[1] >0){
    message("....found corresponding labels for this one")
    fieldsheet <- fieldsheet[isNA,]

    fieldsheet$start_time <- strftime(fieldsheet$start_time, format="%H:%M:%S", tz = 'utc')
    fieldsheet$end_time <- strftime(fieldsheet$end_time, format="%H:%M:%S", tz = 'utc')

    fieldsheet$unix_start_corrected <- my_map_incubations$start
    fieldsheet$unix_stop_corrected <- my_map_incubations$stop

    fieldsheet$timestamp_start <- as.POSIXct(fieldsheet$unix_start_corrected, tz = "UTC", origin = "1970-01-01")
    fieldsheet$timestamp_stop <- as.POSIXct(fieldsheet$unix_stop_corrected, tz = "UTC", origin = "1970-01-01")

    fieldsheet$start_time_corrected <- strftime(fieldsheet$timestamp_start, format="%H:%M:%S", tz = 'utc')
    fieldsheet$end_time_corrected <- strftime(fieldsheet$timestamp_stop, format="%H:%M:%S", tz = 'utc')

    setwd(paste0(path2processed,"/corrected_fieldsheets"))
    f_name <- paste0(subsite,"-Fieldsheet-GHG_corrected.csv")
    write.csv(file = f_name, x = fieldsheet, row.names = F)

  }


}





