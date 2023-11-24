
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

# ---- Directories ----
dropbox_root <- "C:/Users/Camille Minaudo/Dropbox/RESTORE4Cs - Fieldwork/Data"
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
  fieldsheet_temp$date <- as.Date( fieldsheet_temp$date, tryFormats = c("%d.%m.%Y", "%d/%m/%Y"))
  fieldsheet_temp$subsite <- gsub(pattern = "-Fieldsheet-GHG.xlsx",replacement = "",x = basename(f))
  fieldsheet_temp$water_depth <- as.numeric(fieldsheet_temp$water_depth)

  if(isF){
    isF <- F
    fieldsheet <- fieldsheet_temp
  } else {
    fieldsheet <- rbind(fieldsheet, fieldsheet_temp)
  }
}

fieldsheet <- fieldsheet[!is.na(fieldsheet$longitude),]

fieldsheet$unix_start <- get_unix_times(mydate = fieldsheet$date, mytime = fieldsheet$start_time)
fieldsheet$unix_stop <- get_unix_times(mydate = fieldsheet$date, mytime = fieldsheet$end_time)

fieldsheet$start_time <- strftime(fieldsheet$start_time, format="%H:%M:%S", tz = 'utc')
fieldsheet$end_time <- strftime(fieldsheet$end_time, format="%H:%M:%S", tz = 'utc')

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



