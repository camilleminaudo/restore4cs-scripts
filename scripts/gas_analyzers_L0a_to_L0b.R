
# ---
# Authors: Camille Minaudo, Benjamin Misteli
# Project: "RESTORE4Cs"
# date: "Oct 2023"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
# This script loads a raw measurement file (data level L0a) from one of the gas 
# analyzers used in the project, and transform it into a unified harmonized csv 
# file (data level L0b) allowing for visualization of the data and further data processing.


rm(list = ls()) # clear workspace
cat("/014") # clear console


# packages
library(tidyverse)
library(readxl)
library(lubridate)
library(zoo)
library(ggplot2)
library(grid)
library(egg)


# Directories
dropbox_root <- "C:/Users/Camille Minaudo/Dropbox/RESTORE4Cs - Fieldwork/Data"
datapath <- paste0(dropbox_root,"/GHG/RAW data")
fieldsheetpath <- paste0(dropbox_root,"/GHG/Fieldsheets")

setwd(datapath)

# --- SETTINGS
analyser <- "Licor"
subsite_ID <- "S1-CU-A2"
file_to_read <- "S1-CU-A2.data"
who_runs_this <- "Camille Minaudo"

# --- Some functions
read_Licor <- function(file){
  message(paste0("reading ",file_to_read," with script for ",analyser," gas analyser"))
  filename <- paste(datapath,directory_analyser,file_to_read, sep = "/")
  data_raw <- read_lines(filename)
  prefix <- substr(data_raw, start = 1,stop = 5) # isolte first 5 characters for each line
  
  # find line corresponding to headers
  headers <- unlist(strsplit(data_raw[which(prefix == "DATAH")], "\t"))
  units <- unlist(strsplit(data_raw[which(prefix == "DATAU")], "\t"))
  
  data <- read.delim(filename, sep = "\t", header = F, skip = which(prefix == "DATAU"), na.strings = "nan")
  names(data) <- headers
  
  my_data <- data.frame(date = data$DATE,
                        UTCtime = data$TIME,
                        unixtime = data$SECONDS,
                        H2O = data$H2O,
                        CO2 = data$CO2,
                        CH4 = data$CH4/1000, #ppm
                        Press = data$CAVITY_P,
                        label = data$REMARK)
  
  return(my_data)
}


# Read gas analyser's file
if(analyser == "Licor"){
  directory_analyser <- "RAW Data Licor-7810"
  my_data <- read_Licor(file = paste(datapath,directory_analyser,file_to_read, sep = "/"))
} else if(analyser == "Los Gatos"){
  directory_analyser <- "RAW Data Los Gatos"
} else if (analyser == "Picarro"){
  directory_analyser <- "RAW Data Picarro"
}


# select only rows with no CO2 NA values
# my_data <- my_data[!is.na(my_data$CO2),]


# Read corresponding Fieldsheet
fieldsheet_temp <- readxl::read_xlsx(paste0(fieldsheetpath,"/",subsite_ID,"-Fieldsheet-GHG.xlsx"), 
                                     col_names = T)
fieldsheet <- readxl::read_xlsx(paste0(fieldsheetpath,"/",subsite_ID,"-Fieldsheet-GHG.xlsx"), 
                                skip = 2, col_names = F)
names(fieldsheet) <- names(fieldsheet_temp)


get_unix_times <- function(mydate, mytime){
  timestamp <- paste(hour(mytime),minute(mytime),0,sep = ":")
  unix_time <- as.numeric(as.POSIXct(paste(as.Date(mydate),timestamp, 
                                           sep = " "), tz = 'UTC'))
  return(unix_time)
}

fieldsheet$unix_start_time <- get_unix_times(mydate = fieldsheet$date, mytime = fieldsheet$start_time)
fieldsheet$unix_end_time <- get_unix_times(mydate = fieldsheet$date, mytime = fieldsheet$end_time)

# head(fieldsheet)

for (i in  seq(1,length(fieldsheet$pilot_site))){ # for each incubation, proceed with...
  
  my_sel <- my_data[my_data$unixtime>= fieldsheet$unix_start_time[i] & my_data$unixtime<= fieldsheet$unix_end_time[i],]
  
  pCO2 <- ggplot(my_sel, aes(unixtime, CO2))+geom_line()+
    ylab("CO2 [ppm]")+
    ggtitle(paste0(subsite_ID," ",fieldsheet$date[i],", plot ",fieldsheet$plot_id[i],", incub ", i))+
    theme_article()
  pCH4 <- ggplot(my_sel, aes(unixtime, CH4))+geom_line()+
    ylab("CH4 [ppm]")+
    ggtitle(paste(fieldsheet$chamber_type[i], fieldsheet$strata[i], fieldsheet$transparent_dark[i], sep=", "))+
    theme_article()
  
  p <- ggarrange(pCO2,pCH4, nrow = 1)
  
}

