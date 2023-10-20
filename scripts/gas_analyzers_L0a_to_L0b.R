
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
loggerspath <- paste0(datapath,"/RAW Data Logger")

setwd(datapath)

# --- SETTINGS
analyser <- "Licor"
site_ID <- "S1-CU"
subsite_ID <- "S1-CU-P2"
file_to_read <- "S1-CU-P2.data"
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
path2file <- paste0(fieldsheetpath,"/",site_ID,"/",subsite_ID,"-Fieldsheet-GHG.xlsx")
fieldsheet_temp <- readxl::read_xlsx(path2file,
                                     col_names = T)
fieldsheet <- readxl::read_xlsx(path2file,
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


# Read corresponding Loggers data
SN_logger_float <- unique(fieldsheet$logger_floating_chamber)
SN_logger_tube <- unique(fieldsheet$logger_transparent_chamber)

if(SN_logger_float != "NA"){
  is_data_logger_float = T
  data_logger_float <- readxl::read_xlsx(paste0(loggerspath,"/",site_ID,"-",SN_logger_float,".xlsx"),col_names = T)
  data_logger_float$unixtime <- as.numeric(data_logger_float$`Date/hour (UTC)`)
} else {
  message("no data logger linked to the floating chamber!")
  is_data_logger_float = F}

if(SN_logger_tube != "NA"){
  is_data_logger_tube = T
  data_logger_tube <- readxl::read_xlsx(paste0(loggerspath,"/",site_ID,"-",SN_logger_tube,".xlsx"),col_names = T)
  data_logger_tube$unixtime <- as.numeric(data_logger_tube$`Date/hour (UTC)`)
} else {
  is_data_logger_tube = F
  message("no data logger linked to the transparent tube chamber!")
  if(SN_logger_float != "NA"){
    message(".../ using data from floating chamber instead...")
    data_logger_tube = data_logger_float
    is_data_logger_tube = T
  }
}


for (i in  seq(1,length(fieldsheet$pilot_site))){ # for each incubation, proceed with...

  my_sel <- my_data[my_data$unixtime>= (fieldsheet$unix_start_time[i]-30) & my_data$unixtime<= (fieldsheet$unix_end_time[i]+60),]

  my_label <- unique(my_sel$label)
  my_label <- my_label[which(my_label != "")]

  pCO2 <- ggplot(my_sel, aes(unixtime, CO2))+geom_line()+
    geom_line(data = my_data[my_data$label == my_label,],
              aes(unixtime, CO2), colour = "red", alpha = 0.5, linewidth = 1.5)+
    ylab("CO2 [ppm]")+
    ggtitle(paste0(subsite_ID," ",fieldsheet$date[i],", plot ",fieldsheet$plot_id[i],", incub ", i))+
    theme_article()
  pCH4 <- ggplot(my_sel, aes(unixtime, CH4))+geom_line()+
    geom_line(data = my_data[my_data$label == my_label,],
              aes(unixtime, CH4), colour = "red", alpha = 0.5, linewidth = 1.5)+
    ylab("CH4 [ppm]")+
    ggtitle(paste(fieldsheet$chamber_type[i], fieldsheet$strata[i], fieldsheet$transparent_dark[i], sep=", "))+
    theme_article()

  p <- ggarrange(pCO2,pCH4, nrow = 1)


  my_sel <- my_data[my_data$label == my_label,]

  if (fieldsheet$chamber_type[i] == "floating"){
    my_sel$temperature <- approx(data_logger_float$unixtime, data_logger_float$`Ch:1 - Temperature   (°C)`, xout = my_sel$unixtime)$y
  } else if (fieldsheet$chamber_type[i] == "tunbe"){
    my_sel$temperature <- approx(data_logger_tube$unixtime, data_logger_tube$`Ch:1 - Temperature   (°C)`, xout = my_sel$unixtime)$y
  } else {
    warning("chamber type not correct!")
  }


  head(my_sel)

}

