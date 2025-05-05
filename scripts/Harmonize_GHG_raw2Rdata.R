
# ---
# Authors: Miguel Cabrera, Camille Minaudo
# Project: "RESTORE4Cs"
# date: "April 2025"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
# This script imports and harmonizes raw data (ghg timeseries) from Picarro, Los Gatos and Li-COR gas analyzers and saves them as sampling-named Rdata files in RData_path

#Data fixes during import and harmonization: 
#Timezone of Lithuanian licor data is corrected with this script before saving as Rdata. (data from Lithuanian licor used during S4 is set to local Riga/Europe timezone, this script harmonizes to UTC timezone)
#Units for picarro H2O data are corrected, goflux default import function assumes H2O is in mmol/mol, but actual units of our instrument are % (i.e. 100* mol/mol). We correct by multiplying the concentration by 10 before saving to Rdata.


#Last run date CHECK!-----
#The Rdata files in dropbox were last updated with this script on: DATE OF LAST UPDATE
#still not runned, rdata in dropbox as of 5/5/2025 has correct tz but picarro H2O units are not correct. 



rm(list = ls()) # clear workspace
cat("/014") # clear console


# ---- Directories ----
dropbox_root <- "C:/Users/Miguel/Dropbox/RESTORE4Cs - Fieldwork/Data" # You have to make sure this is pointing to the right folder on your local machine

datapath <- paste0(dropbox_root,"/GHG/RAW data")
RData_path <- paste0(dropbox_root,"/GHG/Processed data/RData/")


# --- install goFlux if needed ---
#library(devtools)
#install_github("Qepanna/goFlux")

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
require(msm)
require(data.table)
require(tools)
require(pbapply)

repo_root <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
files.sources = list.files(path = paste0(repo_root,"/functions"), full.names = T)
for (f in files.sources){source(f)}


#Import loop------


#Import, correct and store measurements to RData ---
{
  data_folders <- list.dirs(datapath, full.names = T, recursive = T)[-1]
  i <- grep(pattern = "S", x = data_folders) # selecting the files corresponding to the selected sampling campaign
  data_folders <- data_folders[i]
  
  r <- grep(pattern = "RData",x=data_folders)
  if(length(r)>0){data_folders <- data_folders[-r]}
  
  message("Here is the list of data folders in here:")
  print(data_folders)
  
  
  # Import and store data for for Picarro and LosGatos data
  raw2RData_P_LG <- function(data_folders, instrument, instrumentID, date.format, prec){
    r <- grep(pattern = instrument, x=data_folders)
    for (data_folder in data_folders[r]){
      setwd(data_folder)
      subsite = basename(data_folder)
      message(paste0("processing folder ",basename(data_folder)))
      import2RData(path = data_folder, instrument = instrumentID,
                   date.format = date.format, timezone = 'UTC', keep_all = FALSE,
                   prec = prec)
      
      # load all these R.Data into a single dataframe
      file_list <- list.files(path = paste(data_folder,"/RData",sep=""), full.names = T)
      z <- grep(pattern = instrumentID, x=file_list)
      file_list <- file_list[z]
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
      
      #Correction of units for H2O in case of picarro data: default function "import2RData" assumes mmol/mol, but our instrument reports H2O in mol % (100* mol/mol)
      if(instrument=="Picarro"){
       mydata$H2O_ppm<- mydata$H2O_ppm*10 
      }
      
      # save this dataframe as a new RData file
      setwd(RData_path)
      save(mydata, file = paste0(subsite,"_",instrument,".RData"))
    }
  }
  
  raw2RData_P_LG(data_folders, instrument = "Picarro", instrumentID = "G4301", date.format = "ymd", prec=c(0.025, 0.1, 10))
  raw2RData_P_LG(data_folders, instrument = "Los Gatos", instrumentID = "UGGA", date.format = "mdy", prec =  c(0.2, 1.4, 50))
  
  
  
  # Import and store data for LiCOR data
  fieldsheet_Licor <- fieldsheet[fieldsheet$gas_analyzer=="LI-COR",]
  list_subsites_Licor <- unique(fieldsheet_Licor$subsite)
  
  
  
  r_licor <- grep(pattern = "Licor",x=data_folders)
  for (data_folder in data_folders[r_licor]){
    setwd(data_folder)
    message(paste0("processing folder ",basename(data_folder)))
    if (grepl("S4-CU",data_folder)){
      import2RData(path = data_folder, instrument = "LI-7810",
                   date.format = "ymd", timezone = 'Europe/Riga', keep_all = FALSE, 
                   prec = c(3.5, 0.6, 45))
    } else {
      import2RData(path = data_folder, instrument = "LI-7810",
                   date.format = "ymd", timezone = 'UTC', keep_all = FALSE, 
                   prec = c(3.5, 0.6, 45))
    }
    
    
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
    
    #Correct lithuanian licor tz (licor in local time)
    if (grepl("S4-CU",data_folder)){
      message("S4-CU Licor not in UTC time!! Now corrected to UTC")  
      mydata_imp$POSIX.time<- with_tz(mydata_imp$POSIX.time,tzone = "UTC")
      mydata_imp$DATE<- as.character(as_date(mydata_imp$POSIX.time))
      mydata_imp$TIME<- strftime(mydata_imp$POSIX.time, format="%H:%M:%S",tz = "utc")
    }
    
    # get rid of possible duplicated data
    is_duplicate <- duplicated(mydata_imp$POSIX.time)
    mydata_imp <- mydata_imp[!is_duplicate,]
    
    
    # r_site <- grep(pattern = basename(data_folder), x=list_subsites_Licor)
    
    # create separate folders for each subsite where data actually exists
    for (subsite in list_subsites_Licor){
      corresponding_date <- as.Date(unique(fieldsheet_Licor$date[fieldsheet_Licor$subsite == subsite]))
      
      # check if we already have this data somewhere in the clean file
      n_lines_d <- dim(mydata_imp[which(as.Date(mydata_imp$POSIX.time) == corresponding_date),])[1]
      if(n_lines_d > 0){
        message(paste0("... there is some data to store for subsite ",subsite))
        mydata <- mydata_imp[as.Date(mydata_imp$POSIX.time) == corresponding_date,]
        
        setwd(RData_path)
        # save this dataframe as a new RData file
        save(mydata, file = paste0(subsite,"_","LI-7810",".RData"))
      }
    }
  }
}
