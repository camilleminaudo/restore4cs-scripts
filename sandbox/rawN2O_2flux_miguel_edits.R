

# ---
# Authors: MIGUEL CABRERA
# Project: "RESTORE4Cs"
# date: "Feb 2025"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
# This script computes fluxes for all incubations with N2O data. 
# Table with N2O concentrations in chambers and atmosphere for every UniqueID_notime is needed.
#Incubation definition: incubations for N2O are defined as the duration calculated directly from fieldsheet times. 

#Auxtable extraction adapted from raw2flux

rm(list = ls()) # clear workspace
cat("/014") # clear console


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
# library(goFlux)
# require(dplyr)
# require(purrr)
require(msm)
# require(data.table)
require(tools)
# require(pbapply)

repo_root <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
files.sources = list.files(path = paste0(repo_root,"/functions"), full.names = T)
for (f in files.sources){source(f)}


# ---- Directories ----
dropbox_root <- "C:/Users/Miguel/Dropbox/RESTORE4Cs - Fieldwork/Data" # You have to make sure this is pointing to the write folder on your local machine
datapath <- paste0(dropbox_root,"/GHG/N2O_fluxes")
fieldsheetpath <- paste0(dropbox_root,"/GHG/Fieldsheets")
loggerspath <- paste0(dropbox_root,"/GHG/RAW data/RAW Data Logger")

results_path <- datapath

sampling<- "S4"


#_________________####


#FIELDSHEET PREP --------


#Data pre-processing and harmonization: 
# ---- List GHG chamber fieldsheets in Dropbox and read them ---
# list filenames
myfieldsheets_list <- list.files(fieldsheetpath, pattern = "Fieldsheet-GHG.xlsx", all.files = T, full.names = T, recursive = T)
i <- grep(pattern = sampling, x = myfieldsheets_list) # selecting the files corresponding to the selected sampling campaign
myfieldsheets_list <- myfieldsheets_list[i]
# Read all fieldsheets and put them in a single dataframe
fieldsheet <- read_GHG_fieldsheets(myfieldsheets_list)


#Create UniqueID_notime column to join
fieldsheet$UniqueID_notime<- fieldsheet %>% 
  mutate(UniqueID_notime=tolower(paste(subsite, plot_id,substr(strata,1,1) ,substr(transparent_dark,1,1), sep = "-"))) %>% pull(UniqueID_notime)
  
  
#Check: UniqueID_notime is unique?
length(unique(fieldsheet$UniqueID_notime))==dim(fieldsheet)[1]
  
# recalculating start and stop in proper formats
fieldsheet$timestamp_start <- as.POSIXct(fieldsheet$unix_start, tz = "UTC", origin = "1970-01-01")
fieldsheet$timestamp_stop <- as.POSIXct(fieldsheet$unix_stop, tz = "UTC", origin = "1970-01-01")

fieldsheet$start_time <- strftime(fieldsheet$timestamp_start, format="%H:%M:%S", tz = 'utc')
fieldsheet$end_time <- strftime(fieldsheet$timestamp_stop, format="%H:%M:%S", tz = 'utc')

fieldsheet <- fieldsheet[order(fieldsheet$subsite),]


#_________________####
#FLUX CALCULATION####

# ----- 1. Auxtable prep -----

# This section loops over each subsite and extracts the necessary chamber parameters for flux calculation: T, P, V, A
subsites <- unique(fieldsheet$subsite)

isF_incub <- T
isFsubsite <- T
auxfile <- NULL
for (subsite in subsites){
  message("Now processing ",subsite)
  
  corresp_fs <- fieldsheet[fieldsheet$subsite == subsite,]
  gs <- first(corresp_fs$gas_analyzer)
  
  ## ----1.1. Import looger data----
  # read corresponding temperature logger file and take median temperature during incubation (in raw2flux is only initial initial temperature)
  SN_logger_float <- first(corresp_fs$logger_floating_chamber)
  SN_logger_tube <- first(corresp_fs$logger_transparent_chamber)
  site_ID <- str_sub(subsite, start = 1, end = 5)
  
  # finding out if corresponding file exists, and its extension
  dir_exists_loggdata <- dir.exists(paste0(loggerspath,"/",site_ID,"/"))
  if(dir_exists_loggdata){
    f <- list.files(paste0(loggerspath,"/",site_ID,"/"), full.names = T)
    r <- grep(pattern = ".hobo", x = f)
    if(length(r)>0){f <- f[-r]}
    r <- grep(pattern = ".txt", x = f)
    if(length(r)>0){f <- f[-r]}
    f_ext <- file_ext(f)
    i_f_float <- grep(pattern = SN_logger_float, x = f)[1]
    i_f_tube <- grep(pattern = SN_logger_tube, x = f)[1]
  }
  
  if(!is.na(SN_logger_float) & !is.na(i_f_float)){
    is_data_logger_float = T
    message("...reading corresponding temperature logger file for floating chamber")
    if(f_ext[i_f_float]=="xlsx"){
      data_logger_float <- readxl::read_xlsx(f[i_f_float],col_names = T)
    } else if(f_ext[i_f_float]=="csv"){
      data_logger_float <- read.csv(f[i_f_float], header = T)
    }
    data_logger_float <- data_logger_float[,seq(1,3)]
    names(data_logger_float) <- c("sn","datetime","temperature")
    data_logger_float$sn <- SN_logger_float
    if(is.character(data_logger_float$datetime)){
      data_logger_float$datetime <- as.POSIXct(data_logger_float$datetime, tz = 'utc', tryFormats = c("%m/%d/%y %r", "%d/%m/%Y %H:%M"))
    }
    data_logger_float$unixtime <- as.numeric(data_logger_float$datetime)
  } else {
    message("===> no data logger linked to the floating chamber!")
    is_data_logger_float = F}
  if(dim(data_logger_float)[1]<10){
    message("===> not enough data could be linked to the floating chamber!")
    is_data_logger_float = F
  }
  
  if(!is.na(SN_logger_tube) & !is.na(i_f_tube)){
    is_data_logger_tube = T
    message("...reading corresponding temperature logger file for tube chamber")
    if(f_ext[i_f_tube]=="xlsx"){
      data_logger_tube <- readxl::read_xlsx(f[i_f_tube],col_names = T)
    } else if(f_ext[i_f_tube]=="csv"){
      data_logger_tube <- read.csv(f[i_f_tube], header = T, fill = T)
    }
    data_logger_tube <- data_logger_tube[,seq(1,4)]
    names(data_logger_tube) <- c("sn","datetime","temperature","light")
    data_logger_tube$sn <- SN_logger_tube
    if(is.character(data_logger_tube$datetime)){
      if(length(grep(pattern = "AM", x = first(data_logger_tube$datetime)))>0 | length(grep(pattern = "PM", x = first(data_logger_tube$datetime)))>0){
        data_logger_tube$datetime <- as.POSIXct(data_logger_tube$datetime, tz = 'UTC', format = c("%m/%d/%y %r"))
      } else {
        data_logger_tube$datetime <- as.POSIXct(data_logger_tube$datetime, tz = 'UTC', format = "%m/%d/%Y %H:%M")
      }
    }
    data_logger_tube$unixtime <- as.numeric(data_logger_tube$datetime)
  } else {
    is_data_logger_tube = F
    message("===> no data logger linked to the tube chamber!")
  }
  if(dim(data_logger_tube)[1]<10){
    message("===> not enough data could be linked to the tube chamber!")
    is_data_logger_tube = F
  }
  
  
  
  #----1.2. Per-incubation aux variables----
  # Compute Temperature, Area and Volume from fieldsheet info
  
    #For each row in corresp_fs (i.e. for each incubation in subsite-fieldsheet)
    for (incub in seq_along(corresp_fs$plot_id)){
      
        myTemp <- 15 #ºC (a default temperature to run the flux calculation if no temperature is found)
        myPcham <- 100.1 #kPa (default atmospheric pressure is assumed for all chambers)
        
        #Parameters for floating chambers:
        if (corresp_fs$chamber_type[incub] == "floating"){
          if(is_data_logger_float){
            #Calculate median temperature of data_logger_float during duration of floating incubation.
            myTemp <- median(approx(data_logger_float$unixtime, data_logger_float$temperature, xout = as.numeric(corresp_fs$unix_start[incub]:corresp_fs$unix_stop[incub]))$y )
          }
          myArea = 14365.4439 # cm2
          myVtot = 115/2 # L
          
          #Parameters for tube chambers:
        } else if (corresp_fs$chamber_type[incub] == "tube"){
          if(is_data_logger_tube){
            #Calculate median temperature of data_logger_tube during duration of tube incubation.
            myTemp <- median(approx(data_logger_tube$unixtime, data_logger_tube$temperature, xout = as.numeric(corresp_fs$unix_start[incub]:corresp_fs$unix_stop[incub]))$y )
          }
          myArea = pi*12.1**2 # cm2
          myVtot = myArea*as.numeric(corresp_fs$chamber_height_cm[incub])*1e-3 # L
        } else {
          warning("chamber type not correct!")
        }

        
        # --- 1.3. Create auxfile table ----
        # An auxfile table, made of fieldsheet, adding important variables. The auxfile
        # requires start.time and UniqueID.
        # start.time must be in the format "%Y-%m-%d %H:%M:%S"
        
        auxfile_tmp <- data.frame(subsite = subsite,
                                  UniqueID_notime = corresp_fs$UniqueID_notime[incub],
                                  gas_analiser = gs,
                                  start.time = as.POSIXct((corresp_fs$unix_start[incub]), tz = "UTC"),
                                  duration = (corresp_fs$unix_stop[incub]) - (corresp_fs$unix_start[incub]),
                                  water_depth = corresp_fs$water_depth[incub],
                                  Area = myArea,# cm2
                                  Vtot = myVtot,# L
                                  Tcham = myTemp,# ºC
                                  Pcham = myPcham,#kPa
                                  strata = corresp_fs$strata[incub],
                                  chamberType = corresp_fs$chamber_type[incub],
                                  lightCondition = corresp_fs$transparent_dark[incub])
        if(is.null(auxfile)){
          auxfile <- auxfile_tmp
        } else {
          auxfile <- rbind(auxfile, auxfile_tmp)
        }
        
    }
}

rm(f,f_ext,data_logger_float,data_logger_tube,corresp_fs, auxfile_tmp,i,i_f_float,i_f_tube,incub, gs, is_data_logger_float, is_data_logger_tube, isF_incub, isFsubsite, myArea, myPcham, myTemp, myVtot, dir_exists_loggdata, r, site_ID, SN_logger_float, SN_logger_tube, subsite,subsites, myfieldsheets_list)

   
# ---- 4. N2O flux calculation-------
#---TO ADAPT!! to N2O####

#__________________####


#SAVE RESULTS####







