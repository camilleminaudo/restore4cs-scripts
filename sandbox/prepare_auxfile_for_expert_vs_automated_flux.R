

# ---
# Authors: Camille Minaudo
# Project: "RESTORE4Cs"
# date: "June 2024"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
# This script 


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


#################################

# You have to make sure this is pointing to the write folder on your local machine
dropbox_root <- "C:/Users/Camille Minaudo/Dropbox/RESTORE4Cs - Fieldwork/Data" 

############################


# ---- Directories ----
datapath <- paste0(dropbox_root,"/GHG/RAW data")
fieldsheetpath <- paste0(dropbox_root,"/GHG/Fieldsheets")
loggerspath <- paste0(datapath,"/RAW Data Logger")
RData_path <- paste0(dropbox_root,"/GHG/Processed data/RData/")
plots_path <- paste0(dropbox_root,"/GHG/GHG_expert_vs_automated/plots/")
results_path <- paste0(dropbox_root,"/GHG/GHG_expert_vs_automated/results/")

# ----- Data pre-processing and harmonization -----

# ---- List GHG chamber fieldsheets in Dropbox and read them ---
# list filenames
myfieldsheets_list <- list.files(fieldsheetpath, pattern = "Fieldsheet-GHG.xlsx", all.files = T, full.names = T, recursive = T)

# #S4 data is not ready yet
# i <- grep(pattern = "S4", x = myfieldsheets_list) # selecting the files corresponding to the selected sampling campaign
# myfieldsheets_list <- myfieldsheets_list[-i]

# Read all fieldsheets and put them in a single dataframe
fieldsheet <- read_GHG_fieldsheets(myfieldsheets_list)



# ---- Correct fieldsheets in the case of Picarro data ---

read_map_incubations <- function(path2folder){
  
  my_maps_filenames <- list.files(path2folder, pattern = ".csv", all.files = T, full.names = T, recursive = T)
  i <- grep(pattern = "Picarro", x = my_maps_filenames) # selecting the files corresponding to the Picarro only
  my_maps_filenames <- my_maps_filenames[i]
  
  # i <- grep(pattern = "S4-", x = my_maps_filenames)
  # my_maps_filenames <- my_maps_filenames[-i]
  
  isF <- T
  for(my_maps_filename in my_maps_filenames){
    map_incubations_temp <- read.csv(file = my_maps_filename,
                                     header = T, fill = T)
    if(dim(map_incubations_temp)[1]>0){
      map_incubations_temp <- map_incubations_temp[map_incubations_temp$Species ==  "CH4",]
      map_incubations_temp$subsite <- basename(dirname(my_maps_filename))
      if(isF){
        isF <- F
        map_incubations <- map_incubations_temp
      } else {
        map_incubations <- rbind(map_incubations, map_incubations_temp)
      }
    }
    
  }
  
  map_incubations <- data.frame(subsite =map_incubations$subsite,
                                plot = map_incubations$Comment,
                                time_code = map_incubations$Time.Code,
                                start = map_incubations$start_fit,
                                stop = map_incubations$end_fit)
  
  return(map_incubations)
}


# read all the csv files in data_folder, and group into a single one
map_incubations <- suppressWarnings({read_map_incubations(path2folder = datapath)})
map_incubations <- map_incubations[order(map_incubations$start),]

# check the closest incubation in map_incubations for each row in fieldsheet.
# if more than 4 minutes apart, we consider the row in map_incubations out of sampling
fieldsheet_Picarro <- fieldsheet[fieldsheet$gas_analyzer=="Picarro",]
corresponding_row <- NA*fieldsheet_Picarro$plot_id


for (i in seq_along(fieldsheet_Picarro$plot_id)){
  ind <- which.min(abs(fieldsheet_Picarro$unix_start[i] - map_incubations$start))
  if(length(ind)>0){
    if(abs(fieldsheet_Picarro$unix_start[i] - map_incubations$start[ind])<5*60){
      corresponding_row[i] <- ind
    }
  }
}

# finding corresponding row in case of NA in corresponding_row
if(sum(is.na(corresponding_row))>0){
  ind_NAs <- which(is.na(corresponding_row))
  # it is the most probable that the missing info is in-between two rows with all the data we need
  interp <- approx(x = seq_along(fieldsheet_Picarro$plot_id), y = corresponding_row, xout = ind_NAs, rule = 2)$y
  is_integer <- (interp - floor(interp)) == 0
  corresponding_row[ind_NAs][is_integer] <- interp[is_integer]
}
if(sum(is.na(corresponding_row))>0){
  ind_NAs <- which(is.na(corresponding_row))
  if(length(ind_NAs)>0){
    warning("Could not find a corresponding incubation map for the following incubations:")
    fieldsheet_Picarro$uniqID[ind_NAs]
    
    # removing rows from fieldsheet_Picarro with missing correspondance
    fieldsheet_Picarro <- fieldsheet_Picarro[-ind_NAs,]
    corresponding_row <- corresponding_row[-ind_NAs]
  }
}

# replacing unix_startand unix_stop with new values
fieldsheet_Picarro$unix_start <- map_incubations$start[corresponding_row]
fieldsheet_Picarro$unix_stop <- map_incubations$stop[corresponding_row]






# ---- Correct fieldsheets in the case of LiCOR data, whenever available ---

# load incubation map (run get_exact_incubation_times_Licor.R to update it)
map_incubations <- read.csv( file = paste0(dropbox_root,"/GHG/RAW data/RAW Data Licor-7810/map_incubations.csv"))
# only select realistic start/stop
map_incubations <- map_incubations[which(map_incubations$stop-map_incubations$start < 15*60),]

fieldsheet_Licor <- fieldsheet[fieldsheet$gas_analyzer=="LI-COR",]
corresponding_row <- unix_start_corr <- unix_stop_corr <- NA*fieldsheet_Licor$unix_start


for (i in seq_along(fieldsheet_Licor$plot_id)){
  ind <- which.min(abs(fieldsheet_Licor$unix_start[i] - map_incubations$start))
  if(length(ind)>0){
    if(abs(fieldsheet_Licor$unix_start[i] - map_incubations$start[ind])<3*60){
      corresponding_row[i] <- ind
      unix_start_corr[i] <- map_incubations$start[ind]
      unix_stop_corr[i] <- map_incubations$stop[ind]
    }
  }
}
ind_noNAs <- which(!is.na(corresponding_row))
duration_fieldsheet <- fieldsheet_Licor$unix_stop - fieldsheet_Licor$unix_start
fieldsheet_Licor$unix_start[ind_noNAs] <- unix_start_corr[ind_noNAs]
fieldsheet_Licor$unix_stop <- fieldsheet_Licor$unix_start + duration_fieldsheet



# ---- Merging fieldsheets into a single one ---
fieldsheet_LosGatos <- fieldsheet[fieldsheet$gas_analyzer=="Los Gatos",]

fieldsheet <- rbind(fieldsheet_Licor, fieldsheet_LosGatos, fieldsheet_Picarro)

# taking a little margin to avoid critical moments of manipulation with the chamber and manual gas sampling
fieldsheet$unix_start <- fieldsheet$unix_start+30
fieldsheet$unix_stop <- fieldsheet$unix_stop-30

# recalculating start and stop in proper formats
fieldsheet$timestamp_start <- as.POSIXct(fieldsheet$unix_start, tz = "UTC", origin = "1970-01-01")
fieldsheet$timestamp_stop <- as.POSIXct(fieldsheet$unix_stop, tz = "UTC", origin = "1970-01-01")

fieldsheet$start_time <- strftime(fieldsheet$timestamp_start, format="%H:%M:%S", tz = 'utc')
fieldsheet$end_time <- strftime(fieldsheet$timestamp_stop, format="%H:%M:%S", tz = 'utc')

fieldsheet <- fieldsheet[order(fieldsheet$subsite),]



# ----- Data selection and building auxilliary table -----

# Here, we only focus on open water measurements with the floating chamber
fieldsheet <- fieldsheet[which(fieldsheet$chamber_type=="floating" & fieldsheet$strata=="open water" & fieldsheet$water_depth>0),]


# For each subsite in fieldsheet, go through each incubation and compute co2 and ch4 fluxes
subsites <- unique(fieldsheet$subsite)

isF_incub <- T
isFsubsite <- T
auxfile <- NULL

for (subsite in subsites){
  message("Now processing ",subsite)
  
  corresp_fs <- fieldsheet[fieldsheet$subsite == subsite,]
  gs <- first(corresp_fs$gas_analyzer)
  
  
  # read corresponding temperature logger file and keep initial temperature
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
  
  
  
  setwd(RData_path)
  if(gs== "LI-COR"){
    gs_suffix <- "LI-7810"
  } else {
    gs_suffix <- gs
  }
  
  
  if(file.exists(paste0(subsite,"_",gs_suffix,".RData"))){
    load(file = paste0(subsite,"_",gs_suffix,".RData"))
    
    
    for (incub in seq_along(corresp_fs$plot_id)){
      my_incub <- mydata[as.numeric(mydata$POSIX.time)> corresp_fs$unix_start[incub] &
                           as.numeric(mydata$POSIX.time)< corresp_fs$unix_stop[incub],]
      my_incub <- my_incub[!is.na(my_incub$CO2dry_ppm),]
      
      if (dim(my_incub)[1]>0){
        
        # Compute Temperature, Area and Volume from fieldsheet info
        myTemp <- 15 # a default temperature to run the flux calculation...
        if (corresp_fs$chamber_type[incub] == "floating"){
          if(is_data_logger_float){
            myTemp <- median(approx(data_logger_float$unixtime, data_logger_float$temperature, xout = as.numeric(my_incub$POSIX.time))$y )
          }
          myArea = 14365.4439 # cm2
          myVtot = 115/2 # L
        } else if (corresp_fs$chamber_type[incub] == "tube"){
          if(is_data_logger_tube){
            myTemp <- median(approx(data_logger_tube$unixtime, data_logger_tube$temperature, xout = as.numeric(my_incub$POSIX.time))$y )
          }
          myArea = pi*12.1**2 # cm2
          myVtot = myArea*as.numeric(corresp_fs$chamber_height_cm[incub])*1e-3 # L
        } else {
          warning("chamber type not correct!")
        }
        myPcham <- 100.1 #kPa
        
        # --- Create auxfile table ----
        # An auxfile table, made of fieldsheet, adding important variables. The auxfile
        # requires start.time and UniqueID.
        # start.time must be in the format "%Y-%m-%d %H:%M:%S"
        
        auxfile_tmp <- data.frame(subsite = subsite,
                                  UniqueID = corresp_fs$uniqID[incub],
                                  gas_analiser = gs,
                                  start.time = as.POSIXct((corresp_fs$unix_start[incub]), tz = "UTC"),
                                  duration = (corresp_fs$unix_stop[incub]) - (corresp_fs$unix_start[incub]),
                                  water_depth = corresp_fs$water_depth[incub],
                                  Area = myArea,
                                  Vtot = myVtot,
                                  Tcham = myTemp,
                                  Pcham = myPcham,
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
    
  } else {
    message("---> Could not find corresponding ",gs," data")
  }
  
}

dim(auxfile)

# we only keep incubations longer than 100 secs
auxfile <- auxfile[auxfile$duration>100,]
# we only keep incubations where chamber dimensions are known
auxfile <- auxfile[!is.na(auxfile$Vtot),] # in case chamber height is not specified in the fieldsheet...
auxfile <- auxfile[!is.na(auxfile$Area),]
# we only keep incubations with known Temperature at start
auxfile <- auxfile[!is.na(auxfile$Tcham),]


dim(auxfile)

# saving fluxes estimates
setwd(dirname(results_path))

write.csv(x = auxfile, file = "auxfile.csv", 
          row.names = F)







