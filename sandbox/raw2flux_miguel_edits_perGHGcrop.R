

# ---
# Authors: Camille Minaudo edits by MIGUEL CABRERA (April 2025)
# Project: "RESTORE4Cs"
# date: "Oct 2023"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
# This script computes fluxes for all incubations present in the fieldsheet for a given sampling season
# User has to specify which sampling has to be processed, and all the rest is done automatically.

#TO-DO------
#EDITS Miguel: 
  #DONE: Comments and formating sections
  #DONE: Improvement Picarro starstop corrections
  #DONE: Licor map_injections corrected
  #DONE: Per-GHG cropping of start and stop after visual inspection of every incubation.
  #DONE: Use meteo-stations instead of data-loggers for flux.term (section Auxtable prep)

#TO-DO: Update ebullition method with camille's new function (Section Diffusion vs ebullition)
#TO-DO: set criteria per-GHG for best.flux, & adjust k.mult in goflux (if needed)



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


#---USER OPTIONS ----- 
#sampling is the pattern to select in the fieldsheet filenames, every matching fieldsheet will have their incubations processed.
sampling <- "S" 

#create/update RDATA?
harmonize2RData <- F

#Save plots from goflux?
doPlot <- T

#Minimum duration of incubation (seconds) to calculate flux for
minimum_duration_for_flux<- 100

#Calculate and plot high-leverage points for LM? (For quality control, plots Cook's values for incubations that contain data higher than the typical threshold of  4 / n.obs)
flag_lm_leverage<- F


# ---- Directories ----
dropbox_root <- "C:/Users/Miguel/Dropbox/RESTORE4Cs - Fieldwork/Data" # You have to make sure this is pointing to the write folder on your local machine
datapath <- paste0(dropbox_root,"/GHG/RAW data")
fieldsheetpath <- paste0(dropbox_root,"/GHG/Fieldsheets")
loggerspath <- paste0(datapath,"/RAW Data Logger")
RData_path <- paste0(dropbox_root,"/GHG/Processed data/RData/")
quality_path<- paste0(dropbox_root, "/GHG/Working data/Incubation_quality/") #quality assessment, crop decisions and flags after inspection.
meteo_path<- paste0(dropbox_root, "/Meteo/Formated_data/")

# results_path<- paste0(dropbox_root,"/GHG/Processed data/computed_flux/")
#Final run in local pc
results_path<- "C:/Users/Miguel/Dropbox/TEST_quality_raw2flux/Per_GHG_cropped_fluxes/"
plots_path <- paste0(results_path,"level_incubation/")

#test_newRdata --------
#The results are the same after updating Rdata  with correct picarro h2o units
# results_path<- "C:/Users/Miguel/Dropbox/testing_aquaGHG/"
# plots_path <- paste0(results_path,"level_incubation/")


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



## ---- Correct Picarro fieldsheets ----

#MIGUEL: i.e. substitute start-stoptimes of fieldsheets with those recorded in picarro flux aplication, if available
read_map_incubations <- function(path2folder, sampling){
  
  my_maps_filenames <- list.files(path2folder, pattern = ".csv", all.files = T, full.names = T, recursive = T)
  i <- grep(pattern = "Picarro", x = my_maps_filenames) # selecting the files corresponding to the Picarro only
  my_maps_filenames <- my_maps_filenames[i]
  i <- grep(pattern = sampling, x = my_maps_filenames) # selecting the files corresponding to the selected sampling campaign
  my_maps_filenames <- my_maps_filenames[i]
  
  
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
map_incubations <- suppressWarnings({read_map_incubations(path2folder = datapath, sampling = sampling)})
map_incubations <- map_incubations[order(map_incubations$start),]

# check the closest incubation in map_incubations for each row in fieldsheet.
# if more than 4 minutes apart, we consider the row in map_incubations out of sampling
fieldsheet_Picarro <- fieldsheet[fieldsheet$gas_analyzer=="Picarro",]
corresponding_row <- NA*fieldsheet_Picarro$plot_id


for (i in seq_along(fieldsheet_Picarro$plot_id)){
  ind <- which.min(abs(fieldsheet_Picarro$unix_start[i] - map_incubations$start))
  if(length(ind)>0){
    if(abs(fieldsheet_Picarro$unix_start[i] - map_incubations$start[ind])<4*60){
      corresponding_row[i] <- ind
    }
  }
}

#Check if the 5minute lookup worked for all
if(sum(is.na(corresponding_row))==0){
  message("All fielsheet incubations could be assigned to corrected start-stop times from picarro software")
}else if (sum(is.na(corresponding_row))>0){
  message(paste("looking in a plus or minus 5-minute window did not work for", sum(is.na(corresponding_row)), "incubations:"))
  print(fieldsheet_Picarro[which(is.na(corresponding_row)),]$uniqID)
}
#the following incubations should not be found (no map_incubation for them):
#"s1-ri-a1-18-o-d-14:28" "s1-ri-a1-19-o-d-14:40"

# finding corresponding row in case of NA in corresponding_row
if(sum(is.na(corresponding_row))>0){
  ind_NAs <- which(is.na(corresponding_row))
  # it is the most probable that the missing info is in-between two rows with all the data we need
  interp <- approx(x = seq_along(fieldsheet_Picarro$plot_id), y = corresponding_row, xout = ind_NAs)$y
  is_integer <- (interp - floor(interp)) == 0
  corresponding_row[ind_NAs][is_integer] <- interp[is_integer]
}

if(sum(is.na(corresponding_row))>0){
  ind_NAs <- which(is.na(corresponding_row))
  message("Could not find a corresponding incubation map for the following incubations:")
  print(fieldsheet_Picarro$uniqID[ind_NAs])
  
}
# replacing unix_start and unix_stop with new values
fieldsheet_Picarro$unix_start <- map_incubations$start[corresponding_row]
fieldsheet_Picarro$unix_stop <- map_incubations$stop[corresponding_row]

#The above chunk causes incubations without match in map_incubations to lose their unix_start and unix_stop, re-calculate them
fieldsheet_Picarro<- fieldsheet_Picarro %>% 
  mutate(unix_start=if_else(is.na(unix_start), as.numeric(as.POSIXct(paste(date, start_time),tz = "UTC")), unix_start),
         unix_stop=if_else(is.na(unix_stop), as.numeric(as.POSIXct(paste(date, end_time),tz = "UTC")), unix_stop))



#Check that unix_start and Unix_stop are not duplicated after picarro time-matching: 
{duplicated_starts_picarro<- fieldsheet_Picarro %>% 
    filter(unix_start %in% fieldsheet_Picarro[duplicated(fieldsheet_Picarro$unix_start,incomparables = F),]$unix_start)
  duplicated_stops_picarro<-fieldsheet_Picarro %>% 
    filter(unix_stop %in% fieldsheet_Picarro[duplicated(fieldsheet_Picarro$unix_stop,incomparables = F),]$unix_stop)
  
  #Check potential duplicate incubations created by the approach, code below.
  if(sum(dim(duplicated_starts_picarro)[1],dim(duplicated_stops_picarro)[1])>0){
    message(paste("CAUTION: picarro correction duplicates incubations, check duplicates and correct fieldsheet times"))
    print(duplicated_starts_picarro$uniqID)
  }else{
    message(paste("All Picarro fieldsheet start-stop times are now corrected"))
    rm(duplicated_starts_picarro, duplicated_stops_picarro)}
  }

## ---- Correct LiCOR fieldsheets ----

#Correct Licor start-stop times whenever remarks are available
# load incubation map: updated with last incubations and corrected to get all non-duplicated remarks (removing also remarks that cause wrong-assignments)
map_incubations <- read.csv( file = paste0(dropbox_root,"/GHG/RAW data/RAW Data Licor-7810/map_incubations_touse.csv"))

# only select realistic start/stop
map_incubations <- map_incubations[which(map_incubations$stop-map_incubations$start < 15*60),]

fieldsheet_Licor <- fieldsheet[fieldsheet$gas_analyzer=="LI-COR",]
corresponding_row <- unix_start_corr <- unix_stop_corr <- NA*fieldsheet_Licor$unix_start

#Lookup for corrected start-time correspondence (3 minute window) 
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


#Check how many were found:
if(sum(!is.na(ind_noNAs))>0){
  fieldsheet_Licormapped<- fieldsheet_Licor
  fieldsheet_Licormapped$foundinmap<- FALSE
  fieldsheet_Licormapped$foundinmap[ind_noNAs]<- TRUE
  fieldsheet_Licormapped %>% 
    mutate(sampling=substr(subsite,1,5)) %>% 
    group_by(sampling) %>% 
    summarise(foundinmap=sum(foundinmap), total_incubations= n(), percentfound=foundinmap/total_incubations*100)
  message(paste("Licor start-stop times corrected for",sum(!is.na(ind_noNAs)),"incubations" ))
  
}


#Correct start-stop times based on start remark (if found) + duration fieldsheets
duration_fieldsheet <- fieldsheet_Licor$unix_stop - fieldsheet_Licor$unix_start
fieldsheet_Licor$unix_start[ind_noNAs] <- unix_start_corr[ind_noNAs]
fieldsheet_Licor$unix_stop <- fieldsheet_Licor$unix_start + duration_fieldsheet


## ---- Merging fieldsheets----
#Import los GATOS fieldsheets: start and stop times directly from fieldsheets (no map_incubations)
fieldsheet_LosGatos <- fieldsheet[fieldsheet$gas_analyzer=="Los Gatos",]

#Combine all fieldsheets with corrected times based on map_incubations when available
fieldsheet <- rbind(fieldsheet_Licor, fieldsheet_LosGatos, fieldsheet_Picarro)


#This is legacy, the custom-cropping implemented in auxfile per-gas makes this section redundant.
# fieldsheet$unix_start <- fieldsheet$unix_start+margin_s_after_start
# fieldsheet$unix_stop <- fieldsheet$unix_stop-margin_s_before_end


# recalculating start and stop in proper formats
fieldsheet$timestamp_start <- as.POSIXct(fieldsheet$unix_start, tz = "UTC", origin = "1970-01-01")
fieldsheet$timestamp_stop <- as.POSIXct(fieldsheet$unix_stop, tz = "UTC", origin = "1970-01-01")

fieldsheet$start_time <- strftime(fieldsheet$timestamp_start, format="%H:%M:%S", tz = 'utc')
fieldsheet$end_time <- strftime(fieldsheet$timestamp_stop, format="%H:%M:%S", tz = 'utc')

fieldsheet <- fieldsheet[order(fieldsheet$subsite),]



#_________________####
#RDATA SAVE (optional) ####
#Licor import of S4-CU corrected (Instrument not in UTC time)
#Import and store measurements to RData ---
if(harmonize2RData){
  data_folders <- list.dirs(datapath, full.names = T, recursive = T)[-1]
  i <- grep(pattern = sampling, x = data_folders) # selecting the files corresponding to the selected sampling campaign
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
      
      setwd(RData_path)
      
      # save this dataframe as a new RData file
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
      message("S4-CU Licor not in UTC time!!")  
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

#_________________####
#FLUX CALCULATION####


#----0. Load crop & Meteo-----

#Import table with decisions and cropping amounts: all cropping is based on last run of 5-s margin script, so we have to re-add this 5 seconds to the cropping margins
croping<- read_xlsx(paste0(quality_path, "Inspection_table_allincubations_tocrop_perGHG.xlsx"),na = "NA")

#Simplify and add the additional 5s of cropping to all crop-start and crop-end: Inspection and cropping decisions were taken on fluxes calculated based on 5-s margin. WE re-add this 5s of margin here for all incubations (custom-cropped or not) 
croping_simple<- croping %>%
  select(UniqueID, 
         co2_start_crop_s, co2_end_crop_s,
         ch4_start_crop_s, ch4_end_crop_s) %>%
  mutate(h2o_start_crop_s=5, h2o_end_crop_s=5) %>% #add 5s crop for all h2o
  mutate(co2_start_crop_s=case_when(is.na(co2_start_crop_s)~5,
                                    TRUE~co2_start_crop_s+5),
         co2_end_crop_s=case_when(is.na(co2_end_crop_s)~5,
                                  TRUE~co2_end_crop_s+5),
         ch4_start_crop_s=case_when(is.na(ch4_start_crop_s)~5,
                                    TRUE~ch4_start_crop_s+5),
         ch4_end_crop_s=case_when(is.na(ch4_end_crop_s)~5,
                                  TRUE~ch4_end_crop_s+5))# add 5s to all crop decissions for co2 and ch4

#croping_simple is used to adapt the auxfile for each gas in the section Goflux calculation

#Load data for all meteostations and correspondence of subsite-meteostation (used to subset the correct data)

meteocorresp <- read_csv(paste0(meteo_path, "correspondence_subsite-meteostation.csv"),show_col_types = F)

meteodata <- read_csv(paste0(meteo_path,"allmeteostations_onlysamplingdays.csv"),show_col_types = F)


# ----- 1. Auxtable prep -----

# This section loops over each subsite and performs the flux calculation for CO2 and CH4
subsites <- unique(fieldsheet$subsite)

isF_incub <- T
isFsubsite <- T
for (subsite in subsites){
  message("Now processing ",subsite)
  
  corresp_fs <- fieldsheet[fieldsheet$subsite == subsite,]
  gs <- first(corresp_fs$gas_analyzer)
  
  #Load the appropriate meteo-data for this subsite:
  meteostation<- meteocorresp[meteocorresp$site_subsite==substr(subsite,4,8),]$station_id
  subsite_meteo<- meteodata %>% filter(station_id==meteostation)
  
  setwd(RData_path)
  if(gs== "LI-COR"){
    gs_suffix <- "LI-7810"
  } else {
    gs_suffix <- gs
  }
  
  
  ##----1.1. Per-incubation aux variables----
  auxfile <- NULL
  
  if(file.exists(paste0(subsite,"_",gs_suffix,".RData"))){
    load(file = paste0(subsite,"_",gs_suffix,".RData"))
    
    #Correct the units for picarro H2O,import function assumes mmol/mol, actual units of our instrument are % (i.e. 100* mol/mol). We correct the Rdata by multiplying the concentration by 10
    # if (gs=="Picarro"){mydata$H2O_ppm<- mydata$H2O_ppm*10}
    #The correction of units is now already implemented in Rdata as of 6/5/2025
    
    #Loop over each incubation: 
    for (incub in seq_along(corresp_fs$plot_id)){
      my_incub <- mydata[as.numeric(mydata$POSIX.time)> corresp_fs$unix_start[incub] &
                           as.numeric(mydata$POSIX.time)< corresp_fs$unix_stop[incub],]
      my_incub <- my_incub[!is.na(my_incub$CO2dry_ppm),]
      
      if (dim(my_incub)[1]>0){
        
        # Compute Temperature (ºC) and pressure (kPa) for the initial moment of incub using the meteo data (using approx to interpolate from hourly data):
        myTemp <- 15 # a default temperature to run the flux calculation...
        myPcham <- 100.1 #kPa, a default atm pressure to run the flux calculation
        
        myTemp <- first(approx(as.numeric(subsite_meteo$datetime_utc), subsite_meteo$temp_c, xout=as.numeric(my_incub$POSIX.time))$y)
        
        myPcham<- first(approx(as.numeric(subsite_meteo$datetime_utc), subsite_meteo$Patm_hPa/10, xout=as.numeric(my_incub$POSIX.time))$y)
        
        #Get Area and Volume of chamber from fieldsheet data
        if (corresp_fs$chamber_type[incub] == "floating"){
          #Floating chamber used in Restore4Cs is half sphere of 38cm diameter. radius=19cm
          myArea<- pi*19^2 #cm2
          myVtot<- (((4/3)*pi*(19^3))/1000 )/2# L
          
        } else if (corresp_fs$chamber_type[incub] == "tube"){
          #Tube chamber area is multiplied by height to get the Volume
          myArea = pi*12.1**2 # cm2
          myVtot = myArea*as.numeric(corresp_fs$chamber_height_cm[incub])*1e-3 # L
        } else {
          warning("chamber type not correct!")
        }
        
        # --- 1.2. Create auxfile table ----
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
  
# ---- 2. Checks-------

  if (length(auxfile)>1){
    
    auxfile$Tcham[is.na(auxfile$Tcham)] <- mean(auxfile$Tcham, na.rm = T)
    
    # we only keep incubations longer than X secs (after cropping Xs start and Xs end): see USER OPTIONS
    auxfile <- auxfile[auxfile$duration>minimum_duration_for_flux,]
    
    # we only keep incubations where chamber dimensions are known
    auxfile <- auxfile[!is.na(auxfile$Vtot),] # in case chamber height is not specified in the fieldsheet...
    auxfile <- auxfile[!is.na(auxfile$Area),]
  
    
#3. Per-GHG crop------    
    
    #Create 3 separate auxfiles, one per GHG, with modified start.time and duration based on cropping decissions per GHG (h2o is only cropped by 5s start and end: inspection runned with these margings)
    
    auxfile_co2<- auxfile %>% 
      left_join(croping_simple %>% select(UniqueID,co2_start_crop_s,co2_end_crop_s), by="UniqueID") %>% 
      mutate(start.time=start.time+co2_start_crop_s,
             duration=duration-(co2_start_crop_s+co2_end_crop_s)) %>% 
      select(-c(co2_start_crop_s, co2_end_crop_s))
    
    auxfile_ch4<- auxfile %>% 
      left_join(croping_simple %>% select(UniqueID,ch4_start_crop_s,ch4_end_crop_s), by="UniqueID") %>% 
      mutate(start.time=start.time+ch4_start_crop_s,
             duration=duration-(ch4_start_crop_s+ch4_end_crop_s)) %>% 
      select(-c(ch4_start_crop_s, ch4_end_crop_s))
    
    auxfile_h2o<- auxfile %>% 
      left_join(croping_simple %>% select(UniqueID,h2o_start_crop_s,h2o_end_crop_s), by="UniqueID") %>% 
      mutate(start.time=start.time+h2o_start_crop_s,
             duration=duration-(h2o_start_crop_s+h2o_end_crop_s)) %>% 
      select(-c(h2o_start_crop_s, h2o_end_crop_s))
    
    
    #Define the measurements window of observation (one per GHG)
    #suppressMessagesto avoid the "Do not loop through more than 20 measurements...." message
    suppressMessages({
      mydata_ow_co2 <- obs.win(inputfile = mydata, auxfile = auxfile_co2,
                         obs.length = auxfile_co2$duration, shoulder = 0)
    
      mydata_ow_ch4 <- obs.win(inputfile = mydata, auxfile = auxfile_ch4,
                             obs.length = auxfile_ch4$duration, shoulder = 0)
    
      mydata_ow_h2o <- obs.win(inputfile = mydata, auxfile = auxfile_h2o,
                             obs.length = auxfile_h2o$duration, shoulder = 0)
    })
    # Join mydata_ow with info on start end incubation  (one per GHG)
    mydata_auto_co2 <- lapply(seq_along(mydata_ow_co2), join_auxfile_with_data.loop, flux.unique = mydata_ow_co2) %>%
      map_df(., ~as.data.frame(.x)) %>%
      # Add instrument precision for each gas
      mutate(CO2_prec = first(mydata$CO2_prec), CH4_prec = first(mydata$CH4_prec), 
             N2O_prec = first(mydata$N2O_prec), H2O_prec = first(mydata$H2O_prec))
    
    mydata_auto_ch4 <- lapply(seq_along(mydata_ow_ch4), join_auxfile_with_data.loop, flux.unique = mydata_ow_ch4) %>%
      map_df(., ~as.data.frame(.x)) %>%
      # Add instrument precision for each gas
      mutate(CO2_prec = first(mydata$CO2_prec), CH4_prec = first(mydata$CH4_prec), 
             N2O_prec = first(mydata$N2O_prec), H2O_prec = first(mydata$H2O_prec))
    
    mydata_auto_h2o <- lapply(seq_along(mydata_ow_h2o), join_auxfile_with_data.loop, flux.unique = mydata_ow_h2o) %>%
      map_df(., ~as.data.frame(.x)) %>%
      # Add instrument precision for each gas
      mutate(CO2_prec = first(mydata$CO2_prec), CH4_prec = first(mydata$CH4_prec), 
             N2O_prec = first(mydata$N2O_prec), H2O_prec = first(mydata$H2O_prec))

    #This chunk does not do anything (all info from auxfile is already in mydata_auto_co2)    
    # # Additional auxiliary data required for flux calculation. (one per GHG) 
    # mydata_auto_co2 <- mydata_auto_co2 %>%
    #   left_join(auxfile_co2 %>% select(UniqueID, Area, Vtot, Tcham, Pcham))
    # 
    # mydata_auto_ch4 <- mydata_auto_ch4 %>%
    #   left_join(auxfile_ch4 %>% select(UniqueID, Area, Vtot, Tcham, Pcham))
    # 
    # # mydata_auto_h2o <- mydata_auto_h2o %>%
    # #   left_join(auxfile_h2o %>% select(UniqueID, Area, Vtot, Tcham, Pcham))
    # 
    # # Add instrument precision for each gas  (one per GHG)
    # mydata_auto_co2 <- mydata_auto_co2 %>%
    #   mutate(CO2_prec = first(mydata$CO2_prec), CH4_prec = first(mydata$CH4_prec), 
    #          N2O_prec = first(mydata$N2O_prec), H2O_prec = first(mydata$H2O_prec))
    # 
    # mydata_auto_ch4 <- mydata_auto_ch4 %>%
    #   mutate(CO2_prec = first(mydata$CO2_prec), CH4_prec = first(mydata$CH4_prec), 
    #          N2O_prec = first(mydata$N2O_prec), H2O_prec = first(mydata$H2O_prec))
    # 
    # mydata_auto_h2o <- mydata_auto_h2o %>%
    #   mutate(CO2_prec = first(mydata$CO2_prec), CH4_prec = first(mydata$CH4_prec), 
    #          N2O_prec = first(mydata$N2O_prec), H2O_prec = first(mydata$H2O_prec))
    # 
    
#4. GOflux call -----
    
    
    # Calculate fluxes
    #K.mult sets the multiplier for maximum Kappa (curvature parameter of H.M. model), defaults to 1 (no multiplier). Allow k.mult of 1.1 to discard HM models that exceed K.max
    
    CO2_results_auto <- goFlux(dataframe = mydata_auto_co2, gastype = "CO2dry_ppm",k.mult = 1.1)
    CH4_results_auto <- goFlux(dataframe = mydata_auto_ch4, gastype = "CH4dry_ppb",k.mult = 1.1)
    H2O_results_auto <- goFlux(dataframe = mydata_auto_h2o, gastype = "H2O_ppm", k.mult = 1.1)
    
    
    
    # Use best.flux to select the best flux estimates (LM or HM)
    # based on a list of criteria

    ##BEST FLUX criteria-----
    #Possibility to set different criteria for each gas
    
    #Possible selection criteria: c("g.factor", "kappa", "MAE", "RMSE", "AICc", "SE" )
    
      #g.factor: arbitrary limit for HM/LM, defaults to 2 (if non-linearity is expected, this is not a good selection criteria). If threshold is exceeded, LM is chosen. Decission: set threshold afterwards
    
      #kappa: compares the ratio k.HM/k.max to a threshold (default=1), is a measure of how extreme is the curvature fitted by HM against the theoretical maximum curvature (determined by LM.flux, MDF  and duration). IF the threshold is exceeded, LM is chosen.  Decision: USE, check distribution of kappa and, if needed set a thrsehold lower than 1 so that kappa is used for selection.
    
      #modelfit: equally weighted comparison for all of c("MAE", "RMSE", "AICc", "SE") included in criteria. Selection of model that performs best across most of the criteria included, in case of tie, HM is chosen.
    
    #Warnings: "MDF", "nb.obs", "p-value", "intercept". Each of these included is issued as a warning, MDF (flux below detection), nb.obs (warn if less than 60s incubation), p-value (only for LM, return ~the sigificance of slope), intercept (allows to warn for starting concentrations outside thresholds)
    
    
    
    #Common for all 3 gasses: only RMSE as measure of residuals (MAE conveys the same info)
    
    criteria_co2 <- criteria_ch4 <- criteria_h2o<- c("RMSE", "AICc","SE","kappa", "MDF")
    
    CO2_flux_res_auto <- best.flux(CO2_results_auto, criteria_co2)
    CH4_flux_res_auto <- best.flux(CH4_results_auto, criteria_ch4)
    H2O_flux_res_auto <- best.flux(H2O_results_auto, criteria_h2o)
    
    
    
    
    #----FLAG high-leverage points-----
    
    if (flag_lm_leverage==T){
    # Function to generate gas vs elapsed time and Cook's Distance vs elapsed time plots
    generate_combined_plots <- function(df, gas) {
      # Create the formula dynamically
      formula_str <- paste0(gas, " ~ unixtime")
      
      # Generate the plots
      results_list <- df %>%
        group_by(UniqueID) %>%
        group_split() %>%
        lapply(function(group_df) {
          n <- nrow(group_df)
          if (n < 3) return(NULL)  # Skip groups too small to assess
          
          # Calculate the elapsed time for each UniqueID group
          group_df$elapsed_time <- group_df$unixtime - min(group_df$unixtime)
          
          # Fit the linear model once
          model <- lm(as.formula(formula_str), data = group_df)
          
          # Extract Cook's Distance, Gas, and elapsed time
          cook <- cooks.distance(model)
          gas_values <- group_df[[gas]]
          elapsed_time <- group_df$elapsed_time
          
          # Define the Cook's Distance threshold
          cook_thresh <- 4 / n  # Adjust as needed for stricter or looser thresholds
          
          # Check if any Cook's Distance exceeds the threshold
          if (any(cook > cook_thresh)) {
            
            long_df<- group_df %>% 
              select(UniqueID, elapsed_time) %>% 
              mutate(cook_values=cook, gas_values=gas_values) %>% 
              mutate(cookexceed=cook_values>cook_thresh) %>% 
              pivot_longer(cols = c(cook_values, gas_values),names_to = "names",values_to = "values")
            
            # Create a combined plot using facets
            combined_plot <- ggplot(long_df, aes(x=elapsed_time, y=values,col=cookexceed)) +
              # Gas vs elapsed_time plot
              geom_point() +
              labs(x = "Elapsed Time (seconds)", 
                   title = paste(gas,  unique(group_df$UniqueID))) +
              theme_minimal() +
              facet_wrap(~names, ncol = 1, scales = "free_y")
            
            combined_plot$plot_env$UniqueID <- unique(group_df$UniqueID)
            # Return the combined plot for each UniqueID
            return(combined_plot)
          } else {
            return(NULL)  # Skip if no Cook's Distance exceeds the threshold
          }
        })
      
      # Filter out NULL results (in case no Cook's Distance exceeded the threshold)
      results_list <- results_list[!sapply(results_list, is.null)]
      
      return(results_list)
    }
    
    
    #Create quality plots for CO2
    Cooksplots_co2 <- generate_combined_plots(mydata_auto_co2, gas = "CO2dry_ppm")

    #Save quality plots for co2
    myfilename_cooks <- paste("Cooks",subsite, as.character(as.Date(last(auxfile$start.time))),sep="_")
    
    pdf(file = paste0(plots_path,myfilename_cooks,".pdf"))
    invisible(lapply(Cooksplots_co2, print))
    dev.off()
    
    }
        
    #----5. Plot Goflux (optional)----
    if(doPlot){
      # Plots results
      # Make a list of plots of all measurements, for each gastype
      CO2_flux_plots <- flux.plot(CO2_flux_res_auto, mydata_auto_co2, "CO2dry_ppm")
      CH4_flux_plots <- flux.plot(CH4_flux_res_auto, mydata_auto_ch4, "CH4dry_ppb")
      H2O_flux_plots <- flux.plot(H2O_flux_res_auto, mydata_auto_h2o, "H2O_ppm")
      
      # Combine plot lists into one list
      flux_plot.ls <- c(CO2_flux_plots, CH4_flux_plots, H2O_flux_plots)
      
      # Save plots to pdf in result_subfolder
      myfilename <- paste(subsite, as.character(as.Date(last(auxfile$start.time))),sep="_")
      flux2pdf(flux_plot.ls, outfile = paste0(plots_path,myfilename,".pdf"))
      
    }
    
    #6. Diffusion vs ebullition----
    # estimate ch4 diffusion and ebullition components
    CH4_res_meth1 <- CH4_flux_res_auto  #keep best flux from goflux
    CH4_res_meth1$total_estimated <- NA
    CH4_res_meth1$ebullition <- NA
    CH4_res_meth1$diffusion <- NA
    
    #Only for plots with water, apply estimate of ebullition and diffusion
    for (i in which(auxfile_ch4$water_depth>0)){
      if(auxfile_ch4$water_depth[i]>0){
        my_incub <- mydata[as.numeric(mydata$POSIX.time)> auxfile_ch4$start.time[i] &
                             as.numeric(mydata$POSIX.time)< auxfile_ch4$start.time[i]+auxfile_ch4$duration[i],]
        my_incub <- my_incub[!is.na(my_incub$CO2dry_ppm),]
        # calling dedicated function
        df_ebull <- separate_ebullition_from_diffusion(my_incub, UniqueID = auxfile_ch4$UniqueID[i], doPlot=F)
        # computing fluxes
        H2O_mol = my_incub$H2O_ppm / (1000*1000)
        myfluxterm <- flux.term(auxfile_ch4$Vtot[i], auxfile_ch4$Pcham[i], auxfile_ch4$Area[i],
                                auxfile_ch4$Tcham[i], first(H2O_mol))
        CH4_flux_total <- df_ebull$delta_ch4/df_ebull$duration*myfluxterm # nmol/m2/s
        CH4_flux_diff <- df_ebull$avg_diff_slope*myfluxterm # nmol/m2/s
        CH4_flux_ebull <- CH4_flux_total - CH4_flux_diff
      } else {
        #For plots without water, set ebullition to 0 
        CH4_flux_total <- CH4_flux_ebull <- CH4_res_meth1$best.flux[which(CH4_res_meth1$UniqueID==auxfile_ch4$UniqueID[i])]
        CH4_flux_ebull <- 0
      }
      CH4_res_meth1$total_estimated[which(CH4_res_meth1$UniqueID==auxfile_ch4$UniqueID[i])] <- CH4_flux_total
      CH4_res_meth1$ebullition[which(CH4_res_meth1$UniqueID==auxfile_ch4$UniqueID[i])] <- CH4_flux_ebull
      CH4_res_meth1$diffusion[which(CH4_res_meth1$UniqueID==auxfile_ch4$UniqueID[i])] <- CH4_flux_diff
    }
    
    CH4_res_meth1$ebullition[which(CH4_res_meth1$ebullition<0)] <- 0
    
  
    
    #__________________####
    #SAVE RESULTS####
    #---- 1. Save incubation-level ----
    
    #Save full results in incubation-level path
    myfilenameCO2 <- paste(subsite,"co2_fluxes", as.character(as.Date(last(auxfile$start.time))),sep="_")
    myfilenameCH4 <- paste(subsite,"ch4_fluxes", as.character(as.Date(last(auxfile$start.time))),sep="_")
    write.csv(x = CO2_flux_res_auto, file = paste0(plots_path,myfilenameCO2,".csv"), row.names = F)
    write.csv(x = CH4_res_meth1, file = paste0(plots_path,myfilenameCH4,".csv"), row.names = F)
    
    
    #---- 2. Save Season simple ----
    
    #Join simplified results for both gases in results_path
    table_results <- auxfile %>%
      left_join(CO2_flux_res_auto %>% select(UniqueID, LM.flux, HM.flux, best.flux, model, quality.check),by = "UniqueID") %>%
      rename(CO2_LM.flux = LM.flux, CO2_HM.flux = HM.flux, CO2_best.flux = best.flux, CO2_best.model = model,
             CO2_quality.check = quality.check) %>%
      left_join(CH4_res_meth1 %>% select(UniqueID, LM.flux, HM.flux, best.flux, best.flux, model, quality.check, 
                                         diffusion, ebullition),by = "UniqueID") %>%
      rename(CH4_LM.flux = LM.flux, CH4_HM.flux = HM.flux, CH4_best.flux = best.flux, CH4_best.model = model, 
             CH4_quality.check = quality.check, CH4_diffusive_flux = diffusion, CH4_ebullitive_flux = ebullition)
    
    if (isFsubsite){
      isFsubsite <- F
      table_results_all <- table_results
    } else {
      table_results_all <- rbind(table_results_all, table_results)
    }
  } else {
    message("----- auxfile could not be built for this susbite. It seems data from ",gs," is missing.")
  }
  
}

# setwd(results_path)
myfilename <- paste(sampling,"fluxes",min(as.Date(table_results_all$start.time)),"to",
                    max(as.Date(table_results_all$start.time)), sep = "_")
write.csv(x = table_results_all, file = paste0(results_path, myfilename,".csv"), row.names = F)

s1<-table_results_all %>% filter(grepl("S1", subsite))
s2<-table_results_all %>% filter(grepl("S2", subsite))
s3<-table_results_all %>% filter(grepl("S3", subsite))
s4<-table_results_all %>% filter(grepl("S4", subsite))

myfilename <- paste("S1","fluxes",min(as.Date(s1$start.time)),"to",
                    max(as.Date(s1$start.time)), sep = "_")
write.csv(x = s1, file = paste0(results_path, myfilename,".csv"), row.names = F)

myfilename <- paste("S2","fluxes",min(as.Date(s2$start.time)),"to",
                    max(as.Date(s2$start.time)), sep = "_")
write.csv(x = s2, file = paste0(results_path, myfilename,".csv"), row.names = F)

myfilename <- paste("S3","fluxes",min(as.Date(s3$start.time)),"to",
                    max(as.Date(s3$start.time)), sep = "_")
write.csv(x = s3, file = paste0(results_path, myfilename,".csv"), row.names = F)

myfilename <- paste("S4","fluxes",min(as.Date(s4$start.time)),"to",
                    max(as.Date(s4$start.time)), sep = "_")
write.csv(x = s4, file = paste0(results_path, myfilename,".csv"), row.names = F)

