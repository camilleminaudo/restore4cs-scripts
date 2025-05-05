

# ---
# Authors: Camille Minaudo edits by MIGUEL CABRERA (April 2025)
# Project: "RESTORE4Cs"
# date: "Oct 2023"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
# This script is adapted to produce plots for CO2, CH4 and H2O (using goflux) for S4 vegetated and bare plots with transparent and dark incubations. This will be used to check ventilation events between the transparent and dark incubation and decide for which dark chambers we can calculate N2O fluxes (discard flux N2O flux from dark if no ventilation took place between transparent and dark chamber). We cannot justify an incubation transparent+dark, from a design standpoint. Confirm ventilation or discard dark N2O fluxes. 

#This script was adapted from the typical raw2flux procedure, modifying the fieldsheets to obtain only transparent-dark incubation pairs (same plot, 2 incubations) for S4 sampling campaign. We produce auxfiles from the modified fieldsheets (uniqueID from transparent incubation is kept, duration is calcualted as dark_end minus transparent start, without margins)

#THIS sript does not produce csv files only pdfs to check ventilations.

#For plot S4-VA-A2-5 bare, no ventilation plot  (transparent incubation performed after dark one), but checked manually that chamber was ventilated in between incubations. 

#After inspection of plots, ventilation event are recorded as logical "ventilated" (T/F) in Dark_incubation_ventilationcheck.xlsx inside N2O_fluxes folder for every "UniqueID_notime"

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
sampling <- "S4" 

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


results_path<- "C:/Users/Miguel/Dropbox/TEST_quality_raw2flux/Ventilation_test/"
plots_path <- results_path



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
  
  # removing rows from fieldsheet_Picarro with missing correspondance
  fieldsheet_Picarro <- fieldsheet_Picarro[-ind_NAs,]
  corresponding_row <- corresponding_row[-ind_NAs]
}

# replacing unix_start and unix_stop with new values
fieldsheet_Picarro$unix_start <- map_incubations$start[corresponding_row]
fieldsheet_Picarro$unix_stop <- map_incubations$stop[corresponding_row]


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

# recalculating start and stop in proper formats
fieldsheet$timestamp_start <- as.POSIXct(fieldsheet$unix_start, tz = "UTC", origin = "1970-01-01")
fieldsheet$timestamp_stop <- as.POSIXct(fieldsheet$unix_stop, tz = "UTC", origin = "1970-01-01")

fieldsheet$start_time <- strftime(fieldsheet$timestamp_start, format="%H:%M:%S", tz = 'utc')
fieldsheet$end_time <- strftime(fieldsheet$timestamp_stop, format="%H:%M:%S", tz = 'utc')

fieldsheet <- fieldsheet[order(fieldsheet$subsite),]


#-----ADAPT for vent----

#filter for plots with transparent and dark 
fieldsheet_td<- fieldsheet %>% 
  mutate(plotuniqueID= paste0(subsite,"-",plot_id)) %>% 
  group_by(plotuniqueID) %>% 
  filter(n()==2) %>%  #get only plots with 2 incubations 
#Subsitute transparent unix_stop wiht dark unxistop
  mutate(dark_unixstop= if_else(transparent_dark=="dark",unix_stop,NA_real_),
         unix_stop=mean(dark_unixstop,na.rm=T)) %>% 
  #Remove dark incubations from fieldsheet
  filter(transparent_dark=="transparent")


fieldsheet<- fieldsheet_td

#_________________####
#FLUX CALCULATION####


#----0. Load crop & Meteo-----

#NO crop
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
    
    #Correct the units for picarro H2O,import function assumes mmol/mol, actual units of our instrument are % (i.e. 100* mol/mol). We correct the Rdata by multiplying the concentration by 10. 
    if (gs=="Picarro"){mydata$H2O_ppm<- mydata$H2O_ppm*10}
    
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
    #DO NOT CROP anything, just split auxfile into one per GHG (to avoid modifying code too much)
    auxfile_co2<- auxfile
    auxfile_ch4<- auxfile
    auxfile_h2o<- auxfile
    
    #DO NOT RUN:
    # auxfile_co2<- auxfile %>% 
    #   left_join(croping_simple %>% select(UniqueID,co2_start_crop_s,co2_end_crop_s), by="UniqueID") %>% 
    #   mutate(start.time=start.time+co2_start_crop_s,
    #          duration=duration-(co2_start_crop_s+co2_end_crop_s)) %>% 
    #   select(-c(co2_start_crop_s, co2_end_crop_s))
    # 
    # auxfile_ch4<- auxfile %>% 
    #   left_join(croping_simple %>% select(UniqueID,ch4_start_crop_s,ch4_end_crop_s), by="UniqueID") %>% 
    #   mutate(start.time=start.time+ch4_start_crop_s,
    #          duration=duration-(ch4_start_crop_s+ch4_end_crop_s)) %>% 
    #   select(-c(ch4_start_crop_s, ch4_end_crop_s))
    # 
    # auxfile_h2o<- auxfile %>% 
    #   left_join(croping_simple %>% select(UniqueID,h2o_start_crop_s,h2o_end_crop_s), by="UniqueID") %>% 
    #   mutate(start.time=start.time+h2o_start_crop_s,
    #          duration=duration-(h2o_start_crop_s+h2o_end_crop_s)) %>% 
    #   select(-c(h2o_start_crop_s, h2o_end_crop_s))
    # 
    
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
    
      #g.factor: arbitrary limit for HM/LM, defaults to 2 (if non-linearity is expected, this is not a good selection criteria). If threshold is exceeded, LM is chosen. Decission: DO NOT USE
    
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
    

}
}
