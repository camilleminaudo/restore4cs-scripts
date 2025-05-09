
# ---
# Authors: Miguel Cabrera, Camille Minaudo
# Project: "RESTORE4Cs"
# date: "April 2023"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
# This script creates auxfiles needed to compute in-situ fluxes of CO2 and CH4 with the aquaGHG package from Camille https://github.com/camilleminaudo/aquaGHG

#Two auxfiles (one for CO2, one for CH4) are created using:
# Fieldsheets + map_incubations to set correct start-stop incubations
# Per-GHG cropping decisions logged in excel file Inspection_table_allincubations_tocrop_perGHG.xlsx
# Meteo-data harmonized for each subsite visit. 

#Outputs two csv files (CO2 and CH4) with the following columns: 

#subsite
#UniqueID
#gas_analiser 
#start.time (posixct utc)
#duration (seconds)
#water_depth (cm)
#Area (cm2)
#Vtot (L)
#Tcham (ºC)
#Pcham (kPa)

#strata
#chamberType
#lightCondition

#The columns Start.time and duration might be different for each gas due to ghg-dedicated cropping of artefacts. 



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




# ---- Directories ----
dropbox_root <- "C:/Users/Miguel/Dropbox/RESTORE4Cs - Fieldwork/Data" # You have to make sure this is pointing to the write folder on your local machine

fieldsheetpath <- paste0(dropbox_root,"/GHG/Fieldsheets")#Path with fieldsheets
datapath <- paste0(dropbox_root,"/GHG/RAW data")#Path with raw-data and map_incubations
RData_path <- paste0(dropbox_root,"/GHG/Processed data/RData/")#Path with harmonized Rdata GHG
quality_path<- paste0(dropbox_root, "/GHG/Working data/Incubation_quality/") #Path with quality assessment, crop decisions and flags for each GHG after inspection of each incubation.
meteo_path<- paste0(dropbox_root, "/Meteo/Formated_data/") #Path with meteo-data
auxfile_path<- paste0(dropbox_root,"/GHG/Working data/Auxfiles_4aquaGHG/") #path to save auxfiles


#1. Load and correct fieldsheets-----
minimum_duration_for_flux <- 100 #Minimum duration of incubation to be included in auxfile
sampling<- "S" ##pattern to match with all samplings

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
#the following incubations should not be found (not map_incubation for them):
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
  message(paste("Could not find remarks to correct times for ", dim(fieldsheet_Licor)[1]-sum(!is.na(ind_noNAs)), "incubations" ))
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

#delete intermediate fieldsheets
rm(fieldsheet_Licor, fieldsheet_Licormapped, fieldsheet_LosGatos, fieldsheet_Picarro, map_incubations)




#----2. Load crop & Meteo-----

#Import table with decisions and cropping amounts: all cropping is based on last run of 5-s margin script, so we have to re-add this 5 seconds to the cropping margins
croping<- read_xlsx(paste0(quality_path, "Inspection_table_allincubations_tocrop_perGHG.xlsx"),na = "NA")

#Simplify and add the additional 5s of cropping to all crop-start and crop-end: Inspection and cropping decisions were taken on fluxes calculated based on 5-s margin. WE re-add this 5s of margin here for all incubations (custom-cropped or not) 
croping_simple<- croping %>%
  select(UniqueID, 
         co2_start_crop_s, co2_end_crop_s,
         ch4_start_crop_s, ch4_end_crop_s) %>%
  #add 5s crop for all h2o
  mutate(h2o_start_crop_s=5, h2o_end_crop_s=5) %>% 
  # add 5s to all crop decissions for co2 and ch4
  mutate(co2_start_crop_s=case_when(is.na(co2_start_crop_s)~5,
                                    TRUE~co2_start_crop_s+5),
         co2_end_crop_s=case_when(is.na(co2_end_crop_s)~5,
                                  TRUE~co2_end_crop_s+5),
         ch4_start_crop_s=case_when(is.na(ch4_start_crop_s)~5,
                                    TRUE~ch4_start_crop_s+5),
         ch4_end_crop_s=case_when(is.na(ch4_end_crop_s)~5,
                                  TRUE~ch4_end_crop_s+5))

#croping_simple is used to adapt the auxfile for each gas in the section Goflux calculation
rm(croping)


#Load data for all meteostations and correspondence of subsite-meteostation (used to subset the correct data)
meteocorresp <- read_csv(paste0(meteo_path, "correspondence_subsite-meteostation.csv"),show_col_types = F)

meteodata <- read_csv(paste0(meteo_path,"allmeteostations_onlysamplingdays.csv"),show_col_types = F)




# ----- 3. Auxtable prep -----

# This section loops over each subsite and produces an auxfile with necesary info for flux calculation
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
  
  #Check to only create auxfile for incubations with actual data
  if(file.exists(paste0(subsite,"_",gs_suffix,".RData"))){
    load(file = paste0(subsite,"_",gs_suffix,".RData"))
    
    #Loop over each incubation: 
    for (incub in seq_along(corresp_fs$plot_id)){
      my_incub <- mydata[as.numeric(mydata$POSIX.time)> corresp_fs$unix_start[incub] &
                           as.numeric(mydata$POSIX.time)< corresp_fs$unix_stop[incub],]
      my_incub <- my_incub[!is.na(my_incub$CO2dry_ppm),]
      
      if (dim(my_incub)[1]>0){
        
        # Compute Temperature (ºC) and pressure (kPa) for the initial moment of incub using the meteo data (using approx to interpolate from hourly data).
        #Temperature and Pressure are obtained at the start.time of chamber deployment, irregardless of potential subsequent cropping of GHG timeseries, which is the correct approach (flux term depends on the deployment of chamber and it is the same for the whole incubation)
        
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
    }#End of incub loop
    
  } else {
    message("---> Could not find corresponding ",gs," data")
  }
  
  
  
  if (length(auxfile)>1){
    
    auxfile$Tcham[is.na(auxfile$Tcham)] <- mean(auxfile$Tcham, na.rm = T)
    
    # we only keep incubations longer than X secs (after cropping Xs start and Xs end): see USER OPTIONS
    auxfile <- auxfile[auxfile$duration>minimum_duration_for_flux,]
    
    # we only keep incubations where chamber dimensions are known
    auxfile <- auxfile[!is.na(auxfile$Vtot),] # in case chamber height is not specified in the fieldsheet...
    auxfile <- auxfile[!is.na(auxfile$Area),]
    
    
    #3. Per-GHG crop------    
    
    #Create 2 separate auxfiles, one per GHG, with modified start.time and duration based on cropping decissions per GHG.
    
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
    
  }#End of check for length auxfile  
  
  #Join subsite-auxfiles into a single auxfile for each gas
  
  if(subsite==subsites[1]){
    all_auxfile_co2<- auxfile_co2
    all_auxfile_ch4<- auxfile_ch4
  }else{
    all_auxfile_co2<- rbind(all_auxfile_co2, auxfile_co2)
    all_auxfile_ch4<- rbind(all_auxfile_ch4, auxfile_ch4)
  }
}#end of subsite loop



#Check fiedlsheet not in auxfile
mising<- fieldsheet %>% filter(!uniqID%in%all_auxfile_ch4$UniqueID)
mising %>% select(uniqID,comments)
#Good, only missing incubations without actual gas-analyzer data

#4. Save AUXFILES----

write.csv(all_auxfile_co2, file = paste0(auxfile_path,"co2_auxfile.csv"), row.names = F)
write.csv(all_auxfile_ch4, file = paste0(auxfile_path,"ch4_auxfile.csv"), row.names = F)



