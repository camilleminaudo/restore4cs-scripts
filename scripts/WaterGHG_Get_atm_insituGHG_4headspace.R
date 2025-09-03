#GET ATM GHG concentrations for HEADSPACES
#Description----
# ---
# Authors: Miguel Cabrera Brufau
# Project: "RESTORE4Cs"
# date: "August 2025"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
#MOTIVATION:The GC-derived co2 and ch4 concentrations seem a bit too high, we will use the in-situ measured atmospheric concentrations to correct the calibration factors of the GC. We also need the in-situ CO2 and CH4 atm concentrations for the S4 campaign, as the atm exetainers were only analyzed through the N2O licor. 

#FUNCTIONING: This script extracts, for every subsite visit, a single average CO2 and CH4 concentration measured in the field with the portable analyzers. This will be used to correct the GC cal factors and obtain the headspace GHG concentrations. We will create a loop that extracts the non-incubation concentrations. Then apply a conservative filter of extreme datapoints before calculating a croped average CO2 and CH4 atmospheric concentrations. 

#Steps: 

#DONE: 0. Load functions and incubation auxfiles. 
#DONE: 1. Create subsite-level auxfile (to load full day of data with aquaGHG function). Use first start -Xmin and last end +Xmin. 

# LOOPED: for every subsite-sampling
#DONE: 1. Import rawghgdata
#DONE: 2. Define datastart, dataend
#DONE: 3. Filter rawghgdata for useful atm concentrations (take in-between incubation time with 30s margins)
#DONE: 4. ADD data to table: subsite, datetime, Co2, Ch4
#ENDLOOP: 

#TODO: Define data-driven filter (Co2 within 350-550, Ch4 within 1.90-2.50 for example, define hard-limits after brief exploration of data)

#TODO: Save Atm concentration CO2&CH4 for every subsite-sampling. 



rm(list = ls()) # clear workspace
cat("/014") # clear console



#Load PKGs and functions-----
library(lubridate)
library(tidyr)
library(dplyr)
library(pbapply)
library(ggplot2)
library(egg)
library(goFlux)
library(purrr)
library(ggnewscale)
library(stringr)
library(data.table)

#Load repo functions: 
repo_root_r4cs <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
files.sources = list.files(path = paste0(repo_root_r4cs,"/functions"), full.names = T)
for (f in files.sources){source(f)}

#Unclear aquaGHG is needed:
# make sure you clone or download the latest version of aquaGHG: https://github.com/camilleminaudo/aquaGHG/tree/main 
# files.sources = list.files(path = paste0(dirname(repo_root_r4cs),"/aquaGHG/R/"), full.names = T)
# for (f in files.sources){source(f)}

#Load function
load_incubation <- function(auxfile_i, RData_path){
  
  message("Loading data for ",auxfile_i$UniqueID)
  gas <- unique(auxfile_i$gas_analiser)
  
  setwd(RData_path)
  if(gas== "LI-COR"){
    gs_suffix <- "LI-7810"
  } else {
    gs_suffix <- gas
  }
  load(file = paste0(auxfile_i$subsite,"_",gs_suffix,".RData"))
  mydata <- mydata[,c("POSIX.time", 
                      "CO2dry_ppm", "CH4dry_ppb", "H2O_ppm",
                      "CO2_prec",   "CH4_prec",   "H2O_prec"  )]
  
  mydata <- mydata[which(mydata$POSIX.time>=auxfile_i$start.time & 
                           mydata$POSIX.time<=auxfile_i$start.time+auxfile_i$duration),]
  if(dim(mydata)[1]>0){
    mydata$Etime <- as.numeric(mydata$POSIX.time) - min(as.numeric(mydata$POSIX.time))
    mydata$UniqueID <- auxfile_i$UniqueID
  } else {
    warning(paste0("For ", auxfile_i$UniqueID, ", no measurements were found!"))
  }
  return(mydata)
}



# ---- Directories ----
dropbox_root <- "C:/Users/Miguel/Dropbox/RESTORE4Cs - Fieldwork/Data" # You have to make sure this is pointing to the write folder on your local machine

fieldsheetpath <- paste0(dropbox_root,"/GHG/Fieldsheets")#Path with fieldsheets
datapath <- paste0(dropbox_root,"/GHG/RAW data")#Path with raw-data and map_incubations
RData_path <- paste0(dropbox_root,"/GHG/Processed data/RData/")#Path with harmonized Rdata GHG
headspace_path<- paste0(dropbox_root,"/GHG/Headspaces/")

# quality_path<- paste0(dropbox_root, "/GHG/Working data/Incubation_quality/") #Path with quality assessment, crop decisions and flags for each GHG after inspection of each incubation.
# meteo_path<- paste0(dropbox_root, "/Meteo/Formated_data/") #Path with meteo-data
# auxfile_path<- paste0(dropbox_root,"/GHG/Working data/Auxfiles_4aquaGHG/") #path to save auxfiles





#Create incub auxfile ------
#Create full incubation auxfile from fieldsheet data and remarks: NOT CROPPED, we want the negative of the times recorded. 
#Only needed UniqueID, start.time, duration (no meteo, no cropping, no chamber geometry)

#Adapted from createauxfile without any of the cropping or discarding of incubations. 

{
#1. Load and correct fieldsheets-

sampling<- "S" ##pattern to match with all samplings

#Data pre-processing and harmonization: 
# ---- List GHG chamber fieldsheets in Dropbox and read them ---
# list filenames
myfieldsheets_list <- list.files(fieldsheetpath, pattern = "Fieldsheet-GHG.xlsx", all.files = T, full.names = T, recursive = T)
i <- grep(pattern = sampling, x = myfieldsheets_list) # selecting the files corresponding to the selected sampling campaign
myfieldsheets_list <- myfieldsheets_list[i]
# Read all fieldsheets and put them in a single dataframe
fieldsheet <- read_GHG_fieldsheets(myfieldsheets_list)



## ---- Correct Picarro fieldsheets

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


## ---- Correct LiCOR fieldsheets

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


## ---- Merging fieldsheets

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

# ----- 2. Cropping (Omit)

# ----- 3. Auxtable prep

# This section loops over each subsite and produces an auxfile with necesary info
subsites <- unique(fieldsheet$subsite)

isF_incub <- T
isFsubsite <- T
for (subsite in subsites){
  message("Now processing ",subsite)
  
  corresp_fs <- fieldsheet[fieldsheet$subsite == subsite,]
  gs <- first(corresp_fs$gas_analyzer)
  
  setwd(RData_path)
  if(gs== "LI-COR"){
    gs_suffix <- "LI-7810"
  } else {
    gs_suffix <- gs
  }
  
  
  ##----1.1. Per-incubation aux variables
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
        
        # --- 1.2. Create auxfile table
        # An auxfile table, made of fieldsheet. The auxfile
        # requires start.time, duration and UniqueID.
        # start.time must be in the format "%Y-%m-%d %H:%M:%S"
        
        auxfile_tmp <- data.frame(subsite = subsite,
                                  UniqueID = corresp_fs$uniqID[incub],
                                  gas_analiser = gs,
                                  start.time = as.POSIXct((corresp_fs$unix_start[incub]), tz = "UTC"),
                                  duration = (corresp_fs$unix_stop[incub]) - (corresp_fs$unix_start[incub]))
        
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
  
  #Join subsite-auxfiles into a single auxfile
  
  if(subsite==subsites[1]){
    all_auxfile<- auxfile
  }else{
    all_auxfile<- rbind(all_auxfile, auxfile)
  }
}#end of subsite loop


rm(auxfile, auxfile_tmp,corresp_fs, my_incub, mydata,corresponding_row, ind_NAs, ind_noNAs, ind, incub, is_integer, isF_incub, isFsubsite,subsite,unix_start_corr, unix_stop_corr, sampling, subsites, gs, gs_suffix, i, interp)

}

#Create subsite auxfile ----

#HERE WE HAVE all_auxfile with every single incubation start and duration.

#Remove high-tide sampled day before (RI-R2)
wrong_hightides<- c("s1-ri-r2-1-o-d-13:38","s1-ri-r2-2-o-d-13:52","s1-ri-r2-3-o-d-14:05",
                    "s2-ri-r2-1-o-d-12:32","s2-ri-r2-2-o-d-12:47","s2-ri-r2-3-o-d-13:26",
                    "s3-ri-r2-1-o-d-13:44","s3-ri-r2-2-o-d-13:58","s3-ri-r2-3-o-d-14:12",
                    "s4-ri-r2-1-o-d-10:32","s4-ri-r2-2-o-d-10:45","s4-ri-r2-3-o-d-11:00"
)

##Create subsite-auxfile: with 10min margin before and after first and last incubation.
subsite_auxfile<- all_auxfile %>% 
  filter(!UniqueID%in%wrong_hightides) %>%
  group_by(subsite) %>% 
  summarise(subsite=first(subsite),
            UniqueID=subsite,
            gas_analiser=first(gas_analiser),
            first_start.time=first(start.time),
            last_start.time=last(start.time),
            last_duration=last(duration),
            duration_cor= as.numeric(last_start.time+last_duration-first_start.time, units = "secs"),
            #Add margins (15 minutes before sampling start and 15 minutes after it ends)
            start.time=first_start.time- 900,
            duration=as.numeric(duration_cor + 1800)
) %>% select(subsite, UniqueID,gas_analiser, start.time, duration)


#Extract non-incub data------
#LOOP to extract all non-incubation data
for(s in subsite_auxfile$subsite){

#Load data for UniqueID
#Load auxfiles for UniqueID
auxfile_i <- subsite_auxfile[which(subsite_auxfile$UniqueID==s),]
mydata_i <- load_incubation(auxfile_i, RData_path)


#Create filter for records: create time vector of every second of incubations 
incub_timestamps<- all_auxfile %>% 
  filter(subsite==auxfile_i$subsite) %>% 
  mutate(start.time=start.time-5, #5 seconds before start of incubation
         duration=duration + 5) %>% #5 second margin after end of incubation 
  rowwise() %>% 
  mutate(timestamps=list(as.numeric(start.time+seq(0,duration, by=1), units="secs"))) %>% 
  ungroup() %>% 
  summarise(subsite=first(subsite),
            incub_times= list(c(unlist(timestamps)))) %>% pull(incub_times) %>% unlist() %>% round(digits = 0)


mydata_i_notincub<- mydata_i %>% 
  mutate(timestamp=round(as.numeric(POSIX.time), 0)) %>% 
  filter(!timestamp%in%incub_timestamps) %>% 
  select(UniqueID, POSIX.time, timestamp, CO2dry_ppm, CH4dry_ppb) %>%
  rename(subsite=UniqueID)

#Progressively Rbind atm data
if(s==subsite_auxfile$subsite[1]){
  atm_data<- mydata_i_notincub
} else {
  atm_data<- rbind(atm_data, mydata_i_notincub)
}

}


#Filter and summarise -----
#Summarise atm data with 10% trim (removing highest and lowest 5%)

allsummary_atm_data<- atm_data %>%
  rename(co2_ppm=CO2dry_ppm, ch4_ppb=CH4dry_ppb) %>% 
  pivot_longer(cols = c(co2_ppm, ch4_ppb), names_to = "ghgspecies",values_to = "conc") %>% 
  group_by(subsite, ghgspecies) %>%
  filter(!is.na(conc)) %>% 
  summarise(avg=mean(conc),
            desv=sd(conc), 
            median=median(conc),
            min_=min(conc),
            max_=max(conc),
            nobs=n())




summary_atm_data<- atm_data %>%
  rename(co2_ppm=CO2dry_ppm, ch4_ppb=CH4dry_ppb) %>% 
  pivot_longer(cols = c(co2_ppm, ch4_ppb), names_to = "ghgspecies",values_to = "conc") %>% 
  group_by(subsite, ghgspecies) %>%
  filter(!is.na(conc)) %>% 
  filter(!(ghgspecies=="co2_ppm"&(conc>440|conc<380))) %>% #Data-limits for atm CO2 range
  filter(!(ghgspecies=="ch4_ppb"&(conc>2500|conc<1900))) %>% #Data-limits for atm CH4 range
  # filter(between(conc, quantile(conc, 0.05), quantile(conc,0.95))) %>% 
  summarise(avg=mean(conc),
            desv=sd(conc), 
            median=median(conc),
            min_=min(conc),
            max_=max(conc),
            nobs=n())
            
#Save atm conc -------

write.csv(summary_atm_data, file = paste0(headspace_path,"in_situ_atm_co2&ch4.csv"),row.names = F)








