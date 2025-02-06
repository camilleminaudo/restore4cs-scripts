

# ---
# Authors: Camille Minaudo edits by MIGUEL CABRERA
# Project: "RESTORE4Cs"
# date: "Oct 2023"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
# This script computes fluxes for all incubations present in the fieldsheet for a given sampling season
# Use has to sepcify which sampling has to be processed, and all the rest is done automatically.
#EDITS: 
  #Artefact detector
  #Reduction of "safety margins" to 5s afterstar and 5s beforend



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

#Incubation definition: margins (to be cropped of every incubation) and minimum duration threshold (avoids flux calculation of any incubation shorter than threshold, after removing the margins.
margin_s_after_start<- 5
margin_s_before_end<- 5
minimum_duration_for_flux<- 100

#Shoulder goflux is the oposite of margins: number of seconds to look outside of incubation start-end suplied.
shoulder_goflux<- 0 
#Artefact definition: how many seconds (minimum) of constant slope consistute an incubation with artefact?
artefact_duration_threshold<- 5

# ---- Directories ----
dropbox_root <- "C:/Users/Miguel/Dropbox/RESTORE4Cs - Fieldwork/Data" # You have to make sure this is pointing to the write folder on your local machine
datapath <- paste0(dropbox_root,"/GHG/RAW data")
fieldsheetpath <- paste0(dropbox_root,"/GHG/Fieldsheets")
loggerspath <- paste0(datapath,"/RAW Data Logger")
RData_path <- paste0(dropbox_root,"/GHG/Processed data/RData/")

# results_path<- paste0(dropbox_root,"/GHG/Processed data/computed_flux/")

# results_path <- paste0(dropbox_root, "/GHG/Working data/Incubation_quality/") #Modified folder to save data
results_path <- "C:/Users/Miguel/Dropbox/TEST_quality_raw2flux/"
plots_path <- paste0(results_path,"level_incubation/")



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
    if(abs(fieldsheet_Picarro$unix_start[i] - map_incubations$start[ind])<5*60){
      corresponding_row[i] <- ind
    }
  }
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
  fieldsheet_Picarro$uniqID[ind_NAs]
  # removing rows from fieldsheet_Picarro with missing correspondance
  fieldsheet_Picarro <- fieldsheet_Picarro[-ind_NAs,]
  corresponding_row <- corresponding_row[-ind_NAs]
}

# replacing unix_start and unix_stop with new values
fieldsheet_Picarro$unix_start <- map_incubations$start[corresponding_row]
fieldsheet_Picarro$unix_stop <- map_incubations$stop[corresponding_row]






## ---- Correct LiCOR fieldsheets ----

#Correct Licor start-stop times whenever remarks are available
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



## ---- Merging fieldsheets----
#Import los GATOS fieldsheets: start and stop times directly from fieldsheets (no map_incubations)
fieldsheet_LosGatos <- fieldsheet[fieldsheet$gas_analyzer=="Los Gatos",]

#Combine all fieldsheets with corrected times based on map_incubations when available
fieldsheet <- rbind(fieldsheet_Licor, fieldsheet_LosGatos, fieldsheet_Picarro)


## ----Margins incubations-----
# taking a little margin to avoid critical moments of manipulation with the chamber and manual gas sampling
fieldsheet$unix_start <- fieldsheet$unix_start+margin_s_after_start
fieldsheet$unix_stop <- fieldsheet$unix_stop-margin_s_before_end

# recalculating start and stop in proper formats
fieldsheet$timestamp_start <- as.POSIXct(fieldsheet$unix_start, tz = "UTC", origin = "1970-01-01")
fieldsheet$timestamp_stop <- as.POSIXct(fieldsheet$unix_stop, tz = "UTC", origin = "1970-01-01")

fieldsheet$start_time <- strftime(fieldsheet$timestamp_start, format="%H:%M:%S", tz = 'utc')
fieldsheet$end_time <- strftime(fieldsheet$timestamp_stop, format="%H:%M:%S", tz = 'utc')

fieldsheet <- fieldsheet[order(fieldsheet$subsite),]





#_________________####
#RDATA SAVE (optional) ####
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
    
    import2RData(path = data_folder, instrument = "LI-7810",
                 date.format = "ymd", timezone = 'UTC', keep_all = FALSE, 
                 prec = c(3.5, 0.6, 45))
    
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
    
    # get read of possible duplicated data
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

# ----- 1. Auxtable prep -----

# This section loops over each subsite and performs the flux calculation for CO2 and CH4
subsites <- unique(fieldsheet$subsite)

isF_incub <- T
isFsubsite <- T
for (subsite in subsites){
  message("Now processing ",subsite)
  
  corresp_fs <- fieldsheet[fieldsheet$subsite == subsite,]
  gs <- first(corresp_fs$gas_analyzer)
  
  ## ----1.1. Import looger data----
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
  
  #----1.2. Per-incubation aux variables----
  auxfile <- NULL
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
        
        # --- 1.3. Create auxfile table ----
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
  
# ---- 4. Goflux calculation-------

  if (length(auxfile)>1){
    
    auxfile$Tcham[is.na(auxfile$Tcham)] <- mean(auxfile$Tcham, na.rm = T)
    
    # we only keep incubations longer than X secs (after cropping Xs start and Xs end): see USER OPTIONS
    auxfile <- auxfile[auxfile$duration>minimum_duration_for_flux,]
    
    # we only keep incubations where chamber dimensions are known
    auxfile <- auxfile[!is.na(auxfile$Vtot),] # in case chamber height is not specified in the fieldsheet...
    auxfile <- auxfile[!is.na(auxfile$Area),]
    
    # Define the measurements' window of observation
    # auxfile <- auxfile
    mydata_ow <- obs.win(inputfile = mydata, auxfile = auxfile,
                         obs.length = auxfile$duration, shoulder = shoulder_goflux)
    
    # Join mydata_ow with info on start end incubation
    mydata_auto <- lapply(seq_along(mydata_ow), join_auxfile_with_data.loop, flux.unique = mydata_ow) %>%
      map_df(., ~as.data.frame(.x))
    
    # Additional auxiliary data required for flux calculation.
    mydata_auto <- mydata_auto %>%
      left_join(auxfile %>% select(UniqueID, Area, Vtot, Tcham, Pcham))
    
    # Add instrument precision for each gas
    mydata_auto <- mydata_auto %>%
      mutate(CO2_prec = first(mydata$CO2_prec), CH4_prec = first(mydata$CH4_prec), 
             N2O_prec = first(mydata$N2O_prec), H2O_prec = first(mydata$H2O_prec))
    
    # Calculate fluxes
    CO2_results_auto <- goFlux(dataframe = mydata_auto, gastype = "CO2dry_ppm")
    H2O_results_auto <- goFlux(dataframe = mydata_auto, "H2O_ppm")
    CH4_results_auto <- goFlux(dataframe = mydata_auto, "CH4dry_ppb")
    
    # Use best.flux to select the best flux estimates (LM or HM)
    # based on a list of criteria
    criteria <- c("g.factor", "kappa", "MDF", "R2", "SE.rel")
    
    CO2_flux_res_auto <- best.flux(CO2_results_auto, criteria)
    H2O_flux_res_auto <- best.flux(H2O_results_auto, criteria)
    CH4_flux_res_auto <- best.flux(CH4_results_auto, criteria)
    
    
    #----5. Plot Goflux (optional)----
    if(doPlot){
      # Plots results
      # Make a list of plots of all measurements, for each gastype
      CO2_flux_plots <- flux.plot(CO2_flux_res_auto, mydata_auto, "CO2dry_ppm")
      H2O_flux_plots <- flux.plot(H2O_flux_res_auto, mydata_auto, "H2O_ppm")
      CH4_flux_plots <- flux.plot(CH4_flux_res_auto, mydata_auto, "CH4dry_ppb")
      
      # Combine plot lists into one list
      flux_plot.ls <- c(CO2_flux_plots, CH4_flux_plots, H2O_flux_plots)
      
      # Save plots to pdf in result_subfolder
      myfilename <- paste(subsite, as.character(as.Date(first(auxfile$start.time))),sep="_")
      flux2pdf(flux_plot.ls, outfile = paste0(plots_path,myfilename,".pdf"))
      
    }
    
    #6. Diffusion vs ebullition----
    # estimate ch4 diffusion and ebullition components
    CH4_res_meth1 <- CH4_flux_res_auto  #keep best flux from goflux
    CH4_res_meth1$total_estimated <- NA
    CH4_res_meth1$ebullition <- NA
    CH4_res_meth1$diffusion <- NA
    
    #Only for plots with water, apply estimate of ebullition and diffusion
    for (i in which(auxfile$water_depth>0)){
      if(auxfile$water_depth[i]>0){
        my_incub <- mydata[as.numeric(mydata$POSIX.time)> auxfile$start.time[i] &
                             as.numeric(mydata$POSIX.time)< auxfile$start.time[i]+auxfile$duration[i],]
        my_incub <- my_incub[!is.na(my_incub$CO2dry_ppm),]
        # calling dedicated function
        df_ebull <- separate_ebullition_from_diffusion(my_incub, UniqueID = auxfile$UniqueID[i], doPlot=F)
        # computing fluxes
        H2O_mol = my_incub$H2O_ppm / (1000*1000)
        myfluxterm <- flux.term(auxfile$Vtot[i], auxfile$Pcham[i], auxfile$Area[i],
                                auxfile$Tcham[i], first(H2O_mol))
        CH4_flux_total <- df_ebull$delta_ch4/df_ebull$duration*myfluxterm # nmol/m2/s
        CH4_flux_diff <- df_ebull$avg_diff_slope*myfluxterm # nmol/m2/s
        CH4_flux_ebull <- CH4_flux_total - CH4_flux_diff
      } else {
        #For plots without water, set ebullition to 0 
        CH4_flux_total <- CH4_flux_ebull <- CH4_res_meth1$best.flux[which(CH4_res_meth1$UniqueID==auxfile$UniqueID[i])]
        CH4_flux_ebull <- 0
      }
      CH4_res_meth1$total_estimated[which(CH4_res_meth1$UniqueID==auxfile$UniqueID[i])] <- CH4_flux_total
      CH4_res_meth1$ebullition[which(CH4_res_meth1$UniqueID==auxfile$UniqueID[i])] <- CH4_flux_ebull
      CH4_res_meth1$diffusion[which(CH4_res_meth1$UniqueID==auxfile$UniqueID[i])] <- CH4_flux_diff
    }
    
    CH4_res_meth1$ebullition[which(CH4_res_meth1$ebullition<0)] <- 0
    
    #---- 7. Artefact detection ----
    
    #Initialize artefact columns in CO2 and CH4 fluxes datasets
    CO2_flux_res_auto$contains.artefact <- F
    CH4_res_meth1$contains.artefact <- F
    
    for(i in which(!is.na(auxfile$UniqueID))){
      
      #Select GHG data from UniqueID
      my_incub <- mydata[as.numeric(mydata$POSIX.time)> auxfile$start.time[i] &
                           as.numeric(mydata$POSIX.time)< auxfile$start.time[i]+auxfile$duration[i],]
      # my_incub <- my_incub[!is.na(my_incub$CO2dry_ppm),]
      
      #Function looks for constant slope using second order derivative
      #Check for artefact in CO2 series and flag if found
      CO2_flux_res_auto$contains.artefact[which(CO2_flux_res_auto$UniqueID==auxfile$UniqueID[i])] <- is_constant_slope_any_ghg(my_incub=my_incub, POSIX.time = "POSIX.time", check_cols = c("CO2dry_ppm"), duration = artefact_duration_threshold)
      
      #Check for artefact in CH4 series and flag if found
      CH4_res_meth1$contains.artefact[which(CH4_res_meth1$UniqueID==auxfile$UniqueID[i])] <- is_constant_slope_any_ghg(my_incub=my_incub, POSIX.time = "POSIX.time", check_cols = c("CO2dry_ppm"), duration = artefact_duration_threshold)
      
      
    }
    
    
    
    
    #__________________####
    #SAVE RESULTS####
    #---- 1. Save incubation-level ----
    
    #Save full results in incubation-level path
    myfilenameCO2 <- paste(subsite,"co2_fluxes", as.character(as.Date(first(auxfile$start.time))),sep="_")
    myfilenameCH4 <- paste(subsite,"ch4_fluxes", as.character(as.Date(first(auxfile$start.time))),sep="_")
    write.csv(x = CO2_flux_res_auto, file = paste0(plots_path,myfilenameCO2,".csv"), row.names = F)
    write.csv(x = CH4_res_meth1, file = paste0(plots_path,myfilenameCH4,".csv"), row.names = F)
    
    
    #---- 2. Save Season simple ----
    
    #Join simplified results for both gases in results_path
    table_results <- auxfile %>%
      left_join(CO2_flux_res_auto %>% select(UniqueID, LM.flux, HM.flux, best.flux, model, quality.check, contains.artefact)) %>%
      rename(CO2_LM.flux = LM.flux, CO2_HM.flux = HM.flux, CO2_best.flux = best.flux, CO2_best.model = model,
             CO2_quality.check = quality.check, CO2_contains.artefact=contains.artefact) %>%
      left_join(CH4_res_meth1 %>% select(UniqueID, LM.flux, HM.flux, best.flux, best.flux, model, quality.check, 
                                         diffusion, ebullition, contains.artefact)) %>%
      rename(CH4_LM.flux = LM.flux, CH4_HM.flux = HM.flux, CH4_best.flux = best.flux, CH4_best.model = model, 
             CH4_quality.check = quality.check, CH4_diffusive_flux = diffusion, CH4_ebullitive_flux = ebullition,  CH4_contains.artefact= contains.artefact)
    
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


