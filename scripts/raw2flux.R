
# ---
# Authors: Camille Minaudo
# Project: "RESTORE4Cs"
# date: "Oct 2023"
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
sampling <- "S2"
# USER, please specify if you want plots to be saved
harmonize2RData <- F
doPlot <- T
#################################



# ---- Directories ----
dropbox_root <- "C:/Users/Camille Minaudo/Dropbox/RESTORE4Cs - Fieldwork/Data" # You have to make sure this is pointing to the write folder on your local machine
datapath <- paste0(dropbox_root,"/GHG/RAW data")
fieldsheetpath <- paste0(dropbox_root,"/GHG/Fieldsheets")
corrfieldsheetpath <- paste0(dropbox_root,"/GHG/Processed data/corrected_fieldsheets")
loggerspath <- paste0(datapath,"/RAW Data Logger")
plots_path <- paste0(dropbox_root,"/GHG/Processed data/plots_all_incubations/")
RData_path <- paste0(dropbox_root,"/GHG/Processed data/RData/")
results_path <- paste0(dropbox_root,"/GHG/Processed data/computed_flux/")


# ----- Data pre-processing and harmonization -----

# ---- List GHG chamber fieldsheets in Dropbox and read them ---
# list filenames
myfieldsheets_list <- list.files(fieldsheetpath, pattern = "Fieldsheet-GHG.xlsx", all.files = T, full.names = T, recursive = T)
i <- grep(pattern = sampling, x = myfieldsheets_list) # selecting the files corresponding to the selected sampling campaign
myfieldsheets_list <- myfieldsheets_list[i]
# Read all fieldsheets and put them in a single dataframe
fieldsheet <- read_GHG_fieldsheets(myfieldsheets_list)

fieldsheet$uniqID <- tolower(paste(fieldsheet$subsite,fieldsheet$plot_id,substr(fieldsheet$transparent_dark, 1, 1), sep = "-"))



# ---- Correct fieldsheets in the case of Picarro data ---

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

# recalculating start and stop in propre formats
fieldsheet$timestamp_start <- as.POSIXct(fieldsheet$unix_start, tz = "UTC", origin = "1970-01-01")
fieldsheet$timestamp_stop <- as.POSIXct(fieldsheet$unix_stop, tz = "UTC", origin = "1970-01-01")

fieldsheet$start_time <- strftime(fieldsheet$timestamp_start, format="%H:%M:%S", tz = 'utc')
fieldsheet$end_time <- strftime(fieldsheet$timestamp_stop, format="%H:%M:%S", tz = 'utc')

fieldsheet <- fieldsheet[order(fieldsheet$subsite),]



# ---- Import and store measurements to RData ----
if(harmonize2RData){
  data_folders <- list.dirs(datapath, full.names = T, recursive = T)[-1]
  i <- grep(pattern = sampling, x = data_folders) # selecting the files corresponding to the selected sampling campaign
  data_folders <- data_folders[i]
  
  r <- grep(pattern = "RData",x=data_folders)
  if(length(r)>0){data_folders <- data_folders[-r]}
  
  message("Here is the list of data folders in here:")
  print(data_folders)
  
  
  # Import and store data for for Picarro and LosGatos data
  
  raw2RData_P_LG <- function(data_folders, instrument, instrumentID, date.format){
    r <- grep(pattern = instrument, x=data_folders)
    for (data_folder in data_folders[r]){
      setwd(data_folder)
      subsite = basename(data_folder)
      message(paste0("processing folder ",basename(data_folder)))
      import2RData(path = data_folder, instrument = instrumentID,
                   date.format = date.format, timezone = 'UTC')
      
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
      mydata <- mydata_imp[!is_duplicate,]
      
      setwd(RData_path)
      
      # save this dataframe as a new RData file
      save(mydata, file = paste0(subsite,"_",instrument,".RData"))
    }
  }
  
  raw2RData_P_LG(data_folders, instrument = "Picarro", instrumentID = "G2508", date.format = "ymd")
  raw2RData_P_LG(data_folders, instrument = "Los Gatos", instrumentID = "LGR", date.format = "mdy")
  
  # Import and store data for LiCOR data
  fieldsheet_Licor <- fieldsheet[fieldsheet$gas_analyzer=="LI-COR",]
  list_subsites_Licor <- unique(fieldsheet_Licor$subsite)
  
  r_licor <- grep(pattern = "Licor",x=data_folders)
  for (data_folder in data_folders[r_licor]){
    setwd(data_folder)
    message(paste0("processing folder ",basename(data_folder)))
    
    import2RData(path = data_folder, instrument = "LI-7810",
                 date.format = "ymd", timezone = 'UTC')
    
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
      corresponding_date <- unique(fieldsheet_Licor$date[fieldsheet_Licor$subsite == subsite])
      
      # check if we already have this data somewhere in the clean file
      n_lines_d <- length(mydata_imp$DATE[mydata_imp$DATE == corresponding_date])
      if(n_lines_d > 0){
        message(paste0("... there is some data to store for subsite ",subsite))
        mydata <- mydata_imp[mydata_imp$DATE == corresponding_date,]
        
        setwd(RData_path)
        # save this dataframe as a new RData file
        save(mydata, file = paste0(subsite,"_","LI-7810",".RData"))
      }
    }
  }
}




# ----- Flux calculation -----

# For each subsite in fieldsheet, go through each incubation and compute co2 and ch4 fluxes
if (sampling =="S2"){
  subsites <- unique(fieldsheet$subsite)[-c(13,15)] # remove [-c(13,15)] as soon as Benj fixed the issue with LosGatos S2-DA data
} else {
  subsites <- unique(fieldsheet$subsite)
}
isF_incub <- T
isFsubsite <- T
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
  
  
  
  setwd(RData_path)
  if(gs== "LI-COR"){
    gs_suffix <- "LI-7810"
  } else {
    gs_suffix <- gs
  }
  
  
  
  load(file = paste0(subsite,"_",gs_suffix,".RData"))
  
  auxfile <- NULL
  
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
      
      UniqueID = paste(subsite, seq_along(corresp_fs$pilot_site),corresp_fs$strata,substr(corresp_fs$transparent_dark, 1, 1),sep = "-")[incub]
      
      auxfile_tmp <- data.frame(subsite = subsite,
                                UniqueID = gsub(" ", "", UniqueID, fixed = TRUE),
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
  
  auxfile$Tcham[is.na(auxfile$Tcham)] <- mean(auxfile$Tcham, na.rm = T)
  
  # we only keep incubations longer than 100 secs
  auxfile <- auxfile[auxfile$duration>100,]
  # we only keep incubations where chamber dimensions are known
  auxfile <- auxfile[!is.na(auxfile$Vtot),] # in case chamber height is not specified in the fieldsheet...
  
  # Define the measurements' window of observation
  # auxfile <- auxfile
  mydata_ow <- obs.win(inputfile = mydata, auxfile = auxfile,
                       obs.length = auxfile$duration, shoulder = 2)
  
  # Join mydata_ow with info on start end incubation
  mydata_auto <- lapply(seq_along(mydata_ow), join_auxfile_with_data.loop, flux.unique = mydata_ow) %>%
    map_df(., ~as.data.frame(.x))
  
  # Additional auxiliary data required for flux calculation.
  mydata_auto <- mydata_auto %>%
    left_join(auxfile %>% select(UniqueID, Area, Vtot, Tcham, Pcham))
  
  # Add instrument precision for each gas
  prec = c(3.5, 0.6, 0.4, 45, 45)
  mydata_auto <- mydata_auto %>%
    mutate(CO2_prec = prec[1], CH4_prec = prec[2], N2O_prec = prec[3],
           H2O_prec = prec[4])
  
  # Calculate fluxes
  CO2_results_auto <- goFlux(dataframe = mydata_auto, gastype = "CO2dry_ppm")
  H2O_results_auto <- goFlux(mydata_auto, "H2O_ppm")
  CH4_results_auto <- goFlux(mydata_auto, "CH4dry_ppb")
  
  # Use best.flux to select the best flux estimates (LM or HM)
  # based on a list of criteria
  criteria <- c("g.factor", "kappa", "MDF", "R2", "SE.rel")
  
  CO2_flux_res_auto <- best.flux(CO2_results_auto, criteria)
  H2O_flux_res_auto <- best.flux(H2O_results_auto, criteria)
  CH4_flux_res_auto <- best.flux(CH4_results_auto, criteria)
  
  if(doPlot){
    # Plots results
    # Make a list of plots of all measurements, for each gastype
    CO2_flux_plots <- flux.plot(CO2_flux_res_auto, mydata_auto, "CO2dry_ppm")
    H2O_flux_plots <- flux.plot(H2O_flux_res_auto, mydata_auto, "H2O_ppm")
    CH4_flux_plots <- flux.plot(CH4_flux_res_auto, mydata_auto, "CH4dry_ppb")
    
    # Combine plot lists into one list
    flux_plot.ls <- c(CO2_flux_plots, CH4_flux_plots, H2O_flux_plots)
    
    # Save plots to pdf
    myfilename <- paste(subsite, as.character(as.Date(first(auxfile$start.time))),sep="_")
    flux2pdf(flux_plot.ls, outfile = paste0(results_path,"/level_incubation/",myfilename,".pdf"))
    
  }
  
  
  # estimate ch4 diffusion and ebullition components---------| METHOD 1 |-----------
  CH4_res_meth1 <- CH4_flux_res_auto
  CH4_res_meth1$total_estimated <- NA
  CH4_res_meth1$ebullition <- NA
  CH4_res_meth1$diffusion <- NA
  
  
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
      CH4_flux_total <- CH4_flux_ebull <- CH4_res_meth1$best.flux[which(CH4_res_meth1$UniqueID==auxfile$UniqueID[i])]
      CH4_flux_ebull <- 0
    }
    CH4_res_meth1$total_estimated[which(CH4_res_meth1$UniqueID==auxfile$UniqueID[i])] <- CH4_flux_total
    CH4_res_meth1$ebullition[which(CH4_res_meth1$UniqueID==auxfile$UniqueID[i])] <- CH4_flux_ebull
    CH4_res_meth1$diffusion[which(CH4_res_meth1$UniqueID==auxfile$UniqueID[i])] <- CH4_flux_diff
  }
  
  CH4_res_meth1$ebullition[which(CH4_res_meth1$ebullition<0)] <- 0
  
  
  setwd(paste0(results_path,"/level_incubation"))
  myfilenameCO2 <- paste(subsite,"co2_fluxes", as.character(as.Date(first(auxfile$start.time))),sep="_")
  myfilenameCH4 <- paste(subsite,"ch4_fluxes", as.character(as.Date(first(auxfile$start.time))),sep="_")
  write.csv(x = CO2_flux_res_auto, file = paste0(myfilenameCO2,".csv"), row.names = F)
  write.csv(x = CH4_res_meth1, file = paste0(myfilenameCH4,".csv"), row.names = F)
  
  
  table_results <- auxfile %>%
    left_join(CO2_flux_res_auto %>% select(UniqueID, LM.flux, HM.flux, best.flux, model, quality.check)) %>%
    rename(CO2_LM.flux = LM.flux, CO2_HM.flux = HM.flux, CO2_best.flux = best.flux, CO2_best.model = model,
           CO2_quality.check = quality.check) %>%
    left_join(CH4_res_meth1 %>% select(UniqueID, LM.flux, HM.flux, best.flux, best.flux, model, quality.check, 
                                       diffusion, ebullition)) %>%
    rename(CH4_LM.flux = LM.flux, CH4_HM.flux = HM.flux, CH4_best.flux = best.flux, CH4_best.model = model, 
           CH4_quality.check = quality.check, CH4_diffusive_flux = diffusion, CH4_ebullitive_flux = ebullition)
  
  if (isFsubsite){
    isFsubsite <- F
    table_results_all <- table_results
  } else {
    table_results_all <- rbind(table_results_all, table_results)
  }
}

setwd(results_path)
myfilename <- paste(sampling,"fluxes",min(as.Date(table_results_all$start.time)),"to",
                    max(as.Date(table_results_all$start.time)), sep = "_")
write.csv(x = table_results_all, file = paste0(myfilename,".csv"), row.names = F)



# ---- Some plots ----


table_results_all$campaign_site <- substr(table_results_all$subsite,start = 1, stop = 5)
table_results_all$subsite_short <- substr(table_results_all$subsite,start = 7, stop = 8)


plt_CO2 <- ggplot(table_results_all, aes(subsite_short, CO2_LM.flux,
                                         fill = lightCondition))+
  geom_hline(yintercept = 0)+
  geom_boxplot(alpha=0.5)+
  # geom_jitter(width = 0.2, aes(colour = strata), size=2)+
  theme_article()+
  xlab("subsite")+
  ylab("CO2 flux mmol/m2/s")+
  ggtitle("CO2 flux")+
  scale_fill_viridis_d(begin = 0.2, end = 0.9)+
  scale_colour_viridis_d(begin = 0.2, end = 0.9, option = "C")+
  facet_grid(.~campaign_site)


plt_CH4diff <- ggplot(table_results_all, aes(subsite_short, CH4_LM.flux, 
                                             fill = lightCondition))+
  geom_hline(yintercept = 0)+
  geom_boxplot(alpha=0.5)+
  # geom_jitter(width = 0.2, aes(colour = strata), size=2)+
  theme_article()+
  xlab("subsite")+
  ylab("CH4 flux nmol/m2/s")+
  ggtitle("CH4 flux")+
  scale_fill_viridis_d(begin = 0.2, end = 0.9)+
  scale_colour_viridis_d(begin = 0.2, end = 0.9, option = "C")+
  facet_grid(.~campaign_site)

plt_all <- ggarrange(plt_CO2, plt_CH4diff, ncol = 1)

ggsave(plot = plt_all, filename = paste0(myfilename,".jpg"), path = results_path, 
       width = 10, height = 8, dpi = 300, units = 'in')



for (cs in unique(table_results_all$campaign_site)){
  
  table_results_cs <- table_results_all[table_results_all$campaign_site == cs,]
  
  plt_CO2 <- ggplot(table_results_cs, aes(subsite_short, CO2_LM.flux,
                                          fill = lightCondition))+
    geom_hline(yintercept = 0)+
    geom_boxplot(alpha=0.5, width=0.5)+
    geom_jitter(width = 0.1, size=2, alpha=0.5, aes(shape = strata))+
    theme_article()+
    xlab("subsite")+
    ylab("CO2 flux mmol/m2/s")+
    ggtitle(paste0(cs, ""))+
    scale_fill_viridis_d(begin = 0.2, end = 0.9)+
    scale_colour_viridis_d(begin = 0.2, end = 0.9, option = "C")+
    # facet_grid(lightCondition~.)+
    theme(legend.position ='none')
  
  
  plt_CH4diff <- ggplot(table_results_cs, aes(subsite_short, CH4_LM.flux, 
                                              fill = lightCondition))+
    geom_hline(yintercept = 0)+
    geom_boxplot(alpha=0.5, width=0.5)+
    geom_jitter(width = 0.1, size=2, alpha=0.5, aes(shape = strata))+
    theme_article()+
    xlab("subsite")+
    ylab("CH4 flux nmol/m2/s")+
    # ggtitle("CH4 flux")+
    scale_fill_viridis_d(begin = 0.2, end = 0.9)+
    scale_colour_viridis_d(begin = 0.2, end = 0.9, option = "C")
  
  plt_cs <- ggarrange(plt_CO2, plt_CH4diff, ncol = 2)
  
  
  setwd(results_path)
  myfilename <- paste(cs,"fluxes",min(as.Date(table_results_all$start.time)),"to",
                      max(as.Date(table_results_all$start.time)), sep = "_")
  
  ggsave(plot = plt_cs, filename = paste0(myfilename,".jpg"), path = results_path, 
         width = 10, height = 5, dpi = 300, units = 'in')
}


# 
# 
# 
# ggplot(table_results_all, aes(CO2_LM.flux, CO2_HM.flux, colour = strata))+
#   geom_abline(slope = 1, intercept = 0)+
#   geom_point(size=2, alpha=0.5)+
#   theme_article()+
#   scale_colour_viridis_d(begin = 0.2, end = 0.9)#+facet_grid(.~subsite)
