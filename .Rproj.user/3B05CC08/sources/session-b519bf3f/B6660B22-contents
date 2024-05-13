
# ---
# Authors: Camille Minaudo
# Project: "RESTORE4Cs"
# date: "Janv 2024"
# https://github.com/camilleminaudo/picarroprocessing
# ---

# --- Description of this script
# This script transforms raw data of the Picarro gas scouter analyser into fluxes


rm(list = ls()) # clear workspace
cat("/014") # clear console


# ---- packages ----
library(tidyverse)
library(readxl)
library(lubridate)
library(zoo)
library(ggplot2)
library(grid)
library(egg)
library(GoFluxYourself)
require(dplyr)
require(purrr)
require(msm)
require(pbapply)


repo_root <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
files.sources = list.files(path = paste0(repo_root,"/functions"), full.names = T)
for (f in files.sources){source(f)}



# ---- Directories ----
datapath <- paste0(repo_root,"/data")
rawdatapath <- paste0(datapath,"/raw_data/RAW Data Picarro")
fieldsheetpath <- paste0(datapath,"/fieldsheets")
# corrfieldsheetpath <- paste0(datapath,"/GHG/Processed data/corrected_fieldsheets")
loggerspath <- paste0(datapath,"/raw_data/RAW Data Logger")
plots_path <- paste0(repo_root,"/results/plots")
results_path <- paste0(repo_root,"/results/processed")



# ---- List GHG chamber fieldsheets in data folder and read them ---
# list filenames
myfieldsheets_list <- list.files(fieldsheetpath, pattern = "-Fieldsheet-GHG.xlsx", all.files = T, full.names = T, recursive = T)
# Read all fieldsheets and put them in a single dataframe
fieldsheet <- read_GHG_fieldsheets(myfieldsheets_list)




# ---- Import and store measurements to RData ----
data_folders <- list.dirs(rawdatapath, full.names = T, recursive = T)[-1]
r <- grep(pattern = "RData",x=data_folders)
if(length(r)>0){data_folders <- data_folders[-r]}


message("Here is the list of data folders in here:")
print(data_folders)


for (data_folder in data_folders){
  setwd(data_folder)
  subsite = basename(data_folder)
  message(paste0("processing folder ",basename(data_folder)))
  import2RData(path = data_folder, instrument = "G2508",
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
  mydata <- mydata_imp[!is_duplicate,]

  # creating a folder where to put this data
  dir.create(file.path(results_path, subsite))
  setwd(file.path(results_path, subsite))

  # save this dataframe as a new RData file
  save(mydata, file = paste0("data_",subsite,".RData"))
}



# ---- Correct fieldsheet  times based on "incubation map" ----

read_map_incubations <- function(path2folder){

  my_maps_filenames <- list.files(path2folder, pattern = ".csv", all.files = T, full.names = T, recursive = T)
  isF <- T
  for(my_maps_filename in my_maps_filenames){

    map_incubations_temp <- read.csv(file = my_maps_filename,
                                     header = T)
    map_incubations_temp <- map_incubations_temp[map_incubations_temp$Species ==  "CH4",]
    map_incubations_temp$subsite <- basename(dirname(my_maps_filename))

    if(isF){
      isF <- F
      map_incubations <- map_incubations_temp
    } else {
      map_incubations <- rbind(map_incubations, map_incubations_temp)
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
map_incubations <- read_map_incubations(path2folder = rawdatapath)


# check the closest incubation in map_incubations for each row in fieldsheet.
# if more than 3 minutes apart, we consider the row in map_incubations out of sampling
corresponding_row <- NA*fieldsheet$plot_id

for (i in seq_along(fieldsheet$plot_id)){
  ind <- which.min(abs(fieldsheet$unix_start[i] - map_incubations$start))
  if(abs(fieldsheet$unix_start[i] - map_incubations$start[ind])<3*60){corresponding_row[i] <- ind}
}
corresponding_row <- corresponding_row[!is.na(corresponding_row)]
map_incubations <- map_incubations[corresponding_row, ]

if(dim(map_incubations)[1] != dim(fieldsheet)[1]){
  fieldsheet <- fieldsheet[seq_along(corresponding_row),]
}

fieldsheet$start_time <- format(as.POSIXct(fieldsheet$start_time, tz = 'utc', format = "%T"), "%H:%M:%S")
fieldsheet$end_time <- format(as.POSIXct(fieldsheet$end_time, tz = 'utc', format = "%T"), "%H:%M:%S")

fieldsheet$unix_start_corrected <- map_incubations$start
fieldsheet$unix_stop_corrected <- map_incubations$stop

fieldsheet$timestamp_start <- as.POSIXct(fieldsheet$unix_start_corrected, tz = "UTC", origin = "1970-01-01")
fieldsheet$timestamp_stop <- as.POSIXct(fieldsheet$unix_stop_corrected, tz = "UTC", origin = "1970-01-01")

fieldsheet$start_time_corrected <- strftime(fieldsheet$timestamp_start, format="%H:%M:%S", tz = 'utc')
fieldsheet$end_time_corrected <- strftime(fieldsheet$timestamp_stop, format="%H:%M:%S", tz = 'utc')

# saving this is likely unnecessary
# setwd(paste0(results_path,"/corrected_fieldsheets"))
# f_name <- paste0("Fieldsheet-GHG_corrected.csv")
# write.csv(file = f_name, x = fieldsheet, row.names = F)


# ---- Process incubations timeseries and compute fluxes ----


# ---- Go through each incubation in fieldsheet and compute linear model for co2 and ch4
subsites <- unique(fieldsheet$subsite)
isF_incub <- T
isFsubsite <- T
for (subsite in subsites){
  message("Now processing ",subsite)
  
  corresp_fs <- fieldsheet[fieldsheet$subsite == subsite,]
  gs <- first(corresp_fs$gas_analyzer)
  corresp_fs <- fieldsheet[fieldsheet$subsite == subsite,]
  corresp_fs$unix_start <- corresp_fs$unix_start_corrected
  corresp_fs$unix_stop <- corresp_fs$unix_stop_corrected
  corresp_fs <- corresp_fs[!is.na(corresp_fs$unix_start),]
  
  
  # read corresponding temperature logger file and keep initial temperature
  # --- Read corresponding Loggers data ----
  SN_logger_float <- first(corresp_fs$logger_floating_chamber)
  SN_logger_tube <- first(corresp_fs$logger_transparent_chamber)
  site_ID <- str_sub(subsite, start = 1, end = 5)
  
  # finding out if corresponding file exists, and its extension
  require(tools)
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
  
  
  path2data <- paste0(results_path,"/",subsite)
  if(dir.exists(path2data)){
    
    setwd(path2data)
    load(file = paste0("data_",subsite,".RData"))
    
    
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
        
        UniqueID = paste(subsite, seq_along(corresp_fs$pilot_site),corresp_fs$strata,corresp_fs$transparent_dark,sep = "-")[incub]
        
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
    
    # draw a few incubations randomly
    # draw <- sample(seq_along(auxfile$subsite), 3)
    
    # Define the measurements' window of observation
    auxfile <- auxfile
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
    
    # Plots results
    # Make a list of plots of all measurements, for each gastype
    CO2_flux_plots <- flux.plot(CO2_flux_res_auto, mydata_auto, "CO2dry_ppm")
    H2O_flux_plots <- flux.plot(H2O_flux_res_auto, mydata_auto, "H2O_ppm")
    CH4_flux_plots <- flux.plot(CH4_flux_res_auto, mydata_auto, "CH4dry_ppb")
    
    # estimate ch4 diffusion and ebullition components---------| METHOD 1 |-----------
    CH4_res_meth1 <- CH4_flux_res_auto
    CH4_res_meth1$total_estimated <- NA
    CH4_res_meth1$ebullition <- NA
    CH4_res_meth1$diffusion <- NA
    
    
    for (i in seq_along(auxfile)){
      if(auxfile$water_depth[i]>0){
        my_incub <- mydata[as.numeric(mydata$POSIX.time)> auxfile$start.time[i] &
                             as.numeric(mydata$POSIX.time)< auxfile$start.time[i]+auxfile$duration[i],]
        my_incub <- my_incub[!is.na(my_incub$CO2dry_ppm),]
        # calling dedicated function
        df_ebull <- separate_ebullition_from_diffusion(my_incub = my_incub, UniqueID = auxfile$UniqueID[i])
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
    
    
    setwd(path2data)
    write.csv(x = CO2_flux_res_auto, file = paste0(subsite,"_co2_fluxes.csv"), row.names = F)
    write.csv(x = CH4_res_meth1, file = paste0(subsite,"_ch4_fluxes.csv"), row.names = F)
    
    
    # Combine plot lists into one list
    flux_plot.ls <- c(CO2_flux_plots, CH4_flux_plots, H2O_flux_plots)
    
    # Save plots to pdf
    setwd(plots_path)
    myfilename <- paste(subsite, as.character(as.Date(first(auxfile$start.time))),sep="_")
    # myfilename <- gsub("-","", myfilename)
    flux2pdf(flux_plot.ls, outfile = paste0(plots_path,"/",myfilename,".pdf"))
    
    
    table_results <- auxfile %>%
      left_join(CO2_flux_res_auto %>% select(UniqueID, best.flux, model, quality.check)) %>%
      rename(CO2_flux = best.flux, CO2_model = model, CO2_quality.check = quality.check) %>%
      left_join(CH4_res_meth1 %>% select(UniqueID, best.flux, model, quality.check, diffusion, ebullition)) %>%
      rename(CH4_flux = best.flux, CH4_model = model, CH4_quality.check = quality.check, CH4_diffusive_flux = diffusion, CH4_ebullitive_flux = ebullition)
    
    
    if (isFsubsite){
      isFsubsite <- F
      table_results_all <- table_results
    } else {
      table_results_all <- rbind(table_results_all, table_results)
    }
  } # closing if(dir.exists(path2data)){
}



# ---- Some plots ----


plt_CO2 <- ggplot(table_results_all, aes(strata, CO2_flux, fill = lightCondition))+
  geom_hline(yintercept = 0)+
  geom_boxplot(alpha=0.2)+geom_jitter(width = 0.2)+
  theme_article()+
  ylab("CO2 flux mmol/m2/s")+
  ggtitle("CO2 flux")+
  scale_fill_viridis_d(begin = 0.2, end = 0.9)+facet_grid(.~subsite)


plt_CH4diff <- ggplot(table_results_all, aes(strata, CH4_flux, fill = lightCondition))+
  geom_hline(yintercept = 0)+
  geom_boxplot(alpha=0.2)+geom_jitter(width = 0.2)+
  theme_article()+
  ylab("CH4 flux nmol/m2/s")+
  ggtitle("CH4 flux")+
  scale_fill_viridis_d(begin = 0.2, end = 0.9)+facet_grid(.~subsite)


ggarrange(plt_CO2, plt_CH4diff, ncol = 1)



