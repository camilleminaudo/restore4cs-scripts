# ---
# Authors: Camille Minaudo
# Project: "RESTORE4Cs"
# date: "Oct 2023"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
# This script lists all incubations corresponding to a suite of criterias and create a single RData file out of it
# This data will later be used for automatic versus expert flux assessment.

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
# SPECIFY HERE YOUR NAME
username <- "Camille"

nb_draw <- 5

############################



# ---- Directories ----
dropbox_root <- "C:/Users/Camille Minaudo/Dropbox/RESTORE4Cs - Fieldwork/Data" # You have to make sure this is pointing to the write folder on your local machine
datapath <- paste0(dropbox_root,"/GHG/RAW data")
fieldsheetpath <- paste0(dropbox_root,"/GHG/Fieldsheets")
loggerspath <- paste0(datapath,"/RAW Data Logger")
plots_path <- paste0(dropbox_root,"/GHG/Processed data/plots_all_incubations/")
RData_path <- paste0(dropbox_root,"/GHG/Processed data/RData/")
results_path <- paste0(dropbox_root,"/GHG/Processed data/computed_flux/")


# ----- Data pre-processing and harmonization -----

# ---- List GHG chamber fieldsheets in Dropbox and read them ---
# list filenames
myfieldsheets_list <- list.files(fieldsheetpath, pattern = "Fieldsheet-GHG.xlsx", all.files = T, full.names = T, recursive = T)

#S3 data is not ready yet
i <- grep(pattern = "S3", x = myfieldsheets_list) # selecting the files corresponding to the selected sampling campaign
myfieldsheets_list <- myfieldsheets_list[-i]
# Read all fieldsheets and put them in a single dataframe
fieldsheet <- read_GHG_fieldsheets(myfieldsheets_list)



# ---- Correct fieldsheets in the case of Picarro data ---

read_map_incubations <- function(path2folder){
  
  my_maps_filenames <- list.files(path2folder, pattern = ".csv", all.files = T, full.names = T, recursive = T)
  i <- grep(pattern = "Picarro", x = my_maps_filenames) # selecting the files corresponding to the Picarro only
  my_maps_filenames <- my_maps_filenames[i]
  # i <- grep(pattern = sampling, x = my_maps_filenames) # selecting the files corresponding to the selected sampling campaign
  # my_maps_filenames <- my_maps_filenames[i]
  
  
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
  interp <- approx(x = seq_along(fieldsheet_Picarro$plot_id), y = corresponding_row, xout = ind_NAs)$y
  is_integer <- (interp - floor(interp)) == 0
  corresponding_row[ind_NAs][is_integer] <- interp[is_integer]
}
if(sum(is.na(corresponding_row))>0){
  ind_NAs <- which(is.na(corresponding_row))
  message("Could not find a corresponding incubation map for the following incubations:")
  fieldsheet_Picarro$UniqueID[ind_NAs]
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



# ----- Data selection and random draw of incubations -----

# Here, we only focus on open water measurements with the floating chamber
fieldsheet <- fieldsheet[which(fieldsheet$chamber_type=="floating" & fieldsheet$strata=="open water" & fieldsheet$water_depth>0),]

# draw a few incubations randomly
draw <- sample(seq_along(fieldsheet$subsite), nb_draw)

table_draw <- data.frame(draw = draw,
                         subsite = fieldsheet$subsite[draw],
                         uniquID = fieldsheet$uniqID[draw])


# ---- Go through each incubation selected and build auxfile table ----
auxfile <- NULL
isF_incub <- T
isFsubsite <- T
for (i in seq_along(table_draw$draw)){
  message("Now processing draw ",i, " out of ",nb_draw,", corresponding to incubation ",table_draw$uniquID[i])
  
  corresp_fs <- fieldsheet[fieldsheet$uniqID == table_draw$uniquID[i],]
  subsite <- corresp_fs$subsite
  gs <- corresp_fs$gas_analyzer
  
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
  
  # load corresponding measurements
  load(file = paste0(subsite,"_",gs_suffix,".RData"))
  
  # isolate incubation
  my_incub <- mydata[as.numeric(mydata$POSIX.time)> corresp_fs$unix_start &
                       as.numeric(mydata$POSIX.time)< corresp_fs$unix_stop,]
  my_incub <- my_incub[!is.na(my_incub$CO2dry_ppm),]
  
  if (dim(my_incub)[1]>0){
    
    # Compute Temperature, Area and Volume from fieldsheet info
    myTemp <- 15 # a default temperature to run the flux calculation...
    if (corresp_fs$chamber_type == "floating"){
      if(is_data_logger_float){
        myTemp <- median(approx(data_logger_float$unixtime, data_logger_float$temperature, xout = as.numeric(my_incub$POSIX.time))$y )
      }
      myArea = 14365.4439 # cm2
      myVtot = 115/2 # L
    } else if (corresp_fs$chamber_type == "tube"){
      if(is_data_logger_tube){
        myTemp <- median(approx(data_logger_tube$unixtime, data_logger_tube$temperature, xout = as.numeric(my_incub$POSIX.time))$y )
      }
      myArea = pi*12.1**2 # cm2
      myVtot = myArea*as.numeric(corresp_fs$chamber_height_cm)*1e-3 # L
    } else {
      warning("chamber type not correct!")
    }
    myPcham <- 100.1 #kPa
    
    # --- Create auxfile table ----
    # An auxfile table, made of fieldsheet, adding important variables. The auxfile
    # requires start.time and UniqueID.
    # start.time must be in the format "%Y-%m-%d %H:%M:%S"
    
    auxfile_tmp <- data.frame(subsite = subsite,
                              UniqueID = corresp_fs$uniqID,
                              gas_analiser = gs,
                              start.time = as.POSIXct((corresp_fs$unix_start), tz = "UTC"),
                              duration = (corresp_fs$unix_stop) - (corresp_fs$unix_start),
                              water_depth = corresp_fs$water_depth,
                              Area = myArea,
                              Vtot = myVtot,
                              Tcham = myTemp,
                              Pcham = myPcham,
                              strata = corresp_fs$strata,
                              chamberType = corresp_fs$chamber_type,
                              lightCondition = corresp_fs$transparent_dark)
    
    if(is.null(auxfile)){
      auxfile <- auxfile_tmp
      mydata_all <- mydata
    } else {
      auxfile <- rbind(auxfile, auxfile_tmp)
      mydata_all <- rbind(mydata_all, mydata)
    }
    rm(mydata)
    
  }
}




# ---- Process incubations timeseries and compute fluxes ----

# Define the measurements' window of observation
myauxfile <- auxfile
mydata_ow <- obs.win(inputfile = mydata_all, auxfile = myauxfile,
                     obs.length = myauxfile$duration, shoulder = 2)


# ----------- Compute fluxes after manual selection of CO2 data

# Manually identify measurements by clicking on the start and end points
mydata_manID <- lapply(seq_along(mydata_ow), click.peak.loop,
                       flux.unique = mydata_ow,
                       gastype = "CO2dry_ppm",
                       plot.lim = c(200,1000)) %>%
  map_df(., ~as.data.frame(.x))


# Additional auxiliary data required for flux calculation.
mydata_manID <- mydata_manID %>%
  left_join(auxfile %>% select(username, UniqueID, Area, Vtot, Tcham, Pcham))

# Add instrument precision for each gas
prec = c(3.5, 0.6, 0.4, 45, 45)
mydata_manID <- mydata_manID %>%
  mutate(CO2_prec = prec[1], CH4_prec = prec[2], N2O_prec = prec[3],
         H2O_prec = prec[4])


# Calculate fluxes
CO2_results_manID <- goFlux(mydata_manID, "CO2dry_ppm")
H2O_results_manID <- goFlux(mydata_manID, "H2O_ppm")
CH4_results_manID <- goFlux(mydata_manID, "CH4dry_ppb")

# Use best.flux to select the best flux estimates (LM or HM)
# based on a list of criteria
criteria <- c("g.factor", "kappa", "MDF", "R2", "SE.rel")

CO2_flux_res_manID <- best.flux(CO2_results_manID, criteria)
H2O_flux_res_manID <- best.flux(H2O_results_manID, criteria)
CH4_flux_res_manID <- best.flux(CH4_results_manID, criteria)

# CO2_flux_res_manID <- CO2_flux_res_manID %>%
#   left_join(myauxfile %>% select(UniqueID, strata, chamberType, lightCondition))


# Plots results
# Make a list of plots of all measurements, for each gastype
# CO2_flux_plots <- flux.plot(CO2_flux_res_manID, mydata_manID, "CO2dry_ppm")
# H2O_flux_plots <- flux.plot(H2O_flux_res_manID, mydata_manID, "H2O_ppm")
# CH4_flux_plots <- flux.plot(CH4_flux_res_manID, mydata_manID, "CH4dry_ppb")

# ----------- Compute fluxes blindly without any manual selection

mydata_ow <- obs.win(inputfile = mydata_all, auxfile = myauxfile,
                     obs.length = myauxfile$duration, shoulder = 2)

# Join mydata_ow with info on start end incubation
mydata_auto <- lapply(seq_along(mydata_ow), join_auxfile_with_data.loop, flux.unique = mydata_ow) %>%
  map_df(., ~as.data.frame(.x))

# Additional auxiliary data required for flux calculation.
mydata_auto <- mydata_auto %>%
  left_join(myauxfile %>% select(username, UniqueID, Area, Vtot, Tcham, Pcham))

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


# Use best.flux to select the best flux estimates (LM or HM)
# based on a list of criteria
criteria <- c("g.factor", "kappa", "MDF", "R2", "SE.rel")

CO2_flux_res_auto <- best.flux(CO2_results_auto, criteria)
H2O_flux_res_auto <- best.flux(H2O_results_auto, criteria)
CH4_flux_res_auto <- best.flux(CH4_results_auto, criteria)


CO2_flux_res_auto$flux_method <- "Blind"
CO2_flux_res_manID$flux_method <- "Expert"

CH4_flux_res_auto$flux_method <- "Blind"
CH4_flux_res_manID$flux_method <- "Expert"




# ----------- Estimate CH4 diffusion/ebullition contributions
# method 1: density of prob. of first derivative
# method 2: manual selection of linear pattern in CH4

# estimate ch4 diffusion and ebullition components---------| METHOD 1 |-----------
CH4_res_meth1 <- CH4_flux_res_auto
CH4_res_meth1$total_estimated <- NA
CH4_res_meth1$ebullition <- NA
CH4_res_meth1$diffusion <- NA


for (i in seq_along(auxfile$subsite)){
  if(auxfile$water_depth[i]>0){
    my_incub <- mydata_all[as.numeric(mydata_all$POSIX.time)> auxfile$start.time[i] &
                             as.numeric(mydata_all$POSIX.time)< auxfile$start.time[i]+auxfile$duration[i],]
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
    CH4_flux_total <- CH4_flux_diff <- CH4_res_meth1$best.flux[which(CH4_res_meth1$UniqueID==auxfile$UniqueID[i])]
    CH4_flux_ebull <- 0
  }
  CH4_res_meth1$total_estimated[which(CH4_res_meth1$UniqueID==auxfile$UniqueID[i])] <- CH4_flux_total
  CH4_res_meth1$ebullition[which(CH4_res_meth1$UniqueID==auxfile$UniqueID[i])] <- CH4_flux_ebull
  CH4_res_meth1$diffusion[which(CH4_res_meth1$UniqueID==auxfile$UniqueID[i])] <- CH4_flux_diff
}






# Estimate ch4 diffusion and ebullition components ---------| METHOD 2 |-----------


# Manually identify diffusive (more or less linear) CH4 behaviors by clicking on the start and end points
myCH4_diffusion <- lapply(seq_along(mydata_ow), click.peak.loop,
                          flux.unique = mydata_ow,
                          gastype = "CH4dry_ppb",
                          plot.lim = c(1900,max(fieldsheet$final_ch4)*1000)) %>%
  map_df(., ~as.data.frame(.x))


myCH4_diffusion <- myCH4_diffusion %>%
  mutate(CO2_prec = prec[1], CH4_prec = prec[2], N2O_prec = prec[3],
         H2O_prec = prec[4])


# Calculate fluxes for CH4
CH4_results_diffusion <- goFlux(myCH4_diffusion, "CH4dry_ppb")

# Use best.flux to select the best flux estimates (LM or HM)
# based on a list of criteria
criteria <- c("g.factor", "kappa", "MDF", "R2", "SE.rel")

CH4_res_diff <- best.flux(CH4_results_diffusion, criteria)


CH4_res_meth2 <- CH4_flux_res_manID
CH4_res_meth2$total_estimated <- NA
CH4_res_meth2$ebullition <- NA
CH4_res_meth2$diffusion <- CH4_res_diff$best.flux


# Estimating ebullition component
for (id in unique(CH4_res_meth2$UniqueID)){
  i <- which(CH4_res_meth2$UniqueID == id)
  
  CH4_final <- CH4_flux_res_manID$Ct[i]
  CH4_initial <-  CH4_flux_res_manID$C0[i]
  incubation_time <- myauxfile$duration[which(myauxfile$UniqueID == id)]
  CH4_res_meth2$total_estimated[i] <- (CH4_final-CH4_initial)/incubation_time*CH4_flux_res_manID$flux.term[i] # nmol/m2/s
  CH4_res_meth2$ebullition[i] <- CH4_res_meth2$total_estimated[i]  - CH4_res_meth2$diffusion[i] # total flux - diffusive term
}


# compare method 1 and 2 to estimate ebullitive contribution
CH4_res_meth1$flux_method <- "dydt"
CH4_res_meth2$flux_method <- "Expert"


CO2_flux_res_auto$variable <- "CO2"
CO2_flux_res_manID$variable <- "CO2"
CH4_flux_res_auto$variable <- "CH4"
CH4_flux_res_manID$variable <- "CH4"


table_results <- rbind(CO2_flux_res_auto, CO2_flux_res_manID, CH4_flux_res_auto, CH4_flux_res_manID)
table_results$username <- username

table_results_ebull <- rbind(CH4_res_meth1, CH4_res_meth2)
table_results_ebull$username <- username

# saving fluxes estimates
setwd(results_path)

append_if_exists <- function(filename, data){
  if(file.exists(filename)){
    data <- rbind(read.csv(filename), data)
  } else {
    write.csv(x = data, file = filename, 
              row.names = F)
  }
  return(data)
}


filename <- paste0("BLIND_vs_EXPERT_co2_ch4_fluxes_",Sys.Date(),".csv")
table_results <- append_if_exists(filename, data = table_results)


filename <- paste0("BLIND_vs_EXPERT_ch4_ebullition_",Sys.Date(),".csv")
table_results_ebull <- append_if_exists(filename, data = table_results_ebull)



#----- some plots for the incubations processed in the current session -----



p_auto_vs_manual <- ggplot(data = table_results)+
  geom_abline(slope = 0,intercept = 0, color = 'black')+
  # geom_segment(data = data.frame(UniqueID = CO2_flux_res_auto$UniqueID,
  #                                meth1 = CO2_flux_res_auto$best.flux,
  #                                meth2 = CO2_flux_res_manID$best.flux), aes(x=UniqueID, xend=UniqueID, y = meth1, yend = meth2), linewidth=1, alpha = 0.5)+
  # aes(x=reorder(X_Variable, -Y_Variable, FUN=mean), y=Y_Variable)
  geom_point(aes(reorder(UniqueID, -best.flux, FUN=mean), best.flux, colour = flux_method), size=4, alpha = 0.5)+
  ylab("flux [(mmolCO2 or nmolCH4)/m2/s]")+
  theme_bw()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_colour_viridis_d(begin = 0.1, end = 0.9, option = "F")+facet_grid(.~variable, scales = 'free')+coord_flip()

p_auto_vs_manual
ggsave(filename = paste0("BLIND_vs_EXPERT_co2_ch4_fluxes_",Sys.Date(),".jpeg"), 
       plot = p_auto_vs_manual, path = results_path, width = 10, height = 10, units = 'in', dpi = 300)


p_auto_vs_manual_ch4_ebullition <- ggplot(data = table_results_ebull)+
  geom_abline(slope = 0,intercept = 0, color = 'lightgrey')+
  geom_segment(data = data.frame(UniqueID = CH4_res_meth1$UniqueID,
                                 meth1 = CH4_res_meth1$ebullition,
                                 meth2 = CH4_res_meth2$ebullition), aes(x=UniqueID, xend=UniqueID, y = meth1, yend = meth2), linewidth=1, alpha = 0.5)+
  geom_point(aes(reorder(UniqueID, -ebullition, FUN=mean), ebullition, colour = flux_method), size=4, alpha = 0.5)+
  
  ylab("ebullition component [nmol/m2/s]")+
  theme_bw()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_colour_viridis_d(begin = 0.1, end = 0.9, option = "F")+coord_flip()

p_auto_vs_manual_ch4_ebullition
ggsave(filename = paste0("BLIND_vs_EXPERT_ch4_ebullition_",Sys.Date(),".jpeg"), 
       plot = p_auto_vs_manual_ch4_ebullition, path = results_path, width = 10, height = 10, units = 'in', dpi = 300)




table_results_sprd_CO2  <- table_results[table_results$variable=="CO2",c("UniqueID","flux_method","best.flux")] %>%
  pivot_wider(names_from = flux_method, values_from = c(best.flux))

median((table_results_sprd_CO2$Blind-table_results_sprd_CO2$Expert)/table_results_sprd_CO2$Expert*100, na.rm = T)

ggplot(data = table_results_sprd_CO2)+
  geom_abline(slope = 0,intercept = 0, color = 'lightgrey')+
  # geom_segment(data = data.frame(UniqueID = CO2_flux_res_auto$UniqueID,
  #                                meth1 = CO2_flux_res_auto$best.flux,
  #                                meth2 = CO2_flux_res_manID$best.flux), aes(x=UniqueID, xend=UniqueID, y = meth1, yend = meth2), linewidth=1, alpha = 0.5)+
  geom_point(aes(UniqueID, (Blind-Expert)/Expert*100), size=4, alpha = 0.5)+
  ylab("CO2 flux relative difference [%]")+
  theme_bw()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(table_results_sprd_CO2, aes((Blind-Expert)/Expert*100))+geom_density()


table_results_sprd_CH4  <- table_results[table_results$variable=="CH4",c("UniqueID","flux_method","best.flux")] %>%
  pivot_wider(names_from = flux_method, values_from = c(best.flux))

median((table_results_sprd_CH4$Blind-table_results_sprd_CH4$Expert)/table_results_sprd_CH4$Expert*100, na.rm = T)

ggplot(data = table_results_sprd_CH4)+
  geom_abline(slope = 0,intercept = 0, color = 'lightgrey')+
  # geom_segment(data = data.frame(UniqueID = CH4_flux_res_auto$UniqueID,
  #                                meth1 = CH4_flux_res_auto$best.flux,
  #                                meth2 = CH4_flux_res_manID$best.flux), aes(x=UniqueID, xend=UniqueID, y = meth1, yend = meth2), linewidth=1, alpha = 0.5)+
  geom_point(aes(UniqueID, (Blind-Expert)/Expert*100), size=4, alpha = 0.5)+
  ylab("CH4 flux relative difference [%]")+
  theme_bw()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))




#----- some plots for the all incubations processed so far -----


list_f <- list.files(path = results_path, pattern = "BLIND_vs_EXPERT_co2_ch4_fluxes_", full.names = T)
list_f <- list_f[grep(pattern = "csv", x = list_f)]

isF <- T
for (f in list_f){
  table_results.tmp <- read.csv(f)
  if(isF){
    isF <- F
    table_results_all <- table_results.tmp
  } else {
    table_results_all <- rbind(table_results_all,table_results.tmp)
  }
}

p_auto_vs_manual <- ggplot(data = table_results_all)+
  geom_abline(slope = 0,intercept = 0, color = 'black')+
  # geom_segment(data = data.frame(UniqueID = CO2_flux_res_auto$UniqueID,
  #                                meth1 = CO2_flux_res_auto$best.flux,
  #                                meth2 = CO2_flux_res_manID$best.flux), aes(x=UniqueID, xend=UniqueID, y = meth1, yend = meth2), linewidth=1, alpha = 0.5)+
  # aes(x=reorder(X_Variable, -Y_Variable, FUN=mean), y=Y_Variable)
  geom_point(aes(reorder(UniqueID, -best.flux, FUN=mean), best.flux, colour = flux_method), size=4, alpha = 0.5)+
  ylab("flux [(mmolCO2 or nmolCH4)/m2/s]")+
  theme_bw()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_colour_viridis_d(begin = 0.1, end = 0.9, option = "F")+facet_grid(.~variable, scales = 'free')+coord_flip()

p_auto_vs_manual

ggsave(filename = "BLIND_vs_EXPERT_co2_ch4_all.jpeg", plot = p_auto_vs_manual, path = results_path, width = 10, height = 10, units = 'in', dpi = 300)




list_f <- list.files(path = results_path, pattern = "BLIND_vs_EXPERT_ch4_ebullition_", full.names = T)
list_f <- list_f[grep(pattern = "csv", x = list_f)]

isF <- T
for (f in list_f){
  table_results.tmp <- read.csv(f)
  if(isF){
    isF <- F
    table_results_all <- table_results.tmp
  } else {
    table_results_all <- rbind(table_results_all,table_results.tmp)
  }
}




p_auto_vs_manual_ch4_ebullition <- ggplot(data = table_results_ebull)+
  geom_abline(slope = 0,intercept = 0, color = 'lightgrey')+
  geom_segment(data = data.frame(UniqueID = CH4_res_meth1$UniqueID,
                                 meth1 = CH4_res_meth1$ebullition,
                                 meth2 = CH4_res_meth2$ebullition), aes(x=UniqueID, xend=UniqueID, y = meth1, yend = meth2), linewidth=1, alpha = 0.5)+
  geom_point(aes(reorder(UniqueID, -ebullition, FUN=mean), ebullition, colour = flux_method), size=4, alpha = 0.5)+
  
  ylab("ebullition component [nmol/m2/s]")+
  theme_bw()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_colour_viridis_d(begin = 0.1, end = 0.9, option = "F")+coord_flip()

p_auto_vs_manual_ch4_ebullition
ggsave(filename = "BLIND_vs_EXPERT_ch4_ebullition_all.jpeg", plot = p_auto_vs_manual_ch4_ebullition, path = results_path, width = 10, height = 6, units = 'in', dpi = 300)










