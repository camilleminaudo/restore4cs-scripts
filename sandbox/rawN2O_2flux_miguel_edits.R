

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
        
        #Parameters for floating chambers:
        if (corresp_fs$chamber_type[incub] == "floating"){
          if(is_data_logger_float){
            #Calculate median temperature of data_logger_float during duration of floating incubation.
            myTemp <- median(approx(data_logger_float$unixtime, data_logger_float$temperature, xout = as.numeric(corresp_fs$unix_start[incub]:corresp_fs$unix_stop[incub]))$y )
          }
          
          #Floating chamber used in Restore4Cs is half sphere of 38cm diameter. radius=19cm
          myArea<- pi*19^2 #cm2
          myVtot<- (((4/3)*pi*(19^3))/1000 )/2# L

          #OLD VALUES, Unclear where they came from, wrong!          
          # myArea = 14365.4439 # cm2
          # myVtot = 115/2 # L
          
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

      if(is.na(myTemp)){myTemp <- 25} #ºC (a default temperature to run the flux calculation if no temperature is found)
      myPcham <- 100.1 #kPa (default atmospheric pressure is assumed for all chambers)
      
        
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

   
# ---- 2. N2O flux calculation-------

#Overall approach:
#1st calculate deltaN2O in nmols (using chamber volume and Ideal Gas Law)
#2nd calculate change in nmols per second 
#3rd calculate areal flux in nmols s^-1 m^-2  

#Using Ideal Gas Law: PV = nRT, where:
#P: Pressure (atm)
#V: Volume (L)
#n: number of mols
#T: Temperature (ºK)
#R: Gas Constant,  R = 0.08206 (L·atm·K-1·mol-1)
R_constant<- 0.08206

#Re-arranging: n = PV/RT

#to calculate mols of gas in the chamber at t0 and tf, we will 1st calculate the Volume of n2o (Vn2o, in L) using the ppm and total volume of chamber. Then apply the ideal gas law to get the mols at the given Tcham and Pressure.

#the goflux package uses a slightly different function, with P in KPa, and incorporating a water-vapor correction. Unclear why, we will stick to the more simple Ideal Gas Law. 

# 4.1. Import N2O_ppm and transform units.

#Import atm and tf N2O concentration in ppm, identified with "UniqueID_notime", atm_N2Oppm, tf_N2Oppm 
n2o<- read.csv(paste0(results_path,"/S4_restore4cs_N2Oppm_atmchambers.csv"))


#FLUX CAlCULATION: #Auxfile has Pcham in kPa, Area in cm2 and Tcham in ºcelsius, Vtot is already in L
n2o_flux<- n2o %>% 
  #Join with auxfile
  merge.data.frame(auxfile, by="UniqueID_notime", all.x = T) %>% 
  #Transform Units:
  mutate(Vtot, #Vtot is already in L
         Pcham_atm= 1, # set Pcham to 1 atm
         Tcham_k= Tcham+273.15, #Transform Tcham to Kelvin
         duration,#Duration of incubation (in seconds) calculated from fieldsheet times
         Area_m2=Area*1e-4) %>% # divide by 10,000 to transform cm2 to m2  
  #calculate Vn2o, at t0 and tf
  mutate(V_n2o_t0 = (atm_avgN2Oppm/1e6)*Vtot, 
         V_n2o_tf = (tf_avgN2Oppm/1e6)*Vtot) %>% 
  #calculate nmol_n2o at t0 and tf (apply mol = PV/RT and then multiply by 1e9 to get nmol)
  mutate(nmol_n2o_t0 = 1e9* ( (Pcham_atm*V_n2o_t0)/(R_constant*Tcham_k) ),
         nmol_n2o_tf = 1e9* ( (Pcham_atm*V_n2o_tf)/(R_constant*Tcham_k) )) %>% 
  #calculate delta_nmol_n2o (as nmol difference tf-t0)
  mutate(delta_nmol_n2o = nmol_n2o_tf-nmol_n2o_t0) %>% 
  #calculate absolute flux: nmol_n2o_per_second
  mutate(nmol_n2o_per_second=delta_nmol_n2o/duration) %>% 
  #calculate areal flux: nmol_n2o_per_second_m2
  mutate(nmol_n2o_per_second_m2=nmol_n2o_per_second/Area_m2) %>% 


#CHUNK to calculate the significance of the fluxes (i.e. if the concentration at tf is significantly different from atmospheric.)
mutate(
  # Calculate standard deviations from cv (cv = SD / Mean, so SD = Mean * cv)
  tf_sd = tf_avgN2Oppm * tf_cv,
  atm_sd = atm_avgN2Oppm * atm_cv,
  
  # Calculate standard error for the difference
  se_diff = sqrt((tf_sd^2 / tf_n) + (atm_sd^2 / atm_n)),
  
  # Calculate the t-statistic for the difference
  t_stat = (tf_avgN2Oppm - atm_avgN2Oppm) / se_diff,
  
  # Calculate the degrees of freedom using the formula for unequal variances
  df_ = (((tf_sd^2 / tf_n) + (atm_sd^2 / atm_n))^2) /
    (((tf_sd^2 / tf_n)^2 / (tf_n - 1)) + ((atm_sd^2 / atm_n)^2 / (atm_n - 1))),
  
  # Calculate the p-value using the t-distribution
  p_value = 2 * pt(-abs(t_stat), df_)
)





# ---- 3. N2O flux filtering-------


#Add significance and clean up table (remove unncecessary variables.)
#FILTER also fluxes of which we cannot be sure (plotincubation==2, unclear or incomplete ventilation between transparent and dark incubation)
n2o_flux_simple<- n2o_flux %>% 
  separate(UniqueID_notime, into = c("sampling","pilotsite","subsite","plot","d1","d2"),remove = F) %>% 
  mutate(siteID=paste(pilotsite,subsite, sep = "-")) %>% 
  filter(plotincubation==1) %>% 
  select(UniqueID_notime, sampling, pilotsite, subsite, siteID, start.time, duration, water_depth, strata, chamberType,lightCondition, 
         nmol_n2o_per_second_m2, p_value)
  

message(paste("All Fluxes: Out of ", dim(n2o_flux)[1]), " N2O fluxes, we have ", sum(n2o_flux$p_value<0.05), " significant fluxes (",  round(sum(n2o_flux$p_value<0.05)/dim(n2o_flux)[1]*100,2), "%), at the p<0.05 level")


message(paste("Only reliable fluxes (clear ventilation): Out of ", dim(n2o_flux_simple)[1]), " N2O fluxes, we have ", sum(n2o_flux$p_value<0.05), " significant fluxes (",  round(sum(n2o_flux$p_value<0.05)/dim(n2o_flux_simple)[1]*100,2), "%), at the p<0.05 level")


#__________________####



#SAVE RESULTS####
write.csv(n2o_flux, file = paste0(datapath,"/S4_restore4cs_N2O_EXTRADATA_n2ofluxes.csv.csv"),row.names = F)

write.csv(n2o_flux_simple, file=paste0(datapath,"/S4_restore4cs_N2O_arealflux_nmol_s-1_m-2.csv.csv"),row.names = F)



#Miscelanea: 

ggplot(n2o_flux_simple, aes(y=nmol_n2o_per_second_m2,  fill=p_value<0.05))+
  geom_histogram(position = "identity", alpha=0.3)

