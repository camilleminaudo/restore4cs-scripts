

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

#Auxtable extraction adapted from raw2flux_miguel_edits_perGHGcrop (in order to use meteo data)

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
meteo_path<- paste0(dropbox_root, "/Meteo/Formated_data/")
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

#Load data for all meteostations and correspondence of subsite-meteostation (used to subset the correct data)
meteocorresp <- read_csv(paste0(meteo_path, "correspondence_subsite-meteostation.csv"),show_col_types = F)

meteodata <- read_csv(paste0(meteo_path,"allmeteostations_onlysamplingdays.csv"),show_col_types = F)



# This section loops over each subsite and extracts the necessary chamber parameters for flux calculation: T, P, V, A and saves all in a dataframe called auxfile, identified with UniqueID_notime
subsites <- unique(fieldsheet$subsite)

isF_incub <- T
isFsubsite <- T
auxfile <- NULL

for (subsite in subsites){
  message("Now processing ",subsite)
  
  corresp_fs <- fieldsheet[fieldsheet$subsite == subsite,]
  
  #Load the appropriate meteo-data for this subsite:
  meteostation<- meteocorresp[meteocorresp$site_subsite==substr(subsite,4,8),]$station_id
  subsite_meteo<- meteodata %>% filter(station_id==meteostation)
  
  #Loop over each incubation: 
  for (incub in seq_along(corresp_fs$plot_id)){
    
    #Create my_incub_unix (vector of unixtime for duration of incubation) 
    my_incub_unix<- seq(corresp_fs$unix_start[incub],
                        corresp_fs$unix_stop[incub],
                        by=1)
    
    #Default temperature and pressure in case they are not available in meteo data
      myTemp <- 15 #ºC 
      myPcham <- 1 #Atm
    
      # Compute Temperature (ºC) and pressure (Atm) for the initial moment of incub from the meteo data (using approx to interpolate from hourly data):
      
      myTemp <- first(approx(as.numeric(subsite_meteo$datetime_utc), subsite_meteo$temp_c, xout=my_incub_unix)$y)
      
      myPcham<- first(approx(as.numeric(subsite_meteo$datetime_utc), subsite_meteo$Patm_hPa/1013.25, xout=my_incub_unix)$y)
      
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
      
      
      auxfile_tmp <- data.frame(subsite = subsite,
                                # UniqueID = corresp_fs$uniqID[incub],
                                UniqueID_notime= corresp_fs$UniqueID_notime[incub],
                                start.time = as.POSIXct((corresp_fs$unix_start[incub]), tz = "UTC"),
                                duration = (corresp_fs$unix_stop[incub]) - (corresp_fs$unix_start[incub]),
                                water_depth = corresp_fs$water_depth[incub],#cm
                                Area = myArea,#cm2
                                Vtot = myVtot,#L
                                Tcham = myTemp,#ºC
                                Pcham = myPcham,#atm
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

rm(corresp_fs, auxfile_tmp, isF_incub, isFsubsite, myArea, myPcham, myTemp, myVtot, subsite, subsites, myfieldsheets_list)

   
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



#Import N2O_ppm and transform units.

#Import atm and tf N2O concentration in ppm, identified with "UniqueID_notime", atm_N2Oppm, tf_N2Oppm 
n2o<- read.csv(paste0(results_path,"/S4_restore4cs_N2Oppm_atmchambers.csv"))


#FLUX CAlCULATION: #Auxfile has Pcham in kPa, Area in cm2 and Tcham in ºcelsius, Vtot is already in L
n2o_flux<- n2o %>% 
  #Join with auxfile
  merge.data.frame(auxfile, by="UniqueID_notime", all.x = T) %>% 
  #Transform Units:
  mutate(Vtot, #Vtot is already in L
         Pcham_atm=Pcham, # Pcham already in atm units
         Tcham_k= Tcham+273.15, #Transform Tcham to Kelvin
         duration,#Duration of incubation (in seconds) calculated from fieldsheet times
         Area_m2=Area*1e-4) %>% # divide by 10,000 to transform cm2 to m2  
  #calculate Vn2o, at t0 and tf in litres
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
  # Calculate standard deviations from cv (cv = SD / Mean, so SD = Mean * cv) in nmol
  tf_sd_nmol = nmol_n2o_tf * tf_cv,
  atm_sd_nmol = nmol_n2o_t0 * atm_cv,
  
  # Calculate standard error for the difference
  se_diff_nmol = sqrt((tf_sd_nmol^2 / tf_n) + (atm_sd_nmol^2 / atm_n)),
  
  # Calculate the t-statistic for the difference
  t_stat = delta_nmol_n2o/se_diff_nmol,
  
  # Calculate the degrees of freedom using the formula for unequal variances
  df_ = (((tf_sd_nmol^2 / tf_n) + (atm_sd_nmol^2 / atm_n))^2) /
    (((tf_sd_nmol^2 / tf_n)^2 / (tf_n - 1)) + ((atm_sd_nmol^2 / atm_n)^2 / (atm_n - 1))),
  
  # Calculate the p-value using the t-distribution
  p_value = 2 * pt(-abs(t_stat), df_)
) %>% 
  #get se_diff_nmol as nmol_n2o_per_second_m2 (same units)
  mutate(SE_nmol_n2o_per_second_m2=(se_diff_nmol/duration)/Area_m2)




# ---- 3. N2O flux filtering-------


#Add significance and clean up table (remove unncecessary variables), keep SE of flux and pvalue
#FILTER also fluxes of which we did not ventilate the chamber (unclear or incomplete ventilation between transparent and dark incubation logged after inspection of CO2 and CH4 timeseries)

#Load ventilation check: 
vent_check<- read_xlsx(paste0(results_path, "/Dark_incubation_ventilationcheck.xlsx"),col_types = c("text","numeric","logical"))

not_ventilated<- vent_check %>% 
  filter(ventilated==F) %>% 
  pull(UniqueID_notime)



n2o_flux_simple<- n2o_flux %>% 
  separate(UniqueID_notime, into = c("sampling","pilotsite","subsite","plot","d1","d2"),remove = F) %>% 
  mutate(siteID=paste(pilotsite,subsite, sep = "-")) %>% 
  filter(!UniqueID_notime%in%not_ventilated) %>% 
  #rename variables: 
  rename(N2Oflux_nmol_per_second_per_m2=nmol_n2o_per_second_m2,
         N2Oflux_se=SE_nmol_n2o_per_second_m2,
         N2Oflux_pvalue=p_value) %>% 
  select(UniqueID_notime, sampling, pilotsite, subsite, siteID, start.time, duration, water_depth, strata, chamberType,lightCondition, 
         N2Oflux_nmol_per_second_per_m2, N2Oflux_se, N2Oflux_pvalue)
  



#Check dataset significance
message(paste("All Fluxes: Out of ", dim(n2o_flux)[1]), " N2O fluxes, we have ", sum(n2o_flux$p_value<0.05), " significant fluxes (",  round(sum(n2o_flux$p_value<0.05)/dim(n2o_flux)[1]*100,2), "%), at the p<0.05 level")


message(paste("Only reliable fluxes (clear ventilation): Out of ", dim(n2o_flux_simple)[1]), " reliable N2O incubations, we have ", sum(n2o_flux_simple$N2Oflux_pvalue<0.05), " significant fluxes (",  round(sum(n2o_flux_simple$N2Oflux_pvalue<0.05)/dim(n2o_flux_simple)[1]*100,2), "%), at the p<0.05 level")


#__________________####



#SAVE RESULTS####
write.csv(n2o_flux, file = paste0(datapath,"/S4_restore4cs_N2O_EXTRADATA_n2ofluxes.csv"),row.names = F)

write.csv(n2o_flux_simple, file=paste0(datapath,"/S4_restore4cs_N2O_arealflux_nmol_s-1_m-2.csv"),row.names = F)



#Miscelanea: 

ggplot(n2o_flux_simple, aes(y=N2Oflux_nmol_per_second_per_m2,  fill=N2Oflux_pvalue<0.05))+
  geom_histogram(position = "identity", alpha=0.3)

