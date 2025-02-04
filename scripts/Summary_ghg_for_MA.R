
# ---
# Authors: Miguel Cabrera
# Project: "RESTORE4Cs"
# date: "Feb 2025"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

rm(list = ls()) # clear workspace
cat("/014") # clear console

# --- Description of this script
# This script makes a summary of in-situ GHG CO2 and CH4 fluxes to be used in Meta-Analysis. Main decisions for fluxes selection and temporal integration are recorded throughout the script. 

#Input: per-season and per-subsite CSV files with computed fluxes directly produced by raw2flux.R script
#output: single csv with per-subsite global averages for CO2 and CH4 (mean+sd+N) (with stamp of date created)
# Final units: netCO2: umol/m2/s , netCH4: nmol/m2/s

#Filtering criteria for fluxes fit: 
#We will discard fluxes derived from models with MAE (Mean Absolute Error, arithmetic mean of the absolute residuals) larger than: 

MAE_threshold_CO2<- 10
MAE_threshold_CH4<- 100

#Threshold for MAE taken from Camille's decision in analyseFlux.R script. 
#Possibility to adapt filter based on R2 of models from go-flux


#Save today in order to paste daystamp to output csv-file:
dayoflastrun <- today()



# ---- packages ----
library(tidyverse)
library(readxl)
library(lubridate)
library(zoo)
library(ggplot2)
library(grid)
library(egg)
library(tools)
library(suncalc)


repo_root <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
files.sources = list.files(path = paste0(repo_root,"/functions"), full.names = T)
for (f in files.sources){source(f)}


# ---- Directories ----
dropbox_root <- "C:/Users/Miguel/Dropbox/RESTORE4Cs - Fieldwork/Data" # You have to make sure this is pointing to the write folder on your local machine
results_path <- paste0(dropbox_root,"/GHG/Processed data/computed_flux/")
export_path<- paste0(dropbox_root,"/GHG/Processed data/computed_flux/Summaries for MA/")
fieldsheetpath <- paste0(dropbox_root,"/GHG/Fieldsheets")




# ---- Import fluxes per season ----

#Read all computed fluxes per season
listf<- list.files(path = results_path, pattern = "^S[0-9]{1}", full.names = T)

table_results_all <- NULL
for (f in listf){
  table_results_all <- rbind(table_results_all, 
                             read.csv(file = f, header = T))
}


#Get fit details for every flux:
get_full_detailed_table <- function(table_results_all, variable){
  
  listf <- list.files(path = paste0(results_path,"level_incubation"), pattern = ".csv", all.files = T, full.names = T, recursive = F)
  mytable <- NULL
  for (f in listf[grep(pattern = variable, x = listf)]){
    mytable.tmp <- read.csv(file = f, header = T)
    mytable <- rbind(mytable, mytable.tmp)
  }
  mytable$sampling <- str_sub(mytable$UniqueID, start = 1, 2)
  mytable$pilotsite <- str_sub(mytable$UniqueID, start = 4, 5)
  mytable$subsite <- str_sub(mytable$UniqueID, start = 7, 8)
  mytable$siteID <- str_sub(mytable$UniqueID, start = 4, 8)
  
  ind_match <- match(mytable$UniqueID, table_results_all$UniqueID)
  
  
  mytable <- cbind(mytable, table_results_all[ind_match,c("gas_analiser","start.time","duration","water_depth","strata","chamberType","lightCondition")])
  mytable <- mytable[!is.na(ind_match),]
  return(mytable)
}

table_co2 <- get_full_detailed_table(table_results_all, variable = "co2")
table_ch4 <- get_full_detailed_table(table_results_all, "ch4")


#Get model MAE for each gas:
table_co2_MAE<- table_co2 %>% 
  select(UniqueID, HM.MAE) %>% 
  rename(HM.MAE.co2=HM.MAE)

table_ch4_MAE<- table_ch4 %>% 
  select(UniqueID, HM.MAE) %>% 
  rename(HM.MAE.ch4=HM.MAE)

#Add model MAE to full table: 
table_results_all<- merge.data.frame(table_results_all, table_ch4_MAE,by="UniqueID") %>% 
  merge.data.frame(table_co2_MAE,by="UniqueID" )




#Select and filter Fluxes####

#Selection of fluxes CO2: take CO2 flux from best.flux

#Selection of fluxes CH4:
#Choose best.flux CH4 as diffusion whenever water_depth==0 (i.e. NA for ebullition) or when ebullition==0  (if no ebullition or negative ebullition, goflux should perform better than Camille's method). When ebullition exists or is possitive, chose Camille diffussion estimate as CH4 diffusive flux and keep the ebullitive flux. 

# Summary for CH4 selection:
#Ebullition: if ebullition ==0 | ebullition== NA --> set ebullition to 0
#Diffusion: if ebullition !=0 and != NA --> set diffusion to diffusion
#if ebullition ==0 | ebullition==NA --> set difussion as best.flux from goflux method.
#total flux = diffussion + ebullition



#Filtering of fluxes: 

#for CO2: set to NA fluxes with HM.MAE > MAE_threshold_CO2

#for CH4: set to NA the fluxes from go.flux (water_depth==0|ebullition==0) with HM.MAE > MAE_threshold_CH4
#No filtering for CH4 fluxes coming from camilles method (ebullition!=NA, ebullition!=0). This means we will not filter out fluxes from camilles method if they produce a positive value for ebullition.


#Implementation of selection and filtering decissions: 
#Set unreliable fluxes to NA, and combine CH4 fluxes into diffusive and ebulitive flux, then into CH4_total
table_results_good<- table_results_all %>% 
  select(UniqueID,
         start.time,
         CO2_best.flux,HM.MAE.co2,
         CH4_best.flux, HM.MAE.ch4,
         CH4_diffusive_flux,CH4_ebullitive_flux) %>% 
  mutate(CO2_best.flux=case_when(HM.MAE.co2>MAE_threshold_CO2~NA_real_,#Set CO2 best flux with MAE >10 to NA (Unreliable)
                                 TRUE~CO2_best.flux),
         CH4_best.flux=case_when(HM.MAE.ch4>MAE_threshold_CH4~NA_real_,#Set CH4 best flux with MAE >100 to NA (Unreliable)
                                 TRUE~CH4_best.flux),
         CH4_ebullitive_flux=case_when(is.na(CH4_ebullitive_flux)~0, #substitute NAs ebullition with zero
                                       TRUE~CH4_ebullitive_flux),
         CH4_diffusive_flux=case_when(CH4_ebullitive_flux==0~CH4_best.flux, #when Camille method does not detect ebullition (ebullition ==0 or ebullition ==NA), take best.flux as diffusion
                                      TRUE~CH4_diffusive_flux),
         CH4_ebullitive_flux=case_when(is.na(CH4_diffusive_flux)~NA_real_, #set ebullition to NA, when CH4 does not produce a reliable diffusive flux 
                                       TRUE~CH4_ebullitive_flux)) %>% 
  #Keep only fluxes data (CO2_best.flux, CH4_total), start.time and UniqueID.
  mutate(CH4_total=CH4_diffusive_flux+CH4_ebullitive_flux) %>% 
  select(UniqueID,
         start.time,
         CH4_total,
         CO2_best.flux)


#Remove unnecesary data frames:
rm(table_ch4_MAE, table_co2_MAE, table_ch4, table_co2)


#How many fluxes have we removed with the above filtering?
sum(is.na(table_results_good$CO2_best.flux))
sum(is.na(table_results_good$CH4_total))
#very few fluxes removed, specially for CH4 (due to lack of filtering criteria for Camille's method ebullitionVSdiffusion)



#Import fieldsheet meta-data#----

#Read all GHGfieldsheets (using function from Camille)
fieldsheets<- list.files(path= fieldsheetpath, recursive = T, pattern = "GHG.xlsx",full.names = T)
field<- read_GHG_fieldsheets(fieldsheets)


#Drop unwanted variables and rename key UniqueID:
field2<- field %>% 
  select(-c(person_sampling,gas_analyzer,logger_floating_chamber,logger_transparent_chamber,start_time,initial_co2,initial_ch4,end_time,final_ch4,final_co2,chamber_height_cm, unix_start,unix_stop)) %>% 
  rename(UniqueID=uniqID)


# ---- Join fluxes and meta-data ----
#Key for joining flux data and fieldsheet: 
summary(table_results_good$UniqueID%in%field$uniqID)

#fluxes without fieldsheet details: arising from correction in fieldsheets
table_results_good[which(!table_results_good$UniqueID%in%field$uniqID),]

#fieldsheets without flux calculated (WHY?)
incub_noflux<-field[which(!field$uniqID%in%table_results_good$UniqueID),] %>% 
  mutate(duration=unix_stop-unix_start)

#17 incubations without a flux
#Reasons and checks in Incubations_without_flux_fixed.xlsx (same folder as LifeWatch data example)


#Check for duplicate incubations in fieldsheets (incubations performed 2 times with the same identity except for startime and uniqueID:  
duplicate_incubations<-field2 %>% 
  dplyr::summarise(n = dplyr::n(), .by = c(pilot_site, subsite, date, plot_id, chamber_type, strata, longitude, latitude, water_depth,
                                           transparent_dark)) |>dplyr::filter(n > 1L) %>% 
  merge.data.frame(field2, by=c("pilot_site", "subsite", "date", "plot_id", "chamber_type", "strata", "longitude", "latitude", "water_depth","transparent_dark"), all = F) %>% pull(UniqueID)
duplicate_incubations


#Check for missing transparent_dark
field2 %>% filter(is.na(transparent_dark))



#merge fluxes with fieldsheets dropping the unmatched observations, problematic values and duplicated incubations (no problematic values exists)
ghg_formated<- field2 %>% merge.data.frame(table_results_good, by="UniqueID",all.y = T) %>%
  filter(!is.na(transparent_dark)) %>% 
  filter(!UniqueID%in%duplicate_incubations) %>% 
  select(-UniqueID) %>% 
  pivot_wider(names_from = transparent_dark,values_from = c(start.time, CO2_best.flux, CH4_total,comments)) %>% 
  rowwise() %>% 
  mutate(plot_startime=min(start.time_dark, start.time_transparent, na.rm=T)) %>% 
  ungroup() %>% 
  mutate(plotcode=paste0(subsite,"-",plot_id,"-",toupper(substr(strata,1,1))),#Add first capitalized letter of strata to plotcode
         campaign=str_extract(subsite, pattern="S[0-9]{1}"),
         season=case_when(campaign=="S1"~"fall",
                          campaign=="S2"~"winter",
                          campaign=="S3"~"spring",
                          campaign=="S4"~"summer"),
         subsite=str_extract(subsite, pattern="[A-Z]{2}-[A-Z]{1}[0-9]{1}"),
         sampling=paste(campaign, subsite,sep="-"),
         strata=case_when(strata=="vegetated"&water_depth>0~"vegetated water",#Re-classify strata using water depth (vegetated water when water_depth>0cm, otherwise vegetated land)
                          strata=="vegetated"&water_depth==0~"vegetated land",
                          TRUE~strata)) %>% 
  select(season, campaign,pilot_site, subsite, sampling, date, latitude, longitude, water_depth, strata, plotcode, plot_startime, CO2_best.flux_dark, CO2_best.flux_transparent, CH4_total_dark,CH4_total_transparent)




#Integrate day-night ####

#Integrate day-night per strata

#Calculate daylight hours using latitude and date (hours of light)

# get median coordinates per subsite (including all seasons). Use of median avoids chambers with errors in coordinates afecting the calculation
samplings <- ghg_formated %>% 
  select(subsite, latitude, longitude) %>% 
  group_by(subsite) %>% 
  summarise(subsite_latitude=median(latitude,na.rm=T),
            subsite_longitude=median(longitude,na.rm=T)) %>% 
  ungroup()


#Use dplyr and suncalc to calculate daylight hours for each subsite-date combination.
sampling_daylight <- samplings %>%
  merge.data.frame(y=ghg_formated %>% select(subsite, date), by="subsite", all=T) %>% 
  distinct() %>% 
  mutate(date=as.Date(date)) %>% 
  rowwise() %>%
  mutate(
    sunlight_times = list(getSunlightTimes(date = date, lat = subsite_latitude, lon = subsite_longitude, keep = c("sunrise","sunset"))),
    daylight_duration = as.numeric(difftime(sunlight_times$sunset, sunlight_times$sunrise, units = "hours"))
  ) %>%
  ungroup() %>%
  select(subsite, date, subsite_latitude, subsite_longitude, daylight_duration)

print(sampling_daylight %>% 
        group_by(subsite) %>% 
        summarise(n=n()), n=50)

sampling_daylight %>% filter(subsite=="RI-R2")
#ISSUE: RI-R2 has been sampled in different days in S1 (november2023), S2 (january2024) and S3 (march2024)
# filter(!(subsite=="RI-R2"&date%in%c("2023-11-15","2024-01-24","2024-04-10"))) %>% #remove 1-day-early samplings of RI-R2 for daylighthour calculations





#INTEGRATE to net_GHG exchange per plot (integrate dark and transparent with daylight_duration for vegetated plots, keep dark bare and dark OW fluxes as representative of whole day.
net_ghg_per_plot<- ghg_formated %>% 
  mutate(date=as.Date(date)) %>% 
  merge.data.frame(sampling_daylight, by=c("subsite","date"), all = T) %>% 
  select(season, pilot_site, subsite, sampling, subsite_latitude, subsite_longitude, strata, date, daylight_duration, plotcode, latitude, longitude,  
         CO2_best.flux_dark, CO2_best.flux_transparent,
         CH4_total_dark, CH4_total_transparent) %>% 
  #Integrate transparent and dark with daylight_duration for vegetated fluxes(fluxlight*lighthours+fluxdark*darkhours)/24: same units for net flux (umol/m2/s), but integrated for the whole day. 
  #CO2:
  mutate(net_CO2=case_when(grepl("vegetated",strata)~((CO2_best.flux_transparent*daylight_duration)+(CO2_best.flux_dark*(24-daylight_duration)))/24,
                           #For Bare plots, take dark flux as representative for the whole day
                           grepl("bare",strata)~CO2_best.flux_dark,
                           #For open water plots, use dark flux for whole day
                           grepl("open",strata)~CO2_best.flux_dark),
         
         #CH4: same approach than for CO2 (scaled for vegetated plots, dark fluxes for the rest)
         net_CH4=case_when(grepl("vegetated",strata)~((CH4_total_transparent*daylight_duration)+(CH4_total_dark*(24-daylight_duration)))/24,
                           #For Bare plots, use dark flux for whole day
                           grepl("bare",strata)~CH4_total_dark,
                           #For open water plots, use dark flux for whole day
                           grepl("open",strata)~CH4_total_dark)) %>% 
  #Discard _transparent/_dark ghg fluxes
  select(-c(CO2_best.flux_dark, CO2_best.flux_transparent,CH4_total_dark,CH4_total_transparent)) 


#Summary per subsite####

#Summarise day-integrated fluxes without strata or seasonal groupings, save in export_path.
#Keep median lat/long per subsite 

net_ghg_summary_persubsite<- net_ghg_per_plot %>% 
  select(-c(season, pilot_site, sampling, strata, daylight_duration, plotcode,latitude, longitude,date)) %>% 
  group_by(subsite ,subsite_latitude,subsite_longitude) %>% 
  summarise(avg_netCO2=mean(net_CO2, na.rm=T),
            sd_netCO2=sd(net_CO2, na.rm=T),
            n_netCO2=sum(!is.na(net_CO2)),
            avg_netCH4=mean(net_CH4, na.rm=T),
            sd_netCH4=sd(net_CH4, na.rm=T),
            n_netCH4=sum(!is.na(net_CH4)),
            .groups = "drop"
  )


#Save net_ghg per subsite_with explicit date of last update.
write.csv(net_ghg_summary_persubsite,file = paste0(export_path, "net_ghg_summary_persubsite_V",dayoflastrun, ".csv"),row.names = F)


