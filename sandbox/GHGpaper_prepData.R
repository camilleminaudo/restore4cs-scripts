#GHGpaper_prepchamberdata.R


#Author: Miguel Cabrera-Brufau
#Date: September 2025
#Project: Restore4cs



#Description----
#This scrip is used to prepare the datasets that will be used for statistical analysis and figures of the main GHG paper that will describe the effects of wetland restoration on the fluxes of GHGs using net daily GHG exchange rates of appropriate incubations derived from in-situ measurements.

#MAIN STEPS of DATA PREPARATION: 
#1. Import all fluxes and combine in single table (UniqueID, ghg_best, ghg_se, ghg_model) 
#2. Import harmonized field observations. 
#3. Discard non-appropriate fluxes (rising/receding tide, after veg cut, apropriate strata per subsite only)
#4. Calculate persampling-strata areal composition (using distribution of appropriate chamber deployments (plotcode strata) as aproximation, including deployments without valid flux, to estimate strata composition in subsite). These per-sampling strata proportion composition will be overriden by Remote-sensing derived composition when they become available. 
#5. Integrate day-night fluxes,calculate global warming potential c_co2eq (ch4+co2)
#6. Check completeness fluxes (and lost fluxes due to integration approach for vegetation ch4 and n2o)
#7. Save datasets for further analysis. 


#Inputs: 
#Per-UniqueID: co2, ch4 and n2o best fluxes (and se flux)
#Per-UniqueID and per-plotcode harmonized field observations (created/Updated by  script Harmonize_Fieldsheet-veg_perPLOT_and_perUniqueID). 
#Strata-representativity decissions: common and relevant strata for comparisons.
#Per-Core_id co2, ch4 and n2o fluxes (and se flux)


#Outputs: 
#ChamberData4paper.csv: containing valid daily NEE of co2, ch4 and combined GWP (co2+ch4). N2O fluxes will not be used (only 1 season, often non-significant).
#Stratacomposition_in-situ.csv: containing the strata areal composition derived from chamber deployments. 
#CoresData4paper.csv: formatted dailyflux co2, ch4, n2o, GWPco2andch4, GWPtotal (co2+ch4+n2o) obtained from incubations of cores. 

rm(list = ls()) # clear workspace
cat("/014") # clear console


# ---- Packages ----
#For data-handling & plotting
library(tidyverse)
library(readxl)
library(lubridate)
library(zoo)
# library(ggplot2)
# library(ggpubr)
# library(grid)
# library(egg)
require(purrr)
require(data.table)
require(tools)
library(hms)
library(suncalc)
# library(cowplot) #Alligned axis
# library(ggforce)


#Repo-functions
repo_root <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
files.sources = list.files(path = paste0(repo_root,"/functions"), full.names = T)
for (f in files.sources){source(f)}




# ---- Directories ----
dropbox_root <- "C:/Users/Miguel/Dropbox/RESTORE4Cs - Fieldwork/Data" # You have to make sure this is pointing to the write folder on your local machine

#Path to aquaGHG best.flux results:
results_path <- paste0(dropbox_root,"/GHG/Processed data/computed_flux/")
#Path to N2O flux results:
n2o_results_path <- paste0(dropbox_root,"/GHG/N2O_fluxes/")

#Path to save datasets for main paper: 
paper_path<- paste0(dropbox_root,"/GHG/Main_paper_results/")


#Path to cores dataset: 
cores_path<- paste0(dropbox_root,"/Cores/Cores_flux/")


#__________________------
#Chamber DATA PREP------

#1. Import inst. fluxes-----

#All datasets must be identified at the UniqueID level
co2_all<- read.csv(paste0(results_path,"co2_bestflux.csv")) #umol/m2/s
ch4_all<- read.csv(paste0(results_path,"ch4_bestflux.csv")) #nmol/m2/s

#Check N2O identification scheme: UniqueID_notime.
n2o_all<- read.csv(paste0(n2o_results_path,"S4_restore4cs_N2O_arealflux_nmol_s-1_m-2.csv"))


#Select best fluxes and format: UniqueID, ghg_best, ghg_se
#CO2
co2_best<- co2_all %>% 
  mutate(co2_best=case_when(best.model=="LM"~LM.flux,
                            best.model=="HM"~HM.flux,
                            best.model=="total.flux"~total.flux,
                            best.model=="None  appropriate"~NA_real_),
         co2_se=case_when(best.model=="LM"~LM.flux.se,
                          best.model=="HM"~HM.flux.se,
                          best.model=="total.flux"~total.flux.se,
                          best.model=="None  appropriate"~NA_real_),
         co2_model=best.model) %>% 
  dplyr::select(UniqueID, co2_best, co2_se,co2_model)

#CH4
ch4_best<- ch4_all %>% 
  mutate(ch4_best=case_when(best.model=="LM"~LM.flux,
                            best.model=="HM"~HM.flux,
                            best.model=="total.flux"~total.flux,
                            best.model=="None  appropriate"~NA_real_),
         ch4_se=case_when(best.model=="LM"~LM.flux.se,
                          best.model=="HM"~HM.flux.se,
                          best.model=="total.flux"~total.flux.se,
                          best.model=="None  appropriate"~NA_real_),
         ch4_model=best.model) %>% 
  dplyr::select(UniqueID, ch4_best, ch4_se,ch4_model)

#N2O
n2o_best<- n2o_all %>%
  mutate(n2o_best=N2Oflux_nmol_per_second_per_m2,
         n2o_se=N2Oflux_se) %>%
  dplyr::select(UniqueID_notime, n2o_best, n2o_se)


#JOIN all ghg best-fluxes per UniqueID
ghg_best<- co2_best %>% 
  separate(UniqueID, into = c("c1","c2","c3","c4","c5","c6","c7"), sep = "-", remove = F) %>% 
  mutate(UniqueID_notime=paste(c1,c2,c3,c4,c5,c6,sep = "-")) %>% 
  dplyr::select(-c(c1,c2,c3,c4,c5,c6,c7)) %>% 
  merge.data.frame(ch4_best, by="UniqueID", all=T) %>% 
  merge.data.frame(n2o_best, by="UniqueID_notime", all=T) %>%
  dplyr::select(-UniqueID_notime)

rm(co2_all, co2_best, ch4_all, ch4_best, n2o_all, n2o_best)

#2. Import harm. field Obs.------

#Harmonized field observations per UniqueID
fieldobs_uniqueID<-read.csv(paste0(results_path,"UniqueID_harmonized_field_obs.csv"))

#Harmonized field observations per plotcode 
fieldobs_plotcode<- read.csv(paste0(results_path,"plotcode_harmonized_field_obs.csv"))

#3. Filter valid incubations-----

#Remove fieldobs from non-valid/non-appropriate incubations: Incubations that are not representative of subsite definitions, incubations that are not representative of natural conditions

# 1. Remove rising-tide and receding-tide incubations. (non appropriate fluxes for comparisons)
tidal_IDs<- fieldobs_uniqueID %>% filter(grepl("rising-tide|receding-tide", field_observations)) %>% pull(UniqueID)

# 2. Remove "bare after vegetation removal" (non-appropriate for comparisons)  
vegcut_IDs<- fieldobs_uniqueID %>% filter(grepl("after vegetation removal", field_observations)) %>% pull(UniqueID)


# 3. Remove incubations in non-appropriate strata (according to subsite definition)
#3.1. Camargue (CA): NOTHING TO REMOVE, all incubations are in appropriate strata
CA_wrongstrata_IDs<- c()

#3.2. Curonian lagoon (CU): inconsistent bare representation (only sampled consistently in restored sites although all sites had bare area, non-appropriate for comparisons). We will only consider vegetated and open water plots
CU_wrongstrata_IDs<- fieldobs_uniqueID %>% filter(casepilot=="CU") %>% filter(strata=="bare") %>% pull(UniqueID)

#3.3. Danube delta (DA): inconsistent strata in some samplings: exclude vegetated of S3-DA-A2 (not representative), exclude bare of S1-DA-R2 (only sampled in first visit)
DA_wrongstrata_IDs<- fieldobs_uniqueID %>% filter(casepilot=="DA") %>% 
  filter((sampling=="S3-DA-A2"&strata=="vegetated")|(sampling=="S1-DA-R2"&strata=="bare")) %>% 
  pull(UniqueID)

#3.4. Dutch delta (DU): NOTHING TO REMOVE, all incubations are in appropriate strata
DU_wrongstrata_IDs<- c()

#3.5. Ria d'Aveiro (RI): rising-tide and aftervegcut already removed. Inconsistent incubations: tidal-pool incubations in S3-RI-A1 (only sampled in S3), bare in RI-R1 and in RI-R2 (restored sites are 100% vegetated by definition)
RI_wrongstrata_IDs<- fieldobs_uniqueID %>% filter(casepilot=="RI") %>% 
  filter((sampling=="S3-RI-A1"&grepl("tidal-pool", field_observations))|(subsite=="RI-R1"&strata=="bare")|(subsite=="RI-R2"&strata=="bare")) %>% pull(UniqueID)

#3.6. Valencian wetland (VA): NOTHING TO REMOVE, all incubations are in appropriate strata
VA_wrongstrata_IDs<-c()


#ONLY valid incubation details
valid_fieldobs_uniqueID<- fieldobs_uniqueID %>% 
  filter(!UniqueID%in%c(tidal_IDs,vegcut_IDs,CA_wrongstrata_IDs,CU_wrongstrata_IDs,DA_wrongstrata_IDs,DU_wrongstrata_IDs,RI_wrongstrata_IDs,VA_wrongstrata_IDs))

rm(tidal_IDs,vegcut_IDs,CA_wrongstrata_IDs,CU_wrongstrata_IDs,DA_wrongstrata_IDs,DU_wrongstrata_IDs,RI_wrongstrata_IDs,VA_wrongstrata_IDs)

#4. Strata composition ----

#Construct per-sampling strata composition from chamber deployment: this will be replaced by RS-derived strata composition (when it becomes available)

#Dataset format: season, casepilot,subsite,sampling,strata, strata_prop(0-1)

#strata_prop is calculated for each sampling based on strata distribution of chambers non-excluded by design (those still in valid_plotcodes)
valid_plotcodes<- valid_fieldobs_uniqueID %>% 
  dplyr::select(plotcode, season, casepilot,subsite,sampling, strata) %>% distinct()

#Calculation procedure might vary depending on case pilot, see below: 

#CA: use chamber distribution directly.
ca_stratacomp<- valid_plotcodes %>% filter(casepilot=="CA") %>% 
  #get chambers per strata
  group_by(season,casepilot,subsite,sampling,strata) %>%
  summarise(n_deployed=n()) %>%
  #Calculate proportion of strata
  group_by(season,casepilot,subsite,sampling) %>% 
  mutate(totalchambers=sum(n_deployed),
         strata_prop=n_deployed/totalchambers) %>% 
  dplyr::select(season, casepilot,subsite,sampling,strata,strata_prop) %>%
  ungroup()


#CU: inconsistent boundaries, no seasonal change. All subsites will get the same strata composition (calculated with all avaliable chambers, not grouped)
cu_stratacomp<- valid_plotcodes %>% filter(casepilot=="CU") %>% 
  #Obtain average strata_prop (from all available chambers, no grouping)
  group_by(casepilot,strata) %>%
  summarise(n_deployed=n()) %>% 
  ungroup() %>% 
  mutate(totalchambers=sum(n_deployed),
         strata_prop=n_deployed/totalchambers) %>% 
  dplyr::select(casepilot,strata,strata_prop) %>%
  ungroup() %>% 
  #expand common stratacomp to all subiste samplings
  merge.data.frame(valid_plotcodes %>% filter(casepilot=="CU") %>% dplyr::select(season,casepilot,subsite,sampling,strata) %>% distinct(), by=c("casepilot","strata"),all = T)%>% 
  dplyr::select(season, casepilot,subsite,sampling,strata,strata_prop)


#DA: same strata composition for all preserved and restored subsites (P1,P2,R1,R2), using all chambers, not grouped. DA-A1 and DA-A2 have seasonal variability and distinct strata-composition, use chamber distribution directly. 

da_a_stratacomp<- valid_plotcodes %>% filter(subsite%in%c("DA-A1","DA-A2")) %>% 
  #get chambers per strata
  group_by(season,casepilot,subsite,sampling,strata) %>%
  summarise(n_deployed=n()) %>%
  #Calculate proportion of strata
  group_by(season,casepilot,subsite,sampling) %>% 
  mutate(totalchambers=sum(n_deployed),
         strata_prop=n_deployed/totalchambers) %>% 
  dplyr::select(season, casepilot,subsite,sampling,strata,strata_prop) %>%
  ungroup()

da_pr_stratacomp<- valid_plotcodes %>% filter(subsite%in%c("DA-P1","DA-P2","DA-R1","DA-R2")) %>% 
  #Calculate general strata composition for all chamber deployments
  group_by(casepilot,strata) %>%
  summarise(n_deployed=n()) %>% 
  ungroup() %>% 
  mutate(totalchambers=sum(n_deployed),
         strata_prop=n_deployed/totalchambers) %>% 
  dplyr::select(casepilot,strata,strata_prop) %>%
  ungroup() %>% 
  #expand common strata composition to each subsite
  merge.data.frame(valid_plotcodes %>% filter(subsite%in%c("DA-P1","DA-P2","DA-R1","DA-R2")) %>% dplyr::select(season,casepilot,subsite,sampling,strata) %>% distinct(), by=c("casepilot","strata"),all = T)%>% 
  dplyr::select(season, casepilot,subsite,sampling,strata,strata_prop)


#DU: use chamber distribution directly.
du_stratacomp<- valid_plotcodes %>% filter(casepilot=="DU") %>% 
  #get chambers per strata
  group_by(season,casepilot,subsite,sampling,strata) %>%
  summarise(n_deployed=n()) %>%
  #Calculate proportion of strata
  group_by(season,casepilot,subsite,sampling) %>% 
  mutate(totalchambers=sum(n_deployed),
         strata_prop=n_deployed/totalchambers) %>% 
  dplyr::select(season, casepilot,subsite,sampling,strata,strata_prop) %>%
  ungroup()


#RI: use chamber distribution directly (after removal of non-valid incubations)
ri_stratacomp<- valid_plotcodes %>% filter(casepilot=="RI") %>% 
  #get chambers per strata
  group_by(season,casepilot,subsite,sampling,strata) %>%
  summarise(n_deployed=n()) %>%
  #Calculate proportion of strata
  group_by(season,casepilot,subsite,sampling) %>% 
  mutate(totalchambers=sum(n_deployed),
         strata_prop=n_deployed/totalchambers) %>% 
  dplyr::select(season, casepilot,subsite,sampling,strata,strata_prop) %>%
  ungroup()

#Bare incubations in S3-RI-P1? (benj + antonio), nothing weird, not after veg-cut, taken as true seasonal variation of subsite. BUT only present in S3 sampling... Impact for model?


#VA: use chamber distribution directly with the exception of VA-R1 (no seasonal change): all samplings in subsite VA-R1 will get the same strata composition (average of all seasons)
#VA-R1:
va_r1_stratacomp<- valid_plotcodes %>% filter(subsite%in%c("VA-R1")) %>% 
  #get chambers per strata
  group_by(casepilot,subsite,strata) %>%
  summarise(n_deployed=n()) %>%
  #Calculate proportion of strata
  group_by(casepilot,subsite) %>% 
  mutate(totalchambers=sum(n_deployed),
         strata_prop=n_deployed/totalchambers) %>% 
  dplyr::select(casepilot,subsite,strata,strata_prop) %>%
  ungroup() %>% 
  #expand common strata composition to each sampling
  merge.data.frame(valid_plotcodes %>% filter(subsite%in%c("VA-R1")) %>% dplyr::select(season,casepilot,subsite,sampling,strata) %>% distinct(), by=c("casepilot","subsite","strata"),all = T)%>% 
  dplyr::select(season, casepilot,subsite,sampling,strata,strata_prop)

#Rest of VA: 
va_notR1_stratacomp<- valid_plotcodes %>% filter(subsite%in%c("VA-A1","VA-A2","VA-P1","VA-P2","VA-R2")) %>% 
  #get chambers per strata
  group_by(season,casepilot,subsite,sampling,strata) %>%
  summarise(n_deployed=n()) %>%
  #Calculate proportion of strata
  group_by(season,casepilot,subsite,sampling) %>% 
  mutate(totalchambers=sum(n_deployed),
         strata_prop=n_deployed/totalchambers) %>% 
  dplyr::select(season, casepilot,subsite,sampling,strata,strata_prop) %>%
  ungroup()


#Join and format stratacomp:
stratacomp<- rbind(ca_stratacomp,cu_stratacomp,da_a_stratacomp,da_pr_stratacomp,du_stratacomp,ri_stratacomp,va_notR1_stratacomp,va_r1_stratacomp) %>% 
  #make 0 strata proportion explicit (pivot_wider-->pivot_longer)
  pivot_wider(names_from = strata,values_from=strata_prop ,values_fill = 0) %>% 
  pivot_longer(cols = c(bare,vegetated,`open water`), names_to = "strata",values_to = "strata_prop")

rm(ca_stratacomp,cu_stratacomp,da_a_stratacomp,da_pr_stratacomp,du_stratacomp,ri_stratacomp,va_notR1_stratacomp,va_r1_stratacomp)





#WEIRD patterns/doubts strata_composition: 
#RI preserved, some bare plots (very few), RI-P2: only 1 bare plot in first season. RI-P1: 4 bare plots, all during S3 sampling.



#5. Integrate Daily fluxes ------
#From ghg_best (uniqueID) to ghg_daybest (plotcode)


##5.1. Get daylight hours----

#daylight hours are defined as time from sunrise to sunset (using package suncalc, median lat/long of each subsite and date of sampling)

# get median coordinates per subsite (including all seasons)
samplings <- valid_fieldobs_uniqueID %>% 
  dplyr::select(subsite, latitude, longitude) %>% 
  group_by(subsite) %>% 
  summarise(subsite_latitude=median(latitude,na.rm=T),
            subsite_longitude=median(longitude,na.rm=T)) %>% 
  ungroup()

#Use dplyr and suncalc to calculate daylight hours for each subsite-date combination (each sampling).
sampling_daylight <- samplings %>%
  merge.data.frame(y=valid_fieldobs_uniqueID %>% dplyr::select(subsite, date,sampling), by="subsite", all=T) %>% 
  distinct() %>% 
  mutate(date=as.Date(date)) %>% 
  rowwise() %>%
  mutate(
    sunlight_times = list(getSunlightTimes(date = date, lat = subsite_latitude, lon = subsite_longitude, keep = c("sunrise","sunset"))),
    daylight_duration = as.numeric(difftime(sunlight_times$sunset, sunlight_times$sunrise, units = "hours"))
  ) %>%
  ungroup() %>%
  mutate(hours_day=daylight_duration,
         hours_night= 24-hours_day) %>% 
  dplyr::select(sampling, hours_day, hours_night)




## 5.2. Scale to daily flux-----

#Daily integration will follow the previously agreed decisions (common for all gas species): 
#Vegetated fluxes will be systematically scaled with daylight hours (T+D)
#Open water fluxes will be directly used (only D)
#Bare fluxes will be non-systematically scaled according to data availability: 
#Dark flux only--> dark flux = daily flux
#Transparent only --> transparent flux = daily flux
#T+D available--> scaling with daylight hours (T+D) = daily flux

#NOTE: Inconsistency in scaling of bare fluxes responds to exploration of data, where light condition has seen not as important as plot-to-plot variability and wet/dry sediment, and where discarding transparent bare fluxes will severely limit the representativity of this stratum. 

#get field details of valid plotcodes
valid_plotcode_obs<- fieldobs_plotcode %>%
  filter(plotcode%in%valid_plotcodes$plotcode) 

#Get ghg fluxes for transparent & dark by plotcode (only for valid UniqueIDs)
valid_ghg_transpdark<- valid_fieldobs_uniqueID %>%
  merge.data.frame(ghg_best, by="UniqueID", all.x = T) %>% 
  dplyr::select(plotcode, co2_best, ch4_best,n2o_best, transparent_dark) %>% 
  rename(co2=co2_best, ch4=ch4_best, n2o=n2o_best) %>% 
  pivot_longer(cols=c(co2,ch4,n2o), names_to = "ghg", values_to = "flux") %>% 
  pivot_wider(names_from = transparent_dark, values_from = flux)


#Integrate to dailyflux
daily_ghg_plotcode<- valid_plotcode_obs %>% 
  #Add daylight hours
  merge.data.frame(sampling_daylight, by=c("sampling"), all = T) %>% 
  dplyr::select(plotcode, season, casepilot, subsite, sampling, plot_num, date, plot_start_time, latitude, longitude, gas_analyzer, strata, water_depth, field_observations, ABG_veg_description, ABG_veg_biomass_gpersquaremeter, hours_day, hours_night) %>% 
  #Add valid_ghg_transpdark
  merge.data.frame(valid_ghg_transpdark, by="plotcode", all.x = T) %>% 
  #Integrate dailyflux: change units to per day (from per s to per hour, then scale with daylight_duration in hours)
  #For all vegetated plots, scale flux with daylight
  mutate(dailyflux=case_when(strata=="vegetated"~(dark*60*60*hours_night)+(transparent*60*60*hours_day),
                             #For bare with both fluxes, scale with daylight
                             (strata=="bare"&!is.na(dark*transparent))~(dark*60*60*hours_night)+(transparent*60*60*hours_day),
                             #For bare with NA in transparent, use dark for all day.
                             (strata=="bare"&is.na(transparent))~dark*60*60*24,
                             #For bare with NA in dark, use transparent for all day
                             (strata=="bare"&is.na(dark))~transparent*60*60*24,
                             #For open water, use always  dark for all day
                             strata=="open water"~dark*60*60*24)) %>% 
  #Discard transparent & dark ghg fluxes
  dplyr::select(-c(dark, transparent)) %>% 
  pivot_wider(names_from = ghg, values_from = dailyflux) %>% 
  #Transform molar units for daily flux: 
  #CO2: umol per m2 per day --> mol per m2 per day
  #CH4: nmol per m2 per day --> mmol per m2 per day
  #N2O: nmol per m2 per day --> mmol per m2 per day
  mutate(co2_mol=co2*1e-6,
         ch4_mmol=ch4*1e-6,
         n2o_mmol=n2o*1e-6,
         #Global warming potential GWP100: 
         gwp_co2=co2_mol*44.0095, #molar flux* molar mass-->g/m2/day
         gwp_ch4=(ch4_mmol*1e-3)*16.04246*27,#molar flux * molar mass* GWP100 factor 27
         gwp_n2o=(n2o_mmol*1e-3)*44.0128*273, #molar flux * molar mass 44* GWP100factor 273
         gwp_co2andch4= gwp_co2+gwp_ch4, #Sum of gwp of co2 and ch4 only
         gwp_allghg=gwp_co2+gwp_ch4+gwp_n2o) %>%  #Sum of gwp of all GHGs (co2,ch4,n2o)
  dplyr::select(-c(co2, ch4, n2o)) %>% 
  #Add modeling variables: 
  mutate(status=substr(x = sampling, start = 7, stop=7),
         status=case_when(status=="A"~"Altered",
                          status=="P"~"Preserved",
                          status=="R"~"Restored"))

#Delete objects already used
rm(samplings, sampling_daylight, fieldobs_plotcode, fieldobs_uniqueID, valid_plotcode_obs, valid_plotcodes)






#6. Final format -----

#Produce final dataset to be used in models. Implementing strata-weigths and formating to simplify filtering for modeling. 

#We will use molar fluxes for individual ghgspecies and CO2eq (gCo2eq /m2/d) for GWP of combined ghgspecies. Add column with "units" to be able to differenciate. 


#We format the dataset to be easily used across our modeling approaches: 
#categories are taken as factors, we pivot longer the daily fluxes (response variable), we add the calculation of sample_weight based on strata proportion and sample proportion (non-NA) across strata (corrects sampling bias and ensures equal influence across samplings in the models)

data4models <- daily_ghg_plotcode %>% 
  dplyr::select(plotcode, season, casepilot, status, subsite, sampling, strata, 
                co2_mol, #mol per m2 per day
                ch4_mmol, #mmol per m2 per day
                n2o_mmol, #mmol per m2 per day
                gwp_co2andch4, #gCO2eq per m2 per day
                gwp_allghg) %>% #gCO2eq per m2 per day
  merge.data.frame(stratacomp, by = c("season","casepilot","subsite","sampling","strata"), all.x = TRUE) %>% 
  pivot_longer(
    cols = c(co2_mol, ch4_mmol, n2o_mmol, gwp_co2andch4, gwp_allghg),
    names_to = "ghgspecies",
    values_to = "dailyflux"
  ) %>% 
  #Add unitflux variable
  mutate(unitflux=case_when(ghgspecies=="co2_mol"~"mol per m2 per day",
                            ghgspecies=="ch4_mmol"~"mmol per m2 per day",
                            ghgspecies=="n2o_mmol"~"mmol per m2 per day",
                            ghgspecies=="gwp_co2andch4"~"gCO2eq per m2 per day",
                            ghgspecies=="gwp_allghg"~"gCO2eq per m2 per day")) %>% 
  #Remove units from ghgspecies:
  mutate(ghgspecies=case_when(ghgspecies=="co2_mol"~"co2",
                              ghgspecies=="ch4_mmol"~"ch4",
                              ghgspecies=="n2o_mmol"~"n2o",
                              ghgspecies=="gwp_co2andch4"~"gwp_co2andch4",
                              ghgspecies=="gwp_allghg"~"gwp_allghg")) %>% 
  mutate(across(c(season, casepilot, status, subsite, sampling, strata), as.factor)) %>% 
  
  # Calculate non-NA counts per stratum × sampling × ghg species
  group_by(sampling, strata, ghgspecies) %>% 
  mutate(nonNA_n_perstrata = sum(!is.na(dailyflux))) %>% 
  ungroup() %>%
  
  # Calculate total non-NA count per sampling × ghg species
  group_by(sampling, ghgspecies) %>% 
  mutate(
    nonNA_n_persampling = sum(!is.na(dailyflux)),
    sample_prop = nonNA_n_perstrata / nonNA_n_persampling,
    sample_weight_raw = strata_prop / sample_prop
  ) %>%
  ungroup() %>%
  # Rescale weights to sum to N (non-NA count) per group: i.e. obtain 
  group_by(sampling, ghgspecies) %>%
  mutate(sample_weight = ifelse(
    is.na(dailyflux), NA,
    sample_weight_raw * sum(!is.na(dailyflux)) / sum(sample_weight_raw[!is.na(dailyflux)])
  )) %>%
  ungroup() %>% 
  # Drop intermediate vars
  dplyr::select(-c(sample_prop, nonNA_n_persampling, nonNA_n_perstrata, sample_weight_raw, strata_prop))


#With the above weighting, we can account for sampling bias across strata at each sampling, but we are not equalizing the overall influence of each sampling. This is the most appropriate approach, given that glmmTMB treats weights as precission weights (all weights assumed as 1 when not specified otherwise), so the total sum of weights has to equal N within each sampling. If we try to equalize the influence of each sampling (i.e. overriding the potential effect of imbalanced N across samplings), the weights do not work in the intended way and the output from the models is not comparable (with and without sample weights). However, AIC is still not perfectly comparable between weighted and unweighted models, we should complement AIC with Predictive performance, Parameter stability, Residual diagnostics and Marignal&Conditional R2. 

#Check: do all samplings with valid (non-na) observations have a sum of sample_weights that equals their n?
data4models %>% 
  filter(!is.na(dailyflux)) %>% 
  group_by(sampling, ghgspecies) %>% 
  summarise(sumweights=sum(sample_weight), nobs=n(),.groups = "drop") %>% 
  filter(round(sumweights,digits = 1)!=nobs)#Filter for samplings*ghg combos without weight!=nobs

#OK, all samplings have weights thta match N. 

#7.Save datasets--------

#Save formated and valid in-situ data to be used in paper. 
write.csv(x = data4models, file = paste0(paper_path,"ChamberData4paper.csv"),row.names = F)

#Save stratacomposition derived from appropriate chamber deployments and sampling design decissions. 
write.csv(x = stratacomp, file = paste0(paper_path,"Stratacomposition_in-situ.csv"),row.names = F)


#__________________------

#Cores DATA PREP------


#1. Import cores fluxes -----
coreflux_wide<- read.csv(paste0(cores_path,"All_core_GHGfluxes_and_CH4stocks.csv"))


head(coreflux_wide)

#2. Format core fluxes-----
#Format cores to paper-format:
Coresdata4models<- coreflux_wide %>% 
  dplyr::select(core_id,CO2_flux_mmol_per_m2_per_d,CH4_flux_micromol_per_m2_per_d,N2O_flux_micromol_per_m2_per_d) %>% 
  #ADD relevant qualitative variables for stats: 
  separate(core_id, into = c("season", "casepilot","status_num","corenum"),remove = F,sep = "-") %>% 
  mutate(subsite=paste(casepilot, status_num, sep = "-"),
         status_letter=substr(status_num,1,1),
         status=case_when(status_letter=="A"~"Altered",
                          status_letter=="P"~"Preserved",
                          status_letter=="R"~"Restored"),
         #Add status_difalter (to be able to differentiate DA alterations)
         status_difalter=case_when(subsite=="DA-A1"~"Altered 1",
                                   subsite=="DA-A2"~"Altered 2",
                                   TRUE~status),
         sampling=paste(season, subsite, sep = "-")) %>% 
  dplyr::select(-c(status_letter, status_num, corenum)) %>% 
  #HOMOGENIZE UNITS: co2 (mmol per m2 per day), ch4 (micromol per m2 per day), n2o (micromol per m2 per day), GWPtotal ( g of CO2eq per m2 per day), GWPco2andch4 (g of co2eq per m2 per day)
  mutate(
    #Obtain mol per m2 per day: 
    co2_mol=CO2_flux_mmol_per_m2_per_d*1e-3,
    ch4_mol=CH4_flux_micromol_per_m2_per_d*1e-6,
    n2o_mol=N2O_flux_micromol_per_m2_per_d*1e-6,
    #Calculate Global warming potential (GWP100):  g of CO2eq per m2 per day
    gwp_co2=co2_mol*44.0095, #molar flux * molar mass* GWP100factor(1)
    gwp_ch4=ch4_mol*16.04246*27,#molar flux * molar mass* GWP100factor(27)
    gwp_n2o=n2o_mol*44.0128*273, #molar flux * molar mass * GWP100factor(273)
    GWPco2andch4_gCO2.per.m2.per.day= gwp_co2+gwp_ch4, #Sum of gwp of co2 and ch4 only
    GWPtotal_gC2eq.per.m2.per.day=gwp_co2+gwp_ch4+gwp_n2o) %>%  #Sum of gwp of all GHGs (co2,ch4,n2o)
  rename(corecode=core_id, co2_mmol.per.m2.per.day=CO2_flux_mmol_per_m2_per_d,
         ch4_micromol.per.m2.per.day=CH4_flux_micromol_per_m2_per_d,
         n2o_micromol.per.m2.per.day=N2O_flux_micromol_per_m2_per_d,
  ) %>% 
  #Organize variables:
  dplyr::select(season,casepilot, subsite, sampling,corecode,status,status_difalter,
                co2_mmol.per.m2.per.day,
                ch4_micromol.per.m2.per.day,
                n2o_micromol.per.m2.per.day,
                GWPco2andch4_gCO2.per.m2.per.day,
                GWPtotal_gC2eq.per.m2.per.day) %>%
  #long format:
  pivot_longer(cols = c(co2_mmol.per.m2.per.day,
                        ch4_micromol.per.m2.per.day,
                        n2o_micromol.per.m2.per.day,
                        GWPco2andch4_gCO2.per.m2.per.day,
                        GWPtotal_gC2eq.per.m2.per.day), names_to = c("ghgspecies", "unitflux"), names_sep = "_", values_to = "dailyflux") %>% 
  mutate(unitflux=gsub("\\."," ",unitflux))


#3. Save cores dataset-----
#core4paper as formated csv in paper_path
write.csv(Coresdata4models, paste0(paper_path, "CoresData4paper.csv"),row.names = F)

#__________________------
