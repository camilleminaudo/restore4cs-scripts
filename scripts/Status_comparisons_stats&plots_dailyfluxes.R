#Stats and plots: effects of restoration of coastal wetlands on GHG fluxes

#Author: Miguel Cabrera-Brufau
#Date: June 2025
#Project: Restore4cs


#Description----
#This scrip is used to explore the main effects of conservation status (altered, preserved, restored). Using net daily GHG exchange rates of apropriate incubations. 

#Main steps: 

#DATA PREPARATION: 
#1. Import all fluxes and combine in single table (UniqueID, ghg_best, ghg_se, ghg_model) 
#2. Import harmonized field observations. 
#3. Discard non-appropriate fluxes (rising/receding tide, after veg cut, apropriate strata per subsite only)
#4. Calculate persampling-strata areal composition (using distribution of appropriate chamber deployments (plotcode strata) as aproximation, including deployments without valid flux, to estimate strata composition in subsite). These per-sampling strata proportion composition will be overriden by Remote-sensing derived composition when they become available. 
#5. Integrate day-night fluxes,calculate global warming potential c_co2eq (ch4+co2) and ghg_co2eq (ch4+n2o+co2)

#6. Check completeness fluxes (and lost fluxes due to integration approach for vegetation ch4 and n2o)

#7. Outlier identification& NA handling

#STATS & PLOTS
#PLots and LMM structure (general and per-case pilot): following decissions  


#Inputs: 
#Per-UniqueID: co2, ch4 and n2o best fluxes (and se flux)
#Per-UniqueID and per-plotcode harmonized field observations (created/Updated by  script Harmonize_Fieldsheet-veg_perPLOT_and_perUniqueID). 
#Strata-representativity decissions: common and relevant strata for comparisons


rm(list = ls()) # clear workspace
cat("/014") # clear console


# ---- Packages ----
library(tidyverse)
library(readxl)
library(lubridate)
library(zoo)
library(ggplot2)
library(ggpubr)
library(grid)
library(egg)
require(dplyr)
require(purrr)
require(data.table)
require(tools)
library(hms)
library(suncalc)

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

#Path to save plots from exploration: 
# plots_path <- paste0(dropbox_root,"/GHG/Processed data/computed_flux/plots")

#__________________------
#DATA PREP------

#1. Import inst. fluxes-----
#1. Import bestflux results

#All datasets must be identified at the UniqueID level
co2_all<- read.csv(paste0(results_path,"co2_bestflux.csv")) #umol/m2/s
ch4_all<- read.csv(paste0(results_path,"ch4_bestflux.csv")) #nmol/m2/s

#Check N2O identification scheme: UniqueID_notime
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
  select(UniqueID, co2_best, co2_se,co2_model)

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
  select(UniqueID, ch4_best, ch4_se,ch4_model)

#N2O
n2o_best<- n2o_all %>% 
  mutate(n2o_best=N2Oflux_nmol_per_second_per_m2,
         n2o_se=N2Oflux_se) %>% 
  select(UniqueID_notime, n2o_best, n2o_se)


#JOIN all ghg best-fluxes per UniqueID
ghg_best<- co2_best %>% 
  separate(UniqueID, into = c("c1","c2","c3","c4","c5","c6","c7"), sep = "-", remove = F) %>% 
  mutate(UniqueID_notime=paste(c1,c2,c3,c4,c5,c6,sep = "-")) %>% 
  select(-c(c1,c2,c3,c4,c5,c6,c7)) %>% 
  merge.data.frame(ch4_best, by="UniqueID", all=T) %>% 
  merge.data.frame(n2o_best, by="UniqueID_notime", all=T) %>% 
  select(-UniqueID_notime)

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

#3.2. Curonian lagoon (CU): inconsistent bare representation (only consistent in restored sites, non-appropriate). We will only consider vegetated and open water plots
CU_wrongstrata_IDs<- fieldobs_uniqueID %>% filter(casepilot=="CU") %>% filter(strata=="bare") %>% pull(UniqueID)

#3.3. Danube delta (DA): inconsistent strata in some samplings: exclude vegetated of S3-DA-A2 (not representative), exclude bare of S1-DA-R2 (only sampled in first visit)
DA_wrongstrata_IDs<- fieldobs_uniqueID %>% filter(casepilot=="DA") %>% 
  filter((sampling=="S3-DA-A2"&strata=="vegetated")|(sampling=="S1-DA-R2"&strata=="bare")) %>% 
  pull(UniqueID)

#3.4. Dutch delta (DU): NOTHING TO REMOVE, all incubations are in appropriate strata
DU_wrongstrata_IDs<- c()

#3.5. Ria d'Aveiro (RI): rising-tide and aftervegcut already removed. Inconsistent incubations: tidal-pool incubations in S3-RI-A1, bare in RI-R1 and in RI-R2
RI_wrongstrata_IDs<- fieldobs_uniqueID %>% filter(casepilot=="RI") %>% 
  filter((sampling=="S3-RI-A1"&grepl("tidal-pool", field_observations))|(subsite=="RI-R1"&strata=="bare")|(subsite=="RI-R2"&strata=="bare")) %>% pull(UniqueID)

#3.6. Valencian wetland (VA): NOTHING TO REMOVE, all incubations are in appropriate strata
VA_wrongstrata_IDs<-c()


#ONLY valid incubation details
valid_fieldobs_uniqueID<- fieldobs_uniqueID %>% 
  filter(!UniqueID%in%c(tidal_IDs,vegcut_IDs,CA_wrongstrata_IDs,CU_wrongstrata_IDs,DA_wrongstrata_IDs,DU_wrongstrata_IDs,RI_wrongstrata_IDs,VA_wrongstrata_IDs))

rm(tidal_IDs,vegcut_IDs,CA_wrongstrata_IDs,CU_wrongstrata_IDs,DA_wrongstrata_IDs,DU_wrongstrata_IDs,RI_wrongstrata_IDs,VA_wrongstrata_IDs)

#4. Strata composition ----

###to-check: RI stratacomp----
#WEIRD patterns/doubts strata_composition: 
#RI preserved, some bare plots (very few), RI-P2: only 1 bare plot in first season. RI-P1: 4 bare plots, all during S3 sampling.



#Construct Dummy per-sampling strata composition: this will be replaced by RS-derived strata composition (when it becomes available)

#Dataset format: season, casepilot,subsite,sampling,strata, strata_prop(0-1)

#strata_prop is calculated for each sampling based on strata distribution of chambers non-excluded by design (those still in valid_plotcodes)
valid_plotcodes<- valid_fieldobs_uniqueID %>% 
  select(plotcode, season, casepilot,subsite,sampling, strata) %>% distinct()

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
  select(season, casepilot,subsite,sampling,strata,strata_prop) %>%
  ungroup()


#CU: inconsistent boundaries, no seasonal change. All subsites will get the same strata composition (calculated with all avaliable chambers, not grouped)
cu_stratacomp<- valid_plotcodes %>% filter(casepilot=="CU") %>% 
  #Obtain average strata_prop (from all available chambers, no grouping)
  group_by(casepilot,strata) %>%
  summarise(n_deployed=n()) %>% 
  ungroup() %>% 
  mutate(totalchambers=sum(n_deployed),
         strata_prop=n_deployed/totalchambers) %>% 
  select(casepilot,strata,strata_prop) %>%
  ungroup() %>% 
  #expand common stratacomp to all subiste samplings
  merge.data.frame(valid_plotcodes %>% filter(casepilot=="CU") %>% select(season,casepilot,subsite,sampling,strata) %>% distinct(), by=c("casepilot","strata"),all = T)%>% 
  select(season, casepilot,subsite,sampling,strata,strata_prop)


#DA: same strata composition for all preserved and restored subsites (P1,P2,R1,R2), using all chambers, not grouped. DA-A1 and DA-A2 have seasonal variability and distinct strata-composition, use chamber distribution directly. 

da_a_stratacomp<- valid_plotcodes %>% filter(subsite%in%c("DA-A1","DA-A2")) %>% 
  #get chambers per strata
  group_by(season,casepilot,subsite,sampling,strata) %>%
  summarise(n_deployed=n()) %>%
  #Calculate proportion of strata
  group_by(season,casepilot,subsite,sampling) %>% 
  mutate(totalchambers=sum(n_deployed),
         strata_prop=n_deployed/totalchambers) %>% 
  select(season, casepilot,subsite,sampling,strata,strata_prop) %>%
  ungroup()

da_pr_stratacomp<- valid_plotcodes %>% filter(subsite%in%c("DA-P1","DA-P2","DA-R1","DA-R2")) %>% 
  #Calculate general strata composition for all chamber deployments
  group_by(casepilot,strata) %>%
  summarise(n_deployed=n()) %>% 
  ungroup() %>% 
  mutate(totalchambers=sum(n_deployed),
         strata_prop=n_deployed/totalchambers) %>% 
  select(casepilot,strata,strata_prop) %>%
  ungroup() %>% 
  #expand common strata composition to each subsite
  merge.data.frame(valid_plotcodes %>% filter(subsite%in%c("DA-P1","DA-P2","DA-R1","DA-R2")) %>% select(season,casepilot,subsite,sampling,strata) %>% distinct(), by=c("casepilot","strata"),all = T)%>% 
  select(season, casepilot,subsite,sampling,strata,strata_prop)


#DU: use chamber distribution directly.
du_stratacomp<- valid_plotcodes %>% filter(casepilot=="DU") %>% 
  #get chambers per strata
  group_by(season,casepilot,subsite,sampling,strata) %>%
  summarise(n_deployed=n()) %>%
  #Calculate proportion of strata
  group_by(season,casepilot,subsite,sampling) %>% 
  mutate(totalchambers=sum(n_deployed),
         strata_prop=n_deployed/totalchambers) %>% 
  select(season, casepilot,subsite,sampling,strata,strata_prop) %>%
  ungroup()


#RI: use chamber distribution directly 
ri_stratacomp<- valid_plotcodes %>% filter(casepilot=="RI") %>% 
  #get chambers per strata
  group_by(season,casepilot,subsite,sampling,strata) %>%
  summarise(n_deployed=n()) %>%
  #Calculate proportion of strata
  group_by(season,casepilot,subsite,sampling) %>% 
  mutate(totalchambers=sum(n_deployed),
         strata_prop=n_deployed/totalchambers) %>% 
  select(season, casepilot,subsite,sampling,strata,strata_prop) %>%
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
  select(casepilot,subsite,strata,strata_prop) %>%
  ungroup() %>% 
  #expand common strata composition to each sampling
  merge.data.frame(valid_plotcodes %>% filter(subsite%in%c("VA-R1")) %>% select(season,casepilot,subsite,sampling,strata) %>% distinct(), by=c("casepilot","subsite","strata"),all = T)%>% 
  select(season, casepilot,subsite,sampling,strata,strata_prop)

#Rest of VA: 
va_notR1_stratacomp<- valid_plotcodes %>% filter(subsite%in%c("VA-A1","VA-A2","VA-P1","VA-P2","VA-R2")) %>% 
  #get chambers per strata
  group_by(season,casepilot,subsite,sampling,strata) %>%
  summarise(n_deployed=n()) %>%
  #Calculate proportion of strata
  group_by(season,casepilot,subsite,sampling) %>% 
  mutate(totalchambers=sum(n_deployed),
         strata_prop=n_deployed/totalchambers) %>% 
  select(season, casepilot,subsite,sampling,strata,strata_prop) %>%
  ungroup()


#Join and format stratacomp:
stratacomp<- rbind(ca_stratacomp,cu_stratacomp,da_a_stratacomp,da_pr_stratacomp,du_stratacomp,ri_stratacomp,va_notR1_stratacomp,va_r1_stratacomp) %>% 
  #make 0 strata proportion explicit (pivot_wider-->pivot_longer)
  pivot_wider(names_from = strata,values_from=strata_prop ,values_fill = 0) %>% 
  pivot_longer(cols = c(bare,vegetated,`open water`), names_to = "strata",values_to = "strata_prop")

rm(ca_stratacomp,cu_stratacomp,da_a_stratacomp,da_pr_stratacomp,du_stratacomp,ri_stratacomp,va_notR1_stratacomp,va_r1_stratacomp)




##Strata composition plots-----
strata_colors <- c(
  "bare" = "#d95f02",      # brown/orange
  "vegetated" = "#228B22", # green
  "open water" = "#1f78b4" # blue
)

##CA strata_prop plot
stratacomp %>% filter(casepilot == "CA") %>% 
  ggplot(aes(x = season, y = strata_prop*100, fill = strata)) +
  geom_bar(stat = "identity", col="black") +
  facet_wrap(~subsite) +
  scale_fill_manual(values = strata_colors, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format(scale = 1),expand=c(0,0)) +
  labs(y = "Strata Chambers (%)", x = "Season", fill = "Strata") +
  theme_minimal()

##CU strata_prop plot
stratacomp %>% filter(casepilot == "CU") %>% 
  ggplot(aes(x = season, y = strata_prop*100, fill = strata)) +
  geom_bar(stat = "identity", col="black") +
  facet_wrap(~subsite) +
  scale_fill_manual(values = strata_colors, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format(scale = 1),expand=c(0,0)) +
  labs(y = "Strata Chambers (%)", x = "Season", fill = "Strata") +
  theme_minimal()

##DA strata_prop plot
stratacomp %>% filter(casepilot == "DA") %>% 
  ggplot(aes(x = season, y = strata_prop*100, fill = strata)) +
  geom_bar(stat = "identity", col="black") +
  facet_wrap(~subsite) +
  scale_fill_manual(values = strata_colors, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format(scale = 1),expand=c(0,0)) +
  labs(y = "Strata Chambers (%)", x = "Season", fill = "Strata") +
  theme_minimal()

##DU strata_prop plot
stratacomp %>% filter(casepilot == "DU") %>% 
  ggplot(aes(x = season, y = strata_prop*100, fill = strata)) +
  geom_bar(stat = "identity", col="black") +
  facet_wrap(~subsite) +
  scale_fill_manual(values = strata_colors, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format(scale = 1),expand=c(0,0)) +
  labs(y = "Strata Chambers (%)", x = "Season", fill = "Strata") +
  theme_minimal()

##RI strata_prop plot
stratacomp %>% filter(casepilot == "RI") %>% 
  ggplot(aes(x = season, y = strata_prop*100, fill = strata)) +
  geom_bar(stat = "identity", col="black") +
  facet_wrap(~subsite) +
  scale_fill_manual(values = strata_colors, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format(scale = 1),expand=c(0,0)) +
  labs(y = "Strata Chambers (%)", x = "Season", fill = "Strata") +
  theme_minimal()


##VA strata_prop plot
stratacomp %>% filter(casepilot == "VA") %>% 
  ggplot(aes(x = season, y = strata_prop*100, fill = strata)) +
  geom_bar(stat = "identity", col="black") +
  facet_wrap(~subsite) +
  scale_fill_manual(values = strata_colors, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format(scale = 1),expand=c(0,0)) +
  labs(y = "Strata Chambers (%)", x = "Season", fill = "Strata") +
  theme_minimal()



#WEIRD patterns/doubts strata_composition: 
#RI preserved, some bare plots (very few), RI-P2: only 1 bare plot in first season. RI-P1: 4 bare plots, all during S3 sampling.



#5. Integrate Daily fluxes ------
#From ghg_best (uniqueID) to ghg_daybest (plotcode)


##5.1. Get daylight hours----

#daylight hours are defined as time from sunrise to sunset (using package suncalc, median lat/long of each subsite and date of sampling)

# get median coordinates per subsite (including all seasons)
samplings <- valid_fieldobs_uniqueID %>% 
  select(subsite, latitude, longitude) %>% 
  group_by(subsite) %>% 
  summarise(subsite_latitude=median(latitude,na.rm=T),
            subsite_longitude=median(longitude,na.rm=T)) %>% 
  ungroup()

#Use dplyr and suncalc to calculate daylight hours for each subsite-date combination (each sampling).
sampling_daylight <- samplings %>%
  merge.data.frame(y=valid_fieldobs_uniqueID %>% select(subsite, date,sampling), by="subsite", all=T) %>% 
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
  select(sampling, hours_day, hours_night)




## 5.2. Scale to daily flux-----

#Daily integration will follow the previously agreed decisions (common for all gas species): 
  #Vegetated fluxes will be systematically scaled with daylight hours (T+D)
  #Open water fluxes will be directly used (only D)
  #Bare fluxes will be non-systematically scaled according to data availability: 
      #Dark flux only--> dark flux = daily flux
      #Transparent only --> transparent flux = daily flux
      #T+D available--> scaling with daylight hours (T+D) = daily flux

#NOTE: Inconsistency in scaling of bare fluxes responds to exploration of data, where light condition has seen not as important as plot-to-plot variability and wet/dry sediment, and where discarding transparent bare fluxes will severely limit the representativity of this stratum. 

#start point: vaild_fieldobs  to retrieve fluxes, then 

#get field details of valid plotcodes
valid_plotcode_obs<- fieldobs_plotcode %>%
  filter(plotcode%in%valid_plotcodes$plotcode) 

#Get ghg fluxes for transparent & dark by plotcode (only for valid UniqueIDs)
valid_ghg_transpdark<- valid_fieldobs_uniqueID %>%
  merge.data.frame(ghg_best, by="UniqueID", all.x = T) %>% 
  select(plotcode, co2_best, ch4_best, n2o_best, transparent_dark) %>% 
  rename(co2=co2_best, ch4=ch4_best, n2o=n2o_best) %>% 
  pivot_longer(cols=c(co2,ch4,n2o), names_to = "ghg", values_to = "flux") %>% 
  pivot_wider(names_from = transparent_dark, values_from = flux)


#Integrate to dailyflux
daily_ghg_plotcode<- valid_plotcode_obs %>% 
  #Add daylight hours
  merge.data.frame(sampling_daylight, by=c("sampling"), all = T) %>% 
  select(plotcode, season, casepilot, subsite, sampling, plot_num, date, plot_start_time, latitude, longitude, gas_analyzer, strata, water_depth, field_observations, ABG_veg_description, ABG_veg_biomass_gpersquaremeter, hours_day, hours_night) %>% 
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
select(-c(dark, transparent)) %>% 
  pivot_wider(names_from = ghg, values_from = dailyflux ) 
  #POR AQUI------
#UNITS AND GLOBAL Warming potential 
#Units for daily flux (up to now): 
  #CO2: umol per m2 per day
  #CH4: nmol per m2 per day
  #N2O: nmol per m2 per day

  mutate(c_co2eq=)



#Delete objects already used
rm(samplings, sampling_daylight, fieldobs_plotcode, fieldobs_uniqueID, valid_plotcode_obs, valid_plotcodes)





#POR AQUI------

#DOUBT: propagate error as well to daily? We cannot use the error in the models


