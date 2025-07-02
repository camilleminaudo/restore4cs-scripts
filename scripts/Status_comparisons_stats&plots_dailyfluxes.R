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
#1. Template model and decisions: data-transformations, outlier idetification, strata-scaling implementation


#PLots and LMM structure (general and per-case pilot): following decissions  


#Inputs: 
#Per-UniqueID: co2, ch4 and n2o best fluxes (and se flux)
#Per-UniqueID and per-plotcode harmonized field observations (created/Updated by  script Harmonize_Fieldsheet-veg_perPLOT_and_perUniqueID). 
#Strata-representativity decissions: common and relevant strata for comparisons


rm(list = ls()) # clear workspace
cat("/014") # clear console


# ---- Packages ----
#For data-handling & plotting
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
#For modelling: 
library(lmerTest) #similar to lme4 but allows to estimate significance of predictors (and interactions)
library(bestNormalize)
library(emmeans)
library(DHARMa)
library(performance)
library(ggResidpanel)
library(tibble)
library(partR2)#for breaking down variance by fixed effect
library(car)  # for leveneTest
library(caret)  # for cross-validation
library(MuMIn)  # for AICc if needed

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
  pivot_wider(names_from = ghg, values_from = dailyflux) %>% 
  #Transform molar units for daily flux: 
  #CO2: umol per m2 per day --> mol per m2 per day
  #CH4: nmol per m2 per day --> mmol per m2 per day
  #N2O: nmol per m2 per day --> mmol per m2 per day
  mutate(co2_mol=co2*1e-6,
         ch4_mmol=ch4*1e-6,
         n2o_mmol=n2o*1e-6,
         #Global warming potential GWP100: 
         gwp_co2=co2_mol*44, #molar flux* molar mass-->g/m2/day
         gwp_ch4=(ch4_mmol*1e-3)*16*27,#molar flux * molar mass* GWP100 factor 27
         gwp_n2o=(n2o_mmol*1e-3)*44*273, #molar flux * molar mass 44* GWP100factor 273
         gwp_co2andch4= gwp_co2+gwp_ch4, #gwp of co2 and ch4 only
         gwp_allghg=gwp_co2+gwp_ch4+gwp_n2o) %>%  #gwp of all GHGs
  select(-c(co2, ch4, n2o)) %>% 
  #Add modeling variables: 
  mutate(status=substr(x = sampling, start = 7, stop=7),
         status=case_when(status=="A"~"Altered",
                          status=="P"~"Preserved",
                          status=="R"~"Restored"))

#Delete objects already used
rm(samplings, sampling_daylight, fieldobs_plotcode, fieldobs_uniqueID, valid_plotcode_obs, valid_plotcodes)



#6. Final format -----
  
#Produce final dataset to be used in models. Implementing strata-weigths and formating to simplify filtering for modeling. 

#DOUBT: propagate error as well to daily? We cannot use the error in the models
#DOUBT: how to use instantaneous fluxes if we do not see effect in daily? (not critical)

#We will be using always gwp units (gCo2eq /m2/d), this way we have a more straightforward interpretation and we avoid multiple units. 

#We format the dataset to be easily used across our modeling approaches: 
#categories are taken as factors, we pivot longer the daily fluxes (response variable), we add the calculation of sample_weight based on strata proportion and sample proportion (non-NA) across strata (corrects sampling bias and ensures equal influence across samplings in the models)

data4models <- daily_ghg_plotcode %>% 
  select(plotcode, season, casepilot, status, subsite, sampling, strata, 
         gwp_co2, gwp_ch4, gwp_n2o, gwp_co2andch4, gwp_allghg) %>% 
  merge.data.frame(stratacomp, by = c("season","casepilot","subsite","sampling","strata"), all.x = TRUE) %>% 
  pivot_longer(
    cols = c(gwp_co2, gwp_ch4, gwp_n2o, gwp_co2andch4, gwp_allghg),
    names_to = "ghgspecies",
    values_to = "dailyflux"
  ) %>% 
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
  
  # Normalize weights only over non-NA rows
  group_by(sampling, ghgspecies) %>%
  mutate(sample_weight = ifelse(
    is.na(dailyflux), NA, 
    sample_weight_raw / sum(sample_weight_raw[!is.na(dailyflux)])
  )) %>%
  ungroup() %>%
  # Drop intermediate vars
  select(-c(sample_prop, nonNA_n_persampling, nonNA_n_perstrata, sample_weight_raw, strata_prop))


#Check: do all samplings with valid (non-na observations) share the same overall sum of sample_weights?
data4models %>% 
  filter(!is.na(dailyflux)) %>% 
  group_by(sampling, ghgspecies) %>% 
  summarise(sumweights=sum(sample_weight),.groups = "drop") %>% ungroup() %>%  pull(sumweights) %>% round(., digits = 10)%>% unique()

#Yes, always same overall weight for each sampling




#__________________------
#MODELS------


#Decisions: how to "log-transform" data? we have cero and negative values. 
#RUN each model with various transformations: 
  #DIRECT (not transformed):dailyflux
  #sign-log transformation: log_dailyflux  
  #sign-srqt transformation: root_dailyflux
  #Yeo Johnson: Yeo Johnson Transformation depends on the distribution of the input data, so I must perform it directly on each subset of dataset to be used in a model (cannot do it directly in the whole dataset). Also, Box-cox  cannot be used (doesnt work with negative values). ISSUE: Yeo Johnson sometimes swaps sign of very low-magnitude values. IT CANNOT BE USED
  



#A. FULL models-----
#one per ghgspecies, using data from all casepilots. 

#A.1. CO2 ------
## ========== 0. subset & trans ==========
data_sub <- data4models %>% 
  filter(ghgspecies == "gwp_co2", !is.na(dailyflux)) %>% 
  mutate(
    # Signed log transform (retain sign):retaining the original sign and adding +1 to avoid sign change and cero values. To back-transform original data (for interpretation of model outputs): dailyflux ≈ sign(log_dailyflux) * (exp(abs(log_dailyflux)) - 1)
    dailyflux_logsign = sign(dailyflux) * log(abs(dailyflux) + 1),
    
    # Signed square root transform (retain sign): backtransformation: dailyflux ≈ sign(dailyflux_rootsign)* (dailyflux_rootsign^2) 
    dailyflux_rootsign = sign(dailyflux) * sqrt(abs(dailyflux))
  )

#Define back-transformation functions: 
# For signed log transformation
inv_signedlog <- function(x) sign(x) * (exp(abs(x)) - 1)

# For signed root transformation
inv_rootsign <- function(x) sign(x) * (x^2)




## ========== 1. Untransformed Model ==========
model_untrans <- lmer(dailyflux ~ status * season * casepilot + (1 | subsite), 
                      data = data_sub, weights = sample_weight)

## ========== 2. Signed Log Transform Model ==========
model_logsign <- lmer(dailyflux_logsign ~ status * season * casepilot + (1 | subsite),
                        data = data_sub, weights = sample_weight)

## ========== 3. Signed Root Transform Model ==========
model_rootsign <- lmer(dailyflux_rootsign ~ status * season * casepilot + (1 | subsite),
                 data = data_sub, weights = sample_weight)

## ========== 4. Evaluate Model requirements ==========
check_model_diagnostics <- function(model, model_name, data) {

  
  cat("----- Diagnostics for", model_name, "-----\n")
  
  # 1. Residual diagnostics panel (Q-Q, residuals vs fitted, histogram, etc.)
  print(ggResidpanel::resid_panel(model, qqbands = TRUE, plots = c("qq", "resid", "hist", "fitted")))
  
  # 1b. Residuals vs key predictors (check for heteroscedasticity)
  res <- residuals(model)
  fitted_vals <- fitted(model)
  
  # Residuals vs status
  p1 <- ggplot(data.frame(status = data$status, resid = res), aes(x = status, y = resid)) +
    geom_boxplot() +
    labs(title = paste("Residuals by Status:", model_name), y = "Residuals")
  print(p1)
  
  # Residuals vs casepilot
  p2 <- ggplot(data.frame(casepilot = data$casepilot, resid = res), aes(x = casepilot, y = resid)) +
    geom_boxplot() +
    labs(title = paste("Residuals by Casepilot:", model_name), y = "Residuals")
  print(p2)
  
  # Residuals vs season
  p3 <- ggplot(data.frame(season = data$season, resid = res), aes(x = season, y = resid)) +
    geom_boxplot() +
    labs(title = paste("Residuals by Season:", model_name), y = "Residuals")
  print(p3)
  
  # 1c. Levene's Test for homogeneity of variance
  cat("Levene's Test for residual variance by status:\n")
  print(leveneTest(res ~ data$status))
  
  cat("Levene's Test for residual variance by casepilot:\n")
  print(leveneTest(res ~ data$casepilot))
  
  cat("Levene's Test for residual variance by season:\n")
  print(leveneTest(res ~ data$season))
  
  # 2. Normality check: Shapiro-Wilk test
  cat("Shapiro-Wilk test for residual normality:\n")
  print(shapiro.test(res))
  
  # 3. R² values
  r2 <- performance::r2_nakagawa(model)
  
  # 4. AIC and BIC
  aic_val <- AIC(model)
  bic_val <- BIC(model)
  
  # 5. Predictive performance (10-fold CV RMSE)
  ctrl <- trainControl(method = "cv", number = 10)
  cv_model <- train(
    x = model.matrix(model)[, -1],
    y = getME(model, "y"),
    method = "lm",
    trControl = ctrl
  )
  rmse_cv <- cv_model$results$RMSE
  
  return(list(
    r2 = r2,
    AIC = aic_val,
    BIC = bic_val,
    RMSE_CV = rmse_cv
  ))
}


##__4.1. Visual inspection----
# Run diagnostics for each model.
#Look for straight line at 45º for QQplot, normal-looking residual histogram, No clear pattern or funnel shape in Residual plot
diag_untrans <- check_model_diagnostics(model_untrans, "Untransformed", data_sub)
#Untransformed Observations: clear non-normality of residuals, heterocedasticity at season and casepilot level (not at status level)
diag_logsign <- check_model_diagnostics(model_logsign, "Signed Log", data_sub)
#Logtransformed Observations: Residuals show reasonable normaility. Heterocedasticity at season and casepilot level (not at status level), OK heterocedasticity is expected. BEST
diag_rootsign <- check_model_diagnostics(model_rootsign, "Signed Root", data_sub)
#Roottransformed Observations: clear non-normality of residuals, Heterocedasticity at season and casepilot level (not at status level), OK heterocedasticity is expected. 


##__4.2. Quantitative comparison----
# Compile results 
results <- tibble(
  model = c("Untransformed", "Signed Log", "Signed Root"),
  marginal_r2 = c(diag_untrans$r2$R2_marginal, diag_logsign$r2$R2_marginal, diag_rootsign$r2$R2_marginal),
  conditional_r2 = c(diag_untrans$r2$R2_conditional, diag_logsign$r2$R2_conditional, diag_rootsign$r2$R2_conditional),
  AIC = c(diag_untrans$AIC, diag_logsign$AIC, diag_rootsign$AIC),
  BIC = c(diag_untrans$BIC, diag_logsign$BIC, diag_rootsign$BIC),
  RMSE_CV = c(diag_untrans$RMSE_CV, diag_logsign$RMSE_CV, diag_rootsign$RMSE_CV)
)

# Rank models by each metric (lower AIC, BIC, RMSE are better; higher R² is better)
results <- results %>%
  mutate(
    rank_AIC = rank(AIC),
    rank_BIC = rank(BIC),
    rank_RMSE = rank(RMSE_CV),
    rank_marginal_r2 = rank(-marginal_r2),
    rank_conditional_r2 = rank(-conditional_r2),
    total_score = rank_AIC + rank_BIC + rank_RMSE + rank_marginal_r2 + rank_conditional_r2
  ) %>%
  arrange(total_score)

print(results)

# Identify best model
best_model_name <- results$model[1]
cat("Best model based on combined criteria:", best_model_name, "\n")


## ========== 5. Interpret / Visualize Effects with Best Model ==========

#COPY best model selection: 
best_model<- model_logsign

#Set back-transformation function accordingly: inv_signedlog or inv_rootsign
inv_function<- inv_signedlog

## ---- 5.1. Overall summary ----
summary(best_model)

# This gives you:
# Random effect variance
# Residual variance
# Fixed effect estimates (also for interactions)
# p-values for each main effect and interaction (Satterthwaite's method by default, to change the method you can use: lmerTest::lmerControl(calc.derivs = FALSE, optimizer = "bobyqa")



## ---- 5.2. Explained Variance ----

performance::r2_nakagawa(best_model)
# Marginal R²: variance explained by fixed effects
# Conditional R²: variance explained by fixed + random effects
# 1 - Condinional R²: unnexplained variance


#VARIANCE decomposition (%) (REDUNDANT APPROACH)

#What every component means: 
#var.fixed= explained by fixed effects
#var.random= explained by random effects (duplicated as var.intercept subsite)
#var.residual= unexplained residual variance

# Extract variance components
var_components <- get_variance(model_logsign)
# Drop duplicates: get_variance produces duplicate variances (same var, 2 names assigned)
var_clean <- var_components[c("var.fixed", "var.random", "var.residual")]
# Calculate percentage of total (explained by fixed, explained by random, unnexplained)
round(100 * unlist(var_clean) / sum(unlist(var_clean)), 2)
#Same results as calcualted from maringal and conditional R^2 (REDUNDANT info)



#VARIANCE contribution of individual fixed effects (DOES NOT WORK, Model too complex)
#CHECK---- this requires the model to be calculated with the lme4 package
model_lme4 <- lme4::lmer(dailyflux_logsign ~ status * season * casepilot + (1 | subsite), 
                   data = data_sub, weights = sample_weight)

#TAKES A LONG TIME TO CALCULATE: 20mins with 1000 nboot, issues: "fixed-effect model matrix is rank deficient so dropping 1 column / coefficient" 
partR2(model_lme4,
       partvars = c("status", "season", "casepilot"),
       R2_type = "marginal",
       nboot = 10)
#It fails because: 
  #model too complex (too many interactions, leading to very few points per uniquecombination.)
  #small nboot--> Confident interval fail




## ---- 5.3. EMMs  ----
#Estimated Marginal Means (emms)

# EMMs for status by season and casepilot
emm_all <- emmeans(best_model, ~ status * season * casepilot)
summary(emm_all)
emm_df <- as.data.frame(emm_all)

#Back-transform EMMs
emm_df_bt<- emm_df %>% 
  mutate(across(c(emmean, SE,lower.CL,upper.CL), inv_function))
  
#Plot EMMs 
ggplot(emm_df_bt, aes(x = season, y = emmean, fill = status)) +
  geom_col(position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE),
                position = position_dodge(width = 0.8), width = 0.2) +
  facet_wrap(~ casepilot,scales="free") +
  labs(
    title = "Estimated Marginal Means by Status, Season, and Casepilot",
    y = expression(Back-transformed~EMM~(g~CO2[eq]~d^-1~m^-2)),
    x = "Season",
    col = "Status"
  ) +
  theme_minimal()

## ---- 5.4. Post-hoc on EMMs  ----
# EMM comparison for status by season and casepilot comparisons
pairs(emm_all, by = c("season", "casepilot"))


###POR AQUI---------
#CHECK data for apparent outliers: seems DA S1 altered emmeans are too high in comparison to rest (check magnitude of fluxes and sample_weight's ). REVIEW Jorge's script and use the same approach for detecting outliers (doubt: using transformed or untransformed data?)

#For tomorrow: provide to co-pilot info on database, model structure, main outputs and comparisons with model_no3way. Explore patterns and reach a decission for model complexity. 

#Why status and season are not significant but casepilot and interactions appear strongly significant  (would this serve as justification for casepilot-specific models)
#3-way interaction appears highly significant, 

#Model too complex? 
#check number of observations when broken down by status*season*casepilot (and aditionally subsite)

data4models %>% 
  filter(ghgspecies=="gwp_co2") %>% 
  group_by(status,season, casepilot, subsite) %>% 
  summarise(nobs=sum(!is.na(dailyflux))) %>% 
  ungroup() %>% 
  summarise(min_n=min(nobs),
            max_n=max(nobs),
            median_n=median(nobs),
            mean_n=mean(nobs),
            sd_n=sd(nobs))

data4models %>% 
  filter(ghgspecies=="gwp_co2") %>% 
  group_by(status,season, casepilot) %>% 
  summarise(nobs=sum(!is.na(dailyflux))) %>% 
  ungroup() %>% 
  summarise(min_n=min(nobs),
            max_n=max(nobs),
            median_n=median(nobs),
            mean_n=mean(nobs),
            sd_n=sd(nobs))

#OUTPUT FROM CO-PILOT: Is the Model Too Complex?
#   ✅ Arguments in Favor of Your Model
# You expect interactions: status effects may differ by season and casepilot.
# You have enough observations per group (≥17) to estimate effects, though with caution.
# You're using weights to correct for sampling imbalance, which helps.
# ⚠️ Potential Concerns
# 3-way interactions are hard to interpret and may overfit if not strongly supported by the data.
# If many interaction terms are non-significant, they may be redundant.
# Multicollinearity or sparse cells can inflate standard errors and reduce power.

# How to Evaluate Model Complexity
# 1. Check significance of interaction terms:
car::Anova(best_model, type = 3) # Type III ANOVA (if using lmerTest)

# 2. Compare nested models
# Fit simpler models and compare using AIC/BIC:
model_no3way <- lmer(dailyflux_logsign ~ status + season + casepilot + status:season + status:casepilot + season:casepilot + (1 | subsite),
                     data=data_sub,
                     weights=sample_weight)
anova(model_no3way, best_model)

#Current more complex model significantly improves the fit, how much added explanatory power? is it worth it?


#Recommendation:
# Start with the full model (as you’ve done).
# Test whether the 3-way interaction is significant and improves model fit.
# If not, drop it and retain only 2-way interactions.
# If even those are weak, consider additive models or casepilot-specific models.






#Exploratory plots-----

#GWP_CO2
daily_ghg_plotcode %>%
  mutate(subsite2=substr(subsite,4,5)) %>% 
  group_by(casepilot, status, subsite,subsite2) %>%
  filter(!is.na(gwp_co2)) %>% 
  ggplot(aes(x=subsite2, y=sign(gwp_co2)*log(abs(gwp_co2)+1), fill=status, group=paste0(subsite)))+
  geom_boxplot()+
  facet_wrap(~casepilot,scales="free")



#boxplots using per-plot daily fluxes 
#CA: 

daily_ghg_plotcode %>%
  filter(casepilot=="CA") %>% 
  filter(!is.na(co2_mol)) %>% 
  ggplot(aes(x=subsite,y=co2_mol, fill= status, group=sampling))+
  geom_boxplot()

daily_ghg_plotcode %>%
  filter(casepilot=="CA") %>% 
  filter(!is.na(co2_mol)) %>% 
  ggplot(aes(x=co2_mol))+geom_histogram()


daily_ghg_plotcode %>%
  filter(casepilot=="CA") %>% 
  filter(!is.na(ch4_mmol)) %>% 
  ggplot(aes(x=subsite,y=sign(ch4_mmol)*log10(abs(ch4_mmol)), fill= status, group=sampling))+
  geom_boxplot()

daily_ghg_plotcode %>%
  filter(casepilot=="CA") %>% 
  filter(!is.na(n2o_mmol)) %>% 
  ggplot(aes(x=subsite,y=n2o_mmol, fill= status, group=sampling))+
  geom_boxplot()

daily_ghg_plotcode %>%
  filter(casepilot=="CA") %>% 
  filter(!is.na(c_co2eq)) %>% 
  ggplot(aes(x=subsite,y=c_co2eq, fill= status, group=sampling))+
  geom_boxplot()

daily_ghg_plotcode %>%
  filter(casepilot=="CA") %>% 
  filter(!is.na(ghg_co2eq)) %>% 
  ggplot(aes(x=subsite,y=ghg_co2eq, fill= status, group=sampling))+
  geom_boxplot()

daily_ghg_plotcode %>%
  mutate(subsite2=substr(subsite,4,5)) %>% 
  filter(!is.na(gwp_ch4*gwp_co2)) %>% 
  pivot_longer(cols = c(gwp_ch4,  gwp_co2),names_to = "ghg",values_to = "gwp") %>% 
  ggplot(aes(x=subsite2, y=gwp, fill=ghg, group=paste0(ghg,subsite)))+
  geom_boxplot(outliers = F)+
  facet_wrap(~casepilot,scales="free")



#CU: 

daily_ghg_plotcode %>%
  filter(casepilot=="CU") %>% 
  filter(!is.na(co2_mol)) %>% 
  ggplot(aes(x=subsite,y=co2_mol, fill= status, group=sampling))+
  geom_boxplot()

daily_ghg_plotcode %>%
  filter(casepilot=="CU") %>% 
  filter(!is.na(ch4_mmol)) %>% 
  ggplot(aes(x=subsite,y=ch4_mmol, fill= status, group=sampling))+
  geom_boxplot()+
  scale_y_log10()

daily_ghg_plotcode %>%
  filter(casepilot=="CU") %>% 
  filter(!is.na(n2o_mmol)) %>% 
  ggplot(aes(x=subsite,y=n2o_mmol, fill= status, group=sampling))+
  geom_boxplot()

daily_ghg_plotcode %>%
  filter(casepilot=="CU") %>% 
  filter(!is.na(c_co2eq)) %>% 
  ggplot(aes(x=subsite,y=c_co2eq, fill= status, group=sampling))+
  geom_boxplot()

daily_ghg_plotcode %>%
  filter(casepilot=="CU") %>% 
  filter(!is.na(ghg_co2eq)) %>% 
  ggplot(aes(x=subsite,y=ghg_co2eq, fill= status, group=sampling))+
  geom_boxplot()


#DA: 

daily_ghg_plotcode %>%
  filter(casepilot=="DA") %>% 
  filter(!is.na(co2_mol)) %>% 
  ggplot(aes(x=subsite,y=co2_mol, fill= status, group=sampling))+
  geom_boxplot()

daily_ghg_plotcode %>%
  filter(casepilot=="DA") %>% 
  filter(!is.na(ch4_mmol)) %>% 
  ggplot(aes(x=subsite,y=ch4_mmol, fill= status, group=sampling))+
  geom_boxplot()+
  scale_y_log10()

daily_ghg_plotcode %>%
  filter(casepilot=="DA") %>% 
  filter(!is.na(n2o_mmol)) %>% 
  ggplot(aes(x=subsite,y=n2o_mmol, fill= status, group=sampling))+
  geom_boxplot()

daily_ghg_plotcode %>%
  filter(casepilot=="DA") %>% 
  filter(!is.na(c_co2eq)) %>% 
  ggplot(aes(x=subsite,y=c_co2eq, fill= status, group=sampling))+
  geom_boxplot()

daily_ghg_plotcode %>%
  filter(casepilot=="DA") %>% 
  filter(!is.na(ghg_co2eq)) %>% 
  ggplot(aes(x=subsite,y=ghg_co2eq, fill= status, group=sampling))+
  geom_boxplot()



#RI: 

daily_ghg_plotcode %>%
  filter(casepilot=="RI") %>% 
  filter(!is.na(co2_mol)) %>% 
  ggplot(aes(x=subsite,y=co2_mol, fill= status, group=sampling))+
  geom_boxplot()

daily_ghg_plotcode %>%
  filter(casepilot=="RI") %>% 
  filter(!is.na(ch4_mmol)) %>% 
  ggplot(aes(x=subsite,y=ch4_mmol, fill= status, group=sampling))+
  geom_boxplot()+
  scale_y_log10()

daily_ghg_plotcode %>%
  filter(casepilot=="RI") %>% 
  filter(!is.na(n2o_mmol)) %>% 
  ggplot(aes(x=subsite,y=n2o_mmol, fill= status, group=sampling))+
  geom_boxplot()

daily_ghg_plotcode %>%
  filter(casepilot=="RI") %>% 
  filter(!is.na(c_co2eq)) %>% 
  ggplot(aes(x=subsite,y=c_co2eq, fill= status, group=sampling))+
  geom_boxplot()

daily_ghg_plotcode %>%
  filter(casepilot=="RI") %>% 
  filter(!is.na(ghg_co2eq)) %>% 
  ggplot(aes(x=subsite,y=ghg_co2eq, fill= status, group=sampling))+
  geom_boxplot()


#VA: 

daily_ghg_plotcode %>%
  filter(casepilot=="VA") %>% 
  filter(!is.na(co2_mol)) %>% 
  ggplot(aes(x=subsite,y=co2_mol, fill= status, group=sampling))+
  geom_boxplot()

daily_ghg_plotcode %>%
  filter(casepilot=="VA") %>% 
  filter(!is.na(ch4_mmol)) %>% 
  ggplot(aes(x=subsite,y=ch4_mmol, fill= status, group=sampling))+
  geom_boxplot()+
  scale_y_log10()

daily_ghg_plotcode %>%
  filter(casepilot=="VA") %>% 
  filter(!is.na(n2o_mmol)) %>% 
  ggplot(aes(x=subsite,y=n2o_mmol, fill= status, group=sampling))+
  geom_boxplot()

daily_ghg_plotcode %>%
  filter(casepilot=="VA") %>% 
  filter(!is.na(c_co2eq)) %>% 
  ggplot(aes(x=subsite,y=c_co2eq, fill= status, group=sampling))+
  geom_boxplot()

daily_ghg_plotcode %>%
  filter(casepilot=="VA") %>% 
  filter(!is.na(ghg_co2eq)) %>% 
  ggplot(aes(x=subsite,y=ghg_co2eq, fill= status, group=sampling))+
  geom_boxplot()



