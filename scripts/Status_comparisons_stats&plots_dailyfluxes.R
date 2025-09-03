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
library(outliers) #Outlier exploration
library(rstatix) #Outlier exploration
library(emmeans)
library(glmmTMB)
library(DHARMa)
library(performance)
library(ggResidpanel)
library(tibble)
library(partR2)#for breaking down variance by fixed effect
library(car)  # for leveneTest
library(caret)  # for cross-validation
library(MuMIn)  # for AICc if needed
library(multcomp)# for emmeans post-hoc comparisons


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
plots_path <- paste0(dropbox_root,"/GHG/Processed data/computed_flux/Stats and plots SEFS14")

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


#RI: use chamber distribution directly 
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




##Strata composition plots-----
strata_colors <- c(
  "bare" = "#d95f02",      # brown/orange
  "vegetated" = "#228B22", # green
  "open water" = "#1f78b4" # blue
)

##CA strata_prop plot
strataplot_ca<- stratacomp %>% filter(casepilot == "CA") %>% 
  ggplot(aes(x = season, y = strata_prop*100, fill = strata)) +
  geom_bar(stat = "identity", col="black") +
  facet_wrap(~subsite) +
  scale_fill_manual(values = strata_colors, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format(scale = 1),expand=c(0,0)) +
  labs(y = "Strata Chambers (%)", x = "Season", fill = "Strata") +
  theme_minimal(base_line_size = 0.5)

ggsave(plot = strataplot_ca,filename = "Strata_comp_CA.png",
       width = 115,height = 100,units = "mm",device = "png",dpi = 400,
       path = plots_path)

##CU strata_prop plot
strataplot_cu<- stratacomp %>% filter(casepilot == "CU") %>% 
  ggplot(aes(x = season, y = strata_prop*100, fill = strata)) +
  geom_bar(stat = "identity", col="black") +
  facet_wrap(~subsite) +
  scale_fill_manual(values = strata_colors, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format(scale = 1),expand=c(0,0)) +
  labs(y = "Strata Chambers (%)", x = "Season", fill = "Strata") +
  theme_minimal()

ggsave(plot = strataplot_cu,filename = "Strata_comp_CU.png",
       width = 115,height = 100,units = "mm",device = "png",dpi = 400,
       path = plots_path)

##DA strata_prop plot
strataplot_da<- stratacomp %>% filter(casepilot == "DA") %>% 
  ggplot(aes(x = season, y = strata_prop*100, fill = strata)) +
  geom_bar(stat = "identity", col="black") +
  facet_wrap(~subsite) +
  scale_fill_manual(values = strata_colors, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format(scale = 1),expand=c(0,0)) +
  labs(y = "Strata Chambers (%)", x = "Season", fill = "Strata") +
  theme_minimal()

ggsave(plot = strataplot_da,filename = "Strata_comp_DA.png",
       width = 115,height = 100,units = "mm",device = "png",dpi = 400,
       path = plots_path)

##DU strata_prop plot
strataplot_du<- stratacomp %>% filter(casepilot == "DU") %>% 
  ggplot(aes(x = season, y = strata_prop*100, fill = strata)) +
  geom_bar(stat = "identity", col="black") +
  facet_wrap(~subsite) +
  scale_fill_manual(values = strata_colors, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format(scale = 1),expand=c(0,0)) +
  labs(y = "Strata Chambers (%)", x = "Season", fill = "Strata") +
  theme_minimal()

ggsave(plot = strataplot_du,filename = "Strata_comp_DU.png",
       width = 115,height = 100,units = "mm",device = "png",dpi = 400,
       path = plots_path)


##RI strata_prop plot
strataplot_ri<- stratacomp %>% filter(casepilot == "RI") %>% 
  ggplot(aes(x = season, y = strata_prop*100, fill = strata)) +
  geom_bar(stat = "identity", col="black") +
  facet_wrap(~subsite) +
  scale_fill_manual(values = strata_colors, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format(scale = 1),expand=c(0,0)) +
  labs(y = "Strata Chambers (%)", x = "Season", fill = "Strata") +
  theme_minimal()

ggsave(plot = strataplot_ri,filename = "Strata_comp_RI.png",
       width = 115,height = 100,units = "mm",device = "png",dpi = 400,
       path = plots_path)


##VA strata_prop plot
strataplot_va<- stratacomp %>% filter(casepilot == "VA") %>% 
  ggplot(aes(x = season, y = strata_prop*100, fill = strata)) +
  geom_bar(stat = "identity", col="black") +
  facet_wrap(~subsite) +
  scale_fill_manual(values = strata_colors, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format(scale = 1),expand=c(0,0)) +
  labs(y = "Strata Chambers (%)", x = "Season", fill = "Strata") +
  theme_minimal()

ggsave(plot = strataplot_va,filename = "Strata_comp_VA.png",
       width = 115,height = 100,units = "mm",device = "png",dpi = 400,
       path = plots_path)



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

#start point: vaild_fieldobs  to retrieve fluxes, then 

#get field details of valid plotcodes
valid_plotcode_obs<- fieldobs_plotcode %>%
  filter(plotcode%in%valid_plotcodes$plotcode) 

#Get ghg fluxes for transparent & dark by plotcode (only for valid UniqueIDs)
valid_ghg_transpdark<- valid_fieldobs_uniqueID %>%
  merge.data.frame(ghg_best, by="UniqueID", all.x = T) %>% 
  dplyr::select(plotcode, co2_best, ch4_best, n2o_best, transparent_dark) %>% 
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

#DOUBT: propagate error as well to daily? We cannot use the error in the models
#DOUBT: how to use instantaneous fluxes if we do not see effect in daily? (not critical)

#We will be using always gwp units (gCo2eq /m2/d), this way we have a more straightforward interpretation and we avoid multiple units. 

#We format the dataset to be easily used across our modeling approaches: 
#categories are taken as factors, we pivot longer the daily fluxes (response variable), we add the calculation of sample_weight based on strata proportion and sample proportion (non-NA) across strata (corrects sampling bias and ensures equal influence across samplings in the models)

data4models <- daily_ghg_plotcode %>% 
  dplyr::select(plotcode, season, casepilot, status, subsite, sampling, strata, 
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
  # Rescale weights to sum to N (non-NA count) per group: i.e. obtain 
  group_by(sampling, ghgspecies) %>%
  mutate(sample_weight = ifelse(
    is.na(dailyflux), NA,
    sample_weight_raw * sum(!is.na(dailyflux)) / sum(sample_weight_raw[!is.na(dailyflux)])
  )) %>%
  ungroup() %>% 
  # Drop intermediate vars
  dplyr::select(-c(sample_prop, nonNA_n_persampling, nonNA_n_perstrata, sample_weight_raw, strata_prop))


  #With the above weighting, we are accounting for sampling bias across strata at each sampling, but we are not equalizing the overall influence of each sampling. This is the most appropriate approach, given that glmmTMB treats weights as precission weights (all weights assumed as 1 when not specified otherwise), so the total sum of weights has to equal N within each sampling. If we try to equalize the influence of each sampling (i.e. overriding the potential effect of imbalanced N across samplings), the weights do not work in the intended way and the output from the models is not comparable (with and without sample weights). However, AIC is still not perfectly comparable between weighted and unweighted models, we should complement AIC with Predictive performance, Parameter stability, Residual diagnostics and Marignal&Conditional R2. 
  
#Check: do all samplings with valid (non-na) observations have a sum of sample_weights that equals their n?
data4models %>% 
  filter(!is.na(dailyflux)) %>% 
  group_by(sampling, ghgspecies) %>% 
  summarise(sumweights=sum(sample_weight), nobs=n(),.groups = "drop") %>% 
  filter(round(sumweights,digits = 1)!=nobs)#Filter for samplings*ghg combos without weight!=nobs
  
#OK, all samplings have weights thta match N. 

#7.Save data images--------

save.image(file = paste0(plots_path,"/Data4models.RData"))

load(paste0(plots_path, "/Data4models.RData"))
#__________________------
#MODELS------
#Good tutorial, focus on continuous predictors + random effects. 
# https://ourcodingclub.github.io/tutorials/mixed-models/

#IMPORTANT! -----
#Evaluate model structure and assumptions (residuals), not significance of effects nor even explained variability (we expect very little variance explained by status in many cases). 


#1st step is to decide the approach based on our data structure and distribution. 

#Fixed decissions: 
  #Need to include subsite as random effect (explicitly account for 1.repeated samplings)
  #Need to include sample_weights to account for strata sampling biass (implicitly, not modelled). 


#We want to reach an approach that is valid for all combinations of casepilots (all vs per-casepilot models) and GHG (co2, ch4, n2o, Co2+ch4, co2+ch4+n2o). Try first with LMMs, if not appropriate, fall back to GLMMs (potentially different approach for broad groups: different aproach for all vs per-casepilot or different approach for different GHGs). Ideally one model type for all combos. 


#checks for data (in all combinations for models)

#1. Normality of response variable (not essential but will suggest need for transformation)
#2. Normality of residuals (with the various potential transformations). Moderately important, severe non-normality can bias p-values and confidence intervals. 
#3. Homocedasticity of residuals across factors. Important for valid standard errors and p-values.
#4. Linearity. Residuals vs fitted should not have patterns. Important 
#5. Independence of observations, accounted for by random effect subsite

#Alternatives to LMMs:
  #GLMM: can deal with non-normal data (different distribution families) and heterocedasticity across predictor factors 

#We will evaluate the modelling approach for each GHG species sequentially, First with all casepilots together, then per-casepilot: 

#CO2: evaluate single general model
#Co2: evaluate per-casepilot models

#IN the end, I expect to need to use Generalized Linear Mixed Models (glmmTMB), as it is likely that we face:   #1. NON-normal distributions even after transformations (especially for CH4 and N2O )
  #2. Heterocedasticity across seasons (necessary for SE and p-value accuracy)
#There might be some outliers in CO2 (extreme positive and negative fluxes, check their origin)

#0.Contrast options-----
#NOTES on contrasts: in R by default, contrast is "contr.treatment" wich uses the first level of each factor as the "reference" level and all subsequent as "treatment" levels. This implies that, with default contrast type, when looking at the main effects of my model, what is shown for each main effect is whether there is a significant effect of 1 factor (eg. status) at the reference level of all other factors (i.e. season). What we want is to asses whether there is an overall average effect of status across all levels of season. For this purpose we need to set contrasts to "contr. sum", and always call Anova (model, type="III").

#Set contrasts to contr.sum (for the whole R session)
options(contrasts = c("contr.sum", "contr.poly"))


#0.Transformations------
#Any data-transformation must maintain the flux sign (emission vs capture), it has essential ecological information, at least for CO2. CH4 often follows a logarithmic distribution, but we might have to do a slight shift to get rid of very small negative values, then model with Gamma family.

#TO TRY after SEFS14---------
#Notes on pseudolog (~signed log) transformation: https://stratosida.github.io/regression-regrets/Pseudo_log_explainer.html#finding-a-parameter-that-best-achieves-normality
#Our type of CO2 distribution, is the result of substracting two log-normal distributions (emission-uptake), which leads to a mostly simetrical distribution with heavy tails and mostly centered around 0. 
# pseudo_log <-
#   function(x, sigma = 1, base = 10)
#     asinh(x / (2 * sigma)) / log(base)
# 
# inv_pseudo_log <-
#   function(x, sigma = 1, base = 10)
#     2 * sigma * sinh(x * log(base))
# #THis requires to optimice sigma for each dataset, which can make the distribution slimmer, fatter or even introduce bimodality. Only for modelling, not to use for visual representation. 
# #CA-CO2: 
# #example using CA.  DOES NOT WORK APPROPRIATELY, likely due to negatives (it should, adjust code.)
# x<- sort(data4models %>% filter(casepilot=="DA", ghgspecies=="gwp_co2") %>% filter(!is.na(dailyflux)) %>% pull(dailyflux)*10)
# #Creates pool of sigmas to chose best
# sigmas <- 2 ** seq(-10, 10, 0.2)
# 
# #Calculates correlation between normal distribution and my data untransformed. 
# origcor <- cor(qnorm((1:length(x) - 0.5) / length(x)), x)
# 
# #Recalculates correlation with normal distribution for all options of sigma.
# ncorx <-
#   sapply(sigmas, function(X)
#     cor(qnorm((1:length(
#       x
#     ) - 0.5) / length(x)), pseudo_log(x, X)))
# 
# cat("Optimal sigma: ")
# (optsigma <- sigmas[ncorx == max(ncorx)])
# 
# plot(log2(sigmas), ncorx, ylab = "Correlation with normal deviates", ylim =
#        c(0, 1))
# points(log2(sigmas)[ncorx == max(ncorx)], max(ncorx), pch = 19)
# abline(h = origcor)
# legend(
#   "bottomright",
#   lty = c(1, NA, NA),
#   pch = c(NA, 1, 19),
#   legend = c("Original", "Pseudo-log", "Pseudo-log with optimal sigma")
# )
# box()
# hist(pseudo_log(x, optsigma))


#Transformation functions (with inv_transformations for regaining og scale to interpret)
  #Signed log transformation:
  logsign <- function(x) sign(x) * log(abs(x) + 1)
  inv_logsign <- function(x) sign(x) * (exp(abs(x)) - 1)
  #Signed root transformation:
  rootsign <- function(x) sign(x) * sqrt(abs(x))
  inv_rootsign <- function(x) sign(x) * (x^2)
  # #Signed arcsinh transformation:
  arcsinhsign <- function(x) sign(x) * asinh(abs(x))
  inv_arcsinhsign <- function(x) sign(x) * sinh(abs(x))
  
  

  #Non-appropriate transformations: Box-Cox (non-negatives only), Yeo Johnson (sometimes swaps sign of very low-magnitude values, cannot be used separately for positive and negatives, no sign.YeoJohnson possible).
  
  #Yeo-Johnson transformation (sometimes swaps signs of low-magnitude values)
  # Load required package
  library(bestNormalize)
  
  # Yeo-Johnson transformation
  yeojohnson_transform <- function(x) {
    # Estimate lambda using bestNormalize
    lambda <- yeojohnson(x)$lambda
    yeojohnson(x, lambda = lambda)
  }
  
  # Inverse Yeo-Johnson transformation
  inv_yeojohnson <- function(y, lambda) {
    # Inverse based on the transformation formula
    ifelse(y >= 0,
           if (lambda == 0) exp(y) - 1 else (y * lambda + 1)^(1 / lambda) - 1,
           if (lambda == 2) -exp(-y) + 1 else -((-y * (2 - lambda) + 1)^(1 / (2 - lambda)) - 1))
  }
  
  # Wrapper to store lambda and transformed values
  yeojohnson_wrapper <- function(x) {
    result <- yeojohnson(x)
    attr(result$x.t, "lambda") <- result$lambda
    result$x.t
  }
  
  
  #SELECT AND INSPECT BEST TRANSFORMATION: 

# Lista de transformaciones
  transformations <- list(
    "Original" = identity,
    "Signed Log" = logsign,
    "Signed Sqrt" = rootsign,
    "Signed arcsin" = arcsinhsign,
    "Yeo-Johnson" = yeojohnson_wrapper
  )
  

# Evaluar normalidad con Shapiro-Wilk
evaluate_transformations <- function(x) {
  results <- lapply(names(transformations), function(name) {
    x_t <- transformations[[name]](x)
    shapiro <- tryCatch(shapiro.test(x_t), error = function(e) NULL)
    data.frame(
      Transformation = name,
      Shapiro_W = if (!is.null(shapiro)) shapiro$statistic else NA,
      Shapiro_p = if (!is.null(shapiro)) shapiro$p.value else NA
    )
  })
  bind_rows(results)
}

# Visualización: histogramas + densidad
plot_histograms <- function(x) {
  df <- lapply(names(transformations), function(name) {
    data.frame(
      value = transformations[[name]](x),
      transformation = name
    )
  }) %>% bind_rows()
  
  ggplot(df, aes(x = value)) +
    geom_histogram(aes(y = ..density..), bins = 30, fill = "skyblue", color = "black") +
    geom_density(color = "red", linewidth = 1) +
    facet_wrap(~transformation, scales = "free") +
    theme_minimal() +
    labs(title = "Histograms")
}

# Visualización: Q-Q plots
plot_qq <- function(x) {
  df <- lapply(names(transformations), function(name) {
    data.frame(
      value = transformations[[name]](x),
      transformation = name
    )
  }) %>% bind_rows()
  
  ggplot(df, aes(sample = value)) +
    stat_qq() +
    stat_qq_line() +
    facet_wrap(~transformation, scales = "free") +
    theme_minimal() +
    labs(title = "Q-Q Plots")
}

#For co2: sign-log (likely )
#For ch4: sign-log, likely glmm needed with family =student_t (allows positives and negatives).
#For GPW(co2+ch4): sign-log transformation. 
#For N2O: To be determined
#For GWP_all: To be determined

#We will repeat the checks for all combinations of casepilots and GHGs. 

  #Save results in summary table to be able to compare and reach common decission(s)
  

  #1st try lmer assumptions for GENERAL models (flux~status*casepilot*season)
  
  #Check data distributions and potential transformations for all casepilots together. 
  

#0. Custom formulas-------
#Here formulas to access relevant results from models. Will be used to summarise information.  

#We Cannot calculate R2 (conditional & marginal) for models with dispformula (or gamma family, unclear for t_family distributions), instead calculate pseudo-R2 (log-Likelyhood based) to get the "improvement in model fit relative to a null model (intercept + random effects), capturing the combined explanatory contribution of fixed effects". The null model must use the same data, random effect, dispformula and family distribution as the best model. This pseudoR2 is calculated as: 
#FUNCTION to get pseudoR2 for best_model:
pseudoR2_from_null <- function(best_model, data, group = "subsite", type = "mcfadden") {
  # Extract response variable
  response <- as.character(formula(best_model)[[2]])
  
  # Build null formula: response ~ 1 + (1 | group)
  null_formula <- as.formula(paste(response, "~ 1 + (1 |", group, ")"))
  
  # Extract other model components
  disp_formula <- best_model$call$dispformula
  ziformula <- best_model$call$ziformula
  family <- best_model$modelInfo$family
  
  # Fit null model
  null_model <- glmmTMB(
    formula = null_formula,
    data = data,
    dispformula = disp_formula,
    ziformula = ziformula,
    family = family
  )
  
  # Compute pseudo-R²
  logLik_null <- logLik(null_model)
  logLik_best <- logLik(best_model)
  
  if (type == "mcfadden") {
    r2 <- 1 - (as.numeric(logLik_best) / as.numeric(logLik_null))
  } else if (type == "nagelkerke") {
    n <- nobs(best_model)
    r2 <- (1 - exp((2 / n) * (logLik_null - logLik_best))) /
      (1 - exp((2 / n) * logLik_null))
  } else if (type == "tjur") {
    if (family$family != "binomial") stop("Tjur's R² only applies to binomial models.")
    p1 <- predict(best_model, type = "response")[data[[response]] == 1]
    p0 <- predict(best_model, type = "response")[data[[response]] == 0]
    r2 <- mean(p1) - mean(p0)
  } else {
    stop("Unsupported pseudo-R² type. Choose from 'mcfadden', 'nagelkerke', or 'tjur'.")
  }
  
  return(r2)
}



#Formula to get basic information about model structure from named list of models:
summarize_model_info <- function(model_list) {
  library(glmmTMB)
  library(performance)
  library(dplyr)
  library(purrr)
  
  map_df(names(model_list), function(cp) {
    model <- model_list[[cp]]
    
    r2_vals<- NULL
    # Try to extract R2 values
    r2_vals <- tryCatch({
      MuMIn::r.squaredGLMM(model) #Using MuMin to be able to get marginal R2 even with singularities
    }, error = function(e) {
      warning(paste("R2 extraction failed for", cp))
      return(c(NA_real_, NA_real_))
      
    })
    
    tibble(
      casepilot = cp,
      formula = deparse(formula(model)),
      dispformula = deparse(model$modelInfo$allForm$dispformula),
      family = model$modelInfo$family$family,
      R2_marginal = r2_vals[1,1],
      R2_conditional = r2_vals[1,2],
      is_singular= check_singularity(model)
    )
  })
}


#Summary of Model assumptions: for best_model_list, Requires named lists: best_model_list, data_list.

#We want to structure the results in a table: model structure,  distribution family, dispformula, 
#add pvalue for: normality of residuals, heteroscedasticity (levene test).

summarize_dharma_diagnostics <- function(model_list, data_list) {
  library(DHARMa)
  library(dplyr)
  library(purrr)
  
  map_df(names(model_list), function(cp) {
    model <- model_list[[cp]]
    data  <- data_list[[cp]]
    
    # Simulate residuals
    sim_res <- tryCatch({
      simulateResiduals(fittedModel = model, plot = FALSE)
    }, error = function(e) {
      warning(paste("DHARMa simulation failed for", cp))
      return(NULL)
    })
    
    if (is.null(sim_res)) {
      return(tibble(
        casepilot = cp,
        uniformity_pval = NA_real_,
        dispersion_pval = NA_real_,
        hetero_status_pval = NA_real_,
        hetero_season_pval = NA_real_
      ))
    }
    
    # Run tests
    tibble(
      casepilot = cp,
      formula = deparse(formula(model)),
      dispformula = deparse(model$modelInfo$allForm$dispformula),
      family = model$modelInfo$family$family,
      uniformity_pval = testUniformity(sim_res)$p.value,
      dispersion_pval = testDispersion(sim_res)$p.value,
      hetero_status_pval = testCategorical(sim_res,catPred =data$status,  plot = F)$homogeneity$`Pr(>F)`[1],
      hetero_season_pval = testCategorical(sim_res,catPred =data$season,  plot = F)$homogeneity$`Pr(>F)`[1]
    )
  })
}


#Function to summarise ANOVA results: get main effects significance for list of best_models
summarize_anova_results <- function(model_list) {
  library(car)
  library(glmmTMB)
  library(dplyr)
  library(purrr)
  
  map_df(names(model_list), function(cp) {
    model <- model_list[[cp]]
    
    # Run Type III ANOVA
    anova_res <- tryCatch({
      car::Anova(model, type = "III")
    }, error = function(e) {
      warning(paste("ANOVA failed for", cp))
      return(NULL)
    })
    
    # Extract p-values safely
    get_pval <- function(term) {
      if (is.null(anova_res)) return(NA_real_)
      if (term %in% rownames(anova_res)) return(anova_res[term, "Pr(>Chisq)"])
      return(NA_real_)
    }
    
    tibble(
      casepilot = cp,
      model_call = deparse(formula(model)),
      distribution_family = model$modelInfo$family$family,
      dispformula = deparse(model$modelInfo$allForm$dispformula),
      intercept_pval = get_pval("(Intercept)"),
      status_pval = get_pval("status"),
      season_pval = get_pval("season"),
      interaction_pval = get_pval("status:season")
    )
  })
}




#CO2_______ ------


#Check normality on CO2: datasets:
co2_table_trans<- tibble()
for (cp in unique(data4models$casepilot)) {
  for (ghg in c("gwp_co2")) { #"gwp_ch4","gwp_n2o","gwp_co2andch4"
    x<- data4models %>% filter(casepilot==cp&ghgspecies==ghg) %>% filter(!is.na(dailyflux)) %>% pull(dailyflux)
    table_trans<- evaluate_transformations(x)
    table_trans$casepilot<- cp
    table_trans$ghgspecies<- ghg
    co2_table_trans<- bind_rows(co2_table_trans, table_trans)
    print(plot_histograms(x)+ggtitle(paste(cp, ghg)))
    print(plot_qq(x)+ggtitle(paste(cp, ghg)))
  }
}

#get most-normal transformation for each case
co2_best_trans<- co2_table_trans %>% 
  group_by(ghgspecies,casepilot) %>% 
  filter(Shapiro_p==max(Shapiro_p)) %>% 
  arrange(ghgspecies,casepilot)
co2_best_trans

  
#CASEPILOT models-------
#STEPS: 
  #1. select best transformation (using lmer, no sampleweights)
  #2. Optimize model without sampleweights (lmer as first step --> lme (normal but heterocedasticity) or glmm (non-normal, heterocedasitcity or not))
  #3. Optimize model with sampleweights(glmm)



#1. Best common transformation ------ 
#Focus on normality of residuals, homocedasticity of residuals, model fit, for multiple transformations.

# Initialize summary table
co2_lmer_results <- tibble()

# Loop over all combinations of casepilot and CO2 using lmer models and residuals to check. 
for (cp in unique(data4models$casepilot)) {
  for (ghg in c("gwp_co2")) {
    
    cat("\n\n==================== CASEPILOT:", cp, "| GHG:", ghg, "====================\n")
    
    data_sub <- data4models %>% filter(grepl(cp, casepilot)& ghgspecies == ghg) %>% 
      filter(!is.na(dailyflux))
    
    #Transformations
    data_sub <- data_sub %>%
      mutate(
        dailyflux_logsign = logsign(dailyflux),
        dailyflux_rootsign = rootsign(dailyflux),
        dailyflux_arcsinhsign = arcsinhsign(dailyflux)
      )
    
    # Step 3: Fit lmer models
    model_untrans <- tryCatch(lmer(dailyflux ~ status * season + (1 | subsite), data = data_sub), error = function(e) NULL)
    model_logsign <- tryCatch(lmer(dailyflux_logsign ~ status * season + (1 | subsite), data = data_sub), error = function(e) NULL)
    model_rootsign <- tryCatch(lmer(dailyflux_rootsign ~ status * season + (1 | subsite), data = data_sub), error = function(e) NULL)
    model_arcsinhsign<- tryCatch(lmer(dailyflux_arcsinhsign ~ status * season + (1 | subsite), data = data_sub), error = function(e) NULL)
    
    
    # Step 4: Diagnostics
    check_model <- function(model, name) {
      if (is.null(model)) return(NULL)
      res <- residuals(model)
      shapiro_res <- shapiro.test(res)$p.value
      levene_status <- tryCatch(leveneTest(res ~ data_sub$status)$"Pr(>F)"[1], error = function(e) NA)
      levene_season <- tryCatch(leveneTest(res ~ data_sub$season)$"Pr(>F)"[1], error = function(e) NA)
      r2 <- tryCatch(performance::r2_nakagawa(model), error = function(e) NULL)
      aic_val <- tryCatch(AIC(model), error = function(e) NA)
      bic_val <- tryCatch(BIC(model), error = function(e) NA)
      rmse_cv <- tryCatch({
        ctrl <- trainControl(method = "cv", number = 10)
        cv_model <- train(
          x = model.matrix(model)[, -1],
          y = getME(model, "y"),
          method = "lm",
          trControl = ctrl
        )
        cv_model$results$RMSE
      }, error = function(e) NA)
      
      tibble(
        casepilot = cp,
        gas = ghg,
        model = name,
        shapiro_resid_p = shapiro_res,
        levene_status_p = levene_status,
        levene_season_p = levene_season,
        marginal_r2 = if (!is.null(r2)) r2$R2_marginal else NA,
        conditional_r2 = if (!is.null(r2)) r2$R2_conditional else NA,
        AIC = aic_val,
        BIC = bic_val,
        RMSE_CV = rmse_cv
      )
    }
    
    res_untrans <- check_model(model_untrans, "Untransformed")
    res_logsign <- check_model(model_logsign, "Signed Log")
    res_rootsign <- check_model(model_rootsign, "Signed Root")
    res_arcsinhsign <- check_model(model_arcsinhsign, "Signed Arcsinh")
    
    co2_lmer_results <- bind_rows(co2_lmer_results, res_untrans, res_logsign, res_rootsign, res_arcsinhsign)
  }
}

co2_lmer_results

#CO2 untransfromed data fails for DU and VA (singularities)
co2_lmer_results %>% 
  filter(is.na(conditional_r2))

#Normality of residuals only achieved for some transformations of RI and CA data
co2_lmer_results %>% 
  filter(shapiro_resid_p>0.05)

#Homocedasticity for both factors is never achieved... 
co2_lmer_results %>% 
  filter(levene_season_p>0.05&levene_status_p>0.05)

#Best transformation is consistently Signed Log:
co2_lmer_results %>% 
  group_by(casepilot,gas) %>% 
  filter(AIC==min(AIC))

#Arcsin might improve a bit normality, but not enough to justify change. 
co2_lmer_results %>% 
  group_by(casepilot,gas) %>% 
  filter(shapiro_resid_p==max(shapiro_resid_p))

#CONCLUSION: Lmer is not appropriate for modeling casepilots individually as is. We need to account for heterocedasticity. Heterocedasticity: We could use nlme::lme() and weights (syntax:  "weights = varIdent(form = ~1 | season)" ). BUT we have some non-normal residuals--> glmm for all (allows different family distributions and heterocedasticity)  

#We will use logsign-transformed data, it is the best transformation according to earlier tests with lmer (maximizes normality of residuals and minimizes heterocedasticity). 



#2. Optimice GLMM for each casepilot-------

# 3 out of 6 casepilots show strong non-normality... most show heterocedasticity for season. 

#Which distribution family is appropriate for each casepilot? Do we have heterocedasticity?
#I must account for heterocedasticity while testing the distribution family (it will change the residuals). 

#Normal residuals: CA, RI --> GAUSIAN family 
#Near-normality of residuals: DU --> GAUSIAN family
#Non-normality of residuals: CU, DA, VA ---> TBD (gausian vs t_family)

#Notes on singularity: to have a model with singularity means that some of the random effect terms likely have variance estimates of (almost) zero, or there is redundancy in the random structure. IN our case, the random effect is solely the subsite (to account for subsite-specific differences and repeated samplings across seasons). Our random effect structure is the most simple possible (one random effect for intercept), and we expect this random effect to account for variance when the 2 subsites of each status are fundamentally different. With this structure, we should not be too concerned about singuarity. The inclusion of the subsite as a random effect might not be statistically necesary (variance~0) but it is still conceptually justified.  



## 2.1. CA model ----
#CA: lmer showed normal residuals, heterocedasticity for season (not for status)
co2_lmer_results %>% filter(casepilot=="CA")
#approach: compare fit between gausian models, without heterocedasticity, with for season only and with for season and status. 

ca_co2<- data4models %>% filter(casepilot=="CA"&ghgspecies=="gwp_co2") %>% filter(!is.na(dailyflux)) %>% 
  mutate(dailyflux_logsign=logsign(dailyflux))

#Check effect of logsign tranformation: 
ca_co2 %>% ggplot(aes(x=dailyflux_logsign))+
  geom_histogram()+
  facet_grid(status~season)
#logsign transformation breaks distributions, even when the complete datasets improves normality: 
ca_co2 %>% ggplot(aes(x=dailyflux_logsign))+
  geom_histogram()


#Check data distribution for heterocedasticity across season and status: 
ca_co2 %>% ggplot(aes(x=dailyflux))+
  geom_histogram()+
  facet_grid(status~season)
#Yes, seems heterocedasticity across season and status.

ca_gaus_og<- glmmTMB(formula = dailyflux~status*season + (1|subsite),
                     data = ca_co2,
                     family = gaussian(),
                     dispformula = ~1 #default, homocedasticity assumed
)
res0<- simulateResiduals(ca_gaus_og, n=1000)

plot(res0)
plotResiduals(res0, ca_co2$season) #check heterocedasticity
plotResiduals(res0, ca_co2$status) #check heterocedasticity

ca_t_og<-glmmTMB(formula = dailyflux~status*season + (1|subsite),
                 data = ca_co2,
                 family = t_family,
                 dispformula = ~1 #default, homocedasticity assumed
)
res0_t<- simulateResiduals(ca_t_og, n=1000)
plot(res0_t)
plotQQunif(res0_t)
plotResiduals(res0_t) #PLot general residuals
plotResiduals(res0_t, ca_co2$season) #check heterocedasticity
 plotResiduals(res0_t, ca_co2$status) #check heterocedasticity
testDispersion(res0_t)
#Good QQ-plot, Heterocedasticity for status, strong Dispersion


ca_t_og1<- ca_t_og<-glmmTMB(formula = dailyflux~status*season + (1|subsite),
                            data = ca_co2,
                            family = t_family,
                            dispformula = ~season 
)
res0_t1<- simulateResiduals(ca_t_og1)
plot(res0_t1)
plotQQunif(res0_t1)
plotResiduals(res0_t1) #PLot general residuals
plotResiduals(res0_t1, ca_co2$season) #check heterocedasticity
plotResiduals(res0_t1, ca_co2$status) #check heterocedasticity
testDispersion(res0_t1)
performance::check_overdispersion(res0_t2)
#Good qq, but status heterocedasticity and non-uniformity

ca_t_og2<- ca_t_og<-glmmTMB(formula = dailyflux~status*season + (1|subsite),
                            data = ca_co2,
                            family = t_family,
                            dispformula = ~season + status
)
res0_t2<- simulateResiduals(ca_t_og2)
plot(res0_t2)
plotQQunif(res0_t2)
plotResiduals(res0_t2) #PLot general residuals
plotResiduals(res0_t2, ca_co2$season) #check heterocedasticity
plotResiduals(res0_t2, ca_co2$status) #check heterocedasticity
testDispersion(res0_t2)
performance::check_overdispersion(res0_t2)
check_convergence(ca_t_og2)
#Very good, only caveat: underdispersion (inflated SE, more conservative p-values, not super critical, might indicate overfitting or wrong family distribution)

anova(ca_t_og, ca_t_og1,ca_t_og2)



ca_glmm_gaus0<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite),
                        data = ca_co2,
                        family = gaussian(),
                        dispformula = ~1 #default, homocedasticity assumed
                        )

ca_glmm_gaus1<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite),
                        data = ca_co2,
                        family = gaussian(),
                        dispformula = ~season #heterocedasticity across season levels (not status level)
)
  
ca_glmm_gaus2<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite),
                        data = ca_co2,
                        family = gaussian(),
                        dispformula = ~season+status #heterocedasticity across season and status (no across interaction)
)
  
#Get residuals for all options
res0<- simulateResiduals(ca_glmm_gaus0)
res1<- simulateResiduals(ca_glmm_gaus1)
res2<- simulateResiduals(ca_glmm_gaus2)

#Explore residuals of model options
plot(res0) #check normality, 
plotResiduals(res0, ca_co2$season) #check heterocedasticity
plotResiduals(res0, ca_co2$status) #check heterocedasticity
testDispersion(res0) #checks wether residual variance is higher or lower than expected, n.s.--> no issue. 
testOutliers(res0) #Checks if there are more extreme residuals than expected
hist(res0$scaledResiduals)


plot(res1) #check normality, 
plotResiduals(res1, ca_co2$season) #check heterocedasticity
plotResiduals(res1, ca_co2$status) #check heterocedasticity
testDispersion(res1) #checks wether residual variance is higher or lower than expected, n.s.--> no issue. 
testOutliers(res1) #Checks if there are more extreme residuals than expected
hist(res1$scaledResiduals)
check_singularity(res1)


plot(res2) #check normality, 
plotResiduals(res2, ca_co2$season) #check heterocedasticity
plotResiduals(res2, ca_co2$status) #check heterocedasticity
testDispersion(res2) #checks wether residual variance is higher or lower than expected, n.s.--> no issue. 
testOutliers(res2) #Checks if there are more extreme residuals than expected
hist(res2$scaledResiduals)

#COMPARE model options
anova(ca_glmm_gaus0,ca_glmm_gaus1,ca_glmm_gaus2)

#CA best model: gaussian, allowing heterocedasticity across seasons. 
ca_best_co2<- ca_glmm_gaus1
Anova(ca_best_co2,type = 3) 

#Cannot calculate R2 (conditional & marginal) for models with heterocedasticity,
#Get PseudoR2
pseudoR2_ca_co2<-pseudoR2_from_null(ca_best_co2, ca_co2)
pseudoR2_ca_co2

#R2 of homocedastic model
ca_best_homoc_co2<- ca_glmm_gaus0
r2(ca_best_homoc_co2)


#calculate pseudoR2 from obs vs predicted
actual<- ca_co2$dailyflux_logsign
predicted<- predict(ca_best_co2, type="response")

plot(actual, predicted)
pseudo_r2_ca_co2_Efron<- 1 - sum((actual-predicted)^2)/sum((actual-mean(actual))^2)






rm(ca_glmm_gaus0,ca_glmm_gaus1,ca_glmm_gaus2)

## 2.2. CU model ----
#CU: lmer showed non-normal residuals, heterocedasticity for season and status, and extremely low r2 (~0.06)
co2_lmer_results %>% filter(casepilot=="CU")
#approach: compare fit between gausian and tstudent models without heterocedasticity, with for season only and with for season and status (not for interaction). 

#CU: non-normal residuals
cu_co2<- data4models %>% filter(casepilot=="CU"&ghgspecies=="gwp_co2") %>% filter(!is.na(dailyflux)) %>% 
  mutate(dailyflux_logsign=logsign(dailyflux))

#Check distribution of data: Two near-normal groups (OW and V)
plot_histograms(cu_co2$dailyflux)

cu_glmm_gaus0<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite), 
                       data = cu_co2,
                       family = gaussian(),
                       dispformula = ~1)

cu_glmm_gaus1<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite), 
                        data = cu_co2,
                        family = gaussian(),
                        dispformula = ~season)

cu_glmm_gaus2<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite), 
                        data = cu_co2,
                        family = gaussian(),
                        dispformula = ~season + status)

cu_glmm_stu0<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite), 
                        data = cu_co2,
                        family = t_family,
                        dispformula = ~1)

cu_glmm_stu1<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite), 
                        data = cu_co2,
                        family = t_family,
                        dispformula = ~season)

cu_glmm_stu2<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite), 
                        data = cu_co2,
                        family = t_family,
                        dispformula = ~season + status)

resgaus0<- simulateResiduals(cu_glmm_gaus0, n=1000)
resgaus1<- simulateResiduals(cu_glmm_gaus1, n=1000)
resgaus2<- simulateResiduals(cu_glmm_gaus2, n=1000)
restu0<- simulateResiduals(cu_glmm_stu0, n=1000)
restu1<- simulateResiduals(cu_glmm_stu1, n=1000)
restu2<- simulateResiduals(cu_glmm_stu2, n=1000)

plot(resgaus0) #check normality, 
plotResiduals(resgaus0, cu_co2$season) #check heterocedasticity
plotResiduals(resgaus0, cu_co2$status) #check heterocedasticity
plot(resgaus1) #check normality, 
plotResiduals(resgaus1, cu_co2$season) #check heterocedasticity
plotResiduals(resgaus1, cu_co2$status) #check heterocedasticity

plot(restu0) #check normality, 
plotResiduals(restu0, cu_co2$season) #check heterocedasticity
plotResiduals(restu0, cu_co2$status) #check heterocedasticity
plot(restu1) #check normality, 
plotResiduals(restu1, cu_co2$season) #check heterocedasticity
plotResiduals(restu1, cu_co2$status) #check heterocedasticity

plot(restu2) #check normality, 
plotResiduals(restu2, cu_co2$season) #check heterocedasticity
plotResiduals(restu2, cu_co2$status) #check heterocedasticity


testDispersion(restu1) #checks wether residual variance is higher or lower than expected, n.s.--> no issue. 
testOutliers(restu2) #Checks if there are more extreme residuals than expected

#ANOVA to decide: 
anova(cu_glmm_gaus0,cu_glmm_gaus1,cu_glmm_gaus2,cu_glmm_stu0,cu_glmm_stu1,cu_glmm_stu2)
#Student distribution is consistently better 
anova(cu_glmm_stu0,cu_glmm_stu1,cu_glmm_stu2)
#Big improvemnt with season heterocedasticity, limited improvement with additional status heterocedasticity. 
anova(cu_glmm_stu1,cu_glmm_stu2)


#CU best model: student, allowing heterocedasticity across seasons only. 

cu_best_co2<- cu_glmm_stu2
Anova(cu_best_co2,type = 3) 

#Cannot calculate R2 (conditional & marginal) for models with heterocedasticity,
#Get PseudoR2
pseudoR2_cu_co2<-pseudoR2_from_null(cu_best_co2, cu_co2)
pseudoR2_cu_co2


#R2 of homocedastic model
cu_best_homoc_co2<- cu_glmm_stu0
r2(cu_best_homoc_co2)
MuMIn::r.squaredGLMM(cu_best_homoc_co2) #random effect R2 is 0 (singularity), R2m==R2c

#CU model has singularity (regardless of family distribution used)
performance::check_singularity(cu_glmm_stu0)
performance::check_singularity(cu_glmm_gaus0)

#DOUBT singularity CU:
#This means that the two subsites of each status behave similarly (they are good "replicates"). Even when is not statistically needed, we include it for consistency and for conceptual reasons. Singularity does not impact interpretation but will limit some functionalities (r2).  

#calculate pseudoR2 from obs vs predicted
actual<- cu_co2$dailyflux_logsign
predicted<- predict(cu_best_co2, type="response")

plot(actual, predicted)
pseudo_r2_cu_co2_Efron<- 1 - sum((actual-predicted)^2)/sum((actual-mean(actual))^2)




rm(cu_glmm_gaus0,cu_glmm_gaus1,cu_glmm_gaus2,cu_glmm_stu0,cu_glmm_stu1,cu_glmm_stu2)


## 2.3. DA model ----
#DA: lmer showed non-normal residuals, heterocedasticity for season (and status),
co2_lmer_results %>% filter(casepilot=="DA")
#approach: compare fit between gausian and tstudent models without heterocedasticity, with for season only and with for season and status (not for interaction). 

#DA: non-normal residuals
da_co2<- data4models %>% filter(casepilot=="DA"&ghgspecies=="gwp_co2") %>% filter(!is.na(dailyflux)) %>% 
  mutate(dailyflux_logsign=logsign(dailyflux))

plot_histograms(da_co2$dailyflux)

#calculate all model options
da_glmm_gaus0<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite), 
                        data = da_co2,
                        family = gaussian(),
                        dispformula = ~1)

da_glmm_gaus1<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite), 
                        data = da_co2,
                        family = gaussian(),
                        dispformula = ~season)

da_glmm_gaus2<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite), 
                        data = da_co2,
                        family = gaussian(),
                        dispformula = ~season + status)

da_glmm_stu0<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite), 
                       data = da_co2,
                       family = t_family,
                       dispformula = ~1)

da_glmm_stu1<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite), 
                       data = da_co2,
                       family = t_family,
                       dispformula = ~season)

da_glmm_stu2<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite), 
                       data = da_co2,
                       family = t_family,
                       dispformula = ~season + status)

#Get residuals for all model options
resgaus0<- simulateResiduals(da_glmm_gaus0, n=1000)
resgaus1<- simulateResiduals(da_glmm_gaus1, n=1000)
resgaus2<- simulateResiduals(da_glmm_gaus2, n=1000)
restu0<- simulateResiduals(da_glmm_stu0, n=1000)
restu1<- simulateResiduals(da_glmm_stu1, n=1000)
restu2<- simulateResiduals(da_glmm_stu2, n=1000)

plot(resgaus0) #check normality, 
plotResiduals(resgaus0, da_co2$season) #check heterocedasticity
plotResiduals(resgaus0, da_co2$status) #check heterocedasticity
plot(resgaus1) #check normality, 
plotResiduals(resgaus1, da_co2$season) #check heterocedasticity
plotResiduals(resgaus1, da_co2$status) #check heterocedasticity

plot(restu0) #check normality, 
plotResiduals(restu0, da_co2$season) #check heterocedasticity
plotResiduals(restu0, da_co2$status) #check heterocedasticity
plot(restu1) #check normality, 
plotResiduals(restu1, da_co2$season) #check heterocedasticity
plotResiduals(restu1, da_co2$status) #check heterocedasticity

plot(restu2) #check normality, 
plotResiduals(restu2, da_co2$season) #check heterocedasticity
plotResiduals(restu2, da_co2$status) #check heterocedasticity


testDispersion(restu1) #checks wether residual variance is higher or lower than expected, n.s.--> no issue. 
testOutliers(restu1) #Checks if there are more extreme residuals than expected

#ANOVA to decide: 
anova(da_glmm_gaus0,da_glmm_gaus1,da_glmm_gaus2,da_glmm_stu0,da_glmm_stu1,da_glmm_stu2)
#Student distribution is better always
anova(da_glmm_stu0,da_glmm_stu1,da_glmm_stu2)
#Big improvemnt with season heterocedasticity, limited improvement with additional status heterocedasticity. 
anova(da_glmm_stu1,da_glmm_stu2)


#DA best model: student, allowing heterocedasticity across seasons only. 

da_best_co2<- da_glmm_stu1
Anova(da_best_co2,type = 3)

#Cannot calculate R2 (conditional & marginal) for models with heterocedasticity,
#Get PseudoR2
pseudoR2_da_co2<-pseudoR2_from_null(da_best_co2, da_co2)
pseudoR2_da_co2

#R2 of homocedastic model
da_best_homoc_co2<- da_glmm_stu0
r2(da_best_homoc_co2)
MuMIn::r.squaredGLMM(da_best_homoc_co2) #random effect R2 is 0 (singularity), R2m==R2c
#DA model has singularity (regardless of family distribution used)
performance::check_singularity(da_glmm_stu0)
performance::check_singularity(da_glmm_gaus0)

#DOUBT singularity DA:
#This means that the two subsites of each status behave similarly (they are good "replicates"). Even when is not statistically needed, we include it for consistency and for conceptual reasons. Singularity does not impact interpretation but will limit some functionalities (r2). 
#We expected A1 and A2 subsites to have different CO2 profiles and this to be absorved by the random effect, but apparently not so much (subsite explains no variance)

rm(da_glmm_gaus0,da_glmm_gaus1,da_glmm_gaus2,da_glmm_stu0,da_glmm_stu1,da_glmm_stu2)

## 2.4. DU model ----
#DU: lmer showed near-normal residuals, heterocedasticity for season (not for status)
co2_lmer_results %>% filter(casepilot=="DU")
#approach: compare fit between gausian models: without heterocedasticity, with for season only and with for season and status. 

du_co2<- data4models %>% filter(casepilot=="DU"&ghgspecies=="gwp_co2") %>% filter(!is.na(dailyflux)) %>% 
  mutate(dailyflux_logsign=logsign(dailyflux))

plot_histograms(du_co2$dailyflux)


du_glmm_gaus0<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite),
                        data = du_co2,
                        family = gaussian(),
                        dispformula = ~1 #default, homocedasticity assumed
)

du_glmm_gaus1<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite),
                        data = du_co2,
                        family = gaussian(),
                        dispformula = ~season #heterocedasticity across season levels (not status level)
)

du_glmm_gaus2<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite),
                        data = du_co2,
                        family = gaussian(),
                        dispformula = ~season+status #heterocedasticity across season and status (no across interaction)
)

#Get residuals for all options
res0<- simulateResiduals(du_glmm_gaus0, n=1000)
res1<- simulateResiduals(du_glmm_gaus1, n=1000)
res2<- simulateResiduals(du_glmm_gaus2, n=1000)

#Explore residuals of model options
plot(res0) #check normality, 
plotResiduals(res0, du_co2$season) #check heterocedasticity
plotResiduals(res0, du_co2$status) #check heterocedasticity

plot(res1) #check normality, 
plotResiduals(res1, du_co2$season) #check heterocedasticity
plotResiduals(res1, du_co2$status) #check heterocedasticity
testDispersion(res1) #checks wether residual variance is higher or lower than expected, n.s.--> no issue. 
testOutliers(res1) #Checks if there are more extreme residuals than expected
hist(res1$scaledResiduals)


#COMPARE model options: Heterocedasticity for season improves model performance. , but improves AIC significantly
anova(du_glmm_gaus0,du_glmm_gaus1,du_glmm_gaus2)

#DU best model: gaussian, allowing heterocedasticity across seasons. 
du_best_co2<- du_glmm_gaus1
Anova(du_best_co2,type = 3)

#Cannot calculate R2 (conditional & marginal) for models with heterocedasticity,
#Get PseudoR2
pseudoR2_du_co2<-pseudoR2_from_null(du_best_co2, du_co2)
pseudoR2_du_co2

#R2 of homocedastic model
du_best_homoc_co2<-du_glmm_gaus0
r2(du_best_homoc_co2)

rm(du_glmm_gaus0,du_glmm_gaus1,du_glmm_gaus2)

## 2.5. RI model ----
#RI: lmer showed normal residuals, heterocedasticity for status (but not for season)
co2_lmer_results %>% filter(casepilot=="RI")
#approach: compare fit between gausian models: without heterocedasticity, with for season only and with for season and status. 
##___TO_TRY-------
#Compare logsign transformed (gaussian) and untransformed (t-family) models, use the one with best fit (AIC focus only)

#Subset Co2 data for casepilot
ri_co2<- data4models %>% filter(casepilot=="RI"&ghgspecies=="gwp_co2") %>% filter(!is.na(dailyflux)) %>% 
  mutate(dailyflux_logsign=logsign(dailyflux)) %>% 
  filter(between(dailyflux, -7.5,2.5 )) #Remove 3 outliers

plot_histograms(ri_co2$dailyflux)

ri_co2 %>% ggplot(aes(x=status,y=dailyflux,col=status))+
  geom_violin()+
  ggforce::geom_sina()

#Fit model options
ri_glmm_gaus0<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite),
                        data = ri_co2,
                        family = gaussian(),
                        dispformula = ~1 #default, homocedasticity assumed
)

ri_glmm_gaus1<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite),
                        data = ri_co2,
                        family = gaussian(),
                        dispformula = ~status #heterocedasticity across status levels (not seasons level)
)

ri_glmm_gaus2<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite),
                        data = ri_co2,
                        family = gaussian(),
                        dispformula = ~season+status #heterocedasticity across season and status (no across interaction)
)

#Get residuals for all options
res0<- simulateResiduals(ri_glmm_gaus0)
res1<- simulateResiduals(ri_glmm_gaus1)
res2<- simulateResiduals(ri_glmm_gaus2)

#Explore residuals of model options
plot(res0) #check normality, 
plotResiduals(res0, ri_co2$season) #check heterocedasticity
plotResiduals(res0, ri_co2$status) #check heterocedasticity
testDispersion(res0) #checks wether residual variance is higher or lower than expected, n.s.--> no issue. 
testOutliers(res0) #Checks if there are more extreme residuals than expected
hist(res0$scaledResiduals)


plot(res1) #check normality, 
plotResiduals(res1, ri_co2$season) #check heterocedasticity
plotResiduals(res1, ri_co2$status) #check heterocedasticity
testDispersion(res1) #checks wether residual variance is higher or lower than expected, n.s.--> no issue. 
testOutliers(res1) #Checks if there are more extreme residuals than expected
hist(res1$scaledResiduals)


#BEst model is accounting for heterocedasticity across status, but very little improvement 
anova(ri_glmm_gaus0,ri_glmm_gaus1,ri_glmm_gaus2)


#RI best model: gaussian, assuming homocedasticity
ri_best_co2<- ri_glmm_gaus1
Anova(ri_best_co2,type = 3)

#Explained variance: calculate R2 (conditional & marginal)
r2(ri_best_co2)

#Calculate pseudoR2 (for consistency across models, even when ri_best_co2 has real R2)
#Get PseudoR2
pseudoR2_ri_co2<-pseudoR2_from_null(ri_best_co2, ri_co2)
pseudoR2_ri_co2

#best model is homocedastic:
ri_best_homoc_co2<- ri_glmm_gaus0
r2(ri_best_homoc_co2)
rm(ri_glmm_gaus0,ri_glmm_gaus1,ri_glmm_gaus2)

## 2.6. VA model ----
#VA: lmer showed non-normal residuals, heterocedasticity for season (and status?),
co2_lmer_results %>% filter(casepilot=="VA")
#approach: compare fit between gausian and tstudent models without heterocedasticity, with for season only and with for season and status (not for interaction). 

#VA: non-normal residuals
va_co2<- data4models %>% filter(casepilot=="VA"&ghgspecies=="gwp_co2") %>% filter(!is.na(dailyflux)) %>% 
  mutate(dailyflux_logsign=logsign(dailyflux))

plot_qq(va_co2$dailyflux)
plot_histograms(va_co2$dailyflux) #BIMODAL distribution (sligthly)

#calculate all model options
va_glmm_gaus0<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite), 
                        data = va_co2,
                        family = gaussian(),
                        dispformula = ~1)

va_glmm_gaus1<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite), 
                        data = va_co2,
                        family = gaussian(),
                        dispformula = ~season)

va_glmm_gaus2<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite), 
                        data = va_co2,
                        family = gaussian(),
                        dispformula = ~season + status)

va_glmm_stu0<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite), 
                       data = va_co2,
                       family = t_family,
                       dispformula = ~1)

va_glmm_stu1<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite), 
                       data = va_co2,
                       family = t_family,
                       dispformula = ~season)

va_glmm_stu2<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite), 
                       data = va_co2,
                       family = t_family,
                       dispformula = ~season + status)

#Get residuals for all model options
resgaus0<- simulateResiduals(va_glmm_gaus0, n=1000)
resgaus1<- simulateResiduals(va_glmm_gaus1, n=1000)
resgaus2<- simulateResiduals(va_glmm_gaus2, n=1000)
restu0<- simulateResiduals(va_glmm_stu0, n=1000)
restu1<- simulateResiduals(va_glmm_stu1, n=1000)
restu2<- simulateResiduals(va_glmm_stu2, n=1000)

plot(resgaus0) #check normality, 
plotResiduals(resgaus0, va_co2$season) #check heterocedasticity
plotResiduals(resgaus0, va_co2$status) #check heterocedasticity
plot(resgaus1) #check normality, 
plotResiduals(resgaus1, va_co2$season) #check heterocedasticity
plotResiduals(resgaus1, va_co2$status) #check heterocedasticity

plot(restu0) #check normality, 
plotResiduals(restu0, va_co2$season) #check heterocedasticity
plotResiduals(restu0, va_co2$status) #check heterocedasticity
plot(restu1) #check normality, 
plotResiduals(restu1, va_co2$season) #check heterocedasticity
plotResiduals(restu1, va_co2$status) #check heterocedasticity

plot(restu2) #check normality, 
plotResiduals(restu2, va_co2$season) #check heterocedasticity
plotResiduals(restu2, va_co2$status) #check heterocedasticity


testDispersion(resgaus1) #checks wether residual variance is higher or lower than expected, n.s.--> no issue. 
testOutliers(resgaus1) #Checks if there are more extreme residuals than expected
hist(resgaus1$scaledResiduals)

#ANOVA to decide: 
anova(va_glmm_gaus0,va_glmm_gaus1,va_glmm_gaus2,va_glmm_stu1,va_glmm_stu2)
#Student distribution does not improve fit
anova(va_glmm_gaus0,va_glmm_gaus1,va_glmm_gaus2)
#Big improvemnt with season heterocedasticity, limited improvement with additional status heterocedasticity. 
anova(va_glmm_gaus1,va_glmm_gaus2)


#VA best model: gaus, allowing heterocedasticity across seasons only. 

va_best_co2<- va_glmm_gaus1
Anova(va_best_co2,type = 3)

#Cannot calculate R2 (conditional & marginal) for models with heterocedasticity,
#Get PseudoR2
pseudoR2_va_co2<-pseudoR2_from_null(va_best_co2, va_co2)
pseudoR2_va_co2

##___DOUBT low R2 for VA------
#R2 of homocedastic model
va_best_homoc_co2<- va_glmm_gaus0
r2(va_best_homoc_co2)
#Low R2, no status effect detected

#VA model does not have singularity, 
performance::check_singularity(va_glmm_gaus0)
performance::check_singularity(va_glmm_gaus1)

rm(va_glmm_gaus0,va_glmm_gaus1,va_glmm_gaus2,va_glmm_stu1,va_glmm_stu2)

#3.(skip) Optimize with weights--------

#Here we will create glmms taking into account strata sample_weights. 
#Then compare models with and without sample_weights to see the effect of weights and get the overall best model for each casepilot. 

# DOUBT: should I directly implement the identified best_model option for each casepilot adding sample_weight OR optimize from the beginning (allowing different heterocedasticity structures)??:

  #Applying weights to the already optimized structure can tell us directly if sample weights improve fit. 
  #Optimizing the models from the beginning with sample weights will allow us to check whether the weights influence the overall structure selection BUT we will not be able to directly compare with anova() if models differ in structure (i.e. if they have different fixed/random effects, different dispformula, or different distribution families) 

#We know that we should account for sample weights in the end if/when we get accurate strata composition and if they show substantial and consistent sampling bias. 

## 3.1. CA model (weighted)

#Do sample weights improve the optimized model? 
#Recal best model structure: 
summary(ca_best)

#build the same model with weights
ca_best_weighted<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite),
                           data = ca_co2,
                           family = gaussian(),
                           dispformula = ~season, #heterocedasticity across season levels (not status level)
                           weights=sample_weight)

#Compare fit between weigthed and unweighted model
anova(ca_best, ca_best_weighted)
#inclusion of sample_weight does not significantly improve model. Unclear effect on Df but weights cause no improvement on AIC...

Anova(ca_best,type = 3)
Anova(ca_best_weighted,type = 3)
#Weights also do not change significance of main effects. 

#Do weights lead to different structure? #Run the optimization process to select best model (with/without heterocedasticity) Not essential to check, omit. 


## 3.2. CU model (weighted)

#Do sample weights improve the optimized model? 
#Recal best model structure: 
summary(cu_best)

#build the same model with weights
cu_best_weighted<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite),
                           data = cu_co2,
                           family = t_family,
                           dispformula = ~season, #heterocedasticity across season levels (not status level)
                           weights=sample_weight)

#Compare fit between weigthed and unweighted model
anova(cu_best, cu_best_weighted)
#inclusion of sample_weight does significantly improve model. Unclear effect on Df but weights cause  improvement on AIC...

Anova(cu_best,type = 3)
Anova(cu_best_weighted,type = 3)
#Weights make main effects sigltly more significant.
summary(cu_best_weighted)

#Do weights lead to different structure? #Run the optimization process to select best model (with/without heterocedasticity) Not essential to check, omit. 



## 3.3. DA model (weighted)

#Do sample weights improve the optimized model? 
#Recal best model structure: 
summary(da_best)

#build the same model with weights
da_best_weighted<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite),
                           data = da_co2,
                           family = t_family,
                           dispformula = ~season, #heterocedasticity across season levels (not status level)
                           weights=sample_weight)

#Compare fit between weigthed and unweighted model
anova(da_best, da_best_weighted)
#inclusion of sample_weight does NOT significantly improve model. Unclear effect on Df....

Anova(da_best,type = 3)
Anova(da_best_weighted, type = 3)
#Weights make main effects sigltly less significant.

#Do weights lead to different structure? #Run the optimization process to select best model (with/without heterocedasticity) Not essential to check, omit. 



## 3.4. DU model (weighted)

#Do sample weights improve the optimized model? 
#Recal best model structure: 
summary(du_best)

#build the same model with weights
du_best_weighted<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite),
                           data = du_co2,
                           family = gaussian(),
                           dispformula = ~season, #heterocedasticity across season levels (not status level)
                           weights=sample_weight)

#Compare fit between weigthed and unweighted model
anova(du_best, du_best_weighted)
#inclusion of sample_weight does NOT significantly improve model. Unclear effect on Df....

#No detectable effect of status
Anova(du_best,type = 3)
Anova(du_best_weighted, type = 3)

#Do weights lead to different structure? #Run the optimization process to select best model (with/without heterocedasticity) Not essential to check, omit. 

## 3.5. RI model (weighted)

#Do sample weights improve the optimized model? 
#Recal best model structure: 
summary(ri_best)

#build the same model with weights
ri_best_weighted<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite),
                           data = ri_co2,
                           family = gaussian(),
                           # dispformula = ~season, #HOMOCEDASTICITY assumed
                           weights=sample_weight)

#Compare fit between weigthed and unweighted model
anova(ri_best, ri_best_weighted)
#inclusion of sample_weight does NOT significantly improve model. Unclear effect on Df....


#Same significance of main effects
Anova(ri_best,type = 3)
Anova(ri_best_weighted,type = 3)

#Do weights lead to different structure? #Run the optimization process to select best model (with/without heterocedasticity) Not essential to check, omit. 


## 3.6. VA model (weighted)

#Do sample weights improve the optimized model? 
#Recal best model structure: 
summary(va_best)

#build the same model with weights
va_best_weighted<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite),
                           data = va_co2,
                           family = gaussian(),
                           dispformula = ~season, #HOMOCEDASTICITY assumed
                           weights=sample_weight)

#Compare fit between weigthed and unweighted model
anova(va_best, va_best_weighted)
#inclusion of sample_weight does NOT significantly improve model. Unclear effect on Df....


#Same significance of main effects
Anova(va_best,type = 3)
Anova(va_best_weighted, type = 3)

#Do weights lead to different structure? #Run the optimization process to select best model (with/without heterocedasticity) Not essential to check, omit. 



#4. Effects of best-models-------

#Use best models without strata sample_weights for the moment. 

#From last step: only model in curonian lagoon benefits from sample weights, this is expected as we forced the strata distribution to be fixed across different samplings (most appropriate option to get overall effect). It is the only subsite where strata distribution has to be especially accounted for in comparisons. The rest of casepilots might benefit from weights once we have the actual distribution of them. Good thing is: results do not change much, so our incorporation of weights is not too critical (if criticized, it could be ommited and resuts would hold) 


#BASED ON THE ABOVE OPTIMIZATIONS, provide list of models to evaluate/interpret

#Provide list of homocedastic model variants (to estimate R2) THIS WILL BE PROBLEMATIC FOR MODELS WITH SINGULARITIES... 
homoced_model_list_co2<- list(
  "CA" = ca_best_homoc_co2,
  "CU" = cu_best_homoc_co2, 
  "DA" = da_best_homoc_co2,
  "DU" = du_best_homoc_co2,
  "RI" = ri_best_homoc_co2, #Same as best_model_co2, Only casepilot for which homocedastic model is best
  "VA" = va_best_homoc_co2)

#List of best models for each casepilot (to get residual tests and Significances)
best_model_list_co2<- list(
  "CA" = ca_best_co2,
  "CU" = cu_best_co2, #Only model for which sample_weight improved fit (omit weights for the moment)
  "DA" = da_best_co2,
  "DU" = du_best_co2,
  "RI" = ri_best_co2,
  "VA" = va_best_co2)

#List of dataframes used to build each model (needed to calculate residuals)
model_data_list_co2<- list(
  "CA" = ca_co2,
  "CU" = cu_co2, 
  "DA" = da_co2,
  "DU" = du_co2,
  "RI" = ri_co2,
  "VA" = va_co2)

## 4.1. Overall fit -----
#1st how good is the model--> overall R2 for homocedastic model variant 
co2_model_fit_summary<- summarize_model_info(homoced_model_list_co2)%>% 
  mutate(ghgspecies="gwp_co2")

co2_model_fit_summary
## 4.2. Residual checks---
#Get DHARMa diagnostics for list of bestmodels 

co2_resid_diag_summary <- summarize_dharma_diagnostics(best_model_list_co2, model_data_list_co2)%>% 
  mutate(ghgspecies="gwp_co2")

print(co2_resid_diag_summary)

#Residual tests are ok for most models with 2 exceptions: 
  #dispersion for CU model (this model is always singular, and sample weights are not accounted for in these tests, likely ok) 
  #heterocedasticity for status and season for VA model (this model is very bad: extremely low r2, and even acounting for season heterocedasticity, residuals are still not good, we cannot obtain a interpretable model for VA co2)


##4.2. Significance of effects-----
#Effect of status will first be assessed via Anova (best_model,type=3): is there an effect averaging across all seasons? yes/no, is the effect dependent on season? 


#Obtain significance for main effects of best_models
co2_results_anova<-summarize_anova_results(best_model_list_co2) %>% 
  mutate(ghgspecies="gwp_co2")
print(co2_results_anova)

#Intercept refers to whether the grand mean (on the transformed scale) is significantly different than cero (Is the overall net exchange of CO2 balanced between uptake and emission)
#Status refers to whether the status (averaged across seasons) has a significant effect on flux. 
#Season refers to whether the season (averaged across status) has a significant effect on flux.
#Interaction refers to whether the effect of status depends on season — i.e., the difference between status levels changes across seasons.



#Look for significance at either status or status:season---> there is an effect (that may or may not be evident when averaging across seasons)
#IF status:season significant---> effect depends on season. 

#Does status have an effect (significance of status or interaction)
co2_results_anova %>% filter(status_pval<0.05|interaction_pval<0.05)

# CA, CU , DA and RI show significant effect of interaction (effect of status dependent on season)

#DU and VA show no effect (main or interaction of status)

#Save in dropbox tables with the models summary outputs: 

write.csv(x = co2_model_fit_summary, file = paste0(plots_path,"/CO2models_homocedastic_R2.csv"),row.names = F)
write.csv(x = co2_resid_diag_summary, file = paste0(plots_path,"/CO2models_tests_residuals.csv"),row.names = F)
write.csv(x = co2_results_anova, file = paste0(plots_path,"/CO2models_effects_signififcance.csv"),row.names = F)

##4.3. Effect direction and magnitude----

#Use model scale for statistics, back-transform for interpretation.  

#NOTE on bootstrapping: there is the possibility to bootstrap the models (i.e. recalculate the model many times with re-sampled data, and calculate emmeans for each model, to be able to produce more robust Confidence Intervals for them, but this is not usually done and computationally intensive. We would have to take into account the model structure when re-sampling the data to ensure balanced resamplings (all factors represented by data). WE will not do this (for the moment).

#Some guides to plot and use emmeans: https://rcompanion.org/handbook/G_06.html


#Then emmeans comparisons to get direction of effect and magnitude (averaging over season interaction if any)
plot(emmeans(ri_best_co2, ~status), comparisons=TRUE) #Unclear interpretation of arrows, not the best way of presenting it. 


#Notes on non-gaussian models: upper and lower confidence limits are given as asymp.LCL/asymp.UCL (in our context we can treat them equally as the ones from gaussian models). Additionally, non-gaussian models will have Inf df (degrees of freedom are calculated differently for T-family models), not an issue. 

# Initialize list for storing comparisons
comparisons_list <- list()

# Loop to extract all comparisons for every casepilot best model
for (casepilot_name in names(best_model_list_co2)) {
  ghgspecies<- "gwp_co2"
  cp_model <- best_model_list_co2[[casepilot_name]]
  
  # Status comparisons:
  status_emmeans <- emmeans(cp_model, ~status)
  status_comparisons <- cld(status_emmeans, 
                            alpha = 0.05,
                            Letters = letters,
                            adjust = "sidak") %>%
    mutate(comparison = "status",
           ghgspecies = ghgspecies,
           casepilot = casepilot_name,
           season = NA)%>%    
    #Rename confidence limits for T_family models (if this is the case)
    rename(lower.CL = any_of("asymp.LCL"),
           upper.CL = any_of("asymp.UCL"))
  
  # Season comparisons: 
  season_emmeans <- emmeans(cp_model, ~season)
  season_comparisons <- cld(season_emmeans, 
                            alpha = 0.05,
                            Letters = letters,
                            adjust = "sidak") %>%
    mutate(comparison = "season",
           ghgspecies = ghgspecies,
           casepilot = casepilot_name,
           status = NA)%>%
    #Rename confidence limits for T_family models (if this is the case)
    rename(lower.CL = any_of("asymp.LCL"),
           upper.CL = any_of("asymp.UCL"))
  
  
  # Interaction comparisons: 
  interaction_emmeans <- emmeans(cp_model, ~status:season)
  interaction_comparisons <- cld(interaction_emmeans, 
                                 alpha = 0.05,
                                 Letters = letters,
                                 adjust = "sidak") %>%
    mutate(comparison = "interaction",
           ghgspecies = ghgspecies,
           casepilot = casepilot_name)%>%
    #Rename confidence limits for T_family models (if this is the case)
    rename(lower.CL = any_of("asymp.LCL"),
           upper.CL = any_of("asymp.UCL"))
  
  
  # Join and format CP comparisons
  all_comparisons <- status_comparisons %>%
    full_join(season_comparisons) %>%
    full_join(interaction_comparisons) %>%
    rename(group_letter = .group) %>%
    mutate(
      comparison = factor(comparison, levels = c("status", "season", "interaction")),
      season = factor(season, levels = c("S1", "S2", "S3", "S4")),
      status = factor(status, levels = c("Altered", "Preserved", "Restored")),
      group_letter = gsub(" ", "", group_letter),
      emmean_gwp = inv_logsign(emmean),
      SE_gwp = inv_logsign(SE),
      lower.CL_gwp = inv_logsign(lower.CL),
      upper.CL_gwp = inv_logsign(upper.CL)
    ) %>%
      dplyr::select(any_of(c(
      "ghgspecies", "casepilot", "comparison", "status", "season", "df",
      "emmean", "SE", "lower.CL", "upper.CL", "group_letter",
      "emmean_gwp", "SE_gwp", "lower.CL_gwp", "upper.CL_gwp"
    ))) %>%
    arrange(ghgspecies, casepilot, comparison, status, season)
  
  # Store casepilot results
  comparisons_list[[casepilot_name]] <- all_comparisons
}

# Unir todos los resultados en un solo data.frame
co2_casepilot_comparisons <- bind_rows(comparisons_list)


#Save all comparisons:  
write.csv(x = co2_casepilot_comparisons, file = paste0(plots_path,"/CO2models_emmeans_posthoc.csv"),row.names = F)


#Plot casepilot-per-casepilot, emmeans +- SE/CL (in gwp scale)

##Plots Status effect-----

#Camargue: 
plot_co2_status_CA <- co2_casepilot_comparisons %>%
  filter(casepilot == "CA", comparison == "status") %>%
  ggplot(., aes(x = status, y = emmean_gwp, label = group_letter)) +
      geom_errorbar(aes(ymin = emmean_gwp - SE_gwp, ymax = emmean_gwp + SE_gwp),
                    linewidth = 0.2, width = 0.2) +
      geom_point(shape = 15, size = 4, aes(col=status)) +
      geom_text(nudge_x = 0.2, color = "black") +
      geom_hline(yintercept=0, linetype="dashed")+
      theme_bw() +
      theme(axis.title = element_text(face = "bold"),
            axis.text = element_text(face = "bold"),
            plot.caption = element_text(hjust = 0)) +
      labs(title = expression(Status~effect~on~net~CO[2]~flux),
           subtitle = paste0("Camargue"), 
           y=expression(EMM~(g~CO[2~eq]~d^-1~m^-2)),
           x="Conservation status",
           col=paste0("Status"))

plot_co2_status_CA
ggsave(plot=plot_co2_status_CA, filename="CO2_status_plot_CA.png",
       device = "png",path = plots_path,dpi = 400,
       width = 115, height = 100,units = "mm")


#Curonian Lagoon: 
plot_co2_status_CU<- co2_casepilot_comparisons %>%
  filter(casepilot == "CU", comparison == "status") %>%
  ggplot(., aes(x = status, y = emmean_gwp, label = group_letter)) +
      geom_errorbar(aes(ymin = emmean_gwp - SE_gwp, ymax = emmean_gwp + SE_gwp),
                    linewidth = 0.2, width = 0.2) +
      geom_point(shape = 15, size = 4, aes(col=status)) +
      geom_text(nudge_x = 0.2, color = "black") +
      geom_hline(yintercept=0, linetype="dashed")+
      theme_bw() +
      theme(axis.title = element_text(face = "bold"),
            axis.text = element_text(face = "bold"),
            plot.caption = element_text(hjust = 0)) +
      labs(title = expression(Status~effect~on~net~CO[2]~flux),
           subtitle = paste0("Curonian lagoon"), 
           y=expression(EMM~(g~CO[2~eq]~d^-1~m^-2)),
           x="Conservation status",
           col=paste0("Status"))

plot_co2_status_CU
ggsave(plot=plot_co2_status_CU, filename="CO2_status_plot_CU.png",
       device = "png",path = plots_path,dpi = 400,
       width = 115, height = 100,units = "mm")

#Danube delta: 
plot_co2_status_DA<- co2_casepilot_comparisons %>%
  filter(casepilot == "DA", comparison == "status") %>%
  ggplot(., aes(x = status, y = emmean_gwp, label = group_letter)) +
      geom_errorbar(aes(ymin = emmean_gwp - SE_gwp, ymax = emmean_gwp + SE_gwp),
                    linewidth = 0.2, width = 0.2) +
      geom_point(shape = 15, size = 4, aes(col=status)) +
      geom_text(nudge_x = 0.2, color = "black") +
      geom_hline(yintercept=0, linetype="dashed")+
      theme_bw() +
      theme(axis.title = element_text(face = "bold"),
            axis.text = element_text(face = "bold"),
            plot.caption = element_text(hjust = 0)) +
      labs(title = expression(Status~effect~on~net~CO[2]~flux),
           subtitle = paste0("Danube delta"),
           y=expression(EMM~(g~CO[2~eq]~d^-1~m^-2)),
           x="Conservation status",
           col=paste0("Status"))

plot_co2_status_DA
ggsave(plot=plot_co2_status_DA, filename="CO2_status_plot_DA.png",
       device = "png",path = plots_path,dpi = 400,
       width = 115, height = 100,units = "mm")

#Dutch delta: 
plot_co2_status_DU<- co2_casepilot_comparisons %>%
  filter(casepilot == "DU", comparison == "status") %>%
  ggplot(., aes(x = status, y = emmean_gwp, label = group_letter)) +
      geom_errorbar(aes(ymin = emmean_gwp - SE_gwp, ymax = emmean_gwp + SE_gwp),
                    linewidth = 0.2, width = 0.2) +
      geom_point(shape = 15, size = 4, aes(col=status)) +
      geom_text(nudge_x = 0.2, color = "black") +
      geom_hline(yintercept=0, linetype="dashed")+
      theme_bw() +
      theme(axis.title = element_text(face = "bold"),
            axis.text = element_text(face = "bold"),
            plot.caption = element_text(hjust = 0)) +
      labs(title = expression(Status~effect~on~net~CO[2]~flux),
           subtitle = paste0("Dutch delta"),
           y=expression(EMM~(g~CO[2~eq]~d^-1~m^-2)),
           x="Conservation status",
           col=paste0("Status"))

plot_co2_status_DU
ggsave(plot=plot_co2_status_DU, filename="CO2_status_plot_DU.png",
       device = "png",path = plots_path,dpi = 400,
       width = 115, height = 100,units = "mm")


#Ria d'Aveiro: 
plot_co2_status_RI<- co2_casepilot_comparisons %>%
  filter(casepilot == "RI", comparison == "status") %>%
  ggplot(., aes(x = status, y = emmean_gwp, label = group_letter)) +
      geom_errorbar(aes(ymin = emmean_gwp - SE_gwp, ymax = emmean_gwp + SE_gwp),
                    linewidth = 0.2, width = 0.2) +
      geom_point(shape = 15, size = 4, aes(col=status)) +
      geom_text(nudge_x = 0.2, color = "black") +
      geom_hline(yintercept=0, linetype="dashed")+
      theme_bw() +
      theme(axis.title = element_text(face = "bold"),
            axis.text = element_text(face = "bold"),
            plot.caption = element_text(hjust = 0)) +
      labs(title = expression(Status~effect~on~net~CO[2]~flux),
           subtitle = paste0("Ria d'Aveiro"), 
           y=expression(EMM~(g~CO[2~eq]~d^-1~m^-2)),
           x="Conservation status",
           col=paste0("Status"))

plot_co2_status_RI
ggsave(plot=plot_co2_status_RI, filename="CO2_status_plot_RI.png",
       device = "png",path = plots_path,dpi = 400,
       width = 115, height = 100,units = "mm")

#Valencian wetland: 
plot_co2_status_VA<- co2_casepilot_comparisons %>%
  filter(casepilot == "VA", comparison == "status") %>%
  ggplot(., aes(x = status, y = emmean_gwp, label = group_letter)) +
      geom_errorbar(aes(ymin = emmean_gwp - SE_gwp, ymax = emmean_gwp + SE_gwp),
                    linewidth = 0.2, width = 0.2) +
      geom_point(shape = 15, size = 4, aes(col=status)) +
      geom_text(nudge_x = 0.2, color = "black") +
      geom_hline(yintercept=0, linetype="dashed")+
      theme_bw() +
      theme(axis.title = element_text(face = "bold"),
            axis.text = element_text(face = "bold"),
            plot.caption = element_text(hjust = 0)) +
      labs(title = expression(Status~effect~on~net~CO[2]~flux),
           subtitle = paste0("Valencian wetland"), 
           y=expression(EMM~(g~CO[2~eq]~d^-1~m^-2)),
           x="Conservation status",
           col=paste0("Status"))

plot_co2_status_VA
ggsave(plot=plot_co2_status_VA, filename="CO2_status_plot_VA.png",
       device = "png",path = plots_path,dpi = 400,
       width = 115, height = 100,units = "mm")


##Plots Seasonal effect -----
#To do, little interest for overview ppt (seasonal effect is averaged across the three status)

#Example:
#Camargue: 
plot_co2_season_CA<- co2_casepilot_comparisons %>%
  filter(casepilot == "CA", comparison == "season") %>%
  ggplot(., aes(x = season, y = emmean_gwp, label = group_letter)) +
      geom_errorbar(aes(ymin = emmean_gwp - SE_gwp, ymax = emmean_gwp + SE_gwp),
                    linewidth = 0.2, width = 0.2) +
      geom_point(shape = 15, size = 4, aes(col=season)) +
      geom_text(nudge_x = 0.2, color = "black") +
      geom_hline(yintercept=0, linetype="dashed")+
      theme_bw() +
      theme(axis.title = element_text(face = "bold"),
            axis.text = element_text(face = "bold"),
            plot.caption = element_text(hjust = 0)) +
      labs(title = expression(Seasonal~effect~on~net~CO[2]~flux),
           subtitle = paste0("Camargue"), 
           caption = paste0("Seasonal effect is averaged across the three conservation status.\nBoxes indicate the EM mean, error bars indicate SE,\n",              "Letters indicate different groups"),
           y=expression(EMM~(g~CO[2~eq]~d^-1~m^-2)),
           x="Sampling campaign",
           col=paste0("Sampling"))

plot_co2_season_CA


##Plots Interaction effect ----
#To do, unclear interest for overview ppt (maybe interesting for some casepilots)

#Example:
#Camargue: interaction plot
plot_co2_interaction_CA<- co2_casepilot_comparisons %>%
  filter(casepilot == "CA", comparison == "interaction") %>%
  #ggplot inside brackets to be able to access casepilot name with ggtitle or define dodge for all
  { 
    #Define dodge width for separating status.
    dodge <- position_dodge(width = 0.6)
    
    ggplot(., aes(x = season, y = emmean_gwp, label = group_letter, group=status)) +
      geom_errorbar(aes(ymin = emmean_gwp - SE_gwp, ymax = emmean_gwp + SE_gwp),
                    linewidth = 0.2, width = 0.2, position = dodge) +
      geom_point(shape = 15, size = 4, aes(col=status), position = dodge) +
      geom_text(vjust = -2, color = "black", position = dodge) +
      geom_hline(yintercept=0, linetype="dashed")+
      theme_bw() +
      theme(axis.title = element_text(face = "bold"),
            axis.text = element_text(face = "bold"),
            plot.caption = element_text(hjust = 0)) +
      labs(title = expression(Seasonal~effect~on~net~CO[2]~flux),
           subtitle = paste0("Camargue"), 
           caption = paste0("Boxes indicate the EM mean, error bars indicate SE,\n",              "Letters indicate different groups"),
           y=expression(EMM~(g~CO[2~eq]~d^-1~m^-2)),
           x="Sampling campaign",
           col=paste0("Sampling"))
  }

plot_co2_interaction_CA





#______________________-----

#CH4_______ ------

#FOR CH4, we need to use gamma distributions (right-skewed nature of data), therefore we need to shift the data to positive-values only. WE MUST NOT do shift-log ourselves (REDUNDANT, the "log-transformation" is already accounted for when using family = Gamma(link = "log") in glmmTMB).
#This option assumes: skeweed positive errors and that variance increases with mean (no homocedasticity)

#We can try a gaussian family with dailyflux_logsign transformation (assumes normal residuals and homocedasticity).  

#As we need to keep the shift stored for each case-pilot, we will subset and modify the data4models dataframe.
ch4_formodels<-  data4models %>% 
  filter(ghgspecies=="gwp_ch4") %>% 
  filter(!is.na(dailyflux)) 


#HOW TO SELECT THE BEST SHIFT MAGNITUDE FOR EACH DATASET?
#Is the shift appropriate in magnitude? shift of abs(min)*1.1 causes "broken distributions", with the minimum value being very far from the rest on the log-scale, doubt potential issues? 

#SOLUTION CHOSEN: semi-manual selection for each casepilot: trying to minimize skewness while miniminzing gaps in distribution.
#Selection of shift to: 
  # 1. reduce where possible the skeweness of the transformed data (logscale) and 
  # 2. avoid large gaps in distribution (minimize gap from minimum to percentil 5)

library(e1071)  
library(ggplot2)
library(scales)
library(gridExtra)

# Get casepilot dailyflux
data <- ch4_formodels %>% filter(casepilot=="VA") %>% pull(dailyflux)

#Bound the lower-limit of epsilon (at least 10% higher than minimum value)
minepsilon<-log10(abs(min(data)*1.1))

# Define a range of epsilon values 
epsilons <- 10^seq(minepsilon, -1, length.out = 50)

# Compute skewness and gap for each epsilon
results <- data.frame(
  epsilon = epsilons,
  skewness = NA,
  abs_skew = NA,
  gap = NA
)

for (i in seq_along(epsilons)) {
  eps <- epsilons[i]
  log_data <- log(data + eps)
  results$skewness[i] <- skewness(log_data)
  results$abs_skew[i] <- abs(results$skewness[i])
  results$gap[i] <- quantile(log_data, 0.01, na.rm = T) - min(log_data)
}

# Step 1: Find minimum absolute skewness
min_abs_skew <- min(results$abs_skew, na.rm = T)

# Step 2: Filter epsilons based on skewness (variable criteria)
candidate_eps <- subset(results, abs_skew <= 0.5) #for low-skewed distributions (focus on gap)
candidate_eps <- subset(results, abs_skew <= 1.05 * min_abs_skew) #for high-skewed distributions (balance gap and skewness)

# Step 3: Select epsilon with smallest gap
optimal_row <- candidate_eps[which.min(candidate_eps$gap), ]
optimal_epsilon <- optimal_row$epsilon
{
  # Plot skewness vs epsilon
  p1 <- ggplot(results, aes(x = epsilon, y = skewness)) +
    geom_line() +
    geom_point() +
    geom_point(data=optimal_row, col="red")+
    geom_vline(xintercept = optimal_epsilon, linetype = "dashed", color = "red") +
    scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
    labs(title = "Skewness vs Shift ε", x = "Shift ε (log scale)", y = "Skewness")
  
  #Plot skewness vs gap (minimum-percentil 1)
  p2 <- ggplot(results, aes(x = gap, y = skewness)) +
    geom_line() +
    geom_point() +
    geom_point(data=optimal_row, col="red")+
    labs(title = "Skewness vs gap", x = "Gap (min to percentil1)", y = "Skewness")
  
  
  # Plot histogram of log-transformed data with optimal epsilon
  log_data_opt <- log(data + optimal_epsilon)
  df_log <- data.frame(log_y = log_data_opt)
  
  p3<- ggplot(df_log, aes(x = log_y)) +
    geom_histogram(aes(y = ..density..), bins = 30, fill = "skyblue", color = "black") +
    geom_density(color = "red", linewidth = 1) +
    labs(title = paste0("Log-Transformed Data (ε = ", signif(optimal_epsilon, 2), ")"),
         x = "log(y + ε)", y = "Density")
  
  # Display plots
  grid.arrange(p1, p2,p3, ncol = 3)
}

#shift for casepilots: 
#CA: 0.058 #MInimum gap within skewness of 10% of minimum. 
#CU: 0.0047 #Minimum gap within skewness < 0.5
#DA: 0.033 #MInimum gap within skewness  <0.5 
#DU: 0.067 #Mnimum gap within skewness of 1% of minimum
#RI: 0.081 #Chosen to minimize skewness (gaps unavoidable and similar for positive and negative "outliers")
#VA: 0.096 #Minimum gap within skewness of 5% of minimum. 


#store all tranformations in data-frame to allow comparing distributions. 
ch4_formodels<-  data4models %>% 
  filter(ghgspecies=="gwp_ch4") %>% 
  filter(!is.na(dailyflux)) %>% 
  group_by(casepilot) %>% 
  mutate(shift_casepilot=case_when(casepilot=="CA"~0.058,
                                   casepilot=="CU"~0.0047,
                                   casepilot=="DA"~0.033,
                                   casepilot=="DU"~0.067,
                                   casepilot=="RI"~0.081,
                                   casepilot=="VA"~0.096)) %>% #Shifts selected for each casepilot
  ungroup() %>% 
  mutate(dailyflux_logsign=logsign(dailyflux),
         dailyflux_shift= dailyflux + shift_casepilot)
  

#Store shifts in dataframe independently: 
casepilot_ch4_shifts<- ch4_formodels %>% 
  group_by(casepilot) %>% 
  summarise(shift=mean(shift_casepilot))


#dailyflux_logsign (gaussian) vs dailyflux_shift (gamma)
#I will compare the two options and chose the one with lowest AIC (assuming residuals are good for both)


#As a visual guide for distributions: 
ch4_formodels %>% 
  filter(casepilot=="CA") %>% 
  mutate(dailyflux_shiftlog=log(dailyflux_shift)) %>% 
  pivot_longer(cols = c(dailyflux,dailyflux_logsign, dailyflux_shiftlog), names_to = "transformation", values_to = "value") %>% 
ggplot(aes(x=value))+
  geom_histogram(aes(y = ..density..), bins = 30, fill = "skyblue", color = "black") +
  geom_density(color = "red", linewidth = 1) +
  facet_grid(casepilot~transformation, scales = "free") +
  theme_minimal()

#CASEPILOT models-------
#STEPS: 
#1. select best transformation (using lmer, no sampleweights)
#2. Optimize model without sampleweights (lmer as first step --> lme (normal but heterocedasticity) or glmm (non-normal, heterocedasitcity or not))
#3. Optimize model with sampleweights(glmm)


#CH4 distribution notes: we need a gamma distribution (highly skewed data distribution, non-symetric), for that we require positive values always.
#For CH4, it makes sense to transform the data using a shift-log approach ( log(value+shift)), and then use a gamma family distribution. 

# gaussian distributions + logsign transformation might be appropriate only for some cases (those were we do not expect skeweness, potentially aveiro) for most cases this approach distorts the skeweness of the CH4 data, which is a fundamental property of this type of data, NOT APPROPRIATE. 



#1. Best common transformation ------ 
#CH4 data is naturally positively skewed. We need to retain this pattern. 

#Transformations have been selected already, either shifted (with gamma distribution) or try dailyflux_logsign with gaussian (or t_family)

#DEFAULT IS shifted tranformation, we expect non-normal distributions, so we should directly use gamma families (only supported by GLMM): So NO LMER check (it does not tell us anything of use).

#Check LMER with signlog to check whether we achieve normal residuals in any case. If no normality of residuals, default to glmm with gamma familiy and shift transformation. 

#Focus on normality of residuals, homocedasticity of residuals, model fit, for multiple transformations.

# Initialize summary table
ch4_lmer_results <- tibble()

# Loop over all combinations of casepilot and CH4
for (cp in unique(data4models$casepilot)) {
  for (ghg in c("gwp_ch4")) {
    
    cat("\n\n==================== CASEPILOT:", cp, "| GHG:", ghg, "====================\n")
    
    data_sub <- data4models %>% filter(grepl(cp, casepilot)& ghgspecies == ghg) %>% 
      filter(!is.na(dailyflux))
    
    #Transformations
    data_sub <- data_sub %>%
      mutate(
        dailyflux_logsign = logsign(dailyflux),
        dailyflux_rootsign = rootsign(dailyflux),
        dailyflux_arcsinhsign = arcsinhsign(dailyflux)
      )
    
    # Step 3: Fit lmer models
    model_untrans <- tryCatch(lmer(dailyflux ~ status * season + (1 | subsite), data = data_sub), error = function(e) NULL)
    model_logsign <- tryCatch(lmer(dailyflux_logsign ~ status * season + (1 | subsite), data = data_sub), error = function(e) NULL)
    model_rootsign <- tryCatch(lmer(dailyflux_rootsign ~ status * season + (1 | subsite), data = data_sub), error = function(e) NULL)
    model_arcsinhsign<- tryCatch(lmer(dailyflux_arcsinhsign ~ status * season + (1 | subsite), data = data_sub), error = function(e) NULL)
    
    
    # Step 4: Diagnostics
    check_model <- function(model, name) {
      if (is.null(model)) return(NULL)
      res <- residuals(model)
      shapiro_res <- shapiro.test(res)$p.value
      levene_status <- tryCatch(leveneTest(res ~ data_sub$status)$"Pr(>F)"[1], error = function(e) NA)
      levene_season <- tryCatch(leveneTest(res ~ data_sub$season)$"Pr(>F)"[1], error = function(e) NA)
      r2 <- tryCatch(performance::r2_nakagawa(model), error = function(e) NULL)
      aic_val <- tryCatch(AIC(model), error = function(e) NA)
      bic_val <- tryCatch(BIC(model), error = function(e) NA)
      rmse_cv <- tryCatch({
        ctrl <- trainControl(method = "cv", number = 10)
        cv_model <- train(
          x = model.matrix(model)[, -1],
          y = getME(model, "y"),
          method = "lm",
          trControl = ctrl
        )
        cv_model$results$RMSE
      }, error = function(e) NA)
      
      tibble(
        casepilot = cp,
        gas = ghg,
        model = name,
        shapiro_resid_p = shapiro_res,
        levene_status_p = levene_status,
        levene_season_p = levene_season,
        marginal_r2 = if (!is.null(r2)) r2$R2_marginal else NA,
        conditional_r2 = if (!is.null(r2)) r2$R2_conditional else NA,
        AIC = aic_val,
        BIC = bic_val,
        RMSE_CV = rmse_cv
      )
    }
    
    res_untrans <- check_model(model_untrans, "Untransformed")
    res_logsign <- check_model(model_logsign, "Signed Log")
    res_rootsign <- check_model(model_rootsign, "Signed Root")
    res_arcsinhsign <- check_model(model_arcsinhsign, "Signed Arcsinh")
    
    ch4_lmer_results <- bind_rows(ch4_lmer_results, res_untrans, res_logsign, res_rootsign, res_arcsinhsign)
  }
}

ch4_lmer_results

#Normality of residuals never achieved
ch4_lmer_results %>% 
  filter(shapiro_resid_p>0.05)

#Homocedasticity for both factors only achieved for VA 
ch4_lmer_results %>% 
  filter(levene_season_p>0.05&levene_status_p>0.05)

#Best transformation is consistently Signed Log:
ch4_lmer_results %>% 
  group_by(casepilot,gas) %>% 
  filter(AIC==min(AIC))


#CONCLUSION: Lmer is not appropriate for modeling casepilots individually as is. We need to account for heterocedasticity and non-normality. We have non-normal residuals--> glmm for all (allows different family distributions and heterocedasticity)  


#2. Optimice GLMM for each casepilot-------

#ALL 6 casepilots show strong non-normality... most show heterocedasticity for season.

#We will compare logsign-gaussian/t_family models with shift-gamma models, potentially enabling heterocedasticity via dispformula. We will compare residual behaviours (via plots) and model fit (via AIC). 


#Notes on singularity: to have a model with singularity means that some of the random effect terms likely have variance estimates of (almost) zero, or there is redundancy in the random structure. IN our case, the random effect is solely the subsite (to account for subsite-specific differences and repeated samplings across seasons). Our random effect structure is the most simple possible (one random effect for intercept, not for slope), and we expect this random effect to account for variance when the 2 subsites of each status are fundamentally different. With this structure, we should not be too concerned about singuarity. The inclusion of the subsite as a random effect might not be statistically necessary (variance~0) but it is still conceptually justified.  



## 2.1. CA model ----
#CA: lmer showed non-normal residuals, heterocedasticity for season (marignal for status)
ch4_lmer_results %>% filter(casepilot=="CA")

#approach: compare fit between gausian, t_family models (logsign data) and gamma (shifted data), without heterocedasticity, with for season only and with for season and status. 

ca_ch4<- ch4_formodels %>% filter(casepilot=="CA") %>% 
  filter(dailyflux>-0.025) #Remove 3 outliers 
  

ca_ch4 %>% 
  filter(casepilot=="CA") %>% 
  mutate(dailyflux_shiftlog=log(dailyflux_shift)) %>% 
  pivot_longer(cols = c(dailyflux_logsign, dailyflux_shiftlog), names_to = "transformation", values_to = "value") %>% 
  ggplot(aes(x=value))+
  geom_histogram(aes(y = ..density..), bins = 50, fill = "skyblue", color = "black") +
  geom_density(color = "red", linewidth = 1) +
  facet_grid(casepilot~transformation, scales = "free") +
  theme_minimal()


ca_glmm_gaus0<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite), 
                        data = ca_ch4,
                        family = gaussian(),
                        dispformula = ~1)

ca_glmm_stu0<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite), 
                       data = ca_ch4,
                       family = t_family,
                       dispformula = ~1) 

ca_glmm_gamm0<- glmmTMB(formula = dailyflux_shift~status*season + (1|subsite), 
                        data = ca_ch4,
                        family = Gamma(link = "log"),
                        dispformula = ~1)


resgaus0<- simulateResiduals(ca_glmm_gaus0, n=1000)
restu0<- simulateResiduals(ca_glmm_stu0, n=1000)
resgamm0<- simulateResiduals(ca_glmm_gamm0, n=1000)
  
#Check gaussian model
plot(resgaus0) #check normality, 
plotResiduals(resgaus0, ca_ch4$season) #check heterocedasticity
plotResiduals(resgaus0, ca_ch4$status) #check heterocedasticity
#All asummptions fail, non-normality, funnel shape, deviation from Q-Q, heterocedasticity for both season and status

#Check t_family model
plot(restu0) #check normality, 
plotResiduals(restu0, ca_ch4$season) #check heterocedasticity
plotResiduals(restu0, ca_ch4$status) #check heterocedasticity
#MODEL does not converge and all asumptions fail

#Check gamma model
plotQQunif(resgamm0) #Check QQplot
plotResiduals(resgamm0) #residuals vs predicted
plot(resgamm0) #check Q-Q, residuals vs fitted, 
plotResiduals(resgamm0, ca_ch4$season) #check heterocedasticity
plotResiduals(resgamm0, ca_ch4$status) #check heterocedasticity
testDispersion(resgamm0) #checks wether residual variance is higher or lower than expected, n.s.--> no issue. 
testOutliers(resgamm0) #Checks if there are more extreme residuals than expected

#Gamma looks best, lets try with heteorcedasticity at season or season*status

ca_glmm_gamm1<- glmmTMB(formula = dailyflux_shift~status*season + (1|subsite), 
                        data = ca_ch4,
                        family = Gamma(link = "log"),
                        dispformula = ~season)

ca_glmm_gamm2<- glmmTMB(formula = dailyflux_shift~status*season + (1|subsite), 
                        data = ca_ch4,
                        family = Gamma(link = "log"),
                        dispformula = ~season + status)

resgamm1<- simulateResiduals(ca_glmm_gamm1, n=1000)
resgamm2<- simulateResiduals(ca_glmm_gamm2, n=1000)

plot(resgamm1) #check Q-Q, residuals vs fitted, 
plotResiduals(resgamm1, ca_ch4$season) #check heterocedasticity
plotResiduals(resgamm1, ca_ch4$status) #check heterocedasticity
testDispersion(resgamm1) #checks wether residual variance is higher or lower than expected, n.s.--> no issue. 
testOutliers(resgamm1) #Checks if there are more extreme residuals than expected

plot(resgamm2) #check normality, 
plotResiduals(resgamm2, ca_ch4$season) #check heterocedasticity
plotResiduals(resgamm2, ca_ch4$status) #check heterocedasticity
testDispersion(resgamm2) #checks wether residual variance is higher or lower than expected, n.s.--> no issue. 
testOutliers(resgamm2) #Checks if there are more extreme residuals than expected


#ANOVA to decide: (across families we can only compare AIC,BIC. rest must be ignored)
anova(ca_glmm_gaus0,ca_glmm_gamm0,ca_glmm_gamm1,ca_glmm_gamm2)
#Gamma distribution is is consistently better 
anova(ca_glmm_gamm0,ca_glmm_gamm1,ca_glmm_gamm2)
#Big improvemnt with season heterocedasticity, limited improvement with additional status heterocedasticity. 
anova(ca_glmm_gamm1,ca_glmm_gamm2)


#CA best model: gamma, allowing heterocedasticity across seasons and status (expected, altered will have more dispersion due to higher impact of ebullition). 
ca_best_ch4<- ca_glmm_gamm2
Anova(ca_best_ch4,type = 3) 

#Cannot calculate R2 (conditional & marginal) for models with heterocedasticity,
#Get PseudoR2
pseudoR2_ca_ch4<-pseudoR2_from_null(ca_best_ch4, ca_ch4)
pseudoR2_ca_ch4

#R2 of homocedastic model
ca_best_homoc_ch4<- ca_glmm_gamm0
r2(ca_best_homoc_ch4)
MuMIn::r.squaredGLMM(ca_best_homoc_ch4) 

#Singularity?
performance::check_singularity(ca_best_ch4)
#NO. 

#Very high R2, but presence of residual outliers....do we worry and try to remove them?

#Remove non-best models
rm(ca_glmm_gaus0,ca_glmm_gamm0,ca_glmm_gamm1,ca_glmm_gamm2)


## 2.2. CU model ----
#CU: lmer showed non-normal residuals, heterocedasticity for season and status
ch4_lmer_results %>% filter(casepilot=="CU")
#approach: compare fit between gausian, t_family models (logsign data) and gamma (shifted data), without heterocedasticity, with for season only and with for season and status. 


cu_ch4<- ch4_formodels %>% filter(casepilot=="CU")


cu_ch4 %>% 
  filter(casepilot=="CU") %>% 
  mutate(dailyflux_shiftlog=log(dailyflux_shift)) %>% 
  pivot_longer(cols = c(dailyflux_logsign, dailyflux_shiftlog), names_to = "transformation", values_to = "value") %>% 
  ggplot(aes(x=value))+
  geom_histogram(aes(y = ..density..), bins = 50, fill = "skyblue", color = "black") +
  geom_density(color = "red", linewidth = 1) +
  facet_grid(casepilot~transformation, scales = "free") +
  theme_minimal()


cu_glmm_gaus0<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite), 
                        data = cu_ch4,
                        family = gaussian(),
                        dispformula = ~1)

cu_glmm_stu0<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite), 
                       data = cu_ch4,
                       family = t_family,
                       dispformula = ~1) 

cu_glmm_gamm0<- glmmTMB(formula = dailyflux_shift~status*season + (1|subsite), 
                        data = cu_ch4,
                        family = Gamma(link = "log"),
                        dispformula = ~1)


resgaus0<- simulateResiduals(cu_glmm_gaus0, n=1000)
restu0<- simulateResiduals(cu_glmm_stu0, n=1000)
resgamm0<- simulateResiduals(cu_glmm_gamm0, n=1000)

#Check gaussian model
plot(resgaus0) #check normality, 
plotResiduals(resgaus0, cu_ch4$season) #check heterocedasticity
plotResiduals(resgaus0, cu_ch4$status) #check heterocedasticity
#All asummptions fail, non-normality, funnel shape, deviation from Q-Q, heterocedasticity for both season and status

#Check t_family model
plot(restu0) #check normality, 
plotResiduals(restu0, cu_ch4$season) #check heterocedasticity
plotResiduals(restu0, cu_ch4$status) #check heterocedasticity
#All asumptions fail

#Check gamma model
summary(cu_glmm_gamm0)
plot(resgamm0) #check Q-Q, residuals vs fitted, 
plotResiduals(resgamm0, cu_ch4$season) #check heterocedasticity
plotResiduals(resgamm0, cu_ch4$status) #check heterocedasticity
testDispersion(resgamm0) #checks wether residual variance is higher or lower than expected, n.s.--> no issue. 
testOutliers(resgamm0) #Checks if there are more extreme residuals than expected

#Gamma looks best, lets try with heteorcedasticity at season or season*status

cu_glmm_gamm1<- glmmTMB(formula = dailyflux_shift~status*season + (1|subsite), 
                        data = cu_ch4,
                        family = Gamma(link = "log"),
                        dispformula = ~season)

cu_glmm_gamm2<- glmmTMB(formula = dailyflux_shift~status*season + (1|subsite), 
                        data = cu_ch4,
                        family = Gamma(link = "log"),
                        dispformula = ~season + status)

resgamm1<- simulateResiduals(cu_glmm_gamm1, n=1000)
resgamm2<- simulateResiduals(cu_glmm_gamm2, n=1000)

plot(resgamm1) #check Q-Q, residuals vs fitted, 
plotResiduals(resgamm1, cu_ch4$season) #check heterocedasticity
plotResiduals(resgamm1, cu_ch4$status) #check heterocedasticity
testDispersion(resgamm1) #checks wether residual variance is higher or lower than expected, n.s.--> no issue. 
testOutliers(resgamm1) #Checks if there are more extreme residuals than expected

plot(resgamm2) #check normality, 
plotQQunif(resgamm2) #Check QQplot
plotResiduals(resgamm2, cu_ch4$season) #check heterocedasticity
plotResiduals(resgamm2, cu_ch4$status) #check heterocedasticity
testDispersion(resgamm2) #checks wether residual variance is higher or lower than expected, n.s.--> no issue. 
testOutliers(resgamm2) #Checks if there are more extreme residuals than expected


#GAUssian and student distribution models do not meet assumptions

#Gamma distribution is is consistently better 
anova(cu_glmm_gamm0,cu_glmm_gamm1,cu_glmm_gamm2)
#Limited improvemnt with season heterocedasticity, limited improvement with additional status heterocedasticity. 


#CU best model: gamma, allowing heterocedasticity across seasons and status (expected, altered will have more dispersion due to higher impact of ebullition). 
cu_best_ch4<- cu_glmm_gamm2
Anova(cu_best_ch4,type = 3) 

#Cannot calculate R2 (conditional & marginal) for models with heterocedasticity,
#Get PseudoR2
pseudoR2_cu_ch4<-pseudoR2_from_null(cu_best_ch4, cu_ch4)
pseudoR2_cu_ch4

#R2 of homocedastic model
cu_best_homoc_ch4<- cu_glmm_gamm0
r2(cu_best_homoc_ch4)
MuMIn::r.squaredGLMM(cu_best_homoc_ch4) 

#Singularity?
performance::check_singularity(cu_best_ch4)
#NO. 

#High R2, but presence of residual outliers....do we worry and try to remove them?
#AIC of gaussian and student models is better but they fail their assumptions. 

#Remove non-best models
rm(cu_glmm_gaus0,cu_glmm_stu0,cu_glmm_gamm0,cu_glmm_gamm1,cu_glmm_gamm2)

## 2.3. DA model ----
#DA: lmer showed non-normal residuals, heterocedasticity for season,
ch4_lmer_results %>% filter(casepilot=="DA")

#approach: compare fit between gausian, t_family models (logsign data) and gamma (shifted data), without heterocedasticity, with for season only and with for season and status. 

ch4_formodels %>% 
  filter(casepilot=="DA") %>% 
  mutate(dailyflux_shiftlog=log(dailyflux_shift)) %>% 
  pivot_longer(cols = c(dailyflux,dailyflux_logsign, dailyflux_shiftlog), names_to = "transformation", values_to = "value") %>% 
  ggplot(aes(x=value))+
  geom_histogram(aes(y = ..density..), bins = 50, fill = "skyblue", color = "black") +
  geom_density(color = "red", linewidth = 1) +
  facet_grid(casepilot~transformation, scales = "free") +
  theme_minimal()



da_ch4<- ch4_formodels %>% filter(casepilot=="DA")

da_glmm_gaus0<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite), 
                        data = da_ch4,
                        family = gaussian(),
                        dispformula = ~1)

da_glmm_stu0<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite), 
                       data = da_ch4,
                       family = t_family,
                       dispformula = ~1) 

da_glmm_gamm0<- glmmTMB(formula = dailyflux_shift~status*season + (1|subsite), 
                        data = da_ch4,
                        family = Gamma(link = "log"),
                        dispformula = ~1)


resgaus0<- simulateResiduals(da_glmm_gaus0, n=1000)
restu0<- simulateResiduals(da_glmm_stu0, n=1000)
resgamm0<- simulateResiduals(da_glmm_gamm0, n=1000)

#Check gaussian model
plot(resgaus0) #check normality, 
plotResiduals(resgaus0, da_ch4$season) #check heterocedasticity
plotResiduals(resgaus0, da_ch4$status) #check heterocedasticity
#All asummptions fail, non-normality, funnel shape, deviation from Q-Q, heterocedasticity for both season and status

#Check t_family model
plot(restu0) #check normality, 
plotResiduals(restu0, da_ch4$season) #check heterocedasticity
plotResiduals(restu0, da_ch4$status) #check heterocedasticity
#All asumptions fail

#Check gamma model
plot(resgamm0) #check Q-Q, residuals vs fitted, 
plotResiduals(resgamm0, da_ch4$season) #check heterocedasticity
plotResiduals(resgamm0, da_ch4$status) #check heterocedasticity
testDispersion(resgamm0) #checks wether residual variance is higher or lower than expected, n.s.--> no issue. 
testOutliers(resgamm0) #Checks if there are more extreme residuals than expected

#Gamma looks best, lets try with heteorcedasticity at season or season*status

da_glmm_gamm1<- glmmTMB(formula = dailyflux_shift~status*season + (1|subsite), 
                        data = da_ch4,
                        family = Gamma(link = "log"),
                        dispformula = ~season)

da_glmm_gamm2<- glmmTMB(formula = dailyflux_shift~status*season + (1|subsite), 
                        data = da_ch4,
                        family = Gamma(link = "log"),
                        dispformula = ~season + status)

resgamm1<- simulateResiduals(da_glmm_gamm1, n=1000)
resgamm2<- simulateResiduals(da_glmm_gamm2, n=1000)

plot(resgamm1) #check Q-Q, residuals vs fitted, 
plotResiduals(resgamm1, da_ch4$season) #check heterocedasticity
plotResiduals(resgamm1, da_ch4$status) #check heterocedasticity
testDispersion(resgamm1) #checks wether residual variance is higher or lower than expected, n.s.--> no issue. 
testOutliers(resgamm1) #Checks if there are more extreme residuals than expected

plot(resgamm2) #check normality, 
plotQQunif(resgamm2) #Check QQplot
plotResiduals(resgamm2, da_ch4$season) #check heterocedasticity
plotResiduals(resgamm2, da_ch4$status) #check heterocedasticity
testDispersion(resgamm2) #checks wether residual variance is higher or lower than expected, n.s.--> no issue. 
testOutliers(resgamm2) #Checks if there are more extreme residuals than expected


#GAUssian and student distribution models do not meet assumptions. Student model with better AIC but weird pattern in residuals. 
anova(da_glmm_gaus0, da_glmm_stu0, da_glmm_gamm0)

#Gamma distribution is is consistently better 
anova(da_glmm_gamm0,da_glmm_gamm1,da_glmm_gamm2)
#consistent improvement with season and status heterocedasticity.


#DA best model: gamma, allowing heterocedasticity across seasons and status (expected, altered will have more dispersion due to higher impact of ebullition). 
da_best_ch4<- da_glmm_gamm2
Anova(da_best_ch4,type = 3) 


#Cannot calculate R2 (conditional & marginal) for models with heterocedasticity,
#Get PseudoR2
pseudoR2_da_ch4<-pseudoR2_from_null(da_best_ch4, da_ch4)
pseudoR2_da_ch4

#R2 of homocedastic model
da_best_homoc_ch4<- da_glmm_gamm0
r2(da_best_homoc_ch4)
MuMIn::r.squaredGLMM(da_best_homoc_ch4) 

#Singularity?
performance::check_singularity(da_best_ch4)
#NO. 

#Extremely high improvement of R2 with subiste as random effect--> Subsite explains a lot of CH4 variability (expected field vs inundated farm)
#Overall no effect of status detected (only through interaction with season).
#Potentially for this casepilot it makes sense to separate altered1 and altered2 (very different conditions)

#remove non-best models:
rm(da_glmm_gaus0, da_glmm_stu0, da_glmm_gamm0,da_glmm_gamm1,da_glmm_gamm2)

## 2.4. DU model ----
#DU: lmer showed near-normal residuals, heterocedasticity for season and status
ch4_lmer_results %>% filter(casepilot=="DU")
#approach: compare fit between gausian, t_family models (logsign data) and gamma (shifted data), without heterocedasticity, with for season only and with for season and status. 


du_ch4<- ch4_formodels %>% filter(casepilot=="DU") %>% 
  filter(dailyflux<0.1)# remove outliers



du_ch4 %>% 
  mutate(dailyflux_shiftlog=log(dailyflux_shift)) %>% 
  pivot_longer(cols = c(dailyflux,dailyflux_logsign, dailyflux_shiftlog), names_to = "transformation", values_to = "value") %>% 
  ggplot(aes(x=value))+
  geom_histogram(aes(y = ..density..), bins = 50, fill = "skyblue", color = "black") +
  geom_density(color = "red", linewidth = 1) +
  facet_grid(casepilot~transformation, scales = "free") +
  theme_minimal()

du_ch4 %>% ggplot(aes(x=dailyflux, col=status))+
  geom_histogram()+
  scale_x_continuous(limits = c(-0.5,0.5))+ #ZOOM TO CHECK FOR OUTLIERS
  facet_wrap(~season)

#Gaussian with og data
du_glmm_gaus0<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite), 
                        data = du_ch4,
                        family = gaussian(),
                        dispformula = ~1)

#T-distribution with og data
du_glmm_stu0<- glmmTMB(formula = dailyflux~status*season + (1|subsite), 
                       data = du_ch4,
                       family = t_family,
                       dispformula = ~1) 

du_glmm_gamm0<- glmmTMB(formula = dailyflux_shift~status*season + (1|subsite), 
                        data = du_ch4,
                        family = Gamma(link = "log"),
                        dispformula = ~1)


resgaus0<- simulateResiduals(du_glmm_gaus0)
restu0<- simulateResiduals(du_glmm_stu0)
resgamm0<- simulateResiduals(du_glmm_gamm0)


#Check gaussian model
plot(resgaus0) #check normality, 
plotResiduals(resgaus0, du_ch4$season) #check heterocedasticity
plotResiduals(resgaus0, du_ch4$status) #check heterocedasticity
testOutliers(resgaus0) #Checks if there are more extreme residuals than expected
#Extreme non-normality. All asummptions fail, non-normality, funnel shape, deviation from Q-Q, heterocedasticity for both season and status

#Check t_family model
plot(restu0) #check normality, 
plotQQunif(restu0) #Check QQplot
plotResiduals(restu0)#check residuals 
plotResiduals(restu0, du_ch4$season) #check heterocedasticity
plotResiduals(restu0, du_ch4$status) #check heterocedasticity
testOutliers(restu0) #Checks if there are more extreme residuals than expected
testDispersion(restu0) #checks wether residual variance is higher or lower than expected, n.s.--> no issue. 
#Normal QQ but very weird pattern in residuals (residual calculation failed?)

#Check gamma model
plot(resgamm0) #check Q-Q, residuals vs fitted, 
plotResiduals(resgamm0, du_ch4$season) #check heterocedasticity
plotResiduals(resgamm0, du_ch4$status) #check heterocedasticity
testDispersion(resgamm0) #checks wether residual variance is higher or lower than expected, n.s.--> no issue. 
testOutliers(resgamm0) #Checks if there are more extreme residuals than expected
#BAD-looking residual plots 

#T-family looks best, lets try with heteorcedasticity at season or season+status

du_glmm_stu1<- glmmTMB(formula = dailyflux~status*season + (1|subsite), 
                        data = du_ch4,
                        family = t_family,
                        dispformula = ~season)

du_glmm_stu2<- glmmTMB(formula = dailyflux~status*season + (1|subsite), 
                        data = du_ch4,
                        family = t_family,
                        dispformula = ~season + status)

restu1<- simulateResiduals(du_glmm_stu1)
restu2<- simulateResiduals(du_glmm_stu2)

plot(restu1) #check Q-Q, residuals vs fitted, 
plotResiduals(restu1, du_ch4$season) #check heterocedasticity
plotResiduals(restu1, du_ch4$status) #check heterocedasticity
testDispersion(restu1) #checks wether residual variance is higher or lower than expected, n.s.--> no issue. 
testOutliers(restu1) #Checks if there are more extreme residuals than expected

plot(restu2) #check Q-Q, residuals vs fitted, 
plotQQunif(restu2) #Check QQplot
plotResiduals(restu2, du_ch4$season) #check heterocedasticity
plotResiduals(restu2, du_ch4$status) #check heterocedasticity
testDispersion(restu2) #checks wether residual variance is higher or lower than expected, n.s.--> no issue. 
testOutliers(restu2) #Checks if there are more extreme residuals than expected



du_glmm_gaus1<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite), 
                       data = du_ch4,
                       family = gaussian(),
                       dispformula = ~season)

du_glmm_gaus2<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite), 
                       data = du_ch4,
                       family = gaussian(),
                       dispformula = ~season + status)

resgaus1<- simulateResiduals(du_glmm_gaus1)
resgaus2<- simulateResiduals(du_glmm_gaus2)

plot(resgaus1) #check Q-Q, residuals vs fitted, 
plotResiduals(resgaus1, du_ch4$season) #check heterocedasticity
plotResiduals(resgaus1, du_ch4$status) #check heterocedasticity
testDispersion(resgaus1) #checks wether residual variance is higher or lower than expected, n.s.--> no issue. 
testOutliers(resgaus1) #Checks if there are more extreme residuals than expected

plot(resgaus2) #check Q-Q, residuals vs fitted, 
plotQQunif(resgaus2) #Check QQplot
plotResiduals(resgaus2, du_ch4$season) #check heterocedasticity
plotResiduals(resgaus2, du_ch4$status) #check heterocedasticity
testDispersion(resgaus2) #checks wether residual variance is higher or lower than expected, n.s.--> no issue. 
testOutliers(resgaus2) #Checks if there are more extreme residuals than expected


anova(du_glmm_gaus0, du_glmm_gaus1, du_glmm_gaus2)

#Student model with better AIC and ok residuals
anova(du_glmm_gaus0, du_glmm_stu0, du_glmm_gamm0)

#Gamma distribution is is consistently better (if we discard student model)
anova(du_glmm_stu0,du_glmm_stu1,du_glmm_stu2)
#consistent improvement with season and status heterocedasticity.


#DA best model: t_distributiion with OG data, allowing heterocedasticity across seasons and status. 
du_best_ch4<- du_glmm_gaus2
Anova(du_best_ch4,type = 3) 

#Cannot calculate R2 (conditional & marginal) for models with heterocedasticity,
#Get PseudoR2
pseudoR2_du_ch4<-pseudoR2_from_null(du_best_ch4, du_ch4)
pseudoR2_du_ch4


#R2 of homocedastic model
du_best_homoc_ch4<- du_glmm_gaus0
r2(du_best_homoc_ch4)
MuMIn::r.squaredGLMM(du_best_homoc_ch4) 

#Singularity?
performance::check_singularity(du_best_ch4)
#YES. 

#VERY LOW R2, Ch4 is not varying by season, status

#remove non-best models:
rm(du_glmm_gaus0,du_glmm_gaus1,du_glmm_gaus2, du_glmm_stu0,du_glmm_stu1,du_glmm_stu2, du_glmm_gamm0)




## 2.5. RI model ----
#RI: lmer showed normal residuals, heterocedasticity for status (but not for season)
ch4_lmer_results %>% filter(casepilot=="RI")

ch4_formodels %>% 
  filter(casepilot=="RI") %>% 
  mutate(dailyflux_shiftlog=log(dailyflux_shift)) %>% 
  pivot_longer(cols = c(dailyflux,dailyflux_logsign, dailyflux_shiftlog), names_to = "transformation", values_to = "value") %>% 
  ggplot(aes(x=value))+
  geom_histogram(aes(y = ..density..), bins = 30, fill = "skyblue", color = "black") +
  geom_density(color = "red", linewidth = 1) +
  facet_grid(casepilot~transformation, scales = "free") +
  theme_minimal()
#Logsign does not change distribution of data for RI, use raw values for gaus and student models
ch4_formodels %>% 
  filter(casepilot=="RI") %>% 
  mutate(dailyflux_shiftlog=log(dailyflux_shift)) %>% 
  pivot_longer(cols = c(dailyflux), names_to = "transformation", values_to = "value") %>% 
  ggplot(aes(x=value))+
  geom_histogram(aes(y = ..density..), bins = 50, fill = "skyblue", color = "black") +
  geom_density(color = "red", linewidth = 1) +
  geom_vline(xintercept = -0.02)+
  geom_vline(xintercept = 0.04)+
  
  facet_grid(casepilot~status) +
  theme_minimal()

#visualize outliers
ch4_formodels %>% 
  filter(casepilot=="RI") %>%  
  mutate(subsite2=substr(subsite,4,5)) %>% 
  ggplot(aes(x=subsite2, y=dailyflux_logsign, fill=status, group=paste0(subsite)))+
  geom_violin()+
  geom_hline(yintercept = -0.025)+
  geom_hline(yintercept = 0.025)+
  facet_wrap(~casepilot,scales="free")+
  scale_y_continuous(limits = c(-0.025,0.0251))

#approach: compare fit between gausian, t_family models (raw data) and gamma (shifted data), without heterocedasticity, with for season only and with for season and status. 



ri_ch4<- ch4_formodels %>% filter(casepilot=="RI") %>% 
  filter(between(dailyflux, -0.025,0.025)) #CH4: limits = c(-0.025,0.025), (removes 8 outliers)

plot_histograms(ri_ch4$dailyflux)
ri_glmm_gaus0<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite), 
                        data = ri_ch4,
                        family = gaussian(),
                        dispformula = ~1)

ri_glmm_stu0<- glmmTMB(formula = dailyflux~status*season + (1|subsite), 
                       data = ri_ch4,
                       family = t_family,
                       dispformula = ~1) 

ri_glmm_gamm0<- glmmTMB(formula = dailyflux_shift~status*season + (1|subsite), 
                        data = ri_ch4,
                        family = Gamma(link = "log"),
                        dispformula = ~1)


resgaus0<- simulateResiduals(ri_glmm_gaus0)
restu0<- simulateResiduals(ri_glmm_stu0)
resgamm0<- simulateResiduals(ri_glmm_gamm0)


#Check gaussian model
plot(resgaus0) #check normality, 
plotResiduals(resgaus0, ri_ch4$season) #check heterocedasticity
plotResiduals(resgaus0, ri_ch4$status) #check heterocedasticity
testOutliers(resgaus0) #Checks if there are more extreme residuals than expected
#Extreme non-normality. All asummptions fail, non-normality, funnel shape, deviation from Q-Q, heterocedasticity for status

#Check t_family model
plot(restu0) #check normality, 
plotQQunif(restu0) #Check QQplot
plotResiduals(restu0)#check residuals 
plotResiduals(restu0, ri_ch4$season) #check heterocedasticity
plotResiduals(restu0, ri_ch4$status) #check heterocedasticity
testOutliers(restu0) #Checks if there are more extreme residuals than expected
testDispersion(restu0) #checks wether residual variance is higher or lower than expected, n.s.--> no issue. 
#Extremely good QQplot and residual plots. Homocedasticity for both season and status

#Check gamma model
plot(resgamm0) #check Q-Q, residuals vs fitted, 
plotResiduals(resgamm0, ri_ch4$season) #check heterocedasticity
plotResiduals(resgamm0, ri_ch4$status) #check heterocedasticity
testDispersion(resgamm0) #checks wether residual variance is higher or lower than expected, n.s.--> no issue. 
testOutliers(resgamm0) #Checks if there are more extreme residuals than expected
#BAD-looking residual plots 

#Student looks great! lets try with heteorcedasticity at season or season*status

ri_glmm_stu1<- glmmTMB(formula = dailyflux~status*season + (1|subsite), 
                        data = ri_ch4,
                        family = t_family,
                        dispformula = ~season)


ri_glmm_stu2<- glmmTMB(formula = dailyflux~status*season + (1|subsite), 
                       data = ri_ch4,
                       family = t_family,
                       dispformula = ~season + status)

ri_glmm_gaus1<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite), 
                       data = ri_ch4,
                       family = gaussian(),
                       dispformula = ~season)


ri_glmm_gaus2<- glmmTMB(formula = dailyflux~status*season + (1|subsite), 
                       data = ri_ch4,
                       family = gaussian(),
                       dispformula = ~season + status)


restu1<- simulateResiduals(ri_glmm_stu1, n=1000)
resg1<- simulateResiduals(ri_glmm_gaus1)
resg2<- simulateResiduals(ri_glmm_gaus2)


plot(restu1) #check Q-Q, residuals vs fitted, 
plotResiduals(restu1, ri_ch4$season) #check heterocedasticity
plotResiduals(restu1, ri_ch4$status) #check heterocedasticity
testDispersion(restu1) #checks wether residual variance is higher or lower than expected, n.s.--> no issue. 
testOutliers(restu1) #Checks if there are more extreme residuals than expected
#Heterocedastic model for season also looks good. 


plot(resg1) #check Q-Q, residuals vs fitted, 
plotResiduals(resg1, ri_ch4$season) #check heterocedasticity
plotResiduals(resg1, ri_ch4$status) #check heterocedasticity
testDispersion(resg1) #checks wether residual variance is higher or lower than expected, n.s.--> no issue. 
testOutliers(resg1) #Checks if there are more extreme residuals than expected
#Heterocedastic model for season also looks good. 



#COMPARE AIC between families
anova(ri_glmm_gaus0, ri_glmm_stu0, ri_glmm_gamm0)
#Student distribution is Much better
anova(ri_glmm_stu0,ri_glmm_stu1)
#Adding heterocedasticity does not improve fit. 

anova(ri_glmm_gaus0,ri_glmm_gaus1,ri_glmm_gaus2)

#RI best model: student family, without allowing for heterocedasticity. 
ri_best_ch4<- ri_glmm_gaus0
Anova(ri_best_ch4,type = 3) 


#Calculate pseudoR2 for consistency (with other models)
#Get PseudoR2
pseudoR2_ri_ch4<-pseudoR2_from_null(ri_best_ch4, ri_ch4)
pseudoR2_ri_ch4



#R2 of homocedastic model
ri_best_homoc_ch4<- ri_glmm_gaus0
r2(ri_best_homoc_ch4)
MuMIn::r.squaredGLMM(ri_best_homoc_ch4) 

#Singularity?
performance::check_singularity(ri_best_ch4)
#YES 


#WEIRD: student model shows extremely clean residuals and very low AIC (better than the rest). It has singularity (ie. subsite does not add info). IT SHOWS strongly significant effects, but R2 is extremely LOW. ?!?!
#this might be caused by outliers contributing a lot to residual variance, but fixed effects explaining relatively low of the total variance. Not to worry, check distribution of data and prediction capacity of model to confirm low R2 is due to outliers. 

#Logsign transformation does not change the distribution of RI data, We could fit the model directly with raw-fluxes. SHould we remove the outliers?
#NO; better to leave them. More conceptually sound and shows the strength of the pattern: "even in the presence of a few weird outliers, we can detect a strong effect of restoration/alteration"



# Extract fitted values (on response scale)
fitted_vals <- predict(ri_best_ch4, type = "response")

# Extract observed values
observed_vals <- ri_best_ch4$frame$dailyflux_logsign  # response column used in model

#Plot obs vs pred (no outliers )
data.frame(obs=observed_vals, pred=fitted_vals) %>% 
  # filter(between(obs, -0.025,0.035)) %>% #exclude outliers
  ggplot(aes(x=pred, y=obs))+
  geom_point()+
  geom_smooth(method="lm")+
  geom_abline(intercept = 0, slope = 1, col="red")

#MODEL IS OK, low R2 is due to outliers. We should not remove them (more robust). 

#remove non-best models:
rm(ri_glmm_gamm0,ri_glmm_gaus0,,ri_glmm_gaus1,ri_glmm_gaus2,ri_glmm_stu0,ri_glmm_stu1,ri_glmm_stu2)

## 2.6. VA model ----
#VA: lmer showed non-normal residuals, heterocedasticity for season (and status?),
ch4_lmer_results %>% filter(casepilot=="VA")
#approach: compare fit between gausian, t_family models (logsign data) and gamma (shifted data), without heterocedasticity, with for season only and with for season and status. 

ch4_formodels %>% 
  filter(casepilot=="VA") %>% 
  mutate(dailyflux_shiftlog=log(dailyflux_shift)) %>% 
  pivot_longer(cols = c(dailyflux,dailyflux_logsign, dailyflux_shiftlog), names_to = "transformation", values_to = "value") %>% 
  ggplot(aes(x=value))+
  geom_histogram( bins = 30, fill = "skyblue", color = "black") +
  # geom_density(color = "red", linewidth = 1) +
  facet_grid(casepilot~transformation, scales = "free") +
  theme_minimal()


va_ch4<- ch4_formodels %>% filter(casepilot=="VA")

va_ch4 %>% 
  mutate(dailyflux_shiftlog=log(dailyflux_shift)) %>%
  ggplot(aes(x=dailyflux_shiftlog, col=status))+
  geom_boxplot()+
  # scale_x_continuous(limits = c(-0.5,0.5))+ #ZOOM TO CHECK FOR OUTLIERS
  facet_wrap(~season)


va_glmm_gaus0<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite), 
                        data = va_ch4,
                        family = gaussian(),
                        dispformula = ~1)

va_glmm_stu0<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite), 
                       data = va_ch4,
                       family = t_family,
                       dispformula = ~1) 

va_glmm_gamm0<- glmmTMB(formula = dailyflux_shift~status*season + (1|subsite), 
                        data = va_ch4,
                        family = Gamma(link = "log"),
                        dispformula = ~1)


resgaus0<- simulateResiduals(va_glmm_gaus0)
restu0<- simulateResiduals(va_glmm_stu0)
resgamm0<- simulateResiduals(va_glmm_gamm0)

#Check gaussian model
plot(resgaus0) #check normality, 
plotResiduals(resgaus0, va_ch4$season) #check heterocedasticity
plotResiduals(resgaus0, va_ch4$status) #check heterocedasticity
#All asummptions fail, non-normality, funnel shape, deviation from Q-Q, heterocedasticity for both season and status

#Check t_family model
plot(restu0) #check normality, 
plotResiduals(restu0, va_ch4$season) #check heterocedasticity
plotResiduals(restu0, va_ch4$status) #check heterocedasticity
#All asumptions fail, not convergence

#Check gamma model
summary(va_glmm_gamm0)
plot(resgamm0) #check Q-Q, residuals vs fitted, 
plotResiduals(resgamm0, va_ch4$season) #check heterocedasticity
plotResiduals(resgamm0, va_ch4$status) #check heterocedasticity
testDispersion(resgamm0) #checks wether residual variance is higher or lower than expected, n.s.--> no issue. 
testOutliers(resgamm0) #Checks if there are more extreme residuals than expected

#Gamma looks best, lets try with heteorcedasticity at season or season*status

va_glmm_gamm1<- glmmTMB(formula = dailyflux_shift~status*season + (1|subsite), 
                        data = va_ch4,
                        family = Gamma(link = "log"),
                        dispformula = ~season)

va_glmm_gamm2<- glmmTMB(formula = dailyflux_shift~status*season + (1|subsite), 
                        data = va_ch4,
                        family = Gamma(link = "log"),
                        dispformula = ~season + status)

resgamm1<- simulateResiduals(va_glmm_gamm1)
resgamm2<- simulateResiduals(va_glmm_gamm2)

plot(resgamm1) #check Q-Q, residuals vs fitted, 
plotResiduals(resgamm1, va_ch4$season) #check heterocedasticity
plotResiduals(resgamm1, va_ch4$status) #check heterocedasticity
testDispersion(resgamm1) #checks wether residual variance is higher or lower than expected, n.s.--> no issue. 
testOutliers(resgamm1) #Checks if there are more extreme residuals than expected

plot(resgamm2) #check normality, 
plotQQunif(resgamm2) #Check QQplot
plotResiduals(resgamm2, va_ch4$season) #check heterocedasticity
plotResiduals(resgamm2, va_ch4$status) #check heterocedasticity
testDispersion(resgamm2) #checks wether residual variance is higher or lower than expected, n.s.--> no issue. 
testOutliers(resgamm2) #Checks if there are more extreme residuals than expected

#GAUssian and student distribution models do not meet assumptions

#Gamma distribution is is consistently better 
anova(va_glmm_gamm0,va_glmm_gamm1,va_glmm_gamm2)
#Consistent improvement by including heterocedasticity 


#VA best model: gamma, allowing heterocedasticity across seasons and status (expected, altered will have more dispersion due to higher impact of ebullition). 
va_best_ch4<- va_glmm_gamm2
Anova(va_best_ch4,type = 3) 


#Cannot calculate R2 (conditional & marginal) for models with heterocedasticity,
#Get PseudoR2
pseudoR2_va_ch4<-pseudoR2_from_null(va_best_ch4, va_ch4)
pseudoR2_va_ch4

#R2 of homocedastic model
va_best_homoc_ch4<- va_glmm_gamm0
r2(va_best_homoc_ch4)
MuMIn::r.squaredGLMM(va_best_homoc_ch4) 

#Singularity?
performance::check_singularity(va_best_ch4)
#NO. 

#Good R2, relatively good residuals, Detected effect for both status and season and interactions 
#AIC of gaussian and student models is better but they fail their assumptions. 

#Remove non-best models: 
rm(va_glmm_gamm0,va_glmm_gamm1, va_glmm_gamm2, va_glmm_gaus0,va_glmm_stu0)

#3. (skip) Optimize with weights --------

#Here we will create glmms taking into account strata sample_weights. 
#Then compare models with and without sample_weights to see the effect of weights and get the overall best model for each casepilot. 

# DOUBT: should I directly implement the identified best_model option for each casepilot adding sample_weight OR optimize from the beginning (allowing different heterocedasticity structures)??:

#Applying weights to the already optimized structure can tell us directly if sample weights improve fit. 
#Optimizing the models from the beginning with sample weights will allow us to check whether the weights influence the overall structure selection BUT we will not be able to directly compare with anova() if models differ in structure (i.e. if they have different fixed/random effects, different dispformula, or different distribution families) 

#We know that we should account for sample weights in the end if/when we get accurate strata composition and if they show substantial and consistent sampling bias. 

## 3.1. CA model (weighted)

#Do sample weights improve the optimized model? 
#Recal best model structure: 
summary(ca_best)

#build the same model with weights
ca_best_weighted<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite),
                           data = ca_co2,
                           family = gaussian(),
                           dispformula = ~season, #heterocedasticity across season levels (not status level)
                           weights=sample_weight)

#Compare fit between weigthed and unweighted model
anova(ca_best, ca_best_weighted)
#inclusion of sample_weight does not significantly improve model. Unclear effect on Df but weights cause no improvement on AIC...

Anova(ca_best,type = 3)
Anova(ca_best_weighted,type = 3)
#Weights also do not change significance of main effects. 

#Do weights lead to different structure? #Run the optimization process to select best model (with/without heterocedasticity) Not essential to check, omit. 


## 3.2. CU model (weighted)

#Do sample weights improve the optimized model? 
#Recal best model structure: 
summary(cu_best)

#build the same model with weights
cu_best_weighted<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite),
                           data = cu_co2,
                           family = t_family,
                           dispformula = ~season, #heterocedasticity across season levels (not status level)
                           weights=sample_weight)

#Compare fit between weigthed and unweighted model
anova(cu_best, cu_best_weighted)
#inclusion of sample_weight does significantly improve model. Unclear effect on Df but weights cause  improvement on AIC...

Anova(cu_best,type = 3)
Anova(cu_best_weighted,type = 3)
#Weights make main effects sigltly more significant.
summary(cu_best_weighted)

#Do weights lead to different structure? #Run the optimization process to select best model (with/without heterocedasticity) Not essential to check, omit. 



## 3.3. DA model (weighted)

#Do sample weights improve the optimized model? 
#Recal best model structure: 
summary(da_best)

#build the same model with weights
da_best_weighted<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite),
                           data = da_co2,
                           family = t_family,
                           dispformula = ~season, #heterocedasticity across season levels (not status level)
                           weights=sample_weight)

#Compare fit between weigthed and unweighted model
anova(da_best, da_best_weighted)
#inclusion of sample_weight does NOT significantly improve model. Unclear effect on Df....

Anova(da_best,type = 3)
Anova(da_best_weighted, type = 3)
#Weights make main effects sigltly less significant.

#Do weights lead to different structure? #Run the optimization process to select best model (with/without heterocedasticity) Not essential to check, omit. 



## 3.4. DU model (weighted)

#Do sample weights improve the optimized model? 
#Recal best model structure: 
summary(du_best)

#build the same model with weights
du_best_weighted<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite),
                           data = du_co2,
                           family = gaussian(),
                           dispformula = ~season, #heterocedasticity across season levels (not status level)
                           weights=sample_weight)

#Compare fit between weigthed and unweighted model
anova(du_best, du_best_weighted)
#inclusion of sample_weight does NOT significantly improve model. Unclear effect on Df....

#No detectable effect of status
Anova(du_best,type = 3)
Anova(du_best_weighted, type = 3)

#Do weights lead to different structure? #Run the optimization process to select best model (with/without heterocedasticity) Not essential to check, omit. 

## 3.5. RI model (weighted)

#Do sample weights improve the optimized model? 
#Recal best model structure: 
summary(ri_best)

#build the same model with weights
ri_best_weighted<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite),
                           data = ri_co2,
                           family = gaussian(),
                           # dispformula = ~season, #HOMOCEDASTICITY assumed
                           weights=sample_weight)

#Compare fit between weigthed and unweighted model
anova(ri_best, ri_best_weighted)
#inclusion of sample_weight does NOT significantly improve model. Unclear effect on Df....


#Same significance of main effects
Anova(ri_best,type = 3)
Anova(ri_best_weighted,type = 3)

#Do weights lead to different structure? #Run the optimization process to select best model (with/without heterocedasticity) Not essential to check, omit. 


## 3.6. VA model (weighted)

#Do sample weights improve the optimized model? 
#Recal best model structure: 
summary(va_best)

#build the same model with weights
va_best_weighted<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite),
                           data = va_co2,
                           family = gaussian(),
                           dispformula = ~season, #HOMOCEDASTICITY assumed
                           weights=sample_weight)

#Compare fit between weigthed and unweighted model
anova(va_best, va_best_weighted)
#inclusion of sample_weight does NOT significantly improve model. Unclear effect on Df....


#Same significance of main effects
Anova(va_best,type = 3)
Anova(va_best_weighted, type = 3)

#Do weights lead to different structure? #Run the optimization process to select best model (with/without heterocedasticity) Not essential to check, omit. 



#4. Effects of best-models-------


#BASED ON THE ABOVE OPTIMIZATIONS, provide list of models to evaluate/interpret

#Provide list of homocedastic model variants (to estimate R2)
homoced_model_list_ch4<- list(
  "CA" = ca_best_homoc_ch4,
  "CU" = cu_best_homoc_ch4, 
  "DA" = da_best_homoc_ch4,
  "DU" = du_best_homoc_ch4,
  "RI" = ri_best_homoc_ch4, #Only casepilot for which homocedastic model is best
  "VA" = va_best_homoc_ch4)

#List of best models for each casepilot (to get residual tests and Significances)
best_model_list_ch4<- list(
  "CA" = ca_best_ch4,
  "CU" = cu_best_ch4,
  "DA" = da_best_ch4,
  "DU" = du_best_ch4,
  "RI" = ri_best_ch4,
  "VA" = va_best_ch4)

#List of dataframes used to build each model (needed to calculate residuals)
model_data_list_ch4<- list(
  "CA" = ca_ch4,
  "CU" = cu_ch4, 
  "DA" = da_ch4,
  "DU" = du_ch4,
  "RI" = ri_ch4,
  "VA" = va_ch4)

## 4.1. Overall fit -----
#1st how good is the model--> overall R2 for homocedastic model variant 

ch4_model_fit_summary<- summarize_model_info(homoced_model_list_ch4) %>% 
  mutate(ghgspecies="gwp_ch4")
ch4_model_fit_summary

## 4.2. Residual checks---
#Summary of Model assumptions: for best_model_list using function created in Co2 section:

ch4_resid_diag_summary <- summarize_dharma_diagnostics(best_model_list_ch4, model_data_list_ch4)
print(ch4_resid_diag_summary)

#CHECK-----
#IS  uniformity needed with gamma family? 

# _________
#Residual tests are ok for most models with 2 exceptions: 
#dispersion for CU model (this model is always singular, and sample weights are not accounted for in these tests, likely ok) 
#heterocedasticity for status and season for VA model (this model is very bad: extremely low r2, and even acounting for season heterocedasticity, residuals are still not good, we cannot obtain a interpretable model for VA co2)


##4.2. Significance of effects-----
#Effect of status will first be assessed via Anova (best_model,type=3): is there an effect averaging across all seasons? yes/no, is the effect dependent on season? 

#Function to summarise ANOVA results created in co2 section

#Obtain significance for main effects of best_models
ch4_results_anova<-summarize_anova_results(best_model_list_ch4) %>% 
  mutate(ghgspecies="gwp_ch4")
print(ch4_results_anova)

#Intercept refers to whether the grand mean (on the transformed scale) is significantly different than cero (Is the overall net exchange of CH4 balanced between uptake and emission)
#Status refers to whether the status (averaged across seasons) has a significant effect on flux. 
#Season refers to whether the season (averaged across status) has a significant effect on flux.
#Interaction refers to whether the effect of status depends on season — i.e., the difference between status levels changes across seasons.



#Look for significance at either status or status:season---> there is an effect (that may or may not be evident when averaging across seasons)
#IF status:season significant---> effect depends on season. 

#Does status have an effect (significance of status or interaction)
ch4_results_anova %>% filter(status_pval<0.05|interaction_pval<0.05)

#ALL except DA show average effect of status 
#ALL casepilots show significant effect of interaction


#Save in dropbox tables with the models summary outputs: 

write.csv(x = ch4_model_fit_summary, file = paste0(plots_path,"/CH4models_homocedastic_R2.csv"),row.names = F)
write.csv(x = ch4_resid_diag_summary, file = paste0(plots_path,"/CH4models_tests_residuals.csv"),row.names = F)
write.csv(x = ch4_results_anova, file = paste0(plots_path,"/CH4models_effects_signififcance.csv"),row.names = F)




##4.3. Effect direction and magnitude----

#Use model scale for statistics, back-transform for interpretation.  

#NOTE on bootstrapping (omit): there is the possibility to bootstrap the models (i.e. recalculate the model many times with re-sampled data, and calculate emmeans for each model, to be able to produce more robust Confidence Intervals for them, but this is not usually done and computationally intensive. We would have to take into account the model structure when re-sampling the data to ensure balanced resamplings (all factors represented by data). WE will not do this (for the moment).

#Some guides to plot and use emmeans: https://rcompanion.org/handbook/G_06.html


#Notes on non-gaussian models: upper and lower confidence limits are given as asymp.LCL/asymp.UCL (in our context we can treat them equally as the ones from gaussian models). Additionally, non-gaussian models will have Inf df (degrees of freedom are calculated differently for T-family models), not an issue. 


# Loop to extract all comparisons for every casepilot best model

#The loop has to take into account which unit was modeled in each case, perform the back-transformation outside the loop. 

# Initialize list for storing comparisons
comparisons_list <- list()

for (casepilot_name in names(best_model_list_ch4)) {
  ghgspecies<- "gwp_ch4"
  cp_model <- best_model_list_ch4[[casepilot_name]]
  model_call<- deparse(formula(cp_model))
  # Status comparisons:
  status_emmeans <- emmeans(cp_model, ~status, type="response")
  status_comparisons <- cld(status_emmeans, 
                            alpha = 0.05,
                            Letters = letters,
                            adjust = "sidak") %>%
    mutate(comparison = "status",
           ghgspecies = ghgspecies,
           casepilot = casepilot_name,
           season = NA)%>%    
    #Rename confidence limits for T_family models (if this is the case)
    rename(lower.CL = any_of("asymp.LCL"),
           upper.CL = any_of("asymp.UCL"))
  
  # Season comparisons: 
  season_emmeans <- emmeans(cp_model, ~season, type="response")
  season_comparisons <- cld(season_emmeans, 
                            alpha = 0.05,
                            Letters = letters,
                            adjust = "sidak") %>%
    mutate(comparison = "season",
           ghgspecies = ghgspecies,
           casepilot = casepilot_name,
           status = NA)%>%
    #Rename confidence limits for T_family models (if this is the case)
    rename(lower.CL = any_of("asymp.LCL"),
           upper.CL = any_of("asymp.UCL"))
  
  
  # Interaction comparisons: 
  interaction_emmeans <- emmeans(cp_model, ~status:season, type="response")
  interaction_comparisons <- cld(interaction_emmeans, 
                                 alpha = 0.05,
                                 Letters = letters,
                                 adjust = "sidak") %>%
    mutate(comparison = "interaction",
           ghgspecies = ghgspecies,
           casepilot = casepilot_name)%>%
    #Rename confidence limits for T_family models (if this is the case)
    rename(lower.CL = any_of("asymp.LCL"),
           upper.CL = any_of("asymp.UCL"))
  
  
  # Join and format CP comparisons
  all_comparisons <- status_comparisons %>%
    full_join(season_comparisons) %>%
    full_join(interaction_comparisons) %>%
    rename(group_letter = .group) %>%
    mutate(
      comparison = factor(comparison, levels = c("status", "season", "interaction")),
      season = factor(season, levels = c("S1", "S2", "S3", "S4")),
      status = factor(status, levels = c("Altered", "Preserved", "Restored")),
      group_letter = gsub(" ", "", group_letter)
    ) %>%
    mutate(model_call=model_call) %>% 
    dplyr::select(any_of(c(
      "ghgspecies", "casepilot", "model_call","comparison", "status", "season", "df",
      "emmean", "response","SE", "lower.CL", "upper.CL", "group_letter"
    ))) %>%
    arrange(ghgspecies, casepilot, comparison, status, season)
  
  # Store casepilot results
  comparisons_list[[casepilot_name]] <- all_comparisons
}

#Join all EMMs in a single data frame and back transform where apropriate to _gwp
ch4_casepilot_comparisons <- bind_rows(comparisons_list) %>% 
  #add shift for each casepilot
  left_join(casepilot_ch4_shifts) %>% 
  #Put response from gamma models into variable emmean (do not do for RI, which has correct emmean)
  mutate(emmean=if_else(is.na(emmean),response,emmean)) %>% 
  #Back transform shifted emmean, lower.CL and upper.CL to original scale (for CA,CU,DA,DU, VA), RI and DU are fitted in original scale (do nothing). As SE only measures spread, it is the same in the shifted and the original scales. 
  mutate(emmean_gwp=case_when(!casepilot%in%c("RI","DU")~emmean-shift,
                              TRUE~emmean),
         SE_gwp=SE,
         lower.CL_gwp=case_when(!casepilot%in%c("RI","DU")~lower.CL-shift,
                            TRUE~lower.CL),
         upper.CL_gwp=case_when(!casepilot%in%c("RI","DU")~upper.CL-shift,
                            TRUE~upper.CL))


#Save all comparisons:  
write.csv(x = ch4_casepilot_comparisons, file = paste0(plots_path,"/CH4models_emmeans_posthoc.csv"),row.names = F)



#Plot casepilot-per-casepilot, emmeans +- SE (in gwp scale)

##Status effect plots-----

#Camargue: 
plot_ch4_status_CA<-ch4_casepilot_comparisons %>%
  filter(casepilot == "CA", comparison == "status") %>%
    ggplot(., aes(x = status, y = emmean_gwp, label = group_letter)) +
      geom_errorbar(aes(ymin = emmean_gwp - SE_gwp, ymax = emmean_gwp + SE_gwp),
                    linewidth = 0.2, width = 0.2) +
      geom_point(shape = 15, size = 4, aes(col=status)) +
      geom_text(nudge_x = 0.2, color = "black") +
      geom_hline(yintercept=0, linetype="dashed")+
      theme_bw() +
      theme(axis.title = element_text(face = "bold"),
            axis.text = element_text(face = "bold"),
            plot.caption = element_text(hjust = 0)) +
      labs(title = expression(Status~effect~on~net~CH[4]~flux),
           subtitle = paste0("Camargue"), 
           y=expression(EMM~(g~CO[2~eq]~d^-1~m^-2)),
           x="Conservation status",
           col=paste0("Status"))
  
plot_ch4_status_CA
ggsave(plot=plot_ch4_status_CA, filename="CH4_status_plot_CA.png",
       device = "png",path = plots_path,dpi = 400,
       width = 115, height = 100,units = "mm")


#Curonian Lagoon: 
plot_ch4_status_CU<- ch4_casepilot_comparisons %>%
  filter(casepilot == "CU", comparison == "status") %>%
    ggplot(., aes(x = status, y = emmean_gwp, label = group_letter)) +
      geom_errorbar(aes(ymin = emmean_gwp - SE_gwp, ymax = emmean_gwp + SE_gwp),
                    linewidth = 0.2, width = 0.2) +
      geom_point(shape = 15, size = 4, aes(col=status)) +
      geom_text(nudge_x = 0.2, color = "black") +
      geom_hline(yintercept=0, linetype="dashed")+
      theme_bw() +
      theme(axis.title = element_text(face = "bold"),
            axis.text = element_text(face = "bold"),
            plot.caption = element_text(hjust = 0)) +
      labs(title = expression(Status~effect~on~net~CH[4]~flux),
           subtitle = paste0("Curonian lagoon"), 
           y=expression(EMM~(g~CO[2~eq]~d^-1~m^-2)),
           x="Conservation status",
           col=paste0("Status"))
  
plot_ch4_status_CU
ggsave(plot=plot_ch4_status_CU, filename="CH4_status_plot_CU.png",
       device = "png",path = plots_path,dpi = 400,
       width = 115, height = 100,units = "mm")

#Danube delta: 
plot_ch4_status_DA<- ch4_casepilot_comparisons %>%
  filter(casepilot == "DA", comparison == "status") %>%
  ggplot(., aes(x = status, y = emmean_gwp, label = group_letter)) +
      geom_errorbar(aes(ymin = emmean_gwp - SE_gwp, ymax = emmean_gwp + SE_gwp),
                    linewidth = 0.2, width = 0.2) +
      geom_point(shape = 15, size = 4, aes(col=status)) +
      geom_text(nudge_x = 0.2, color = "black") +
      geom_hline(yintercept=0, linetype="dashed")+
      theme_bw() +
      theme(axis.title = element_text(face = "bold"),
            axis.text = element_text(face = "bold"),
            plot.caption = element_text(hjust = 0)) +
      labs(title = expression(Status~effect~on~net~CH[4]~flux),
           subtitle = paste0("Danube delta"), 
           y=expression(EMM~(g~CO[2~eq]~d^-1~m^-2)),
           x="Conservation status",
           col=paste0("Status"))

plot_ch4_status_DA

ggsave(plot=plot_ch4_status_DA, filename="CH4_status_plot_DA.png",
       device = "png",path = plots_path,dpi = 400,
       width = 115, height = 100,units = "mm")

#Dutch delta: 
plot_ch4_status_DU<- ch4_casepilot_comparisons %>%
  filter(casepilot == "DU", comparison == "status") %>%
    ggplot(., aes(x = status, y = emmean_gwp, label = group_letter)) +
      geom_errorbar(aes(ymin = emmean_gwp - SE_gwp, ymax = emmean_gwp + SE_gwp),
                    linewidth = 0.2, width = 0.2) +
      geom_point(shape = 15, size = 4, aes(col=status)) +
      geom_text(nudge_x = 0.2, color = "black") +
      geom_hline(yintercept=0, linetype="dashed")+
      theme_bw() +
      theme(axis.title = element_text(face = "bold"),
            axis.text = element_text(face = "bold"),
            plot.caption = element_text(hjust = 0)) +
      labs(title = expression(Status~effect~on~net~CH[4]~flux),
           subtitle = paste0("Dutch delta"),
           y=expression(EMM~(g~CO[2~eq]~d^-1~m^-2)),
           x="Conservation status",
           col=paste0("Status"))

plot_ch4_status_DU

ggsave(plot=plot_ch4_status_DU, filename="CH4_status_plot_DU.png",
       device = "png",path = plots_path,dpi = 400,
       width = 115, height = 100,units = "mm")


#Ria d'Aveiro: 
plot_ch4_status_RI<-ch4_casepilot_comparisons %>%
  filter(casepilot == "RI", comparison == "status") %>%
    ggplot(., aes(x = status, y = emmean_gwp, label = group_letter)) +
      geom_errorbar(aes(ymin = emmean_gwp - SE_gwp, ymax = emmean_gwp + SE_gwp),
                    linewidth = 0.2, width = 0.2) +
      geom_point(shape = 15, size = 4, aes(col=status)) +
      geom_text(nudge_x = 0.2, color = "black") +
      geom_hline(yintercept=0, linetype="dashed")+
      theme_bw() +
      theme(axis.title = element_text(face = "bold"),
            axis.text = element_text(face = "bold"),
            plot.caption = element_text(hjust = 0)) +
      labs(title = expression(Status~effect~on~net~CH[4]~flux),
           subtitle = paste0("Ria d'Aveiro"), 
           y=expression(EMM~(g~CO[2~eq]~d^-1~m^-2)),
           x="Conservation status",
           col=paste0("Status"))

plot_ch4_status_RI
ggsave(plot=plot_ch4_status_RI, filename="CH4_status_plot_RI.png",
       device = "png",path = plots_path,dpi = 400,
       width = 115, height = 100,units = "mm")

#Valencian wetland: 
plot_ch4_status_VA<-ch4_casepilot_comparisons %>%
  filter(casepilot == "VA", comparison == "status") %>%
  ggplot(., aes(x = status, y = emmean_gwp, label = group_letter)) +
      geom_errorbar(aes(ymin = emmean_gwp - SE_gwp, ymax = emmean_gwp + SE_gwp),
                    linewidth = 0.2, width = 0.2) +
      geom_point(shape = 15, size = 4, aes(col=status)) +
      geom_text(nudge_x = 0.2, color = "black") +
      geom_hline(yintercept=0, linetype="dashed")+
      theme_bw() +
      theme(axis.title = element_text(face = "bold"),
            axis.text = element_text(face = "bold"),
            plot.caption = element_text(hjust = 0)) +
      labs(title = expression(Status~effect~on~net~CH[4]~flux),
           subtitle = paste0("Valencian wetland"), 
           y=expression(EMM~(g~CO[2~eq]~d^-1~m^-2)),
           x="Conservation status",
           col=paste0("Status"))

plot_ch4_status_VA
ggsave(plot=plot_ch4_status_VA, filename="CH4_status_plot_VA.png",
       device = "png",path = plots_path,dpi = 400,
       width = 115, height = 100,units = "mm")


grid.arrange(plot_ch4_status_CA,plot_ch4_status_CU,plot_ch4_status_DA, plot_ch4_status_DU, plot_ch4_status_RI, plot_ch4_status_VA, ncol=3)
grid.arrange(plot_ch4_status_CA, plot_ch4_status_CU,p3, ncol = 3)



##Seasonal effect plots-----
#To do, not too much interest for overview ppt (seasonal effect is averaged across the three status)

#Example: CA
#Camargue: seasonal plot
plot_ch4_season_CA<-ch4_casepilot_comparisons %>%
  filter(casepilot == "CA", comparison == "season") %>%
  #ggplot inside brackets to be able to access casepilot name with ggtitle
  { 
    ggplot(., aes(x = season, y = emmean_gwp, label = group_letter)) +
      geom_errorbar(aes(ymin = emmean_gwp - SE_gwp, ymax = emmean_gwp + SE_gwp),
                    linewidth = 0.2, width = 0.2) +
      geom_point(shape = 15, size = 4, aes(col=season)) +
      geom_text(nudge_y = 0.4, color = "black") +
      geom_hline(yintercept=0, linetype="dashed")+
      theme_bw() +
      theme(axis.title = element_text(face = "bold"),
            axis.text = element_text(face = "bold"),
            plot.caption = element_text(hjust = 0)) +
      labs(title = expression(Seasonal~effect~on~net~CH[4]~flux),
           subtitle = paste0("Camargue"), 
           caption = paste0("Seasonal effect is averaged across the three conservation status.\nBoxes indicate the EM mean, error bars indicate SE,\n",              "Letters indicate different groups"),
           y=expression(EMM~(g~CO[2~eq]~d^-1~m^-2)),
           x="Sampling campaign",
           col=paste0("Sampling"))
  }

plot_ch4_season_CA


##Interaction effect plots----
#To do, unclear interest for overview ppt (maybe interesting for some casepilots)

#Example: CA
#Camargue: interaction plot
plot_ch4_interaction_CA<-ch4_casepilot_comparisons %>%
  filter(casepilot == "CA", comparison == "interaction") %>%
  #ggplot inside brackets to be able to access casepilot name with ggtitle
  { 
    #Define dodge width for separating status.
    dodge <- position_dodge(width = 0.6)
    
    ggplot(., aes(x = season, y = emmean_gwp, label = group_letter, group=status)) +
      geom_errorbar(aes(ymin = emmean_gwp - SE_gwp, ymax = emmean_gwp + SE_gwp),
                    linewidth = 0.2, width = 0.2, position = dodge) +
      geom_point(shape = 15, size = 4, aes(col=status), position = dodge) +
      geom_text(vjust = -2, color = "black", position = dodge) +
      geom_hline(yintercept=0, linetype="dashed")+
      theme_bw() +
      theme(axis.title = element_text(face = "bold"),
            axis.text = element_text(face = "bold"),
            plot.caption = element_text(hjust = 0)) +
      labs(title = expression(Seasonal~effect~on~net~CO[2]~flux),
           subtitle = paste0("Camargue"), 
           caption = paste0("Boxes indicate the EM mean, error bars indicate SE,\n",              "Letters indicate different groups"),
           y=expression(EMM~(g~CO[2~eq]~d^-1~m^-2)),
           x="Sampling campaign",
           col=paste0("Sampling"))
  }

plot_ch4_interaction_CA



#______________________-----

#CO2+CH4_______ ------
#IMPORTANT, crossing of NAs for Co2 and Ch4 might cause strata distributions to become biased or non-balanced between subsites/status in different casepilots. Especially critical for systems with a lot of ebullition (which often causes artefacts in CO2 leading to NA fluxes)

#Check normality on CO2+CH4: datasets:
table_all<- tibble()
for (cp in unique(data4models$casepilot)) {
  for (ghg in c("gwp_co2andch4")) { #"gwp_ch4","gwp_n2o","gwp_co2andch4"
    x<- data4models %>% filter(casepilot==cp&ghgspecies==ghg) %>% filter(!is.na(dailyflux)) %>% pull(dailyflux)
    table_trans<- evaluate_transformations(x)
    table_trans$casepilot<- cp
    table_trans$ghgspecies<- ghg
    table_all<- bind_rows(table_all, table_trans)
    print(plot_histograms(x)+ggtitle(paste(cp, ghg)))
    print(plot_qq(x)+ggtitle(paste(cp, ghg)))
  }
}

#get most-normal transformation for each case
best_trans<- table_all %>% 
  group_by(ghgspecies,casepilot) %>% 
  filter(Shapiro_p==max(Shapiro_p)) %>% 
  arrange(ghgspecies,casepilot)
best_trans


#CASEPILOT models-------
#STEPS: 
#1. select best transformation (using lmer, no sampleweights)
#2. Optimize model without sampleweights (lmer as first step --> lme (normal but heterocedasticity) or glmm (non-normal, heterocedasitcity or not))
#3. (skip)Optimize model with sampleweights(glmm)



#1. Best common transformation ------ 
#Focus on normality of residuals, homocedasticity of residuals, model fit, for multiple transformations.

# Initialize summary table
co2andch4_lmer_results <- tibble()

# Loop over all combinations of casepilot and CO2 using lmer models and residuals to check. 
for (cp in unique(data4models$casepilot)) {
  for (ghg in c("gwp_co2andch4")) {
    
    cat("\n\n==================== CASEPILOT:", cp, "| GHG:", ghg, "====================\n")
    
    data_sub <- data4models %>% filter(grepl(cp, casepilot)& ghgspecies == ghg) %>% 
      filter(!is.na(dailyflux))
    
    #Transformations
    data_sub <- data_sub %>%
      mutate(
        dailyflux_logsign = logsign(dailyflux),
        dailyflux_rootsign = rootsign(dailyflux),
        dailyflux_arcsinhsign = arcsinhsign(dailyflux)
      )
    
    # Step 3: Fit lmer models
    model_untrans <- tryCatch(lmer(dailyflux ~ status * season + (1 | subsite), data = data_sub), error = function(e) NULL)
    model_logsign <- tryCatch(lmer(dailyflux_logsign ~ status * season + (1 | subsite), data = data_sub), error = function(e) NULL)
    model_rootsign <- tryCatch(lmer(dailyflux_rootsign ~ status * season + (1 | subsite), data = data_sub), error = function(e) NULL)
    model_arcsinhsign<- tryCatch(lmer(dailyflux_arcsinhsign ~ status * season + (1 | subsite), data = data_sub), error = function(e) NULL)
    
    
    # Step 4: Diagnostics
    check_model <- function(model, name) {
      if (is.null(model)) return(NULL)
      res <- residuals(model)
      shapiro_res <- shapiro.test(res)$p.value
      levene_status <- tryCatch(leveneTest(res ~ data_sub$status)$"Pr(>F)"[1], error = function(e) NA)
      levene_season <- tryCatch(leveneTest(res ~ data_sub$season)$"Pr(>F)"[1], error = function(e) NA)
      r2 <- tryCatch(performance::r2_nakagawa(model), error = function(e) NULL)
      aic_val <- tryCatch(AIC(model), error = function(e) NA)
      bic_val <- tryCatch(BIC(model), error = function(e) NA)
      rmse_cv <- tryCatch({
        ctrl <- trainControl(method = "cv", number = 10)
        cv_model <- train(
          x = model.matrix(model)[, -1],
          y = getME(model, "y"),
          method = "lm",
          trControl = ctrl
        )
        cv_model$results$RMSE
      }, error = function(e) NA)
      
      tibble(
        casepilot = cp,
        gas = ghg,
        model = name,
        shapiro_resid_p = shapiro_res,
        levene_status_p = levene_status,
        levene_season_p = levene_season,
        marginal_r2 = if (!is.null(r2)) r2$R2_marginal else NA,
        conditional_r2 = if (!is.null(r2)) r2$R2_conditional else NA,
        AIC = aic_val,
        BIC = bic_val,
        RMSE_CV = rmse_cv
      )
    }
    
    res_untrans <- check_model(model_untrans, "Untransformed")
    res_logsign <- check_model(model_logsign, "Signed Log")
    res_rootsign <- check_model(model_rootsign, "Signed Root")
    res_arcsinhsign <- check_model(model_arcsinhsign, "Signed Arcsinh")
    
    co2andch4_lmer_results <- bind_rows(co2andch4_lmer_results, res_untrans, res_logsign, res_rootsign, res_arcsinhsign)
  }
}

co2andch4_lmer_results

#Untransfromed data fails for DU and VA (singularities)
co2andch4_lmer_results %>% 
  filter(is.na(conditional_r2))

#Normality of residuals never achieved
co2andch4_lmer_results %>% 
  filter(shapiro_resid_p>0.05)

#Homocedasticity for both factors is only for untransformed RI 
co2andch4_lmer_results %>% 
  filter(levene_season_p>0.05&levene_status_p>0.05)

#Best transformation is consistently Signed Log:
co2andch4_lmer_results %>% 
  group_by(casepilot,gas) %>% 
  filter(AIC==min(AIC))


#CONCLUSION: Lmer is not appropriate for modeling casepilots individually as is. We need to account for heterocedasticity. Heterocedasticity: We could use nlme::lme() and weights (syntax:  "weights = varIdent(form = ~1 | season)" ). BUT we have some non-normal residuals--> glmm for all (allows different family distributions and heterocedasticity)  


#We will use logsign-transformed data, it is the best transformation according to earlier tests with lmer (maximizes normality of residuals and minimizes heterocedasticity). 


#2. Optimice GLMM for each casepilot-------

#casepilots show strong non-normality... most show heterocedasticity for season. 

#Which distribution family is appropriate for each casepilot? Do we have heterocedasticity?
#I must account for heterocedasticity while testing the distribution family (it will change the residuals). 

#Notes on singularity: to have a model with singularity means that some of the random effect terms likely have variance estimates of (almost) zero, or there is redundancy in the random structure. IN our case, the random effect is solely the subsite (to account for subsite-specific differences and repeated samplings across seasons). Our random effect structure is the most simple possible (one random effect for intercept), and we expect this random effect to account for variance when the 2 subsites of each status are fundamentally different. With this structure, we should not be too concerned about singuarity. The inclusion of the subsite as a random effect might not be statistically necesary (variance~0) but it is still conceptually justified.  



## 2.1. CA model ----
#CA: lmer showed non-normal residuals, heterocedasticity for season (not for status)
co2andch4_lmer_results %>% filter(casepilot=="CA")
#approach: compare fit between gausian models, without heterocedasticity, with for season only and with for season and status. 

ca_co2andch4<- data4models %>% filter(casepilot=="CA"&ghgspecies=="gwp_co2andch4") %>% filter(!is.na(dailyflux)) %>% 
  mutate(dailyflux_logsign=logsign(dailyflux))

#Check data distribution: 
plot_histograms(ca_co2andch4$dailyflux)


ca_glmm_gaus0<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite),
                        data = ca_co2andch4,
                        family = gaussian(),
                        dispformula = ~1 #default, homocedasticity assumed
)

ca_glmm_gaus1<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite),
                        data = ca_co2andch4,
                        family = gaussian(),
                        dispformula = ~season #heterocedasticity across season levels (not status level)
)

ca_glmm_gaus2<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite),
                        data = ca_co2andch4,
                        family = gaussian(),
                        dispformula = ~season+status #heterocedasticity across season and status (no across interaction)
)

#Get residuals for all options
res0<- simulateResiduals(ca_glmm_gaus0)
res1<- simulateResiduals(ca_glmm_gaus1)
res2<- simulateResiduals(ca_glmm_gaus2)

#Explore residuals of model options
plot(res0) #check normality, 
plotQQunif(res0)
plotResiduals(res0, ca_co2andch4$season) #check heterocedasticity
plotResiduals(res0, ca_co2andch4$status) #check heterocedasticity
testDispersion(res0) #checks wether residual variance is higher or lower than expected, n.s.--> no issue. 
testOutliers(res0) #Checks if there are more extreme residuals than expected
#Very good residuals, heterocedasticity at season


plot(res1) #check normality, 
plotQQunif(res0)
plotResiduals(res1, ca_co2andch4$season) #check heterocedasticity
plotResiduals(res1, ca_co2andch4$status) #check heterocedasticity
testDispersion(res1) #checks wether residual variance is higher or lower than expected, n.s.--> no issue. 
testOutliers(res1) #Checks if there are more extreme residuals than expected
#Very good residuals, all tests pass

plot(res2) #check normality, 
plotResiduals(res2, ca_co2andch4$season) #check heterocedasticity
plotResiduals(res2, ca_co2andch4$status) #check heterocedasticity
testDispersion(res2) #checks wether residual variance is higher or lower than expected, n.s.--> no issue. 
testOutliers(res2) #Checks if there are more extreme residuals than expected
hist(res2$scaledResiduals)
#Very good residuals, all tests pass

#COMPARE model options
anova(ca_glmm_gaus0,ca_glmm_gaus1,ca_glmm_gaus2)

#CA best model: gaussian, allowing heterocedasticity across seasons. 
ca_best_co2andch4<- ca_glmm_gaus1
Anova(ca_best_co2andch4,type = 3) 

#Cannot calculate R2 (conditional & marginal) for models with heterocedasticity,
#Get PseudoR2
pseudoR2_ca_co2andch4<-pseudoR2_from_null(ca_best_co2andch4, ca_co2andch4)
pseudoR2_ca_co2andch4

#R2s of homocedastic model
ca_best_homoc_co2andch4<- ca_glmm_gaus0
r2(ca_best_homoc_co2andch4)


rm(ca_glmm_gaus0,ca_glmm_gaus1,ca_glmm_gaus2)

## 2.2. CU model ----
#CU: lmer showed non-normal residuals, heterocedasticity for season, and extremely low r2 (~0.06)
co2andch4_lmer_results %>% filter(casepilot=="CU")

#approach: compare fit between gausian and tstudent models without heterocedasticity, with for season only and with for season and status (not for interaction). 

#CU: non-normal residuals
cu_co2andch4<- data4models %>% filter(casepilot=="CU"&ghgspecies=="gwp_co2andch4") %>% filter(!is.na(dailyflux)) %>% 
  mutate(dailyflux_logsign=logsign(dailyflux))

#Check distribution of data: Bimodal distribution...Two near-normal groups (probably OW and V)
plot_histograms(cu_co2andch4$dailyflux)

cu_glmm_gaus0<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite), 
                        data = cu_co2andch4,
                        family = gaussian(),
                        dispformula = ~1)

cu_glmm_gaus1<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite), 
                        data = cu_co2andch4,
                        family = gaussian(),
                        dispformula = ~season)

cu_glmm_gaus2<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite), 
                        data = cu_co2andch4,
                        family = gaussian(),
                        dispformula = ~season + status)

cu_glmm_stu0<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite), 
                       data = cu_co2andch4,
                       family = t_family,
                       dispformula = ~1)

cu_glmm_stu1<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite), 
                       data = cu_co2andch4,
                       family = t_family,
                       dispformula = ~season)

cu_glmm_stu2<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite), 
                       data = cu_co2andch4,
                       family = t_family,
                       dispformula = ~season + status)

resgaus0<- simulateResiduals(cu_glmm_gaus0, n=1000)
resgaus1<- simulateResiduals(cu_glmm_gaus1, n=1000)
resgaus2<- simulateResiduals(cu_glmm_gaus2, n=1000)
restu0<- simulateResiduals(cu_glmm_stu0, n=1000)
restu1<- simulateResiduals(cu_glmm_stu1, n=1000)
restu2<- simulateResiduals(cu_glmm_stu2, n=1000)

plot(resgaus0) #check normality, 
plotQQunif(resgaus0)#plot QQ
plotResiduals(resgaus0, cu_co2andch4$season) #check heterocedasticity
plotResiduals(resgaus0, cu_co2andch4$status) #check heterocedasticity
#Non-normal residuals

plot(resgaus1) #check normality, 
plotResiduals(resgaus1, cu_co2andch4$season) #check heterocedasticity
plotResiduals(resgaus1, cu_co2andch4$status) #check heterocedasticity
#Non-normal

plot(restu0) #check normality, 
plotQQunif(restu0)#plot QQ
plotResiduals(restu0, cu_co2andch4$season) #check heterocedasticity
plotResiduals(restu0, cu_co2andch4$status) #check heterocedasticity
#Good residuals but heterocedasticity

plot(restu1) #check normality, 
plotQQunif(restu1)#plot QQ
plotResiduals(restu1, cu_co2andch4$season) #check heterocedasticity
plotResiduals(restu1, cu_co2andch4$status) #check heterocedasticity
#Good residuals, still heterocedasticity for status

plot(restu2) #check normality, 
plotQQunif(restu2)#plot QQ
plotResiduals(restu2, cu_co2andch4$season) #check heterocedasticity
plotResiduals(restu2, cu_co2andch4$status) #check heterocedasticity
#ALL ASSUMPITONS MET

testDispersion(restu2) #checks wether residual variance is higher or lower than expected, n.s.--> no issue. 
testOutliers(restu2) #Checks if there are more extreme residuals than expected

#ANOVA to decide: 
anova(cu_glmm_gaus0,cu_glmm_gaus1,cu_glmm_gaus2,cu_glmm_stu0,cu_glmm_stu1,cu_glmm_stu2)
#Student distribution is consistently better 
anova(cu_glmm_stu0,cu_glmm_stu1,cu_glmm_stu2)
#Big improvemnt with season heterocedasticity, best with season and status.
anova(cu_glmm_stu1,cu_glmm_stu2)


#CU best model: student, allowing heterocedasticity across seasons and status. 

cu_best_co2andch4<- cu_glmm_stu2
Anova(cu_best_co2andch4,type = 3) 


#Cannot calculate R2 (conditional & marginal) for models with heterocedasticity,
#Get PseudoR2
pseudoR2_cu_co2andch4<-pseudoR2_from_null(cu_best_co2andch4, cu_co2andch4)
pseudoR2_cu_co2andch4

#R2 of homocedastic model
cu_best_homoc_co2andch4<- cu_glmm_stu0
r2(cu_best_homoc_co2andch4)
MuMIn::r.squaredGLMM(cu_best_homoc_co2andch4)

#NO singularity:
performance::check_singularity(cu_best_co2andch4)
performance::check_singularity(cu_best_homoc_co2andch4)

rm(cu_glmm_gaus0,cu_glmm_gaus1,cu_glmm_gaus2,cu_glmm_stu0,cu_glmm_stu1,cu_glmm_stu2)




## 2.3. DA model ----
#DA: lmer showed non-normal residuals, heterocedasticity for season (and status),
co2andch4_lmer_results %>% filter(casepilot=="DA")
#approach: compare fit between gausian and tstudent models without heterocedasticity, with for season only and with for season and status (not for interaction). 

#DA: non-normal residuals
da_co2andch4<- data4models %>% filter(casepilot=="DA"&ghgspecies=="gwp_co2andch4") %>% filter(!is.na(dailyflux)) %>% 
  mutate(dailyflux_logsign=logsign(dailyflux))

plot_histograms(da_co2andch4$dailyflux) #Bimodal distribution for transformations

#calculate all model options
da_glmm_gaus0<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite), 
                        data = da_co2andch4,
                        family = gaussian(),
                        dispformula = ~1)

da_glmm_gaus1<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite), 
                        data = da_co2andch4,
                        family = gaussian(),
                        dispformula = ~season)

da_glmm_gaus2<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite), 
                        data = da_co2andch4,
                        family = gaussian(),
                        dispformula = ~season + status)

da_glmm_stu0<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite), 
                       data = da_co2andch4,
                       family = t_family,
                       dispformula = ~1)

da_glmm_stu1<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite), 
                       data = da_co2andch4,
                       family = t_family,
                       dispformula = ~season)

da_glmm_stu2<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite), 
                       data = da_co2andch4,
                       family = t_family,
                       dispformula = ~season + status)

#Get residuals for all model options
resgaus0<- simulateResiduals(da_glmm_gaus0, n=1000)
resgaus1<- simulateResiduals(da_glmm_gaus1, n=1000)
resgaus2<- simulateResiduals(da_glmm_gaus2, n=1000)
restu0<- simulateResiduals(da_glmm_stu0, n=1000)
restu1<- simulateResiduals(da_glmm_stu1, n=1000)
restu2<- simulateResiduals(da_glmm_stu2, n=1000)

plot(resgaus0) #check normality, 
plotResiduals(resgaus0, da_co2andch4$season) #check heterocedasticity
plotResiduals(resgaus0, da_co2andch4$status) #check heterocedasticity
#NOn-normal

plot(resgaus1) #check normality, 
plotResiduals(resgaus1, da_co2andch4$season) #check heterocedasticity
plotResiduals(resgaus1, da_co2andch4$status) #check heterocedasticity
#Non-normal

plot(restu0) #check normality, 
plotResiduals(restu0, da_co2andch4$season) #check heterocedasticity
plotResiduals(restu0, da_co2andch4$status) #check heterocedasticity

plot(restu1) #check normality, 
plotResiduals(restu1, da_co2andch4$season) #check heterocedasticity
plotResiduals(restu1, da_co2andch4$status) #check heterocedasticity
#All asumptions met

plot(restu2) #check normality, 
plotResiduals(restu2, da_co2andch4$season) #check heterocedasticity
plotResiduals(restu2, da_co2andch4$status) #check heterocedasticity


testDispersion(restu1) #checks wether residual variance is higher or lower than expected, n.s.--> no issue. 
testOutliers(restu1) #Checks if there are more extreme residuals than expected

#ANOVA to decide: 
anova(da_glmm_gaus0,da_glmm_gaus1,da_glmm_gaus2,da_glmm_stu0,da_glmm_stu1,da_glmm_stu2)
#Student distribution is better always
anova(da_glmm_stu0,da_glmm_stu1,da_glmm_stu2)
#Big improvemnt with season heterocedasticity, No improvement with additional status heterocedasticity. 
anova(da_glmm_stu1,da_glmm_stu2)


#DA best model: student, allowing heterocedasticity across seasons only. 

da_best_co2andch4<- da_glmm_stu1
Anova(da_best_co2andch4,type = 3)

#Cannot calculate R2 (conditional & marginal) for models with heterocedasticity,
#Get PseudoR2
pseudoR2_da_co2andch4<-pseudoR2_from_null(da_best_co2andch4, da_co2andch4)
pseudoR2_da_co2andch4


#R2 of homocedastic model
da_best_homoc_co2andch4<- da_glmm_stu0
r2(da_best_homoc_co2andch4)
MuMIn::r.squaredGLMM(da_best_homoc_co2andch4) #random effect R2 is 0 (singularity), R2m==R2c

#CHECK singularity
performance::check_singularity(da_glmm_stu0) #Singularity for T-family
performance::check_singularity(da_glmm_gaus0) #Not for gaussian

#DOUBT singularity DA:
#This means that the two subsites of each status behave similarly (they are good "replicates"). Even when is not statistically needed, we include it for consistency and for conceptual reasons. Singularity does not impact interpretation but will limit some functionalities (r2). 
#We expected A1 and A2 subsites to have different CO2andch4 profiles and this to be absorved by the random effect, but apparently not so much (subsite explains no variance)
rm(da_glmm_gaus0,da_glmm_gaus1,da_glmm_gaus2,da_glmm_stu0,da_glmm_stu1,da_glmm_stu2)

## 2.4. DU model ----
#DU: lmer showed near-normal residuals, heterocedasticity for season (not for status)
co2andch4_lmer_results %>% filter(casepilot=="DU")
#approach: compare fit between gausian models: without heterocedasticity, with for season only and with for season and status. 

du_co2andch4<- data4models %>% filter(casepilot=="DU"&ghgspecies=="gwp_co2andch4") %>% filter(!is.na(dailyflux)) %>% 
  mutate(dailyflux_logsign=logsign(dailyflux))

plot_histograms(du_co2andch4$dailyflux)


du_glmm_gaus0<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite),
                        data = du_co2andch4,
                        family = gaussian(),
                        dispformula = ~1 #default, homocedasticity assumed
)

du_glmm_gaus1<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite),
                        data = du_co2andch4,
                        family = gaussian(),
                        dispformula = ~season #heterocedasticity across season levels (not status level)
)

du_glmm_gaus2<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite),
                        data = du_co2andch4,
                        family = gaussian(),
                        dispformula = ~season+status #heterocedasticity across season and status (no across interaction)
)

#Get residuals for all options
res0<- simulateResiduals(du_glmm_gaus0, n=1000)
res1<- simulateResiduals(du_glmm_gaus1, n=1000)
res2<- simulateResiduals(du_glmm_gaus2, n=1000)

#Explore residuals of model options
plot(res0) #check normality, 
plotResiduals(res0, du_co2andch4$season) #check heterocedasticity
plotResiduals(res0, du_co2andch4$status) #check heterocedasticity

plot(res1) #check normality, 
plotResiduals(res1, du_co2andch4$season) #check heterocedasticity
plotResiduals(res1, du_co2andch4$status) #check heterocedasticity
testDispersion(res1) #checks wether residual variance is higher or lower than expected, n.s.--> no issue. 
testOutliers(res1) #Checks if there are more extreme residuals than expected
hist(res1$scaledResiduals)


#COMPARE model options: Heterocedasticity for season improves model performance. , but improves AIC significantly
anova(du_glmm_gaus0,du_glmm_gaus1,du_glmm_gaus2)

#DU best model: gaussian, allowing heterocedasticity across seasons. 
du_best_co2andch4<- du_glmm_gaus1
Anova(du_best_co2andch4,type = 3)

#Cannot calculate R2 (conditional & marginal) for models with heterocedasticity,
#Get PseudoR2
pseudoR2_du_co2andch4<-pseudoR2_from_null(du_best_co2andch4, du_co2andch4)
pseudoR2_du_co2andch4

#R2 of homocedastic model
du_best_homoc_co2andch4<-du_glmm_gaus0
r2(du_best_homoc_co2andch4)

#Remove non-best models
rm(du_glmm_gaus0,du_glmm_gaus1,du_glmm_gaus2)


## 2.5. RI model ----
#RI: lmer showed normal residuals, heterocedasticity for status (but not for season)
co2andch4_lmer_results %>% filter(casepilot=="RI")
#approach: compare fit between gausian models: without heterocedasticity, with for season only and with for season and status. 
##___TO_TRY-------
#Compare logsign transformed (gaussian) and untransformed (t-family) models, use the one with best fit (AIC focus only)
#CO2+CH4: limits = c(-7.5,2.5), (removes 3 outliers)

#Subset Co2 data for casepilot
ri_co2andch4<- data4models %>% filter(casepilot=="RI"&ghgspecies=="gwp_co2andch4") %>% filter(!is.na(dailyflux)) %>% 
  mutate(dailyflux_logsign=logsign(dailyflux)) %>% 
  filter(between(dailyflux, -7.5,2.5))

plot_histograms(ri_co2andch4$dailyflux)

#Fit model options
ri_glmm_gaus0<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite),
                        data = ri_co2andch4,
                        family = gaussian(),
                        dispformula = ~1 #default, homocedasticity assumed
)

ri_glmm_gaus1<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite),
                        data = ri_co2andch4,
                        family = gaussian(),
                        dispformula = ~status #heterocedasticity across status levels (not seasons level)
)

ri_glmm_gaus2<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite),
                        data = ri_co2andch4,
                        family = gaussian(),
                        dispformula = ~season+status #heterocedasticity across season and status (no across interaction)
)

#Get residuals for all options
res0<- simulateResiduals(ri_glmm_gaus0, n=1000)
res1<- simulateResiduals(ri_glmm_gaus1, n=1000)
res2<- simulateResiduals(ri_glmm_gaus2, n=1000)

#Explore residuals of model options
plot(res0) #check normality, 
plotResiduals(res0, ri_co2andch4$season) #check heterocedasticity
plotResiduals(res0, ri_co2andch4$status) #check heterocedasticity
testDispersion(res0) #checks wether residual variance is higher or lower than expected, n.s.--> no issue. 
#Normal but heterocedasticity for status


plot(res1) #check normality, 
plotResiduals(res1, ri_co2andch4$season) #check heterocedasticity
plotResiduals(res1, ri_co2andch4$status) #check heterocedasticity
testDispersion(res1) #checks wether residual variance is higher or lower than expected, n.s.--> no issue. 
testOutliers(res1) #Checks if there are more extreme residuals than expected
hist(res1$scaledResiduals)
#ALL ASUMPTIONS MET

#COMPARE model options: Slight heterocedasticity across status, no significant improvement when modelled.
anova(ri_glmm_gaus0,ri_glmm_gaus1,ri_glmm_gaus2)

#RI best model: gaussian, assuming homocedasticity
ri_best_co2andch4<- ri_glmm_gaus1
Anova(ri_best_co2andch4,type = 3)

#Set homoc model as best model 
ri_best_homoc_co2andch4<- ri_glmm_gaus0

#Explained variance: calculate R2 (conditional & marginal)
r2(ri_best_homoc_co2andch4)


#For consistency, calculate pseudo-R2 (log-Likelyhood based):
#Get PseudoR2
pseudoR2_ri_co2andch4<-pseudoR2_from_null(ri_best_co2andch4, ri_co2andch4)
pseudoR2_ri_co2andch4


#Remove non-best models
rm(ri_glmm_gaus0,ri_glmm_gaus1,ri_glmm_gaus2)

## 2.6. VA model ----
#VA: lmer showed non-normal residuals, heterocedasticity for season (and status?),
co2andch4_lmer_results %>% filter(casepilot=="VA")
#approach: compare fit between gausian and tstudent models without heterocedasticity, with for season only and with for season and status (not for interaction). 



#VA: non-normal residuals
va_co2andch4<- data4models %>% filter(casepilot=="VA"&ghgspecies=="gwp_co2andch4") %>% filter(!is.na(dailyflux)) %>% 
  mutate(dailyflux_logsign=logsign(dailyflux))

plot_qq(va_co2andch4$dailyflux)
plot_histograms(va_co2andch4$dailyflux) #BIMODAL distribution (sligthly)

#calculate all model options
va_glmm_gaus0<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite), 
                        data = va_co2andch4,
                        family = gaussian(),
                        dispformula = ~1)

va_glmm_gaus1<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite), 
                        data = va_co2andch4,
                        family = gaussian(),
                        dispformula = ~season)

va_glmm_gaus2<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite), 
                        data = va_co2andch4,
                        family = gaussian(),
                        dispformula = ~season + status)

va_glmm_stu0<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite), 
                       data = va_co2andch4,
                       family = t_family,
                       dispformula = ~1)

va_glmm_stu1<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite), 
                       data = va_co2andch4,
                       family = t_family,
                       dispformula = ~season)

va_glmm_stu2<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite), 
                       data = va_co2andch4,
                       family = t_family,
                       dispformula = ~season + status)

#Get residuals for all model options
resgaus0<- simulateResiduals(va_glmm_gaus0, n=1000)
resgaus1<- simulateResiduals(va_glmm_gaus1, n=1000)
resgaus2<- simulateResiduals(va_glmm_gaus2, n=1000)
restu0<- simulateResiduals(va_glmm_stu0, n=1000)
restu1<- simulateResiduals(va_glmm_stu1, n=1000)
restu2<- simulateResiduals(va_glmm_stu2, n=1000)

plot(resgaus0) #check normality, 
plotResiduals(resgaus0, va_co2andch4$season) #check heterocedasticity
plotResiduals(resgaus0, va_co2andch4$status) #check heterocedasticity
#Non-normal

plot(resgaus1) #check normality, 
plotResiduals(resgaus1, va_co2andch4$season) #check heterocedasticity
plotResiduals(resgaus1, va_co2andch4$status) #check heterocedasticity
#Better normality, but heterocedasticity

plot(resgaus2) #check normality, 
plotResiduals(resgaus2, va_co2andch4$season) #check heterocedasticity
plotResiduals(resgaus2, va_co2andch4$status) #check heterocedasticity
#Better normality, but heterocedasticity remains

plot(restu0) #check normality, 
plotResiduals(restu0, va_co2andch4$season) #check heterocedasticity
plotResiduals(restu0, va_co2andch4$status) #check heterocedasticity
#MODEL DID NOT CONVERGE

plot(restu1) #check normality, 
plotResiduals(restu1, va_co2andch4$season) #check heterocedasticity
plotResiduals(restu1, va_co2andch4$status) #check heterocedasticity
#Better residuals, heterocedasticity

plot(restu2) #check normality, 
plotResiduals(restu2, va_co2andch4$season) #check heterocedasticity
plotResiduals(restu2, va_co2andch4$status) #check heterocedasticity


testDispersion(resgaus1) #checks wether residual variance is higher or lower than expected, n.s.--> no issue. 
testOutliers(resgaus1) #Checks if there are more extreme residuals than expected
hist(resgaus1$scaledResiduals)

#ANOVA to decide: 
anova(va_glmm_gaus0,va_glmm_gaus1,va_glmm_gaus2,va_glmm_stu0,va_glmm_stu1,va_glmm_stu2)
#Student distribution does not improve fit
anova(va_glmm_gaus0,va_glmm_gaus1,va_glmm_gaus2)
#Big improvemnt with season heterocedasticity, limited improvement with additional status heterocedasticity. 
anova(va_glmm_gaus1,va_glmm_gaus2)


#VA best model: gaus, allowing heterocedasticity across seasons only. 

va_best_co2andch4<- va_glmm_gaus1
Anova(va_best_co2andch4,type = 3)

#Cannot calculate R2 (conditional & marginal) for models with heterocedasticity,
#Get PseudoR2
pseudoR2_va_co2andch4<-pseudoR2_from_null(va_best_co2andch4, va_co2andch4)
pseudoR2_va_co2andch4


##___DOUBT low R2 for VA------
#R2 of homocedastic model
va_best_homoc_co2andch4<- va_glmm_gaus0
r2(va_best_homoc_co2andch4)
#Low R2, no status effect detected

#Check singularity: 
performance::check_singularity(va_glmm_gaus0)#No singularity
performance::check_singularity(va_glmm_gaus1)#No singularity

#Remove non-best models
rm(va_glmm_gaus0,va_glmm_gaus1,va_glmm_gaus2,va_glmm_stu0,va_glmm_stu1,va_glmm_stu2)


#3.(skip) Optimize with weights--------
#TO ADAPT _co2andch4
#Here we will create glmms taking into account strata sample_weights. 
#Then compare models with and without sample_weights to see the effect of weights and get the overall best model for each casepilot. 

# DOUBT: should I directly implement the identified best_model option for each casepilot adding sample_weight OR optimize from the beginning (allowing different heterocedasticity structures)??:

#Applying weights to the already optimized structure can tell us directly if sample weights improve fit. 
#Optimizing the models from the beginning with sample weights will allow us to check whether the weights influence the overall structure selection BUT we will not be able to directly compare with anova() if models differ in structure (i.e. if they have different fixed/random effects, different dispformula, or different distribution families) 

#We know that we should account for sample weights in the end if/when we get accurate strata composition and if they show substantial and consistent sampling bias. 

## 3.1. CA model (weighted)

#Do sample weights improve the optimized model? 
#Recal best model structure: 
summary(ca_best)

#build the same model with weights
ca_best_weighted<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite),
                           data = ca_co2,
                           family = gaussian(),
                           dispformula = ~season, #heterocedasticity across season levels (not status level)
                           weights=sample_weight)

#Compare fit between weigthed and unweighted model
anova(ca_best, ca_best_weighted)
#inclusion of sample_weight does not significantly improve model. Unclear effect on Df but weights cause no improvement on AIC...

Anova(ca_best,type = 3)
Anova(ca_best_weighted,type = 3)
#Weights also do not change significance of main effects. 

#Do weights lead to different structure? #Run the optimization process to select best model (with/without heterocedasticity) Not essential to check, omit. 


## 3.2. CU model (weighted)

#Do sample weights improve the optimized model? 
#Recal best model structure: 
summary(cu_best)

#build the same model with weights
cu_best_weighted<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite),
                           data = cu_co2,
                           family = t_family,
                           dispformula = ~season, #heterocedasticity across season levels (not status level)
                           weights=sample_weight)

#Compare fit between weigthed and unweighted model
anova(cu_best, cu_best_weighted)
#inclusion of sample_weight does significantly improve model. Unclear effect on Df but weights cause  improvement on AIC...

Anova(cu_best,type = 3)
Anova(cu_best_weighted,type = 3)
#Weights make main effects sigltly more significant.
summary(cu_best_weighted)

#Do weights lead to different structure? #Run the optimization process to select best model (with/without heterocedasticity) Not essential to check, omit. 



## 3.3. DA model (weighted)

#Do sample weights improve the optimized model? 
#Recal best model structure: 
summary(da_best)

#build the same model with weights
da_best_weighted<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite),
                           data = da_co2,
                           family = t_family,
                           dispformula = ~season, #heterocedasticity across season levels (not status level)
                           weights=sample_weight)

#Compare fit between weigthed and unweighted model
anova(da_best, da_best_weighted)
#inclusion of sample_weight does NOT significantly improve model. Unclear effect on Df....

Anova(da_best,type = 3)
Anova(da_best_weighted, type = 3)
#Weights make main effects sigltly less significant.

#Do weights lead to different structure? #Run the optimization process to select best model (with/without heterocedasticity) Not essential to check, omit. 



## 3.4. DU model (weighted)

#Do sample weights improve the optimized model? 
#Recal best model structure: 
summary(du_best)

#build the same model with weights
du_best_weighted<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite),
                           data = du_co2,
                           family = gaussian(),
                           dispformula = ~season, #heterocedasticity across season levels (not status level)
                           weights=sample_weight)

#Compare fit between weigthed and unweighted model
anova(du_best, du_best_weighted)
#inclusion of sample_weight does NOT significantly improve model. Unclear effect on Df....

#No detectable effect of status
Anova(du_best,type = 3)
Anova(du_best_weighted, type = 3)

#Do weights lead to different structure? #Run the optimization process to select best model (with/without heterocedasticity) Not essential to check, omit. 

## 3.5. RI model (weighted)

#Do sample weights improve the optimized model? 
#Recal best model structure: 
summary(ri_best)

#build the same model with weights
ri_best_weighted<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite),
                           data = ri_co2,
                           family = gaussian(),
                           # dispformula = ~season, #HOMOCEDASTICITY assumed
                           weights=sample_weight)

#Compare fit between weigthed and unweighted model
anova(ri_best, ri_best_weighted)
#inclusion of sample_weight does NOT significantly improve model. Unclear effect on Df....


#Same significance of main effects
Anova(ri_best,type = 3)
Anova(ri_best_weighted,type = 3)

#Do weights lead to different structure? #Run the optimization process to select best model (with/without heterocedasticity) Not essential to check, omit. 


## 3.6. VA model (weighted)

#Do sample weights improve the optimized model? 
#Recal best model structure: 
summary(va_best)

#build the same model with weights
va_best_weighted<- glmmTMB(formula = dailyflux_logsign~status*season + (1|subsite),
                           data = va_co2,
                           family = gaussian(),
                           dispformula = ~season, #HOMOCEDASTICITY assumed
                           weights=sample_weight)

#Compare fit between weigthed and unweighted model
anova(va_best, va_best_weighted)
#inclusion of sample_weight does NOT significantly improve model. Unclear effect on Df....


#Same significance of main effects
Anova(va_best,type = 3)
Anova(va_best_weighted, type = 3)

#Do weights lead to different structure? #Run the optimization process to select best model (with/without heterocedasticity) Not essential to check, omit. 



#4. Effects of best-models-------

#Use best models without strata sample_weights for the moment. 

#From last step: only model in curonian lagoon benefits from sample weights, this is expected as we forced the strata distribution to be fixed across different samplings (most appropriate option to get overall effect). It is the only subsite where strata distribution has to be especially accounted for in comparisons. The rest of casepilots might benefit from weights once we have the actual distribution of them. Good thing is: results do not change much, so our incorporation of weights is not too critical (if criticized, it could be ommited and resuts would hold) 


#BASED ON THE ABOVE OPTIMIZATIONS, provide list of models to evaluate/interpret

#Provide list of homocedastic model variants (to estimate R2) THIS WILL BE PROBLEMATIC FOR MODELS WITH SINGULARITIES... 
homoced_model_list_co2andch4<- list(
  "CA" = ca_best_homoc_co2andch4,
  "CU" = cu_best_homoc_co2andch4, 
  "DA" = da_best_homoc_co2andch4,
  "DU" = du_best_homoc_co2andch4,
  "RI" = ri_best_homoc_co2andch4, #Only casepilot for which homocedastic model is best
  "VA" = va_best_homoc_co2andch4)

#List of best models for each casepilot (to get residual tests and Significances)
best_model_list_co2andch4<- list(
  "CA" = ca_best_co2andch4,
  "CU" = cu_best_co2andch4, #Only model for which sample_weight improved fit
  "DA" = da_best_co2andch4,
  "DU" = du_best_co2andch4,
  "RI" = ri_best_co2andch4,
  "VA" = va_best_co2andch4)

#List of dataframes used to build each model (needed to calculate residuals)
model_data_list_co2andch4<- list(
  "CA" = ca_co2andch4,
  "CU" = cu_co2andch4, 
  "DA" = da_co2andch4,
  "DU" = du_co2andch4,
  "RI" = ri_co2andch4,
  "VA" = va_co2andch4)

## 4.1. Overall fit -----
#1st how good is the model--> overall R2 for homocedastic model variant 

co2andch4_model_fit_summary<- summarize_model_info(homoced_model_list_co2andch4)%>% 
  mutate(ghgspecies="gwp_co2andch4")


## 4.2. Residual checks---
#Summary of Model assumptions: for best_model_list
#We want to structure the results in a table: model structure,  distribution family, dispformula, 
#homoc_marginal_r2, homoc_conditional_r2. add pvalue for: normality of residuals, heteroscedasticity (levene test).

#Get summary dharma diagnostics (using function created in co2 section)
co2andch4_resid_diag_summary <- summarize_dharma_diagnostics(best_model_list_co2andch4, model_data_list_co2andch4)%>% 
  mutate(ghgspecies="gwp_co2andch4")

print(co2andch4_resid_diag_summary)
#Residual tests are ok for most models with 2 exceptions: 
#Uniformity for DU
#Heterocedasticity for VA


##4.2. Significance of effects-----
#Effect of status will first be assessed via Anova (best_model,type=3): is there an effect averaging across all seasons? yes/no, is the effect dependent on season? 

#Obtain significance for main effects of best_models using function built in co2 section
co2andch4_results_anova<-summarize_anova_results(best_model_list_co2andch4) %>% 
  mutate(ghgspecies="gwp_co2andch4")
print(co2andch4_results_anova)

#Intercept refers to whether the grand mean (on the transformed scale) is significantly different than cero (Is the overall net exchange of CO2 balanced between uptake and emission)
#Status refers to whether the status (averaged across seasons) has a significant effect on flux. 
#Season refers to whether the season (averaged across status) has a significant effect on flux.
#Interaction refers to whether the effect of status depends on season — i.e., the difference between status levels changes across seasons.



#Look for significance at either status or status:season---> there is an effect (that may or may not be evident when averaging across seasons)
#IF status:season significant---> effect depends on season. 

#Does status have an effect (significance of status or interaction)
co2andch4_results_anova %>% filter(status_pval<0.05|interaction_pval<0.05)

# CA, DA and RI show significant effect of interaction (effect of status dependent on season)

#CU, DU and VA show no effect (main or interaction of status)

#Save in dropbox tables with the models summary outputs: 

write.csv(x = co2andch4_model_fit_summary, file = paste0(plots_path,"/CO2plusCH4models_homocedastic_R2.csv"),row.names = F)
write.csv(x = co2andch4_resid_diag_summary, file = paste0(plots_path,"/CO2plusCH4models_tests_residuals.csv"),row.names = F)
write.csv(x = co2andch4_results_anova, file = paste0(plots_path,"/CO2plusCH4models_effects_signififcance.csv"),row.names = F)

##4.3. Effect direction and magnitude----

#Use model scale for statistics, back-transform for interpretation.  

#NOTE on bootstrapping: there is the possibility to bootstrap the models (i.e. recalculate the model many times with re-sampled data, and calculate emmeans for each model, to be able to produce more robust Confidence Intervals for them, but this is not usually done and computationally intensive. We would have to take into account the model structure when re-sampling the data to ensure balanced resamplings (all factors represented by data). WE will not do this (for the moment).

#Some guides to plot and use emmeans: https://rcompanion.org/handbook/G_06.html


#Then emmeans comparisons to get direction of effect and magnitude (averaging over season interaction if any)
plot(emmeans(ri_best_co2andch4, ~status), comparisons=TRUE) #Unclear interpretation of arrows, not the best way of presenting it. 


#Notes on non-gaussian models: upper and lower confidence limits are given as asymp.LCL/asymp.UCL (in our context we can treat them equally as the ones from gaussian models). Additionally, non-gaussian models will have Inf df (degrees of freedom are calculated differently for T-family models), not an issue. 

# Initialize list for storing comparisons
comparisons_list <- list()

# Loop to extract all comparisons for every casepilot best model
for (casepilot_name in names(best_model_list_co2andch4)) {
  ghgspecies<- "gwp_co2andch4"
  cp_model <- best_model_list_co2andch4[[casepilot_name]]
  
  # Status comparisons:
  status_emmeans <- emmeans(cp_model, ~status)
  status_comparisons <- cld(status_emmeans, 
                            alpha = 0.05,
                            Letters = letters,
                            adjust = "sidak") %>%
    mutate(comparison = "status",
           ghgspecies = ghgspecies,
           casepilot = casepilot_name,
           season = NA)%>%    
    #Rename confidence limits for T_family models (if this is the case)
    rename(lower.CL = any_of("asymp.LCL"),
           upper.CL = any_of("asymp.UCL"))
  
  # Season comparisons: 
  season_emmeans <- emmeans(cp_model, ~season)
  season_comparisons <- cld(season_emmeans, 
                            alpha = 0.05,
                            Letters = letters,
                            adjust = "sidak") %>%
    mutate(comparison = "season",
           ghgspecies = ghgspecies,
           casepilot = casepilot_name,
           status = NA)%>%
    #Rename confidence limits for T_family models (if this is the case)
    rename(lower.CL = any_of("asymp.LCL"),
           upper.CL = any_of("asymp.UCL"))
  
  
  # Interaction comparisons: 
  interaction_emmeans <- emmeans(cp_model, ~status:season)
  interaction_comparisons <- cld(interaction_emmeans, 
                                 alpha = 0.05,
                                 Letters = letters,
                                 adjust = "sidak") %>%
    mutate(comparison = "interaction",
           ghgspecies = ghgspecies,
           casepilot = casepilot_name)%>%
    #Rename confidence limits for T_family models (if this is the case)
    rename(lower.CL = any_of("asymp.LCL"),
           upper.CL = any_of("asymp.UCL"))
  
  
  # Join and format CP comparisons
  all_comparisons <- status_comparisons %>%
    full_join(season_comparisons) %>%
    full_join(interaction_comparisons) %>%
    rename(group_letter = .group) %>%
    mutate(
      comparison = factor(comparison, levels = c("status", "season", "interaction")),
      season = factor(season, levels = c("S1", "S2", "S3", "S4")),
      status = factor(status, levels = c("Altered", "Preserved", "Restored")),
      group_letter = gsub(" ", "", group_letter),
      emmean_gwp = inv_logsign(emmean),
      SE_gwp = inv_logsign(SE),
      lower.CL_gwp = inv_logsign(lower.CL),
      upper.CL_gwp = inv_logsign(upper.CL)
    ) %>%
      dplyr::select(any_of(c(
      "ghgspecies", "casepilot", "comparison", "status", "season", "df",
      "emmean", "SE", "lower.CL", "upper.CL", "group_letter",
      "emmean_gwp", "SE_gwp", "lower.CL_gwp", "upper.CL_gwp"
    ))) %>%
    arrange(ghgspecies, casepilot, comparison, status, season)
  
  # Store casepilot results
  comparisons_list[[casepilot_name]] <- all_comparisons
}

# Unir todos los resultados en un solo data.frame
co2andch4_casepilot_comparisons <- bind_rows(comparisons_list)


#Save all comparisons:  
write.csv(x = co2andch4_casepilot_comparisons, file = paste0(plots_path,"/CO2plusCH4models_emmeans_posthoc.csv"),row.names = F)


#Plot casepilot-per-casepilot, emmeans +- SE/CL (in gwp scale)

##Plots Status effect-----

#Camargue: 
plot_co2andch4_status_CA <- co2andch4_casepilot_comparisons %>%
  filter(casepilot == "CA", comparison == "status") %>%
  ggplot(., aes(x = status, y = emmean_gwp, label = group_letter)) +
  geom_errorbar(aes(ymin = emmean_gwp - SE_gwp, ymax = emmean_gwp + SE_gwp),
                linewidth = 0.2, width = 0.2) +
  geom_point(shape = 15, size = 4, aes(col=status)) +
  geom_text(nudge_x = 0.2, color = "black") +
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw() +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  labs(title = expression(Status~effect~on~NEE~(CO[2]+CH[4])),
       subtitle = paste0("Camargue"), 
       y=expression(EMM~(g~CO[2~eq]~d^-1~m^-2)),
       x="Conservation status",
       col=paste0("Status"))

plot_co2andch4_status_CA
ggsave(plot=plot_co2andch4_status_CA, filename="CO2plusCH4_status_plot_CA.png",
       device = "png",path = plots_path,dpi = 400,
       width = 115, height = 100,units = "mm")


#Curonian Lagoon: 
plot_co2andch4_status_CU<- co2andch4_casepilot_comparisons %>%
  filter(casepilot == "CU", comparison == "status") %>%
  ggplot(., aes(x = status, y = emmean_gwp, label = group_letter)) +
  geom_errorbar(aes(ymin = emmean_gwp - SE_gwp, ymax = emmean_gwp + SE_gwp),
                linewidth = 0.2, width = 0.2) +
  geom_point(shape = 15, size = 4, aes(col=status)) +
  geom_text(nudge_x = 0.2, color = "black") +
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw() +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  labs(title = expression(Status~effect~on~NEE~(CO[2]+CH[4])),
       subtitle = paste0("Curonian lagoon"), 
       y=expression(EMM~(g~CO[2~eq]~d^-1~m^-2)),
       x="Conservation status",
       col=paste0("Status"))

plot_co2andch4_status_CU
ggsave(plot=plot_co2andch4_status_CU, filename="CO2plusCH4_status_plot_CU.png",
       device = "png",path = plots_path,dpi = 400,
       width = 115, height = 100,units = "mm")

#Danube delta: 
plot_co2andch4_status_DA<- co2andch4_casepilot_comparisons %>%
  filter(casepilot == "DA", comparison == "status") %>%
  ggplot(., aes(x = status, y = emmean_gwp, label = group_letter)) +
  geom_errorbar(aes(ymin = emmean_gwp - SE_gwp, ymax = emmean_gwp + SE_gwp),
                linewidth = 0.2, width = 0.2) +
  geom_point(shape = 15, size = 4, aes(col=status)) +
  geom_text(nudge_x = 0.2, color = "black") +
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw() +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  labs(title = expression(Status~effect~on~NEE~(CO[2]+CH[4])),
       subtitle = paste0("Danube delta"),
       y=expression(EMM~(g~CO[2~eq]~d^-1~m^-2)),
       x="Conservation status",
       col=paste0("Status"))

plot_co2andch4_status_DA
ggsave(plot=plot_co2andch4_status_DA, filename="CO2plusCH4_status_plot_DA.png",
       device = "png",path = plots_path,dpi = 400,
       width = 115, height = 100,units = "mm")

#Dutch delta: 
plot_co2andch4_status_DU<- co2andch4_casepilot_comparisons %>%
  filter(casepilot == "DU", comparison == "status") %>%
  ggplot(., aes(x = status, y = emmean_gwp, label = group_letter)) +
  geom_errorbar(aes(ymin = emmean_gwp - SE_gwp, ymax = emmean_gwp + SE_gwp),
                linewidth = 0.2, width = 0.2) +
  geom_point(shape = 15, size = 4, aes(col=status)) +
  geom_text(nudge_x = 0.2, color = "black") +
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw() +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  labs(title = expression(Status~effect~on~NEE~(CO[2]+CH[4])),
       subtitle = paste0("Dutch delta"),
       y=expression(EMM~(g~CO[2~eq]~d^-1~m^-2)),
       x="Conservation status",
       col=paste0("Status"))

plot_co2andch4_status_DU
ggsave(plot=plot_co2andch4_status_DU, filename="CO2plusCH4_status_plot_DU.png",
       device = "png",path = plots_path,dpi = 400,
       width = 115, height = 100,units = "mm")


#Ria d'Aveiro: 
plot_co2andch4_status_RI<- co2andch4_casepilot_comparisons %>%
  filter(casepilot == "RI", comparison == "status") %>%
  ggplot(., aes(x = status, y = emmean_gwp, label = group_letter)) +
  geom_errorbar(aes(ymin = emmean_gwp - SE_gwp, ymax = emmean_gwp + SE_gwp),
                linewidth = 0.2, width = 0.2) +
  geom_point(shape = 15, size = 4, aes(col=status)) +
  geom_text(nudge_x = 0.2, color = "black") +
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw() +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  labs(title = expression(Status~effect~on~NEE~(CO[2]+CH[4])),
       subtitle = paste0("Ria d'Aveiro"), 
       y=expression(EMM~(g~CO[2~eq]~d^-1~m^-2)),
       x="Conservation status",
       col=paste0("Status"))

plot_co2andch4_status_RI
ggsave(plot=plot_co2andch4_status_RI, filename="CO2plusCH4_status_plot_RI.png",
       device = "png",path = plots_path,dpi = 400,
       width = 115, height = 100,units = "mm")

#Valencian wetland: 
plot_co2andch4_status_VA<- co2andch4_casepilot_comparisons %>%
  filter(casepilot == "VA", comparison == "status") %>%
  ggplot(., aes(x = status, y = emmean_gwp, label = group_letter)) +
  geom_errorbar(aes(ymin = emmean_gwp - SE_gwp, ymax = emmean_gwp + SE_gwp),
                linewidth = 0.2, width = 0.2) +
  geom_point(shape = 15, size = 4, aes(col=status)) +
  geom_text(nudge_x = 0.2, color = "black") +
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw() +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  labs(title = expression(Status~effect~on~NEE~(CO[2]+CH[4])),
       subtitle = paste0("Valencian wetland"), 
       y=expression(EMM~(g~CO[2~eq]~d^-1~m^-2)),
       x="Conservation status",
       col=paste0("Status"))

plot_co2andch4_status_VA
ggsave(plot=plot_co2andch4_status_VA, filename="CO2plusCH4_status_plot_VA.png",
       device = "png",path = plots_path,dpi = 400,
       width = 115, height = 100,units = "mm")


##Plots Seasonal effect -----
#To do, little interest for overview ppt (seasonal effect is averaged across the three status)

#Example:
#Camargue: 
plot_co2_season_CA<- co2_casepilot_comparisons %>%
  filter(casepilot == "CA", comparison == "season") %>%
  ggplot(., aes(x = season, y = emmean_gwp, label = group_letter)) +
  geom_errorbar(aes(ymin = emmean_gwp - SE_gwp, ymax = emmean_gwp + SE_gwp),
                linewidth = 0.2, width = 0.2) +
  geom_point(shape = 15, size = 4, aes(col=season)) +
  geom_text(nudge_x = 0.2, color = "black") +
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw() +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  labs(title = expression(Seasonal~effect~on~NEE~(CO[2]+CH[4])),
       subtitle = paste0("Camargue"), 
       caption = paste0("Seasonal effect is averaged across the three conservation status.\nBoxes indicate the EM mean, error bars indicate SE,\n",              "Letters indicate different groups"),
       y=expression(EMM~(g~CO[2~eq]~d^-1~m^-2)),
       x="Sampling campaign",
       col=paste0("Sampling"))

plot_co2_season_CA


##Plots Interaction effect ----
#To do, unclear interest for overview ppt (maybe interesting for some casepilots)

#Example:
#Camargue: interaction plot
plot_co2_interaction_CA<- co2_casepilot_comparisons %>%
  filter(casepilot == "CA", comparison == "interaction") %>%
  #ggplot inside brackets to be able to access casepilot name with ggtitle or define dodge for all
  { 
    #Define dodge width for separating status.
    dodge <- position_dodge(width = 0.6)
    
    ggplot(., aes(x = season, y = emmean_gwp, label = group_letter, group=status)) +
      geom_errorbar(aes(ymin = emmean_gwp - SE_gwp, ymax = emmean_gwp + SE_gwp),
                    linewidth = 0.2, width = 0.2, position = dodge) +
      geom_point(shape = 15, size = 4, aes(col=status), position = dodge) +
      geom_text(vjust = -2, color = "black", position = dodge) +
      geom_hline(yintercept=0, linetype="dashed")+
      theme_bw() +
      theme(axis.title = element_text(face = "bold"),
            axis.text = element_text(face = "bold"),
            plot.caption = element_text(hjust = 0)) +
      labs(title = expression(Seasonal~effect~on~NEE~(CO[2]+CH[4])),
           subtitle = paste0("Camargue"), 
           caption = paste0("Boxes indicate the EM mean, error bars indicate SE,\n",              "Letters indicate different groups"),
           y=expression(EMM~(g~CO[2~eq]~d^-1~m^-2)),
           x="Sampling campaign",
           col=paste0("Sampling"))
  }

plot_co2_interaction_CA





#______________________-----
#SUMMARY across GHGs-------
#CHECK-----

#R2 interpretation notes:
#As the mayority of our GLMMs are not homocedastic, we do not have a R2 (cannot calculate). We rely on the marginal (fixed effects) and conditional (fixed + random effects) R2 calculated for the homocedastic variant of each model to provide an approximated partitioning of the explanatory power of fixed vs random effects. We can use the marginal+conditional R2 for this "partitioning" for all family distributions used, HOWEVER for gamma family the R2 does not indicate variance explained of my data (R2=0.5 != "model explained 50% of total variance in data"), but instead model-derived variance on the link scale, we can only use the 2 R2s as partition of variance explained (although not exact). For gaussian family R2 exactly matches the %of variance explained interpretation, for t_family is approximate due to different impact of tails. 

#We could Add to this a Pseudo-R2 for the actual model used (in cases of heterocedasticity), to provide the overall model fit. This pseudo-R2 is calculated as 1- (deviance_model / null_deviance). Is not an aditive property (total deviance!= residual + explained)
  #A pseudoR2 of 0.47 can be interpreted as: 
# “The model reduced deviance by 47% relative to a null model, indicating substantial improvement in model fit.”
# “The pseudo-R² of 0.47 indicates that the fitted model explains 47% of the deviance (i.e., discrepancy from a saturated model) compared to an intercept-only model.”
# “The model accounts for 47% of the explainable structure in the data relative to a null model.”
#“Out of the total deviance (i.e., mismatch between model and a perfect fit), 47% was removed by including the predictors of our model.”

#To calculate the "null_model" to compare deviance against my full model, I must build it with the exact same: data, family distribution, dispformula and random effects. 


#1. Table model structure-----

#compile modeled variable (dailyflux,dailyflux_logsing, dailyflux_shift), family (gaussian,t, gamma), dispformula

modelstr_summary<- rbind(co2_model_fit_summary, ch4_model_fit_summary, co2andch4_model_fit_summary) %>%
  #get modelled variable (to infer transformation)
  mutate(response_var=gsub(" ~ status \\* season \\+ \\(1 \\| subsite\\)", "", formula)) %>% 
  dplyr::select(-c(R2_marginal, R2_conditional,is_singular,formula))%>% 
  pivot_longer(cols = c(response_var, family, dispformula),names_to = "parameter",values_to = "value") %>%
  pivot_wider(names_from = ghgspecies, values_from = value)
modelstr_summary


#2. Table  Main effects + R2 -----

#Compile Co2, Ch4 and Co2+CH4 main effects in a single table to present one variant with p_values, other with significance symbols (***,**,*,.) instead.

#Columns: casepilot, Parameter (status,season,status:season, R2m, R2c), CO2, CH4, CO2+CH4
r2_summary<- rbind(co2_model_fit_summary, ch4_model_fit_summary, co2andch4_model_fit_summary) %>% 
  dplyr::select(-c(formula, family, dispformula, is_singular)) %>% 
  pivot_longer(cols = c(R2_marginal, R2_conditional),names_to = "parameter",values_to = "value") %>% 
  mutate(value=round(value, digits = 3)) %>% 
  pivot_wider(names_from = ghgspecies, values_from = value)
r2_summary


#Compile PseudoR2 (Log-lik approach), values stored in environment:
# Get all objects in the environment that match the naming pattern
ps_r2_objects <- ls(pattern = "^pseudoR2_")

ps_r2_objects <- Filter(function(x) is.numeric(get(x)) && length(get(x)) == 1, ps_r2_objects) #remove function pseudoR2_from_null

# Initialize an empty list to store rows
r2_list <- lapply(ps_r2_objects, function(obj_name) {
  # Extract casepilot and ghg from the object name
  parts <- strsplit(obj_name, "_")[[1]]
  casepilot <- toupper(parts[2])
  ghgspecies <- paste0("gwp_",parts[3])
  
  # Get the value of the object
  r2_value <- get(obj_name)
  
  # Return as a named list
  data.frame(casepilot = casepilot, ghgspecies = ghgspecies, pseudoR2 = r2_value)
})

# Combine all rows into a single dataframe formated to merge with summary_r2
pseudo_R2_df <- do.call(rbind, r2_list) %>% 
  pivot_longer(cols = pseudoR2, names_to = "parameter", values_to = "value") %>% 
  mutate(value=round(value, digits = 3))%>% 
  pivot_wider(names_from = ghgspecies, values_from = value)

#Remove innecesary objects
rm(ps_r2_objects)
rm(list = ls(pattern = "^pseudoR2_"))


#main effects table with numeric pvalues, 
effects_summary<-  rbind(co2_results_anova, ch4_results_anova, co2andch4_results_anova) %>% 
  dplyr::select(-c(model_call, distribution_family, dispformula,intercept_pval)) %>%
  rename(status=status_pval, season=season_pval, interaction=interaction_pval) %>% 
  pivot_longer(cols=c(status,season,interaction),names_to = "parameter",values_to = "value") %>% 
  # mutate(value=round(value, digits = 3)) %>% 
  pivot_wider(names_from = ghgspecies, values_from = value) 

effects_summary


#main effects effects using symbols for significance

#Function to get symbols
pval_to_symbol <- function(p) {
  ifelse(p < 0.001, "***",
         ifelse(p < 0.01, "**",
                ifelse(p < 0.05, "*",
                       ifelse(p < 0.1, ".", "ns")
                )
         )
  )
}

#Transform pvalues to symbols
effects_summary_symbols<- effects_summary %>% 
  mutate(across(c(gwp_co2,gwp_ch4,gwp_co2andch4),pval_to_symbol))


#Combine tables R2 and Main effects
sig_and_r2_numeric_summary<-rbind(r2_summary,effects_summary,pseudo_R2_df) %>% 
  mutate(parameter=factor(parameter, levels = c("status","season", "interaction", "R2_marginal","R2_conditional","pseudoR2"))) %>% 
  arrange(casepilot,parameter)

sig_and_r2_symbol_summary <- rbind(r2_summary,pseudo_R2_df) %>% 
  mutate(across(c(gwp_co2,gwp_ch4,gwp_co2andch4),as.character)) %>% 
  rbind(effects_summary_symbols) %>% 
  mutate(parameter=factor(parameter, levels = c("status","season", "interaction", "R2_marginal","R2_conditional","pseudoR2"))) %>% 
  arrange(casepilot,parameter)

sig_and_r2_numeric_summary
sig_and_r2_symbol_summary

#Save Summary tables
#Model structures:
write.csv(modelstr_summary, file = paste0(plots_path,"/Summary_model_structures.csv"),row.names = F)

#Effects and R2 (numeric)
write.csv(sig_and_r2_numeric_summary, file = paste0(plots_path,"/Summary_effectSig_and_R2.csv"),row.names = F)

#Effects and R2 (Symbols)
write.csv(sig_and_r2_symbol_summary, file = paste0(plots_path,"/Summary_effectSigSymbols_and_R2.csv"))


# 3.Table all EMMEANS-----

emmeans_all<- merge.data.frame(co2_casepilot_comparisons, ch4_casepilot_comparisons, all = T) %>% 
  merge.data.frame(co2andch4_casepilot_comparisons,all=T)
#____________________--------
#Check effect of transformations------

# Load libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)


# 1. Create a dummy skewed dataset (centered at 0)
n <- 1000
pos_skew <- rexp(n / 2, rate = 0.1)    # long right tail
neg_skew <- -rexp(n / 2, rate = 0.1)   # long left tail

data <- c(pos_skew, neg_skew)
data <- data - mean(data)              # center at 0

# 2. Define the sign-log transformation
sign_log_transform <- function(x) {
  sign(x) * log(abs(x) + 1)
}

# 3. Define scales
scales <- c(1, 365, 10000)

# 4. Build a tidy dataframe with original and transformed values
df <- tibble(
  id = 1:length(data),
  original = data/2
) %>%
  crossing(scale = scales) %>%
  mutate(
    scaled_value = original * scale,
    transformed = sign_log_transform(scaled_value),
    scale_label = paste0("×", scale)
  )

# 5. Reshape for comparison plots
long_df <- df %>%
  dplyr::select(id, scale_label, scaled_value, transformed) %>%
  pivot_longer(cols = c(scaled_value, transformed),
               names_to = "type",
               values_to = "value") %>%
  mutate(type = dplyr::recode(type,
                       "scaled_value" = "Original",
                       "transformed" = "Transformed"))

# 6. Plot: Distribution before and after transformation
ggplot(long_df, aes(x = value, fill = type)) +
  geom_histogram(bins = 80, alpha = 0.6, position = "identity") +
  facet_wrap(type ~ scale_label, scales = "free") +
  scale_fill_manual(values = c("gray60", "steelblue")) +
  theme_minimal() +
  labs(title = "Effect of sign-log Transformation Across Scales",
       x = "Value", y = "Count", fill = "Data Type")


ggplot(long_df, aes(x = type, y= value)) +
  geom_point()+
  geom_line(aes(group=id))+
  facet_wrap(~scale_label, scales="free")
  geom_histogram(bins = 80, alpha = 0.6, position = "identity") +
  facet_wrap(type ~ scale_label, scales = "free") +
  scale_fill_manual(values = c("gray60", "steelblue")) +
  theme_minimal() +
  labs(title = "Effect of sign-log Transformation Across Scales",
       x = "Value", y = "Count", fill = "Data Type")



#____________________--------
#SEFS PLOTS-----------
  
  #CHECK PLOT TUTORIAL-----------
  # https://r-graph-gallery.com/web-violinplot-with-ggstatsplot.html
  
# save.image(paste0(plots_path,"/model_results_last.Rdata"))
  
load(paste0(plots_path,"/model_results_last.Rdata"))
      
#3 panel plots, CO2, CH4, CO2+CH4
#Boxplots + violin plot (sina dots)
#Letters from GLMM post-hoc groups
#EMMean from GLMM post-hoc groups
#Log10-scales when needed.
   
library(ggforce)
  
  log10sign<-  function(x) sign(x) * log10(abs(x) + 1)
  #CO2:  needs log10scale for visual clarity. 
  log10breaks_real<- c(-300,-100,-30,-10,-3,-1,0,1,3,10,30,100,300)
  log10breaks_trans<- log10sign(log10breaks_real)
  
  
#1. Camargue-------

  ##___co2--------
  #CO2:  needs log10scale for visual clarity. 
  
  
  ca_co2<- data4models %>% filter(casepilot=="CA", ghgspecies=="gwp_co2") %>% filter(!is.na(dailyflux))
  ca_co2_emmean<- emmeans_all %>% filter(casepilot=="CA",comparison=="status",ghgspecies=="gwp_co2")
  
ca_co2_sefsplot<-
ggplot(ca_co2, aes(x=status, y=log10sign(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(aes(col=status),position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  # geom_errorbar(data=ca_co2_emmean, aes(x=status,
  #                                       ymin=log10sign(lower.CL_gwp),
  #                                       ymax=log10sign(upper.CL_gwp)),inherit.aes = F,
  #               width=0.1)+ #ADD EMMEANS CI
  geom_point(data=ca_co2_emmean, aes(x=status,y=log10sign(emmean_gwp)),shape=23,size=3,fill = "black")+ #ADD EMMEANS
  geom_text(data=ca_co2_emmean, aes(x=status, label = group_letter, y=log10sign(100)), vjust=0.9, color = "black", inherit.aes = F, size=6, fontface="bold") +
  scale_y_continuous(minor_breaks = NULL,
    # limits = c(log10sign(-100),log10sign(100)),
    breaks = log10breaks_trans,
    labels = log10breaks_real

  ) +
  theme_bw() +
  theme(
        axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs(y= expression(CO[2]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Conservation status (n = ",nrow(ca_co2),")"),
       fill=paste0("Status"))

ca_co2_sefsplot

##___ch4 ------
# needs log10scale for visual clarity. 

ca_ch4<- ch4_formodels %>% filter(casepilot=="CA") %>% 
  filter(dailyflux>-0.025) %>% #Remove 3 outliers  
  filter(!is.na(dailyflux))
ca_ch4_emmean<- emmeans_all %>% 
  filter(casepilot=="CA",comparison=="status",ghgspecies=="gwp_ch4")

ca_ch4_sefsplot<-
ggplot(ca_ch4, aes(x=status, y=log10sign(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(aes(col=status),position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  # geom_errorbar(data=ca_ch4_emmean, aes(x=status,
  #                                       ymin=log10sign(lower.CL_gwp),
  #                                       ymax=log10sign(upper.CL_gwp)),inherit.aes = F,
  #               width=0.1)+ #ADD EMMEANS CI
  geom_point(data=ca_ch4_emmean, aes(x=status,y=log10sign(emmean_gwp)),shape=23,size=3,fill = "black")+ #ADD EMMEANS
  geom_text(data=ca_ch4_emmean, aes(x=status, label = group_letter, y=log10sign(300)), vjust=0.9, color = "black", inherit.aes = F, size=6, fontface="bold") +
  scale_y_continuous(minor_breaks = NULL,
    # limits = c(log10sign(-100),log10sign(100)),
    breaks = log10breaks_trans,
    labels = log10breaks_real
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs(y= expression(CH[4]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Conservation status (n = ",nrow(ca_ch4),")"),
       fill=paste0("Status"))

ca_ch4_sefsplot

##___co2+ch4-------
#: same trans than CO2. 

ca_co2andch4<- data4models %>% 
  filter(casepilot=="CA", ghgspecies=="gwp_co2andch4") %>% filter(!is.na(dailyflux))
ca_co2andch4_emmean<- emmeans_all %>% filter(casepilot=="CA",comparison=="status",ghgspecies=="gwp_co2andch4")

ca_co2andch4_sefsplot<- 
ggplot(ca_co2andch4, aes(x=status, y=log10sign(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(aes(col=status),position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  # geom_errorbar(data=ca_co2andch4_emmean, aes(x=status,
  #                                       ymin=log10sign(lower.CL_gwp),
  #                                       ymax=log10sign(upper.CL_gwp)),inherit.aes = F,
  #               width=0.1)+ #ADD EMMEANS CI
  geom_point(data=ca_co2andch4_emmean, aes(x=status,y=log10sign(emmean_gwp)),shape=23,size=3,fill = "black")+ #ADD EMMEANS
  geom_text(data=ca_co2andch4_emmean, aes(x=status, label = group_letter, y=log10sign(280)), vjust=0.9, color = "black", inherit.aes = F, size=6, fontface="bold") +
  scale_y_continuous(minor_breaks = NULL,
    # limits = c(log10sign(-100),log10sign(100)),
    breaks = log10breaks_trans,
    labels = log10breaks_real
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs(y= expression(CO[2]+CH[4]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Conservation status (n = ",nrow(ca_co2andch4),")"),
       fill=paste0("Status"))

ca_co2andch4_sefsplot

ca_COMBO_sefsplot <- grid.arrange(ca_co2_sefsplot,ca_ch4_sefsplot,ca_co2andch4_sefsplot, ncol=3)
ggsave(plot = ca_COMBO_sefsplot, 
       filename = "SEFSPLOTS_CA_3ghg_status_violin.png",path = plots_path,device = "png",
       dpi = 400,
       units = "mm",
       height = 100, width = 350)


#2. Curonian-------

##___co2--------
#CO2:  needs log10scale for visual clarity. 


cu_co2<- data4models %>% filter(casepilot=="CU", ghgspecies=="gwp_co2") %>% filter(!is.na(dailyflux))
cu_co2_emmean<- emmeans_all %>% filter(casepilot=="CU",comparison=="status",ghgspecies=="gwp_co2")

cu_co2_sefsplot<-
  ggplot(cu_co2, aes(x=status, y=log10sign(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(aes(col=status),position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  # geom_errorbar(data=cu_co2_emmean, aes(x=status,
  #                                       ymin=log10sign(lower.CL_gwp),
  #                                       ymax=log10sign(upper.CL_gwp)),inherit.aes = F,
  #               width=0.1)+ #ADD EMMEANS CI
  geom_point(data=cu_co2_emmean, aes(x=status,y=log10sign(emmean_gwp)),shape=23,size=3,fill = "black")+ #ADD EMMEANS
  geom_text(data=cu_co2_emmean, aes(x=status, label = group_letter, y=log10sign(100)), vjust=0.9, color = "black", inherit.aes = F, size=6, fontface="bold") +
  scale_y_continuous(minor_breaks = NULL,
    # limits = c(log10sign(-100),log10sign(100)),
    breaks = log10breaks_trans,
    labels = log10breaks_real
    
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs(y= expression(CO[2]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Conservation status (n = ",nrow(cu_co2),")"),
       fill=paste0("Status"))

cu_co2_sefsplot


##___ch4 ------
# needs log10scale for visual clarity. 

cu_ch4<- data4models %>% 
  filter(casepilot=="CU", ghgspecies=="gwp_ch4") %>% filter(!is.na(dailyflux))
cu_ch4_emmean<- emmeans_all %>% 
  filter(casepilot=="CU",comparison=="status",ghgspecies=="gwp_ch4")

cu_ch4_sefsplot<-
  ggplot(cu_ch4, aes(x=status, y=log10sign(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(aes(col=status),position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  geom_point(data=cu_ch4_emmean, aes(x=status,y=log10sign(emmean_gwp)),shape=23,size=3,fill = "black")+ #ADD EMMEANS
  geom_text(data=cu_ch4_emmean, aes(x=status, label = group_letter, y=log10sign(300)), vjust=0.9, color = "black", inherit.aes = F, size=6, fontface="bold") +
  scale_y_continuous(minor_breaks = NULL,
    # limits = c(log10sign(-100),log10sign(100)),
    breaks = log10breaks_trans,
    labels = log10breaks_real
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs(y= expression(CH[4]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Conservation status (n = ",nrow(cu_ch4),")"),
       fill=paste0("Status"))

cu_ch4_sefsplot

##___co2+ch4-------
#: same trans than CO2. 

cu_co2andch4<- data4models %>% 
  filter(casepilot=="CU", ghgspecies=="gwp_co2andch4") %>% filter(!is.na(dailyflux))
cu_co2andch4_emmean<- emmeans_all %>% filter(casepilot=="CU",comparison=="status",ghgspecies=="gwp_co2andch4")

cu_co2andch4_sefsplot<- 
  ggplot(cu_co2andch4, aes(x=status, y=log10sign(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(aes(col=status),position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  # geom_errorbar(data=cu_co2andch4_emmean, aes(x=status,
  #                                       ymin=log10sign(lower.CL_gwp),
  #                                       ymax=log10sign(upper.CL_gwp)),inherit.aes = F,
  #               width=0.1)+ #ADD EMMEANS CI
  geom_point(data=cu_co2andch4_emmean, aes(x=status,y=log10sign(emmean_gwp)),shape=23,size=3,fill = "black")+ #ADD EMMEANS
  geom_text(data=cu_co2andch4_emmean, aes(x=status, label = group_letter, y=log10sign(250)), vjust=0.9, color = "black", inherit.aes = F, size=6, fontface="bold") +
  scale_y_continuous(minor_breaks = NULL,
    # limits = c(log10sign(-100),log10sign(100)),
    breaks = log10breaks_trans,
    labels = log10breaks_real
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs(y= expression(CO[2]+CH[4]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Conservation status (n = ",nrow(cu_co2andch4),")"),
       fill=paste0("Status"))

cu_co2andch4_sefsplot

cu_COMBO_sefsplot <- grid.arrange(cu_co2_sefsplot,cu_ch4_sefsplot,cu_co2andch4_sefsplot, ncol=3)
ggsave(plot = cu_COMBO_sefsplot, 
       filename = "SEFSPLOTS_CU_3ghg_status_violin.png",path = plots_path,device = "png",
       dpi = 400,
       units = "mm",
       height = 100, width = 350)



#3. Danube-------

##___co2--------
#CO2:  needs log10scale for visual clarity. 


da_co2<- data4models %>% filter(casepilot=="DA", ghgspecies=="gwp_co2") %>% filter(!is.na(dailyflux))
da_co2_emmean<- emmeans_all %>% filter(casepilot=="DA",comparison=="status",ghgspecies=="gwp_co2")

da_co2_sefsplot<-
  ggplot(da_co2, aes(x=status, y=log10sign(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(aes(col=status),position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  # geom_errorbar(data=da_co2_emmean, aes(x=status,
  #                                       ymin=log10sign(lower.CL_gwp),
  #                                       ymax=log10sign(upper.CL_gwp)),inherit.aes = F,
  #               width=0.1)+ #ADD EMMEANS CI
  geom_point(data=da_co2_emmean, aes(x=status,y=log10sign(emmean_gwp)),shape=23,size=3,fill = "black")+ #ADD EMMEANS
  geom_text(data=da_co2_emmean, aes(x=status, label = group_letter, y=log10sign(200)), vjust=0.9, color = "black", inherit.aes = F, size=6, fontface="bold") +
  scale_y_continuous(minor_breaks = NULL,
    # limits = c(log10sign(-100),log10sign(100)),
    breaks = log10breaks_trans,
    labels = log10breaks_real
    
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs(y= expression(CO[2]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Conservation status (n = ",nrow(da_co2),")"),
       fill=paste0("Status"))

da_co2_sefsplot


##___ch4 ------
# needs log10scale for visual clarity. 

da_ch4<- data4models %>% 
  filter(casepilot=="DA", ghgspecies=="gwp_ch4") %>% filter(!is.na(dailyflux))
da_ch4_emmean<- emmeans_all %>% 
  filter(casepilot=="DA",comparison=="status",ghgspecies=="gwp_ch4")

da_ch4_sefsplot<-
  ggplot(da_ch4, aes(x=status, y=log10sign(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(aes(col=status),position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  # geom_errorbar(data=da_ch4_emmean, aes(x=status,
  #                                       ymin=log10sign(lower.CL_gwp),
  #                                       ymax=log10sign(upper.CL_gwp)),inherit.aes = F,
  #               width=0.1)+ #ADD EMMEANS CI
  geom_point(data=da_ch4_emmean, aes(x=status,y=log10sign(emmean_gwp)),shape=23,size=3,fill = "black")+ #ADD EMMEANS
  geom_text(data=da_ch4_emmean, aes(x=status, label = group_letter, y=log10sign(700)), vjust=0.9, color = "black", inherit.aes = F, size=6, fontface="bold") +
  scale_y_continuous(minor_breaks = NULL,
    # limits = c(log10sign(-100),log10sign(100)),
    breaks = log10breaks_trans,
    labels = log10breaks_real
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs(y= expression(CH[4]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Conservation status (n = ",nrow(da_ch4),")"),
       fill=paste0("Status"))

da_ch4_sefsplot

##___co2+ch4-------
#: same trans than CO2. 

da_co2andch4<- data4models %>% 
  filter(casepilot=="DA", ghgspecies=="gwp_co2andch4") %>% filter(!is.na(dailyflux))
da_co2andch4_emmean<- emmeans_all %>% filter(casepilot=="DA",comparison=="status",ghgspecies=="gwp_co2andch4")

da_co2andch4_sefsplot<- 
  ggplot(da_co2andch4, aes(x=status, y=log10sign(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(aes(col=status),position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  # geom_errorbar(data=da_co2andch4_emmean, aes(x=status,
  #                                       ymin=log10sign(lower.CL_gwp),
  #                                       ymax=log10sign(upper.CL_gwp)),inherit.aes = F,
  #               width=0.1)+ #ADD EMMEANS CI
  geom_point(data=da_co2andch4_emmean, aes(x=status,y=log10sign(emmean_gwp)),shape=23,size=3,fill = "black")+ #ADD EMMEANS
  geom_text(data=da_co2andch4_emmean, aes(x=status, label = group_letter, y=log10sign(600)), vjust=0.9, color = "black", inherit.aes = F, size=6, fontface="bold") +
  scale_y_continuous(minor_breaks = NULL,
    # limits = c(log10sign(-100),log10sign(100)),
    breaks = log10breaks_trans,
    labels = log10breaks_real
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs(y= expression(CO[2]+CH[4]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Conservation status (n = ",nrow(da_co2andch4),")"),
       fill=paste0("Status"))

da_co2andch4_sefsplot

da_COMBO_sefsplot <- grid.arrange(da_co2_sefsplot,da_ch4_sefsplot,da_co2andch4_sefsplot, ncol=3)
ggsave(plot = da_COMBO_sefsplot, 
       filename = "SEFSPLOTS_DA_3ghg_status_violin.png",path = plots_path,device = "png",
       dpi = 400,
       units = "mm",
       height = 100, width = 350)



#4. Dutch delta-------
#WE might try non-transformed CH4 removing a couple of outliers. 

##___co2--------
#CO2:  needs log10scale for visual clarity. 


du_co2<- data4models %>% 
  filter(casepilot=="DU", ghgspecies=="gwp_co2") %>% filter(!is.na(dailyflux))
du_co2_emmean<- emmeans_all %>% 
  filter(casepilot=="DU",comparison=="status",ghgspecies=="gwp_co2")

du_co2_sefsplot<-
  ggplot(du_co2, aes(x=status, y=log10sign(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(aes(col=status),position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5,,scale = "width")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  geom_point(data=du_co2_emmean, aes(x=status,y=log10sign(emmean_gwp)),shape=23,size=3,fill = "black")+ #ADD EMMEANS
  geom_text(data=du_co2_emmean, aes(x=status, label = group_letter, y=log10sign(100)), vjust=0.9, color = "black", inherit.aes = F,fontface = "bold",size=6 ) +
  scale_y_continuous(minor_breaks = NULL,
    # limits = c(log10sign(-100),log10sign(100)),
    breaks = log10breaks_trans,
    labels = log10breaks_real
    
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs(y= expression(CO[2]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Conservation status (n = ",nrow(du_co2),")"),
       fill=paste0("Status"))

du_co2_sefsplot


##___ch4 ------
# 

du_ch4<- data4models %>% 
  filter(casepilot=="DU", ghgspecies=="gwp_ch4") %>% filter(!is.na(dailyflux)) %>% 
  filter(dailyflux<0.1)# remove  outliers

du_ch4_emmean<- emmeans_all %>% 
  filter(casepilot=="DU",comparison=="status",ghgspecies=="gwp_ch4")

du_ch4_sefsplot<-
  ggplot(du_ch4, aes(x=status, y=log10sign(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(aes(col=status),position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5,,scale = "width")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  geom_point(data=du_ch4_emmean, aes(x=status,y=log10sign(emmean_gwp)),shape=23,size=3,fill = "black")+ #ADD EMMEANS
  geom_text(data=du_ch4_emmean, aes(x=status, label = group_letter, y=log10sign(0.2)), vjust=0.9, color = "black", inherit.aes = F,fontface = "bold",size=6 ) +
  scale_y_continuous(minor_breaks = NULL,
    limits = c(log10sign(-0.1),log10sign(0.3)),
    breaks = log10sign(c(-0.1,-0.03,-0.01,0,0.01,0.03,0.1,0.3)),
    labels =c(-0.1,-0.03,-0.01,0,0.01,0.03,0.1,0.3)
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs(y= expression(CH[4]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Conservation status (n = ",nrow(du_ch4),")"),
       fill=paste0("Status"))

du_ch4_sefsplot

##___co2+ch4-------
#: same trans than CO2. 

du_co2andch4<- data4models %>% 
  filter(casepilot=="DU", ghgspecies=="gwp_co2andch4") %>% filter(!is.na(dailyflux))
du_co2andch4_emmean<- emmeans_all %>%
  filter(casepilot=="DU",comparison=="status",ghgspecies=="gwp_co2andch4")

du_co2andch4_sefsplot<- 
  ggplot(du_co2andch4, aes(x=status, y=log10sign(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(aes(col=status),position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5,scale = "width")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  geom_point(data=du_co2andch4_emmean, aes(x=status,y=log10sign(emmean_gwp)),shape=23,size=3,fill = "black")+ #ADD EMMEANS
  geom_text(data=du_co2andch4_emmean, aes(x=status, label = group_letter, y=log10sign(90)), vjust=0.9, color = "black", inherit.aes = F, size=6, fontface="bold") +
  scale_y_continuous(minor_breaks = NULL,
    # limits = c(log10sign(-100),log10sign(100)),
    breaks = log10breaks_trans,
    labels = log10breaks_real
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs(y= expression(CO[2]+CH[4]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Conservation status (n = ",nrow(du_co2andch4),")"),
       fill=paste0("Status"))

du_co2andch4_sefsplot

du_COMBO_sefsplot <- grid.arrange(du_co2_sefsplot,du_ch4_sefsplot,du_co2andch4_sefsplot, ncol=3)
ggsave(plot = du_COMBO_sefsplot, 
       filename = "SEFSPLOTS_DU_3ghg_status_violin.png",path = plots_path,device = "png",
       dpi = 400,
       units = "mm",
       height = 100, width = 350)



#5. Ria de Aveiro-------

##___co2 (log)--------

ri_co2<- data4models %>% 
  filter(casepilot=="RI", ghgspecies=="gwp_co2") %>% filter(!is.na(dailyflux))%>% 
  filter(between(dailyflux, -7.5,2.5 )) #Remove 3 outliers

ri_co2_emmean<- emmeans_all %>% 
  filter(casepilot=="RI",comparison=="status",ghgspecies=="gwp_co2")

ri_co2_sefsplot<-
  ggplot(ri_co2, aes(x=status, y=log10sign(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(aes(col=status),position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  geom_point(data=ri_co2_emmean, aes(x=status,y=log10sign(emmean_gwp)),shape=23,size=3,fill = "black")+ #ADD EMMEANS
  geom_text(data=ri_co2_emmean, aes(x=status, label = group_letter, y=log10sign(3)), vjust=0.9, color = "black", inherit.aes = F, fontface = "bold",size = 6) +
  scale_y_continuous(minor_breaks = NULL,
    limits = c(log10sign(-10),log10sign(10)),
    breaks = log10breaks_trans,
    labels = log10breaks_real
    
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs(y= expression(CO[2]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Conservation status (n = ",nrow(ri_co2),")"),
       fill=paste0("Status"))

ri_co2_sefsplot


##___ch4 (log) ------


ri_ch4<- data4models %>% 
  filter(casepilot=="RI", ghgspecies=="gwp_ch4") %>% filter(!is.na(dailyflux))%>% 
  filter(between(dailyflux, -0.025,0.025)) #CH4: limits = c(-0.025,0.025), (removes 8 outliers)

ri_ch4_emmean<- emmeans_all %>% 
  filter(casepilot=="RI",comparison=="status",ghgspecies=="gwp_ch4")

ri_ch4_sefsplot<-
  ggplot(ri_ch4, aes(x=status, y=log10sign(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(aes(col=status),position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  # geom_errorbar(data=ri_ch4_emmean, aes(x=status,
  #                                       ymin=log10sign(lower.CL_gwp),
  #                                       ymax=log10sign(upper.CL_gwp)),inherit.aes = F,
  #               width=0.1)+ #ADD EMMEANS CI
  geom_point(data=ri_ch4_emmean, aes(x=status,y=log10sign(emmean_gwp)),shape=23,size=3,fill = "black")+ #ADD EMMEANS
  geom_text(data=ri_ch4_emmean, aes(x=status, label = group_letter, y=log10sign(0.03)), vjust=0.9, color = "black", inherit.aes = F, fontface = "bold",size = 6) +
  scale_y_continuous(minor_breaks = NULL,
    limits = c(log10sign(-0.03),log10sign(0.03)),
    breaks = log10sign(c(-0.03,-0.01,0,0.01,0.03))
,
    labels = c(-0.03,-0.01,0,0.01,0.03)
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs( 
       y= expression(CH[4]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Conservation status (n = ",nrow(ri_ch4),")"),
       fill=paste0("Status"))
ri_ch4_sefsplot

##___co2+ch4 (log)-------
#: same trans than CO2. 

ri_co2andch4<- data4models %>% 
  filter(casepilot=="RI", ghgspecies=="gwp_co2andch4") %>% filter(!is.na(dailyflux))%>% 
  filter(between(dailyflux, -7.5,2.5 )) #Remove 3 outliers

ri_co2andch4_emmean<- emmeans_all %>%
  filter(casepilot=="RI",comparison=="status",ghgspecies=="gwp_co2andch4")

ri_co2andch4_sefsplot<- 
  ggplot(ri_co2andch4, aes(x=status, y=log10sign(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(aes(col=status),position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  geom_point(data=ri_co2andch4_emmean, aes(x=status,y=log10sign(emmean_gwp)),shape=23,size=3,fill = "black")+ #ADD EMMEANS
  geom_text(data=ri_co2andch4_emmean, aes(x=status, label = group_letter, y=log10sign(3)), vjust=0.9, color = "black", inherit.aes = F, size = 6,fontface = "bold") +
  scale_y_continuous(minor_breaks = NULL,
    limits = c(log10sign(-10),log10sign(10)),
    breaks = log10breaks_trans,
    labels = log10breaks_real
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs(
       y= expression(CO[2]+CH[4]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Conservation status (n = ",nrow(ri_ch4),")"),
       fill=paste0("Status"))

ri_co2andch4_sefsplot

ri_COMBO_sefsplot <- grid.arrange(ri_co2_sefsplot,ri_ch4_sefsplot,ri_co2andch4_sefsplot, ncol=3)
ggsave(plot = ri_COMBO_sefsplot, 
       filename = "SEFSPLOTS_RI_3ghg_status_violin.png",path = plots_path,device = "png",
       dpi = 400,
       units = "mm",
       height = 100, width = 350)

##___co2 (linear)--------

ri_co2<- data4models %>% 
  filter(casepilot=="RI", ghgspecies=="gwp_co2") %>% filter(!is.na(dailyflux))%>% 
  filter(between(dailyflux, -7.5,2.5 )) #Remove 3 outliers

ri_co2_emmean<- emmeans_all %>% 
  filter(casepilot=="RI",comparison=="status",ghgspecies=="gwp_co2")

ri_co2_sefsplot_linear<-
  ggplot(ri_co2, aes(x=status, y=(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(aes(col=status),position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  geom_point(data=ri_co2_emmean, aes(x=status,y=(emmean_gwp)),shape=23,size=3,fill = "black")+ #ADD EMMEANS
  geom_text(data=ri_co2_emmean, aes(x=status, label = group_letter, y=(3)), vjust=0.9, color = "black", inherit.aes = F, fontface = "bold",size = 6) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs(y= expression(CO[2]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Conservation status (n = ",nrow(ri_co2),")"),
       fill=paste0("Status"))

ri_co2_sefsplot_linear


##___ch4 (linear)------


ri_ch4<- data4models %>% 
  filter(casepilot=="RI", ghgspecies=="gwp_ch4") %>% filter(!is.na(dailyflux))%>% 
  filter(between(dailyflux, -0.025,0.025)) #CH4: limits = c(-0.025,0.025), (removes 8 outliers)

ri_ch4_emmean<- emmeans_all %>% 
  filter(casepilot=="RI",comparison=="status",ghgspecies=="gwp_ch4")

ri_ch4_sefsplot_linear<-
  ggplot(ri_ch4, aes(x=status, y=(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(aes(col=status),position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  geom_point(data=ri_ch4_emmean, aes(x=status,y=(emmean_gwp)),shape=23,size=3,fill = "black")+ #ADD EMMEANS
  geom_text(data=ri_ch4_emmean, aes(x=status, label = group_letter, y=(0.03)), vjust=0.9, color = "black", inherit.aes = F, fontface = "bold",size = 6) +

  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs( 
    y= expression(CH[4]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
    x=paste0("Conservation status (n = ",nrow(ri_ch4),")"),
    fill=paste0("Status"))
ri_ch4_sefsplot_linear

##___co2+ch4 (linear)-------
#: same trans than CO2. 

ri_co2andch4<- data4models %>% 
  filter(casepilot=="RI", ghgspecies=="gwp_co2andch4") %>% filter(!is.na(dailyflux))%>% 
  filter(between(dailyflux, -7.5,2.5 )) #Remove 3 outliers

ri_co2andch4_emmean<- emmeans_all %>%
  filter(casepilot=="RI",comparison=="status",ghgspecies=="gwp_co2andch4")

ri_co2andch4_sefsplot_linear<- 
  ggplot(ri_co2andch4, aes(x=status, y=(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(aes(col=status),position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  geom_point(data=ri_co2andch4_emmean, aes(x=status,y=(emmean_gwp)),shape=23,size=3,fill = "black")+ #ADD EMMEANS
  geom_text(data=ri_co2andch4_emmean, aes(x=status, label = group_letter, y=(3)), vjust=0.9, color = "black", inherit.aes = F, size = 6,fontface = "bold") +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs(
    y= expression(CO[2]+CH[4]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
    x=paste0("Conservation status (n = ",nrow(ri_ch4),")"),
    fill=paste0("Status"))

ri_co2andch4_sefsplot_linear

ri_COMBO_sefsplot_linear <- grid.arrange(ri_co2_sefsplot_linear,ri_ch4_sefsplot_linear,ri_co2andch4_sefsplot_linear, ncol=3)
ggsave(plot = ri_COMBO_sefsplot_linear, 
       filename = "SEFSPLOTS_RI_linear_3ghg_status_violin.png",path = plots_path,device = "png",
       dpi = 400,
       units = "mm",
       height = 100, width = 350)




#6. Valencia-------

##___co2--------
#CO2:  needs log10scale for visual clarity. 


va_co2<- data4models %>% filter(casepilot=="VA", ghgspecies=="gwp_co2") %>% filter(!is.na(dailyflux))
va_co2_emmean<- emmeans_all %>% filter(casepilot=="VA",comparison=="status",ghgspecies=="gwp_co2")

va_co2_sefsplot<-
  ggplot(va_co2, aes(x=status, y=log10sign(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(aes(col=status),position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  # geom_errorbar(data=va_co2_emmean, aes(x=status,
  #                                       ymin=log10sign(lower.CL_gwp),
  #                                       ymax=log10sign(upper.CL_gwp)),inherit.aes = F,
  #               width=0.1)+ #ADD EMMEANS CI
  geom_point(data=va_co2_emmean, aes(x=status,y=log10sign(emmean_gwp)),shape=23,size=3,fill = "black")+ #ADD EMMEANS
  geom_text(data=va_co2_emmean, aes(x=status, label = group_letter, y=log10sign(120)), vjust=0.9, color = "black", inherit.aes = F, nudge_x = 0.1, size = 6,fontface = "bold") +
  scale_y_continuous(minor_breaks = NULL,
    # limits = c(log10sign(-100),log10sign(100)),
    breaks = log10breaks_trans,
    labels = log10breaks_real
    
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs(y= expression(CO[2]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Conservation status (n = ",nrow(va_co2),")"),
       fill=paste0("Status"))

va_co2_sefsplot


##___ch4 ------
# needs log10scale for visual clarity. 

va_ch4<- data4models %>% 
  filter(casepilot=="VA", ghgspecies=="gwp_ch4") %>% filter(!is.na(dailyflux))
va_ch4_emmean<- emmeans_all %>% 
  filter(casepilot=="VA",comparison=="status",ghgspecies=="gwp_ch4")

va_ch4_sefsplot<-
  ggplot(va_ch4, aes(x=status, y=log10sign(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(aes(col=status),position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  # geom_errorbar(data=va_ch4_emmean, aes(x=status,
  #                                       ymin=log10sign(lower.CL_gwp),
  #                                       ymax=log10sign(upper.CL_gwp)),inherit.aes = F,
  #               width=0.1)+ #ADD EMMEANS CI
  geom_point(data=va_ch4_emmean, aes(x=status,y=log10sign(emmean_gwp)),shape=23,size=3,fill = "black")+ #ADD EMMEANS
  geom_text(data=va_ch4_emmean, aes(x=status, label = group_letter, y=log10sign(30)), vjust=0.9, color = "black", inherit.aes = F, nudge_x = 0.1, size = 6,fontface = "bold") +
  scale_y_continuous(minor_breaks = NULL,
    # limits = c(log10sign(-100),log10sign(100)),
    breaks = log10breaks_trans,
    labels = log10breaks_real
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs(y= expression(CH[4]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Conservation status (n = ",nrow(va_ch4),")"),
       fill=paste0("Status"))

va_ch4_sefsplot

##___co2+ch4-------
#: same trans than CO2. 

va_co2andch4<- data4models %>% 
  filter(casepilot=="VA", ghgspecies=="gwp_co2andch4") %>% filter(!is.na(dailyflux))
va_co2andch4_emmean<- emmeans_all %>% filter(casepilot=="VA",comparison=="status",ghgspecies=="gwp_co2andch4")

va_co2andch4_sefsplot<- 
  ggplot(va_co2andch4, aes(x=status, y=log10sign(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(aes(col=status),position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  geom_point(data=va_co2andch4_emmean, aes(x=status,y=log10sign(emmean_gwp)),shape=23,size=3,fill = "black")+ #ADD EMMEANS
  geom_text(data=va_co2andch4_emmean, aes(x=status, label = group_letter, y=log10sign(120)), vjust=0.9, color = "black", inherit.aes = F, nudge_x = 0.1, size = 6,fontface = "bold") +
  scale_y_continuous(minor_breaks = NULL,
    # limits = c(log10sign(-100),log10sign(100)),
    breaks = log10breaks_trans,
    labels = log10breaks_real
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs(y= expression(CO[2]+CH[4]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Conservation status (n = ",nrow(va_co2andch4),")"),
       fill=paste0("Status"))

va_co2andch4_sefsplot

va_COMBO_sefsplot <- grid.arrange(va_co2_sefsplot,va_ch4_sefsplot,va_co2andch4_sefsplot, ncol=3)
ggsave(plot = va_COMBO_sefsplot, 
       filename = "SEFSPLOTS_VA_3ghg_status_violin.png",path = plots_path,device = "png",
       dpi = 400,
       units = "mm",
       height = 100, width = 350)



#Preserved across CP--------


#___Boxplot + violin -------
data4models %>% 
  filter(status=="Preserved", ghgspecies%in%c("gwp_co2")) %>% 
  mutate(casepilot=factor(casepilot, levels = c("DU","RI","CA","VA","DA","CU"))) %>% 
  ggplot(aes(x=casepilot, y=log10sign(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(col ="#00BA38",position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, fill="#00BA38",size = 0.7)+
  scale_y_continuous(minor_breaks = NULL,
                     # limits = c(log10sign(-100),log10sign(100)),
                     breaks = log10breaks_trans,
                     labels = log10breaks_real
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
  guides(color = "none", fill="none")+
  labs(title=expression(CO[2]~GWP) ,y= expression(CO[2]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Casepilot"),
       fill=paste0("Status"))


ggsave( 
       filename = "Preserved only_CO2_gwp.png",path = plots_path,device = "png",
       dpi = 400,
       units = "mm",
       height = 100, width = 125)




data4models %>% 
  filter(status=="Preserved", ghgspecies%in%c("gwp_ch4")) %>% 
  mutate(casepilot=factor(casepilot, levels = c("DU","RI","CA","VA","DA","CU"))) %>% 
  ggplot(aes(x=casepilot, y=log10sign(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(col ="#00BA38",position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, fill="#00BA38",size = 0.7)+
  scale_y_continuous(minor_breaks = NULL,
                     # limits = c(log10sign(-100),log10sign(100)),
                     breaks = log10breaks_trans,
                     labels = log10breaks_real
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
  guides(color = "none", fill="none")+
  labs(title=expression(CH[4]~GWP) ,y= expression(CH[4]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Casepilot"),
       fill=paste0("Status"))


ggsave( 
  filename = "Preserved only_CH4_gwp.png",path = plots_path,device = "png",
  dpi = 400,
  units = "mm",
  height = 100, width = 125)


#___Boxplot only-------

data4models %>% 
  filter(status=="Preserved", ghgspecies%in%c("gwp_co2")) %>% 
  mutate(casepilot=factor(casepilot, levels = c("DU","RI","CA","VA","DA","CU"))) %>% 
  ggplot(aes(x=casepilot, y=log10sign(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  # geom_violin(col ="#00BA38",position = position_dodge(0.9),linewidth = 1,scale = "width")+
  # geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, fill="#00BA38",size = 0.7)+
  scale_y_continuous(minor_breaks = NULL,
                     # limits = c(log10sign(-100),log10sign(100)),
                     breaks = log10breaks_trans,
                     labels = log10breaks_real
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
  guides(color = "none", fill="none")+
  labs(title=expression(CO[2]~GWP) ,y= expression(CO[2]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Casepilot"),
       fill=paste0("Status"))


ggsave( 
  filename = "Preserved only_just_boxplot_CO2_gwp.png",path = plots_path,device = "png",
  dpi = 400,
  units = "mm",
  height = 100, width = 125)




data4models %>% 
  filter(status=="Preserved", ghgspecies%in%c("gwp_ch4")) %>% 
  mutate(casepilot=factor(casepilot, levels = c("DU","RI","CA","VA","DA","CU"))) %>% 
  ggplot(aes(x=casepilot, y=log10sign(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  # geom_violin(col ="#00BA38",position = position_dodge(0.9),linewidth = 1,scale = "width")+
  # geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, fill="#00BA38",size = 0.7)+
  scale_y_continuous(minor_breaks = NULL,
                     # limits = c(log10sign(-100),log10sign(100)),
                     breaks = log10breaks_trans,
                     labels = log10breaks_real
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
  guides(color = "none", fill="none")+
  labs(title=expression(CH[4]~GWP) ,y= expression(CH[4]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Casepilot"),
       fill=paste0("Status"))


ggsave( 
  filename = "Preserved only_just boxplot_CH4_gwp.png",path = plots_path,device = "png",
  dpi = 400,
  units = "mm",
  height = 100, width = 125)




#____________________--------


#(skip) GLOBAL models -----
#Below are the initial tests to produce global models for each GHG, as they become very complex (3-way interactions, different data tranformations,...) We will not focus on global models but on the CASEPILOT models (from previous sections of script)

#_________--------

#CO2------


#one per ghgspecies, using data from all casepilots. Excluding sample_weight (for the moment)

## ========== 0. subset & trans ==========
data_co2_all <- data4models %>% 
  filter(ghgspecies == "gwp_co2", !is.na(dailyflux)) %>% 
  mutate(
    #logsign transformation
    dailyflux_logsign =logsign(dailyflux),
    #logroot transformation 
    dailyflux_rootsign = rootsign(dailyflux)
  )

#Outlier exploration: 
#Identify outliers (with data already transformed)
co2log_outliers<-data_co2_all %>% 
  group_by(status,season,casepilot,subsite,strata) %>% 
  mutate(is.extreme=is_extreme(dailyflux_logsign)) %>% 
  filter(is.extreme==T) %>% 
  mutate(subsite2=substr(subsite,4,5))

#visualize outliers
data_co2_all %>% 
  mutate(subsite2=substr(subsite,4,5)) %>% 
  ggplot(aes(x=subsite2, y=dailyflux_logsign, fill=status, group=paste0(subsite)))+
  geom_boxplot()+
  geom_point(data = co2log_outliers, col="red")+
  facet_wrap(~casepilot,scales="free")

#Extreme outliers (3*IQR, by status,season,casepilot, subsite, strata) could be removed, if removed, sample_weights should be adjusted accordingly to retain strata-representativity (re-normalize). 
#SHOULD WE REMOVE OUTLIERS? 
#DECISSION: NO. Our "outliers" represent true variability in field conditions. Even when grouped by strata,  data-driven outlier removal should NOT be done (as we cannot account for small scale variability in our dataset) 


## ========== 1. Untransformed Model ==========
model_untrans <- lmer(dailyflux ~ status * season * casepilot + (1 | subsite), 
                      data = data_co2_all)

## ========== 2. Signed Log Transform Model ==========
model_logsign <- lmer(dailyflux_logsign ~ status * season * casepilot + (1 | subsite),
                      data = data_co2_all)

## ========== 3. Signed Root Transform Model ==========
model_rootsign <- lmer(dailyflux_rootsign ~ status * season * casepilot + (1 | subsite),
                       data = data_co2_all)

## ========== 4. Evaluate Model requirements ==========
#This function produces diagnostics for the LMM models, either the general ones (including casepilot as factor) or the per-casepilot LMMs
check_model_diagnostics <- function(model, model_name, data) {
  cat("----- Diagnostics for", model_name, "-----\n")
  
  print(ggResidpanel::resid_panel(model, qqbands = TRUE, plots = c("qq", "resid", "hist", "fitted")))
  
  res <- residuals(model)
  
  # Residuals vs status
  if ("status" %in% names(data)) {
    p1 <- ggplot(data.frame(status = data$status, resid = res), aes(x = status, y = resid)) +
      geom_boxplot() +
      labs(title = paste("Residuals by Status:", model_name), y = "Residuals")
    print(p1)
  }
  
  # Residuals vs casepilot
  if ("casepilot" %in% names(data) && nlevels(factor(data$casepilot)) > 1) {
    p2 <- ggplot(data.frame(casepilot = data$casepilot, resid = res), aes(x = casepilot, y = resid)) +
      geom_boxplot() +
      labs(title = paste("Residuals by Casepilot:", model_name), y = "Residuals")
    print(p2)
  }
  
  # Residuals vs season
  if ("season" %in% names(data)) {
    p3 <- ggplot(data.frame(season = data$season, resid = res), aes(x = season, y = resid)) +
      geom_boxplot() +
      labs(title = paste("Residuals by Season:", model_name), y = "Residuals")
    print(p3)
  }
  
  # Levene's Tests
  if ("status" %in% names(data) && nlevels(factor(data$status)) > 1) {
    cat("Levene's Test for residual variance by status:\n")
    print(leveneTest(res ~ data$status))
  }
  
  if ("casepilot" %in% names(data) && nlevels(factor(data$casepilot)) > 1) {
    cat("Levene's Test for residual variance by casepilot:\n")
    print(leveneTest(res ~ data$casepilot))
  }
  
  if ("season" %in% names(data) && nlevels(factor(data$season)) > 1) {
    cat("Levene's Test for residual variance by season:\n")
    print(leveneTest(res ~ data$season))
  }
  
  # Shapiro-Wilk test
  cat("Shapiro-Wilk test for residual normality:\n")
  print(shapiro.test(res))
  
  # R², AIC, BIC
  r2 <- performance::r2_nakagawa(model)
  aic_val <- AIC(model)
  bic_val <- BIC(model)
  
  # 10-fold CV RMSE
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
diag_untrans <- check_model_diagnostics(model_untrans, "Untransformed", data_co2_all)
#Untransformed Observations: clear non-normality of residuals, heterocedasticity at season and casepilot level (not at status level)
diag_logsign <- check_model_diagnostics(model_logsign, "Signed Log", data_co2_all)
#Logtransformed Observations: Residuals show reasonable normality (visual). Heterocedasticity at season and casepilot level (not at status level), OK heterocedasticity is expected. BEST
diag_rootsign <- check_model_diagnostics(model_rootsign, "Signed Root", data_co2_all)
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

#Set back-transformation function accordingly: inv_logsign or inv_rootsign
inv_function<- inv_logsign

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


## ---- 5.3. EMMs  ----
#Estimated Marginal Means (emms)


#EMMs of status, for each combination of season and case-pilot: Within each season and casepilot, how do the statuses compare?
emm_status<-emmeans(best_model, ~ status | season * casepilot)
summary(emm_status)

#The above call of emmeans is the option to be used for post-hoc comparisons of status effects within ecological and seasonal contexts, results are equivalent to emm_all <- emmeans(best_model,~ status * season * casepilot) but grouped differently, so that calling pairs(emm_status) with default options will already provide the meaningful comparisons. If I was to call pairs(emm_all) it would return every pairwise comparison across all combinations fo the 3 factors, which is meaningless. 

#Back-transform EMMs
emm_status_bt<- as.data.frame(emm_status) %>% 
  mutate(across(c(emmean, SE,lower.CL,upper.CL), inv_function))


#Plot EMMs 
ggplot(emm_status_bt, aes(x = season, y = emmean, fill = status)) +
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

#Compare EEMs vs weigthed averages (calcualte SE of weighted average)
weighted_avgs_og <- data_co2_all %>%
  group_by(status, casepilot, season) %>%
  summarise(
    weighted_mean = sum(dailyflux * sample_weight) / sum(sample_weight),
    sum_w = sum(sample_weight),
    sum_w2 = sum(sample_weight^2),
    # Temporarily compute mean for variance calculation
    mean_w = sum(dailyflux * sample_weight) / sum(sample_weight),
    # Weighted variance
    var_w = sum(sample_weight * (dailyflux - mean_w)^2) / sum(sample_weight),
    # Effective sample size
    n_eff = (sum_w)^2 / sum_w2,
    # Standard error
    se = sqrt(var_w) / sqrt(n_eff),
    .groups = "drop"
  ) %>%
  # Keep only relevant output columns
  dplyr::select(status, casepilot, season, weighted_mean, se)

#Plot weighted averages (original scale) 
ggplot(weighted_avgs_og, aes(x = season, y = weighted_mean, fill = status)) +
  geom_col(position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = weighted_mean - se, ymax = weighted_mean + se),
                position = position_dodge(width = 0.8), width = 0.2) +
  facet_wrap(~ casepilot,scales="free") +
  labs(
    title = "Weighted avg by Status, Season, and Casepilot",
    y = expression(Average~CO[2]~dailyflux~(g~CO2[eq]~d^-1~m^-2)),
    x = "Season",
    col = "Status"
  ) +
  theme_minimal()

#PLot back transformed EMMs vs Weighted averages (original scale)
emm_status_bt %>% 
  left_join(weighted_avgs_og) %>% 
  ggplot(aes(x=weighted_mean,y=emmean))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  ggtitle("BT EMMs vs Weighted avg of untransformed\n(same groups)")
#There is a large deviation because the EMMs are calculated based on transformed values (as they should), and the back transformed average is not the same as the average of back transformed (original scale)

#JUST FOR CONSISTENCY: if we calculate the weighted averages on the transformed values and then back transform the averages, the results are equal to those from the model emmeans .



## ---- 5.4. Post-hoc on EMMs  ----
#All statistics, including post-hoc must be performed in the same scale as the model (i.e. sign-log transformed), only  back-transform for plots, interpretation and reporting effect sizes. 

# EMM comparison for status by season and casepilot comparisons
pairs(emm_status, adjust="tukey")

#Only significant differences for a few cases of DANUBE DELTA: 
pairs(emm_status, adjust="tukey") %>% as.data.frame() %>% 
  filter(p.value<0.05)


## ---- 5.5. Model complexity  ----

#Model too complex? 

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
                     data=data_co2_all)
anova(model_no3way, best_model)
#Current more complex model significantly improves the fit, how much added explanatory power? is it worth it?


#BASED ON Type III ANOVA Results: 
# ✅ Interpretation of the Type III ANOVA Results
# The main effect of status is not significant, but its interactions with season and casepilot are.
# The main effect of casepilot is highly significant, which aligns with your expectation that these are ecologically distinct systems.
# 🔍 What this suggests:
#   The effect of status is not consistent across all casepilots or seasons — it depends on context.
# The strong casepilot effect may indeed obscure or dilute the status effect when analyzed in a pooled model, especially if the direction or magnitude of the status effect varies across casepilots.


#THEREFORE: Modeling Each Casepilot Separately is the better option.

# Interactions are significant and complex (as in your case).
# We expect heterogeneity in how status affects emissions across ecosystems.
# We want to avoid overfitting or interpretational complexity from 3-way interactions.
# By modeling each casepilot individually:

# We simplify the model (no 3-way interactions).
# We gain clarity on how status behaves in each ecological context.
# We can later compare or summarize these effects across casepilots (e.g., via forest plots or meta-analytic summaries).

#_________----
#CH4---------


#SHIFT would have to be selected for the whole dataset (i.e. apply biggest shift to all casepilots)
#TO BE ADAPTED: 

## ========== 0. subset & trans ==========
data_co2_all <- data4models %>% 
  filter(ghgspecies == "gwp_co2", !is.na(dailyflux)) %>% 
  mutate(
    #logsign transformation
    dailyflux_logsign =logsign(dailyflux),
    #logroot transformation 
    dailyflux_rootsign = rootsign(dailyflux)
  )

#Outlier exploration: 
#Identify outliers (with data already transformed)
co2log_outliers<-data_co2_all %>% 
  group_by(status,season,casepilot,subsite,strata) %>% 
  mutate(is.extreme=is_extreme(dailyflux_logsign)) %>% 
  filter(is.extreme==T) %>% 
  mutate(subsite2=substr(subsite,4,5))

#visualize outliers
data_co2_all %>% 
  mutate(subsite2=substr(subsite,4,5)) %>% 
  ggplot(aes(x=subsite2, y=dailyflux_logsign, fill=status, group=paste0(subsite)))+
  geom_boxplot()+
  geom_point(data = co2log_outliers, col="red")+
  facet_wrap(~casepilot,scales="free")

#Extreme outliers (3*IQR, by status,season,casepilot, subsite, strata) could be removed, if removed, sample_weights should be adjusted accordingly to retain strata-representativity (re-normalize). 
#SHOULD WE REMOVE OUTLIERS? 
#DECISSION: NO. Our "outliers" represent true variability in field conditions. Even when grouped by strata,  data-driven outlier removal should NOT be done (as we cannot account for small scale variability in our dataset) 


## ========== 1. Untransformed Model ==========
model_untrans <- lmer(dailyflux ~ status * season * casepilot + (1 | subsite), 
                      data = data_co2_all)

## ========== 2. Signed Log Transform Model ==========
model_logsign <- lmer(dailyflux_logsign ~ status * season * casepilot + (1 | subsite),
                      data = data_co2_all)

## ========== 3. Signed Root Transform Model ==========
model_rootsign <- lmer(dailyflux_rootsign ~ status * season * casepilot + (1 | subsite),
                       data = data_co2_all)

## ========== 4. Evaluate Model requirements ==========
#This function produces diagnostics for the LMM models, either the general ones (including casepilot as factor) or the per-casepilot LMMs
check_model_diagnostics <- function(model, model_name, data) {
  cat("----- Diagnostics for", model_name, "-----\n")
  
  print(ggResidpanel::resid_panel(model, qqbands = TRUE, plots = c("qq", "resid", "hist", "fitted")))
  
  res <- residuals(model)
  
  # Residuals vs status
  if ("status" %in% names(data)) {
    p1 <- ggplot(data.frame(status = data$status, resid = res), aes(x = status, y = resid)) +
      geom_boxplot() +
      labs(title = paste("Residuals by Status:", model_name), y = "Residuals")
    print(p1)
  }
  
  # Residuals vs casepilot
  if ("casepilot" %in% names(data) && nlevels(factor(data$casepilot)) > 1) {
    p2 <- ggplot(data.frame(casepilot = data$casepilot, resid = res), aes(x = casepilot, y = resid)) +
      geom_boxplot() +
      labs(title = paste("Residuals by Casepilot:", model_name), y = "Residuals")
    print(p2)
  }
  
  # Residuals vs season
  if ("season" %in% names(data)) {
    p3 <- ggplot(data.frame(season = data$season, resid = res), aes(x = season, y = resid)) +
      geom_boxplot() +
      labs(title = paste("Residuals by Season:", model_name), y = "Residuals")
    print(p3)
  }
  
  # Levene's Tests
  if ("status" %in% names(data) && nlevels(factor(data$status)) > 1) {
    cat("Levene's Test for residual variance by status:\n")
    print(leveneTest(res ~ data$status))
  }
  
  if ("casepilot" %in% names(data) && nlevels(factor(data$casepilot)) > 1) {
    cat("Levene's Test for residual variance by casepilot:\n")
    print(leveneTest(res ~ data$casepilot))
  }
  
  if ("season" %in% names(data) && nlevels(factor(data$season)) > 1) {
    cat("Levene's Test for residual variance by season:\n")
    print(leveneTest(res ~ data$season))
  }
  
  # Shapiro-Wilk test
  cat("Shapiro-Wilk test for residual normality:\n")
  print(shapiro.test(res))
  
  # R², AIC, BIC
  r2 <- performance::r2_nakagawa(model)
  aic_val <- AIC(model)
  bic_val <- BIC(model)
  
  # 10-fold CV RMSE
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
diag_untrans <- check_model_diagnostics(model_untrans, "Untransformed", data_co2_all)
#Untransformed Observations: clear non-normality of residuals, heterocedasticity at season and casepilot level (not at status level)
diag_logsign <- check_model_diagnostics(model_logsign, "Signed Log", data_co2_all)
#Logtransformed Observations: Residuals show reasonable normality (visual). Heterocedasticity at season and casepilot level (not at status level), OK heterocedasticity is expected. BEST
diag_rootsign <- check_model_diagnostics(model_rootsign, "Signed Root", data_co2_all)
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

#Set back-transformation function accordingly: inv_logsign or inv_rootsign
inv_function<- inv_logsign

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


## ---- 5.3. EMMs  ----
#Estimated Marginal Means (emms)


#EMMs of status, for each combination of season and case-pilot: Within each season and casepilot, how do the statuses compare?
emm_status<-emmeans(best_model, ~ status | season * casepilot)
summary(emm_status)

#The above call of emmeans is the option to be used for post-hoc comparisons of status effects within ecological and seasonal contexts, results are equivalent to emm_all <- emmeans(best_model,~ status * season * casepilot) but grouped differently, so that calling pairs(emm_status) with default options will already provide the meaningful comparisons. If I was to call pairs(emm_all) it would return every pairwise comparison across all combinations fo the 3 factors, which is meaningless. 

#Back-transform EMMs
emm_status_bt<- as.data.frame(emm_status) %>% 
  mutate(across(c(emmean, SE,lower.CL,upper.CL), inv_function))


#Plot EMMs 
ggplot(emm_status_bt, aes(x = season, y = emmean, fill = status)) +
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

#Compare EEMs vs weigthed averages (calcualte SE of weighted average)
weighted_avgs_og <- data_co2_all %>%
  group_by(status, casepilot, season) %>%
  summarise(
    weighted_mean = sum(dailyflux * sample_weight) / sum(sample_weight),
    sum_w = sum(sample_weight),
    sum_w2 = sum(sample_weight^2),
    # Temporarily compute mean for variance calculation
    mean_w = sum(dailyflux * sample_weight) / sum(sample_weight),
    # Weighted variance
    var_w = sum(sample_weight * (dailyflux - mean_w)^2) / sum(sample_weight),
    # Effective sample size
    n_eff = (sum_w)^2 / sum_w2,
    # Standard error
    se = sqrt(var_w) / sqrt(n_eff),
    .groups = "drop"
  ) %>%
  # Keep only relevant output columns
  dplyr::select(status, casepilot, season, weighted_mean, se)

#Plot weighted averages (original scale) 
ggplot(weighted_avgs_og, aes(x = season, y = weighted_mean, fill = status)) +
  geom_col(position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = weighted_mean - se, ymax = weighted_mean + se),
                position = position_dodge(width = 0.8), width = 0.2) +
  facet_wrap(~ casepilot,scales="free") +
  labs(
    title = "Weighted avg by Status, Season, and Casepilot",
    y = expression(Average~CO[2]~dailyflux~(g~CO2[eq]~d^-1~m^-2)),
    x = "Season",
    col = "Status"
  ) +
  theme_minimal()

#PLot back transformed EMMs vs Weighted averages (original scale)
emm_status_bt %>% 
  left_join(weighted_avgs_og) %>% 
  ggplot(aes(x=weighted_mean,y=emmean))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  ggtitle("BT EMMs vs Weighted avg of untransformed\n(same groups)")
#There is a large deviation because the EMMs are calculated based on transformed values (as they should), and the back transformed average is not the same as the average of back transformed (original scale)

#JUST FOR CONSISTENCY: if we calculate the weighted averages on the transformed values and then back transform the averages, the results are equal to those from the model emmeans .



## ---- 5.4. Post-hoc on EMMs  ----
#All statistics, including post-hoc must be performed in the same scale as the model (i.e. sign-log transformed), only  back-transform for plots, interpretation and reporting effect sizes. 

# EMM comparison for status by season and casepilot comparisons
pairs(emm_status, adjust="tukey")

#Only significant differences for a few cases of DANUBE DELTA: 
pairs(emm_status, adjust="tukey") %>% as.data.frame() %>% 
  filter(p.value<0.05)


## ---- 5.5. Model complexity  ----

#Model too complex? 

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
                     data=data_co2_all)
anova(model_no3way, best_model)
#Current more complex model significantly improves the fit, how much added explanatory power? is it worth it?


#BASED ON Type III ANOVA Results: 
# ✅ Interpretation of the Type III ANOVA Results
# The main effect of status is not significant, but its interactions with season and casepilot are.
# The main effect of casepilot is highly significant, which aligns with your expectation that these are ecologically distinct systems.
# 🔍 What this suggests:
#   The effect of status is not consistent across all casepilots or seasons — it depends on context.
# The strong casepilot effect may indeed obscure or dilute the status effect when analyzed in a pooled model, especially if the direction or magnitude of the status effect varies across casepilots.


#THEREFORE: Modeling Each Casepilot Separately is the better option.

# Interactions are significant and complex (as in your case).
# We expect heterogeneity in how status affects emissions across ecosystems.
# We want to avoid overfitting or interpretational complexity from 3-way interactions.
# By modeling each casepilot individually:

# We simplify the model (no 3-way interactions).
# We gain clarity on how status behaves in each ecological context.
# We can later compare or summarize these effects across casepilots (e.g., via forest plots or meta-analytic summaries).





#___________________------
#MISCELANEA------------
#TEST_MODELDECISSION_1-------
#CHUNK generated with Co-pilot:

# ================================:
# GHG Flux Modeling Decision Tree
# ================================:

# Load required packages
library(tidyverse)
library(lme4)
library(glmmTMB)
library(bestNormalize)
library(car)
library(performance)
library(caret)
library(DHARMa)
library(ggResidpanel)


# 1. Subset Data

# Replace with your casepilot and gas of interest
data_sub <- data %>% filter(casepilot == "CA", gas == "CO2")


# 2. Initial Distribution Checks

# Histogram and Q-Q plot
ggplot(data_sub, aes(x = dailyflux)) + 
  geom_histogram(bins = 30) + 
  ggtitle("Histogram of Raw Flux")

qqnorm(data_sub$dailyflux)
qqline(data_sub$dailyflux)

# Shapiro-Wilk test
shapiro.test(data_sub$dailyflux)


# 3. Try Transformations

# Signed log
data_sub <- data_sub %>%
  mutate(dailyflux_logsign = sign(dailyflux) * log10(abs(dailyflux) + 1),
         dailyflux_rootsign = sign(dailyflux) * sqrt(abs(dailyflux)))

# Box-Cox (only for positive values)
if (all(data_sub$dailyflux > 0)) {
  bc <- boxcox(dailyflux ~ 1, data = data_sub)
}


# 4. Fit lmer Models (with weights)

model_untrans <- lmer(dailyflux ~ status * season + (1 | subsite), 
                      data = data_sub, weights = sample_weight)

model_logsign <- lmer(dailyflux_logsign ~ status * season + (1 | subsite), 
                      data = data_sub, weights = sample_weight)

model_rootsign <- lmer(dailyflux_rootsign ~ status * season + (1 | subsite), 
                       data = data_sub, weights = sample_weight)


# 5. Residual Diagnostics Function

check_model_diagnostics <- function(model, model_name, data) {
  cat("----- Diagnostics for", model_name, "-----\n")
  print(ggResidpanel::resid_panel(model, qqbands = TRUE, plots = c("qq", "resid", "hist", "fitted")))
  
  res <- residuals(model)
  
  # Levene's Test
  if ("status" %in% names(data) && nlevels(factor(data$status)) > 1) {
    cat("Levene's Test by status:\n")
    print(leveneTest(res ~ data$status))
  }
  if ("season" %in% names(data) && nlevels(factor(data$season)) > 1) {
    cat("Levene's Test by season:\n")
    print(leveneTest(res ~ data$season))
  }
  
  # Shapiro-Wilk
  cat("Shapiro-Wilk test:\n")
  print(shapiro.test(res))
  
  # R2, AIC, BIC
  r2 <- performance::r2_nakagawa(model)
  aic_val <- AIC(model)
  bic_val <- BIC(model)
  
  # RMSE via 10-fold CV
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


# 6. Run Diagnostics

diag_untrans <- check_model_diagnostics(model_untrans, "Untransformed", data_sub)
diag_logsign <- check_model_diagnostics(model_logsign, "Signed Log", data_sub)
diag_rootsign <- check_model_diagnostics(model_rootsign, "Signed Root", data_sub)


# 7. Compare Models

results <- tibble(
  model = c("Untransformed", "Signed Log", "Signed Root"),
  marginal_r2 = c(diag_untrans$r2$R2_marginal, diag_logsign$r2$R2_marginal, diag_rootsign$r2$R2_marginal),
  conditional_r2 = c(diag_untrans$r2$R2_conditional, diag_logsign$r2$R2_conditional, diag_rootsign$r2$R2_conditional),
  AIC = c(diag_untrans$AIC, diag_logsign$AIC, diag_rootsign$AIC),
  BIC = c(diag_untrans$BIC, diag_logsign$BIC, diag_rootsign$BIC),
  RMSE_CV = c(diag_untrans$RMSE_CV, diag_logsign$RMSE_CV, diag_rootsign$RMSE_CV)
) %>%
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

# 8. If All Models Violate Assumptions → Try glmmTMB

# Choose family based on distribution
family_choice <- gaussian()  # or Gamma(link = "log") for CH4/N2O

model_glmmTMB <- glmmTMB(
  dailyflux ~ status * season + (1 | subsite),
  data = data_sub,
  weights = sample_weight,
  family = family_choice,
  dispformula = ~ season  # Optional: model heteroscedasticity
)

# Residual diagnostics
sim_output <- simulateResiduals(model_glmmTMB)
plot(sim_output)

# Model performance
print(performance::r2(model_glmmTMB))


#TEST2_MODELDECISSIONTREE---------
# GHG Flux Modeling Decision Tree Template
# Author: Copilot
# Description: Systematic modeling workflow with assumption checks and summary table

library(tidyverse)
library(lme4)
library(glmmTMB)
library(DHARMa)
library(ggResidpanel)
library(car)
library(bestNormalize)
library(performance)
library(caret)

# Initialize summary table
summary_table <- tibble(
  casepilot = character(),
  gas = character(),
  normality_data = character(),
  normality_residuals = character(),
  homoscedasticity_status = character(),
  homoscedasticity_season = character(),
  linearity = character(),
  independence = "assumed via random effects",
  recommended_model = character(),
  notes = character()
)

# Define transformation functions
signed_log <- function(x) sign(x) * log10(abs(x) + 1)
signed_root <- function(x) sign(x) * sqrt(abs(x))

# Example: Loop through one casepilot × gas combination
# Replace with your actual loop or manual selection
casepilot_id <- "CA"
gas_type <- "CO2"
data_sub <- data %>% filter(casepilot == casepilot_id, gas == gas_type)

# Step 1: Explore distribution of raw data
hist(data_sub$dailyflux, main = "Histogram of Raw Flux", xlab = "dailyflux")
qqnorm(data_sub$dailyflux); qqline(data_sub$dailyflux)
shapiro_raw <- shapiro.test(data_sub$dailyflux)
normality_data <- ifelse(shapiro_raw$p.value > 0.05, "yes", "no")

# Step 2: Apply transformations
data_sub <- data_sub %>%
  mutate(
    dailyflux_logsign = signed_log(dailyflux),
    dailyflux_rootsign = signed_root(dailyflux)
  )

# Step 3: Fit lmer models
model_untrans <- lmer(dailyflux ~ status * season + (1 | subsite), data = data_sub, weights = sample_weight)
model_logsign <- lmer(dailyflux_logsign ~ status * season + (1 | subsite), data = data_sub, weights = sample_weight)
model_rootsign <- lmer(dailyflux_rootsign ~ status * season + (1 | subsite), data = data_sub, weights = sample_weight)

# Step 4: Residual diagnostics
check_residuals <- function(model, data, model_name) {
  res <- residuals(model)
  fitted_vals <- fitted(model)
  qq <- shapiro.test(res)
  normality <- ifelse(qq$p.value > 0.05, "yes", "no")
  
  levene_status <- tryCatch({
    leveneTest(res ~ data$status)
  }, error = function(e) NULL)
  
  levene_season <- tryCatch({
    leveneTest(res ~ data$season)
  }, error = function(e) NULL)
  
  hom_status <- ifelse(!is.null(levene_status) && levene_status$`Pr(>F)`[1] > 0.05, "yes", "no")
  hom_season <- ifelse(!is.null(levene_season) && levene_season$`Pr(>F)`[1] > 0.05, "yes", "no")
  
  linearity <- ifelse(any(grepl("pattern", capture.output(plot(res ~ fitted_vals)))), "no", "yes")
  
  list(
    normality = normality,
    hom_status = hom_status,
    hom_season = hom_season,
    linearity = linearity
  )
}

diag_untrans <- check_residuals(model_untrans, data_sub, "Untransformed")
diag_logsign <- check_residuals(model_logsign, data_sub, "Signed Log")
diag_rootsign <- check_residuals(model_rootsign, data_sub, "Signed Root")

# Step 5: Choose best model based on diagnostics
model_scores <- tibble(
  model = c("untrans", "logsign", "rootsign"),
  normality = c(diag_untrans$normality, diag_logsign$normality, diag_rootsign$normality),
  hom_status = c(diag_untrans$hom_status, diag_logsign$hom_status, diag_rootsign$hom_status),
  hom_season = c(diag_untrans$hom_season, diag_logsign$hom_season, diag_rootsign$hom_season),
  linearity = c(diag_untrans$linearity, diag_logsign$linearity, diag_rootsign$linearity)
)

# Select best model (all "yes" preferred)
model_scores <- model_scores %>%
  mutate(score = rowSums(across(normality:linearity, ~ . == "yes"))) %>%
  arrange(desc(score))

best_model <- model_scores$model[1]
recommended_model <- ifelse(model_scores$score[1] == 4, paste("lmer:", best_model), "glmmTMB")

# Step 6: Append to summary table
summary_table <- summary_table %>%
  add_row(
    casepilot = casepilot_id,
    gas = gas_type,
    normality_data = normality_data,
    normality_residuals = model_scores$normality[1],
    homoscedasticity_status = model_scores$hom_status[1],
    homoscedasticity_season = model_scores$hom_season[1],
    linearity = model_scores$linearity[1],
    recommended_model = recommended_model,
    notes = ifelse(recommended_model == "glmmTMB", "fallback due to assumption violations", "assumptions met")
  )

print(summary_table)

#TEST3_MODELDECISSIONTREEE---------
#SImilar to test2 but looped.
#Need to adapt GLMM distribution family depending on data distribution


library(tidyverse)
library(lme4)
library(glmmTMB)
library(DHARMa)
library(ggResidpanel)
library(car)
library(caret)
library(bestNormalize)
library(performance)

# Initialize summary table
summary_results <- tibble()

# Define inverse transformations
inv_logsign <- function(x) sign(x) * (10^abs(x) - 1)
inv_rootsign <- function(x) sign(x) * (abs(x)^2)

# Loop over all combinations of casepilot and GHG
for (cp in unique(data4models$casepilot)) {
  for (ghg in unique(data4models$ghgspecies)) {
    
    cat("\n\n==================== CASEPILOT:", cp, "| GHG:", ghg, "====================\n")
    
    data_sub <- data4models %>% filter(grepl(cp, casepilot)& ghgspecies == ghg) %>% filter(!is.na(dailyflux))
    
    # Skip if too few observations
    if (nrow(data_sub) < 20) next
    
    # Step 1: Distribution checks
    response <- data_sub$dailyflux
    shapiro_raw <- shapiro.test(response)
    skewness <- e1071::skewness(response)
    
    # Step 2: Transformations
    data_sub <- data_sub %>%
      mutate(
        dailyflux_logsign = sign(dailyflux) * log10(abs(dailyflux) + 1),
        dailyflux_rootsign = sign(dailyflux) * sqrt(abs(dailyflux))
      )
    
    # Step 3: Fit lmer models
    model_untrans <- tryCatch(lmer(dailyflux ~ status * season + (1 | subsite), data = data_sub, weights = sample_weight), error = function(e) NULL)
    model_logsign <- tryCatch(lmer(dailyflux_logsign ~ status * season + (1 | subsite), data = data_sub, weights = sample_weight), error = function(e) NULL)
    model_rootsign <- tryCatch(lmer(dailyflux_rootsign ~ status * season + (1 | subsite), data = data_sub, weights = sample_weight), error = function(e) NULL)
    
    # Step 4: Diagnostics
    check_model <- function(model, name) {
      if (is.null(model)) return(NULL)
      res <- residuals(model)
      shapiro_res <- shapiro.test(res)$p.value
      levene_status <- tryCatch(leveneTest(res ~ data_sub$status)$"Pr(>F)"[1], error = function(e) NA)
      levene_season <- tryCatch(leveneTest(res ~ data_sub$season)$"Pr(>F)"[1], error = function(e) NA)
      r2 <- tryCatch(performance::r2_nakagawa(model), error = function(e) NULL)
      aic_val <- tryCatch(AIC(model), error = function(e) NA)
      bic_val <- tryCatch(BIC(model), error = function(e) NA)
      rmse_cv <- tryCatch({
        ctrl <- trainControl(method = "cv", number = 10)
        cv_model <- train(
          x = model.matrix(model)[, -1],
          y = getME(model, "y"),
          method = "lm",
          trControl = ctrl
        )
        cv_model$results$RMSE
      }, error = function(e) NA)
      
      tibble(
        casepilot = cp,
        gas = ghg,
        model = name,
        shapiro_resid_p = shapiro_res,
        levene_status_p = levene_status,
        levene_season_p = levene_season,
        marginal_r2 = if (!is.null(r2)) r2$R2_marginal else NA,
        conditional_r2 = if (!is.null(r2)) r2$R2_conditional else NA,
        AIC = aic_val,
        BIC = bic_val,
        RMSE_CV = rmse_cv
      )
    }
    
    res_untrans <- check_model(model_untrans, "Untransformed")
    res_logsign <- check_model(model_logsign, "Signed Log")
    res_rootsign <- check_model(model_rootsign, "Signed Root")
  
    summary_results <- bind_rows(summary_results, res_untrans, res_logsign, res_rootsign)
  }
}



#CH4---------

#Check distribution families that best fit the actual data: no transformation will produce normal-distributed data. 

# Load required packages
library(ggplot2)
library(MASS)
library(statmod)  # for dinvgauss


data4models %>% filter(ghgspecies=="gwp_ch4") %>%
  filter(!is.na(dailyflux)) %>%
ggplot(aes(x=dailyflux, col=dailyflux>0))+
  geom_histogram()+
  facet_wrap(~dailyflux>0, scales="free")

y<-data4models %>% filter(ghgspecies=="gwp_ch4") %>%
  filter(!is.na(dailyflux)) %>%
  mutate(dailyflux=dailyflux+(abs(min(dailyflux)*1.1))) %>% 
  pull(dailyflux)

min(y)
hist(y,breaks = 70)

# Estimate parameters for each distribution
params_norm <- list(mean = mean(y), sd = sd(y))
params_gamma <- fitdistr(y, "gamma")$estimate
params_lnorm <- fitdistr(y, "lognormal")$estimate

# Plot histogram with overlaid densities
ggplot(data.frame(y), aes(x = y)) +
  geom_histogram(aes(y = ..density..), bins = 50, fill = "gray90", color = "black") +
  stat_function(fun = dnorm, args = params_norm, color = "blue", size = 1.1, linetype = "solid", aes(linetype = "Normal")) +
  stat_function(fun = dgamma, args = list(shape = params_gamma["shape"], rate = params_gamma["rate"]),
                color = "red", size = 1.1, linetype = "dashed", aes(linetype = "Gamma")) +
  stat_function(fun = dlnorm, args = list(meanlog = params_lnorm["meanlog"], sdlog = params_lnorm["sdlog"]),
                color = "green4", size = 1.1, linetype = "dotdash", aes(linetype = "Lognormal")) +
  scale_linetype_manual(name = "Distributions",
                        values = c("Normal" = "solid",
                                   "Gamma" = "dashed",
                                   "Lognormal" = "dotdash")) +
  labs(title = "Histogram with Fitted Continuous Distributions", x = "y", y = "Density") +
  theme_minimal()






#Exploratory plots-----
library(ggforce)
# install.packages("ggforce")
#CA:

#Orignal scale (non-tranformed):
filter(casepilot == "RI", comparison == "status") %>%
  ggplot(., aes(x = status, y = emmean_gwp, label = group_letter)) +
  geom_errorbar(aes(ymin = emmean_gwp - SE_gwp, ymax = emmean_gwp + SE_gwp),
                linewidth = 0.2, width = 0.2) +
  geom_point(shape = 15, size = 4, aes(col=status)) +
  geom_text(nudge_x = 0.2, color = "black") +
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw() +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  labs(title = expression(Status~effect~on~NEE~(CO[2]+CH[4])),
       subtitle = paste0("Ria d'Aveiro"), 
       y=expression(EMM~(g~CO[2~eq]~d^-1~m^-2)),
       x="Conservation status",
       col=paste0("Status"))

data4models %>% 
  filter(casepilot=="RI") %>% 
  filter(!is.na(dailyflux)) %>% 
  mutate(dailyflux_logsign=logsign(dailyflux)) %>% 
  filter(ghgspecies%in%c("gwp_co2andch4")) %>% 
  mutate(ghgspecies=factor(ghgspecies, levels = c("gwp_co2","gwp_ch4","gwp_co2andch4"))) %>% 
  ggplot(aes(x=status, y=dailyflux, fill=status))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_boxplot(outliers = F)+
  theme_bw() +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0))+
  labs(title = expression(Status~effect~on~NEE~(CO[2]+CH[4])),
       subtitle = paste0("Ria d'Aveiro"), 
       y=expression(NEE~(g~CO[2~eq]~d^-1~m^-2)),
       x="Conservation status",
       col=paste0("Status"))

ggsave(filename="CO2plusCH4_boxplot_status_RI.png",
       device = "png",path = plots_path,dpi = 400,
       width = 115, height = 100,units = "mm")

#Log-sign scale: 

#set appropriate breaks:
breaks_real <- c(-64, -32, -16, -8, -4, -2, -1, 0, 1, 2, 4, 8, 16, 32, 64, 128)
breaks_trans <- logsign(breaks_real)

data4models %>% 
  filter(casepilot=="RI") %>% 
  filter(!is.na(dailyflux)) %>% 
  mutate(dailyflux_logsign=logsign(dailyflux)) %>% 
  filter(ghgspecies%in%c("gwp_co2","gwp_ch4","gwp_co2andch4")) %>% 
  mutate(ghgspecies=factor(ghgspecies, levels = c("gwp_co2","gwp_ch4","gwp_co2andch4"))) %>% 
  ggplot(aes(x=status, y=dailyflux_logsign, fill=status))+
  geom_boxplot()+
  scale_y_continuous(minor_breaks = NULL,
    name = "Daily flux (log-sign scale)",
    breaks = breaks_trans,  # choose appropriate range
    labels = breaks_real  # inverse labels to interpret values
  ) +
  facet_wrap(~ghgspecies, scales="free")


#USING geom sina:
data4models %>% 
  filter(casepilot=="RI") %>% 
  filter(!is.na(dailyflux)) %>% 
  mutate(dailyflux_logsign=logsign(dailyflux)) %>% 
  filter(ghgspecies%in%c("gwp_co2","gwp_ch4","gwp_co2andch4")) %>% 
  mutate(ghgspecies=factor(ghgspecies, levels = c("gwp_co2","gwp_ch4","gwp_co2andch4"))) %>% 
  ggplot(aes(x=status, y=dailyflux_logsign,col=status))+
  geom_hline(yintercept = logsign(0))+
  geom_violin(scale="width")+
  geom_sina(aes(col=status),scale=F)+
  scale_y_continuous(minor_breaks = NULL,
    name = "Daily flux (log-sign scale)",
    breaks = breaks_trans,  # choose appropriate range
    labels = breaks_real  # inverse labels to interpret values
  ) +
  facet_wrap(~ghgspecies, scales="free")


data4models %>% 
  filter(casepilot=="CA") %>% 
  filter(!is.na(dailyflux)) %>% 
  mutate(dailyflux_logsign=logsign(dailyflux)) %>% 
  filter(ghgspecies%in%c("gwp_co2","gwp_ch4","gwp_co2andch4")) %>% 
  mutate(ghgspecies=factor(ghgspecies, levels = c("gwp_co2","gwp_ch4","gwp_co2andch4"))) %>% 
  ggplot(aes(x=status, y=dailyflux_logsign))+
  geom_boxplot(aes(col=status))+
  geom_sina(aes(col=status),scale=F)+
  scale_y_continuous(minor_breaks = NULL,
    name = "Daily flux (log-sign scale)",
    breaks = breaks_trans,  # choose appropriate range
    labels = breaks_real  # inverse labels to interpret values
  ) +
  facet_wrap(~ghgspecies, scales="free")



#GWP_CO2: all seasons combined
daily_ghg_plotcode %>%
  mutate(subsite2=substr(subsite,4,5)) %>% 
  group_by(casepilot, status, subsite,subsite2) %>%
  filter(!is.na(gwp_co2)) %>% 
  ggplot(aes(x=subsite2, y=sign(gwp_co2)*log(abs(gwp_co2)+1), fill=status, group=paste0(subsite)))+
  geom_boxplot()+
  facet_wrap(~casepilot,scales="free")

#GWP_CH4: all seasons combined
daily_ghg_plotcode %>%
  mutate(subsite2=substr(subsite,4,5)) %>% 
  group_by(casepilot, status, subsite,subsite2) %>%
  filter(!is.na(gwp_co2)) %>% 
  ggplot(aes(x=subsite2, y=sign(gwp_ch4)*log(abs(gwp_ch4)+1), fill=status, group=paste0(subsite)))+
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
  filter(!is.na(gwp_co2andch4)) %>% 
  ggplot(aes(x=subsite,y=gwp_co2andch4, fill= status, group=sampling))+
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
  filter(!is.na(gwp_co2andch4)) %>% 
  ggplot(aes(x=subsite,y=gwp_co2andch4, fill= status, group=sampling))+
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





daily_ghg_plotcode %>%
  filter(casepilot=="RI") %>% 
  filter(!is.na(gwp_co2)) %>% 
  ggplot(aes(x=subsite,y=gwp_co2, fill= status, group=sampling))+
  geom_boxplot()+
  geom_violin()

valid_ghg_transpdark %>% filter(grepl("RI",plotcode)) %>% 
  filter(ghg=="co2") %>% 
  separate(plotcode, into = c("season", "casepilot", "subsite", "plotnum"),sep = "-",remove = F) %>% 
  mutate(status=substr(subsite,1,1)) %>% 
  ggplot(aes(x=subsite, y=dark, fill = status))+
  geom_boxplot()+
  geom_point()

valid_ghg_transpdark %>% filter(grepl("RI",plotcode)) %>% 
  filter(ghg=="co2") %>% 
  separate(plotcode, into = c("season", "casepilot", "subsite", "plotnum"),sep = "-",remove = F) %>% 
  mutate(status=substr(subsite,1,1)) %>% 
  ggplot(aes(x=subsite, y=transparent, fill = status))+
  geom_boxplot()+
  geom_point()
