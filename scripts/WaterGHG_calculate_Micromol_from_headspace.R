#Water GHG from headspace

#Description----
# ---
# Authors: Miguel Cabrera Brufau
# Project: "RESTORE4Cs"
# date: "September 2025"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
#MOTIVATION: We need to obtain from atmospheric and equilibrated headspaces concentrations the original ghg concentration in the water used for equilibration, in Molar units. 


#FUNCTIONING: This script uses the atmospheric and headspace exetainer ghg concentrations to calulate the molar GHG concentration in the water. It requires ppm-ghg concentrations (atm and headspace) as well as the equilibration water temperature. 


#STEPS: 

#TODO:  Import, format and filter headspace and atmospheric GHG data. (from GC and from Licor-injections) into a single data frame of format:  subsite_ID, exetainer_ID, sample identity (atm|headspace), co2, ch4, n2o=NA, method:"GC-derived"|"Licor-derived"

#TODO: Inspect and remove weird atmospheric data, impute with similar average atm if needed. 
#TODO: Import accessory data (water temperatures from waterteam)
#TODO: Use Jorge functions to calcuate water-GHG concentrations. 
#TODO: Save water-GHG concentrations with plot identifyer (UniqueID from in-situGHG dataset)


rm(list = ls())

# ---- Directories ----
dropbox_root <- "C:/Users/Miguel/Dropbox/RESTORE4Cs - Fieldwork/Data/"

fieldsheet_path <- paste0(dropbox_root,"GHG/Fieldsheets")
headspace_path<- paste0(dropbox_root,"GHG/Headspaces/")
gc_path<- paste0(dropbox_root,"GHG/GC_data/")


#---Packages and functions -----

library(tidyverse)
library(readxl)
library(broom)
library(ggpmisc)
library(egg)

repo_root <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
files.sources = list.files(path = paste0(repo_root,"/functions"), full.names = T)
for (f in files.sources){source(f)}

#Headspace functions are stored in script within repo functions folder, also in:  
# source("https://raw.githubusercontent.com/JorgeMonPe/Functions_headspace_calculation/main/Functions_gas_concentration.R")


#1. Import and format-----


##GC-data----

gc_dat<- read.csv(paste0(gc_path, "allsamples_ppm_noblankcorrection.csv")) %>% 
  filter(exe_type%in%c("atmosphere", "headspace"))

names(gc_dat)

##Licor-data----

licor_dat<- read.csv(paste0(headspace_path,"S4_Headspace&atm_co2_ch4_n2o.csv"))

names(licor_dat)

##In-situ atm data-----
#To fill-in co2 and ch4 concentrations from S4 (not measured with CO2&CH4 licor, atmospheric exetainers were only analyzed with N2O-licor)

atm_insitu_dat<- read.csv(paste0(headspace_path, "in_situ_atm_co2&ch4.csv")) %>% 
  select(subsite,ghgspecies, avg) %>% 
  rename(subsite_ID=subsite) %>% 
  pivot_wider(names_from = ghgspecies,values_from = avg) %>% 
  mutate(ppm_ch4_insitu=ch4_ppb/1000) %>% 
  rename(ppm_co2_insitu=co2_ppm) %>% 
  select(subsite_ID, ppm_co2_insitu, ppm_ch4_insitu)
#DO NOT USE in-situ data for GC-data vials. (calfactor has bias and this will impact the ghgwater calculation)


##Water-temp data----
#Import water-temperature per Headspace
watertemp<- read.csv(paste0(headspace_path, "headspace_temperatures.csv")) %>% 
  rename(subsite_ID=subsiteID, exetainer_ID=exetainerID, temp_c=temp) %>% 
  select(subsite_ID, exetainer_ID, temp_c)


#2. Atm Summary ----
#Inspect and remove (impute if needed) weird atmospheric values. Produce single average atm GHGconc per subsiteID
#From GC data and from Licor data 
atm_gc<- gc_dat %>% 
  filter(exe_type=="atmosphere")

head(atm_gc)

atm_gc %>% 
  filter(campaign=="S3") %>% 
  filter(!exetainer_ID%in%discard_ch4_atm) %>% 
  ggplot(aes(x=subsite_ID, y=ppm_ch4))+
  geom_boxplot()+
  geom_point()+
  geom_label(aes(label=exetainer_ID))

atm_licor<- licor_dat %>% 
  filter(exe_type=="atmosphere")

atm_licor %>% 
  filter(campaign=="S4") %>% 
  filter(!exetainer_ID%in%discard_n2o_atm) %>% 
  ggplot(aes(x=subsite_ID, y=ppm_n2o))+
  geom_boxplot()+
  geom_point()+
  geom_label(aes(label=exetainer_ID))



#After inspection: fill the vectors with atmospheric samples to be discarded per each GHG 
discard_co2_atm<- c("S1-CU-P2-1", "S1-CU-R2-12","S2-CA-P1-1","S2-CA-R2-1","S2-RI-P1-1","S3-CA-P1-28","S3-CA-P2-26","S3-CA-A2-26","S3-CU-R1-11","S3-VA-R1-1","S4-CU-P2-1", "S4-CU-P2-2","S4-DU-R2-12")

discard_ch4_atm<- c("S1-CA-R1-22","S1-DA-A2-16","S1-DA-R2-3","S1-VA-A2-3","S1-VA-R1-12","S1-RI-R1-27","S2-CA-A2-18","S2-CU-A2-9","S2-RI-P1-14","S2-DA-R1-1","S3-CA-R1-1","S3-DA-A1-1","S3-DA-R1-1","S3-CU-A2-20","S3-CU-R1-11","S3-VA-P1-21","S3-CU-P2-8","S4-CU-P2-1", "S4-CU-P2-2")

discard_n2o_atm<- c()


#Join all exetainer analysis (and discard wrong values)
all_atm<- atm_gc %>% 
  select(-obs) %>% 
  mutate(ppm_n2o=NA_real_) %>% 
  rbind(atm_licor) %>% #Join GC data with atmospheric data
  #Set discard gas atm samples to NA
  mutate(ppm_co2=case_when(exetainer_ID%in%discard_co2_atm~NA_real_,
                           TRUE~ppm_co2),
         ppm_ch4=case_when(exetainer_ID%in%discard_ch4_atm~NA_real_,
                           TRUE~ppm_ch4),
         ppm_n2o=case_when(exetainer_ID%in%discard_n2o_atm~NA_real_,
                           TRUE~ppm_n2o))



#Obtain per-subsite atm average from exetainers 
atm_subsite_average<- all_atm %>% 
  group_by(subsite_ID) %>% 
  summarise(atm_ppm_co2=mean(ppm_co2, na.rm=T),
            atm_ppm_ch4=mean(ppm_ch4, na.rm=T),
            atm_ppm_n2o=mean(ppm_n2o, na.rm=T)) %>% 
  #Impute GC co2 average to subsite S1-CU-P2 (no data from vials)
  mutate(atm_ppm_co2=case_when(subsite_ID=="S1-CU-P2"~mean(atm_ppm_co2,na.rm=T),
                               TRUE~atm_ppm_co2)) %>% 
  #Join insitu co2 and ch4 to use for Licor vials (S4)
  full_join(atm_insitu_dat, by="subsite_ID") %>% 
  #Use insitu co2 and ch4 for licor S4 vials
  mutate(atm_ppm_co2=case_when(grepl("S4",subsite_ID)~ppm_co2_insitu,
                               TRUE~atm_ppm_co2),
         atm_ppm_ch4=case_when(grepl("S4", subsite_ID)~ppm_ch4_insitu,
                               TRUE~atm_ppm_ch4)) %>% 
  select(subsite_ID, atm_ppm_co2, atm_ppm_ch4, atm_ppm_n2o)


#3. Join all data----
#Join needed data into single data frame (one row per headspace)

gc_hs<- gc_dat %>% 
  filter(exe_type=="headspace") %>% 
  mutate(ppm_n2o=NA_real_) %>% 
  select(exetainer_ID, subsite_ID, exe_type, ppm_co2, ppm_ch4, ppm_n2o)


licor_hs<- licor_dat %>% 
  filter(exe_type=="headspace") %>% 
  select(exetainer_ID, subsite_ID, exe_type, ppm_co2, ppm_ch4, ppm_n2o)


#Exetainers to discard:
discard_exetainers<- c("S1-CU-P1-20", #Sample lost (errorGC)
                       "S1-CU-A2-21", #Sample lost (no vial)
                       "S1-DA-R1-7", #Sample lost (no vial)
                       "S1-DA-P1-11", #Sample lost (errorGC)
                       "S1-RI-A1-17", #Sample lost (no vial)
                       "S1-RI-A1-29", #Sample lost (no vial)
                       "S2-RI-A1-17", #Sample lost (no vial)
                       "S2-DA-A1-2" #NO temperature, only sample in this subsite, remove from dataset
)


#Join all necessary data for ghgWATER calculation: 
all_hs<- gc_hs %>% 
  rbind(licor_hs) %>% 
  left_join(atm_subsite_average, by="subsite_ID") %>% 
  left_join(watertemp,  by=c("exetainer_ID", "subsite_ID"))%>% 
  filter(!exetainer_ID%in%c(discard_exetainers))
  

#4. Calculate waterGHG-----
#Use Functions from Jorge to calculate waterGHG for every headspace. 
waterGHG<- all_hs %>% 
  mutate(co2_water_micromolar=nGHG_water_uM(Specie="CO2", GHG_ppmv = ppm_co2, Vol_H2O_HS = 30,Vol_air_HS = 30,T_Celsius = temp_c,GHG_atm_ppmv = atm_ppm_co2),
         ch4_water_micromolar=nGHG_water_uM(Specie="CH4", GHG_ppmv = ppm_ch4, Vol_H2O_HS = 30,Vol_air_HS = 30,T_Celsius = temp_c,GHG_atm_ppmv = atm_ppm_ch4),
         n2o_water_micromolar=nGHG_water_uM(Specie="CO2", GHG_ppmv = ppm_n2o, Vol_H2O_HS = 30,Vol_air_HS = 30,T_Celsius = temp_c,GHG_atm_ppmv = atm_ppm_n2o)) %>% 
  select(subsite_ID, exetainer_ID, exe_type,co2_water_micromolar,ch4_water_micromolar,n2o_water_micromolar)
  


waterGHG %>% 
  separate(subsite_ID, into = c("season","site","Subsite"),sep = "-",remove = F) %>% 
  mutate(subsite=paste(site,Subsite, sep = "-")) %>% 
  ggplot(aes(x=subsite, y=ch4_water_micromolar))+
  geom_point()+
  geom_boxplot()+
  facet_wrap(.~site, scales="free")


#5. Save data-----

#Final format and save waterGHG data with in-situGHG plot identifyers.

write.csv(waterGHG, file = paste0(headspace_path, "WaterGHG_from_headspace.csv"),row.names = F)


#____________________#----


#Testing of functions------

#How much does the water-concentration change per change of cal factor of (atm and hs ppm)

atm_ch4_conc<- 2.050

testdf<- data.frame(temp_c= rep(25,15),
                    vol_air=rep(30,15),
                    vol_water=rep(30,15),
                    factor=c(0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3)) %>% 
  mutate(atm_ppm=2.05*factor,
         headsp_ppm=atm_ppm*5)


testdf_results<- testdf %>% 
  mutate(micromolar_ch4_water=nGHG_water_uM(Specie = "CH4",GHG_ppmv = headsp_ppm,Vol_H2O_HS = vol_water, Vol_air_HS = vol_air,T_Celsius = temp_c, Patm_eq = 1, GHG_atm_ppmv = atm_ppm))


ggplot(testdf_results, aes(x=factor, y=micromolar_ch4_water))+
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "n")))+
  geom_point()

#SEEMS linear, no "compensation".  the calculated water concentration is directly dependent on the calibration factor of the source atmospheric and headspace ppm values





