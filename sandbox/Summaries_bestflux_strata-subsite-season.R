
# ---
# Authors: Miguel Cabrera
# Project: "RESTORE4Cs"
# date: "May 2025"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
#This script uses best.flux simplified csv files re-calculated from aquaGHG approach to produce per-strata, per-subsite

#Antonio necesita S3 per-strata&per-subsite fluxes



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


# ---- Directories ----
dropbox_root <- "C:/Users/Miguel/Dropbox/RESTORE4Cs - Fieldwork/Data" # You have to make sure this is pointing to the write folder on your local machine

data_path<- paste0(dropbox_root, "/GHG/Processed data/computed_flux/Updated_bestflux_aquaGHG/")


#Import bestflux
ch4<- read.csv(file = paste0(data_path, "ch4_bestflux.csv"))
co2<- read.csv(file = paste0(data_path, "co2_bestflux.csv"))



#Per strata summaries CH4 (average, no lightcondition taken into account)

ch4_plusbasic<- ch4 %>% 
  select(UniqueID, best.flux, best_model,best_model_flags) %>% 
  separate(UniqueID, into = c("season", "case_pilot", "subsite_code","plotnum", "strata", "light", "utchourminute"), sep = "-", remove = F) %>% 
  mutate(sampling=paste(season, case_pilot, subsite_code, sep = "-"))




ch4_strata<- ch4_plusbasic %>% 
  select(season, case_pilot, subsite_code, sampling, strata, best.flux) %>% 
  group_by(season, case_pilot, subsite_code, sampling, strata) %>% 
  summarise(avg_flux=mean(best.flux, na.rm=T),
            sd_flux=sd(best.flux, na.rm=T), 
            n_flux=sum(!is.na(best.flux)),
            n_discard=sum(is.na(best.flux)),
            min_flux=min(best.flux, na.rm = T),
            max_flux=max(best.flux, na.rm = T),
            median_flux=median(best.flux, na.rm = T))



write.csv(ch4_strata,file = paste0(data_path, "CH4avgflux_nmolCH4persquaremeterpersecond_summary_per_strata.csv"), row.names = F)

ch4 %>% 
  group_by(best_model) %>% 
  summarise(n=n())

co2 %>% 
  group_by(best_model) %>% 
  summarise(n=n())


