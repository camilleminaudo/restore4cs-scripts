
# ---
# Authors: Miguel Cabrera
# Project: "RESTORE4Cs"
# date: "Dec 2023"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
# This script makes summaries of CO2 and CH4 fluxes calculated in the raw2flux for the Curonian Lagoon Open Water incubations. Using best.fluxes derived from goflux package. 

#Input: per-season CSV files with computed fluxes

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

# repo_root <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
# files.sources = list.files(path = paste0(repo_root,"/functions"), full.names = T)
# for (f in files.sources){source(f)}


# ---- Directories ----
dropbox_root <- "C:/Users/Miguel/Dropbox/RESTORE4Cs - Fieldwork/Data" # You have to make sure this is pointing to the write folder on your local machine
# datapath <- paste0(dropbox_root,"/GHG/RAW data")
# fieldsheetpath <- paste0(dropbox_root,"/GHG/Fieldsheets")
# loggerspath <- paste0(datapath,"/RAW Data Logger")
# RData_path <- paste0(dropbox_root,"/GHG/Processed data/RData/")
results_path <- paste0(dropbox_root,"/GHG/Processed data/computed_flux/")
plots_path <- paste0(dropbox_root,"/GHG/Processed data/plots_all_incubations/")


#Read all computed fluxes per season
season_csvs<- list.files(path = results_path, pattern = "^S")
dat<- data.frame()
for (i in season_csvs){
  s<- read.csv(paste0(results_path,i))
  
  dat<- rbind(dat,s)
}
rm(s, i)



#List UniqueID of incubations with weird patterns based on exploration of incubation plots.
weird_incubations_uniqueID<- c("s2-cu-p2-15-o-d-11:39",
                               "s1-cu-a2-2-o-d-06:52", 
                               "s3-cu-a2-15-o-d-13:07", 
                               "s3-cu-a1-1-o-d-06:39", 
                               "s1-cu-a1-2-o-d-07:28",
                               "s4-cu-r2-6-o-d-10:33")

#Filter open water from Curonian Lagoon, use only best.flux estimates, (transform CH4 from nmol to umol m-2 s-1), calculate emissions in CO2eq (mgCO2eq m-2 h-1), applying GWP of 28 for CH4 (1gCH4 = 28g CO2eq)

# calculate fluxes in mg m-2 h-1

cu<- dat %>% 
  separate(subsite, into = c("season", "site", "subsite_only"), remove = F) %>% 
  filter(site=="CU"&strata=="open water") %>% 
  filter(!UniqueID%in%wrong_incubations_uniqueID) %>% #Exclude weird-looking incubations
  mutate(CH4_best.flux=CH4_best.flux/1000, #get CH4 flux in umol m-2 s-1
         GWP_ch4=28*((CH4_best.flux*16.043)*1e-3*60*60), #28 * mgCH4 emited m-2 g-1 ~mgCO2eq m-2 h-1
         GWP_co2=((CO2_best.flux*44.009)*1e-3*60*60),#mgCO2 emited m-2 h-1
         GWP_total=GWP_ch4+GWP_co2) %>% #mgCO2eq m-2 h-1
  select(-c(CO2_LM.flux,CO2_HM.flux,CH4_diffusive_flux,CH4_ebullitive_flux,CH4_LM.flux,CH4_HM.flux,Area,Vtot,Tcham,Pcham,lightCondition,chamberType,duration))


write.csv(cu, file = paste0(results_path,"summaries_curonian_ow/Curonian_All_openwater_S1-S4_CO2_CH4_umolm-2_s-1_AND_gwpCO2eq_mgm-2h-1.csv"),row.names = F)

#Obtain fluxes in long format and in umol m-2 s-1
fluxes_long<- cu %>% 
  # filter(!UniqueID%in%wrong_incubations_uniqueID) %>% 
  select(season,subsite_only,UniqueID, CO2_best.flux, CH4_best.flux,GWP_ch4,GWP_co2,GWP_total) %>% 
  pivot_longer(cols = c(CO2_best.flux, CH4_best.flux,GWP_ch4,GWP_co2,GWP_total), names_to = "fluxtype", values_to = "flux")
  

ggplot(subset(fluxes_long, fluxtype%in%c("GWP_ch4","GWP_co2","GWP_total")), aes(x=season, y=flux))+
  geom_boxplot(outlier.size = 0.6)+
  facet_grid(fluxtype~subsite_only,scales="free")+
  theme_bw()+
  scale_y_continuous(name=expression(CO[2] * " eq flux (mg CO"[2] * " m"^-2 * " h"^-1 * ")"))

ggsave(filename = paste0(results_path,"summaries_curonian_ow/CU_ow_seasonal_GWP.jpg"),
       device = "jpg",units = "cm",width = 15,height = 15,)



ggplot(subset(fluxes_long, fluxtype%in%c("GWP_ch4","GWP_co2","GWP_total")), aes(x=subsite_only, y=flux))+
  geom_boxplot()+
  facet_grid(fluxtype~.,scales="free")+
  theme_bw()+
  scale_y_continuous(name=expression(CO[2] * " eq flux (mg CO"[2] * " m"^-2 * " h"^-1 * ")"))

ggsave(filename = paste0(results_path,"summaries_curonian_ow/CU_ow_year_GWP.jpg"),
       device = "jpg",units = "cm",width = 8,height = 10,)


#Summary per subsite and season
cu_sum<- fluxes_long %>% 
  filter(fluxtype%in%c("CO2_best.flux","CH4_best.flux","GWP_ch4","GWP_co2","GWP_total")) %>% 
  group_by(season,subsite_only,fluxtype) %>% 
  summarise(
    mean_value = mean(flux, na.rm = TRUE),
    sd_value = sd(flux, na.rm = TRUE),
    nobs = sum(!is.na(flux)),
    se = sd_value / sqrt(nobs),
    ci95lower = mean_value - qt(0.975, df = nobs - 1) * se,
    ci95upper = mean_value + qt(0.975, df = nobs - 1) * se
  )

write.csv(cu_sum, file = paste0(results_path,"summaries_curonian_ow/Curonian_seasonal_summary_openwater_S1-S4_CO2_CH4_umolm-2_s-1_AND_gwpCO2eq_mgm-2h-1.csv"),row.names = F)


ggplot(cu_sum, aes(x=season, y=mean_value))+
  geom_point()+
  geom_errorbar(aes(ymax=mean_value+se,ymin=mean_value-se),width = 0.5)+
  facet_grid(fluxtype~subsite_only,scales="free")





#Summary per subsite (all incubations together, not weighted by number of observations per season)
cu_sum_y<- fluxes_long %>% 
  filter(fluxtype%in%c("CO2_best.flux","CH4_best.flux","GWP_ch4","GWP_co2","GWP_total")) %>% 
  group_by(subsite_only,fluxtype) %>% 
  summarise(
    mean_value = mean(flux, na.rm = TRUE),
    sd_value = sd(flux, na.rm = TRUE),
    nobs = n(),
    se = sd_value / sqrt(nobs),
    ci95lower = mean_value - qt(0.975, df = nobs - 1) * se,
    ci95upper = mean_value + qt(0.975, df = nobs - 1) * se
  )

write.csv(cu_sum_y, file = paste0(results_path,"summaries_curonian_ow/Curonian_yearly_summary_openwater_S1-S4_CO2_CH4_umolm-2_s-1_AND_gwpCO2eq_mgm-2h-1.csv"),row.names = F)

ggplot(cu_sum_y, aes(x=subsite_only, y=mean_value))+
  geom_point()+
  geom_errorbar(aes(ymax=mean_value+se,ymin=mean_value-se),width = 0.3)+
  facet_grid(fluxtype~.,scales="free")


















#Summary all year using CO2_best.flux, CH4_diff_plus_bubling and GWP.
#Adjust weigths for equal influence of each season (irrespective of different Nobs per season). 

#CHECK PROCEDURE!!! Not sure of the approach taken here:
cu_sum_y_notsure<- fluxes_long %>%
  filter(fluxtype%in%c("CO2_best.flux","CH4_best.flux","GWP_ch4","GWP_co2","GWP_total")) %>% 
  group_by(season,subsite_only,fluxtype) %>% 
  summarise(
    nobs = n(),  # Number of observations per season
    mean_value = mean(flux, na.rm = TRUE),  # Mean of the variable
    sd_value = sd(flux, na.rm = TRUE),  # Standard deviation
    se_value = sd_value / sqrt(nobs),  # Standard error
    var_value = var(flux, na.rm = TRUE)  # Variance for weighted calculations
  ) %>%
  group_by(subsite_only, fluxtype) %>% 
  summarise(
    # Weighted mean calculation
    weighted_mean = sum(mean_value * nobs) / sum(nobs),
    
    # Weighted standard deviation calculation (variance first)
    weighted_variance = sum((nobs - 1) * var_value) / (sum(nobs) - length(unique(season))),
    weighted_sd = sqrt(weighted_variance),
    
    # Weighted standard error
    weighted_se = weighted_sd / sqrt(sum(nobs)),
    
    # 95% Confidence Interval
    ci95lower = weighted_mean - qt(0.975, df = sum(nobs) - length(unique(season))) * weighted_se,
    ci95upper = weighted_mean + qt(0.975, df = sum(nobs) - length(unique(season))) * weighted_se,
    
    #Total N
    nobs_total=sum(nobs)
    )


ggplot(cu_sum_y_notsure, aes(x=subsite_only, y=weighted_mean))+
  geom_point()+
  geom_errorbar(aes(ymax=ci95upper,ymin=ci95lower),width = 0.5)+
  facet_grid(fluxtype~.,scales="free")

