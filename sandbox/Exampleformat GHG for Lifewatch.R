
# ---
# Authors: Miguel Cabrera
# Project: "RESTORE4Cs"
# date: "Dec 2023"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
# This script re-formats the fluxes calculated in raw2flux script for the preeliminary data structure to send to LifeWatch Italy

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

repo_root <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
files.sources = list.files(path = paste0(repo_root,"/functions"), full.names = T)
for (f in files.sources){source(f)}


# ---- Directories ----
dropbox_root <- "C:/Users/Miguel/Dropbox/RESTORE4Cs - Fieldwork/Data" # You have to make sure this is pointing to the write folder on your local machine
# datapath <- paste0(dropbox_root,"/GHG/RAW data")
fieldsheetpath <- paste0(dropbox_root,"/GHG/Fieldsheets")
# loggerspath <- paste0(datapath,"/RAW Data Logger")
# RData_path <- paste0(dropbox_root,"/GHG/Processed data/RData/")
results_path <- paste0(dropbox_root,"/GHG/Processed data/computed_flux/")
vegetation_path<- paste0(dropbox_root,"/Vegetation/")
# plots_path <- paste0(dropbox_root,"/GHG/Processed data/plots_all_incubations/")



# Final variables and units: 
#season: summer,spring,autum, winter
#case_pilot: CA,DA,CU,RI,VA,DU
#subsite: casepilot_A1/A2/P1/P2/R1/R2
#date:
#plot_latitude
#plot_longitude
#plot_waterdepth: in cm
#plot_strata: open water, bare, vegetated land, vegetated water
#plot_vegetation_biomass
#plot vegetation description
#plot_uniqueID: subsite_plotnumber_o/b/vl/vw
#plot_starttime: first start.datetime of plot
#transparent_light_intensity: lux of transparent incubation
#dark_light_intensity: lux of dark incubation
#transparent_temperature: ºC of transparent
#dark_temperature: ºC of dark
#transparent_co2flux: umol/m2/s
#dark_co2flux: umol/m2/s
#transparent_diff_ch4flux: nmol/m2/s
#dark_diff_ch4flux: nmol/m2/s
#transparent_ebull_ch4flux: nmol/m2/s
#dark_ebull_ch4flux: nmol/m2/s
#transparent_n2oflux: nmol/m2/s
#dark_n2oflux: nmol/m2/s


# ---- Import fluxes per season ----

#Read all computed fluxes per season
season_csvs<- list.files(path = results_path, pattern = "^S")
dat<- data.frame()
for (i in season_csvs){
  s<- read.csv(paste0(results_path,i))
  
  dat<- rbind(dat,s)
}
rm(s, i)


str(dat)


#Keep only fluxes data, start.time and UniqueID.
#For CO2 fluxes, take only CO2_best.flux
#For CH4 fluxes, take CH4_diffusive_flux and CH4_ebullitive_flux
#We will have to adapt the raw2flux script so that CH4 has always data for diffusion and ebullition (if no ebullition, ebullition==0)

dat2<- dat %>% 
  select(UniqueID,
         start.time,
         CO2_best.flux,
         #CO2_LM.flux,CO2_HM.flux,CO2_best.flux,CO2_best.model,CO2_quality.check,
         #CH4_LM.flux,CH4_HM.flux,CH4_best.flux,CH4_best.model, CH4_quality.check,
         CH4_diffusive_flux,CH4_ebullitive_flux)



#Import fieldsheet meta-data#----

#Read all GHGfieldsheets (using function from Camille)
fieldsheets<- list.files(path= fieldsheetpath, recursive = T, pattern = "GHG.xlsx",full.names = T)
field<- read_GHG_fieldsheets(fieldsheets)

#Drop unwanted variables and rename key UniqueID:
field2<- field %>% 
  select(-c(person_sampling,gas_analyzer,logger_floating_chamber,logger_transparent_chamber,start_time,initial_co2,initial_ch4,end_time,final_ch4,final_co2,chamber_height_cm, unix_start,unix_stop)) %>% 
  rename(UniqueID=uniqID)


# ---- Format and join ----
#Key for joining flux data and fieldsheet: 
summary(dat$UniqueID%in%field$uniqID)




#Check for duplicate incubations in fieldsheets (incubations performed 2 times with the same identity except for startime and uniqueID:  
duplicate_incubations<-field2 %>% 
  dplyr::summarise(n = dplyr::n(), .by = c(pilot_site, subsite, date, plot_id, chamber_type, strata, longitude, latitude, water_depth,
                                           transparent_dark)) |>dplyr::filter(n > 1L) %>% 
  merge.data.frame(field2, by=c("pilot_site", "subsite", "date", "plot_id", "chamber_type", "strata", "longitude", "latitude", "water_depth","transparent_dark"), all = F) %>% pull(UniqueID)


#merge fluxes with fieldsheets dropping the unmatched observations, problematic values and duplicated incubations
a<- field2 %>% merge.data.frame(dat2, by="UniqueID",all = F) %>%
  filter(!is.na(transparent_dark)) %>% 
  filter(!UniqueID%in%duplicate_incubations) %>% 
  select(-UniqueID) %>% 
  pivot_wider(names_from = transparent_dark,values_from = c(start.time, CO2_best.flux, CH4_diffusive_flux,CH4_ebullitive_flux,comments)) %>% 
  mutate(plot_startime=min(start.time_dark, start.time_transparent, na.rm=T),
         plotcode=paste(subsite, plot_id,sep = "-"),
         campaign=str_extract(subsite, pattern="S[0-9]{1}"),
         season=case_when(campaign=="S1"~"fall",
                          campaign=="S2"~"winter",
                          campaign=="S3"~"spring",
                          campaign=="S4"~"summer"),
         subsite=str_extract(subsite, pattern="[A-Z]{2}-[A-Z]{1}[0-9]{1}"),
         sampling=paste(campaign, subsite,sep="-"),
         strata=case_when(strata=="vegetated"&water_depth>=10~"vegetated water",#Re-classify strata using water depth (vegetated water when water_depth>=10cm, otherwise vegetated land)
                          strata=="vegetated"&water_depth<10~"vegetated land",
                          TRUE~strata))





# ---- Data loggers ----
#We need Temp and Light intensity for dark and transparent conditions. 
#extract from data loggers files using info from field: logger_floating_chamber, logger_transparent_chamber, unixstart, unixstop, uniqID
#In raw2flux, temperature extraction already implemented(CHECK and adapt for Light intensity)



# ---- Vegetation ID & DW ----

#####REVISAR######
#Revisar codigo de vegetacion, quitar cosas innecesarias y  commentar codigo


#Vegetation description and biomassDW: 
#1. Format vegetation files to Identify based on Season-site-subsite-plot information. Final variables: samplecode, DW (g), comments_veg.DONE

#Import pattern: RESTORE4Cs_vegetation.........xlsx (rbind all). Take headers and drop row 1 with col-description
#Headers: c("person_measuring","season","pilot_site", "subsite", "plot_id", "dry weigth", "comments")
#season = S1/S2/S3/S4
#pilot_site= CA/CU/DA/DU/VA/RI
#subsite= Preserved 1/2, Altered 1/2, Restored 1/2
#plot_id= integer
#dry_weight=numeric
#comment= text

#Some issues identified in vegetation excels (duplicate labels, no labels, ...)
#Cross-referenced with GHG fieldsheets to correct obvious mistakes.

#List all vegetation files
veg_files<- list.files(path = vegetation_path, pattern = "^RESTORE4Cs_vegetation_",full.names = T,recursive = T)
#Read all computed fluxes per season
veg<- read_vegetation_fieldsheets(veg_files)

#Duplicated plant DWs
veg %>% 
  dplyr::summarise(n = dplyr::n(), .by = c(plotcode)) |>dplyr::filter(n > 1L) %>% 
  merge.data.frame(veg, by=c("plotcode"), all = F)


#2. Extract comments_dark and comments_transparent from fieldsheets (vegetated plots)
vegetated_plots<- field2 %>% 
  filter(!UniqueID%in%c("s1-va-a1-12-v-t-14:58","s1-va-a1-12-v-d-15:06")) %>% #remove repeated incubation
  filter(strata=="vegetated") %>%#Select only vegetated plots
  mutate(plotcode=paste(subsite, plot_id,sep="-")) %>% 
  select(plotcode, strata,transparent_dark, comments) %>% 
  pivot_wider(values_from = comments, names_from = transparent_dark) %>% 
  distinct()

vegetated_plots_veg<- vegetated_plots %>% 
  merge.data.frame(veg, by=c("plotcode"),all = T)


#vegetation DW data without vegetated plotGHG: 14 unmatched vegetation samples
veg[!veg$plotcode%in%vegetated_plots$plotcode,]

#Vegetated plotsGHG without vegetationDW 
vegetated_plots[!vegetated_plots$plotcode%in%veg$plotcode,]



#3. Join and combine comments to get vegetation description (sp/genus/family/), calculate ABG_biomass using tube area (pi*(12.1cm radius)^2)
vegetation_dw_final<- vegetated_plots %>% 
  merge.data.frame(veg, by=c("plotcode"),all.x = T) %>% 
  mutate(ABG_biomass=dry_weight/(pi*12.1^2)) %>% 
  select(plotcode, strata, ABG_biomass, transparent, dark, comments) %>% 
  rename(comment1=transparent, comment2=dark, comment3=comments)




