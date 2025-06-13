#Inspect & harmonize fieldsheet-veg metadata for database

#Author: Miguel Cabrera-Brufau
#Date: June 2025
#Project: Restore4cs


#Description----
#This script is used to harmonize the fieldsheet and vegetation information into useful tables to be used in the database and data-exploration. 

#Two csv files are created: 

  #One with harmonized details per UniqueID (including all incubations performed)
  #One with harmonized details per plotcode (per chamber deployment, discarding incubations performed "after vegetation cut")


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
require(dplyr)
require(purrr)
require(data.table)
require(tools)
library(hms)

repo_root <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
files.sources = list.files(path = paste0(repo_root,"/functions"), full.names = T)
for (f in files.sources){source(f)}


# ---- Directories ----
dropbox_root <- "C:/Users/Miguel/Dropbox/RESTORE4Cs - Fieldwork/Data" # You have to make sure this is pointing to the write folder on your local machine
#Path to Co2 and Ch4 auxfiles:
auxfile_path<- paste0(dropbox_root,"/GHG/Working data/Auxfiles_4aquaGHG/") 
#Path to aquaGHG best.flux results:
results_path <- paste0(dropbox_root,"/GHG/Processed data/computed_flux/")
#Path to save plots from exploration: 
plots_path <- paste0(dropbox_root,"/GHG/Processed data/computed_flux/plots")
#Path with fieldsheets:
fieldsheet_path <- paste0(dropbox_root,"/GHG/Fieldsheets")
#Path with vegetation data: 
vegetation_path<- paste0(dropbox_root, "/Vegetation/")


#----Import fluxes-----
#1. Import bestflux results
co2_best<- read.csv(paste0(results_path,"co2_bestflux.csv"))
ch4_best<- read.csv(paste0(results_path,"ch4_bestflux.csv"))



#Import meta-data------

  ##Fieldsheets----
#Read all GHGfieldsheets (using function from Camille)
fieldsheets<- list.files(path= fieldsheet_path, recursive = T, pattern = "GHG.xlsx",full.names = T)
field<- read_GHG_fieldsheets(fieldsheets)

  ##Vegetation data-----
#Import harmnonized vegetation
veg_data<- read.csv(paste0(vegetation_path,"RESTORE4Cs_finalvegetation_ABG_biomass_description.csv"))

#Join and keep only relevant info:
field_veg<- field %>% 
  rename(UniqueID=uniqID) %>% 
  separate(UniqueID, into = c("c1","c2","c3","c4","c5","c6","c7"), sep = "-", remove = F) %>% 
  mutate(plotcode=toupper(paste(c1,c2,c3,c4,sep = "-"))) %>% 
  merge.data.frame(veg_data %>% select(plotcode, vegetation_description, ABG_biomass_gpersquaremeter), by="plotcode", all=T) %>% 
  rename(sampling=subsite) %>% 
  mutate(season=substr(sampling, 1,2), 
         subsite=substr(sampling,4,8))
  
field_veg2<- field_veg %>% 
  #Remove not useful details
  select(-c(c1,c2,c3,c4,c5,c6,c7, 
            logger_floating_chamber,logger_transparent_chamber,person_sampling,
            chamber_type, chamber_height_cm, 
            unix_start,initial_co2, initial_ch4, 
            end_time, unix_stop,final_co2,final_ch4)) %>% 
  #Remove project-specific redundant info:
  select(-c(season, pilot_site, subsite))


#Exclude: After vegetation cut----
#get plotcodes with more than 1 strata (resulting from V--> B after veg removal)
inconsist_plotcodes<- field_veg2 %>% 
  select(plotcode, strata) %>% 
  distinct() %>%
  group_by(plotcode) %>% 
  mutate(nstrata=n()) %>% 
  filter(nstrata>1) %>% pull(plotcode)

#get details of plotcodes with issues
issue_plots<- field_veg2 %>% 
  filter(plotcode%in%inconsist_plotcodes)

#48 plots with bare after vegetation cut. Corresponding to 7 subsite visits (samplings): DA-A1 seasons 1-3 (not in S4), and some other subsite visits (not-consistent) 
issue_plots %>% select(sampling) %>% distinct()

#DECISSION: fluxes after such manipulation are not representative of real conditions. It's of little interest to other researchers. We will not include BARE after Cut in the final database or use them for the paper. 

#Harmonize veg-cut comment  to be able to filter out these incubations: "after vegetation removal" 
issue_plots %>% filter(strata=="bare") %>% group_by(sampling,comments) %>% summarise(n_incubations=n())

after_veg_cut_IDs<- field_veg2 %>% filter(grepl(pattern = "after vegetation removal",x = comments)) %>% pull(UniqueID)

  
  #Harmonize comments ------

#we need to decide on a way to include the relevant observations, Harmonizing comments into general observations not related to vegetation: Check what would those be:

  #"Rising tide (x-xcm depth)/Receding tide. water_depth==avg(min-max cm): DONE
  #"Submerged vegetation" (no further especification): DONE
  #"water_depth is lower estimate": DONE
  #"Tidal pool": DONE
  #"Microbial/algae mat visible": check for "green" "mat" "phytobentos", "microbial" "colour" "color" + pictures. Keywords and fieldsheets reviewed, DONE. 
  #dry/wet sediment | dry/wet soil | dry/wet sand: water_depth>0 + comments + pictures. Keywords and fieldsheets reviewed. Missing harmonization to unique obs: Dry/Wet soil/sediment/sand (single word for all? NO),DONE 
  
  
#Tide harmonization: 
  #check that water_depth is average of min/max
  #Check High-tide vs receding tide vs rising tide --> Only rising tide vs receding tide
  
tides<-field_veg2 %>% 
    # filter(grepl("tide", comments)) %>% 
    select(UniqueID,water_depth, comments) %>% 
    mutate(tide_behaviour=str_extract(pattern = "Rising-tide|Receding-tide", comments),
           depth_range=str_extract(pattern="([0-9]{1,2}-[0-9]{1,2}cm depth)", comments),
           tide_behaviour=if_else(is.na(tide_behaviour),"",tide_behaviour),
           depth_range=if_else(is.na(depth_range),"",depth_range),
           tide_presence=grepl("tide",comments),
           harmonized_tide_obs=paste(tide_behaviour,depth_range))
  
#Any non-tide with info captured in tide-behaviour or depth range?
  tides %>% filter(!tide_presence) %>% filter(!is.na(tide_behaviour)) %>% select(harmonized_tide_obs) %>% distinct()
  tides %>% filter(!tide_presence) %>% filter(!is.na(depth_range))%>% select(harmonized_tide_obs) %>% distinct()

  #Check unique combinations of tide_obs
  tides %>% filter(tide_presence) %>% select(harmonized_tide_obs) %>% distinct()  

  
#Submerged vegetation obs harmonization: 
sub.veg<- field_veg2 %>% 
  filter(strata=="open water") %>% 
  filter(!is.na(comments)) %>% 
  select(UniqueID, comments) %>% 
  mutate(harmonized_subveg_obs=if_else(grepl("submerged",comments, ignore.case = T), "Submerged vegetation",""))

#Check, any non-open_water with "submerged vegetation"?
field_veg2 %>% 
  mutate(harmonized_subveg_obs=if_else(grepl("submerged",comments, ignore.case = T), "Submerged vegetation","")) %>% 
           filter(strata!="open water") %>% 
  filter(harmonized_subveg_obs!="")
         

#Water depth obs  harmonization
larger_wd<- field_veg2 %>% 
  filter(!is.na(comments)) %>% 
  select(UniqueID, water_depth, comments) %>% 
  mutate(lager_water_depth_obs=if_else(str_detect(pattern = ">", comments),"water_depth is lower estimate", ""))

#Checked all comments containing ">" refer to water_depth greater than measure capability.


tidalpools<- field_veg2 %>%
  filter(!is.na(comments)) %>% 
  filter(water_depth>0) %>% 
  select(UniqueID, water_depth, comments) %>% 
  mutate(pool_obs=str_extract(pattern="Tidal-pool", comments))

#Check DU open waters== tidal pools? Only A1 and P1, based on GPS and pictures
field_veg2 %>% 
  filter(grepl("du",UniqueID)) %>% 
  filter(strata=="open water") %>% 
  select(UniqueID, comments)


#BARE obs harmonization 
#wet/"dry sediment"
#"Vissible biofilm"

bare_obs<- field_veg2 %>% 
  select(UniqueID,strata,water_depth,comments) %>% 
  mutate(moisture=str_extract(pattern = "wet soil|dry soil|dry sediment|wet sediment|wet sand|dry sand", comments),
         obs_moisture=if_else(is.na(moisture),"",moisture),
         biofilm=str_extract(pattern = "Biofilm|biofilm|phyto|microbial",comments),
         
         obs_biofilm=if_else(is.na(biofilm),"","microbial/algae mat visible"))

#Some biofilm/moisture comments in non-bare incubations: OK
bare_obs %>% 
  filter((!is.na(moisture)|(obs_biofilm!=""))&strata!="bare")


bare_obs %>% select(obs_moisture) %>% distinct()
bare_obs %>% select(obs_biofilm) %>% distinct()
bare_obs %>% filter(is.na(obs_biofilm))


#Join harmonized obs----
#Produce final meta-data: per UniqueID (with all incubations) per plotcode (without after veg-cut incubations).

##Details per UniqueID-----

#combine all harmonized obs,
meta_harm_uniqueID<-field_veg2 %>% 
  #Tide obs:
  mutate(tide_behaviour=tolower(str_extract(pattern = "Rising-tide|Receding-tide", comments)),
         depth_range=str_extract(pattern="([0-9]{1,2}-[0-9]{1,2}cm depth)", comments),
         obs_tide_behaviour=if_else(is.na(tide_behaviour),"",tide_behaviour),
         obs_depth_range=if_else(is.na(depth_range),"",depth_range),
         tide_presence=grepl("tide",comments),
         obs_tide=paste(tide_behaviour,depth_range)) %>% 
  select(-c(tide_behaviour,depth_range, tide_presence, obs_tide)) %>% 
  #Lower-estimate water depth obs: 
  mutate(obs_lager_depth=if_else(grepl(pattern = ">", comments),"water_depth is lower estimate", "")) %>% 
#Submerged vegetation obs harmonization: 
  mutate(obs_subveg=if_else(grepl("submerged",comments, ignore.case = T), "submerged vegetation","")) %>% 
  #Tidal-pools obs:
  mutate(obs_tidalpool=if_else(grepl("Tidal-pool",comments, ignore.case = T), "tidal-pool","")) %>% 
  #Wet sediment & pressence of microbial mat obs: 
    mutate(obs_moisture=str_extract(pattern = "wet soil|dry soil|dry sediment|wet sediment|wet sand|dry sand", comments),
           obs_moisture=if_else(is.na(obs_moisture),"",obs_moisture),
           obs_biofilm=if_else(grepl(pattern = "Biofilm|biofilm|phyto|microbial",comments),"microbial/algae mat visible","")) %>% 
  #COMBINE harmonized obs into a single string: either obs_combined or ""
  rowwise() %>% 
  mutate(obs_combined = paste(c_across(c(obs_tide_behaviour,obs_depth_range,obs_tidalpool, obs_lager_depth, obs_subveg,obs_moisture, obs_biofilm))[c_across(c(obs_tide_behaviour,obs_depth_range,obs_tidalpool, obs_lager_depth, obs_subveg,obs_moisture, obs_biofilm)) != ""], collapse = "; ")) %>%
  ungroup() %>% 
  mutate(obs_combined = sub("^; ", "", obs_combined)) %>% 
  #Remove individual obs: 
  select(-c(obs_tide_behaviour,obs_depth_range,obs_tidalpool, obs_lager_depth, obs_subveg,obs_moisture, obs_biofilm)) %>% 
  #ADAPT vegetation data: for non-vegetated plots: AGB=0, aboveground_vegetation_description="none"
  mutate(ABG_veg_description=if_else(strata=="vegetated",vegetation_description,"no aboveground vegetation"),
         ABG_veg_biomass_gpersquaremeter=if_else(ABG_veg_description=="no aboveground vegetation",0,ABG_biomass_gpersquaremeter)) %>% 
  select(-c(vegetation_description,ABG_biomass_gpersquaremeter))
  
#check that no grepl is assigning obs to non-appropriate comments. DONE
combine_to_review<- meta_harm_uniqueID %>% 
  filter(obs_combined!="") %>% 
  select(UniqueID,obs_combined, comments)


#FINAL HARMONIZED sampling details per UniqueID (all incubations)
uniqueIDs_harmonized_field_obs<- meta_harm_uniqueID %>% 
  #Rename plot_num and field_observations
  rename(plot_num=plot_id, field_observations=obs_combined) %>% 
  select(UniqueID, plotcode, sampling, plot_num, date, start_time, latitude,longitude, 
         gas_analyzer, strata, water_depth, field_observations, ABG_veg_description, ABG_veg_biomass_gpersquaremeter,transparent_dark)



##Details per plotcode-----

#FIRST: CHECK if there is any plotcode with inconsistent information (non-matching details for transparent and dark)
inconsistentplots<- meta_harm_uniqueID %>% 
  #remove incubations "after vegetation cut"  (they will be different strata and vegetation_data)
  filter(!UniqueID%in%after_veg_cut_IDs) %>% 
  #count, for each plotcode, how many distinct values exist for each variable (using variables that should have the same value for all incubations of the same plotcode
  group_by(plotcode) %>% 
  transmute(n_rows=n(),
         n_date=n_distinct(date),
         n_gas_analyzer=n_distinct(gas_analyzer),
         n_latitude=n_distinct(latitude),
         n_longitude=n_distinct(longitude),
         n_strata=n_distinct(strata),
         n_water_depth=n_distinct(water_depth),
         n_obs_combined=n_distinct(obs_combined),
         n_ABG_veg_description=n_distinct(ABG_veg_description),
         n_ABG_veg_biomass_gpersquaremeter=n_distinct(ABG_veg_biomass_gpersquaremeter)
  ) %>% 
  pivot_longer(cols = -c(plotcode,n_rows),names_to = "var", values_to = "count_obs") %>% 
  filter(n_rows>1) %>% #filter plotcodes with more than 1 incubation
  filter(count_obs>1) #any variable with more than 1 distinct value per plotcode

#NO inconsistent plotcode (appart from after-veg-cut incubations)
head(inconsistentplots)


#Extra check: is chamber height always the same for each plotcode?
ok_difheights<- c("S4-DA-P2-1","S4-DA-P2-11","S4-DA-P2-4","S4-DA-R1-2","S4-DA-R1-3","S4-DA-R1-4")
#In S4-DA no ring was used and chamber height is not the same for transparent and dark incubation: appropriate. 
field_veg %>% 
  group_by(plotcode) %>% 
  summarise(heights=n_distinct(chamber_height_cm)) %>% 
  filter(heights>1) %>% 
  filter(!plotcode%in%ok_difheights)




#Plotcode observations: 
#COMBINE harmonized observations and details into a single row per plotcode (assigned to transparent and dark UniqueIDs), one row per chamber deployment

plotcode_harmonized_field_obs<- meta_harm_uniqueID %>% 
  #remove incubations "after vegetation cut" 
  filter(!UniqueID%in%after_veg_cut_IDs) %>% 
  #add plot_start_time (first start_time per plotcode)
  mutate(start_time_hms = as_hms(start_time)) %>%
  group_by(plotcode) %>%
  mutate(plot_start_time = as.character(start_time[which.min(start_time_hms)])) %>%
  ungroup() %>%
  select(-c(start_time_hms,start_time,comments)) %>% 
  #Combine transparent_dark with UniqueID: migrate UniqueID identifier to transparent_UniqueID or dark_UniqueID
  pivot_wider(names_from = transparent_dark, values_from = UniqueID) %>% 
  rename(plot_num=plot_id, dark_UniqueID=dark, transparent_UniqueID=transparent, field_observations=obs_combined) %>% 
  #Reorder variables
  select(plotcode, sampling, plot_num, date, plot_start_time, latitude, longitude, gas_analyzer, strata, water_depth, field_observations, ABG_veg_description, ABG_veg_biomass_gpersquaremeter,dark_UniqueID,transparent_UniqueID) %>% 
  #arrange obs
  arrange(sampling, plot_num)
  


#Save fielddetails:
  #Per UniqueID
write.csv(x = uniqueIDs_harmonized_field_obs,file = paste0(results_path,"UniqueID_harmonized_field_obs.csv"),row.names = F)

  #Per plotcode
write.csv(x = plotcode_harmonized_field_obs,file = paste0(results_path,"Plotcode_harmonized_field_obs.csv"),row.names = F)







#______________________####
#EXTRA: Other decissions-----

##doubt channel VA-A1-----
#DOUBTs subsite delimitations: 
#VA-A1--> surrounding channel to be considered?

field_veg2 %>% 
  filter(grepl("VA-A1", sampling)) %>% 
  filter(strata=="open water") %>% 
  select(sampling, comments)

#Channel sampled in all seasons, it is the only water-incubation for S4. 
#Without channel, S1, S2 and S3, have >4 open water incubations. 

#DECISSION: fluxes from peripheral channel are not representative of subsite conditions, temporary water-bodies in seasons with water presence are well represented without including the channel samples. We will exclude these from paper (and database)

peripheral_channel_IDs<- field_veg2 %>% filter(grepl("VA-A1", sampling)) %>% filter(strata=="open water") %>% filter(grepl("Peripheral channel", comments)) %>% pull(UniqueID)


##doubt: bare-transparent----
#here we check if filtering out bare-transparent incubations, keeping only dark as per protocol would negatively impact the results making us lose representativity and statistical power.  

#Check subsites where bare is inconsistently sampled (bare-transparent without bare-dark)
bare_transparent_IDs<- field_veg2 %>% 
  filter(strata=="bare") %>% 
  filter(transparent_dark=="transparent") %>% pull(UniqueID)

only_bare_transparent_IDs<- field_veg2 %>% 
  filter(strata=="bare") %>% 
  group_by(plotcode) %>% 
  mutate(inc_perplot=n()) %>% 
  filter(inc_perplot==1&transparent_dark=="transparent") %>% pull(UniqueID)


#29 samplings with at least 1 bare-transparent
samplingswith_bare_t<- field_veg2 %>% filter(UniqueID%in%bare_transparent_IDs) %>% select(sampling) %>% distinct() %>% pull(sampling)

#Check subsite samplings that would lose bare-representativity (bare_n>=3) if we excluded bare-transparent incubations: 

subsite_bare_repre<- field_veg2 %>% 
  filter(sampling %in%samplingswith_bare_t) %>% 
  filter(strata=="bare") %>% 
  group_by(sampling) %>% 
  summarise(bare_plots=n_distinct(plotcode),
            bare_trans=sum(transparent_dark=="transparent"),
            bare_dark=sum(transparent_dark=="dark"),
            bare_total=bare_trans+bare_dark) %>%
  #Exclude RI Restored, no bare in actual subsite
  filter(!grepl("RI-R",sampling)) %>% 
  mutate(obs=case_when(bare_plots==bare_dark~"Good, all plots have dark flux",
                       grepl("RI-A1",sampling)&bare_dark<6~"100% bare site, not enough bare-dark",
                       bare_dark<3~"ISSUE,not enough bare-dark",
                       bare_dark>=3~"Good, enough bare-dark")) %>% 
  arrange(obs, bare_dark)

#GOOD: many samplings with bare-dark in every bare plot (transparent is extra incub):
good_samplings<-c("S1-VA-P2", #OK, 9 bare plots with dark (3 with tranparent)
                  "S1-RI-A1", #OK, 13 bare plots with dark (7 with transparent)
                  "S2-CU-R2", #Irrelevant, bare is not representative anyway
                  "S2-DA-A1", #OK, 7 bare plots with dark (6 with transparent)
                  "S2-DU-A1", #OK, 4 bare plots with dark (4 with transparent)
                  "S2-DU-A2", #OK, 3 bare plots with dark (3 with transparent)
                  "S2-DU-P2", #OK, 3 bare plots with dark (3 with transparent)
                  "S2-VA-R2", #OK, 6 bare plots with dark (1 with transparent)
                  "S3-DU-A1", #OK, 3 bare plots with dark (3 with transparent)
                  "S3-DU-A2", #OK, 3 bare plots with dark (3 with transparent)
                  "S3-RI-A1", #OK, 9 bare plots with dark (9 with transparent)
                  "S3-VA-P2", #OK, 3 bare plots with dark (3 with transparent)
                  "S3-VA-R2", #OK, 5 bare plots with dark (1 with transparent)
                  "S4-DU-A1", #OK, 3 bare plots with dark (3 with transparent)
                  "S4-DU-A2", #OK, 3 bare plots with dark (2 with transparent)
                  "S4-DU-P2", #OK, 3 bare plots with dark (3 with transparent)
                  "S4-VA-A2" #OK, 6 bare plots with dark (1 with transparent)
)

#GOOD enough: 
good_enough_samplings<- c("S3-DU-P2", #2 bare-plots with dark, (2 with transparent), 3rd bare plot lost (field complications)
                          "S1-DA-R2" #3 bare-transparent only BUT NON-issue: bare is not representative of subsite anyways (bare strata only considered in S1 season)
)
subsite_bare_repre %>% 
  filter(!sampling%in%c(good_samplings,good_enough_samplings))


#ISSUES: 

issue_samplings<- c("S2-RI-A1", #S2-RI-A1: 100% bare site, 13 plots (13 with transparent,BUT only first 3 also with dark, logger-confirmed) is it enough representativity to only keep dark? In seasons S1 and S3 also a lot of transparent incubations (extra), in S4 only bare-dark incubations. Most likely transparent after first 3 plots is decision due to lack of differences in flux. CHECK CO2&CH4 fluxes for all incubations. 
                    "S1-DU-A1",#S1-DU-A1: 4 bare-plots, 3 only transparent + 1 transparent&Dark
                    "S1-DU-A2",#S1-DU-A2: 3 bare-plots, 3 only transparent
                    "S1-DU-P2",#S1-DU-P2: 3 bare-plots, 2 only transparent + 1 transparent&Dark
                    #CHECK ALL fluxes transparent-dark S1-DU together
                    "S2-CA-P2" #S2-CA-P2: 9 bare-plots, 7 only transparent + 2 transparent&Dark, with so many bare, most likely decision to do transparent incubations is due to lack of difference in field.
) 


#Explore to decide: check dependence of bare-fluxes on light condition (graphically) for Co2 and Ch4 separately. 

bare_lightdark<- field_veg2 %>% 
  #Filter for incubations bare & samplings with transparent_dark
  filter(strata=="bare") %>% 
  # filter(sampling%in%samplingswith_bare_t) %>% 
  merge.data.frame(co2_best %>% select(UniqueID,best.model,LM.flux,HM.flux) %>% rename(co2_model=best.model,co2_LM.flux=LM.flux, co2_HM.flux=HM.flux), all.x = T) %>% 
  merge.data.frame(ch4_best %>% select(UniqueID,best.model,LM.flux,HM.flux,total.flux) %>% rename(ch4_model=best.model,ch4_LM.flux=LM.flux, ch4_HM.flux=HM.flux, ch4_total.flux=total.flux), all.x = T) %>% 
  #select best flux
  mutate(co2_bestflux=case_when(co2_model=="LM"~co2_LM.flux,
                                co2_model=="HM"~co2_HM.flux,
                                TRUE~NA_real_),
         ch4_bestflux=case_when(ch4_model=="LM"~ch4_LM.flux,
                                ch4_model=="HM"~ch4_HM.flux,
                                ch4_model=="total.flux"~ch4_total.flux,
                                TRUE~NA_real_)) %>% 
  #remove extra cols:
  select(-c(co2_model, co2_LM.flux,co2_HM.flux,
            ch4_model, ch4_LM.flux,ch4_HM.flux, ch4_total.flux)) %>% 
  #filter out incubations without any useful flux
  filter(!(is.na(co2_bestflux*ch4_bestflux)))


#FIRST ISSUE: S2-RI-A1: only 3 bare-dark for 100% bare subsite
#For CO2, only dark should be used. For CH4, decision not critical.
#Comparison transparent/dark CO2 
bare_lightdark %>% 
  filter(grepl("RI-A",sampling)) %>%  
  ggplot(aes(x=factor(plot_id), y=co2_bestflux, col=transparent_dark))+
  geom_point()+
  facet_wrap(~sampling)
#S2-RI-A1: strong & consistent difference bare-transparent vs bare-dark for CO2. Use only dark despite of very low n (still net negative flux for most samplings)


#Comparison transparent/dark CH4 
bare_lightdark %>% 
  filter(grepl("RI-A",sampling)) %>%  
  ggplot(aes(x=factor(plot_id), y=ch4_bestflux, col=transparent_dark))+
  geom_point()+
  facet_wrap(~sampling)
#S2-RI-A1: Consistent but not too strong effect due to light. Decision to include or exclude bare-transparent not critical. 



#SECOND ISSUE: #S2-CA-P2: 60% bare subsite, 9 bare-plots, only 2 with dark (+transparent), 7 just transparent
#No detectable effect of light on fluxes (high variability in general, no-light dependent)
#CO2, similar strata-average flux, could work with n=2 (only dark)
#CH4, similar  strata-average flux, could work with n=2 (only dark)
#We should use all bare fluxes, no matter their light condition.

#Comparison transparent/dark CO2 
bare_lightdark %>% 
  filter(grepl("CA-P2",sampling)) %>%  
  ggplot(aes(x=factor(plot_id), y=co2_bestflux, col=transparent_dark))+
  geom_point()+
  geom_boxplot(aes(x=7))+
  facet_wrap(~sampling)
#S2-CA-P2: no apparent effect of light on bare-CO2 fluxes. 3 last plots with near cero emissions (wet sediment, no light-dependent)

#Comparison transparent/dark CH4 
bare_lightdark %>% 
  filter(grepl("CA-P2",sampling)) %>%  
  ggplot(aes(x=factor(plot_id), y=ch4_bestflux, col=transparent_dark))+
  geom_point()+
  geom_boxplot(aes(x=7))+
  facet_wrap(~sampling)
#S2-CA-P2: no apparent effect of light on bare-CH4.  3 last plots with very-high fluxes  (wet sediment, no light-dependent). 



#ISSUES in DUTCH delta: 
# "S1-DU-A1",#S1-DU-A1: 4 bare-plots, 3 only transparent + 1 transparent&Dark
# "S1-DU-A2",#S1-DU-A2: 3 bare-plots, 3 only transparent
# "S1-DU-P2",#S1-DU-P2: 3 bare-plots, 2 only transparent + 1 transparent&Dark

#WHAT TO DO: 1 subsite losses all bare-representativity if we only use dark incubations. 

bare_lightdark %>% 
  filter(grepl("DU-", sampling))%>%  
  ggplot(aes(x=factor(plot_id), y=co2_bestflux, col=transparent_dark))+
  geom_point()+
  facet_wrap(~sampling)
#S1-DU-A1: similar values for bare-light and bare-transparent, no obvious difference
#S1-DU-A2: Only bare transparent nothing else that could be used. 1 emits, 2 drawdown
#s1-DU-P2: 1 bare-transparent with strong uptake, 2 transparent and 1 dark near cero flux with similar values. 

bare_lightdark %>% 
  filter(grepl("DU-", sampling))%>%  
  # filter(!grepl("S4-DU-P2", sampling)) %>% 
  ggplot(aes(x=factor(plot_id), y=ch4_bestflux, col=transparent_dark))+
  geom_point()+
  facet_wrap(~sampling)
#No differences light-dark for CH4

#DECISSION BARE-TRANSPARENT: only 1 subsite completely loses bare-strata representation, for others losing some bare strata representativity remaining fluxes are still within regular distribution. We will only use bare-dark for the paper. We can still include bare-transparent in database. 

