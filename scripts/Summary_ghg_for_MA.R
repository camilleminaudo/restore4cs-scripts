
# ---
# Authors: Miguel Cabrera
# Project: "RESTORE4Cs"
# date: "June 2025"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

rm(list = ls()) # clear workspace
cat("/014") # clear console

# --- Description of this script
# This script makes a summary of in-situ GHG CO2, CH4 and N2O fluxes to be used in Meta-Analysis. Main decisions for fluxes selection and temporal integration are recorded throughout the script. 

#Inputs: 

#1. best.flux csv files for CO2 and CH4 produced by  Selection criteria for aquaGHG best.flux.R script
#2. N2O fluxes for season 4: S4_restore4cs_N2O_arealflux_nmol_s-1_m-2.csv

#output: single csv with per-subsite global averages for each GHG (mean+sd+N) (with stamp of date created)
# Final units: netCO2: umol/m2/s , netCH4: nmol/m2/s, netN2O: nmol/m2/s


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


#Save today in order to paste daystamp to output csv-file:
dayoflastrun <- today()


# ---- Import fluxes ----

#Read all CO2 and CO2 computed fluxes (Selection of best.flux based on model fit, and ebullition detection for CH4)
co2_best<- read.csv(paste0(results_path,"co2_bestflux.csv"))
ch4_best<- read.csv(paste0(results_path,"ch4_bestflux.csv"))

n2o_best<- read.csv(paste0(dropbox_root,"/GHG/N2O_fluxes/","S4_restore4cs_N2O_arealflux_nmol_s-1_m-2.csv"))



#---Join & format----
#get for each GHG: UniqueID_notime, best.flux, best.flux.se
#Incubations with no best.model selection are assigned NA as flux (this are incubations marked to discard)

#CO2: no filtering, discard decision do not have a best.model
co2_format<- co2_best %>% 
  mutate(co2_flux=case_when(best.model=="LM"~LM.flux,
                            best.model=="HM"~HM.flux,
                            best.model=="total.flux"~total.flux,
                            TRUE~NA_real_),
         co2_se=case_when(best.model=="LM"~LM.flux.se,
                          best.model=="HM"~HM.flux.se,
                          best.model=="total.flux"~total.flux.se,
                          TRUE~NA_real_)) %>%
  separate(UniqueID, into = c("c1","c2","c3","c4","c5","c6","c7"),sep = "-",remove = F) %>% 
  mutate(UniqueID_notime=paste(c1,c2,c3,c4,c5,c6,sep = "-")) %>% 
  select(UniqueID_notime, co2_flux,co2_se)

#CH4:no filtering, discard decision do not have a best.model
ch4_format<- ch4_best %>% 
  mutate(ch4_flux=case_when(best.model=="LM"~LM.flux,
                            best.model=="HM"~HM.flux,
                            best.model=="total.flux"~total.flux,
                            TRUE~NA_real_), 
         ch4_se=case_when(best.model=="LM"~LM.flux.se,
                          best.model=="HM"~HM.flux.se,
                          best.model=="total.flux"~total.flux.se,
                          TRUE~NA_real_)) %>%
  separate(UniqueID, into = c("c1","c2","c3","c4","c5","c6","c7"),sep = "-",remove = F) %>% 
  mutate(UniqueID_notime=paste(c1,c2,c3,c4,c5,c6,sep = "-")) %>% 
  select(UniqueID_notime, ch4_flux,ch4_se)


#N2O: 
#Filtering of N2O fluxes should be based on the precission (relative_SE) of the concentrations leading to the flux, not on the precission/significance of the flux itself. 

#filter out N2O fluxes calculated based on imprecise (CV>5%) final concentrations, marked in N2Oflux_qualityflag column
n2o_format <- n2o_best %>% 
  mutate(n2o_flux=N2Oflux_nmol_per_second_per_m2,
         n2o_se=N2Oflux_se, 
         rel_se=n2o_se/n2o_flux) %>% 
  filter(N2Oflux_qualityflag=="Reliable flux") %>% 
  select(UniqueID_notime, n2o_flux, n2o_se)


#Join GHG fluxes
ghg_format<- co2_format %>% 
  merge.data.frame(ch4_format, by="UniqueID_notime", all = T) %>% 
  merge.data.frame(n2o_format, by="UniqueID_notime", all =T)
  




#----Import fieldsheets----

#Read all GHGfieldsheets (using function from Camille)
fieldsheets<- list.files(path= fieldsheetpath, recursive = T, pattern = "GHG.xlsx",full.names = T)
field<- read_GHG_fieldsheets(fieldsheets)

#Drop unwanted variables and rename key UniqueID:
field2<- field %>% 
  select(-c(person_sampling,gas_analyzer,logger_floating_chamber,logger_transparent_chamber,initial_co2,initial_ch4,end_time,final_ch4,final_co2,chamber_height_cm, unix_start,unix_stop)) %>% 
  rename(UniqueID=uniqID) %>% 
  separate(UniqueID, into = c("c1","c2","c3","c4","c5","c6","c7"),sep = "-",remove = F) %>% 
  mutate(UniqueID_notime=paste(c1,c2,c3,c4,c5,c6,sep = "-")) %>% 
  select(-c(c1,c2,c3,c4,c5,c6,c7))

#Check uniqueID_notime does not lose UNIQUEness
length(unique(field2$UniqueID))==length(unique(field2$UniqueID_notime))



# ---- Join fluxes and meta-data ----
#Key for joining flux data and fieldsheet: 
summary(ghg_format$UniqueID_notime%in%field2$UniqueID_notime)

#fluxes without fieldsheet details: arising from correction in fieldsheets
ghg_format[which(!ghg_format$UniqueID_notime%in%field2$UniqueID_notime),]

#fieldsheets without flux calculated (WHY?)
incub_noflux<-field2[which(!field2$UniqueID_notime%in%ghg_format$UniqueID_notime),]

#7 incubations without a flux (errors with gas analizers)
incub_noflux %>% select(UniqueID_notime, comments)


#Check for duplicate incubations in fieldsheets (incubations performed 2 times with the same identity except for startime and uniqueID:  
duplicate_incubations<-field2 %>% 
  dplyr::summarise(n = dplyr::n(), .by = c(pilot_site, subsite, date, plot_id, chamber_type, strata, longitude, latitude, water_depth,
                                           transparent_dark)) |>dplyr::filter(n > 1L) %>% 
  merge.data.frame(field2, by=c("pilot_site", "subsite", "date", "plot_id", "chamber_type", "strata", "longitude", "latitude", "water_depth","transparent_dark"), all = F) %>% pull(UniqueID)
duplicate_incubations


#Check for missing transparent_dark
field2 %>% filter(is.na(transparent_dark))


#merge fluxes with fieldsheets dropping the unmatched observations (no flux was calculated)

ghg_widelight<- field2 %>% merge.data.frame(ghg_format, by="UniqueID_notime",all.y = T) %>%
  filter(!is.na(transparent_dark)) %>% 
  filter(!UniqueID%in%duplicate_incubations) %>% 
  select(-c(UniqueID, UniqueID_notime,co2_se,ch4_se, n2o_se)) %>% 
  pivot_wider(names_from = transparent_dark,values_from = c(start_time, co2_flux, ch4_flux, n2o_flux, comments)) %>% 
  rowwise() %>% 
  mutate(plot_startime=min(start_time_dark, start_time_transparent, na.rm=T)) %>% 
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
  select(season, campaign,pilot_site, subsite, sampling, date, latitude, longitude, water_depth, strata, plotcode, plot_startime, 
         co2_flux_dark, co2_flux_transparent,
         ch4_flux_dark, ch4_flux_transparent,
         n2o_flux_dark, n2o_flux_transparent)




#Get day-night hours####

#Integrate day-night per strata

#Calculate daylight hours using latitude and date (hours of light)

# get median coordinates per subsite (including all seasons). Use of median avoids chambers with errors in coordinates afecting the calculation
samplings <- ghg_widelight %>% 
  select(subsite, latitude, longitude) %>% 
  group_by(subsite) %>% 
  summarise(subsite_latitude=median(latitude,na.rm=T),
            subsite_longitude=median(longitude,na.rm=T)) %>% 
  ungroup()


#Use dplyr and suncalc to calculate daylight hours for each subsite-date combination.
sampling_daylight <- samplings %>%
  merge.data.frame(y=ghg_widelight %>% select(subsite, date), by="subsite", all=T) %>% 
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

sampling_daylight %>% filter(subsite=="RI-R2") %>% arrange(date)
#ISSUE: RI-R2 has been sampled in different days in S1 (november2023), S2 (january2024), S3 (april2024) and S4 (july2024)
#Not a problem, lighthours are almost the same. Joining with subsite and date as keys. 


#Test Integrate day-night fluxes####
#INTEGRATE to net_GHG exchange per plot (integrate dark and transparent with daylight_duration for vegetated plots, keep dark bare and dark OW fluxes as representative of whole day.

net_ghg_per_plot<- ghg_widelight %>% 
  mutate(date=as.Date(date)) %>% 
  merge.data.frame(sampling_daylight, by=c("subsite","date"), all = T) %>% 
  select(season, pilot_site, subsite, sampling, subsite_latitude, subsite_longitude, strata, date, daylight_duration, plotcode, latitude, longitude,  
         co2_flux_dark, co2_flux_transparent,
         ch4_flux_dark, ch4_flux_transparent,
         n2o_flux_dark, n2o_flux_transparent) %>% 
  #Integrate transparent and dark with daylight_duration for vegetated fluxes(fluxlight*lighthours+fluxdark*darkhours)/24: same units for net flux (umol/m2/s), but integrated for the whole day. 
  #CO2:
  mutate(net_CO2=case_when(grepl("vegetated",strata)~((co2_flux_transparent*daylight_duration)+(co2_flux_dark*(24-daylight_duration)))/24,
                           #For Bare plots, take dark flux as representative for the whole day
                           grepl("bare",strata)~co2_flux_dark,
                           #For open water plots, use dark flux for whole day
                           grepl("open",strata)~co2_flux_dark),
         
         #CH4: same approach than for CO2 (scaled for vegetated plots, dark fluxes for the rest)
         net_CH4=case_when(grepl("vegetated",strata)~((ch4_flux_transparent*daylight_duration)+(ch4_flux_dark*(24-daylight_duration)))/24,
                           #For Bare plots, use dark flux for whole day
                           grepl("bare",strata)~ch4_flux_dark,
                           #For open water plots, use dark flux for whole day
                           grepl("open",strata)~ch4_flux_dark),
         #N2O: same approach than for CO2 (scaled for vegetated plots, dark fluxes for the rest)
         net_N2O=case_when(grepl("vegetated",strata)~((n2o_flux_transparent*daylight_duration)+(n2o_flux_dark*(24-daylight_duration)))/24,
                           #For Bare plots, use dark flux for whole day
                           grepl("bare",strata)~n2o_flux_dark,
                           #For open water plots, use dark flux for whole day
                           grepl("open",strata)~n2o_flux_dark),
  ) %>% 
  #Discard _transparent/_dark ghg fluxes
  select(-c(co2_flux_dark, co2_flux_transparent,
            ch4_flux_dark, ch4_flux_transparent,
            n2o_flux_dark, n2o_flux_transparent)) 


#How many plots end up with a valid net_flux after all filtering (%oftotal)
message(paste0("We have valid net_CO2 fluxes for ",round(sum(!is.na(net_ghg_per_plot$net_CO2))/dim(net_ghg_per_plot)[1]*100,2),"% of the chambers deployed, after filtering and integrating"))

message(paste0("We have valid net_CH4 fluxes for ",round(sum(!is.na(net_ghg_per_plot$net_CH4))/dim(net_ghg_per_plot)[1]*100,2),"% of the chambers deployed, after filtering and integrating"))

message(paste0("We have valid net_N2O fluxes for ",round(sum(!is.na(net_ghg_per_plot$net_N2O))/dim(net_ghg_per_plot)[1]*100,2),"% of the chambers deployed, after filtering and integrating"))


#How many bare fluxes have we discarded due to only having transparent flux?

#How many bare plots with at least 1 flux (transparent or dark)
n_fluxplotbare_co2<-dim(ghg_widelight %>% filter(strata=="bare") %>% 
  filter((is.na(co2_flux_dark)+is.na(co2_flux_transparent)<2)))[1]
n_fluxplotbare_ch4<-dim(ghg_widelight %>% filter(strata=="bare") %>% 
  filter((is.na(ch4_flux_dark)+is.na(ch4_flux_transparent)<2)))[1]
n_fluxplotbare_n2o<-dim(ghg_widelight %>% filter(strata=="bare") %>% 
  filter((is.na(n2o_flux_dark)+is.na(n2o_flux_transparent)<2)))[1]


n_netfluxplotbare_co2<- dim(net_ghg_per_plot %>% filter(strata=="bare") %>% filter(!is.na(net_CO2)))[1]
n_netfluxplotbare_ch4<- dim(net_ghg_per_plot %>% filter(strata=="bare") %>% filter(!is.na(net_CH4)))[1]
n_netfluxplotbare_n2o<- dim(net_ghg_per_plot %>% filter(strata=="bare") %>% filter(!is.na(net_N2O)))[1]

message(paste0("We have valid net_CO2 for ", n_netfluxplotbare_co2, " chambers out of ", n_fluxplotbare_co2, " bare chambers with at least 1 valid flux (transparent or dark): ", round(n_netfluxplotbare_co2/n_fluxplotbare_co2*100), "%"))

message(paste0("We have valid net_CH4 for ", n_netfluxplotbare_ch4, " chambers out of ", n_fluxplotbare_ch4, " bare chambers with at least 1 valid flux (transparent or dark): ", round(n_netfluxplotbare_ch4/n_fluxplotbare_ch4*100), "%"))

message(paste0("We have valid net_N2O for ", n_netfluxplotbare_n2o, " chambers out of ", n_fluxplotbare_n2o, " bare chambers with at least 1 valid flux (transparent or dark): ", round(n_netfluxplotbare_n2o/n_fluxplotbare_n2o*100), "%"))


#Which subsite samplings correspond to discarded net_co2 and net_ch4 fluxes?
ghg_widelight %>% filter(strata=="bare") %>% 
  select(plotcode, sampling,co2_flux_dark,co2_flux_transparent,ch4_flux_dark,ch4_flux_transparent) %>% 
  merge.data.frame(net_ghg_per_plot, by=c("plotcode","sampling")) %>% 
  select(plotcode, sampling,co2_flux_dark,co2_flux_transparent,ch4_flux_dark,ch4_flux_transparent, net_CO2, net_CH4) %>% 
  group_by(sampling) %>% 
  summarise(plots_anyflux_co2=sum((is.na(co2_flux_dark)+is.na(co2_flux_transparent)<2)),
            net_fluxes_co2=sum(!is.na(net_CO2)),
            plots_anyflux_ch4=sum((is.na(ch4_flux_dark)+is.na(ch4_flux_transparent)<2)),
            net_fluxes_ch4=sum(!is.na(net_CH4))) %>% 
  filter(plots_anyflux_co2!=net_fluxes_co2|plots_anyflux_ch4!=net_fluxes_ch4)
  
#IMPORTANT---- 
#taking only dark fluxes from bare strata would almost all fluxes from 2 samplings with bare as dominant strata (S2-CA-P2 and S2-RI-A1), additionally, 4 samplings lose the representation of bare strata. 

#Decision: we will use mean of transparent and dark fluxes for bare plots (this way we balance representativity and effect of light/dark). Not ideal, but we cannot discard 2 whole subsites. Both subsites have some incubations with corresponding light/dark, we assume that the decission to do only transparent bare incubations after that was based on expert criteria in the field. This ensures that we keep strata representativity (taking only 1 flux per plot) while avoiding the loss of data coming from non-sistematic samplings. 


#Integrate day-night (V2)-----

net_ghg_per_plot_v2<- ghg_widelight %>% 
  mutate(date=as.Date(date)) %>% 
  merge.data.frame(sampling_daylight, by=c("subsite","date"), all = T) %>% 
  select(season, pilot_site, subsite, sampling, subsite_latitude, subsite_longitude, strata, date, daylight_duration, plotcode, latitude, longitude,  
         co2_flux_dark, co2_flux_transparent,
         ch4_flux_dark, ch4_flux_transparent,
         n2o_flux_dark, n2o_flux_transparent) %>% 
  #Integrate transparent and dark with daylight_duration for vegetated fluxes(fluxlight*lighthours+fluxdark*darkhours)/24: same units for net flux (umol/m2/s), but integrated for the whole day. 
  #net_GHG comes from day/night integration for all strata
  mutate(net_CO2=((co2_flux_transparent*daylight_duration)+(co2_flux_dark*(24-daylight_duration)))/24,
         net_CH4=((ch4_flux_transparent*daylight_duration)+(ch4_flux_dark*(24-daylight_duration)))/24,
         net_N2O=((n2o_flux_transparent*daylight_duration)+(n2o_flux_dark*(24-daylight_duration)))/24) %>% 
  #mean_GHG comes from averaging any fluxes (transparent or dark) of each plot (row-wise)
  rowwise() %>% 
  mutate(mean_CO2=mean(c_across(c(co2_flux_transparent, co2_flux_dark)), na.rm=T),
         mean_CH4=mean(c_across(c(ch4_flux_transparent,ch4_flux_dark)), na.rm=T),
         mean_N2O=mean(c_across(c(n2o_flux_transparent,n2o_flux_dark)), na.rm=T)) %>% 
  #Override NAs in net_GHG fluxes for strata bare and open water (this way, we are using fluxes dark or transparent as representative of whole day for these strata when only 1 flux is available): this will substitute all open water fluxes and only bare fluxes that do not have matching light&dark fluxes
mutate(net_CO2=case_when(grepl("bare|open", strata)&is.na(net_CO2)~mean_CO2,
                         TRUE~net_CO2),
       net_CH4=case_when(grepl("bare|open", strata)&is.na(net_CH4)~mean_CH4,
                         TRUE~net_CH4),
       net_N2O=case_when(grepl("bare|open", strata)&is.na(net_N2O)~mean_N2O,
                         TRUE~net_N2O))%>% 
  #Discard _transparent/_dark ghg fluxes
  select(-c(co2_flux_dark, co2_flux_transparent,
            ch4_flux_dark, ch4_flux_transparent,
            n2o_flux_dark, n2o_flux_transparent,
            mean_CO2,mean_CH4,mean_N2O)) 


#Summary per subsite####

#Summarise day-integrated fluxes without strata or seasonal groupings, save in export_path.
#Keep median lat/long per subsite 

net_ghg_summary_persubsite<- net_ghg_per_plot_v2 %>% 
  select(-c(season, pilot_site, sampling, strata, daylight_duration, plotcode,latitude, longitude,date)) %>% 
  separate(subsite, into=c("pilot_site","subsite")) %>% 
  mutate(status=substr(subsite,1,1)) %>% 
  group_by(pilot_site, status, subsite ,subsite_latitude,subsite_longitude) %>% 
  summarise(avg_netCO2=mean(net_CO2, na.rm=T),
            sd_netCO2=sd(net_CO2, na.rm=T),
            n_netCO2=sum(!is.na(net_CO2)),
            avg_netCH4=mean(net_CH4, na.rm=T),
            sd_netCH4=sd(net_CH4, na.rm=T),
            n_netCH4=sum(!is.na(net_CH4)),
            avg_netN2O=mean(net_N2O, na.rm=T),
            sd_netN2O=sd(net_N2O, na.rm=T),
            n_netN2O=sum(!is.na(net_N2O)),
            .groups = "drop"
  )



# Summarise day-integrated fluxes using only the central 95% of each variable in each subsite, handling NAs. DO not filter this for N2O, with very few datapoints, this filtering is not appropriate (it removes at minimum 2 datapoints (highest & lowest)), many N2O samplings with very few 5-7 fluxes. 

net_ghg_summary_central95_persubsite<- net_ghg_per_plot_v2 %>%
  select(-c(season, pilot_site, sampling, strata, daylight_duration, plotcode,latitude, longitude,date)) %>% 
  separate(subsite, into=c("pilot_site","subsite")) %>% 
  mutate(status=substr(subsite,1,1)) %>% 
  group_by(pilot_site, status, subsite ,subsite_latitude,subsite_longitude) %>% 
  # For net_CO2
  summarise(
    net_CO2_lower = quantile(net_CO2, 0.025, na.rm = TRUE),
    net_CO2_upper = quantile(net_CO2, 0.975, na.rm = TRUE),
    avg_netCO2_95central = mean(net_CO2[net_CO2 >= net_CO2_lower & net_CO2 <= net_CO2_upper], na.rm = TRUE),
    sd_netCO2_95central = sd(net_CO2[net_CO2 >= net_CO2_lower & net_CO2 <= net_CO2_upper], na.rm = TRUE),
    n_netCO2_95central = sum(!is.na(net_CO2) & net_CO2 >= net_CO2_lower & net_CO2 <= net_CO2_upper),
    
    # For net_CH4
    net_CH4_lower = quantile(net_CH4, 0.025, na.rm = TRUE),
    net_CH4_upper = quantile(net_CH4, 0.975, na.rm = TRUE),
    avg_netCH4_95central = mean(net_CH4[net_CH4 >= net_CH4_lower & net_CH4 <= net_CH4_upper], na.rm = TRUE),
    sd_netCH4_95central = sd(net_CH4[net_CH4 >= net_CH4_lower & net_CH4 <= net_CH4_upper], na.rm = TRUE),
    n_netCH4_95central = sum(!is.na(net_CH4) & net_CH4 >= net_CH4_lower & net_CH4 <= net_CH4_upper),
    .groups = "drop"
  ) %>% 
  select(-c(net_CO2_lower,net_CO2_upper,net_CH4_lower,net_CH4_upper))



#Plot results: 
net_ghg_per_plot_toplot<- net_ghg_per_plot_v2 %>% 
  mutate(subsite_only=substr(subsite,4,5),
         status=substr(subsite_only,1,1))

net_ghg_per_plot_95<- net_ghg_per_plot_toplot %>% 
  group_by(pilot_site,subsite) %>% 
  mutate(upper_netCO2=quantile(net_CO2,0.975,na.rm=T),
         lower_netCO2=quantile(net_CO2,0.025,na.rm=T),
         upper_netCH4=quantile(net_CH4,0.975,na.rm=T),
         lower_netCH4=quantile(net_CH4,0.025,na.rm=T)) %>% 
  mutate(net_CO2=case_when(between(net_CO2,lower_netCO2,upper_netCO2)~net_CO2,
                           TRUE~NA_real_),
         net_CH4=case_when(between(net_CH4,lower_netCH4, upper_netCH4)~net_CH4,
                           TRUE~NA_real_))
  
  #net_CO2: in black the data outside the 95% central distribution
ggplot(net_ghg_per_plot_95, aes(x=subsite_only, y=net_CO2, col=status))+
  geom_boxplot()+
  geom_point(data=net_ghg_per_plot_toplot, col="black")+
  geom_point()+
  facet_grid(pilot_site~., scales="free")

ggplot(net_ghg_per_plot_95, aes(x=subsite_only, y=net_CO2, col=status))+
  geom_boxplot()+
  geom_point()+
  facet_grid(pilot_site~., scales="free")

  #net_CH4: in black the data outside the 95% central distribution
ggplot(net_ghg_per_plot_95, aes(x=subsite_only, y=net_CH4, col=status))+
  geom_boxplot()+
  geom_point(data=net_ghg_per_plot_toplot, col="black")+
  geom_point()+
  facet_grid(pilot_site~., scales="free")

ggplot(net_ghg_per_plot_95, aes(x=subsite_only, y=net_CH4, col=status))+
  geom_boxplot()+
  geom_point()+
  facet_grid(pilot_site~., scales="free")

  #net_N2O: in black the data outside the 95% central distribution
ggplot(net_ghg_per_plot_95, aes(x=subsite_only, y=net_N2O, col=status))+
  geom_boxplot()+
  geom_point(data=net_ghg_per_plot_toplot, col="black")+
  geom_point()+
  facet_grid(pilot_site~., scales="free")


#Save net_ghg per subsite_with explicit date of last update.
write.csv(net_ghg_summary_persubsite,file = paste0(export_path, "net_ghg_summary_persubsite_V",dayoflastrun, ".csv"),row.names = F)

#Save the central 95% summary
write.csv(net_ghg_summary_central95_persubsite,file = paste0(export_path, "net_ghg_summary_central95_persubsite_V",dayoflastrun, ".csv"),row.names = F)



#Check improvement on SD and n discarded for CO2 and CH4 between the two summaries: 
names(net_ghg_summary_persubsite)
net_ghg_summary_persubsite %>% 
  merge.data.frame(net_ghg_summary_central95_persubsite, by=c("pilot_site","status","subsite","subsite_latitude","subsite_longitude")) %>% 
  mutate(sd_rel_reduction_CO2=(sd_netCO2-sd_netCO2_95central)/sd_netCO2,
         sd_rel_reduction_CH4=(sd_netCH4-sd_netCH4_95central)/sd_netCH4,
         n_discard_CO2=n_netCO2-n_netCO2_95central,
         n_discard_CH4=n_netCH4-n_netCH4_95central) %>% 
  summarise(avg_relSD_reduct_CO2=mean(sd_rel_reduction_CO2),
            avg_ndiscard_CO2=mean(n_discard_CO2),
            avg_relSD_reduct_CH4=mean(sd_rel_reduction_CH4),
            avg_ndiscard_CH4=mean(n_discard_CH4),
            avg_n_subsite_CO2=mean(n_netCO2),
            avg_n_subsite_CH4=mean(n_netCH4))




