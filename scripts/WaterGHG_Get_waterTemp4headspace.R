#GET WATER TEMPERATURE FOR HEADSPACEs 

# ---
# Authors: Miguel Cabrera Brufau
# Project: "RESTORE4Cs"
# date: "August 2025"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
# This script extracts the obtains water temperature data from CNR team and assigns water-temperatures to headspace samples based on proximity of water-samples and headspace samples. 

#Steps: 

# DONE: 0. Add functions and copy from other scripts for fieldsheet import 

# DONE: 1. Extract from CNR's data: season, site,subiste, tide, watercode (uniqueID for water),date, time, lat, long, temp, (extra: o2 mg/L, conductivity microS/cm, salinity PSU, pH, )

# DONE: 2. get headspace sampling details, exetainerID, plotID correspondence

# DONE: 3. get headspacee sampling details from chamber fieldsheets 

# DONE: 4. Use proximity in space and local knowledge of water-body fragmentation to assign water temperature from Waterdata to headspace exetainerID. 


rm(list = ls()) # clear workspace
cat("/014") # clear console


# ---- packages ----
library(dplyr)
library(ggplot2)
library(ggpubr)
library(lme4)
library(hms)
library(lubridate)
library(stringr)
library(tidyverse)
library(geosphere)

repo_root <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
files.sources = list.files(path = paste0(repo_root,"/functions"), full.names = T)
for (f in files.sources){source(f)}

# ---- Directories ----
dropbox_root <- "C:/Users/Miguel/Dropbox/RESTORE4Cs - Fieldwork/Data/"
water_path <- paste0(dropbox_root,"Water/")

fieldsheet_path <- paste0(dropbox_root,"GHG/Fieldsheets")


  #1. Get waterdata------
#Import
allwater<- readxl::read_xlsx(paste0(water_path,"Water sampling and filtration_all data.xlsx"), sheet = 1)

str(allwater)


#Format and filter
#Save in a single table the appropriate data. 
allwatertemps<- allwater %>% 
  select(Season, Site, Subsite, Replicate, 
         `Label Sample ID`, #HT or LT denotes tide, some inconsistencies in subsitecombos for HT, in format and order of sub-codes, in separators (-, _)
         `Date dd/mm/yyyy`, 
         `Sampling time (UTC)`, #ISSUE: gets imported as POSIXct with defatult 1899-12-31 date
         `Latitude Y °N (decimal)`,
         `Longitude X °E (decimal)`,
         `Tw (ºC)`,
         Salinity, #PSU
         pH, 
         `DO (mg/L)`,
         `Water depth (m)`
         ) %>% 
  rename_with(.fn = tolower) %>% 
  rename(waterID=`label sample id`,
         water_depth=`water depth (m)`) %>% 
  mutate(hightide=grepl("HT",waterID),
         water_depth=100*water_depth # transform water depth to cm
         ) %>% 
  rename(date=`date dd/mm/yyyy`,
         utctime=`sampling time (utc)`,
         latitude=`latitude y °n (decimal)`,
         longitude= `longitude x °e (decimal)`,
         temp= `tw (ºc)`,
          sal=salinity,
         ph=ph,
         do_mgperlitre=`do (mg/l)`) %>% 
  mutate(timeonly=substr(as.character(utctime), 12, 19),
         datetime=(paste(as.character(date), timeonly)),
         datetime=as.POSIXlt.character(datetime, tz="UTC",format = "%Y-%m-%d %H:%M:%OS" )) %>% 
  filter(!is.na(latitude)) %>% 
  select(-c(timeonly, utctime)) %>% 
  #Remove duplicated details (all replicates taken in the same spot: same data)
  filter(!waterID%in%c("S1-CA-P2-R2","S1-CA-P2-R3"))



##----CHECKED high tides-----
#Check that the headspace samples all correspond to low-tide samples, that way we can more easily correct codes

#NO; there are high and low-tide headspaces in both watertemps and headsp_details
#We need to "duplicate" watertemps for high-tides where the subsite is ambiguous, some subsites are listed as A1/P1/R1, 

#Decission: high-tides of waterteam are not the same as rising-tides of GHGteam. use only low-tide temperatures for matching. Hightide temperatures are very uniform for all subsites, impute with a single value for each RI season. 

hightideswatertemp<- allwatertemps %>%  filter(hightide==T) %>% 
  group_by(site, season) %>% 
  summarise(temp=mean(temp, na.rm=T))
  


#Correct lowtidewater waterID codes . 
watertemps<- allwatertemps %>%  filter(hightide==F) %>% 
  #In water data, Samples of S1-DU-R1 and S1-DU-A2 are swapped (wrong codes), correct subsite content and remake waterID
  mutate(subsite=case_when(waterID%in%c("S1-DU-R1-R1","S1-DU-R1-R2","S1-DU-R1-R3")~"A2",
                           waterID%in%c("S1-DU-A2-R1","S1-DU-A2-R2","S1-DU-A2-R3")~"R1",
                           TRUE~subsite)) %>% 
  mutate(waterID= paste(season, site, subsite, replicate, sep = "-")) %>% 
  mutate(subsite_ID= paste(season, site, subsite, sep = "-")) %>% 
  select(-c(season,site, subsite)) %>% 
  select(subsite_ID, replicate, waterID, hightide, date,datetime, latitude, longitude,water_depth, temp, sal, ph, do_mgperlitre) %>% 
  mutate(date=as.character(date)) %>% 
  rename(subsiteID=subsite_ID)

watertemps

  
  #2. Get headspace sampling----
#Create UniqueIDnotime from headspace fieldsheets,import GHGfieldsheets and filter for headspace samplings. 

#Import fieldsheets exetainers field samples
fieldsheets_exe<-list.files(fieldsheet_path, pattern = "exetainers.xlsx", recursive = T, full.names = T)

for(i in fieldsheets_exe){
  
  a<- readxl::read_xlsx(i, trim_ws = T)
  a<- a[!is.na(a[,1]),]
  
  if(i==fieldsheets_exe[1]){A<- a}else {A<- rbind(A,a)}
}
rm(a,i)

#Are all exetainer_ID values in A different?
length(unique(A$exetainer_ID))==dim(A)[1]


#get only headspaces and atm (separate)
atm<- A %>% filter(`headspace/trapped/atm`=="atmosphere")

#get headspace exetainerID-plotID correspondence
headsp<- A %>% filter(`headspace/trapped/atm`=="headspace") %>% 
  select(subsite_ID, plot_ID, exetainer_ID) %>%
  rename(subsiteID=subsite_ID, plotID= plot_ID, exetainerID=exetainer_ID) %>% 
  mutate(plotID=paste(subsiteID, plotID, sep = "-")) %>% 
  select(-subsiteID)

headsp


#Read all GHGfieldsheets (using function from Camille)
fieldsheets<- list.files(path= fieldsheet_path, recursive = T, pattern = "GHG.xlsx",full.names = T)
field<- read_GHG_fieldsheets(fieldsheets)


#Filter for first plotID (duplicated plotID for transp/dark inubations, but same plot) 
ghgwater<- field %>% 
  mutate(plotID=paste(subsite, plot_id, sep = "-")) %>% 
  #Filter to get only first ocurrence of any plotID
  group_by(plotID) %>% 
  filter(row_number()==1) %>% 
  ungroup()


#Get headspace sampling details (plotID, exetainerID, date, datetime, latitude, longitude, water_depth)
headsp_details<- headsp %>% 
  left_join(ghgwater, by="plotID") %>% 
  mutate(datetime= as.POSIXlt.character(paste(date, start_time), tz="UTC",format = "%Y-%m-%d %H:%M:%OS")) %>%
  mutate(date=as.character(date)) %>% 
  mutate(subsiteID=substr(plotID, 1, 8)) %>% 
  mutate(high_tide=grepl(pattern = "Rising-tide", comments)) %>% 
  select(subsiteID,plotID, exetainerID, date, datetime, latitude, longitude, water_depth, high_tide,comments)
    
headsp_details %>% filter(is.na(latitude))


rm(A, allwater, allwatertemps)

  #3. Correspondence ------

head(watertemps)
head(headsp_details)



# Ensure datetime and date columns are properly formatted
watertemps <- watertemps %>%
  mutate(datetime = ymd_hms(datetime), date = ymd(date)) %>% 
  #Ensure watertemps to match always contain a temperature (avoids matching headspace to NA temperature)
  filter(!is.na(temp))


headsp_details <- headsp_details %>%
  mutate(datetime = ymd_hms(datetime), date = ymd(date))


# Step 1: Left join to retain all headsp_details_lowtide rows
distance_df <- headsp_details %>%
  left_join(watertemps, by = "subsiteID", suffix = c("_headsp", "_water"), relationship = "many-to-many") %>%
  mutate(
    geo_distance_m = if_else(
      !is.na(latitude_water) & !is.na(longitude_water),
      distGeo(
        matrix(c(longitude_water, latitude_water), ncol = 2),
        matrix(c(longitude_headsp, latitude_headsp), ncol = 2)
      ),
      NA_real_
    )
  )

# Step 2: Filter for closest match per exetainerID (if any)
closest_match_df <- distance_df %>%
  group_by(exetainerID) %>%
  slice_min(geo_distance_m, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    date_diff_days = as.numeric(date_headsp - date_water),
    datetime_diff_minutes = as.numeric(difftime(datetime_headsp, datetime_water, units = "mins")),
    depth_diff = abs(water_depth_headsp - water_depth_water)
  )



#Any exetainerID with more than 1 correspondence?
closest_match_df %>% group_by(exetainerID) %>% summarise(n=n()) %>% filter(n>1)

#Any exetainerID without a match?
headsp_details$exetainerID[!headsp_details$exetainerID%in%closest_match_df$exetainerID]



##Check matching-------

#Check farthest matches: 
closest_match_df %>% 
  arrange(desc(geo_distance_m)) %>% 
  select(subsiteID, plotID, exetainerID, geo_distance_m, waterID, temp)

#OK, Errors Already corrected: S1-DU-R1 and S1-DU-A2 were matched with ~23km distance (waterdata had swapped subsite codes)

#Check temperature NAs:
closest_match_df %>% 
  filter(is.na(temp))
#ONE headspace sample taken in a rain-pond in the agricultural field (not sampled by water-team)
#S2-DA-A1 (pond in agricultural field) without any watersample to match:--> discard headspace


#Split high-tides and lowtides:

lowtides_filled<- closest_match_df %>% 
  filter(high_tide==F)

hightides_filled<- closest_match_df %>% 
  filter(high_tide==T)


##Impute High-tides----
#Given that the ghgteam and waterteam do not share the same high-tides locations (and timings), along with high-tide temperatures being highly uniform regardless of subsite, we will overwrite all high-tide headspaces with the average water temperature of high-tide samples for each season. removing the rest of details (DO, sal, ph,....,)
hightides_filled


hightides_filled<-hightides_filled %>% 
  mutate(season=substr(subsiteID, 1,2)) %>% 
  mutate(temp= case_when(season=="S1"~hightideswatertemp %>% filter(season=="S1") %>% pull(temp),
                         season=="S2"~hightideswatertemp %>% filter(season=="S2") %>% pull(temp),
                         season=="S3"~hightideswatertemp %>% filter(season=="S3") %>% pull(temp),
                         season=="S4"~hightideswatertemp %>% filter(season=="S4") %>% pull(temp))) %>% 
  #Remove wrong columns (arising from improper match)
  select(-c(season, replicate, waterID, hightide, date_water, datetime_water, latitude_water, longitude_water, water_depth_water, sal, ph, do_mgperlitre, geo_distance_m, date_diff_days, datetime_diff_minutes, depth_diff))


str(hightides_filled)



#Re-join high and lowtides
filled_headsp_details<- lowtides_filled %>% 
  full_join(hightides_filled) %>% 
  #Remove 1 headspace with NA temp, single sample from DA agricultural field (no temperature available, rainpond only sampled during S2, not of interest)
  filter(!is.na(temp)) %>% 
  select(-c(replicate, hightide, latitude_water, longitude_water, date_water, datetime_water, geo_distance_m, date_diff_days, datetime_diff_minutes, water_depth_water,depth_diff)) #remove waterdata-variables not of interest


 str(filled_headsp_details)

 
   #4. Save details -----
 
 
 
 write.csv(x = filled_headsp_details, file = paste0(dropbox_root,"GHG/Headspaces/","headspace_temperatures.csv"),row.names = F)
 