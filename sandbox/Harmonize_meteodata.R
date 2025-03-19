

# ---
# Authors: Miguel Cabrera
# Project: "RESTORE4Cs"
# date: "March 2025"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
#This script harmonizes the meteo data uploaded by the case-pilots to a format useful for GHG chamber temperature and pressure correction of fluxes. 

# 1. Compile Subsites lat/long centroids.

# 2. Harmonize various meteo-data files, retaining common variables, each station has its own unique station_id

# 3. Produce Map_meteo stations with assignment of each subsite to the closest meteo-station


rm(list = ls()) # clear workspace
cat("/014") # clear console

# install.packages("geosphere")
# ---- packages ----
library(tidyverse)
library(readxl)
library(lubridate)
library(ggplot2)
library(sf)#for geographical calculations
library(geosphere)

repo_root <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
files.sources = list.files(path = paste0(repo_root,"/functions"), full.names = T)
for (f in files.sources){source(f)}

# ---- Directories ----
dropbox_root <- "C:/Users/Miguel/Dropbox/RESTORE4Cs - Fieldwork/Data" # You have to make sure this is pointing to the write folder on your local machine
meteopath <- paste0(dropbox_root,"/Meteo")
fieldsheetpath <- paste0(dropbox_root,"/GHG/Fieldsheets")

#_________________####



#---1. Compile subsites info -----

#Load fieldsheets and calculate: single location per subsite (median lat&long), dates of visit 


# List GHG chamber fieldsheets in Dropbox and read them
myfieldsheets_list <- list.files(fieldsheetpath, pattern = "Fieldsheet-GHG.xlsx", all.files = T, full.names = T, recursive = T)

# Read all fieldsheets and put them in a single dataframe
fieldsheet <- read_GHG_fieldsheets(myfieldsheets_list)


#Get centroids of every subsite as latitude and longitude
centroids_subsites <- fieldsheet %>%
  select(subsite, latitude, longitude) %>% 
  mutate(site_subsite = substr(subsite, 4, 8)) %>%  # drop season from subsite 
  group_by(site_subsite) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%   # WGS 84 CRS
  summarise(centroid = st_centroid(st_union(geometry))) %>% # calculate centroids
  mutate(centroid_coords = st_coordinates(centroid)) %>% # transform into long lat columns
  select(site_subsite, centroid_coords) %>%
  st_drop_geometry() %>%
  mutate(
    latitude = centroid_coords[, 2],  # Extract latitudes
    longitude = centroid_coords[, 1]  # Extract longitudes
  ) %>%
  select(site_subsite, latitude, longitude) 



#Load meteo-stations coordinates:

names_meteo_metadata<- names(read_xlsx(path = paste0(meteopath,"/General Meteo-station details per Case-Pilot.xlsx"),col_names = T, n_max = 0))

meteo_metadata<- read_xlsx(path = paste0(meteopath,"/General Meteo-station details per Case-Pilot.xlsx"),skip=2, col_names = F, n_max = 20)

names(meteo_metadata)<- names_meteo_metadata

meteo_metadata<- meteo_metadata %>% 
  filter(!is.na(Site)) %>% 
  filter(`Air temperature`=="yes" & `Atmospheric pressure`=="yes") %>% 
  select(station_id, Latitude, Longitude) %>% 
  rename(lat_station=Latitude, lon_station=Longitude) %>% 
  distinct()



#Calculate and assign the closest station for each subsite

# Convert the data frames into sf objects
subsites_sf <- st_as_sf(centroids_subsites, coords = c("longitude", "latitude"), crs = 4326)
meteo_sf <- st_as_sf(meteo_metadata, coords = c("lon_station", "lat_station"), crs = 4326)

# Calculate distances and find the closest station for each subsite
subsite_station <- subsites_sf %>%
  st_join(meteo_sf, join = st_nearest_feature) %>% #join based on closest
  mutate(subsite_coords = st_coordinates(geometry)) %>%
  st_drop_geometry() %>%
  mutate(lat_subsite = subsite_coords[, 2], lon_subsite = subsite_coords[, 1]) %>%   #Keep subsite coords as columns
  select(-subsite_coords) %>% 
  merge.data.frame(meteo_metadata, by="station_id") %>% 
  mutate(distance_km = distVincentySphere(cbind(lon_subsite, lat_subsite), cbind(lon_station, lat_station))/1000) %>% #Calculate distance subsite-station
  select(site_subsite, lat_subsite, lon_subsite, station_id, distance_km)
  
rm(meteo_sf, subsites_sf, centroids_subsites,meteo_metadata)

#Save subsite-station correspondence
write.csv(subsite_station, file = paste0(meteopath,"/Formated_data/correspondence_meteostation.csv"),row.names = F)



#---2. Harmonize meteostation data -----

#Common variables: 

# station_id, 
# datetime_utc
# temp_c,
# Patm_hPa,
# precip_mm,
# humidity_percent,
# __________________________globalrad_UNITS_TBD 
# windspeed_ms,
# winddir_degrees

#optional variables: cloudcover_percent (CU)


#CAMARGUE
#Span: 2023-01-01 to 2024-10-28
#Single station for hourly data (includes solar radiation and atm pressure from another station, merged into a single file)

#Solar radiation in Joule per square cm

ca<- read_xlsx(paste0(meteopath,"/Meteo_CA/Meteo_TdV_hourly_2023_2024.xlsx"))
ca<- ca[,c(2,6,7,8,9,12,13,15,16)] #keep only variables of interest

#rename to something useful
names(ca)<- c("station_id","AAAAMMDDHH","precip_1/10mm","windspeed_ms","winddir_degrees","temp_c","humidity_percent","globalrad_joulepersquarecm","Patm_hPa")

ca_Arles<- ca %>% 
  filter(station_id=="ARLES") %>% 
  mutate(datetime_utc=as.POSIXct(as.character(AAAAMMDDHH),format="%Y%m%d%H",tz = "utc"),#Checked in website, times are utc
         station_id="CA_Arles",
         precip_mm=`precip_1/10mm`/10,
         Patm_hPa=as.numeric(Patm_hPa),
         globalrad_joulepersquarecm=as.numeric(globalrad_joulepersquarecm),
         windspeed_ms=as.numeric(windspeed_ms),
         winddir_degrees=as.numeric(winddir_degrees)) %>% 
  select(station_id, datetime_utc, temp_c,Patm_hPa, precip_mm,humidity_percent, globalrad_joulepersquarecm, #STILL TO DETERMINE UNIFIED UNITS FOR IRRADIANCE
         windspeed_ms,winddir_degrees)


rm(ca)


#CURONIAN
#Span: 2023-09-01 to 2024-09-30
#Hourly resolution
#Two stations provided, only 1 with both temp_c and Patm: cu_Nida
#Solar radiation requested but still not available


#Variables:
# observationTimeUtc - time of performed observations (UTC time zone).
# airTemperature - air temperature, °C.
# feelsLikeTemperature - perceived temperature, °C.
# windSpeed - wind speed, m/s.
# windGust - wind gust, m/s. Maximum gust per hour.
# windDirection - wind direction, °. Values: 0 - from the north, 180 - from the south, etc. i.e.
# cloudCover - cloud cover, %. Values: 0 - clear, 100 - cloudy. If cloud cover cannot be determined (for example due to fog), a null value is returned.
# seaLevelPressure - pressure at sea level, hPa.
# relativeHumidity - relative air humidity, %.
# precipitation - amount of precipitation, mm. Amount of precipitation per hour.
# conditionCode - weather conditions, code.


cu_Nida<- read_tsv(file = paste0(meteopath,"/Meteo_CU/Nida_202309-202409.txt")) %>% 
  mutate(station_id="CU_Nida") %>% 
  rename(datetime_utc=observationTimeUtc,temp_c=airTemperature,Patm_hPa=seaLevelPressure,
         precip_mm=precipitation, humidity_percent=relativeHumidity, 
         windspeed_ms=windSpeed, winddir_degrees=windDirection,cloudcover_percent=cloudCover) %>% 
  select(station_id, datetime_utc, temp_c,Patm_hPa, precip_mm,humidity_percent, #globalrad_joulepersquarecm, #NO DATA AVAILABLE YET
         windspeed_ms,winddir_degrees,
         cloudcover_percent)
  


#DANUBE
#3 stations supplied, 2 close to subsites: DA_Mahmudia, 	DA_Tulcea
