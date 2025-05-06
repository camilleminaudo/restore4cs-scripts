

# ---
# Authors: Miguel Cabrera
# Project: "RESTORE4Cs"
# date: "March 2025"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
#This script harmonizes the meteo data uploaded by the case-pilots to a format useful for GHG chamber temperature and pressure correction of fluxes. 

# 1. Compile Subsites lat/long centroids. Produce Map_meteo stations with assignment of each subsite to the closest meteo-station

# 2. Harmonize various meteo-data files, retaining common variables, each station has its own unique station_id

# 3. Combine dataset, check and subset dates. Save the data for sampling days and for whole data-availability
#All subsites hourly except for RI (10-minute resolution)

# 4. Retrieve modelled historical data for each subsite individually from open-meteo.com


#All data is saved in Meteo / Formated_data as csv files





rm(list = ls()) # clear workspace
cat("/014") # clear console

# install.packages("geosphere")
# install.packages("openmeteo")

# ---- packages ----
library(tidyverse)
library(readxl)
library(lubridate)
library(ggplot2)
library(sf)#for geographical calculations
library(geosphere)
library(openmeteo)
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
write.csv(subsite_station, file = paste0(meteopath,"/Formated_data/correspondence_subsite-meteostation.csv"),row.names = F)



#---2. Harmonize meteostation data -----

#Common variables: 

# station_id, 
# datetime_utc
# temp_c,
# Patm_hPa,
# precip_mm,
# humidity_percent,
# globalrad_wperm2
# windspeed_ms,
# winddir_degrees
#cloudcover_percent (CU)

##CA####
#CAMARGUE
#Span: 2023-01-01 to 2024-10-28
#Frequency: hourly
#Timezone: "UTC" checked
#Radiation units: Global rad (Joule /cm2)

#Single station for hourly data (includes solar radiation and atm pressure from another station, merged into a single file, pretend is all from the same station)

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
         winddir_degrees=as.numeric(winddir_degrees),
         cloudcover_percent=NA_real_) %>% 
  mutate(globalrad_Wperm2= globalrad_joulepersquarecm*(10000/3600)) %>%  #multiply 10^4 to transform to meters, divide by 3600 s in one hour. This assumes the joule per cm2 measure was integrated over 1h. 
  select(station_id, datetime_utc, temp_c,Patm_hPa, precip_mm,humidity_percent, globalrad_Wperm2,
         windspeed_ms,winddir_degrees,cloudcover_percent)


rm(ca)




##CU####
#CURONIAN
#Span: 2023-09-01 to 2024-09-30
#Frequency: hourly
#Timezone: "UTC" checked
#Radiation units: NA (data not available)
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
  mutate(station_id="CU_Nida",
         globalrad_Wperm2=NA_real_) %>% 
  rename(datetime_utc=observationTimeUtc,temp_c=airTemperature,Patm_hPa=seaLevelPressure,
         precip_mm=precipitation, humidity_percent=relativeHumidity, 
         windspeed_ms=windSpeed, winddir_degrees=windDirection,cloudcover_percent=cloudCover) %>% 
  select(station_id, datetime_utc, temp_c,Patm_hPa, precip_mm,humidity_percent,
         globalrad_Wperm2, #NO DATA AVAILABLE YET
         windspeed_ms,winddir_degrees,
         cloudcover_percent)
  




##DA####
#DANUBE
#Stations: Two valid
#Span: 2023-09-01 to 2024-09-30
#Frequency: hourly
#Timezone: UTC (Contantin email)
#Radiation units: "global_radiation" KJ/m2 (constantin email)


#3 stations supplied, 2 close to subsites: DA_Mahmudia, 	DA_Tulcea 
#Hourly data

# da_Mahmudia and da Tulcea, both with the same columns and names
#time (YYYY-MM-DD HH:MM:SS),
#air_temp (ºC)
#pressure: hPa
#precipitation (units?), assumed mm/h
#wind_speed (units?), assumed m/s
#wind_direction (degrees, 1-360) substitute 360 N with 0 N
#humidity (percent humidity)
#global_radiation (Units?):Kj/m2 (according to Constantin email) 

#1J/s = 1 W

#Apparently KJ /m2, reported every hour, most likely, KJ/(m2*h), multiply by 1000 and divide by 3600 s/h

da_Mahmudia<- read.csv(paste0(meteopath, "/Meteo_DA/Mahmudia.csv")) %>% 
  rename(datetime_utc=time, temp_c=air_temp, Patm_hPa=pressure, precip_mm=precipitation,windspeed_ms=wind_speed,
         winddir_degrees=wind_direction, humidity_percent=humidity, globalrad_KJperm2=global_radiation) %>% 
  mutate(station_id="DA_Mahmudia", 
         datetime_utc=as.POSIXct(datetime_utc, tz="utc", format="%Y-%m-%d %H:%M:%S"),
         winddir_degrees=ifelse(winddir_degrees==360,0,winddir_degrees),
         cloudcover_percent=NA_real_,
         globalrad_Wperm2=globalrad_KJperm2*1000/3600) %>% 
  select(station_id, datetime_utc, temp_c,Patm_hPa, precip_mm, humidity_percent, windspeed_ms,winddir_degrees, globalrad_Wperm2,cloudcover_percent) %>% 
  #Remove duplicate times if any (keep only first)
  arrange(datetime_utc) %>%
  filter(!duplicated(datetime_utc))


da_Tulcea<- read.csv(paste0(meteopath, "/Meteo_DA/Tulcea.csv")) %>% 
  rename(datetime_utc=time, temp_c=air_temp, Patm_hPa=pressure, precip_mm=precipitation,windspeed_ms=wind_speed,
         winddir_degrees=wind_direction, humidity_percent=humidity, globalrad_KJperm2=global_radiation) %>% 
  mutate(station_id="DA_Tulcea", 
         datetime_utc=as.POSIXct(datetime_utc, tz="utc", format="%Y-%m-%d %H:%M:%S"),
         winddir_degrees=ifelse(winddir_degrees==360,0,winddir_degrees),
         cloudcover_percent=NA_real_,
         globalrad_Wperm2=globalrad_KJperm2*1000/3600) %>% 
  select(station_id, datetime_utc, temp_c,Patm_hPa, precip_mm, humidity_percent, windspeed_ms,winddir_degrees, globalrad_Wperm2,cloudcover_percent) %>% 
  #Remove duplicate times if any (keep only first)
  arrange(datetime_utc) %>%
  filter(!duplicated(datetime_utc))
  
  



##DU####
#DUTCH Delta
#Stations: Two valid, same file
#Span: 2023-09-01 to 2024-10-01
#Frequency: hourly
#Timezone: UTC checked in website
#Radiation units: Global radiation (in J/cm2) 


#Two stations supplied, in the same file. Identified by 
#DU_323, DU_340

#32 lines describing the data at the top of the file:
# Opmerking: door stationsverplaatsingen en veranderingen in waarneemmethodieken zijn deze tijdreeksen van uurwaarden mogelijk inhomogeen! Dat betekent dat deze reeks van gemeten waarden niet geschikt is voor trendanalyse. Voor studies naar klimaatverandering verwijzen we naar de gehomogeniseerde dagreeksen <http://www.knmi.nl/nederland-nu/klimatologie/daggegevens> of de Centraal Nederland Temperatuur <http://www.knmi.nl/kennis-en-datacentrum/achtergrond/centraal-nederland-temperatuur-cnt>.
# 
# SOURCE: ROYAL NETHERLANDS METEOROLOGICAL INSTITUTE (KNMI)
# Comment: These time series are inhomogeneous because of station relocations and changes in observation techniques. As a result these series are not suitable for trend analysis. For climate change studies we refer to the homogenized series of daily data <http://www.knmi.nl/nederland-nu/klimatologie/daggegevens> or the Central Netherlands Temperature <http://www.knmi.nl/kennis-en-datacentrum/achtergrond/centraal-nederland-temperatuur-cnt>.
# 
# STN         LON(east)   LAT(north)  ALT(m)      NAME
# 323         3.884       51.527      1.40        Wilhelminadorp
# 340         4.342       51.449      19.20       Woensdrecht 

#Column names, descriptions and units: 
#"# STN": numeric code of station 
#"YYYYMMDD": 
#"HH": hour of day 1-24
#"DD": Mean wind direction (in degrees) during the 10-minute period preceding the time of observation (360=north; 90=east; 180=south; 270=west; 0=calm 990=variable)
#"FH": Hourly mean wind speed (in 0.1 m/s)
#"FF": Mean wind speed (in 0.1 m/s) during the 10-minute period preceding the time of observation
#"FX": Maximum wind gust (in 0.1 m/s) during the hourly division
#"T": Temperature (in 0.1 degrees Celsius) at 1.50 m at the time of observation
#"T10N": Minimum temperature (in 0.1 degrees Celsius) at 0.1 m in the preceding 6-hour period
#"TD": Dew point temperature (in 0.1 degrees Celsius) at 1.50 m at the time of observation
#"SQ": Sunshine duration (in 0.1 hour) during the hourly division; calculated from global radiation (-1 for <0.05 hour)
#"Q": Global radiation (in J/cm2) during the hourly division
#"DR": Precipitation duration (in 0.1 hour) during the hourly division
#"RH": Hourly precipitation amount (in 0.1 mm) (-1 for <0.05 mm)
#"P": Air pressure (in 0.1 hPa) reduced to mean sea level; at the time of observation
#"VV": Horizontal visibility at the time of observation (0=less than 100m; 1=100-200m; 2=200-300m;...; 49=4900-5000m; 50=5-6km; 56=6-7km; 57=7-8km; ...; 79=29-30km; 80=30-35km; 81=35-40km;...; 89=more than 70km)
#"N": Cloud cover (in octants); at the time of observation (9=sky invisible)
#"U": Relative atmospheric humidity (in percents) at 1.50 m at the time of observation
#"WW": Present weather code (00-99); description for the hourly division. http://bibliotheek.knmi.nl/scholierenpdf/weercodes_Nederland
#"IX": Indicator present weather code (1=manned and recorded (using code from visual observations); 2;3=manned and omitted (no significant weather phenomenon to report; not available); 4=automatically recorded (using code from visual observations); 5;6=automatically omitted (no significant weather phenomenon to report; not available); 7=automatically set (using code from automated observations)
#"M": Fog 0=no occurrence; 1=occurred during the preceding hour and/or at the time of observation
#"R": Rainfall 0=no occurrence; 1=occurred during the preceding hour and/or at the time of observation
#"S": Snow 0=no occurrence; 1=occurred during the preceding hour and/or at the time of observation
#"O": Thunder  0=no occurrence; 1=occurred during the preceding hour and/or at the time of observation
#"Y":  Ice formation 0=no occurrence; 1=occurred during the preceding hour and/or at the time of observation

du<- read_csv(paste0(meteopath, "/Meteo_DU/Dutch_weather.txt"),skip = 32)

du_selection<- c("# STN",#add "DU_" to "# STN" and rename to station_id 
                 "YYYYMMDD", #yearmonthday, ok
                 "HH",# hour of day, 1-24 transform with yearmonthday to datetime_utc
                 "T", #air temp, divide by 10 and rename to temp_c
                 "P", #air pressure, divide by 10 and rename to Patm_hPa
                 "RH",# precipitation, subsitute -1 for 0, and divide by 10, rename to precip_mm
                 "U", #relative humidity in percent, rename to humidity_percent
                 "Q", #Global radiation (in J/cm2), rename to globalrad_joulepersquarecm
                 "FH", #Hourly mean wind speed (in 0.1 m/s), divide by 10 and rename to windspeed_ms
                 "DD", #Mean wind direction (in degrees) during the 10-minute period preceding the time of observation (360=north; 90=east; 180=south; 270=west; 0=calm 990=variable), set 0 and 990 to NA, substitute 360 by 0 and rename to winddir_degrees
                 "N") #Cloud cover (in octants); at the time of observation (9=sky invisible due to fog, set to NA). TRANSFORM TO cloudcover_percent (0=clear, 100= total cover)


#Formatting dutch weather

du_format <- du[, du_selection] %>% 
  mutate(station_id=paste0("DU_",`# STN`)) %>% #Create station_id
  #Adjust and combine YYYYMMDD and HH (1-24) to create datetime_utc
  mutate(HH = if_else(HH == 24, 0, HH),
         YYYYMMDD = if_else(HH == 0, (as.Date(as.character(YYYYMMDD), format = "%Y%m%d") + 1), as.Date(as.character(YYYYMMDD), format = "%Y%m%d")),
         datetime_utc = ymd_h(paste(YYYYMMDD, HH),tz = "utc")) %>% 
  mutate(temp_c=`T`/10, #get T in ºC
         Patm_hPa=P/10, #get P in hPa
         precip_mm=if_else(RH==-1, 0, RH/10),#substitute "<0.5mm" code (-1) with 0, and transform to mm units
         humidity_percent=U,
         windspeed_ms=FH/10, 
         winddir_degrees=case_when(DD%in%c(0,990)~NA_real_, #set 0=calm 990=variable to NA
                                   DD==360~0, #set North from 360 to 0
                                   TRUE~DD),
         cloudcover_percent=ifelse(N==9, NA_real_, N/8*100),#Set 9 (sky not vissible fog/heavysnow) to NA, transform to % 
         globalrad_Wperm2=Q*(10000/3600)) %>% #Transform to W/m2 (multiply to get m2 and divide by 3600s in an hour)
  select(station_id, datetime_utc, temp_c,Patm_hPa, precip_mm, humidity_percent, windspeed_ms,winddir_degrees,
         globalrad_Wperm2, cloudcover_percent)
  


#Separate into stations: 
du_323<- du_format %>% 
  filter(station_id=="DU_323")

du_340<- du_format %>% 
  filter(station_id=="DU_340")

rm(du, du_format, du_selection)





##RI####
#RIA DE AVEIRO
#Stations: one
#Span: 2023-09-01 to 2024-09-30
#Frequency: 10 minutes
#Timezone: RI data is in local time (there were jumps consistent with summertime adjustments in sunrise/sunset data) convert to UTC at import
#Radiation units: "Rsolar"  globalrad_Wperm2, infered from values magnitudes

#RIA DE AVEIRO
#Time resolution: 10minutes

#Variables: 
#"TIME"     "Prec_Tot" "TA02"     "TA10"     "TA20"     "TA30"     "HR02"     "HR10"     "HR20"     "HR30"    "V10m"     "V10m_Max" "V10m_Dir" "V30m"     "V30m_Max" "V30m_Dir" "Rsolar"   "Press"   

#We will use the variables at 10m, 
#Assumed units: 
#TIME : it is tz WET 
#Prec_tot: mm accummulated in last 10 minutes, it will be integrated hourly (i.e. precipitation in the last hour) even when reporting it every 10 minutes. 
#TA10: atmospheric temperature at 10m (ºC)
#HR10: percent relative humidity at 10m (%)
#V10m: windspeed at 10m (m/s)
#V10m_Dir: winddirection at 10m height (degrees)
#Rsolar: globalrad_ ------------------------------------UNITS TO CHECK
#Press: atmospheric pressure (hPa)



ri<- read_xlsx(path=paste0(meteopath,"/Meteo_RI/Dados_09_2023_a_09_2024.xlsx"))

ri_selection<- c("TIME", "Prec_Tot", "TA10", "HR10", "V10m", "V10m_Dir", "Rsolar", "Press")

ri_EMA<- ri[,ri_selection] %>% 
  mutate(station_id="RI_EMA",
         datetime_utc=force_tzs(TIME, tzones="WET",tzone_out = "UTC"), #forze the timezone change to utc (read_xlsx always assumes UTC)
         temp_c=TA10,
         Patm_hPa=Press,
         precip_mm= Prec_Tot +  
           lag(Prec_Tot, 1, default = NA) + 
           lag(Prec_Tot, 2, default = NA) + 
           lag(Prec_Tot, 3, default = NA) + 
           lag(Prec_Tot, 4, default = NA) + 
           lag(Prec_Tot, 5, default = NA), #Integrate precip_mm to each hour (every 10 minutes, the sum of precipitation in last hour is reported)
         humidity_percent=HR10,
         windspeed_ms=V10m,
         winddir_degrees=ifelse(V10m_Dir==360, 0, V10m_Dir),# substitute 360 N for 0 N
         globalrad_Wperm2=Rsolar,#Global radiation is in W/m2
         cloudcover_percent=NA_real_) %>% 
         select(station_id, datetime_utc, temp_c,Patm_hPa, precip_mm, humidity_percent, windspeed_ms,winddir_degrees, globalrad_Wperm2,cloudcover_percent) %>% 
  #Remove duplicate times if any (keep only first)
  arrange(datetime_utc) %>%
  filter(!duplicated(datetime_utc))


####hourly_RI: TO-DO------
#Produce another dataset for RI with hourly resolution (to have an homogeneous dataset across case-pilots)

# ri_EMAh<- ri %>% 
  #Generate truncated datetime_utc only hour resolution, group by truncated_datetime, summarise (mean, na.rm=T), EXCEPTIONS: precip_mm: subsitute NAs with 0 and sum per hour, winddir_degrees: special function to search trigonometric mean (rad transformation so that mean of 0.1 and 359.9 produces 0), Globalrad_units: Investigate first how the global radiation is reported, to see if it makes sense to do mean or sum (if is cumulative radiation E/area or energyflux/area). 
  
rm(ri,ri_selection)         
         




##VA####
#VALENCIA
#Stations: one
#Span: 2023-01-01 to 2024-07-31
#Frequency: hourly
#Timezone: UTC (assumed)
#Radiation units: Solar radiation" W/squaremeter

#VALENCIA: 

#only 1 dataset provided, no precipitation data
va_CEAV<- read_xlsx(paste0(meteopath,"/Meteo_VA/Meteo Marjal Moros - Sagunto CEAV.xlsx"),sheet = "Data", skip = 2)
units_va_CEAV<- as.character(read_xlsx(paste0(meteopath,"/Meteo_VA/Meteo Marjal Moros - Sagunto CEAV.xlsx"),sheet = "Data",n_max = 1)[1,])
names(va_CEAV)<- names(read_xlsx(paste0(meteopath,"/Meteo_VA/Meteo Marjal Moros - Sagunto CEAV.xlsx"),sheet = "Data",n_max = 1))

head(va_CEAV)
print(units_va_CEAV)

va_selection<- c("Date", "HORA", "Wind velocity", "Wind direction", "Temp.", "Relative Hum", "Pressure", "Solar radiation")

#Variables and units
#Date: dd/mm/yyyy
#HORA: 0-23 ----------------CHECK TIMEZONE
#"Wind velocity": m/s 
#"Wind direction": degrees 0-360  (sustitute 360 with 0) 
#"Temp.": ºC
#"Relative Hum" %
#"Pressure": mbar=hPa
#"Solar radiation" W/squaremeter

#Missing: precipitation, cloudcover, 

va_CEAV<- va_CEAV[,va_selection] %>% 
  mutate(station_id="VA_CEAV",
         datetime_utc=ymd_h(paste(Date, HORA),tz = "utc"),
         temp_c=Temp.,
         Patm_hPa=Pressure,
         humidity_percent=`Relative Hum`,
         windspeed_ms=`Wind velocity`,
         winddir_degrees=ifelse(`Wind direction`==360, 0, `Wind direction`),
         globalrad_Wperm2=`Solar radiation`,
         cloudcover_percent=NA_real_,
         precip_mm=NA_real_
         ) %>% 
  select(station_id, datetime_utc, temp_c,Patm_hPa,
         # precip_mm,  # MISSING PRECIPITATION
         humidity_percent, windspeed_ms,winddir_degrees, globalrad_Wperm2,precip_mm,cloudcover_percent)


rm(va_selection, units_va_CEAV)
  

#---3. Check, subset and save -----

#Solar radiation can be reported in Joules per area or watt per area. The difference is that Joules is energy and watt is power (energy per s). 

#We will transform all units to global irradiance: SI: W per m^2 (instantaneous measure, radiative flux/area)
#1W =1J per s

#IF irradiance is reported in jules / area, it is important to know how much time was the irradiance integrated for (usually this will be the frequency of measures). 

#Watt/m^2: global irradiance ranges: Average annual solar radiation arriving at the top of the Earth's atmosphere is roughly 1361 W/m2.[34] The Sun's rays are attenuated as they pass through the atmosphere, leaving maximum normal surface irradiance at approximately 1000 W/m2 at sea level on a clear day. When 1361 W/m2 is arriving above the atmosphere (when the Sun is at the zenith in a cloudless sky), direct sun is about 1050 W/m2, and global radiation on a horizontal surface at ground level is about 1120 W/m2.[35] The latter figure includes radiation scattered or reemitted by the atmosphere and surroundings.

#To transform joulespercm2 to wperm2 we have assumed that the frequency of measures is the integration time of solar energy (1h=3600s)


#JOIN data from all subsites AND CHECK THAT all values are consistent with their reported units
all<- rbind(ca_Arles,cu_Nida, da_Mahmudia, da_Tulcea, du_323, du_340, ri_EMA, va_CEAV)


all %>% 
  mutate(site_id=substr(station_id,1,2)) %>% 
  mutate(variable="Global Radiation (W per m^2)") %>% 
  filter(!is.na(globalrad_Wperm2))%>% 
  group_by(site_id,variable) %>%
  select(globalrad_Wperm2) %>% 
  summarise(avg_=mean(globalrad_Wperm2,na.rm=T),
            sd_=sd(globalrad_Wperm2,na.rm=T),
            max_=max(globalrad_Wperm2,na.rm=T),
            median_=median(globalrad_Wperm2,na.rm=T))


all %>% 
  filter(!is.na(globalrad_Wperm2))%>% 
  mutate(site=substr(station_id,1,2)) %>% 
  filter(globalrad_Wperm2>0) %>% 
  ggplot(aes(x=site, y=globalrad_Wperm2))+
  geom_boxplot()+
  ggtitle("Overview, non-zero Global radiation (W per m^2)")

#DANUBE solar radiation DATA was reported in KJ/m2 (according to Constantin), we inferred that it was integrated over 1 hour and, after unit transformation to W (J/s) per squaremeter, the values match the expected results. 


#Get combination of case_pilot& samplingday to subset full dataset to only dates of sampling:
site_daystokeep <- fieldsheet %>%
  mutate(season_site = substr(subsite, 1, 5)) %>%
  group_by(season_site) %>%
  mutate(start_date = min(date),
         end_date = max(date)) %>%
  select(season_site, start_date, end_date) %>%
  distinct() %>%
  mutate(sampling_days = purrr::map2(start_date, end_date, seq.Date, by = "day")) %>%
  unnest(sampling_days) %>% 
  ungroup() %>% 
  mutate(site=substr(season_site,4,5)) %>% 
  select(site, sampling_days) %>% 
  distinct() %>% 
  mutate(site_daystokeep=paste(site,sampling_days,sep = "_")) %>% 
  pull(site_daystokeep)
  

#Subset dataset for sampling dates at each case_pilot (keeping different station_id per case_pilot if present)
meteo_sampling_days<- all %>% 
  mutate(site_day=paste(substr(station_id, 1,2),as_date(datetime_utc), sep = "_")) %>% 
  filter(site_day%in%site_daystokeep) %>% 
  select(-site_day) %>% 
  arrange(datetime_utc,station_id)

#Check is any datetime_utc duplicated?
meteo_sampling_days %>% 
  group_by(station_id,datetime_utc) %>% 
  summarise(n=n(), .groups = "drop") %>% 
  filter(n>1)

all %>% 
  group_by(station_id,datetime_utc) %>% 
  summarise(n=n(), .groups = "drop") %>% 
  filter(n>1)


#Save meteo_sampling_days
write.csv(meteo_sampling_days, paste0(meteopath, "/Formated_data/allmeteostations_onlysamplingdays.csv"),row.names = F)

#Save meteo_all
write.csv(all, paste0(meteopath, "/Formated_data/allmeteostations_alldays.csv"), row.names = F)

rm(ca_Arles, cu_Nida,da_Mahmudia,da_Tulcea,du_323,du_340, ri_EMA, va_CEAV)



#4. API retrieve open-meteo-----

#Open-meteo provides historical meteorological data from re-analysis and forecast at specific locations with 9km spatial resolution and hourly. 

 
#Set start and end dates to retrieve: full year for all subsites sept23-sept24
start<- "2023-09-01"
end<- "2024-08-31"

#Create a table for API-download from open-meteo.com with site_subsite and location(lat,lon):
subsite_locations<- subsite_station %>% 
  rowwise() %>% 
  mutate(location=list(c(lat_subsite, lon_subsite))) %>% 
  select(site_subsite, location) 

#Create list of variables to retrieve
hourly_variables<- c("temperature_2m",
                     "relative_humidity_2m",
                     "pressure_msl",
                     "surface_pressure",
                     "precipitation",
                     "cloud_cover",
                     "wind_speed_10m",
                     "wind_direction_10m",
                     "diffuse_radiation",
                     "shortwave_radiation",
                     "soil_temperature_0_to_7cm",
                     "soil_moisture_0_to_7cm"
                     )



#This loop does not work (exceeded the Minutely API request limit)
for(i in subsite_locations$site_subsite){
print(i)
  
location=unlist(subsite_locations[which(subsite_locations$site_subsite==i),]$location)

#Retrieve hourly variables
openmeteo_subsite <- weather_history(location = location,
                              start=start,end=end,
                              response_units = list(temperature_unit = "celsius", 
                                                    windspeed_unit = "ms", 
                                                    precipitation_unit = "mm"),
                              hourly = hourly_variables,
                              timezone = "UTC")

#Wait for a bit to avoid exceeding API request limit
Sys.sleep(2.5)

#Add site_subsite identifier (eg. CA-A1, DU-R2,...)
openmeteo_subsite$site_subsite<- i

#Combine openmeteo_substie with rest
if(i==subsite_locations$site_subsite[1]){
  openmeteo_all<- openmeteo_subsite }else{
    openmeteo_all<- rbind(openmeteo_all, openmeteo_subsite)}

rm(openmeteo_subsite)

}

#Data from different subsites is very similar for most cases (VA is the same for all variables except for surface pressure and relative humidity)
openmeteo_all %>% 
  filter(grepl("VA", site_subsite)) %>%
  pivot_longer(cols = -c(datetime, site_subsite), names_to = "var",values_to = "value") %>% 
  group_by(datetime, var) %>% 
  summarise(sd_=sd(value)) %>% 
  filter(sd_>0) %>% 
  ungroup() %>% 
  select(var) %>% 
  distinct()



#Save openmeteo_all
write.csv(x = openmeteo_all, file = paste0(meteopath,"/Formated_data/openmeteo_retrieved_per_subsite.csv"),row.names = F)



#OPENMETEO citation: 
# Zippenfenig, P. (2023). Open-Meteo.com Weather API [Computer software]. Zenodo. https://doi.org/10.5281/ZENODO.7970649
# 
# Hersbach, H., Bell, B., Berrisford, P., Biavati, G., Horányi, A., Muñoz Sabater, J., Nicolas, J., Peubey, C., Radu, R., Rozum, I., Schepers, D., Simmons, A., Soci, C., Dee, D., Thépaut, J-N. (2023). ERA5 hourly data on single levels from 1940 to present [Data set]. ECMWF. https://doi.org/10.24381/cds.adbb2d47
# 
# Muñoz Sabater, J. (2019). ERA5-Land hourly data from 2001 to present [Data set]. ECMWF. https://doi.org/10.24381/CDS.E2161BAC
# 
# Schimanke S., Ridal M., Le Moigne P., Berggren L., Undén P., Randriamampianina R., Andrea U., Bazile E., Bertelsen A., Brousseau P., Dahlgren P., Edvinsson L., El Said A., Glinton M., Hopsch S., Isaksson L., Mladek R., Olsson E., Verrelle A., Wang Z.Q. (2021). CERRA sub-daily regional reanalysis data for Europe on single levels from 1984 to present [Data set]. ECMWF. https://doi.org/10.24381/CDS.622A565A
