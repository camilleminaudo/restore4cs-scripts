

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
         winddir_degrees=as.numeric(winddir_degrees)) %>% 
  select(station_id, datetime_utc, temp_c,Patm_hPa, precip_mm,humidity_percent, globalrad_joulepersquarecm, #STILL TO DETERMINE UNIFIED UNITS FOR IRRADIANCE
         windspeed_ms,winddir_degrees)


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
  mutate(station_id="CU_Nida") %>% 
  rename(datetime_utc=observationTimeUtc,temp_c=airTemperature,Patm_hPa=seaLevelPressure,
         precip_mm=precipitation, humidity_percent=relativeHumidity, 
         windspeed_ms=windSpeed, winddir_degrees=windDirection,cloudcover_percent=cloudCover) %>% 
  select(station_id, datetime_utc, temp_c,Patm_hPa, precip_mm,humidity_percent, #globalrad_joulepersquarecm, #NO DATA AVAILABLE YET
         windspeed_ms,winddir_degrees,
         cloudcover_percent)
  




##DA####
#DANUBE
#Stations: Two valid
#Span: 2023-09-01 to 2024-09-30
#Frequency: hourly
#Timezone: UNKNOWN
#Radiation units: "global_radiation" UNKNOWN


#3 stations supplied, 2 close to subsites: DA_Mahmudia, 	DA_Tulcea 
#Hourly data

# da_Mahmudia and da Tulcea, both with the same columns and names
#time (YYYY-MM-DD HH:MM:SS), assumed utc CONFIRM
#air_temp (ºC)
#pressure: hPa
#precipitation (units?), assumed mm/h
#wind_speed (units?), assumed m/s
#wind_direction (degrees, 1-360) substitute 360 N with 0 N
#humidity (percent humidity)
#global_radiation (Units?) NO UNITS ASSUMED, will check values to infere

da_Mahmudia<- read.csv(paste0(meteopath, "/Meteo_DA/Mahmudia.csv")) %>% 
  rename(datetime_utc=time, temp_c=air_temp, Patm_hPa=pressure, precip_mm=precipitation,windspeed_ms=wind_speed,
         winddir_degrees=wind_direction, humidity_percent=humidity, globalrad_units=global_radiation) %>% 
  mutate(station_id="DA_Mahmudia", 
         datetime_utc=as.POSIXct(datetime_utc, tz="utc", format="%Y-%m-%d %H:%M:%S"),
         winddir_degrees=ifelse(winddir_degrees==360,0,winddir_degrees)) %>% 
  select(station_id, datetime_utc, temp_c,Patm_hPa, precip_mm, humidity_percent, windspeed_ms,winddir_degrees,
         globalrad_units)


da_Tulcea<- read.csv(paste0(meteopath, "/Meteo_DA/Tulcea.csv")) %>% 
  rename(datetime_utc=time, temp_c=air_temp, Patm_hPa=pressure, precip_mm=precipitation,windspeed_ms=wind_speed,
         winddir_degrees=wind_direction, humidity_percent=humidity, globalrad_units=global_radiation) %>% 
  mutate(station_id="DA_Tulcea", 
         datetime_utc=as.POSIXct(datetime_utc, tz="utc", format="%Y-%m-%d %H:%M:%S"),
         winddir_degrees=ifelse(winddir_degrees==360,0,winddir_degrees)) %>% 
  select(station_id, datetime_utc, temp_c,Patm_hPa, precip_mm, humidity_percent, windspeed_ms,winddir_degrees,
         globalrad_units)
  
  



##DU####
#DUTCH Delta
#Stations: Two valid, same file
#Span: 2023-09-01 to 2024-10-01
#Frequency: hourly
#Timezone: UNKNOWN
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
         globalrad_joulepersquarecm=Q,
         windspeed_ms=FH/10, 
         winddir_degrees=case_when(DD%in%c(0,990)~NA_real_, #set 0=calm 990=variable to NA
                                   DD==360~0, #set North from 360 to 0
                                   TRUE~DD),
         cloudcover_percent=ifelse(N==9, NA_real_, N/8*100)) %>% #Set 9 (sky not vissible fog/heavysnow) to NA, transform to % 
  select(station_id, datetime_utc, temp_c,Patm_hPa, precip_mm, humidity_percent, windspeed_ms,winddir_degrees,
         globalrad_joulepersquarecm, cloudcover_percent)
  


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
#Timezone: UNKNOWN
#Radiation units: "Rsolar"  UNKNOWN

#RIA DE AVEIRO
#Time resolution: 10minutes

#Variables: 
#"TIME"     "Prec_Tot" "TA02"     "TA10"     "TA20"     "TA30"     "HR02"     "HR10"     "HR20"     "HR30"    "V10m"     "V10m_Max" "V10m_Dir" "V30m"     "V30m_Max" "V30m_Dir" "Rsolar"   "Press"   

#We will use the variables at 10m, 
#Assumed units: 
#TIME : assumed utc 
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
         datetime_utc=TIME,
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
         winddir_degrees=ifelse(V10m_Dir==360, 0, V10m_Dir),
         globalrad_units=Rsolar) %>% # substitute 360 N for 0 N
         select(station_id, datetime_utc, temp_c,Patm_hPa, precip_mm, humidity_percent, windspeed_ms,winddir_degrees,
                globalrad_units)

str(ri_EMA)

####TO-DO------
#Produce another dataset for RI with hourly resolution (to have an homogeneous dataset across case-pilots)
# ri_EMAh<- ri %>% 
  #Generate truncated datetime_utc only hour resolution, group by truncated_datetime, summarise (mean, na.rm=T), EXCEPTIONS: precip_mm: subsitute NAs with 0 and sum per hour, winddir_degrees: special function to search trigonometric mean (rad transformation so that mean of 0.1 and 359.9 produces 0), Globalrad_units: Investigate first how the global radiation is reported, to see if it makes sense to do mean or sum (if is cumulative radiation E/area or energyflux/area). 
  
rm(ri,ri_selection)         
         




##VA####
#VALENCIA
#Stations: one
#Span: 2023-01-01 to 2024-07-31
#Frequency: hourly
#Timezone: UNKNOWN
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
         globalrad_Wattpersquarmeter=`Solar radiation`
         ) %>% 
  select(station_id, datetime_utc, temp_c,Patm_hPa,
         # precip_mm,  # MISSING PRECIPITATION
         humidity_percent, windspeed_ms,winddir_degrees, globalrad_Wattpersquarmeter)


rm(va_selection, units_va_CEAV)
  

#---3. Check tz and radiation units -----

#CAMARGUE
#Span: 2023-01-01 to 2024-10-28
#Frequency: hourly
#Timezone: "UTC" checked
#Radiation units: Global rad (Joule /cm2) (mean71, max371)

#CURONIAN
#Span: 2023-09-01 to 2024-09-30
#Frequency: hourly
#Timezone: "UTC" checked
#Radiation units: NA (data not available)

#DANUBE
#Stations: Two valid
#Span: 2023-09-01 to 2024-09-30
#Frequency: hourly
#Timezone: UNKNOWN
#Radiation units: "global_radiation" UNKNOWN (mean614, max3470) 
#EMAIL SENT TO CONSTANTIN 20250320


#DUTCH Delta
#Stations: Two valid, same file
#Span: 2023-09-01 to 2024-10-01
#Frequency: hourly
#Timezone: UNKNOWN
#Radiation units: "Global radiation" (in J/cm2)  (mean45, max333)

#RIA DE AVEIRO
#Stations: one
#Span: 2023-09-01 to 2024-09-30
#Frequency: 10 minutes
#Timezone: UNKNOWN
#Radiation units: "Rsolar"  UNKNOWN (mean171, max1152) infered: W/squaremeter

#VALENCIA
#Stations: one
#Span: 2023-01-01 to 2024-07-31
#Frequency: hourly
#Timezone: UNKNOWN
#Radiation units: Solar radiation" W/squaremeter (mean 210, max1070)


#We need to check timezone for valencia, aveiro, dutch delta, danube use radiation data for sunrise and sunset, it shoudl be mostly ok (differences should be apparent in a year of data)





summary(ca_Arles)
summary(da_Tulcea)
summary(du_323)
summary(ri_EMA)
summary(va_CEAV)
names(ca_Arles)


da_Tulcea %>% 
  filter(globalrad_units >0) %>% 
  pull(globalrad_units ) %>% 
  hist(breaks = 30)


threshold<- 0.5

sunrise_sunset <- da_Mahmudia %>%
  mutate(
    date = as.Date(datetime_utc),
    hour = hour(datetime_utc),
    solar_radiation=globalrad_units
  ) %>% 
  group_by(date) %>%
  arrange(datetime_utc) %>%
  mutate(
    sunrise = if_else(lag(solar_radiation) <= threshold & solar_radiation > threshold, datetime_utc, NA),
    sunset = if_else(lag(solar_radiation) > threshold & solar_radiation <= threshold, datetime_utc, NA)
  ) %>%
  summarise(
    sunrise = min(sunrise, na.rm = TRUE),  # Find the first occurrence of sunrise for each day
    sunset = max(sunset, na.rm = TRUE)     # Find the last occurrence of sunset for each day
  )



