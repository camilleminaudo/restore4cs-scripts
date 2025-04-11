#Get_gaiaGPS_sampling_details

# ---
# Authors: Miguel Cabrera
# Project: "RESTORE4Cs"
# date: "Dec 2023"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
# This script uses exported GEOjson files to compile all the sampling details (GaiaID, GaiaTitle, coordinates, time, notes) and associated picture-URLs in a single table per gaia_user.

#Input files: exported GEOjson files from GAIA (can be GAIAfolders, GAIAwaypoints or even GAIAtracks).

#Output files: Table with GAIA sampling details and picture-URLs per user


#IMPORTANT NOTES: 
#Waypoints and folders have to be set to Public Sharing==ON to download the pictures. 

#GAIA Subfolders are lost when exporting a master folder as GEOjson: Master>Subfolder>Waypoint in GAIA becomes Master>Waypoint in GEOJson file. Adapt your Gaia data organization accordingly when possible. (if its too cumbersome to rename all the waypoints, we will retrieve the season-site-subsite part of the plotcode with date and coordinates).


rm(list = ls()) # clear workspace
cat("/014") # clear console


#---Packages----
library(jsonlite)
library(httr)
library(tidyverse)
library(purrr)


#---DIRECTORIES----

# Path to your GAIA folder (1 per user) where we will have GeoJSON files exported from GAIA (folders and potentially single waypoints). Here we will store the compiled table with all sampling information from each user. 

#A name for us to know who took the pictures & coordinates with GAIA
gaia_user<- "Camille"

#dropbox path
dropbox_path<- "C:/Users/Miguel/Dropbox/RESTORE4Cs - Fieldwork/Data/"

#GaiaGPS sampling path 
gaia_path<- paste0(dropbox_path,"GHG/Gaia_GPS/")

#Path to exported gaia files (GEOjson), ONE PER USER
geojson_path<- paste0(gaia_path,gaia_user, "_exported_geojson_files/")

{  
#List all geojson files exported from GAIA
geojson_files<- list.files(path = geojson_path, pattern = ".geojson",full.names = F)


#LOOP over GAIA files####
for (geo in geojson_files){
  
  message(paste("Retrieving ",geo," sampling info and picture URLs" ))
  
  # Read the GeoJSON file
  geojson_data <- fromJSON(paste0(geojson_path, geo))
  
  #2 types of files: Point/Feature and FeatureCollection
  #2 options for pictures (they exsist or not)
  
##Folders ######
#IF file is of type Folder: 
if(geojson_data$type=="FeatureCollection"){

#Import identifiers (title, time, lat,long,notes) from features$properties 
waypoints_data<- geojson_data$features$properties %>% 
  select(id,title, time_created,latitude,longitude, notes)%>% 
  rename(waypoint_id=id)


# Extract id and fullsize_url from geojson file, handle missing entries (for waypoints without picture)
photo_data <- map_dfr(geojson_data$features$properties$photos, ~{
  # Check if the data frame exists and has title and fullsize_url columns (i.e. if the waypoint has a picture attached)
  if (is.null(.) || !all(c("waypoint_id", "fullsize_url") %in% colnames(.))) {
    tibble(waypoint_id = NA, fullsize_url = NA)  # If missing, return a row with NAs
  } else {
    select(., waypoint_id, fullsize_url)  # Otherwise, extract the columns
  }
})

#for tracks only: to avoid issues of script if they are included
if(is.null(geojson_data$features$properties$photos) ){
  photo_data<-tibble(waypoint_id = NA, fullsize_url = NA)  # If missing, return a row with NAs
}
}

#Single waypoints####
#IF file geojson is single waypoint: 
if(geojson_data$type%in%c("Point","Feature")){

#Get waypoint details in any case 
waypoints_data<-   as.data.frame(geojson_data$properties[c("id","title", "time_created","latitude","longitude", "notes")]) %>% rename(waypoint_id=id)

# geojson_data$properties$photos is an empty list for waypoints without pictures, if that happens create empty photo_data dataframe
if(is.list(geojson_data$properties$photos)){
  photo_data<- data.frame(waypoint_id=NA, fullsize_url=NA)
}

# geojson_data$properties$photos is a data frame for waypoints with pictures,if that happens, simply select id and fullsize_url
if(is.data.frame(geojson_data$properties$photos)){
  photo_data <-geojson_data$properties$photos %>% 
    select(waypoint_id, fullsize_url)
}

}


#Merge identifiers and photo URLs (keeping mutliple URLs if they exist)
waypoints_data<- waypoints_data %>% 
  full_join(photo_data, by="waypoint_id",multiple = "all") %>% 
  filter(!is.na(waypoint_id)) %>% #Drop rows without data (created by loop for operational reasons)
  mutate(gaia_file=geo,
         gaia_user=gaia_user)

rm(photo_data)


#Compile table####
#Add file identifier and add file data to compilation of GAIA table
if (geo==geojson_files[1]){
  all_waypoints<- waypoints_data
} else { all_waypoints<- rbind(all_waypoints,waypoints_data)}

}

all_waypoints<- all_waypoints %>% arrange(time_created)  

#SAVE TABLE####
#Save GAIA data (expandoing any previous data)
write.table(
  all_waypoints, 
  file = paste0(gaia_path, "GAIA_table_", gaia_user, ".csv"), 
  sep = ",", 
  row.names = FALSE, 
  col.names = !file.exists(paste0(gaia_path, "GAIA_table_", gaia_user, ".csv")), 
  append = F
)
}

