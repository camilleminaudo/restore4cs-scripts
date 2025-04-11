#GAIA_export

# ---
# Authors: Miguel Cabrera
# Project: "RESTORE4Cs"
# date: "Dec 2023"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
# This script uses exported GEOjson files to download all pictures associated with the chambers and stores the sampling details (plotname, coordinates, time, notes) in a single table.

#Input files: exported GEOjson files from GAIA (can be folders, waypoints or even tracks).

#Output files: Named pictures of chambers and table with sampling details 


#Some notes: 
#GAIA Subfolders are lost when exported to GEOjson: Master>Subfolder>Waypoint in GAIA becomes Master>Waypoint in GEOJson file

#Waypoints and folders have to be set to Public Sharing==ON to download the pictures. 


rm(list = ls()) # clear workspace
cat("/014") # clear console


#---Packages----
library(jsonlite)
library(httr)
library(tidyverse)
library(purrr)


#---DIRECTORIES----

# Path to your GAIA folder (1 per user) where we will have GeoJSON files exported from GAIA (folders and potentially single waypoints). Here we will store the compiled table with all sampling information. 

gaia_user<- "Miguel"

gaia_path<- "C:/Users/Miguel/Documents/GHG_chamber_GAIA_Miguel/"



# Create a folder to save images if it doesn't exist
  if (!dir.exists(paste0(gaia_path,"photos"))) {
    dir.create(paste0(gaia_path,"photos"))
  }

#DO not modify!
geojson_path<- paste0(gaia_path, "exported_geojson_files/")
downloaded_images_path<- paste0(gaia_path, "photos/")
  
  
  
  
  
#List all geojson files exported from GAIA
geojson_files<- list.files(path = geojson_path, pattern = ".geojson",full.names = F)



#LOOP over GAIA files####
for (geo in geojson_files){
  
  message(paste("Retrieving ",geo," sampling info and pictures" ))
  
  # Read the GeoJSON file
  geojson_data <- fromJSON(paste0(geojson_path, geo))
  
  #2 types of files: Point/Feature and FeatureCollection
  #2 options for pictures (they exsist or not)
  
##Folders ######
#IF file is of type Folder: 
if(geojson_data$type=="FeatureCollection"){

#Import identifiers (title, time, lat,long,notes) from features$properties 
waypoints_data<- geojson_data$features$properties %>% 
  select(title, time_created,latitude,longitude, notes)


# Extract title and fullsize_url from geojson file, handle missing entries (for waypoints without picture)
photo_data <- map_dfr(geojson_data$features$properties$photos, ~{
  # Check if the data frame exists and has title and fullsize_url columns (i.e. if the waypoint has a picture attached)
  if (is.null(.) || !all(c("title", "fullsize_url") %in% colnames(.))) {
    tibble(title = NA, fullsize_url = NA)  # If missing, return a row with NAs
  } else {
    select(., title, fullsize_url)  # Otherwise, extract the columns
  }
})

#for tracks only: to avoid issues of script if they are included
if(is.null(geojson_data$features$properties$photos) ){
  photo_data<-tibble(title = NA, fullsize_url = NA)  # If missing, return a row with NAs
}
}

#Single waypoints####
#IF file geojson is single waypoint: 
if(geojson_data$type%in%c("Point","Feature")){

#Get waypoint details in any case 
waypoints_data<-   as.data.frame(geojson_data$properties[c("title", "time_created","latitude","longitude", "notes")])

# geojson_data$properties$photos is an empty list for waypoints without pictures, if that happens create empty photo_data dataframe
if(is.list(geojson_data$properties$photos)){
  photo_data<- data.frame(title=NA, fullsize_url=NA)
}

# geojson_data$properties$photos is a data frame for waypoints with pictures,if that happens, simply select title and fullsize_url
if(is.data.frame(geojson_data$properties$photos)){
  photo_data <-geojson_data$properties$photos %>% 
    select(title, fullsize_url)
}

}


#Merge identifiers and photo URLs (keeping mutliple URLs if they exist)
waypoints_data<- waypoints_data %>% 
  full_join(photo_data, by="title",multiple = "all") %>% 
  filter(!is.na(title)) %>% #Drop rows without data
  group_by(title) %>% 
  mutate(uniquetitle=paste(title, row_number(),sep = "_")) %>% 
  ungroup() 

rm(photo_data)

# Photo Download #####
# Download each full-size photo and save with the title as part of the file name
for (i in 1:nrow(waypoints_data)) {
  photo_url <- waypoints_data$fullsize_url[i]
  photo_title <- waypoints_data$uniquetitle[i]
  
  
  
  #IF photo_url exist, try to download it and print Success/Error 
  if (!is.null(photo_url)&!is.na(photo_url)){
       if (nzchar(photo_url)) { # Check if the URL is not empty
    # Use the title in the file name, sanitize to remove unsafe characters
    image_name <- paste0(gsub("[^[:alnum:]_]", "_", photo_title), ".jpg")
    image_path <- file.path(paste0(downloaded_images_path, image_name))
    
    # Check if the file already exists
    if (file.exists(image_path)) {
      message(paste("Picture", image_name, "already exists"))
      next()  # Skip to the next iteration of the loop
    }
    
    # Download the image
    tryCatch({
      GET(photo_url, write_disk(image_path))
      message(paste("Downloaded:", image_name))
    }, error = function(e) {
      message(paste("Error downloading", photo_url, ":", e$message))
      message(paste("Check that Public sharing is enabled in GAIAgps.com for plot", waypoints_data$title[i]))
    })
  }
  } 
  if ((is.null(photo_url)|is.na(photo_url))){message(paste("No picture URL for waypoint", photo_title))}
}

#Compile table####
#Add file identifier and add file data to compilation of GAIA table
waypoints_data<- waypoints_data %>% 
  mutate(gaia_file=geo,
         gaia_user=gaia_user)

if (geo==geojson_files[1]){
  all_waypoints<- waypoints_data
} else { all_waypoints<- rbind(all_waypoints,waypoints_data)}

}


#SAVE TABLE####
#Save GAIA data (expandoing any previous data)
write.table(
  all_waypoints, 
  file = paste0(gaia_path, "GAIA_table_", gaia_user, ".csv"), 
  sep = ",", 
  row.names = FALSE, 
  col.names = !file.exists(paste0(gaia_path, "GAIA_table_", gaia_user, ".csv")), 
  append = TRUE
)
