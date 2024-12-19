#GAIA_download_pictures

# ---
# Authors: Miguel Cabrera
# Project: "RESTORE4Cs"
# date: "Dec 2023"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
# This script downloads all pictures from GAIAgps naming them according a the unified table with data from all users  produced by Gaia_get_sampling_details script and modified to have a consistent naming of Season-Site-Subsite-Plot_picturenumber

#Input files: Harmonized table with sampling details from all Gaiagps users

#Output files: Pictures named in a consistent way across all sampling of Restore4Cs


rm(list = ls()) # clear workspace
cat("/014") # clear console


#---Packages----
library(jsonlite)
library(httr)
library(tidyverse)
library(purrr)


#---DIRECTORIES----

# Path to your GAIA folder (1 per user) where we will have GeoJSON files exported from GAIA (folders and potentially single waypoints). Here we will store the compiled table with all sampling information. 

#The directory where harmonized table with samplin details is located. Inside, a folder will be created with all pictures.
gaia_path<- "C:/Users/Miguel/Documents/GHG_chamber_GAIA_allusers/"


# Create a folder to save images if it doesn't exist
if (!dir.exists(paste0(gaia_path,"photos"))) {
  dir.create(paste0(gaia_path,"photos"))
}

downloaded_images_path<- paste0(gaia_path, "photos/")

#Load harmonized Sampling details #####
# EDIT when we have all gaia details

waypoints_data<- 
  
  





# Photo Download #####
#EDIT naming convention of pictures based on harmonized table with gaia details


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
        message(paste("Check that Public sharing is enabled in GAIAgps.com for chamber", waypoints_data$title[i]))
      })
    }
  } 
  if ((is.null(photo_url)|is.na(photo_url))){message(paste("No picture URL for waypoint", photo_title))}
}