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

#The directory where harmonized table with sampling details is located. Inside, a folder will be created with all pictures.
gaia_path<- "C:/Users/Miguel/Documents/GHG_chamber_GAIA_tableANDpictures_allusers/"


# Create a folder to save images if it doesn't exist
if (!dir.exists(paste0(gaia_path,"photos"))) {
  dir.create(paste0(gaia_path,"photos"))
}

downloaded_images_path<- paste0(gaia_path, "photos/")

#Load harmonized Sampling details #####

waypoints_data<- read.csv(paste0(gaia_path, "Pictures_urls.csv")) %>% 
  mutate(plotcode_notext=gsub("plot","",plotcode)) %>% 
  arrange(time_created) %>% 
  group_by(plotcode_notext) %>% 
  mutate(uniquetitle=gsub("-","_",paste(plotcode_notext,row_number(),sep="_"))) %>% #Modify uniquetitle for filenames and add subsrits "_1", "_2" for duplicated plotcodes
  ungroup() %>% 
  select(-plotcode_notext) %>% 
  mutate(downloaded=NA,
         download_error=NA)
  
#Check uniqueness of uniquetitle
waypoints_data%>%
  filter(duplicated(uniquetitle) | duplicated(uniquetitle, fromLast = TRUE))

#view duplicated plotcodes and corresponding uniquetitle  
duplicated<-waypoints_data %>%
  filter(duplicated(plotcode) | duplicated(plotcode, fromLast = TRUE)) %>% 
  select(plotcode, uniquetitle)

head(duplicated)


# Photo Download #####
#EDIT naming convention of pictures based on harmonized table with gaia details:DONE


#EDIT loop to get data.frame with unaccessible pictures (title plotcode, gaia_file, gaia_user) to be able to request granting of access to users. We can add it to the waypoints_data, to have also a list of all the pictures downloaded and the ones with issues. Adding variable downloaded (T/F), and errormessage ("text" originating from error in loop)

# Download each full-size photo and save with the plotcode as part of the file name
for (i in 1:nrow(waypoints_data)) {
  photo_url <- waypoints_data$fullsize_url[i]
  photo_title <- waypoints_data$uniquetitle[i]
  
  # IF photo_url exists, try to download it and print Success/Error
  if (!is.null(photo_url) & !is.na(photo_url) & nzchar(photo_url)) {  # Check if the URL is not empty
    # Use the title in the file name, sanitize to remove unsafe characters
    image_name <- paste0(gsub("[^[:alnum:]_]", "_", photo_title), ".jpg")
    image_path <- file.path(paste0(downloaded_images_path, image_name))
    
    # Check if the file already exists
    if (file.exists(image_path)) {
      # Skip to next iteration if file already exists
      # message(paste("Picture", image_name, "already exists"))
      waypoints_data$downloaded[i] <- TRUE
      next()  # Skip to the next iteration of the loop
    }
    
    # Make a GET request to fetch the image
    response <- GET(photo_url)
    
    # Check if the response status is OK (status 200)
    if (response$status_code == 200) {
      # Check if the content length (size) is greater than zero
      content_length <- length(response$content)
      
      if (content_length > 0) {
        # Save the image to the file path if it's not empty
        writeBin(response$content, image_path)
        message(paste("Downloaded:", image_name))
        waypoints_data$downloaded[i] <- TRUE  # Set to TRUE since the download was successful
        waypoints_data$download_error[i] <- NA  # Clear any previous errors if download was successful
      } else {
        # Content is empty, don't create the file, mark as failed
        message(paste("Error: Downloaded content is empty for", photo_url))
        waypoints_data$downloaded[i] <- FALSE
        waypoints_data$download_error[i] <- "Downloaded content is empty"
      }
    } else {
      # If status code is not 200 (failure)
      message(paste("Error: Failed to download image from", photo_url))
      waypoints_data$downloaded[i] <- FALSE
      waypoints_data$download_error[i] <- paste("Failed with status code", response$status_code)
    }
  } else {
    # If photo_url is NULL or empty
    message(paste("No valid URL for waypoint", photo_title))
    waypoints_data$downloaded[i] <- FALSE
    waypoints_data$download_error[i] <- "No valid URL"
  }
}

#Write csv with sucess/failure of download process
write.csv(waypoints_data, file = paste0(gaia_path, "Pictures_urls_download_status.csv"),row.names = F)

#Check which gaia_folder and users failed the download
waypoints_data %>% 
  filter(downloaded!=T) %>% 
  select(uniquetitle, gaia_user, gaia_file, download_error) %>% 
  group_by(gaia_user, gaia_file, download_error) %>% 
  summarise(failed_downloads=n())






# Below, old version of loop (does not work properly)

# 
# for (i in 1:nrow(waypoints_data)) {
#   photo_url <- waypoints_data$fullsize_url[i]
#   photo_title <- waypoints_data$uniquetitle[i]
#   
#   # IF photo_url exists, try to download it and print Success/Error
#   if (!is.null(photo_url) & !is.na(photo_url)) {
#     if (nzchar(photo_url)) {  # Check if the URL is not empty
#       # Use the title in the file name, sanitize to remove unsafe characters
#       image_name <- paste0(gsub("[^[:alnum:]_]", "_", photo_title), ".jpg")
#       image_path <- file.path(paste0(downloaded_images_path, image_name))
#       
#       # Check if the file already exists
#       if (file.exists(image_path)) {
#         # Skip to next iteration if file already exists
#         message(paste("Picture", image_name, "already exists"))
#         waypoints_data$downloaded[i] <- TRUE
#         next()  # Skip to the next iteration of the loop
#       }
#       
#       # Try to download the image
#       tryCatch({
#         # Make a GET request to fetch the image
#         response <- GET(photo_url)
#         
#         # Check if the response status is OK (status 200)
#         if (response$status_code == 200) {
#           # Check if the content length (size) is greater than zero
#           content_length <- length(response$content)
#           if (content_length > 0) {
#             # Save the image to the file path if it's not empty
#             writeBin(response$content, image_path)
#             message(paste("Downloaded:", image_name))
#             waypoints_data$downloaded[i] <- TRUE
#             waypoints_data$download_error[i] <- NA  # Clear any previous errors if download was successful
#           } else {
#             stop("Downloaded content is empty.")
#           }
#         } else {
#           stop("Failed to download the image or server error.")
#         }
#       }, error = function(e) {
#         # Handle the error: log the message and update dataframe
#         message(paste("Error downloading", photo_url, ":", e$message))
#         waypoints_data$downloaded[i] <- FALSE
#         waypoints_data$download_error[i] <- e$message  # Store the error message in the dataframe
#         # Make sure the empty file is not created in case of failure
#         if (file.exists(image_path)) {
#           file.remove(image_path)  # Remove the empty file if it exists
#         }
#         message(paste("Check that Public sharing is enabled in GAIAgps.com for chamber", waypoints_data$title[i], "corresponding to plotcode", photo_title))
#       })
#     }
#   } else {
#     message(paste("No picture URL for waypoint", photo_title))
#     waypoints_data$downloaded[i] <- FALSE  # No URL means failed download
#     waypoints_data$download_error[i] <- "No URL provided"
#   }
# }




