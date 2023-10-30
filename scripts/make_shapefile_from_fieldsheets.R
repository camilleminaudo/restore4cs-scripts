
# ---
# Authors: Camille Minaudo
# Project: "RESTORE4Cs"
# date: "Oct 2023"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
# This script loads all the fieldsheets contained in the Dropbox folder, and
# make a shapefile out of it.


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
library(sp)
library(sf)
library(readxl)


source(paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/get_unix_times.R"))


# ---- Directories ----
dropbox_root <- "C:/Users/Camille Minaudo/Dropbox/RESTORE4Cs - Fieldwork/Data"


# list files in Dropbox
f <- list.files(dropbox_root, pattern = "ieldsheet", all.files = T, full.names = T, recursive = T)
i <- grep(pattern = ".xlsx", x = f)
f <- f[i]
r <- grep(pattern = "template",x=f)
f <- f[-r]
r <- grep(pattern = "Scans",x=f)
f <- f[-r]
r <- grep(pattern = "exetainer",x=f)
myfiles <- f[-r]

myfiles

# load each fieldsheet and extract relevant information
isF <- T
for (f in myfiles){

  # load file
  fieldsheet_temp <- readxl::read_xlsx(f,col_names = T)
  fieldsheet <- readxl::read_xlsx(f,skip = 2, col_names = F)
  names(fieldsheet) <- tolower(names(fieldsheet_temp))
  fieldsheet$date <- as.Date( fieldsheet$date, tryFormats = c("%d.%m.%y", "%d/%m/%y"))

  indSed <- grep(pattern = "Sediment", x=f)
  indGHG <- grep(pattern = "GHG", x=f)

  if(length(indSed) > 0){
    shp_data_temp <- data.frame(long = fieldsheet[[8]],
                                lat = fieldsheet[[9]],
                                variable = "sediment",
                                sampleID = fieldsheet$sample_id,
                                water_depth = fieldsheet[[10]],
                                strata = fieldsheet$strata,
                                date = fieldsheet$date,
                                unix_time_utc = get_unix_times(mydate = fieldsheet$date, mytime = fieldsheet[[7]]))
  } else if(length(indGHG) > 0){
    shp_data_temp <- data.frame(long = fieldsheet$longitude,
                                lat = fieldsheet$latitude,
                                variable = "chamber_measurement",
                                sampleID = paste(gsub(x = basename(f), pattern = "-Fieldsheet-GHG.xlsx", replacement = ""), "-plot",fieldsheet$plot_id,sep = ""),
                                filename = basename(f),
                                water_depth = fieldsheet$water_depth,
                                strata = fieldsheet$strata,
                                date =  fieldsheet$date,
                                unix_time_utc = get_unix_times(mydate = fieldsheet$date, mytime = fieldsheet$start_time))
  } else {
    message("variable name could not be recognized")
  }

  if (isF){
    isF <- F
    shp_data <- shp_data_temp
  } else {
    shp_data <- rbind(shp_data,shp_data_temp)
  }
}

dim(shp_data)



# ------------ Create a GIS file and save it ----------------


my_shp <- st_as_sf(shp_data[!is.na(shp_data$long),],
         coords = c("long", "lat"),
         crs = 4326)

setwd(paste0(dropbox_root,"/GIS"))
st_write(my_shp, dsn = "sampling_points_S1.shp", layer = "sampling_points_S1.shp", driver = "ESRI Shapefile")





library(maps)
world <- sf::st_as_sf(map('world', plot = FALSE, fill = TRUE))

ggplot(data = world) +
  geom_sf(fill= "antiquewhite") +
  geom_sf(data = my_shp, color = "red") +

  # coord_sf(crs = st_crs(3035))+
  coord_sf(xlim = c(-10, 44), ylim = c(24, 56), expand = FALSE)

