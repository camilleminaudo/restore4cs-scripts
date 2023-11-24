
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

source(paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/get_unix_times.R"))

# ---- Settings ----
dropbox_root <- "C:/Users/Camille Minaudo/Dropbox/RESTORE4Cs - Fieldwork/Data/"
sampling <- "S1"

# list files in Dropbox
f <- list.files(dropbox_root, pattern = "Fieldsheet", all.files = T, full.names = T, recursive = T)
i <- grep(pattern = ".xlsx", x = f)
f <- f[i]
i <- grep(pattern = sampling, x = f) # selecting the files corresponding to the selected sampling campaign
f <- f[i]
# r <- grep(pattern = "template",x=f)
# f <- f[-r]
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
                                sampleID = fieldsheet[[5]],
                                filename = basename(f),
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


setwd(dropbox_root)
filename_w <- paste0(dropbox_root,"Water/Water sampling and filtration_all data.xlsx")

# load file
fieldsheet_water <- readxl::read_xlsx(filename_w,
                                      col_names = T, n_max = 3*6*6*4,
                                      sheet = "Water_sampling_master_ONLINE", skip = 8,
                                      col_types = c(rep("text",7),rep("numeric",17),rep("text",4)))

fieldsheet_water$date <- as.Date(as.numeric(fieldsheet_water$`Date dd/mm/yyyy`), origin = "1899-12-30")
myhour <- (floor(fieldsheet_water$`Sampling time (Local)`*24))
myminutes <- (round((fieldsheet_water$`Sampling time (Local)`*24 - floor(fieldsheet_water$`Sampling time (Local)`*24))*60))


shp_data_water <- data.frame(long = fieldsheet_water$`Longitude X °E (decimal)`,
                             lat = fieldsheet_water$`Latitude Y °N (decimal)`,
                             variable = "water",
                             sampleID = fieldsheet_water$`Label Sample ID`,
                             filename = basename(filename_w),
                             water_depth = fieldsheet_water$`Water depth (m)`,
                             strata = "open water",
                             date = fieldsheet_water$date,
                             unix_time_utc = NA)



shp_data <- rbind(shp_data,shp_data_water)

shp_data <- shp_data[!is.na(shp_data$long),]
dim(shp_data)


# ------------ Create a GIS file and save it ----------------

my_shp <- st_as_sf(shp_data[!is.na(shp_data$long),],
                   coords = c("long", "lat"),
                   crs = 4326)


setwd(paste0(dropbox_root,"/GIS"))

myfilename <- "sampling_points_S1"


if(dir.exists(myfilename)){
  require(fs)
  dir_delete(myfilename)
}
#
# # removing files if this layer was already created before
# extensions <- c("shp","dbf","prj","shx")
# if(file.exists(paste0(myfilename,".shp"))){
#   for(extension in extensions){
#     file.remove(paste0(myfilename,extension))
#   }
# }

# writing shp file to GIS directory
st_write(my_shp, dsn = myfilename, layer = myfilename, driver = "ESRI Shapefile")





library(maps)
world <- sf::st_as_sf(map('world', plot = FALSE, fill = TRUE))

ggplot(data = world) +
  geom_sf(fill= "antiquewhite") +
  geom_sf(data = my_shp, color = "red") +

  # coord_sf(crs = st_crs(3035))+
  coord_sf(xlim = c(-10, 44), ylim = c(38, 56), expand = FALSE)

# shp_data$water_depth <- as.numeric(shp_data$water_depth)
# ggplot(shp_data[!is.na(shp_data$water_depth),])+geom_density(aes(water_depth))+theme_article()




# some statistics
message(paste(dim(shp_data[shp_data$variable=="chamber_measurement",])[1],"GHG chamber plots"))
message(paste(dim(shp_data[shp_data$variable=="water",])[1],"water samples"))
message(paste(dim(shp_data[shp_data$variable=="sediment",])[1],"sediment samples"))


ggplot(shp_data, aes(variable))+geom_bar()+theme_bw()




