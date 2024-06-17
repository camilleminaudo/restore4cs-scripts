
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

source(paste0(dirname(dirname(rstudioapi::getSourceEditorContext()$path)),"/functions/get_unix_times.R"))



# ---- Settings ----
dropbox_root <- "C:/Users/Camille Minaudo/Dropbox/RESTORE4Cs - Fieldwork/Data/" #Camille Minaudo
# dropbox_root <- "C:/Users/misteli/Dropbox/RESTORE4Cs - Fieldwork/Data/" #Benjamin Misteli

#################################
sampling <- "S3"
#################################

# ---- Read fieldsheets ----

# list files in Dropbox
f <- list.files(dropbox_root, pattern = "Fieldsheet", all.files = T, full.names = T, recursive = T)
i <- grep(pattern = ".xlsx", x = f)
f <- f[i]
i <- grep(pattern = sampling, x = f) # selecting the files corresponding to the selected sampling campaign
f <- f[i]
# r <- grep(pattern = "template",x=f)
# f <- f[-r]
r <- grep(pattern = "Scans",x=f)
if(length(r)>0){
  f <- f[-r]  
}
r <- grep(pattern = "exetainer",x=f)
myfiles <- f[-r]

myfiles

# load each fieldsheet and extract relevant information
isF <- T
for (f in myfiles){

  # load file
  fieldsheet_temp <- readxl::read_xlsx(f,col_names = T, n_max = 1, .name_repair = "unique_quiet")
  fieldsheet <- readxl::read_xlsx(f,skip = 2, col_names = F, n_max = 100, .name_repair = "unique_quiet")

  if(dim(fieldsheet)[1]>1){
    names(fieldsheet) <- tolower(names(fieldsheet_temp))
    fieldsheet$date <- as.Date(fieldsheet[[3]], tryFormats = c("%d.%m.%Y", "%d/%m/%Y"))
    
    
    indSed <- grep(pattern = "Sediment", x=f)
    indGHG <- grep(pattern = "GHG", x=f)
    
    if(length(indSed) > 0){
      shp_data_temp <- data.frame(long = fieldsheet[[9]],
                                  lat = fieldsheet[[10]],
                                  variable = "sediment",
                                  variable_2 = fieldsheet[[6]],
                                  sampleID = fieldsheet[[5]],
                                  filename = basename(f),
                                  water_depth = fieldsheet[[11]],
                                  strata = fieldsheet$strata,
                                  date = fieldsheet$date,
                                  unix_time_utc = get_unix_times(mydate = fieldsheet$date, mytime = fieldsheet[[8]]))
    } else if(length(indGHG) > 0){
      shp_data_temp <- data.frame(long = fieldsheet[[11]],
                                  lat = fieldsheet[[12]],
                                  variable = "chamber_measurement",
                                  variable_2 = fieldsheet$strata,
                                  sampleID = paste(gsub(x = basename(f), pattern = "-Fieldsheet-GHG.xlsx", replacement = ""), "-plot",fieldsheet[[8]],sep = ""),
                                  filename = basename(f),
                                  water_depth = fieldsheet[[13]],
                                  strata = fieldsheet[[10]],
                                  date =  fieldsheet$date,
                                  unix_time_utc = get_unix_times(mydate = fieldsheet$date, mytime = fieldsheet[[15]]))
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
}



# load fieldsheet for water sampling and extract relevant information

setwd(dropbox_root)
filename_w <- paste0(dropbox_root,"Water/Water sampling and filtration_all data.xlsx")

# load file
fieldsheet_water <- readxl::read_xlsx(filename_w, 
                                      col_names = T, n_max = 3*6*6*4,
                                      sheet = "Water_sampling_master_ONLINE", skip = 8,
                                      col_types = c(rep("text",9),rep("numeric",21),rep("text",4)))

fieldsheet_water <- fieldsheet_water[fieldsheet_water$Survey==sampling,]

fieldsheet_water$date <- as.Date(as.numeric(fieldsheet_water$`Date dd/mm/yyyy`), origin = "1899-12-30")
fieldsheet_water$`Sampling time (Local)` <- as.numeric(fieldsheet_water$`Sampling time (Local)`)
myhour <- (floor(fieldsheet_water$`Sampling time (Local)`*24))
myminutes <- (round((fieldsheet_water$`Sampling time (Local)`*24 - floor(fieldsheet_water$`Sampling time (Local)`*24))*60))


shp_data_water <- data.frame(long = fieldsheet_water$`Longitude X °E (decimal)`,
                             lat = fieldsheet_water$`Latitude Y °N (decimal)`,
                             variable = "water",
                             variable_2 = "water",
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

myfilename <- paste0(sampling,"_sampling_points")


if(dir.exists(myfilename)){
  require(fs)
  dir_delete(myfilename)
}


# writing shp file to GIS directory
st_write(my_shp, dsn = myfilename, layer = myfilename, driver = "ESRI Shapefile", append = F)





library(maps)
world <- sf::st_as_sf(maps::map('world', plot = FALSE, fill = TRUE))

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
message(paste(dim(shp_data[shp_data$variable_2=="sediment sample",])[1],"sediment samples"))
message(paste(dim(shp_data[shp_data$variable_2=="core",])[1],"cores incubated"))


ggplot(shp_data, aes(variable))+geom_bar()+theme_bw()




