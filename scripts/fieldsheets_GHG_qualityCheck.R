# ---
# Authors: Camille Minaudo
# Project: "RESTORE4Cs"
# date: "Oct 2023"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
# This script loads all GHG chamber measurement fieldsheets and runs a simple quality check:
# ==> makes sure that end time is always > start time
#


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
require(dplyr)
require(purrr)


source(paste0(dirname(dirname(rstudioapi::getSourceEditorContext()$path)),"/functions/get_unix_times.R"))
source(paste0(dirname(dirname(rstudioapi::getSourceEditorContext()$path)),"/functions/read_GHG_fieldsheets.R"))



# ---- Directories ----
dropbox_root <- "C:/Users/Camille Minaudo/Dropbox/RESTORE4Cs - Fieldwork/Data"
fieldsheetpath <- paste0(dropbox_root,"/GHG/Fieldsheets")

# ---- List GHG chamber fieldsheets in Dropbox and read them ---
# list filenames
myfieldsheets_list <- list.files(fieldsheetpath, pattern = "Fieldsheet-GHG.xlsx", all.files = T, full.names = T, recursive = T)

#################################
# sampling <- "S3"
# #################################
# 
# i <- grep(pattern = sampling, x = myfieldsheets_list) # selecting the files corresponding to the selected sampling campaign
# myfieldsheets_list <- myfieldsheets_list[i]

# Read all fieldsheets and put them in a single dataframe
fieldsheet <- read_GHG_fieldsheets(myfieldsheets_list)

# --------- rows with stop time < start time
ind_erronous_times <- which(fieldsheet$unix_stop < fieldsheet$unix_start)

message("the following rows show erronous start/stop times")
as.data.frame(fieldsheet[ind_erronous_times,c("pilot_site","subsite","plot_id","start_time","end_time")])


# rows with duration suspicious
ind_duration_neg <- which(fieldsheet$unix_stop-fieldsheet$unix_start < 0)
ind_duration_tooshort <- which(fieldsheet$unix_stop-fieldsheet$unix_start < 3*60)
ind_duration_toolong <- which(fieldsheet$unix_stop-fieldsheet$unix_start > 15*60)

ind_duration <- c(ind_duration_neg, ind_duration_tooshort, ind_duration_toolong)

fieldsheet$duration <- fieldsheet$unix_stop-fieldsheet$unix_start

message("the following rows show suspicious start/stop times")
as.data.frame(fieldsheet[ind_duration,c("pilot_site","subsite","plot_id","start_time","end_time","duration")])




# --------- rows with possible error with co2 values

ind_suspicious_co2_initial <- which(fieldsheet$initial_co2 > 1000)
ind_suspicious_co2_final <- which(fieldsheet$final_co2 > 1000)
ind_suspicious_co2 <- unique(c(ind_suspicious_co2_initial, ind_suspicious_co2_final))

message("the following rows show suspiciously HIGH CO2 levels")
as.data.frame(fieldsheet[ind_suspicious_co2,c("pilot_site","subsite","plot_id","initial_co2","final_co2")])


ind_suspicious_co2_initial <- which(fieldsheet$initial_co2 < 200)
ind_suspicious_co2_final <- which(fieldsheet$final_co2 < 200)
ind_suspicious_co2 <- unique(c(ind_suspicious_co2_initial, ind_suspicious_co2_final))

message("the following rows show suspiciously LOW CO2 levels")
as.data.frame(fieldsheet[ind_suspicious_co2,c("pilot_site","subsite","plot_id","initial_co2","final_co2")])


# --------- rows with possible error with CH4 units (ppb instead of ppm)

ind_suspicious_ch4_initial <- which(fieldsheet$initial_ch4 > 2000)
ind_suspicious_ch4_final <- which(fieldsheet$final_ch4 > 2000)
ind_suspicious_ch4 <- unique(c(ind_suspicious_ch4_initial, ind_suspicious_ch4_final))

message("the following rows show suspiciously HIGH methane levels")
as.data.frame(fieldsheet[ind_suspicious_ch4,c("pilot_site","subsite","plot_id","initial_ch4","final_ch4")])


ind_suspicious_ch4_initial <- which(fieldsheet$initial_ch4 < 1)
ind_suspicious_ch4_final <- which(fieldsheet$final_ch4 < 1)
ind_suspicious_ch4 <- unique(c(ind_suspicious_ch4_initial, ind_suspicious_ch4_final))

message("the following rows show suspiciously LOW methane levels")
as.data.frame(fieldsheet[ind_suspicious_ch4,c("pilot_site","subsite","plot_id","initial_ch4","final_ch4")])



# --------- rows with depth possibly reported in m instead of cm
fieldsheet_not_dry <- fieldsheet[fieldsheet$water_depth != 0,]
ind_suspicious_depth <- which(fieldsheet_not_dry$water_depth < 2)

message("the following rows show suspicious water depth ranges:")
as.data.frame(fieldsheet_not_dry[ind_suspicious_depth,c("pilot_site","subsite","plot_id","water_depth")])
message("which corresponds to the following subsites:")
unique(fieldsheet_not_dry$subsite[ind_suspicious_depth])



# --------- missing logger data info
ind_loggers <- which(is.na(fieldsheet$logger_floating_chamber) | is.na(fieldsheet$logger_transparent_chamber))

message("the following rows do not have logger data info")
unique(fieldsheet$subsite[ind_loggers])


# --------- Possible mistake with chamber type or chamber_height
message(" Possible mistake with chamber type or chamber_height")

ind_suspicious_chamb_type <- which(fieldsheet$chamber_type=="tube"& fieldsheet$chamber_height_cm==0 & fieldsheet$strata=="open water")
as.data.frame(fieldsheet[ind_suspicious_chamb_type,c("pilot_site","subsite","plot_id","chamber_type","strata","chamber_height_cm")])




# --------- Possible mistake with chamber type or chamber_height

my_shp <- st_as_sf(fieldsheet,
                   coords = c("longitude", "latitude"),
                   crs = 4326)
setwd(paste0(dropbox_root,"/GIS"))
# writing shp file to GIS directory
# myfilename <- paste0("GHG_sampling_points")
# st_write(my_shp, dsn = myfilename, layer = myfilename, driver = "ESRI Shapefile", append = F)

# load subsites polygons
setwd(dropbox_root)
shp_sites <- st_read("./GIS/Shape_files_subsites/Rough shapes for export/Subsites_rough_forExport.shp")

# intersect my_shp with shp_sites and find points outside polygons
# find corresponding basins
polygs = st_as_sf(shp_sites)
points = st_as_sf(my_shp)
sf_use_s2(FALSE)

intrsct = st_intersection(points, polygs) # intersect polygons with points, keeping the information from both

ind_match <- match(x = points$uniqID, table = intrsct$uniqID)

message("Wrong coordinates")
points$uniqID[is.na(ind_match)]




