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
library(openxlsx)
library(chron)

source(paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/get_unix_times.R"))

# ---- Settings ----
#dropbox_root <- "C:/Users/Camille Minaudo/Dropbox/RESTORE4Cs - Fieldwork/Data/" #Camille Minaudo
dropbox_root <- "D:/Dropbox/RESTORE4Cs - Fieldwork/Data/Sediment/" #Benjamin Misteli
sampling <- "S1"
setwd(dropbox_root)

# list files in Dropbox
f <- list.files(dropbox_root, pattern = "Fieldsheet", all.files = T, full.names = T, recursive = T)
i <- grep(pattern = ".xlsx", x = f)
f <- f[i]

myfiles <- f
fieldsheets_list <- list()

for (f in myfiles){
  
  # load file
  fieldsheet_temp <- readxl::read_xlsx(f,col_names = T)
  fieldsheet <- readxl::read_xlsx(f,skip = 2, col_names = F, n_max = 100)
  names(fieldsheet) <- tolower(names(fieldsheet_temp))
  fieldsheet$date <- as.Date(fieldsheet$date, tryFormats = c("%d.%m.%Y", "%d/%m/%Y"))
  fieldsheets_list[[f]] <- fieldsheet
}
complete_list <- do.call(rbind, fieldsheets_list)
row.names(complete_list) <- NULL
complete_list_sediments <- complete_list %>%
  filter(type == "sediment sample")

datasheet_wcl_names <- readxl::read_xlsx("Lab data/Datasheet-Sediment.xlsx",col_names = T)
datasheet_wcl <- readxl::read_xlsx("Lab data/Datasheet-Sediment.xlsx",skip = 2,col_names = F, n_max = 100)
names(datasheet_wcl) <- tolower(names(datasheet_wcl_names))

merged_data <- merge(complete_list_sediments, datasheet_wcl, by = "sample_id", all = TRUE)
merged_data <- merged_data %>%
  arrange(desc(package_arrived))

write.xlsx(merged_data, file = "Lab data/test.xlsx")

