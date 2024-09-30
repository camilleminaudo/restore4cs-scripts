# ---
# Authors: Benjamin Misteli
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

source(paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/get_unix_times.R"))

# ---- Settings ----
#dropbox_root <- "C:/Users/Camille Minaudo/Dropbox/RESTORE4Cs - Fieldwork/Data/" #Camille Minaudo
#dropbox_root <- "D:/Dropbox/RESTORE4Cs - Fieldwork/Data/Sediment/" #Benjamin Misteli PC
dropbox_root <- "C:/Users/misteli/Dropbox/RESTORE4Cs - Fieldwork/Data/Sediment" #Benjamin Misteli Laptop
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

write.xlsx(complete_list_sediments, file = "Lab data/fieldsheet_all_summary.xlsx")

fieldsheet_all_summary <- readxl::read_xlsx("Lab data/fieldsheet_all_summary.xlsx",col_names = T)
data_afdm <- readxl::read_xlsx("Lab data/LabData_dry_afdm.xlsx",col_names = T)
data_afdm_fieldsheet <- merge(fieldsheet_all_summary, data_afdm, by = "sample_id", all.x = TRUE)

data_cn <- readxl::read_xlsx("Lab data/LabData_cn.xlsx",col_names = T)
data_afdm_cn_fieldsheet <- merge(data_afdm_fieldsheet, data_cn, by = "bowl_id", all.x = TRUE)

data_afdm_cn_fieldsheet <- data_afdm_cn_fieldsheet %>%
  select(-c(1,3:8, 14:23, 25:27)) %>%
  rename(comments_fieldsheet = comments)%>%
  separate(sample_id, into = c("Season", "Site", "Subsite", "Replicate"), sep = "-", remove = FALSE)

write.xlsx(data_afdm_cn_fieldsheet, file = "Lab data/RESTORE4Cs_Sediment_Data_complete.xlsx")