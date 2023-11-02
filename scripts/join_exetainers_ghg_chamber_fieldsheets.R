
# ---
# Authors: Camille Minaudo
# Project: "RESTORE4Cs"
# date: "Oct 2023"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
# This script combine an exetainer fieldsheet with a GHG fieldsheet to provide
# an indication of the concentrations of CO2 and CH4 expected in the gas samples.
# This should be useful to prepare for Gas Chromatography analysis of the samples.


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


campaign <- "S1"
site <- "CU"
subsites <- c("P1","P2","A1","A2","R1","R2")

dropbox_root <- "C:/Users/Camille Minaudo/Dropbox/RESTORE4Cs - Fieldwork/Data/"

path2output <- "C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/data/RESTORE4Cs/GasChromato"


# ---- Functions ----


# Read corresponding Fieldsheet
read_fieldsheet <- function(path2file){
  fieldsheet_temp <- readxl::read_xlsx(path2file,
                                       col_names = T)
  fieldsheet <- readxl::read_xlsx(path2file,
                                  skip = 2, col_names = F)
  names(fieldsheet) <- names(fieldsheet_temp)

  abnormal_col_name <- grep("transparent or dark", x = names(fieldsheet))
  if(length(abnormal_col_name)>0){
    names(fieldsheet)[grep("transparent or dark", x = names(fieldsheet))] <- "transparent_dark"
  }

  # fieldsheet$unix_start_time <- get_unix_times(mydate = fieldsheet$date, mytime = fieldsheet$start_time)
  # fieldsheet$unix_end_time <- get_unix_times(mydate = fieldsheet$date, mytime = fieldsheet$end_time)
  return(fieldsheet)
}


join_fieldsheets <- function(exetainers_fieldsheetpath, ghg_fieldsheetpath){
  exetainers_fieldsheet <- read_fieldsheet(path2file = exetainers_fieldsheetpath)
  ghg_fieldsheet <- read_fieldsheet(path2file = ghg_fieldsheetpath)


  # provide a unique identification for each chamber incubation
  ghg_fieldsheet$unique_id <- paste("plot",ghg_fieldsheet$plot_id,ghg_fieldsheet$strata,ghg_fieldsheet$transparent_dark, sep = "_")

  exetainers_fieldsheet$unique_id <- paste("plot",exetainers_fieldsheet$PlotID,exetainers_fieldsheet$Strata,exetainers_fieldsheet$transparent_dark, sep = "_")

  exetainers_fieldsheet <- exetainers_fieldsheet %>%
    left_join(ghg_fieldsheet %>% select(unique_id, final_co2, final_ch4, comments))


  exetainers_fieldsheet$final_co2[exetainers_fieldsheet$`headspace/trapped/atm` == "atmosphere"] <- mean(ghg_fieldsheet$initial_co2[ghg_fieldsheet$chamber_type == "floating" |
                                                                                                                                      (ghg_fieldsheet$chamber_type == "tube" & ghg_fieldsheet$transparent_dark == "transparent")])

  exetainers_fieldsheet$final_ch4[exetainers_fieldsheet$`headspace/trapped/atm` == "atmosphere"] <- mean(ghg_fieldsheet$initial_ch4[ghg_fieldsheet$chamber_type == "floating" |
                                                                                                                                      (ghg_fieldsheet$chamber_type == "tube" & ghg_fieldsheet$transparent_dark == "transparent")])

  exetainers_fieldsheet$final_co2[exetainers_fieldsheet$`headspace/trapped/atm` == "headspace"] <- exetainers_fieldsheet$final_ch4[exetainers_fieldsheet$`headspace/trapped/atm` == "headspace"] <- NA

  df <- as.data.frame(exetainers_fieldsheet)


  df <- rename(df, expected_co2_ppm = final_co2, expected_ch4_ppm = final_ch4)

  df <- df[,c("unique_id","date","person_sampling","ExetainerID","headspace/trapped/atm","expected_co2_ppm","expected_ch4_ppm")]

  return(df)
}


# ---- Directories ----

isF <- T
for(subsite in subsites){
  message(paste0("Processing ",campaign,"-",site,"-",subsite))
  exetainers_fieldsheetpath <- paste0(dropbox_root, "GHG/Fieldsheets/",campaign,"-",site,"/",campaign,"-",site,"-",subsite,"-Fieldsheet-exetainers.xlsx")
  ghg_fieldsheetpath <- paste0(dropbox_root, "GHG/Fieldsheets/",campaign,"-",site,"/",campaign,"-",site,"-",subsite,"-Fieldsheet-GHG.xlsx")


  df <- join_fieldsheets(exetainers_fieldsheetpath, ghg_fieldsheetpath)



  setwd(path2output)
  myfilename <- paste0(campaign,"-",site,"-",subsite,"-exetainers_co2_ch4_concentrations.csv")
  write.csv(x = df, file = myfilename, sep = ";", dec = ".", row.names = F, col.names = T)


  if(isF){
    isF <- F
    df_all <- df
  } else {
    df_all <- rbind(df_all, df)
  }
}

setwd(path2output)
myfilename <- paste0(campaign,"-",site,"-all-exetainers_co2_ch4_concentrations.csv")
write.csv(x = df_all, file = myfilename, sep = ";", dec = ".", row.names = F, col.names = T)





