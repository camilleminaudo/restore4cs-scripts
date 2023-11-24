
# ---
# Authors: Camille Minaudo
# Project: "RESTORE4Cs"
# date: "Oct 2023"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
# This script


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
library(GoFluxYourself)
require(dplyr)
require(purrr)

source(paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/get_unix_times.R"))
source(paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/read_GHG_fieldsheets.R"))



# ---- Directories ----
dropbox_root <- "C:/Users/Camille Minaudo/Dropbox/RESTORE4Cs - Fieldwork/Data"
datapath <- paste0(dropbox_root,"/GHG/RAW data")
fieldsheetpath <- paste0(dropbox_root,"/GHG/Fieldsheets")

# ---- List GHG chamber fieldsheets in Dropbox ----
myfieldsheets_list <- list.files(fieldsheetpath, pattern = "Fieldsheet-GHG.xlsx", all.files = T, full.names = T, recursive = T)


# ---- Read all fieldsheets and put them in a single dataframe ----

fieldsheet <- read_GHG_fieldsheets(myfieldsheets_list)


# ---- Go through each incubation in fieldsheet and make a plot, organized by subsite ----

for (subsite in unique(fieldsheet$subsite)){

  corresp_fs <- fieldsheet[fieldsheet$subsite == subsite,]
  gs <- first(corresp_fs$gas_analyzer)

  if(gs == "LI-COR"){
    gs_folder <- "RAW Data Licor-7810"
  } else if (gs == "Los Gatos"){
    gs_folder <- "RAW Data Los Gatos"
  } else if (gs == "Picarro"){
    gs_folder <- "RAW Data Picarro"
  } else{
    warning("------------> gas analyser not properly detected!")
  }

  path2data <- paste0(datapath,"/",gs_folder,"/RData/",subsite)
  setwd(path2data)
  load(file = paste0("data_",subsite,".RData"))





  # # Create a list of dataframe (by UniqueID)
  # data_split <- mydata %>%
  #   right_join(corresp_fs, by = c("UniqueID")) %>% group_by(UniqueID) %>%
  #   group_split()




  for (incub in seq_along(corresp_fs$plot_id)){
    my_incub <- mydata[as.numeric(mydata$POSIX.time)> corresp_fs$unix_start[incub] &
                         as.numeric(mydata$POSIX.time)< corresp_fs$unix_stop[incub],]

    plt_CO2 <- ggplot(my_incub, aes(POSIX.time, CO2dry_ppm))+geom_line()+
      theme_article()+
      xlab("time UTC")+
      ylab("CO2dry [ppm]")+
      ggtitle(paste0(subsite," plot ",
                     corresp_fs$plot_id[incub]," ",corresp_fs$strata[incub]," ",
                     corresp_fs$transparent_dark[incub]))
    plt_CH4 <- ggplot(my_incub, aes(POSIX.time, CH4dry_ppb))+geom_line()+
      theme_article()+
      xlab("time UTC")+
      ylab("CH4dry [ppm]")
    plt_H2O <- ggplot(my_incub, aes(POSIX.time, H2O_ppm))+geom_line()+
      theme_article()+
      xlab("time UTC")+
      ylab("H2O [ppm]")

    plt <- ggarrange(plt_CO2, plt_CH4, plt_H2O, ncol = 1)
    plt
  }







  isFincub <- T
  for (incub in seq_along(corresp_fs$plot_id)){
    my_incub <- mydata[as.numeric(mydata$POSIX.time)> corresp_fs$unix_start[incub] &
                         as.numeric(mydata$POSIX.time)< corresp_fs$unix_stop[incub],]
    my_incub$UniqueID <- paste0(subsite,"-plot",
                                corresp_fs$plot_id[incub],"-",corresp_fs$strata[incub],"-",
                                corresp_fs$transparent_dark[incub])
    if(isFincub){
      isFincub = F
      my_incub_all <- my_incub
    } else {
      my_incub_all <- rbind(my_incub_all, my_incub)
    }
  }





}



