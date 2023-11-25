
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
corrfieldsheetpath <- paste0(dropbox_root,"/GHG/Processed data/corrected_fieldsheets")
plots_path <- paste0(dropbox_root,"/GHG/Processed data/plots_all_incubations/")

# ---- List GHG chamber fieldsheets in Dropbox and read them ---
# list filenames
myfieldsheets_list <- list.files(fieldsheetpath, pattern = "Fieldsheet-GHG.xlsx", all.files = T, full.names = T, recursive = T)
# Read all fieldsheets and put them in a single dataframe
fieldsheet <- read_GHG_fieldsheets(myfieldsheets_list)



# ---- List GHG chamber corrected fieldsheets in Dropbox ---
# list filenames
myfieldsheets_corrected_list <- list.files(corrfieldsheetpath, pattern = "Fieldsheet-GHG_corrected.csv", all.files = T, full.names = T, recursive = T)
subsites_corrected_fs <- gsub(pattern = "-Fieldsheet-GHG_corrected.csv",replacement = "",x = basename(myfieldsheets_corrected_list))


# ---- function to save a list of plots into pdf file ----

gg_save_pdf = function(list, filename) {
  pdf(filename)
  for (p in list) {
    print(p)
  }
  dev.off()
  invisible(NULL)
}

# ---- Go through each incubation in fieldsheet and make a plot, organized by subsite ----
subsites <- unique(fieldsheet$subsite)
for (subsite in subsites[31:36]){
  message("Now processing ",subsite)

  corresp_fs <- fieldsheet[fieldsheet$subsite == subsite,]
  gs <- first(corresp_fs$gas_analyzer)

  # check if there is a corresponding corrected fielsheet for this subsite
  existsCorr_fs <- which(subsites_corrected_fs == subsite)
  if (length(existsCorr_fs) > 0){
    fs_corr <- read.csv(file = myfieldsheets_corrected_list[existsCorr_fs])
    corresp_fs$unix_start <- fs_corr$unix_start_corrected[match(corresp_fs$start_time, fs_corr$start_time)]
    corresp_fs$unix_stop <- fs_corr$unix_stop_corrected[match(corresp_fs$start_time, fs_corr$start_time)]
    corresp_fs <- corresp_fs[!is.na(corresp_fs$unix_start),]
  }

  if(gs == "LI-COR"){
    gs_folder <- "RAW Data Licor-7810"
  } else if (gs == "Los Gatos"){
    gs_folder <- "RAW Data Los Gatos"
  } else if (gs == "Picarro"){
    gs_folder <- "RAW Data Picarro"
  } else{
    warning("------> gas analyser not properly detected!")
  }

  path2data <- paste0(datapath,"/",gs_folder,"/RData/",subsite)
  if(dir.exists(path2data)){
    setwd(path2data)
    load(file = paste0("data_",subsite,".RData"))

    plt_list <- vector('list', length(corresp_fs$plot_id))
    for (incub in seq_along(corresp_fs$plot_id)){
      my_incub <- mydata[as.numeric(mydata$POSIX.time)> corresp_fs$unix_start[incub] &
                           as.numeric(mydata$POSIX.time)< corresp_fs$unix_stop[incub],]
      my_incub <- my_incub[!is.na(my_incub$CO2dry_ppm),]
      if (dim(my_incub)[1]>0){
        plt_CO2 <- ggplot(my_incub, aes(POSIX.time, CO2dry_ppm))+geom_line()+
          theme_article()+
          xlab("time UTC")+
          ylab("CO2dry [ppm]")+
          ggtitle(paste0(subsite," plot ",
                         corresp_fs$plot_id[incub]," ",corresp_fs$strata[incub]," ",
                         corresp_fs$transparent_dark[incub], ", depth = ",corresp_fs$water_depth[incub], " cm"))
        plt_CH4 <- ggplot(my_incub, aes(POSIX.time, CH4dry_ppb))+geom_line()+
          theme_article()+
          xlab("time UTC")+
          ylab("CH4dry [ppm]")
        plt_H2O <- ggplot(my_incub, aes(POSIX.time, H2O_ppm))+geom_line()+
          theme_article()+
          xlab("time UTC")+
          ylab("H2O [ppm]")

        # plt <- ggarrange(plt_CO2, plt_CH4, plt_H2O, ncol = 1)
        plt_list[[incub]] <- ggarrange(plt_CO2, plt_CH4, plt_H2O, ncol = 1)
      }
    }
    # Print pdf
    setwd(plots_path)
    gg_save_pdf(list = plt_list, filename = paste0(subsite,".pdf"))
  }
}



