
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
require(msm)

source(paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/get_unix_times.R"))
source(paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/read_GHG_fieldsheets.R"))
source(paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/flux.term.R"))
source(paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/LM.flux.R"))



# ---- Directories ----
dropbox_root <- "C:/Users/Camille Minaudo/Dropbox/RESTORE4Cs - Fieldwork/Data"
datapath <- paste0(dropbox_root,"/GHG/RAW data")
fieldsheetpath <- paste0(dropbox_root,"/GHG/Fieldsheets")
corrfieldsheetpath <- paste0(dropbox_root,"/GHG/Processed data/corrected_fieldsheets")
loggerspath <- paste0(datapath,"/RAW Data Logger")
plots_path <- paste0(dropbox_root,"/GHG/Processed data/plots_all_incubations/")
results_path <- paste0(dropbox_root,"/GHG/Processed data/computed_flux/")


doPlot <- F

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

# ---- Go through each incubation in fieldsheet and compute linear model for co2 and ch4 ----
subsites <- unique(fieldsheet$subsite)
isF_incub <- T
for (subsite in subsites){
  message("Now processing ",subsite)

  corresp_fs <- fieldsheet[fieldsheet$subsite == subsite,]
  gs <- first(corresp_fs$gas_analyzer)

  # check if there is a corresponding corrected fielsheet for this subsite
  existsCorr_fs <- which(subsites_corrected_fs == subsite)
  if (length(existsCorr_fs) > 0){
    message("using exact incubation start and stop times")
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


  # read corresponding temperature logger file and keep initial temperature

  # --- Read corresponding Loggers data ----
  SN_logger_float <- first(corresp_fs$logger_floating_chamber)
  SN_logger_tube <- first(corresp_fs$logger_transparent_chamber)
  site_ID <- str_sub(subsite, start = 1, end = 5)

  # finding out if corresponding file exists, and its extension
  require(tools)
  dir_exists_loggdata <- dir.exists(paste0(loggerspath,"/",site_ID,"/"))
  if(dir_exists_loggdata){
    f <- list.files(paste0(loggerspath,"/",site_ID,"/"), full.names = T)
    r <- grep(pattern = ".hobo", x = f)
    if(length(r)>0){f <- f[-r]}
    r <- grep(pattern = ".txt", x = f)
    if(length(r)>0){f <- f[-r]}
    f_ext <- file_ext(f)
    i_f_float <- grep(pattern = SN_logger_float, x = f)[1]
    i_f_tube <- grep(pattern = SN_logger_tube, x = f)[1]
  }

  if(!is.na(SN_logger_float) & !is.na(i_f_float)){
    is_data_logger_float = T
    message("...reading corresponding temperature logger file for floating chamber")
    if(f_ext[i_f_float]=="xlsx"){
      data_logger_float <- readxl::read_xlsx(f[i_f_float],col_names = T)
    } else if(f_ext[i_f_float]=="csv"){
      data_logger_float <- read.csv(f[i_f_float], header = T)
    }
    data_logger_float <- data_logger_float[,seq(1,3)]
    names(data_logger_float) <- c("sn","datetime","temperature")
    data_logger_float$sn <- SN_logger_float
    if(is.character(data_logger_float$datetime)){
      data_logger_float$datetime <- as.POSIXct(data_logger_float$datetime, tz = 'utc', tryFormats = c("%m/%d/%y %r", "%d/%m/%Y %H:%M"))
    }
    data_logger_float$unixtime <- as.numeric(data_logger_float$datetime)
  } else {
    message("===> no data logger linked to the floating chamber!")
    is_data_logger_float = F}

  if(!is.na(SN_logger_tube) & !is.na(i_f_tube)){
    is_data_logger_tube = T
    message("...reading corresponding temperature logger file for tube chamber")
    if(f_ext[i_f_tube]=="xlsx"){
      data_logger_tube <- readxl::read_xlsx(f[i_f_tube],col_names = T)
    } else if(f_ext[i_f_tube]=="csv"){
      data_logger_tube <- read.csv(f[i_f_tube], header = T, fill = T)
    }
    data_logger_tube <- data_logger_tube[,seq(1,4)]
    names(data_logger_tube) <- c("sn","datetime","temperature","light")
    data_logger_tube$sn <- SN_logger_tube
    if(is.character(data_logger_tube$datetime)){
      if(length(grep(pattern = "AM", x = first(data_logger_tube$datetime)))>0 | length(grep(pattern = "PM", x = first(data_logger_tube$datetime)))>0){
        data_logger_tube$datetime <- as.POSIXct(data_logger_tube$datetime, tz = 'UTC', format = c("%m/%d/%y %r"))
      } else {
        data_logger_tube$datetime <- as.POSIXct(data_logger_tube$datetime, tz = 'UTC', format = "%m/%d/%Y %H:%M")
      }
    }
    data_logger_tube$unixtime <- as.numeric(data_logger_tube$datetime)
  } else {
    is_data_logger_tube = F
    message("===> no data logger linked to the tube chamber!")
  }


  path2data <- paste0(datapath,"/",gs_folder,"/RData/",subsite)
  if(dir.exists(path2data)){
    setwd(path2data)
    load(file = paste0("data_",subsite,".RData"))

    if(doPlot){plt_list <- vector('list', length(corresp_fs$plot_id))}

    for (incub in seq_along(corresp_fs$plot_id)){
      my_incub <- mydata[as.numeric(mydata$POSIX.time)> corresp_fs$unix_start[incub] &
                           as.numeric(mydata$POSIX.time)< corresp_fs$unix_stop[incub],]
      my_incub <- my_incub[!is.na(my_incub$CO2dry_ppm),]
      if (dim(my_incub)[1]>0){

        my_incub$elapsed_time <- as.numeric(my_incub$POSIX.time - corresp_fs$unix_start[incub])

        if(doPlot){
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

        # Compute Temperature, Area and Volume from fieldsheet info
        myTemp <- 15 # a default temperature to run the flux calculation...
        if (corresp_fs$chamber_type[incub] == "floating"){
          if(is_data_logger_float){
            myTemp <- median(approx(data_logger_float$unixtime, data_logger_float$temperature, xout = as.numeric(my_incub$POSIX.time))$y )
          }
          myArea = 14365.4439 # cm2
          myVtot = 115 # L
        } else if (corresp_fs$chamber_type[incub] == "tube"){
          if(is_data_logger_tube){
            myTemp <- median(approx(data_logger_tube$unixtime, data_logger_tube$temperature, xout = as.numeric(my_incub$POSIX.time))$y )
          }
          myArea = pi*12.1**2 # cm2
          myVtot = myArea*as.numeric(corresp_fs$chamber_height_cm[incub])*1e-3 # L
        } else {
          warning("chamber type not correct!")
        }


        # Fit a simple linear model using functions from GoFluxYourself package
        myfluxterm <- flux.term(V_L = myArea,
                                P_kPa = 100,
                                A_cm2 = myArea,
                                T_C = myTemp,
                                H2O_mol = first(my_incub$H2O_ppm)/(1000*1000))
        lm_co2 <- LM.flux(gas.meas = my_incub$CO2dry_ppm,
                          time.meas = my_incub$elapsed_time,
                          flux.term = myfluxterm)
        lm_ch4 <- LM.flux(gas.meas = my_incub$CH4dry_ppb*1e-3, #in ppm
                          time.meas = my_incub$elapsed_time,
                          flux.term = myfluxterm)

        names(lm_co2) <- paste0(names(lm_co2),"_co2")
        names(lm_ch4) <- paste0(names(lm_ch4),"_ch4")

        flux_table_incub <- cbind(corresp_fs[incub,], lm_co2, lm_ch4)

        if(isF_incub) {
          isF_incub <- F
          flux_table_incub_all <- flux_table_incub
        } else {
          flux_table_incub_all <- rbind(flux_table_incub_all, flux_table_incub)
        }
      }
    }
    if(doPlot){
      # Print pdf
      setwd(plots_path)
      gg_save_pdf(list = plt_list, filename = paste0(subsite,".pdf"))
    }
  }
}

setwd(results_path)
myfilename <- paste("flux_all_incubations_linear_model",min(flux_table_incub_all$date),max(flux_table_incub_all$date), sep = "_")
write.csv(x = flux_table_incub_all, file = paste0(myfilename,".csv"), row.names = F)



# ---- Some plots ----

ggplot(flux_table_incub_all, aes(LM.r2_co2, fill = strata))+geom_density(alpha=0.5)+
  theme_article()

ggplot(flux_table_incub_all, aes(LM.r2_ch4, fill = strata))+geom_density(alpha=0.5)+
  theme_article()


ggplot(flux_table_incub_all[flux_table_incub_all$LM.r2_co2>0.5,], aes(subsite, LM.flux_co2, colour = transparent_dark))+
  geom_hline(yintercept = 0, alpha = 0.2)+
  geom_boxplot()+facet_grid(strata~.)+
  theme_article()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ylab("CO2 flux [µmol/m2/s]")



ggplot(flux_table_incub_all[flux_table_incub_all$LM.r2_ch4>0.5,], aes(subsite, LM.flux_ch4, colour = transparent_dark))+
  geom_hline(yintercept = 0, alpha = 0.2)+
  geom_boxplot()+facet_grid(strata~.)+
  theme_article()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_y_log10()+ylab("CH4 flux [µmol/m2/s]")





ggplot(flux_table_incub_all, aes(pilot_site, LM.flux_co2, colour = transparent_dark))+
  geom_hline(yintercept = 0, alpha = 0.2)+
  geom_boxplot()+facet_grid(strata~.)+
  theme_article()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ylab("CO2 flux [µmol/m2/s]")



ggplot(flux_table_incub_all, aes(pilot_site, LM.flux_ch4, colour = transparent_dark))+
  geom_hline(yintercept = 0, alpha = 0.2)+
  geom_boxplot()+facet_grid(strata~.)+
  theme_article()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_y_log10()+
  ylab("CH4 flux [µmol/m2/s]")

