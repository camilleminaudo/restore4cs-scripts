
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


############################################
   # USER, please specify if you want plots to be saved
doPlot <- F
############################################

# --- install GoFluxYourself if needed ---
#library(devtools)
#install_github("Qepanna/goFlux")

# ---- packages ----
library(tidyverse)
library(readxl)
library(lubridate)
library(zoo)
library(ggplot2)
library(grid)
library(egg)
library(goFlux)
require(dplyr)
require(purrr)
require(msm)


repo_root <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
files.sources = list.files(path = paste0(repo_root,"/functions"), full.names = T)
for (f in files.sources){source(f)}




# ---- Directories ----
dropbox_root <- "C:/Users/misteli/Dropbox/RESTORE4Cs - Fieldwork/Data"
datapath <- paste0(dropbox_root,"/GHG/RAW data")
fieldsheetpath <- paste0(dropbox_root,"/GHG/Fieldsheets")
corrfieldsheetpath <- paste0(dropbox_root,"/GHG/Processed data/corrected_fieldsheets")
loggerspath <- paste0(datapath,"/RAW Data Logger")
plots_path <- paste0(dropbox_root,"/GHG/Processed data/plots_all_incubations/")
results_path <- paste0(dropbox_root,"/GHG/Processed data/computed_flux/")





# ---- List GHG chamber fieldsheets in Dropbox and read them ---
# list filenames
myfieldsheets_list <- list.files(fieldsheetpath, pattern = "Fieldsheet-GHG.xlsx", all.files = T, full.names = T, recursive = T)
# Read all fieldsheets and put them in a single dataframe
fieldsheet <- read_GHG_fieldsheets(myfieldsheets_list)



# ---- List GHG chamber corrected fieldsheets in Dropbox ---
# list filenames
myfieldsheets_corrected_list <- list.files(corrfieldsheetpath, pattern = "Fieldsheet-GHG_corrected.csv", all.files = T, full.names = T, recursive = T)
subsites_corrected_fs <- gsub(pattern = "-Fieldsheet-GHG_corrected.csv",replacement = "",x = basename(myfieldsheets_corrected_list))


# ---- Go through each incubation in fieldsheet and compute linear model for co2 and ch4 ----
subsites <- unique(fieldsheet$subsite)
isF_incub <- T
isFsubsite <- T
for (subsite in subsites[1:36]){
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
    
    auxfile <- NULL
    
    for (incub in seq_along(corresp_fs$plot_id)){
      my_incub <- mydata[as.numeric(mydata$POSIX.time)> corresp_fs$unix_start[incub] &
                           as.numeric(mydata$POSIX.time)< corresp_fs$unix_stop[incub],]
      my_incub <- my_incub[!is.na(my_incub$CO2dry_ppm),]
      
      if (dim(my_incub)[1]>0){
        
        # Compute Temperature, Area and Volume from fieldsheet info
        myTemp <- 15 # a default temperature to run the flux calculation...
        if (corresp_fs$chamber_type[incub] == "floating"){
          if(is_data_logger_float){
            myTemp <- median(approx(data_logger_float$unixtime, data_logger_float$temperature, xout = as.numeric(my_incub$POSIX.time))$y )
          }
          myArea = 14365.4439 # cm2
          myVtot = 115/2 # L
        } else if (corresp_fs$chamber_type[incub] == "tube"){
          if(is_data_logger_tube){
            myTemp <- median(approx(data_logger_tube$unixtime, data_logger_tube$temperature, xout = as.numeric(my_incub$POSIX.time))$y )
          }
          myArea = pi*12.1**2 # cm2
          myVtot = myArea*as.numeric(corresp_fs$chamber_height_cm[incub])*1e-3 # L
        } else {
          warning("chamber type not correct!")
        }
        myPcham <- 100.1 #kPa
        
        # --- Create auxfile table ----
        # An auxfile table, made of fieldsheet, adding important variables. The auxfile
        # requires start.time and UniqueID.
        # start.time must be in the format "%Y-%m-%d %H:%M:%S"
        
        UniqueID = paste(subsite, seq_along(corresp_fs$pilot_site),corresp_fs$strata,corresp_fs$transparent_dark,sep = "-")[incub]
        
        auxfile_tmp <- data.frame(subsite = subsite,
                                  UniqueID = gsub(" ", "", UniqueID, fixed = TRUE),
                                  gas_analiser = gs,
                                  start.time = as.POSIXct((corresp_fs$unix_start[incub]), tz = "UTC"),
                                  duration = (corresp_fs$unix_stop[incub]) - (corresp_fs$unix_start[incub]),
                                  water_depth = corresp_fs$water_depth[incub],
                                  Area = myArea,
                                  Vtot = myVtot,
                                  Tcham = myTemp,
                                  Pcham = myPcham,
                                  strata = corresp_fs$strata[incub],
                                  chamberType = corresp_fs$chamber_type[incub],
                                  lightCondition = corresp_fs$transparent_dark[incub])
        if(is.null(auxfile)){
          auxfile <- auxfile_tmp
        } else {
          auxfile <- rbind(auxfile, auxfile_tmp)
        }
        
      }
    }
    
    auxfile$Tcham[is.na(auxfile$Tcham)] <- mean(auxfile$Tcham, na.rm = T)
    
    auxfile <- auxfile[auxfile$duration>100,]
    auxfile <- auxfile[!is.na(auxfile$Vtot),] # in case chamber height is not specified in the fieldsheet...
    
    # Define the measurements' window of observation
    # auxfile <- auxfile
    mydata_ow <- obs.win(inputfile = mydata, auxfile = auxfile,
                         obs.length = auxfile$duration, shoulder = 2)
    
    # Join mydata_ow with info on start end incubation
    mydata_auto <- lapply(seq_along(mydata_ow), join_auxfile_with_data.loop, flux.unique = mydata_ow) %>%
      map_df(., ~as.data.frame(.x))
    
    # Additional auxiliary data required for flux calculation.
    mydata_auto <- mydata_auto %>%
      left_join(auxfile %>% select(UniqueID, Area, Vtot, Tcham, Pcham))
    
    # Add instrument precision for each gas
    prec = c(3.5, 0.6, 0.4, 45, 45)
    mydata_auto <- mydata_auto %>%
      mutate(CO2_prec = prec[1], CH4_prec = prec[2], N2O_prec = prec[3],
             H2O_prec = prec[4])
    
    # Calculate fluxes
    CO2_results_auto <- goFlux(dataframe = mydata_auto, gastype = "CO2dry_ppm")
    H2O_results_auto <- goFlux(mydata_auto, "H2O_ppm")
    CH4_results_auto <- goFlux(mydata_auto, "CH4dry_ppb")
    
    # Use best.flux to select the best flux estimates (LM or HM)
    # based on a list of criteria
    criteria <- c("g.factor", "kappa", "MDF", "R2", "SE.rel")
    
    CO2_flux_res_auto <- best.flux(CO2_results_auto, criteria)
    H2O_flux_res_auto <- best.flux(H2O_results_auto, criteria)
    CH4_flux_res_auto <- best.flux(CH4_results_auto, criteria)
    
    if(doPlot){
      # Plots results
      # Make a list of plots of all measurements, for each gastype
      CO2_flux_plots <- flux.plot(CO2_flux_res_auto, mydata_auto, "CO2dry_ppm")
      H2O_flux_plots <- flux.plot(H2O_flux_res_auto, mydata_auto, "H2O_ppm")
      CH4_flux_plots <- flux.plot(CH4_flux_res_auto, mydata_auto, "CH4dry_ppb")
      
      # Combine plot lists into one list
      flux_plot.ls <- c(CO2_flux_plots, CH4_flux_plots, H2O_flux_plots)
      
      # Save plots to pdf
      myfilename <- paste(subsite, as.character(as.Date(first(auxfile$start.time))),sep="_")
      flux2pdf(flux_plot.ls, outfile = paste0(results_path,myfilename,".pdf"))
      
    }
    
    
    # estimate ch4 diffusion and ebullition components---------| METHOD 1 |-----------
    CH4_res_meth1 <- CH4_flux_res_auto
    CH4_res_meth1$total_estimated <- NA
    CH4_res_meth1$ebullition <- NA
    CH4_res_meth1$diffusion <- NA
    
    
    for (i in which(auxfile$water_depth>0)){
      if(auxfile$water_depth[i]>0){
        my_incub <- mydata[as.numeric(mydata$POSIX.time)> auxfile$start.time[i] &
                             as.numeric(mydata$POSIX.time)< auxfile$start.time[i]+auxfile$duration[i],]
        my_incub <- my_incub[!is.na(my_incub$CO2dry_ppm),]
        # calling dedicated function
        df_ebull <- separate_ebullition_from_diffusion(my_incub, UniqueID = auxfile$UniqueID[i], doPlot=F)
        # computing fluxes
        H2O_mol = my_incub$H2O_ppm / (1000*1000)
        myfluxterm <- flux.term(auxfile$Vtot[i], auxfile$Pcham[i], auxfile$Area[i],
                                auxfile$Tcham[i], first(H2O_mol))
        CH4_flux_total <- df_ebull$delta_ch4/df_ebull$duration*myfluxterm # nmol/m2/s
        CH4_flux_diff <- df_ebull$avg_diff_slope*myfluxterm # nmol/m2/s
        CH4_flux_ebull <- CH4_flux_total - CH4_flux_diff
      } else {
        CH4_flux_total <- CH4_flux_ebull <- CH4_res_meth1$best.flux[which(CH4_res_meth1$UniqueID==auxfile$UniqueID[i])]
        CH4_flux_ebull <- 0
      }
      CH4_res_meth1$total_estimated[which(CH4_res_meth1$UniqueID==auxfile$UniqueID[i])] <- CH4_flux_total
      CH4_res_meth1$ebullition[which(CH4_res_meth1$UniqueID==auxfile$UniqueID[i])] <- CH4_flux_ebull
      CH4_res_meth1$diffusion[which(CH4_res_meth1$UniqueID==auxfile$UniqueID[i])] <- CH4_flux_diff
    }
    
    CH4_res_meth1$ebullition[which(CH4_res_meth1$ebullition<0)] <- 0
    
    
    setwd(results_path)
    myfilenameCO2 <- paste("co2_fluxes",subsite, as.character(as.Date(first(auxfile$start.time))),sep="_")
    myfilenameCH4 <- paste("ch4_fluxes",subsite, as.character(as.Date(first(auxfile$start.time))),sep="_")
    write.csv(x = CO2_flux_res_auto, file = paste0(myfilenameCO2,".csv"), row.names = F)
    write.csv(x = CH4_res_meth1, file = paste0(myfilenameCH4,".csv"), row.names = F)
    
    
    table_results <- auxfile %>%
      left_join(CO2_flux_res_auto %>% select(UniqueID, LM.flux, HM.flux, best.flux, model, quality.check)) %>%
      rename(CO2_LM.flux = LM.flux, CO2_HM.flux = HM.flux, CO2_best.flux = best.flux, CO2_best.model = model,
             CO2_quality.check = quality.check) %>%
      left_join(CH4_res_meth1 %>% select(UniqueID, LM.flux, HM.flux, best.flux, best.flux, model, quality.check, 
                                         diffusion, ebullition)) %>%
      rename(CH4_LM.flux = LM.flux, CH4_HM.flux = HM.flux, CH4_best.flux = best.flux, CH4_best.model = model, 
             CH4_quality.check = quality.check, CH4_diffusive_flux = diffusion, CH4_ebullitive_flux = ebullition)
    
    if (isFsubsite){
      isFsubsite <- F
      table_results_all <- table_results
    } else {
      table_results_all <- rbind(table_results_all, table_results)
    }
  }
}

setwd(results_path)
myfilename <- paste("000_fluxes_all",min(as.Date(table_results_all$start.time)),
                    max(as.Date(table_results_all$start.time)), sep = "_")
write.csv(x = table_results_all, file = paste0(myfilename,".csv"), row.names = F)



# ---- Some plots ----


table_results_all$campaign_site <- substr(table_results_all$subsite,start = 1, stop = 5)
table_results_all$subsite_short <- substr(table_results_all$subsite,start = 7, stop = 8)


plt_CO2 <- ggplot(table_results_all, aes(subsite_short, CO2_LM.flux,
                                         fill = lightCondition))+
  geom_hline(yintercept = 0)+
  geom_boxplot(alpha=0.5)+
  # geom_jitter(width = 0.2, aes(colour = strata), size=2)+
  theme_article()+
  xlab("subsite")+
  ylab("CO2 flux mmol/m2/s")+
  ggtitle("CO2 flux - Linear Model")+
  scale_fill_viridis_d(begin = 0.2, end = 0.9)+
  scale_colour_viridis_d(begin = 0.2, end = 0.9, option = "C")+
  facet_grid(.~campaign_site)


plt_CH4diff <- ggplot(table_results_all, aes(subsite_short, CH4_LM.flux, 
                                             fill = lightCondition))+
  geom_hline(yintercept = 0)+
  geom_boxplot(alpha=0.5)+
  # geom_jitter(width = 0.2, aes(colour = strata), size=2)+
  theme_article()+
  xlab("subsite")+
  ylab("CH4 flux nmol/m2/s")+
  ggtitle("CH4 flux - Linear Model")+
  scale_fill_viridis_d(begin = 0.2, end = 0.9)+
  scale_colour_viridis_d(begin = 0.2, end = 0.9, option = "C")+
  facet_grid(.~campaign_site)

plt_all <- ggarrange(plt_CO2, plt_CH4diff, ncol = 1)

ggsave(plot = plt_all, filename = paste0(myfilename,".jpg"), path = results_path, 
       width = 10, height = 8, dpi = 300, units = 'in')



for (cs in unique(table_results_all$campaign_site)){
  
  table_results_cs <- table_results_all[table_results_all$campaign_site == cs,]
  
  plt_CO2 <- ggplot(table_results_cs, aes(subsite_short, CO2_LM.flux,
                                           fill = lightCondition))+
    geom_hline(yintercept = 0)+
    geom_boxplot(alpha=0.5)+
    geom_jitter(width = 0.2, size=2, alpha=0.5)+
    theme_article()+
    xlab("subsite")+
    ylab("CO2 flux mmol/m2/s")+
    ggtitle("CO2 flux - Linear Model")+
    scale_fill_viridis_d(begin = 0.2, end = 0.9)+
    scale_colour_viridis_d(begin = 0.2, end = 0.9, option = "C")+
    facet_grid(.~campaign_site)
  
  
  plt_CH4diff <- ggplot(table_results_cs, aes(subsite_short, CH4_LM.flux, 
                                               fill = lightCondition))+
    geom_hline(yintercept = 0)+
    geom_boxplot(alpha=0.5)+
    geom_jitter(width = 0.2, size=2, alpha=0.5)+
    theme_article()+
    xlab("subsite")+
    ylab("CH4 flux nmol/m2/s")+
    ggtitle("CH4 flux - Linear Model")+
    scale_fill_viridis_d(begin = 0.2, end = 0.9)+
    scale_colour_viridis_d(begin = 0.2, end = 0.9, option = "C")+
    facet_grid(.~campaign_site)
  
  plt_cs <- ggarrange(plt_CO2, plt_CH4diff, ncol = 1)
  
  myfilename <- paste("000_fluxes",cs,min(as.Date(table_results_cs$start.time)), sep = "_")
  
  ggsave(plot = plt_cs, filename = paste0(myfilename,".jpg"), path = results_path, 
         width = 6, height = 8, dpi = 300, units = 'in')
  
}




# 
# 
# table_results_water <- table_results_all[table_results_all$strata=="open water",]
# 
# plt_CO2 <- ggplot(table_results_water, aes(subsite_short, CO2_LM.flux,
#                                          fill = lightCondition))+
#   geom_hline(yintercept = 0)+
#   geom_boxplot(alpha=0.5)+
#   geom_jitter(width = 0.2, size=2)+
#   theme_article()+
#   xlab("subsite")+
#   ylab("CO2 flux mmol/m2/s")+
#   ggtitle("CO2 flux")+
#   scale_fill_viridis_d(begin = 0.2, end = 0.9)+
#   scale_colour_viridis_d(begin = 0.2, end = 0.9, option = "C")+
#   facet_grid(.~campaign_site)
# 
# 
# plt_CH4diff <- ggplot(table_results_water, aes(subsite_short, CH4_LM.flux, 
#                                              fill = lightCondition))+
#   geom_hline(yintercept = 0)+
#   geom_boxplot(alpha=0.5)+
#   geom_jitter(width = 0.2, size=2)+
#   theme_article()+
#   xlab("subsite")+
#   ylab("CH4 flux nmol/m2/s")+
#   ggtitle("CH4 flux")+
#   scale_fill_viridis_d(begin = 0.2, end = 0.9)+
#   scale_colour_viridis_d(begin = 0.2, end = 0.9, option = "C")+
#   facet_grid(.~campaign_site)
# 
# ggarrange(plt_CO2, plt_CH4diff, ncol = 1)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ggplot(table_results_all, aes(CO2_LM.flux, CO2_HM.flux, colour = strata))+
#   geom_abline(slope = 1, intercept = 0)+
#   geom_point(size=2, alpha=0.5)+
#   theme_article()+
#   scale_colour_viridis_d(begin = 0.2, end = 0.9)#+facet_grid(.~subsite)
