
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
require(pbapply)


source(paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/click.peak.R"))
source(paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/click.peak.loop.R"))
source(paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/flux.term.R"))
# source(paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/read_Licor.R"))
source(paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/get_unix_times.R"))
source(paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/G2508_import.R"))
source(paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/import2RData.R"))



# ---- Directories ----

datapath <- "C:/Users/Camille Minaudo/Dropbox/RESTORE4Cs - Fieldwork/Data/GHG/RAW data/RAW Data Picarro/S1-CA-P1/"
fieldsheetpath <- datapath
# loggerspath <- paste0(datapath,"/RAW Data Logger")

# path_to_L0b <- paste0(dropbox_root,"/GHG/Processed data")

setwd(datapath)

# ---- SETTINGS ----
# site_ID <- "S1-CU"
# subsite_ID <- "S1-CU-A2"
# file_to_read <- "NOMAD-20231108-094148Z-DataLog_User_Minimal.dat"
who_runs_this <- "Camille Minaudo"


# Read corresponding Fieldsheet
path2file <- paste0(fieldsheetpath,"/","MONISTROL_20231031_085029.csv")
fieldsheet <- read.csv(file = path2file,
                                     header = T)

fieldsheet$unix_start_time <- fieldsheet$start_fit
fieldsheet$unix_end_time <- fieldsheet$end_fit

fieldsheet <- fieldsheet[fieldsheet$Species ==  "CH4",]
# head(fieldsheet)

analyser <- "Picarro" # this could create a problem if
# several gas analyzers are used in
# the same day at the same subsite

# mydata_imp <- G2508_import(inputfile = file_to_read, timezone = "UTC")

setwd(datapath)
import2RData(path = datapath, instrument = "G2508", date.format = "ymd", timezone = 'UTC')

# load all these R.Data
file_list <- list.files(path = paste(datapath,"RData",sep="/"), full.names = T)
isF <- T
for(i in seq_along(file_list)){
  load(file_list[i])
  if(isF){
    isF <- F
    mydata_imp <- data.raw
  } else {
    mydata_imp <- rbind(mydata_imp, data.raw)
  }
  rm(data.raw)
}



# --- Create auxfile table ----
# An auxfile table, made of fieldsheet, adding important variables. The auxfile
# requires start.time and UniqueID.
# start.time must be in the format "%Y-%m-%d %H:%M:%S"
auxfile <- NULL
for (i in c(1,2,3,4)){
# for (i in seq_along(fieldsheet$Time.Code)){

  # my_sel <- mydata_imp[as.numeric(mydata_imp$POSIX.time)>= (fieldsheet$unix_start_time[i]) & as.numeric(mydata_imp$POSIX.time)<= (fieldsheet$unix_end_time[i]),]

  UniqueID = fieldsheet$Comment[i]

  myTcham <- fieldsheet$chamber_temperature[i]
  myPcham <- 100.92503288 #kPa
  myArea = 14365.4439 # cm2
  myVtot = 115 # L

  auxfile_tmp <- data.frame(UniqueID = gsub(" ", "", UniqueID, fixed = TRUE),
                            gas_analiser = analyser,
                            start.time = as.POSIXct((fieldsheet$unix_start_time[i]), tz = "UTC"),
                            duration = (fieldsheet$unix_end_time[i]) - (fieldsheet$unix_start_time[i]),
                            Area = myArea,
                            Vtot = myVtot,
                            Tcham = myTcham,
                            Pcham = myPcham,
                            strata = "NA",
                            chamberType = "NA",
                            lightCondition = "NA")
  if(is.null(auxfile)){
    auxfile <- auxfile_tmp
  } else {
    auxfile <- rbind(auxfile, auxfile_tmp)
  }
}

#----- compute CO2 fluxes -----

# Define the measurements' window of observation
mydata_ow <- obs.win(inputfile = mydata_imp, auxfile = auxfile, gastype = "CO2dry_ppm",
                     obs.length = auxfile$duration, shoulder = 30)
mydata_ow

# Manually identify measurements by clicking on the start and end points
myCO2data_manID <- lapply(seq_along(mydata_ow), click.peak.loop,
                          flux.unique = mydata_ow,
                          gastype = "CO2dry_ppm",
                          plot.lim = c(200,1000)) %>%
  map_df(., ~as.data.frame(.x))

# Additional auxiliary data required for flux calculation.
myCO2data_manID <- myCO2data_manID %>%
  left_join(auxfile %>% select(UniqueID, Area, Vtot, Tcham, Pcham))


# saving CO2 timeseries file (L0B)
# setwd(path_to_L0b)
# myfilename <- paste(subsite_ID, as.character(as.Date(first(myCO2data_manID$DATE))),analyser,"CO2",sep="_")
# write.csv(x = myCO2data_manID, file = paste0(myfilename, ".csv"), sep = ";", dec = ".", row.names = F, col.names = T)


# Calculate fluxes for CO2 and H2O
CO2_results <- goFlux(myCO2data_manID, "CO2dry_ppm")
H2O_results <- goFlux(myCO2data_manID, "H2O_ppm")

# Use best.flux to select the best flux estimates (LM or HM)
# based on a list of criteria
criteria <- c("g.factor", "kappa", "MDF", "R2", "SE.rel")

CO2_flux_res <- best.flux(CO2_results, criteria)
H2O_flux_res <- best.flux(H2O_results, criteria)

CO2_flux_res <- CO2_flux_res %>%
  left_join(auxfile %>% select(UniqueID, strata, chamberType, lightCondition))


# Plots results
# Make a list of plots of all measurements, for each gastype
CO2_flux_plots <- flux.plot(CO2_flux_res, myCO2data_manID, "CO2dry_ppm")
# CO2_flux_plots
H2O_flux_plots <- flux.plot(H2O_flux_res, myCO2data_manID, "H2O_ppm")

# a simple boxplot for CO2 fluxes
ggplot(CO2_flux_res, aes(lightCondition, best.flux, fill = lightCondition))+
  geom_boxplot(alpha=0.2)+geom_jitter(width = 0.2)+
  theme_article()+facet_wrap(.~strata, scales = "free")+
  scale_fill_viridis_d(begin = 0.2, end = 0.9)

ggplot(H2O_flux_res, aes(lightCondition, best.flux, fill = lightCondition))+
  geom_boxplot(alpha=0.2)+geom_jitter(width = 0.2)+
  theme_article()+facet_wrap(.~strata, scales = "free")+
  scale_fill_viridis_d(begin = 0.2, end = 0.9)





#----- computing CH4 fluxes -----


# Define the measurements' window of observation
mydata_ow <- obs.win(inputfile = mydata_imp, auxfile = auxfile,
                     obs.length = auxfile$duration, shoulder = 30)


# Manually identify start/end CH4 measurements by clicking on the start and end points
myCH4data_manID <- lapply(seq_along(mydata_ow), click.peak.loop,
                          flux.unique = mydata_ow,
                          gastype = "CH4dry_ppb",
                          plot.lim = c(1900,4000)) %>%
  map_df(., ~as.data.frame(.x))


# saving CH4 timeseries file (L0B)
# setwd(path_to_L0b)
# myfilename <- paste(subsite_ID, as.character(as.Date(first(myCH4data_manID$DATE))),analyser,"CH4",sep="_")
# myfilename <- paste0(myfilename, ".csv")
# write.csv(x = myCH4data_manID, file = myfilename, sep = ";", dec = ".", row.names = F, col.names = T)


# linking window of start/end time to mydata_ow list
i=0
for (id in unique(myCH4data_manID$UniqueID)){
  i=i+1
  mydata_ow[[i]]$start.time <- unique(myCH4data_manID$start.time_corr[myCH4data_manID$UniqueID == id])
  mydata_ow[[i]]$end.time <- unique(myCH4data_manID$end.time[myCH4data_manID$UniqueID == id])
  mydata_ow[[i]]$duration <- as.numeric(mydata_ow[[i]]$end.time) - as.numeric(mydata_ow[[i]]$start.time)
  mydata_ow[[i]]$c0 <- first(myCH4data_manID$c0[myCH4data_manID$UniqueID == id])

  auxfile$c0[i] <- first(myCH4data_manID$c0[myCH4data_manID$UniqueID == id])
  auxfile$cf[i] <- first(myCH4data_manID$cf[myCH4data_manID$UniqueID == id])
}

# Manually identify diffusive (more or less linear) CH4 behaviors by clicking on the start and end points
myCH4_diffusion <- lapply(seq_along(mydata_ow), click.peak.loop,
                          flux.unique = mydata_ow,
                          gastype = "CH4dry_ppb",
                          plot.lim = c(1900,4000)) %>%
  map_df(., ~as.data.frame(.x))


# Calculate fluxes for CH4
CH4_results_diffusion <- goFlux(myCH4_diffusion, "CH4dry_ppb")

# Use best.flux to select the best flux estimates (LM or HM)
# based on a list of criteria
criteria <- c("g.factor", "kappa", "MDF", "R2", "SE.rel")

CH4_flux_res <- best.flux(CH4_results_diffusion, criteria)

CH4_flux_plots <- flux.plot(CH4_flux_res, myCH4_diffusion, "CH4dry_ppb")
# CH4_flux_plots


table_results_CH4 <- auxfile %>%
  left_join(CH4_flux_res %>% select(UniqueID, best.flux, model))



# Estimating ebullition component
i=0
for (id in unique(myCH4data_manID$UniqueID)){
  i=i+1

  CH4_final <- table_results_CH4$cf[i]
  CH4_initial <-  table_results_CH4$c0[i]
  incubation_time <- first(mydata_ow[[i]]$duration)

  H2O_mol = mydata_ow[[i]]$H2O_ppm / (1000*1000)
  myfluxterm <- flux.term(table_results_CH4$Vtot[i], table_results_CH4$Pcham[i], table_results_CH4$Area[i],
                          table_results_CH4$Tcham[i], first(H2O_mol))

  CH4_flux_total <- (CH4_final-CH4_initial)/incubation_time*myfluxterm # ppb/m2/s


  CH4_ebullition <- CH4_flux_total - CH4_flux_res$best.flux[i] # total flux - diffusive term
  CH4_ebullition[CH4_ebullition<0] <- 0

  table_results_CH4$CH4_ebullition[i] <- CH4_ebullition
}



ggplot(table_results_CH4, aes(lightCondition, best.flux, fill = lightCondition))+
  geom_boxplot(alpha=0.2)+geom_jitter(width = 0.2)+
  theme_article()+facet_wrap(.~strata, scales = "free")+
  ylab("CH4 flux nmol/m2/s")+
  scale_fill_viridis_d(begin = 0.2, end = 0.9)





