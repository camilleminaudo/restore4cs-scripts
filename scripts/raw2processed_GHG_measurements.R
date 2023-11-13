
# ---
# Authors: Camille Minaudo
# Project: "RESTORE4Cs"
# date: "Oct 2023"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
# This script loads raw measurement files (data level L0a) from one of the gas
# analyzers used in the project, transform it into a unified harmonized csv
# file (data level L0b), and computes CO2 and CH4 fluxes


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


source(paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/click.peak.R"))
source(paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/click.peak.loop.R"))
source(paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/flux.term.R"))
source(paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/read_Licor.R"))
source(paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/get_unix_times.R"))
source(paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/G2508_import.R"))
source(paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/import2RData.R"))



# ---- Directories ----
dropbox_root <- "C:/Users/Camille Minaudo/Dropbox/RESTORE4Cs - Fieldwork/Data"

datapath <- paste0(dropbox_root,"/GHG/RAW data")
fieldsheetpath <- paste0(dropbox_root,"/GHG/Fieldsheets")
loggerspath <- paste0(datapath,"/RAW Data Logger")

path_to_L0b <- paste0(dropbox_root,"/GHG/Processed data")

setwd(datapath)

# ---- SETTINGS ----
site_ID <- "S1-VA"
subsite_ID <- "S1-VA-P2"
# file_to_read <- "S1-CU-R2 and Cores.data"
who_runs_this <- "Camille Minaudo"


# Read corresponding Fieldsheet
path2file <- paste0(fieldsheetpath,"/",site_ID,"/",subsite_ID,"-Fieldsheet-GHG.xlsx")
fieldsheet_temp <- readxl::read_xlsx(path2file,
                                     col_names = T)
fieldsheet <- readxl::read_xlsx(path2file,
                                skip = 2, col_names = F)
names(fieldsheet) <- names(fieldsheet_temp)

fieldsheet$unix_start_time <- get_unix_times(mydate = fieldsheet$date, mytime = fieldsheet$start_time)
fieldsheet$unix_end_time <- get_unix_times(mydate = fieldsheet$date, mytime = fieldsheet$end_time)
# head(fieldsheet)

analyser <- first(fieldsheet$gas_analyzer) # this could create a problem if
                                           # several gas analyzers are used in
                                           # the same day at the same subsite

# Read gas analyser's file
if(analyser == "LI-COR"){
  directory_analyser <- "RAW Data Licor-7810"
  # my_data <- read_Licor(file = paste(datapath,directory_analyser,file_to_read, sep = "/"))

  mydata_imp <- LI7810_import(inputfile = paste(datapath,directory_analyser,file_to_read, sep = "/"))

} else if(analyser == "Los Gatos"){
  directory_analyser <- paste(datapath,"RAW Data Los Gatos",subsite_ID,sep="/")
  setwd(directory_analyser)
  import2RData(path = directory_analyser, instrument = "LGR", date.format = "mdy", timezone = 'UTC')

  # load all these R.Data
  file_list <- list.files(path = paste(directory_analyser,"RData",sep="/"), full.names = T)
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

} else if (analyser == "Picarro"){
  directory_analyser <- paste(datapath,"RAW Data Picarro",subsite_ID,sep="/")
  setwd(directory_analyser)
  import2RData(path = directory_analyser, instrument = "G2508", date.format = "ymd", timezone = 'UTC')

  # load all these R.Data
  file_list <- list.files(path = paste(directory_analyser,"RData",sep="/"), full.names = T)
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


}





# --- Read corresponding Loggers data ----
SN_logger_float <- first(fieldsheet$logger_floating_chamber)
SN_logger_tube <- first(fieldsheet$logger_transparent_chamber)

if(SN_logger_float != "NA"){
  is_data_logger_float = T
  data_logger_float <- readxl::read_xlsx(paste0(loggerspath,"/",site_ID,"/",site_ID,"-",SN_logger_float,".xlsx"),col_names = T)
  names(data_logger_float) <- c("sn","datetime","temperature","unixtime")
  data_logger_float$sn <- SN_logger_float
  data_logger_float$unixtime <- as.numeric(data_logger_float[[2]])
} else {
  message("no data logger linked to the floating chamber!")
  is_data_logger_float = F}

if(SN_logger_tube != "NA"){
  is_data_logger_tube = T
  data_logger_tube <- readxl::read_xlsx(paste0(loggerspath,"/",site_ID,"/",site_ID,"-",SN_logger_tube,".xlsx"),col_names = T)
  data_logger_tube <- data_logger_tube[,seq(1,4)]
  names(data_logger_tube) <- c("sn","datetime","temperature","light","unixtime")
  data_logger_tube$sn <- SN_logger_tube
  data_logger_tube$unixtime <- as.numeric(data_logger_tube[[2]])
} else {
  is_data_logger_tube = F
  message("no data logger linked to the transparent tube chamber!")
  if(SN_logger_float != "NA"){
    message(".../ using data from floating chamber instead...")
    data_logger_tube = data_logger_float
    is_data_logger_tube = T
  }
}

# --- Create auxfile table ----
# An auxfile table, made of fieldsheet, adding important variables. The auxfile
# requires start.time and UniqueID.
# start.time must be in the format "%Y-%m-%d %H:%M:%S"
auxfile <- NULL
for (i in 5:10){
# for (i in seq_along(fieldsheet$pilot_site)){

  my_sel <- mydata_imp[as.numeric(mydata_imp$POSIX.time)>= (fieldsheet$unix_start_time[i]+30) & as.numeric(mydata_imp$POSIX.time)<= (fieldsheet$unix_end_time[i]-30),]

  UniqueID = paste("plot",fieldsheet$plot_id[i],"_",fieldsheet$chamber_type[i],"_",fieldsheet$transparent_dark[i],
                   sep = "")



  if (!is_data_logger_float & !is_data_logger_tube){
    myTcham = 15
  } else {
    if (fieldsheet$chamber_type[i] == "floating"){
      my_sel$temperature <- approx(data_logger_float$unixtime, data_logger_float[[3]], xout = as.numeric(my_sel$POSIX.time))$y
    } else if (fieldsheet$chamber_type[i] == "tube"){
      my_sel$temperature <- approx(data_logger_tube$unixtime, data_logger_tube[[3]], xout = as.numeric(my_sel$POSIX.time))$y
    } else {
      warning("chamber type not correct!")
    }
    myTcham = mean(my_sel$temperature)
  }

  if (fieldsheet$chamber_type[i] == "floating"){
    myArea = 14365.4439 # cm2
    myVtot = 115 # L
  } else if (fieldsheet$chamber_type[i] == "tube"){
    myArea = pi*12.1**2 # cm2
    myVtot = myArea*fieldsheet$chamber_height_cm[i]*1e-3 # L
  } else {
    warning("chamber type not correct!")
  }

  auxfile_tmp <- data.frame(UniqueID = UniqueID,
                            gas_analiser = analyser,
                            start.time = as.POSIXct((fieldsheet$unix_start_time[i]-30), tz = "UTC"),
                            duration = (fieldsheet$unix_end_time[i]+30) - (fieldsheet$unix_start_time[i]-30),
                            Area = myArea,
                            Vtot = myVtot,
                            Tcham = myTcham,
                            Pcham = 99.4,
                            strata = fieldsheet$strata[i],
                            chamberType = fieldsheet$chamber_type[i],
                            lightCondition = fieldsheet$transparent_dark[i])
  if(is.null(auxfile)){
    auxfile <- auxfile_tmp
  } else {
    auxfile <- rbind(auxfile, auxfile_tmp)
  }
}


#----- compute CO2 fluxes -----

# Define the measurements' window of observation
mydata_ow <- obs.win(inputfile = mydata_imp, auxfile = auxfile,
                     obs.length = auxfile$duration, shoulder = 30)

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
setwd(path_to_L0b)
myfilename <- paste(subsite_ID, as.character(as.Date(first(myCO2data_manID$DATE))),analyser,"CO2",sep="_")
write.csv(x = myCO2data_manID, file = paste0(myfilename, ".csv"), sep = ";", dec = ".", row.names = F, col.names = T)


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
                       plot.lim = c(1900,max(fieldsheet$final_ch4)*1000)) %>%
  map_df(., ~as.data.frame(.x))


# saving CH4 timeseries file (L0B)
setwd(path_to_L0b)
myfilename <- paste(subsite_ID, as.character(as.Date(first(myCH4data_manID$DATE))),analyser,"CH4",sep="_")
myfilename <- paste0(myfilename, ".csv")
write.csv(x = myCH4data_manID, file = myfilename, sep = ";", dec = ".", row.names = F, col.names = T)


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
                          plot.lim = c(1900,max(fieldsheet$final_ch4)*1000)) %>%
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





#----- joining CO2 and CH4 fluxes estimates into a single table -----


table_results <- auxfile[,-c(which(names(auxfile) =="c0"),which(names(auxfile) =="cf"))] %>%
  left_join(CO2_flux_res %>% select(UniqueID, best.flux, model, quality.check)) %>%
  rename(CO2_flux = best.flux, CO2_model = model, CO2_quality.check = quality.check) %>%
  left_join(CH4_flux_res %>% select(UniqueID, best.flux, model, quality.check)) %>%
  rename(CH4_diffusive_flux = best.flux, CH4_model = model, CH4_quality.check = quality.check) %>%
  left_join(table_results_CH4 %>% select(UniqueID, CH4_ebullition)) %>%
  rename(CH4_ebullition_flux = CH4_ebullition)


plt_CO2 <- ggplot(table_results, aes(lightCondition, CO2_flux, fill = lightCondition))+
  geom_hline(yintercept = 0)+
  geom_boxplot(alpha=0.2)+geom_jitter(width = 0.2)+
  theme_article()+facet_wrap(.~strata, scales = "free")+
  ylab("CO2 flux mmol/m2/s")+
  ggtitle(paste0(subsite_ID,", CO2 flux"))+
  scale_fill_viridis_d(begin = 0.2, end = 0.9)


plt_CH4diff <- ggplot(table_results, aes(lightCondition, CH4_diffusive_flux, fill = lightCondition))+
  geom_hline(yintercept = 0)+
  geom_boxplot(alpha=0.2)+geom_jitter(width = 0.2)+
  theme_article()+facet_wrap(.~strata, scales = "free")+
  ylab("CH4 flux nmol/m2/s")+
  ggtitle(paste0(subsite_ID,", CH4 diffusive flux"))+
  scale_fill_viridis_d(begin = 0.2, end = 0.9)


plt_CH4ebull <- ggplot(table_results, aes(lightCondition, CH4_ebullition_flux, fill = lightCondition))+
  geom_hline(yintercept = 0)+
  geom_boxplot(alpha=0.2)+geom_jitter(width = 0.2)+
  theme_article()+facet_wrap(.~strata, scales = "free")+
  ylab("CH4 flux nmol/m2/s")+
  ggtitle(paste0(subsite_ID,", CH4 ebullition"))+
  scale_fill_viridis_d(begin = 0.2, end = 0.9)

ggarrange(plt_CO2, plt_CH4diff, plt_CH4ebull, ncol = 1)




# saving fluxes estimates
setwd(path_to_L0b)
myfilename <- paste(subsite_ID, as.character(as.Date(first(table_results$start.time))),analyser,"fluxes",sep="_")
myfilename <- paste0(myfilename, ".csv")
write.csv(x = table_results, file = myfilename, sep = ";", dec = ".", row.names = F, col.names = T)



#----- joining CO2 and CH4 fluxes estimates into a single table -----


# Combine plot lists into one list
flux_plot.ls <- c(CO2_flux_plots, CH4_flux_plots, H2O_flux_plots)

# Save plots to pdf
plot_path <- paste0(path_to_L0b,"/plots_",subsite_ID)

dir.create(plot_path)
setwd(plot_path)
myfilename <- paste(subsite_ID, as.character(as.Date(first(table_results$start.time))),analyser,sep="_")
myfilename <- gsub("-","", myfilename)
flux2pdf(flux_plot.ls, outfile = paste0(myfilename,".pdf"))


