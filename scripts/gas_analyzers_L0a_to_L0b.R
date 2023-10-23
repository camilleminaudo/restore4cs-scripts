
# ---
# Authors: Camille Minaudo, Benjamin Misteli
# Project: "RESTORE4Cs"
# date: "Oct 2023"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
# This script loads a raw measurement file (data level L0a) from one of the gas
# analyzers used in the project, and transform it into a unified harmonized csv
# file (data level L0b) allowing for visualization of the data and further data processing.


rm(list = ls()) # clear workspace
cat("/014") # clear console


# packages
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



# Directories
dropbox_root <- "C:/Users/Camille Minaudo/Dropbox/RESTORE4Cs - Fieldwork/Data"
datapath <- paste0(dropbox_root,"/GHG/RAW data")
fieldsheetpath <- paste0(dropbox_root,"/GHG/Fieldsheets")
loggerspath <- paste0(datapath,"/RAW Data Logger")

setwd(datapath)

# --- SETTINGS
analyser <- "Licor"
site_ID <- "S1-CU"
subsite_ID <- "S1-CU-P2"
file_to_read <- "S1-CU-P2.data"
who_runs_this <- "Camille Minaudo"

# --- Some functions
read_Licor <- function(file){
  message(paste0("reading ",file_to_read," with script for ",analyser," gas analyser"))
  filename <- paste(datapath,directory_analyser,file_to_read, sep = "/")
  data_raw <- read_lines(filename)
  prefix <- substr(data_raw, start = 1,stop = 5) # isolte first 5 characters for each line

  # find line corresponding to headers
  headers <- unlist(strsplit(data_raw[which(prefix == "DATAH")], "\t"))
  units <- unlist(strsplit(data_raw[which(prefix == "DATAU")], "\t"))

  data <- read.delim(filename, sep = "\t", header = F, skip = which(prefix == "DATAU"), na.strings = "nan")
  names(data) <- headers

  my_data <- data.frame(date = data$DATE,
                        UTCtime = data$TIME,
                        unixtime = data$SECONDS,
                        H2O = data$H2O,
                        CO2 = data$CO2,
                        CH4 = data$CH4/1000, #ppm
                        Press = data$CAVITY_P,
                        label = data$REMARK)

  return(my_data)
}



# Read gas analyser's file
if(analyser == "Licor"){
  directory_analyser <- "RAW Data Licor-7810"
  my_data <- read_Licor(file = paste(datapath,directory_analyser,file_to_read, sep = "/"))
} else if(analyser == "Los Gatos"){
  directory_analyser <- "RAW Data Los Gatos"
} else if (analyser == "Picarro"){
  directory_analyser <- "RAW Data Picarro"
}


# select only rows with no CO2 NA values
# my_data <- my_data[!is.na(my_data$CO2),]


# Read corresponding Fieldsheet
path2file <- paste0(fieldsheetpath,"/",site_ID,"/",subsite_ID,"-Fieldsheet-GHG.xlsx")
fieldsheet_temp <- readxl::read_xlsx(path2file,
                                     col_names = T)
fieldsheet <- readxl::read_xlsx(path2file,
                                skip = 2, col_names = F)
names(fieldsheet) <- names(fieldsheet_temp)

get_unix_times <- function(mydate, mytime){
  timestamp <- paste(hour(mytime),minute(mytime),0,sep = ":")
  unix_time <- as.numeric(as.POSIXct(paste(as.Date(mydate),timestamp,
                                           sep = " "), tz = 'UTC'))
  return(unix_time)
}

fieldsheet$unix_start_time <- get_unix_times(mydate = fieldsheet$date, mytime = fieldsheet$start_time)
fieldsheet$unix_end_time <- get_unix_times(mydate = fieldsheet$date, mytime = fieldsheet$end_time)
# head(fieldsheet)


# Read corresponding Loggers data
SN_logger_float <- unique(fieldsheet$logger_floating_chamber)
SN_logger_tube <- unique(fieldsheet$logger_transparent_chamber)

if(SN_logger_float != "NA"){
  is_data_logger_float = T
  data_logger_float <- readxl::read_xlsx(paste0(loggerspath,"/",site_ID,"-",SN_logger_float,".xlsx"),col_names = T)
  data_logger_float$unixtime <- as.numeric(data_logger_float$`Date/hour (UTC)`)
} else {
  message("no data logger linked to the floating chamber!")
  is_data_logger_float = F}

if(SN_logger_tube != "NA"){
  is_data_logger_tube = T
  data_logger_tube <- readxl::read_xlsx(paste0(loggerspath,"/",site_ID,"-",SN_logger_tube,".xlsx"),col_names = T)
  data_logger_tube$unixtime <- as.numeric(data_logger_tube$`Date/hour (UTC)`)
} else {
  is_data_logger_tube = F
  message("no data logger linked to the transparent tube chamber!")
  if(SN_logger_float != "NA"){
    message(".../ using data from floating chamber instead...")
    data_logger_tube = data_logger_float
    is_data_logger_tube = T
  }
}

#
# for (i in  seq(1,length(fieldsheet$pilot_site))){ # for each incubation, proceed with...
#
  # my_sel <- my_data[my_data$unixtime>= (fieldsheet$unix_start_time[i]-30) & my_data$unixtime<= (fieldsheet$unix_end_time[i]+60),]
  #
  # my_label <- unique(my_sel$label)
  # my_label <- my_label[which(my_label != "")]
  #
  # pCO2 <- ggplot(my_sel, aes(unixtime, CO2))+geom_line()+
  #   geom_line(data = my_data[my_data$label == my_label,],
  #             aes(unixtime, CO2), colour = "red", alpha = 0.5, linewidth = 1.5)+
  #   ylab("CO2 [ppm]")+
  #   ggtitle(paste0(subsite_ID," ",fieldsheet$date[i],", plot ",fieldsheet$plot_id[i],", incub ", i))+
  #   theme_article()
  # pCH4 <- ggplot(my_sel, aes(unixtime, CH4))+geom_line()+
  #   geom_line(data = my_data[my_data$label == my_label,],
  #             aes(unixtime, CH4), colour = "red", alpha = 0.5, linewidth = 1.5)+
  #   ylab("CH4 [ppm]")+
  #   ggtitle(paste(fieldsheet$chamber_type[i], fieldsheet$strata[i], fieldsheet$transparent_dark[i], sep=", "))+
  #   theme_article()
  #
  # p <- ggarrange(pCO2,pCH4, nrow = 1)
#
# }


# load Li-COR file with GoFluxYourself package

mydata_imp <- LI7810_import(inputfile = paste(datapath,directory_analyser,file_to_read, sep = "/"))



# The auxfile requires start.time and UniqueID
# start.time must be in the format "%Y-%m-%d %H:%M:%S"
auxfile <- NULL
# for (i in 15){
for (i in seq_along(fieldsheet$pilot_site)){

  my_sel <- my_data[my_data$unixtime>= (fieldsheet$unix_start_time[i]) & my_data$unixtime<= (fieldsheet$unix_end_time[i]),]
  my_label <- unique(my_sel$label)
  my_label <- my_label[which(my_label != "")]
  my_sel <- my_data[my_data$label == my_label,]

  if (fieldsheet$chamber_type[i] == "floating"){
    my_sel$temperature <- approx(data_logger_float$unixtime, data_logger_float$`Ch:1 - Temperature   (°C)`, xout = my_sel$unixtime)$y
    myArea = 14365.4439 # cm2
    myVtot = 115 # L
  } else if (fieldsheet$chamber_type[i] == "tube"){
    my_sel$temperature <- approx(data_logger_tube$unixtime, data_logger_tube$`Ch:1 - Temperature   (°C)`, xout = my_sel$unixtime)$y
    myArea = pi*12.1**2 # cm2
    myVtot = myArea*fieldsheet$chamber_height_cm[i]*1e-3 # L
  } else {
    warning("chamber type not correct!")
  }
  auxfile_tmp <- data.frame(UniqueID = paste("plot",fieldsheet$plot_id[i],"_",fieldsheet$chamber_type[i],"_",fieldsheet$transparent_dark[i],
                                             sep = ""),
                            start.time = as.POSIXct((fieldsheet$unix_start_time[i]-30), tz = "UTC"),
                            duration = (fieldsheet$unix_end_time[i]+30) - (fieldsheet$unix_start_time[i]-30),
                            Area = myArea,
                            Vtot = myVtot,
                            Tcham = mean(my_sel$temperature),
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
mydata_manID <- lapply(seq_along(mydata_ow), click.peak.loop,
                     flux.unique = mydata_ow,
                     gastype = "CO2dry_ppm",
                     plot.lim = c(200,1000)) %>%
  map_df(., ~as.data.frame(.x))

# Additional auxiliary data required for flux calculation.
mydata_manID <- mydata_manID %>%
  left_join(auxfile %>% select(UniqueID, Area, Vtot, Tcham, Pcham))


# Calculate fluxes for CO2 and H2O
CO2_results <- goFlux(mydata_manID, "CO2dry_ppm")
H2O_results <- goFlux(mydata_manID, "H2O_ppm")

# Use best.flux to select the best flux estimates (LM or HM)
# based on a list of criteria
criteria <- c("g.factor", "kappa", "MDF", "R2", "SE.rel")

CO2_flux_res <- best.flux(CO2_results, criteria)
H2O_flux_res <- best.flux(H2O_results, criteria)

CO2_flux_res <- CO2_flux_res %>%
  left_join(auxfile %>% select(UniqueID, strata, chamberType, lightCondition))


# Plots results

# Make a list of plots of all measurements, for each gastype
CO2_flux_plots <- flux.plot(CO2_flux_res, mydata_manID, "CO2dry_ppm")
H2O_flux_plots <- flux.plot(H2O_flux_res, mydata_manID, "H2O_ppm")

# Combine plot lists into one list
flux_plot.ls <- c(CO2_flux_plots, H2O_flux_plots)

# Save plots to pdf
# flux2pdf(flux_plot.ls, outfile = "demo.results.pdf")

ggplot(CO2_flux_res, aes(lightCondition, best.flux, fill = lightCondition))+
  geom_boxplot(alpha=0.2)+geom_jitter(width = 0.2)+
  theme_article()+facet_wrap(.~strata, scales = "free")+
  scale_fill_viridis_d(begin = 0.2, end = 0.9)










#----- computing CH4 fluxes -----


# Define the measurements' window of observation
mydata_ow <- obs.win(inputfile = mydata_imp, auxfile = auxfile,
                     obs.length = auxfile$duration, shoulder = 30)


# Manually identify start/end CH4 measurements by clicking on the start and end points
mydata_manID <- lapply(seq_along(mydata_ow), click.peak.loop,
                       flux.unique = mydata_ow,
                       gastype = "CH4dry_ppb",
                       plot.lim = c(1900,4000)) %>%
  map_df(., ~as.data.frame(.x))

# linking window of start/end time to mydata_ow list
i=0
for (id in unique(mydata_manID$UniqueID)){
  i=i+1
  mydata_ow[[i]]$start.time <- unique(mydata_manID$start.time_corr[mydata_manID$UniqueID == id])
  mydata_ow[[i]]$end.time <- unique(mydata_manID$end.time[mydata_manID$UniqueID == id])
  mydata_ow[[i]]$duration <- as.numeric(mydata_ow[[i]]$end.time) - as.numeric(mydata_ow[[i]]$start.time)
  mydata_ow[[i]]$c0 <- first(mydata_manID$c0[mydata_manID$UniqueID == id])

  auxfile$c0[i] <- first(mydata_manID$c0[mydata_manID$UniqueID == id])
  auxfile$cf[i] <- first(mydata_manID$cf[mydata_manID$UniqueID == id])
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
CH4_flux_plots


table_results <- auxfile %>%
  left_join(CH4_flux_res %>% select(UniqueID, best.flux, model))



# Estimating ebullitive component
i=0
for (id in unique(mydata_manID$UniqueID)){
  i=i+1

  CH4_final <- table_results$cf[i]
  CH4_initial <-  table_results$c0[i]
  incubation_time <- first(mydata_ow[[i]]$duration)

  H2O_mol = mydata_ow[[i]]$H2O_ppm / (1000*1000)
  myfluxterm <- flux.term(table_results$Vtot[i], table_results$Pcham[i], table_results$Area[i],
            table_results$Tcham[i], first(H2O_mol))

  CH4_flux_total <- (CH4_final-CH4_initial)/incubation_time*myfluxterm # ppb/m2/s


  CH4_ebullition <- CH4_flux_total - CH4_flux_res$best.flux[i] # total flux - diffusive term
  CH4_ebullition[CH4_ebullition<0] <- 0

  table_results$CH4_ebullition[i] <- CH4_ebullition
}


