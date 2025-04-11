# ---
# Authors: Camille Minaudo
# Project: "RESTORE4Cs"
# date: "Avril 2025"

# ---

# --- Description of this script
# 


# clearing workspace and console
rm(list = ls()) # clear workspace
cat("/014") # clear console


library(lubridate)
library(tidyr)
library(dplyr)
library(pbapply)
library(ggplot2)
library(egg)
library(goFlux)
library(purrr)
library(ggnewscale)
library(stringr)
library(data.table)


# make sure you clone or download the latest version of aquaGHG: https://github.com/camilleminaudo/aquaGHG/tree/main 
# then, correctly specify the path to the corresponding folder and source all scripts from there:
files.sources = list.files(path = "C:/Projects/myGit/aquaGHG/R/", full.names = T)
for (f in files.sources){source(f)}


# ---- Directories ----
dropbox_root <- "C:/Users/Camille Minaudo/Dropbox/RESTORE4Cs - Fieldwork/Data" 

# datapath <- paste0(dropbox_root,"/GHG/RAW data")
RData_path <- paste0(dropbox_root, "/GHG//Processed data/RData") # be careful with the data in the Dropbox. On 24/03/2025, message sent to Miguel about it
results_path <- paste0(dropbox_root,"/GHG/Processed data/computed_flux/")
auxfile_path <- "SET here"

# ---- Loading auxfile ----
setwd(auxfile_path)
auxfile <- read.csv(file = "auxfile.csv")
auxfile <- auxfile %>% mutate(start.time = as.POSIXct(round(start.time,"secs")))
auxfile$obs.length <- floor(auxfile$duration)


# ---- Re-compute fluxes ----

load_incubation <- function(auxfile_i, RData_path){
  
  message("Loading data for ",auxfile_i$UniqueID)
  gas <- unique(auxfile_i$gas_analiser)
  
  setwd(RData_path)
  if(gas== "LI-COR"){
    gs_suffix <- "LI-7810"
  } else {
    gs_suffix <- gas
  }
  load(file = paste0(auxfile_i$subsite,"_",gs_suffix,".RData"))
  mydata <- mydata[,c("POSIX.time", 
                      "CO2dry_ppm", "CH4dry_ppb", "H2O_ppm",
                      "CO2_prec",   "CH4_prec",   "H2O_prec"  )]
  
  mydata <- mydata[which(mydata$POSIX.time>=auxfile_i$start.time & 
                           mydata$POSIX.time<=auxfile_i$start.time+auxfile_i$duration),]
  
  if(dim(mydata)[1]>0){
    mydata$Etime <- as.numeric(mydata$POSIX.time) - min(as.numeric(mydata$POSIX.time))
    mydata$UniqueID <- auxfile_i$UniqueID
  } else {
    warning(paste0("For ", auxfile_i$UniqueID, ", no measurements were found!"))
  }
  
  return(mydata)
}


CO2_flux.auto <- CH4_flux.auto <- df.no_measurements <- NULL
list_ids <- unique(auxfile$UniqueID)
for (k in seq_along(list_ids)){
  i = list_ids[k]
  message(paste0("processing ",i))
  
  auxfile_i <- auxfile[which(auxfile$UniqueID==i),]
  
  if(dim(auxfile_i)[1]==0){
    message(paste0("Could not find corresponding auxfile for ",i))
    df.no_measurements <- rbind(df.no_measurements,
                                data.frame(UniqueID=i,
                                           message = "no auxfile"))
  } else {
    mydata <- load_incubation(auxfile_i, RData_path)
    
    if(dim(mydata)[1]==0){
      df.no_measurements <- rbind(df.no_measurements,
                                  data.frame(UniqueID=i,
                                             message = "no measurements"))
    } else {
      
      # print(plot_incubations(dataframe = mydata))
      
      # message("... calculate CO2 flux with aquaGHG, no manual selection")
      CO2_flux.auto_i <- automaticflux(dataframe = mydata, myauxfile = auxfile_i, shoulder = 0, gastype = "CO2dry_ppm", 
                                       fluxSeparation = FALSE, displayPlots = FALSE, method = "trust.it.all")
      # message("... calculate CH4 flux with aquaGHG, no manual selection")
      CH4_flux.auto_i <- automaticflux(dataframe = mydata, myauxfile = auxfile_i, shoulder = 0, gastype = "CH4dry_ppb", 
                                       fluxSeparation = T, displayPlots = FALSE, method = "trust.it.all")
      
      CO2_flux.auto <- rbind(CO2_flux.auto, CO2_flux.auto_i)
      CH4_flux.auto <- rbind(CH4_flux.auto, CH4_flux.auto_i)
      
    }
    rm(auxfile_i)
    }
}

# Save files

setwd(results_path)
save(list = c("CO2_flux.auto",
              "CH4_flux.auto",
              "df.no_measurements"), file = "results_aquaGHG_recalculated.Rdata")






