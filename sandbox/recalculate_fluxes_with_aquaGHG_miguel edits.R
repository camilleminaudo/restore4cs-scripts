# ---
# Authors: Camille Minaudo, Miguel Cabrera
# Project: "RESTORE4Cs"
# date: "June 2025"

# ---

# --- Description----
#this script uses the functions from aquaGHG to calculate CO2 and CH4 fluxes for all incubations separately. Takes auxfiles created with script Create_CO2andCH4_auxfiles_4aquaGHG.R and Rdata updated with script Harmonize_GHG_raw2Rdata.R (accurate timezones and correct units, last update 6/5/2025)

#separated loops for CO2 and CH4 (each can have different start.time and duration for a given UniqueID)

#Needs for restore4cs: 
#CO2: no need to calulate fluxes with aquaGHG, all results have been correctly produced and model chosen with previous raw2flux scripts. Done here for consistency and to test. 

#CH4: We run independent aquaGHG loops (with and W/O flux separation), irrespective of water presence (there are ebullitive dynamics without water). Improvements needed for flux separation algorithm.

#3 methods for each incubation: LM, HM, total.flux. 
#total.flux has been modified to be calcualted using the mean concentration difference (first vs last 10s, and the elapsed time between these means). THIS ensures an accurate & appropriate estimate (albeit slightly underestimated when linear patterns in windows) associated with a reliable uncertainty estimate propagated SE (need to use SE to be comparable with uncertainties of LM/HM derived flux SE estimates). For deltaconc use Cf - C0 propagate SE, for deltatime use elapsed time between center of initial and final time-windows (duration - 10s).


#TO-DO------
#Improve flux separation: 
  #1. Improve ebullition detection
  #2. Improve selection of difusive chunk
  #3. Set criteria for quality of separation


#BUGs and issues------
# #SOME ISSUES WITH start.time class (rounding causes class change to POSIXt), easy fix by re-transforming to POSIXct

#BUG: total.flux within aquaGHG calculates total flux with the incubation cropped until the start of the identified difusive chunk. Total flux should include the whole incubation, irregardless of which section is identified as the diffusive chunk. In this script, corrected variables are added to the CH4water_flux results:
# full.total.flux (and SD), full.ebullition.flux (using full.total.flux - difusive.flux)


#BUG: flux.plot fails when HM.flux is NA (breaking the loop): 
#Error en flux.plot(flux.results = best.flux_auto, dataframe = mydata_auto, : 
# 'HM.flux' in 'flux.results' must be of class numeric

#FIX: inside automaticflux(), force all outputs of HM model to numeric, even when they are NA.

#IN ADDITION, UniqueID not displayed in plots without flux separation. 
#Pseudofix: manually add plottitle in loop using grid::grid.text()
#I've tried to add print(p + ggtitle) within automaticflux, but apparently p is not of class ggplot
#With flux separation, plots are displayed without issue and with correct title.

#ISSUE/DOUBT: How to modify the criteria for best flux using the wrapper of aquaGHG? (not essential can perform selection after the flux calculation run)


# clearing workspace and console
rm(list = ls()) # clear workspace
cat("/014") # clear console



# ---- Directories ----
dropbox_root <- "C:/Users/Miguel/Dropbox/RESTORE4Cs - Fieldwork/Data" # You have to make sure this is pointing to the write folder on your local machine

#Path to Rdata, udated and fixed by script "Harmonize_GHG_raw2Rdata.R" UTCtimes and correct units (as of 6/5/2025).
RData_path <- paste0(dropbox_root, "/GHG/Processed data/RData") 

#Path to Co2 and Ch4 auxfiles with corrected start.time and duration (cropping implemented)
auxfile_path<- paste0(dropbox_root,"/GHG/Working data/Auxfiles_4aquaGHG/") 

#Set results_path Within R4Cs dropbox
results_path <- paste0(dropbox_root,"/GHG/Processed data/computed_flux/")

plots_path<- paste0(results_path,"aquaGHG_incubation_plots/")
if (!dir.exists(plots_path)) {
  dir.create(plots_path, recursive = TRUE)
}


#Load PKGs and functions-----
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
repo_root_r4cs <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
files.sources = list.files(path = paste0(dirname(repo_root_r4cs),"/aquaGHG/R/"), full.names = T)
for (f in files.sources){source(f)}

#Load function
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



# ---- Loading auxfiles ----

#CO2 auxfile:
co2_auxfile <- read.csv(file = paste0(auxfile_path,"co2_auxfile.csv"))

co2_auxfile <- co2_auxfile %>% 
  mutate(start.time = as.POSIXct(round(as.POSIXct(start.time,tz = "UTC"),"secs")),#rounding transforms to class different to POSIXct, causing error, forze as.POSIXct to fix
         obs.length=floor(duration))

#CH4 auxfile:
ch4_auxfile <- read.csv(file = paste0(auxfile_path,"ch4_auxfile.csv"))

ch4_auxfile <- ch4_auxfile %>% 
  mutate(start.time = as.POSIXct(round(as.POSIXct(start.time,tz = "UTC"),"secs")),#rounding transforms to class different to POSIXct, causing error, forze as.POSIXct to fix
         obs.length=floor(duration))


#-----subset to test (Optional)------

#select only a couple of subsites to test the script
pattern<- "s"

ch4_auxfile<- ch4_auxfile %>% 
  filter(grepl(pattern, UniqueID))

co2_auxfile<- co2_auxfile %>% 
  filter(grepl(pattern, UniqueID))


#-----Initialize objects----
CO2_flux.auto <- CH4all_nosep_flux.auto <- CH4all_fluxsep_flux.auto <- NULL

CO2_df.no_measurements <- CH4all_nosep_df.no_measurements <-  CH4all_fluxsep_df.no_measurements <- NULL

#----CO2 loop (pdf plots)-----
#This takes approx 30 minutes
rm(subsite, k)
for (subsite in unique(co2_auxfile$subsite)){
  message(paste0("processing subsite",subsite))
  
  #Start pdf for plots
  pdf(file = paste0(plots_path,"CO2_",subsite,".pdf"))
  
  #loop over uniqueID of current subsite
  for (k in seq_along(unique(co2_auxfile[which(co2_auxfile$subsite==subsite),]$UniqueID))){
    
    #Select uniqueID and display progress
    i = co2_auxfile[which(co2_auxfile$subsite==subsite),]$UniqueID[k]
    
    incubnum<- which(co2_auxfile$UniqueID==i)
    message(paste0("Processing CO2 incubation ",incubnum, " of ",length(co2_auxfile$UniqueID)," (",round(100*incubnum/length(co2_auxfile$UniqueID),0), " %)"))
    
    #Load auxfiles for UniqueID
    co2_auxfile_i <- co2_auxfile[which(co2_auxfile$UniqueID==i),]
    
    #check that UniqueID has auxfile
    if(dim(co2_auxfile_i)[1]==0){
      message(paste0("Could not find corresponding auxfile for ",i))
      CO2_df.no_measurements <- rbind(CO2_df.no_measurements,
                                           data.frame(UniqueID=i,
                                                      message = "no auxfile"))
    } else {
      #Load data for UniqueID
      mydata_co2 <- load_incubation(co2_auxfile_i, RData_path)
      
      #Check that there is data for UniqueID
      if(dim(mydata_co2)[1]==0){
        CO2_df.no_measurements <- rbind(CO2_df.no_measurements,
                                             data.frame(UniqueID=i,
                                                        message = "no measurements"))
      } else {
        #Calculate flux for UniqueID and plot    
        CO2_flux.auto_i <- automaticflux(dataframe = mydata_co2, myauxfile = co2_auxfile_i, shoulder = 0, gastype = "CO2dry_ppm", 
                                              fluxSeparation = F,
                                              displayPlots = TRUE, 
                                              method = "trust.it.all")
        
        #ADD the uniqueID to the displayed plot (to save it in pdf): 
        grid::grid.text(as.character(i),
                        x = 0.5, y = 0.975,
                        gp = grid::gpar(fontsize = 14, fontface = "bold"))
        
        ##ADD full.total.flux------
        #ADD full.total.flux with whole incubation period: 
        {
          #Calculate total.flux and total.flux.SD for whole incubation. Using chunks copied from aquaGHG:
          
          #smoothing signal:
          # smooth signal from incubation (whole duration)
          mydf <- data.frame(POSIX.time = mydata_co2$POSIX.time,
                             time = as.numeric(mydata_co2$POSIX.time-first(mydata_co2$POSIX.time)),
                             conc = mydata_co2[["CO2dry_ppm"]])
          
          concsmooth <- smooth.spline(x = mydf$time, y = mydf$conc, nknots = round(length(mydf$time)/3), spar = 0.8)
          mydf$concsmooth <- approx(concsmooth$x, concsmooth$y, xout = mydf$time, rule = 2)$y
          
          
          #get initial, final concentration and incubation time (whole duration)
          t.win <- 30
          C0 = min(mydf$concsmooth[mydf$time<t.win])
          Cf = max(mydf$concsmooth[mydf$time>max(mydf$time)-t.win])
          incubation_time = last(mydf$time)
          #get deltaconcs
          deltaconcs = Cf-C0
          #Calculate SD for initial final and deltaconc
          SD_C0 <- sd(mydf$conc[mydf$time<t.win])
          SD_Cf <- sd(mydf$conc[mydf$time>max(mydf$time)-t.win])
          SD_deltaconcs <- sqrt(SD_C0^2+SD_Cf^2)
          
          #Calculate full.total.flux aplying flux.term from aquaGHG
          full.total.flux <- (deltaconcs)/incubation_time*CO2_flux.auto_i$flux.term # umol/m2/s
          #Calculate SD of total flux
          SD_full.total.flux <- abs(full.total.flux) * SD_deltaconcs/deltaconcs
          
          #ADD full.total.flux to flux.auto_i dataframe
          CO2_flux.auto_i$full.total.flux<- full.total.flux
          CO2_flux.auto_i$SD_full.total.flux<- SD_full.total.flux
        }
        
        #ADD mean.total.flux------
        {
          #calculate mean.C0, mean.C0.sd, mean.Cf, mean.Cf.sd, mean.total.flux, mean.total.flux.SD, mean.total.flux.SE
          #THe approach used in aquaGHG (smoothing--> minmax difference) cannot be properly assigned to an uncertainty measure and its a bit too complex to provide a clear definition to readers/reviewers. Here we calculate instead the mean.total.flux (along mean.C0 and mean.Cf),using a 10s time window and raw concentrations. Uncertainty is easily propagated onto SD and SE (SE needed for consistency with LM.flux and HM.flux)
          
          #get average initial, final concentration and incubation time (whole duration, time-window 10s)
          t.win <- 10
          mean.C0 <- mean(mydf$conc[mydf$time<t.win])
          mean.Cf <- mean(mydf$conc[mydf$time>max(mydf$time)-t.win])
          #Correct incubation_time to be calcualted between central distribution of initial and final time windows. (avoids underestimation when we have), duration-(2*(t.win/2))= duration-t.win
          mean.incubation_time <- max(mydf$time)-t.win
          #get deltaconcs using mean C0 and Cf
          mean.deltaconcs <- mean.Cf-mean.C0
          
          #Calculate SDs for initial final and deltaconc
          mean.C0.SD <- sd(mydf$conc[mydf$time<t.win])
          mean.Cf.SD <- sd(mydf$conc[mydf$time>max(mydf$time)-t.win])
          mean.deltaconcs.SD <- sqrt(mean.C0.SD^2+mean.Cf.SD^2)
          # Calculate SE for mean.deltaconc (using t.win as n for each mean conc)
          mean.deltaconcs.SE <- sqrt((mean.C0.SD^2 / t.win) + (mean.Cf.SD^2 / t.win))
          
          #Calculate mean.total.flux (and associated uncertainties) by aplying flux.term from aquaGHG
          mean.total.flux <- (mean.deltaconcs)/mean.incubation_time*CO2_flux.auto_i$flux.term # umol/m2/s
          #Calculate SD of mean total flux
          mean.total.flux.SD <- abs(mean.total.flux) * mean.deltaconcs.SD/mean.deltaconcs
          # Calculate SE of mean total flux (error propagation)
          mean.total.flux.SE <- (CO2_flux.auto_i$flux.term / mean.incubation_time) * mean.deltaconcs.SE
          
          #ADD mean.XX variables to flux.auto_i dataframe
          CO2_flux.auto_i$mean.C0<- mean.C0
          CO2_flux.auto_i$mean.C0.SD<- mean.C0.SD
          CO2_flux.auto_i$mean.Cf<- mean.Cf
          CO2_flux.auto_i$mean.Cf.SD<- mean.Cf.SD
          
          CO2_flux.auto_i$mean.total.flux <- mean.total.flux
          CO2_flux.auto_i$mean.total.flux.SD <- mean.total.flux.SD
          CO2_flux.auto_i$mean.total.flux.SE <- mean.total.flux.SE
          }
        
        #Join flux to rest of dataset
        CO2_flux.auto <- rbind(CO2_flux.auto, CO2_flux.auto_i)
      }
    }
    
    #Remove auxfile for UniqueID (avoids re-usage, if next does not exist)
    rm(co2_auxfile_i)
  }#End of UniqueID loop
  dev.off()
}#end subsite loop


# Save results to Rdata
save(list = c("CO2_flux.auto",
              "CO2_df.no_measurements"), file = paste0(results_path,"CO2_aquaGHG.Rdata"))




#----CH4all_nosep loop (pdf plots)-----
#This takes approx 30 minutes
rm(subsite, k)
for (subsite in unique(ch4_auxfile$subsite)){
  message(paste0("processing subsite",subsite))
  
  #Start pdf for plots
  pdf(file = paste0(plots_path,"CH4all_nosep_",subsite,".pdf"))
  
  #loop over uniqueID of current subsite
  for (k in seq_along(unique(ch4_auxfile[which(ch4_auxfile$subsite==subsite),]$UniqueID))){
    
    #Select uniqueID and display progress
    i = ch4_auxfile[which(ch4_auxfile$subsite==subsite),]$UniqueID[k]
    
    incubnum<- which(ch4_auxfile$UniqueID==i)
    message(paste0("Processing CH4 incubation (all,no sep) ",incubnum, " of ",length(ch4_auxfile$UniqueID)," (",round(100*incubnum/length(ch4_auxfile$UniqueID),0), " %)"))
    
    #Load auxfiles for UniqueID
    ch4_auxfile_i <- ch4_auxfile[which(ch4_auxfile$UniqueID==i),]
    
    #check that UniqueID has auxfile
    if(dim(ch4_auxfile_i)[1]==0){
      message(paste0("Could not find corresponding auxfile for ",i))
      CH4all_nosep_df.no_measurements <- rbind(CH4all_nosep_df.no_measurements,
                                      data.frame(UniqueID=i,
                                                 message = "no auxfile"))
    } else {
      #Load data for UniqueID
      mydata_ch4 <- load_incubation(ch4_auxfile_i, RData_path)
      
      #Check that there is data for UniqueID
      if(dim(mydata_ch4)[1]==0){
        CH4all_nosep_df.no_measurements <- rbind(CH4all_nosep_df.no_measurements,
                                        data.frame(UniqueID=i,
                                                   message = "no measurements"))
      } else {
        #Calculate flux for UniqueID and plot    
        CH4all_nosep_flux.auto_i <- automaticflux(dataframe = mydata_ch4, myauxfile = ch4_auxfile_i, shoulder = 0, gastype = "CH4dry_ppb", 
                                         fluxSeparation = F,
                                         displayPlots = TRUE, 
                                         method = "trust.it.all")
        
        #ADD the uniqueID to the displayed plot (to save it in pdf): 
        grid::grid.text(as.character(i),
                        x = 0.5, y = 0.975,
                        gp = grid::gpar(fontsize = 14, fontface = "bold"))
        
        
        ##ADD full.total.flux------
        #ADD full.total.flux with whole incubation period:
        {
          #Calculate total.flux and total.flux.SD for whole incubation. Using chunks copied from aquaGHG:
          
          #smoothing signal:
          # smooth signal from incubation (whole duration)
          mydf <- data.frame(POSIX.time = mydata_ch4$POSIX.time,
                             time = as.numeric(mydata_ch4$POSIX.time-first(mydata_ch4$POSIX.time)),
                             conc = mydata_ch4[["CH4dry_ppb"]])
          
          concsmooth <- smooth.spline(x = mydf$time, y = mydf$conc, nknots = round(length(mydf$time)/3), spar = 0.8)
          mydf$concsmooth <- approx(concsmooth$x, concsmooth$y, xout = mydf$time, rule = 2)$y
          
          #get initial, final concentration and incubation time (whole duration)
          t.win <- 30
          C0 = min(mydf$concsmooth[mydf$time<t.win])
          Cf = max(mydf$concsmooth[mydf$time>max(mydf$time)-t.win])
          incubation_time = last(mydf$time)
          #get deltaconcs
          deltaconcs = Cf-C0
          #Calculate SD for initial final and deltaconc
          SD_C0 <- sd(mydf$conc[mydf$time<t.win])
          SD_Cf <- sd(mydf$conc[mydf$time>max(mydf$time)-t.win])
          SD_deltaconcs <- sqrt(SD_C0^2+SD_Cf^2)
          
          #Calculate full.total.flux aplying flux.term from aquaGHG
          full.total.flux <- (deltaconcs)/incubation_time*CH4all_nosep_flux.auto_i$flux.term # nmol/m2/s
          #Calculate SD of total flux
          SD_full.total.flux <- abs(full.total.flux) * SD_deltaconcs/deltaconcs
          
          #ADD full.total.flux to flux.auto_i dataframe
          CH4all_nosep_flux.auto_i$full.total.flux<- full.total.flux
          CH4all_nosep_flux.auto_i$SD_full.total.flux<- SD_full.total.flux
        }
        
        #ADD mean.total.flux------
        {
          #calculate mean.C0, mean.C0.sd, mean.Cf, mean.Cf.sd, mean.total.flux, mean.total.flux.SD, mean.total.flux.SE
          #THe approach used in aquaGHG (smoothing--> minmax difference) cannot be properly assigned to an uncertainty measure and its a bit too complex to provide a clear definition to readers/reviewers. Here we calculate instead the mean.total.flux (along mean.C0 and mean.Cf),using a 10s time window and raw concentrations. Uncertainty is easily propagated onto SD and SE (SE needed for consistency with LM.flux and HM.flux)
          
          #get average initial, final concentration and incubation time (whole duration, time-window 10s)
          t.win <- 10
          mean.C0 <- mean(mydf$conc[mydf$time<t.win])
          mean.Cf <- mean(mydf$conc[mydf$time>max(mydf$time)-t.win])
          #Correct incubation_time to be calcualted between central distribution of initial and final time windows. (avoids underestimation when we have), duration-(2*(t.win/2))= duration-t.win
          mean.incubation_time <- max(mydf$time)-t.win
          #get deltaconcs using mean C0 and Cf
          mean.deltaconcs <- mean.Cf-mean.C0
          
          #Calculate SDs for initial final and deltaconc
          mean.C0.SD <- sd(mydf$conc[mydf$time<t.win])
          mean.Cf.SD <- sd(mydf$conc[mydf$time>max(mydf$time)-t.win])
          mean.deltaconcs.SD <- sqrt(mean.C0.SD^2+mean.Cf.SD^2)
          # Calculate SE for mean.deltaconc (using t.win as n for each mean conc)
          mean.deltaconcs.SE <- sqrt((mean.C0.SD^2 / t.win) + (mean.Cf.SD^2 / t.win))
          
          #Calculate mean.total.flux (and associated uncertainties) by aplying flux.term from aquaGHG
          mean.total.flux <- (mean.deltaconcs)/mean.incubation_time*CH4all_nosep_flux.auto_i$flux.term # umol/m2/s
          #Calculate SD of mean total flux
          mean.total.flux.SD <- abs(mean.total.flux) * mean.deltaconcs.SD/mean.deltaconcs
          # Calculate SE of mean total flux (error propagation)
          mean.total.flux.SE <- (CH4all_nosep_flux.auto_i$flux.term / mean.incubation_time) * mean.deltaconcs.SE
          
          #ADD mean.XX variables to flux.auto_i dataframe
          CH4all_nosep_flux.auto_i$mean.C0<- mean.C0
          CH4all_nosep_flux.auto_i$mean.C0.SD<- mean.C0.SD
          CH4all_nosep_flux.auto_i$mean.Cf<- mean.Cf
          CH4all_nosep_flux.auto_i$mean.Cf.SD<- mean.Cf.SD
          
          CH4all_nosep_flux.auto_i$mean.total.flux <- mean.total.flux
          CH4all_nosep_flux.auto_i$mean.total.flux.SD <- mean.total.flux.SD
          CH4all_nosep_flux.auto_i$mean.total.flux.SE <- mean.total.flux.SE
        }
        
        
        #Join flux to rest of dataset
        CH4all_nosep_flux.auto <- rbind(CH4all_nosep_flux.auto, CH4all_nosep_flux.auto_i)
      }
    }
    
    #Remove auxfile for UniqueID (avoids re-usage, if next does not exist)
    rm(ch4_auxfile_i)
  }#End of UniqueID loop
  dev.off()
}#end subsite loop


# Save results to Rdata
save(list = c("CH4all_nosep_flux.auto",
              "CH4all_nosep_df.no_measurements"), file = paste0(results_path,"CH4all_nosep_aquaGHG.Rdata"))



#----CH4all flux separation(pdf plots)-----
#This takes approx XX minutes
rm(subsite, k)
for (subsite in unique(ch4_auxfile$subsite)){
  message(paste0("processing subsite",subsite))
  
  #Start pdf for plots
  pdf(file = paste0(plots_path,"CH4all_fluxsep_",subsite,".pdf"))
  
  #loop over uniqueID of current subsite
  for (k in seq_along(unique(ch4_auxfile[which(ch4_auxfile$subsite==subsite),]$UniqueID))){
    
    #Select uniqueID and display progress
    i = ch4_auxfile[which(ch4_auxfile$subsite==subsite),]$UniqueID[k]

    incubnum<- which(ch4_auxfile$UniqueID==i)
    message(paste0("Processing CH4 flux separation incubation ",incubnum, " of ",length(ch4_auxfile$UniqueID)," (",round(100*incubnum/length(ch4_auxfile$UniqueID),0), " %)"))
    
    #Load auxfiles for UniqueID
    ch4_auxfile_i <- ch4_auxfile[which(ch4_auxfile$UniqueID==i),]
   
    #check that UniqueID has auxfile for ch4
    if(dim(ch4_auxfile_i)[1]==0){
      message(paste0("Could not find corresponding auxfile for ",i))
      CH4all_fluxsep_df.no_measurements <- rbind(CH4all_fluxsep_df.no_measurements,
                                  data.frame(UniqueID=i,
                                             message = "no auxfile ch4"))
    } else {
      #Load ch4 data for UniqueID
      mydata_ch4 <- load_incubation(ch4_auxfile_i, RData_path)
      
      #Check that there is ch4 data for UniqueID
      if(dim(mydata_ch4)[1]==0){
        CH4all_fluxsep_df.no_measurements <- rbind(CH4all_fluxsep_df.no_measurements,
                                    data.frame(UniqueID=i,
                                               message = "no measurements for ch4"))
      } else {
        #Calculate CH4 flux for UniqueID  with flux separation and plot    
        CH4all_fluxsep_flux.auto_i <- automaticflux(dataframe = mydata_ch4, myauxfile = ch4_auxfile_i, shoulder = 0, gastype = "CH4dry_ppb", 
                                         fluxSeparation = T,
                                         displayPlots = TRUE, 
                                         method = "trust.it.all")
        
        ##full.total.flux & full.ebullition.flux------
        #ADD full.total.flux and full.ebullition.flux with whole incubation period:
        {
        #Calculate total.flux and total.flux.SD for whole incubation (to overwrite that coming from automaticflux, which only selects that after the start of identified difusion). Using chunks copied from aquaGHG:
        
        #smoothing signal:
        # smooth signal from incubation (whole duration)
        mydf <- data.frame(POSIX.time = mydata_ch4$POSIX.time,
                           time = as.numeric(mydata_ch4$POSIX.time-first(mydata_ch4$POSIX.time)),
                           conc = mydata_ch4[["CH4dry_ppb"]])
        
        concsmooth <- smooth.spline(x = mydf$time, y = mydf$conc, nknots = round(length(mydf$time)/3), spar = 0.8)
        mydf$concsmooth <- approx(concsmooth$x, concsmooth$y, xout = mydf$time, rule = 2)$y
        
        #get initial, final concentration and incubation time (whole duration)
        t.win <- 30
        C0 = min(mydf$concsmooth[mydf$time<t.win])
        Cf = max(mydf$concsmooth[mydf$time>max(mydf$time)-t.win])
        incubation_time = last(mydf$time)
        #get deltaconcs
        deltaconcs = Cf-C0
        #Calculate SD for initial final and deltaconc
        SD_C0 <- sd(mydf$conc[mydf$time<t.win])
        SD_Cf <- sd(mydf$conc[mydf$time>max(mydf$time)-t.win])
        SD_deltaconcs <- sqrt(SD_C0^2+SD_Cf^2)
        
        #Calculate full.total.flux aplying flux.term from aquaGHG
        full.total.flux <- (deltaconcs)/incubation_time*CH4all_fluxsep_flux.auto_i$flux.term # nmol/m2/s
        full.ebullition.flux <- full.total.flux-CH4all_fluxsep_flux.auto_i$diffusion.flux
        #Calculate SD of total flux
        SD_full.total.flux <- abs(full.total.flux) * SD_deltaconcs/deltaconcs
        
        #Re-calculate SD_full.ebullition.flux using updated SD_full.total.flux and extracted SD_diffusion.flux
        SD_diffusion.flux <- CH4all_fluxsep_flux.auto_i$diffusion.flux.SD
        SD_full.ebullition.flux <- sqrt(SD_diffusion.flux^2+SD_full.total.flux^2)
        
        #ADD to CH4all_fluxsep_flux.auto_i results for full.total flux and full.ebullition.flux
        CH4all_fluxsep_flux.auto_i$full.total.flux<- full.total.flux
        CH4all_fluxsep_flux.auto_i$SD_full.total.flux<- SD_full.total.flux
        CH4all_fluxsep_flux.auto_i$full.ebullition.flux<- full.ebullition.flux
        CH4all_fluxsep_flux.auto_i$SD_full.ebullition.flux<- SD_full.ebullition.flux
        
        }
        
        
        #Join flux to rest of CH4w_sep dataset
        CH4all_fluxsep_flux.auto <- rbind(CH4all_fluxsep_flux.auto, CH4all_fluxsep_flux.auto_i)
      }
    }
    
    #Remove auxfile for UniqueID (avoids re-usage, if next does not exist)
    rm(ch4_auxfile_i)
  }#End of UniqueID loop
  dev.off()
}#end subsite loop


# Save results to Rdata
save(list = c("CH4all_fluxsep_flux.auto",
              "CH4all_fluxsep_df.no_measurements"), file = paste0(results_path,"CH4all_fluxsep_aquaGHG.Rdata"))



#------_____________________----

#---OLD LOOPs (DONT RUN)-----


#----DONT RUN CH4 dry loop (no plots)-------
#This takes approx 6minutes
rm(k,i)
#loop over uniqueID of full auxfile
for (k in seq_along(ch4_dry_auxfile$UniqueID)){
  
  #Select uniqueID and display progress
  i = ch4_dry_auxfile$UniqueID[k]
  message(paste0("Processing CH4dry incubation ",k, " of ",length(ch4_dry_auxfile$UniqueID)," (",round(100*k/length(ch4_dry_auxfile$UniqueID),0), " %)"))
  
  #Load auxfiles for UniqueID
  ch4_dry_auxfile_i <- ch4_dry_auxfile[which(ch4_dry_auxfile$UniqueID==i),]
  
  #check that UniqueID has auxfile for ch4 dry
  if(dim(ch4_dry_auxfile_i)[1]==0){
    message(paste0("Could not find corresponding auxfile for ",i))
    CH4dry_df.no_measurements <- rbind(CH4dry_df.no_measurements,
                                       data.frame(UniqueID=i,
                                                  message = "no auxfile ch4"))
  } else {
    #Load ch4 dry data for UniqueID
    mydata_ch4_dry <- load_incubation(ch4_dry_auxfile_i, RData_path)
    
    #Check that there is data for UniqueID
    if(dim(mydata_ch4_dry)[1]==0){
      CH4dry_df.no_measurements <- rbind(CH4dry_df.no_measurements,
                                         data.frame(UniqueID=i,
                                                    message = "no measurements for ch4"))
    } else {
      #Calculate CH4 dry flux for UniqueID
      CH4dry_flux.auto_i <- automaticflux(dataframe = mydata_ch4_dry, myauxfile = ch4_dry_auxfile_i, shoulder = 0, gastype = "CH4dry_ppb", 
                                          fluxSeparation = FALSE,
                                          displayPlots = FALSE,
                                          method = "trust.it.all")
      
      #Join flux to rest of CH4 dry dataset
      CH4dry_flux.auto <- rbind(CH4dry_flux.auto, CH4dry_flux.auto_i)
    }
  }
  #Remove auxfile for UniqueID (avoids re-usage, if next does not exist)
  rm(ch4_dry_auxfile_i)
  
}#End of UniqueID loop


# Save results to Rdata
save(list = c("CH4dry_flux.auto",
              "CH4dry_df.no_measurements"), file = paste0(results_path,"CH4dry_results_aquaGHG.Rdata"))



#----DONT RUN CO2 loop (No plots)-----
#this takes approx 10 min
rm(k,i)
#loop over uniqueID of full auxfile
for (k in seq_along(co2_auxfile$UniqueID)){
  
  #Select uniqueID and display progress
  i = co2_auxfile$UniqueID[k]
  message(paste0("Processing CO2 incubation ",k, " of ",length(co2_auxfile$UniqueID)," (",round(100*k/length(co2_auxfile$UniqueID),0), " %)"))
  
  #Load auxfiles for UniqueID
  co2_auxfile_i <- co2_auxfile[which(co2_auxfile$UniqueID==i),]
  
  #check that UniqueID has auxfile for co2
  if(dim(co2_auxfile_i)[1]==0){
    message(paste0("Could not find corresponding auxfile for ",i))
    CO2_df.no_measurements <- rbind(CO2_df.no_measurements,
                                    data.frame(UniqueID=i,
                                               message = "no auxfile co2"))
  } else {
    #Load co2 data for UniqueID
    mydata_co2 <- load_incubation(co2_auxfile_i, RData_path)
    
    #Check that there is co2 data for UniqueID
    if(dim(mydata_co2)[1]==0){
      CO2_df.no_measurements <- rbind(CO2_df.no_measurements,
                                      data.frame(UniqueID=i,
                                                 message = "no measurements for co2"))
    } else {
      #Calculate CO2 flux for UniqueID
      CO2_flux.auto_i <- automaticflux(dataframe = mydata_co2, myauxfile = co2_auxfile_i, shoulder = 0, gastype = "CO2dry_ppm", 
                                       fluxSeparation = FALSE,
                                       displayPlots = FALSE,
                                       method = "trust.it.all")
      
      #Join flux to rest of CO2 dataset
      CO2_flux.auto <- rbind(CO2_flux.auto, CO2_flux.auto_i)
    }
  }
  #Remove auxfile for UniqueID (avoids re-usage, if next does not exist)
  rm(co2_auxfile_i)
  
}#End of UniqueID loop





#----- DONT RUN CH4 loop (No plots)------
rm(k,i)
for (k in seq_along(list_ids)){
  i = list_ids[k]
  message(paste0("processing ",i))
  
  #Load auxfiles for UniqueID
  ch4_auxfile_i <- ch4_auxfile[which(ch4_auxfile$UniqueID==i),]
  
  #check that UniqueID has auxfile for ch4
  if(dim(ch4_auxfile_i)[1]==0){
    message(paste0("Could not find corresponding auxfile for ",i))
    df.no_measurements <- rbind(df.no_measurements,
                                data.frame(UniqueID=i,
                                           message = "no auxfile ch4"))
  } else {
    #Load ch4 data for UniqueID
    mydata_ch4 <- load_incubation(ch4_auxfile_i, RData_path)
    
    #Check that there is ch4 data for UniqueID
    if(dim(mydata_ch4)[1]==0){
      df.no_measurements <- rbind(df.no_measurements,
                                  data.frame(UniqueID=i,
                                             message = "no measurements for ch4"))
    } else {
      #Calculate CH4 flux for UniqueID      
      CH4_flux.auto_i <- automaticflux(dataframe = mydata_ch4, myauxfile = ch4_auxfile_i, shoulder = 0, gastype = "CH4dry_ppb", 
                                       fluxSeparation = T,
                                       displayPlots = F, 
                                       method = "trust.it.all")
      
      #Join flux to rest of CH4 dataset
      CH4_flux.auto <- rbind(CH4_flux.auto, CH4_flux.auto_i)
    }
  }
  
  #Remove auxfile for UniqueID (avoids re-usage, if next does not exist)
  rm(ch4_auxfile_i)
}#End of UniqueID loop





#----DONT RUN CO2 loop (with pdf plots)-----
rm(subsite, k)
for (subsite in unique(co2_auxfile$subsite)){
  message(paste0("processing subsite",subsite))
  
  #Start pdf for plots
  pdf(file = paste0(results_path,"co2_",subsite,".pdf"))
  
  #loop over uniqueID of current subsite
  for (k in seq_along(unique(co2_auxfile[which(co2_auxfile$subsite==subsite),]$UniqueID))){
    #Select uniqueID and display it
    i = co2_auxfile[which(co2_auxfile$subsite==subsite),]$UniqueID[k]
    message(paste0("processing ",i))
    
    #Load auxfiles for UniqueID
    co2_auxfile_i <- co2_auxfile[which(co2_auxfile$UniqueID==i),]
    
    #check that UniqueID has auxfile for co2
    if(dim(co2_auxfile_i)[1]==0){
      message(paste0("Could not find corresponding auxfile for ",i))
      df.no_measurements <- rbind(df.no_measurements,
                                  data.frame(UniqueID=i,
                                             message = "no auxfile co2"))
    } else {
      #Load co2 data for UniqueID
      mydata_co2 <- load_incubation(co2_auxfile_i, RData_path)
      
      #Check that there is co2 data for UniqueID
      if(dim(mydata_co2)[1]==0){
        df.no_measurements <- rbind(df.no_measurements,
                                    data.frame(UniqueID=i,
                                               message = "no measurements for co2"))
      } else {
        #Calculate CO2 flux for UniqueID
        CO2_flux.auto_i <- automaticflux(dataframe = mydata_co2, myauxfile = co2_auxfile_i, shoulder = 0, gastype = "CO2dry_ppm", 
                                         fluxSeparation = FALSE,
                                         displayPlots = TRUE,
                                         method = "trust.it.all")
        
        #Join flux to rest of CO2 dataset
        CO2_flux.auto <- rbind(CO2_flux.auto, CO2_flux.auto_i)
        
        #Plot incubation:
        # print(plot_incubations(dataframe = mydata))
        
      }
      
    }
    
    #Remove auxfile for UniqueID (avoids re-usage, if next does not exist)
    rm(co2_auxfile_i)
    
  }#End of UniqueID loop
  
  dev.off()
}#end subsite loop



