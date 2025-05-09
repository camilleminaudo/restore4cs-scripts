# ---
# Authors: Camille Minaudo, MIGUEL edits 
# Project: "RESTORE4Cs"
# date: "Avril 2025"

# ---

# --- Description----
#this script uses the functions from aquaGHG to calculate CO2 and CH4 fluxes for all incubations separately. Takes auxfiles created with script Create_CO2andCH4_auxfiles_4aquaGHG.R and Rdata updated with script Harmonize_GHG_raw2Rdata.R (accurate timezones and correct units, last update 6/5/2025)

#separate loops for CO2 and CH4 (each can have different start.time and duration for a given UniqueID)

#Needs for restore4cs: 
#CO2: no need to calulate fluxes with aquaGHG, all results have been correctly produced and model chosen with previous raw2flux scripts. Done here for consistency and to test. 

#CH4: we need to use aquaGHG to be able to separate difusion and ebullition and decide on method to produce final estimates ("best.flux"). For it we need to plot to evaluate the performance of the approach in multiple UniqueID cases. 


#BUGs and issues------
# #SOME ISSUES WITH start.time class (rounding causes class change to POSIXt), easy fix by re-transforming to POSIXct



#BUG: flux.plot fails when HM.flux is NA (breaking the loop): 
#Error en flux.plot(flux.results = best.flux_auto, dataframe = mydata_auto, : 
# 'HM.flux' in 'flux.results' must be of class numeric
# Además: Aviso:
  # Flux estimate is too close to zero to estimate HM flux in UniqueID s3-ca-a2-3-v-d-08:05. NAs produced. 

#In addition, UniqueID not displayed in plots for co2. (UniqueID is displayed without issue for ch4 with fluxSeparation = T and method = "trust.it.all")


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

#Set testing results_path
results_path<- "C:/Users/Miguel/Dropbox/testing_aquaGHG/"

# results_path <- paste0(dropbox_root,"/GHG/Processed data/computed_flux/")


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


#-----Separate CH4 auxfiles dry/water-----
#Separate ch4_auxfile based on water depth. Incubations with water-depth==0 will be calculated with goflux directly. Incubations with water will be checked for ebullition.

ch4_dry_auxfile<- ch4_auxfile %>% filter(water_depth==0)

ch4_water_auxfile<- ch4_auxfile %>% filter(water_depth>0)



#-----Initialize objects----
CO2_flux.auto <- CH4dry_flux.auto <- CH4water_flux.auto <- NULL

CO2_df.no_measurements <- CH4dry_df.no_measurements <- CH4water_df.no_measurements <- NULL

#----CO2 loop (No plots)-----
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
  

# Save results to Rdata
save(list = c("CO2_flux.auto",
              "CO2_df.no_measurements"), file = paste0(results_path,"CO2_results_aquaGHG_recalculated.Rdata"))


#----CH4 dry loop (no plots)-------
rm(k,i)
#loop over uniqueID of full auxfile
for (k in seq_along(ch4_dry_auxfile$UniqueID)){
  
  #Select uniqueID and display progress
  i = ch4_dry_auxfile$UniqueID[k]
  message(paste0("Processing CH4dry incubation ",k, " of ",length(ch4_dry_auxfile$UniqueID)," (",round(100*k/length(ch4_dry_auxfile$UniqueID),0), " %)"))
  
  #Load auxfiles for UniqueID
  ch4_dry_auxfile_i <- ch4_dry_auxfile[which(ch4_dry_auxfile$UniqueID==i),]
  
  #check that UniqueID has auxfile for ch4 dry
  if(dim(ch4_dry_auxfile)[1]==0){
    message(paste0("Could not find corresponding auxfile for ",i))
    CH4dry_df.no_measurements <- rbind(CH4dry_df.no_measurements,
                                data.frame(UniqueID=i,
                                           message = "no auxfile ch4"))
  } else {
    #Load ch4 dry data for UniqueID
    mydata_ch4_dry <- load_incubation(ch4_dry_auxfile_i, RData_path)
    
    #Check that there is co2 data for UniqueID
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
              "CH4dry_df.no_measurements"), file = paste0(results_path,"CH4dry_results_aquaGHG_recalculated.Rdata"))



#----CH4 water loop (with pdf plots)-----
rm(subsite, k)
for (subsite in unique(ch4_water_auxfile$subsite)){
  message(paste0("processing subsite",subsite))
  
  #Start pdf for plots
  pdf(file = paste0(results_path,"CH4water_",subsite,".pdf"))
  
  #loop over uniqueID of current subsite
  for (k in seq_along(unique(ch4_water_auxfile[which(ch4_water_auxfile$subsite==subsite),]$UniqueID))){
    
    #Select uniqueID and display progress
    i = ch4_water_auxfile[which(ch4_water_auxfile$subsite==subsite),]$UniqueID[k]

    message(paste0("Processing CH4water incubation ",k, " of ",length(ch4_water_auxfile$UniqueID)," (",round(100*k/length(ch4_water_auxfile$UniqueID),0), " %)"))
    
    #Load auxfiles for UniqueID
    ch4_water_auxfile_i <- ch4_water_auxfile[which(ch4_water_auxfile$UniqueID==i),]
    
    #check that UniqueID has auxfile for ch4 (with water)
    if(dim(ch4_water_auxfile_i)[1]==0){
      message(paste0("Could not find corresponding auxfile for ",i))
      CH4water_df.no_measurements <- rbind(CH4water_df.no_measurements,
                                  data.frame(UniqueID=i,
                                             message = "no auxfile ch4"))
    } else {
      #Load ch4 data for UniqueID with water
      mydata_ch4_water <- load_incubation(ch4_water_auxfile_i, RData_path)
      
      #Check that there is ch4 data for UniqueID
      if(dim(mydata_ch4_water)[1]==0){
        CH4water_df.no_measurements <- rbind(CH4water_df.no_measurements,
                                    data.frame(UniqueID=i,
                                               message = "no measurements for ch4"))
      } else {
        #Calculate CH4 flux for UniqueID  with flux separation and plot    
        CH4water_flux.auto_i <- automaticflux(dataframe = mydata_ch4_water, myauxfile = ch4_water_auxfile_i, shoulder = 0, gastype = "CH4dry_ppb", 
                                         fluxSeparation = T,
                                         displayPlots = TRUE, 
                                         method = "trust.it.all")
        
        #Join flux to rest of CH4water dataset
        CH4water_flux.auto <- rbind(CH4water_flux.auto, CH4water_flux.auto_i)
      }
    }
    
    #Remove auxfile for UniqueID (avoids re-usage, if next does not exist)
    rm(ch4_water_auxfile_i)
  }#End of UniqueID loop
  dev.off()
}#end subsite loop


# Save results to Rdata
save(list = c("CH4water_flux.auto",
              "CH4water_df.no_measurements"), file = paste0(results_path,"CH4water_results_aquaGHG_recalculated.Rdata"))



#Ch4 water auxfile 4 review------
#save ch4 auxfile to log mistakes that could be corrected (i.e. cases of obious ebullition from visual inspection) log in them sucess or miss (with short description)

write.csv(x = ch4_water_auxfile, file = paste0(results_path, "ch4_water_auxfile4inspection.csv"), row.names = F)



#------_____________________----


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







#Miscelanea-----


#UniqueID s1-ca-r1-1-o-d had an artefact at the beginning first ~260 s with constant ch4, check that the actual data. 
#There was an issue with an incomplete rawfile in the wrong folder that caused data to be incomplete for this UniqueID, wrong rawfile was deleted and rdata re-imported. Now this is fixed. 

#WE will have to re-run alos the co2 loop to update the flux of this UniqueID. 
ch4_auxfile %>% filter(grepl("s1-ca-r1-1-o-d", UniqueID))
co2_auxfile %>% filter(grepl("s1-ca-r1-1-o-d", UniqueID))
#Duration is 710 in both auxfiles, starting at 8:40:39

a<-load_incubation(auxfile_i = ch4_auxfile %>% filter(grepl("s1-ca-r1-1-o-d", UniqueID)), RData_path = RData_path)
#when loaded, only 450 s of data. starting at 8:45:00
a$POSIX.time[1]
a %>% 
  ggplot(aes(x=Etime, y=CO2dry_ppm))+
  geom_point()

a %>% 
  ggplot(aes(x=POSIX.time, y=CH4dry_ppb))+
  geom_point()

#NO artefacts in loaded data

#Check Rdata file: the Rdatafile starts at 8:45:00 
load(file = paste0(RData_path, "/S1-CA-R1_LI-7810.RData"))
first(mydata$POSIX.time) 

start<- ch4_auxfile %>% filter(grepl("s1-ca-r1-1-o-d", UniqueID)) %>% pull(start.time)
stop<- start+ ch4_auxfile %>% filter(grepl("s1-ca-r1-1-o-d", UniqueID)) %>% pull(duration)

mydataauxfile<- mydata %>% 
  filter(between(POSIX.time,start,stop ))

rawdata<- read.delim(file=)
#Check raw-data: raw data file has non-NA ghg data starting from 08:23:29 (more than 20 minutes before start of incubation). 
#first non-NA POSIX.time in raw file: 2023-10-31	08:23:29, 

#raw-data file has remark starting at 08:40:34 (start of incubation)

#It seems the first part of rawdata file is not imported as Rdata  #WHY???!!

#Check data in rawdata and Rdata is the same (match ghg concentrations and time?)

#YES!! Same data in both,
# 08:45:00	13898.349	  439.44772	  2034.1851 in rawfile
# 08:45:00  13898.35  439.4477  2034.185 in Rdata 
 
#fieldsheet start.time: 8:41 

#Only "filter" when harmonizing Rdata is removing duplicate Posixtimes
#In data.raw (Rdata direct from raw-file) only 4 duplicated times, and not within incubation timeframe. 

#different number of lines in data.raw (Rdata direct from rawfile)

is_duplicate <- duplicated(data.raw$POSIX.time)
data.raw_nodup <- data.raw[!is_duplicate,]

#FOUND a different raw-file with same data,but starting at 8:45:00, in wrong folder that seems to override the data (deleted incomplete file of S1-ca-r1 day from folder S1-CU-DA-DU). 
