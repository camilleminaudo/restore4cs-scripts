
# ---
# Authors: Miguel Cabrera
# Project: "RESTORE4Cs"
# date: "Dec 2023"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
# This script re-formats the fluxes calculated in raw2flux script for the preeliminary data structure to send to LifeWatch Italy

#Input: per-season CSV files with computed fluxes

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
library(tools)

repo_root <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
files.sources = list.files(path = paste0(repo_root,"/functions"), full.names = T)
for (f in files.sources){source(f)}


# ---- Directories ----
dropbox_root <- "C:/Users/Miguel/Dropbox/RESTORE4Cs - Fieldwork/Data" # You have to make sure this is pointing to the write folder on your local machine
# datapath <- paste0(dropbox_root,"/GHG/RAW data")
fieldsheetpath <- paste0(dropbox_root,"/GHG/Fieldsheets")
loggerspath <- paste0(dropbox_root,"/GHG/RAW data/RAW Data Logger")
# RData_path <- paste0(dropbox_root,"/GHG/Processed data/RData/")
results_path <- paste0(dropbox_root,"/GHG/Processed data/computed_flux/")
vegetation_path<- paste0(dropbox_root,"/Vegetation/")
lightplots_path<- paste0(dropbox_root,"/GHG/Processed data/computed_flux/plots/")
# plots_path <- paste0(dropbox_root,"/GHG/Processed data/plots_all_incubations/")



# Final variables and units: 
#season: summer,spring,autum, winter
#case_pilot: CA,DA,CU,RI,VA,DU
#subsite: casepilot_A1/A2/P1/P2/R1/R2
#date:
#plot_latitude
#plot_longitude
#plot_waterdepth: in cm
#plot_strata: open water, bare, vegetated land, vegetated water
#plot_vegetation_biomass
#plot vegetation description
#plot_uniqueID: subsite_plotnumber_o/b/vl/vw
#plot_starttime: first start.datetime of plot
#transparent_light_intensity: lux of transparent incubation
#dark_light_intensity: lux of dark incubation
#transparent_temperature: ºC of transparent
#dark_temperature: ºC of dark
#transparent_co2flux: umol/m2/s
#dark_co2flux: umol/m2/s
#transparent_diff_ch4flux: nmol/m2/s
#dark_diff_ch4flux: nmol/m2/s
#transparent_ebull_ch4flux: nmol/m2/s
#dark_ebull_ch4flux: nmol/m2/s
#transparent_n2oflux: nmol/m2/s
#dark_n2oflux: nmol/m2/s


# ---- Import fluxes per season ----

#Read all computed fluxes per season
season_csvs<- list.files(path = results_path, pattern = "^S")
dat<- data.frame()
for (i in season_csvs){
  s<- read.csv(paste0(results_path,i))
  
  dat<- rbind(dat,s)
}
rm(s, i)


str(dat)


#Keep only fluxes data, start.time and UniqueID.
#For CO2 fluxes, take only CO2_best.flux
#For CH4 fluxes, take CH4_diffusive_flux and CH4_ebullitive_flux
#We will have to adapt the raw2flux script so that CH4 has always data for diffusion and ebullition (if no ebullition, ebullition==0)

dat2<- dat %>% 
  select(UniqueID,
         start.time,
         CO2_best.flux,
         #CO2_LM.flux,CO2_HM.flux,CO2_best.flux,CO2_best.model,CO2_quality.check,
         #CH4_LM.flux,CH4_HM.flux,CH4_best.flux,CH4_best.model, CH4_quality.check,
         CH4_diffusive_flux,CH4_ebullitive_flux)



#Import fieldsheet meta-data#----

#Read all GHGfieldsheets (using function from Camille)
fieldsheets<- list.files(path= fieldsheetpath, recursive = T, pattern = "GHG.xlsx",full.names = T)
field<- read_GHG_fieldsheets(fieldsheets)

#Drop unwanted variables and rename key UniqueID:
field2<- field %>% 
  select(-c(person_sampling,gas_analyzer,logger_floating_chamber,logger_transparent_chamber,start_time,initial_co2,initial_ch4,end_time,final_ch4,final_co2,chamber_height_cm, unix_start,unix_stop)) %>% 
  rename(UniqueID=uniqID)


# ---- Format and join ----
#Key for joining flux data and fieldsheet: 
summary(dat$UniqueID%in%field$uniqID)


#Check for duplicate incubations in fieldsheets (incubations performed 2 times with the same identity except for startime and uniqueID:  
duplicate_incubations<-field2 %>% 
  dplyr::summarise(n = dplyr::n(), .by = c(pilot_site, subsite, date, plot_id, chamber_type, strata, longitude, latitude, water_depth,
                                           transparent_dark)) |>dplyr::filter(n > 1L) %>% 
  merge.data.frame(field2, by=c("pilot_site", "subsite", "date", "plot_id", "chamber_type", "strata", "longitude", "latitude", "water_depth","transparent_dark"), all = F) %>% pull(UniqueID)


#Check for missing transparent_dark
field2 %>% filter(is.na(transparent_dark))


#merge fluxes with fieldsheets dropping the unmatched observations, problematic values and duplicated incubations
ghg_formated<- field2 %>% merge.data.frame(dat2, by="UniqueID",all = F) %>%
  filter(!is.na(transparent_dark)) %>% 
  filter(!UniqueID%in%duplicate_incubations) %>% 
  select(-UniqueID) %>% 
  pivot_wider(names_from = transparent_dark,values_from = c(start.time, CO2_best.flux, CH4_diffusive_flux,CH4_ebullitive_flux,comments)) %>% 
  mutate(plot_startime=min(start.time_dark, start.time_transparent, na.rm=T),
         plotcode=paste(subsite, plot_id,sep = "-"),
         campaign=str_extract(subsite, pattern="S[0-9]{1}"),
         season=case_when(campaign=="S1"~"fall",
                          campaign=="S2"~"winter",
                          campaign=="S3"~"spring",
                          campaign=="S4"~"summer"),
         subsite=str_extract(subsite, pattern="[A-Z]{2}-[A-Z]{1}[0-9]{1}"),
         sampling=paste(campaign, subsite,sep="-"),
         strata=case_when(strata=="vegetated"&water_depth>=10~"vegetated water",#Re-classify strata using water depth (vegetated water when water_depth>=10cm, otherwise vegetated land)
                          strata=="vegetated"&water_depth<10~"vegetated land",
                          TRUE~strata)) %>% 
  select(season, campaign,pilot_site, subsite, sampling, date, latitude, longitude, water_depth, strata, plotcode, plot_startime, CO2_best.flux_dark, CO2_best.flux_transparent, CH4_diffusive_flux_dark, CH4_ebullitive_flux_dark, CH4_diffusive_flux_transparent, CH4_ebullitive_flux_transparent)




# ---- Vegetation ID & DW ----

#Vegetation description and biomassDW: 
#1. Format vegetation files to Identify based on Season-site-subsite-plot information. Final variables: samplecode, DW (g), comments_veg.DONE


#List all vegetation files 
#Import pattern: RESTORE4Cs_vegetation.........xlsx (rbind all). 
veg_files<- list.files(path = vegetation_path, pattern = "^RESTORE4Cs_vegetation_",full.names = T,recursive = T)

#Read all computed fluxes per season
veg<- read_vegetation_fieldsheets(veg_files)

#Duplicated plant DWs
veg %>% 
  dplyr::summarise(n = dplyr::n(), .by = c(plotcode)) |>dplyr::filter(n > 1L) %>% 
  merge.data.frame(veg, by=c("plotcode"), all = F)


#2. Extract from fieldsheets the comments of vegetated plots. Format into comments_dark and comments_transparent
vegetated_plots<- field2 %>% 
  filter(!UniqueID%in%c("s1-va-a1-12-v-t-14:58","s1-va-a1-12-v-d-15:06")) %>% #remove repeated incubation
  filter(strata=="vegetated") %>%#Select only vegetated plots
  mutate(plotcode=paste(subsite, plot_id,sep="-")) %>% 
  select(plotcode, strata,transparent_dark, comments) %>% 
  pivot_wider(values_from = comments, names_from = transparent_dark) %>% 
  distinct()


#Join fieldsheet and vegetation info
vegetated_plots_veg<- vegetated_plots %>% 
  merge.data.frame(veg, by=c("plotcode"),all = T)

#vegetation DW data without vegetated plotGHG: 7 unmatched vegetation samples
veg[!veg$plotcode%in%vegetated_plots$plotcode,]


#S2-RI, S1-CA, S2-CA, S2-CU, S3-VA, S4-CA: unmatched vegetation data was corrected based on sample distribution of chambers and vegetation samples. Modifications are commented in vegetation excel.  
#S2-RI had some correction, but remaining unmatched vegetation data cannot be assigned to vegetated chambers (no obvious correspondence or typos)


#Vegetated plotsGHG without vegetationDW 92 chambers without vegetation data
vegetated_plots[!vegetated_plots$plotcode%in%veg$plotcode,]



#3. Join and combine comments to get vegetation description (sp/genus/family/), calculate ABG_biomass using tube area (pi*(12.1cm radius)^2)
#Unmatched vegetation data (veg data without chamber) are excluded
vegetation_dw_final<- vegetated_plots %>% 
  merge.data.frame(veg, by=c("plotcode"),all.x = T) %>% 
  mutate(ABG_biomass_gpersquaremeter=dry_weight/(pi*0.121^2)) %>%# g/square m (radius in meters) 
  select(plotcode,
         strata, ABG_biomass_gpersquaremeter, transparent, dark, comments) %>% 
  rename(comment1=transparent, comment2=dark, comment3=comments)

write.csv(vegetation_dw_final,file = paste0(vegetation_path, "RESTORE4Cs_ABGbiomass_vegetated_plots.csv"), row.names = F)


rm(veg, vegetated_plots,vegetated_plots_veg,vegetation_dw_final)
#DONE____Manually, harmonize plot vegetation description: DONE, some inconsistencies found.

#Check harmonized vegetation

veg_formated<- read_xlsx(path = paste0(vegetation_path,"RESTORE4Cs_ABGbiomass_and_description_vegetated_plots.xlsx"),na = "NA")

#check unique vegetation descriptors:
unique(veg_formated$vegetation_description)
mean(veg_formated$ABG_biomass_gpersquaremeter,na.rm = T)


#Different levels of specificity: Genus, Family, Common description (eg. Reed, Grass, Macroalgae, Charophyte)
#Many NAs that could be filled by looking at sampling pictures




# ---- Data loggers ----
#We need Temp and Light intensity for dark and transparent conditions. 
#extract from data loggers files using info from field: logger_floating_chamber, logger_transparent_chamber, unixstart, unixstop, uniqID, chamber_type
#In raw2flux, temperature extraction already implemented(CHECK and adapt for Light intensity)


loggers_maps<- field %>% 
  select(subsite, logger_floating_chamber,logger_transparent_chamber,chamber_type, uniqID, unix_start, unix_stop) %>% 
  filter(!uniqID%in%duplicate_incubations)







#Most of the below  code is borrowed from raw2flux script.

subsites <- unique(loggers_maps$subsite)
logger_summary<- data.frame()

#Initialize list of plots to save lightplots
plotslight <- list()

for (subsite in subsites){
  message("Now processing ",subsite)
  
  corresp_fs <- loggers_maps[loggers_maps$subsite == subsite,]
  
  # read corresponding temperature logger file and keep initial temperature
  SN_logger_float <- first(corresp_fs$logger_floating_chamber)
  SN_logger_tube <- first(corresp_fs$logger_transparent_chamber)
  site_ID <- str_sub(subsite, start = 1, end = 5)
  
  # finding out if corresponding file exists, and its extension
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
  
  #Read data from data logger float (sn, datetime, temperature, unixtime)
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
    #Temperature as numeric:
    data_logger_float$temperature<- as.numeric(data_logger_float$temperature)
    if(is.character(data_logger_float$datetime)){
      data_logger_float$datetime <- as.POSIXct(data_logger_float$datetime, tz = 'utc', tryFormats = c("%m/%d/%y %r", "%d/%m/%Y %H:%M"))
    }
    data_logger_float$unixtime <- as.numeric(data_logger_float$datetime)
  } else {
    message("===> no data logger linked to the floating chamber!")
    is_data_logger_float = F}
  
  #Check that we actually have data from logger
  if(dim(data_logger_float)[1]<10){
    message("===> not enough data could be linked to the floating chamber!")
    is_data_logger_float = F
  }
  
  #Read data from data logger tube (sn, datetime, temperature, unixtime)
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
    #Temperature and light as numeric
    data_logger_tube$temperature <- as.numeric(data_logger_tube$temperature)
    data_logger_tube$light <- as.numeric(data_logger_tube$light)
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
  #Check that we actually have data from logger
  if(dim(data_logger_tube)[1]<10){
    message("===> not enough data could be linked to the tube chamber!")
    is_data_logger_tube = F
  }
  
  
  
  for (incub in seq_along(corresp_fs$uniqID)){
    
    #if floating -> subset temp data and calculate descriptive stats for temperature. Add results to loggersdata (subsite, uniqID, descriptive stats)
    
    #If tube-> subset logger data, calculate descriptive stats for temp and light (mean, sd, median). Add results to loggersdata (subsite, uniqID, descriptive stats)
    
    #Get start and stop of incubation, 
    incubstart<- corresp_fs[incub,]$unix_start
    incubstop<- corresp_fs[incub,]$unix_stop
    
  #IF floating: 
    if (corresp_fs$chamber_type[incub] == "floating"){
      if(is_data_logger_float){
        
      incub_log<- data_logger_float[data_logger_float$unixtime>=incubstart&data_logger_float$unixtime<incubstop,]  
      
      #Calculate descriptive stats for temp (fill light variables with NA and 0 obs)
      incubsave<- data.frame(
        "uniqID"=corresp_fs$uniqID[incub],
        "temp_mean"=mean(incub_log$temperature, na.rm=T),
        "temp_sd"=sd(incub_log$temperature, na.rm = T),
        "temp_median"=median(incub_log$temperature,na.rm = T),
        "temp_n"=sum(!is.na(incub_log$temperature)),
        "light_mean"=as.numeric(NA),
        "light_sd"=as.numeric(NA),
        "light_median"=as.numeric(NA),
        "light_n"=sum(!is.na(NA))
      )
      
      }
      #IF tube 
    } else if (corresp_fs$chamber_type[incub] == "tube"){
      if(is_data_logger_tube){
        
        #Subset data from incubation period
        incub_log<- data_logger_tube[data_logger_tube$unixtime>=incubstart&data_logger_tube$unixtime<incubstop,]  
        
        #Calculate descriptive stats for temp and light
        incubsave<- data.frame(
          "uniqID"=corresp_fs$uniqID[incub],
          "temp_mean"=mean(incub_log$temperature, na.rm=T),
          "temp_sd"=sd(incub_log$temperature, na.rm = T),
          "temp_median"=median(incub_log$temperature,na.rm = T),
          "temp_n"=sum(!is.na(incub_log$temperature)),
          "light_mean"=mean(incub_log$light, na.rm=T),
          "light_sd"=sd(incub_log$light, na.rm = T),
          "light_median"=median(incub_log$light,na.rm = T),
          "light_n"=sum(!is.na(incub_log$light))
        )
        
        }
    } else {
      warning("chamber type not correct!")
    }
    logger_summary<- rbind(logger_summary,incubsave)
  }
  
  if(is_data_logger_tube){
  #For each subsite, make a plot with the light sensor data and the corresponding transparent_dark incubations
    
  #add UniqID to data logger tube
    corresp_fs_tube<- corresp_fs %>% filter(chamber_type=="tube")
  data_logger_tube$uniqID <- sapply(data_logger_tube$unixtime, function(unixtime) {
    # Find the label for each unixtime based on the unixstart and unixstop
    uniqID <- corresp_fs_tube$uniqID[corresp_fs_tube$unix_start <= unixtime & corresp_fs_tube$unix_stop >= unixtime]
    if (length(uniqID) > 0) {
      return(uniqID)
    } else {
      return(NA)  # or another value like "No label" if no match
    }
  })
  #add light condition to data logger ("t" or "d")
  data_logger_tube$condition_light<- substr(str_extract(data_logger_tube$uniqID, pattern = "[a-z]{1}-[0-9]{2}"),1,1)
  
  #PLot of the uniqIDs in the light data (with 3 h before and after the incubation times)
  p<-data_logger_tube %>% 
    filter(unixtime>min(corresp_fs_tube$unix_start, na.rm=T)-(60*180)&unixtime<=max(corresp_fs_tube$unix_stop, na.rm=T)+(60*180)) %>% 
    filter(!is.na(light)) %>% 
    ggplot(aes(x=datetime, y=light))+
    #ADDING labels fails (error when printing pdf within loop)
    # geom_label(data = data_logger_tube %>%  filter(!is.na(uniqID)) %>% group_by(uniqID) %>% summarise(datetime=mean(datetime),light=first(light)), aes(x=datetime, label = uniqID, y=0),vjust=0.5,hjust = 1,angle = 90)+
    geom_line()+
    geom_point(data=subset(data_logger_tube, !is.na(condition_light)),aes(col=condition_light))+
    labs(col="Light?")+
    scale_x_datetime(breaks = "hour")+
    # coord_cartesian(ylim = c(-40000, NA)) +  
    #scale y adjustment is what fails!
    # scale_y_continuous(limits = c(-45000,
    #                               max(data_logger_tube %>%
    #                                     filter(unixtime>min(corresp_fs$unix_start, na.rm=T)-(60*180)&unixtime<=max(corresp_fs$unix_stop,na.rm=T)+(60*180)) %>% pull(light))))+
    ggtitle(paste(subsite, "SN:",SN_logger_tube))+
    theme_bw()
  # }else{
  #   p<- ggplot()+ ggtitle(paste(subsite, "SN:",SN_logger_tube))+theme_bw()
  }
  
  # Store each subsite light data plot in the list
  plotslight[[subsite]] <- p
  
}

#save lightplot for every subsite: 
pdf(file = paste0(lightplots_path,"Lightsensor_per_subsite.pdf"),width = 15)  # Open PDF device

# Loop through the list of plots and print each plot
for (plot_name in names(plotslight)) {
  print(plotslight[[plot_name]])
}

dev.off()  # Close the PDF device


####POR AQUI####

#Up to here works, but many problems reading the tube_loggers. 
#Some subsites have duplicated data logger info (eg. light sensor for S1-CU-P2 and S1-CU-P2 have the exact same light data). Check files and check loop (to know that the plot is not re-using the last valid  data_logger_tube data)
#Cannot find a way of consistently plotting the uniqID labels inside each plot (printing fails when there is an issue with the data)

#ALSO probably time from fieldsheets not appropiate (alternative?).

#REVISE import of data loggers for timezone changes. 
#ASK Camille if he adjusted the times of LICOR PICARRO, LOS GATOS. 

names(plotslight)[1]




