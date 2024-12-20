
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
lifewatch_example_path<- paste0(dropbox_root,"/GHG/Processed data/computed_flux/LifeWatch Italy GHG dataset example/")
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
  rowwise() %>% 
  mutate(plot_startime=min(start.time_dark, start.time_transparent, na.rm=T)) %>% 
  ungroup() %>% 
  mutate(plotcode=paste0(subsite,"-",plot_id,"-",toupper(substr(strata,1,1))),#Add first capitalized letter of strata to plotcode
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


write.csv(ghg_formated, file = paste0(lifewatch_example_path, "GHG_per_plot_formated.csv"),row.names = F)


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


veg_formated<- veg_formated %>% 
  select(plotcode,ABG_biomass_gpersquaremeter,vegetation_description) %>% 
  mutate(plotcode=paste0(plotcode,"-V"))#add 1st letter of strata to plotcode

write.csv(veg_formated, file = paste0(lifewatch_example_path, "vegetation_per_plot_formated.csv"),row.names = F)



# ---- Data loggers ----
#We need Temp and Light intensity for dark and transparent conditions. 
#extract from data loggers files using info from field: logger_floating_chamber, logger_transparent_chamber, unixstart, unixstop, uniqID, chamber_type
#In raw2flux, temperature extraction already implemented(CHECK and adapt for Light intensity)

#Get info from fieldsheets
loggers_maps<- field %>% 
  select(subsite, logger_floating_chamber,logger_transparent_chamber,chamber_type, uniqID, unix_start, unix_stop) %>%
  #REMOVE DUPLICATE INCUBATIONS
  filter(!uniqID%in%duplicate_incubations)







#Most of the below  code is borrowed from raw2flux script.

subsites <- unique(loggers_maps$subsite)

# Use expand.grid to create all combinations
designed_subsites <- expand.grid(campaign = c("S1", "S2", "S3", "S4"), pilot_site = c("CA", "CU", "DA", "DU", "VA", "RI"), status = c("A1", "A2", "P1", "P2", "R1", "R2")) %>% mutate(all_subsites=paste(campaign,pilot_site,status,sep = "-")) %>% pull(all_subsites)

#Check no subsite missing or duplicate
unique(subsites %in%designed_subsites)
unique(designed_subsites %in%subsites)
rm(designed_subsites)

#Check NAs in logger data (subsite, SN loggers, chambertype, uniqID, unixstart, unixstop)
loggers_maps_nas <- loggers_maps[apply(loggers_maps, 1, function(x) any(is.na(x))), ]

print(loggers_maps_nas)
#14 incubations with undefined time-period (check fieldsheets and correct)

#Extract data from loggers only for well-defined incubations
complete_loggers_maps<- loggers_maps %>% drop_na()


#Initialize data frame to save logger data
logger_summary<- data.frame()

#Initialize list of plots to save lightplots
plotslight <- list()

for (subsite in subsites){
  message("Now processing ",subsite)
  
  corresp_fs <- complete_loggers_maps[complete_loggers_maps$subsite == subsite,]
  
  #get details from fieldsheet to find data_logger files:
  SN_logger_float <- first(corresp_fs$logger_floating_chamber)
  SN_logger_tube <- first(corresp_fs$logger_transparent_chamber)
  site_ID <- str_sub(subsite, start = 1, end = 5)
  
  # finding out if corresponding data_logger files exist, and their extension
  dir_exists_loggdata <- dir.exists(paste0(loggerspath,"/",site_ID,"/"))
  if(dir_exists_loggdata){
    f <- list.files(paste0(loggerspath,"/",site_ID,"/"), full.names = T)
    #Discard hobo files and txt files as sources for data_logger (only import from .csv or .xlsx)
    r <- grep(pattern = ".hobo", x = f)
    if(length(r)>0){f <- f[-r]}
    r <- grep(pattern = ".txt", x = f)
    if(length(r)>0){f <- f[-r]}
    f_ext <- file_ext(f)
    i_f_float <- grep(pattern = SN_logger_float, x = f)[1]
    i_f_tube <- grep(pattern = SN_logger_tube, x = f)[1]
  }
  
  
  #If the data logger float file is found, read it and get data (sn, datetime, temperature, unixtime)
  if(!is.na(SN_logger_float) & !is.na(i_f_float)){
    is_data_logger_float = T
    # message("...reading corresponding temperature logger file for floating chamber")
    if(f_ext[i_f_float]=="xlsx"){
      data_logger_float <- readxl::read_xlsx(f[i_f_float],col_names = T)
    } else if(f_ext[i_f_float]=="csv"){
      data_logger_float <- read.csv(f[i_f_float], header = T)
    }
    data_logger_float <- data_logger_float[,seq(1,3)]
    names(data_logger_float) <- c("sn","datetime","temperature")
    #add sn instead of rownumber
    data_logger_float$sn <- SN_logger_float
    
    #Format data (as numeric, datetime,..)
    data_logger_float$temperature<- as.numeric(data_logger_float$temperature)
    if(is.character(data_logger_float$datetime)){
      data_logger_float$datetime <- as.POSIXct(data_logger_float$datetime, tz = 'utc', tryFormats = c("%m/%d/%y %r", "%d/%m/%Y %H:%M"))
    }
    data_logger_float$unixtime <- as.numeric(data_logger_float$datetime)
  } else {
    #If we dont find the float file, message and record its absence within loop
    message("===> no data logger file linked to the floating chamber!")
    is_data_logger_float = F}
  
  #Check that we have enough float data from file
  if(dim(data_logger_float)[1]<10){
    message("===> not enough data could be linked to the floating chamber!")
    is_data_logger_float = F
  }
  
  
  #IF the data logger tube is found, read it and get data (sn, datetime, temperature, light, unixtime)
  if(!is.na(SN_logger_tube) & !is.na(i_f_tube)){
    is_data_logger_tube = T
    # message("...reading corresponding temperature logger file for tube chamber")
    
    if(f_ext[i_f_tube]=="xlsx"){
      data_logger_tube <- readxl::read_xlsx(f[i_f_tube],col_names = T)
    } else if(f_ext[i_f_tube]=="csv"){
      data_logger_tube <- read.csv(f[i_f_tube], header = T, fill = T)
    }
    ##ISSUE: error when reading------
    #Next line throws an error if data has less than 4 columns. If this happens, the loop re-uses the last data_logger_tube for all calculations and plots. 
    #Decission: make a check for 4 columns or skip
    data_logger_tube <- data_logger_tube[,seq(1,4)]
    names(data_logger_tube) <- c("sn","datetime","temperature","light")
    #add sn instead of rownumber
    data_logger_tube$sn <- SN_logger_tube
    #Format data (as.numeric, datetime, )
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
    #If we dont find the tube file, message and record its absence within loop
    is_data_logger_tube = F
    message("===> no data logger linked to the tube chamber!")
  }
  #Check that we have enough tube data from file
  if(dim(data_logger_tube)[1]<10){
    message("===> not enough data could be linked to the tube chamber!")
    is_data_logger_tube = F
  }

  
  

  #Loop over incubations of subsite and get 1 row with descriptive statistics for temperature and light. 
  for (incub in seq_along(corresp_fs$uniqID)){
    
    #Get start and stop of incubation, 
    incubstart<- corresp_fs[incub,]$unix_start
    incubstop<- corresp_fs[incub,]$unix_stop
    
  #IF incubation is floating: subset data float and calculate descriptive stats for temperature. Add row to loggers_summary
  #IF incubation is tube: subset data tube and calculate descriptive stats. Add row to loggers_summary  
    if (corresp_fs$chamber_type[incub] == "floating"){
      #IF there is data_logger float
      if(is_data_logger_float){
        
        #Subset sensor data from incubation period
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
        #Add incubsave floating to logger_summary
        logger_summary<- rbind(logger_summary,incubsave)
        }
  #IF incubation is tube:
    } else if (corresp_fs$chamber_type[incub] == "tube"){
      #If there is data_logger float
      if(is_data_logger_tube){
        
        #Subset sensor data from incubation period
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
        #Add incubsave floating to logger_summary
        logger_summary<- rbind(logger_summary,incubsave)
        
        }
    } else {warning("chamber type not correct!")}

  }
  
  
  #For each subsite, check if we have data in data_logger_tube, if so, make a plot with the light sensor data and the corresponding transparent_dark incubations. 
  #If we do not have data, create an empty plot that says so.
  if(is_data_logger_tube){
  
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
    filter(unixtime>min(corresp_fs_tube$unix_start, na.rm=T)-(60*180)&
             unixtime<=max(corresp_fs_tube$unix_stop, na.rm=T)+(60*180)) %>% 
    filter(!is.na(light)) %>% 
    ggplot(aes(x=datetime, y=light))+
    #ADDING labels fails (error when printing pdf within loop)
    # geom_label(data = data_logger_tube %>%  filter(!is.na(uniqID)) %>% group_by(uniqID) %>% summarise(datetime=mean(datetime,na.rm=T)), aes(x=datetime, label = uniqID, y=0),vjust=0.5,hjust = 1,angle = 90)+
    geom_line()+
    geom_point(data=subset(data_logger_tube, !is.na(condition_light)),aes(col=condition_light))+
    labs(col="Light?")+
    scale_x_datetime(breaks = "hour")+
    # coord_cartesian(ylim = c(-40000, NA)) +  
    #scale y adjustment fails!
    # scale_y_continuous(limits = c(-45000,
    #                               max(data_logger_tube %>%
    #                                     filter(unixtime>min(corresp_fs$unix_start, na.rm=T)-(60*180)&unixtime<=max(corresp_fs$unix_stop,na.rm=T)+(60*180)) %>% pull(light))))+
    ggtitle(paste(subsite, "SN:",SN_logger_tube))+
    theme_bw()
  }else{
    p<- ggplot()+ ggtitle(paste(subsite, "File not found for SN:",SN_logger_tube))+theme_bw()
  }
 
   # Store each subsite light data plot in the list
  plotslight[[subsite]] <- p
  
  }
  
#   }
# }

#save lightplot for every subsite: 
pdf(file = paste0(lightplots_path,"Lightsensor_per_subsite.pdf"),width = 15)  # Open PDF device

# Loop through the list of plots and print each plot
for (plot_name in names(plotslight)) {
  print(plotslight[[plot_name]])
}

dev.off()  # Close the PDF device



#Save logger_summary
write.csv(x = logger_summary, file = paste0(lifewatch_example_path,"logger_per_incubation.csv"),row.names = F)


####Data loggers ISSUES: ####

#Detail exploration per sn and subsites based on pdf: Lightsensor_per_subsite.pdf 

#S1-CA: sn10237447 A1, A2, R1 ~5min
#S1-CA: sn10257417 P1, P2,R2 ~5min

#S1-CU: sn21329647 A1,P1, ~5min. For R1 no data exist(confirmed with Benj)
#S1-CU: sn23524567 A2, P2, R2 no data for light intensity exists (confirmed with Benj)

#S1-DA: sn21245425  P2,R1,R2    NOT SYNC: TIMEZONES! "All sensors configured in local time (GMT+3)" note in excel
#S1-DA: sn21329647  P1,A1,A2    NOT SYNC: TIMEZONES! "All sensors configured in local time (GMT+3)" note in excel

#S1-DU: file not found for any (only 1 sensor-data in folder, is a floating SN). Rafa Carballeira was in charge of loggers, contacted on 20241220 for feedback.

#S1-RI: sn10257417 A1,P2,R1  NOT SYNC CHECK SN, DATES, TIMEZONES
#S1-RI: sn10257418  Good

#S1-VA: all good, some light_condition erroneous in fieldsheets

#S2-CA: sn10257417 P1,P2,R2 NOT SYNC: TIMEZONE issue most likely 
#S2-CA: sn10257418 A1, A2, R1 very fewpoints (logger in hour resolution)

#S2-CU: A1, A2 P1, P2, R1, R2 NOT SYNC CHECK SN, DATES, TIMEZONES

#S2-DA: sn10257417 P2,R1,R2: NOT SYNC CHECK SN, DATES, TIMEZONES
#S2-DA: sn10257418 A1,P1, very few points, A2 no data. (logger in hour resolution)

#S2-DU: A1, A2 P1, P2, R1, R2 NOT SYNC CHECK SN, DATES, TIMEZONES

#S2-RI: sn10257417 A1, P2, R1 ~5min
#S2-RI: sn10257418 A2, P1, R2 ~5min

#S2-VA: sn10257417 A1, P2, R1 ~1min
#S2-VA: sn10257418 A2, P1, R2 ~5min, 

#S3-CA: sn10257417 A1, A2, R1 ~5 min,
#S3-CA: sn10257418 P1, P2, R2 ~5min

#S3-CU: sn10257417 A2, P2, R2 ~5 min,
#S3-CU: sn10257418 A1, P1, R1 ~5 min,

#S3-DA: sn10257417 A1,A2,P1 ~5 min,
#S3-DA: sn10257418 R1,R2,P2 ~1 min,

#S3-DU: sn10257417 A1, A2, P2 ~5 min,
#S3-DU: sn10257418 R1, R2, P1 ~1 min,

#S3-RI: sn10257417 A2, P1, R2 ~1h delay TIMEZONE issue
#S3-RI: sn10257418 A1, P2, R1 ~1h delay TIMEZONE issue

#S3-VA: sn10257417 A1, P2, R1 ~1 min,
#S3-VA: sn10257418 A2, P1, R2 ~5 min,

#S4-CA: sn10257417 A1, P1, P2 NOT sync CHECK SN, DATES, TIMEZONES
#S4-CA: sn10257418 A2, R1, R2 NOT sync CHECK SN, DATES, TIMEZONES

#S4-CU: sn10257418 A1, R1, P1 ~1min
#S4-CU: sn10257417 A2, R2, P2 ~3min

#S4-DA: sn10257417 A1, A2, P1 ~5 min
#S4-DA: sn10257418 P2, R1, R2 on time perfect

#S4-DU: sn10257417 A1, A2, ~2 min delay, P2 unclear sync 
#S4-DU: sn10257418 R1, R2, P1 ~2 min 

#S4-RI: sn10257417 A1, P2, R1 ~2 min
#S4-RI: sn10257418 A2, P1, R2 ~1 min
 
#S4-VA: sn10257418 A1, A2, R1  NOT sync CHECK SN, DATES, TIMEZONES
#S4-VA: sn10257417 P1, P2, R2 NOT sync CHECK SN, DATES, TIMEZONES



#Notes from first execution of loop: 
  #subsites with light data duplicated: eg. S1-CU-A1 light data is re-used for S1-CU-A2 (with a different SN_logger and which data-logger file does not have light data). 1st. WHY IS THE DATA RE-USED? 2nd(less important)Why no data in file

  #Subsites with no data: 
    #No file found (is_data_logger_tube==F) check import of data logger files and the files themselves)
    #file exist (is_data_logger_tube==T) but empty light data in timewindow of incubations 

  #Subsites with data but only 6-7 data-points check og files (data in hourly resolution)



#Things to fix in loop: 


#Cannot find a way of consistently plotting the uniqID labels inside each plot (printing fails when there is an issue with the data)

#REVISE import of data loggers for timezone changes. 
#ASK Camille if he adjusted the times of LICOR PICARRO, LOS GATOS. 




###Format data_logger data####
logger_summary<- read.csv(paste0(lifewatch_example_path,"logger_per_incubation.csv"))

logger_summary %>%   dplyr::summarise(n = dplyr::n(), .by = c(uniqID)) |>
  dplyr::filter(n > 1L) 
#The loop mantains uniqueness for uniqID (no duplicates)


#WE HAVE TO INCLUDE a letter Denoting strata_short to plotcode to avoid duplicated records. 

logger_summary_formated<- logger_summary %>% 
  #Drop duplicated incubations (repetitions)
  filter(!uniqID %in%duplicate_incubations) %>% 
  separate(uniqID, into = c("campaign","case_pilot","status","plotnum","strat_short","light_short","startime"),sep = "-",remove = F) %>% 
  mutate(plotcode=toupper(paste(campaign, case_pilot,status,plotnum,strat_short,sep = "-")),#Add 1st letter of strata to plotcode.eg: S1-CA-R2-9-V
           light_condition=case_when(light_short=="d"~"dark",
                                     light_short=="t"~"transparent",
                                     TRUE~NA)) %>% 
  #Drop plots with NA in light_condition:
  filter(!is.na(light_condition)) %>% 
  #Select central measure for temperature and light data: median (for the moment)
  transmute(plotcode, light_condition, temperature=temp_median, light=light_median) %>% 
  pivot_wider(names_from = light_condition, values_from = c(temperature, light))


write.csv(logger_summary_formated, file = paste0(lifewatch_example_path, "logger_per_plot_formated.csv"),row.names = F)





#JOIN 3 per-plot datasets####

#Join GHG, Vegetation and data logger datasets into 1
ghg<- read.csv(paste0(lifewatch_example_path,"GHG_per_plot_formated.csv"))
veg<- read.csv(paste0(lifewatch_example_path,"vegetation_per_plot_formated.csv"))
logger<- read.csv(paste0(lifewatch_example_path,"logger_per_plot_formated.csv"))

#Add vegetation and logger data to GHG fluxes 
all<- ghg %>% 
  merge.data.frame(., veg, by="plotcode",all.x = T) %>% 
  merge.data.frame(., logger, by="plotcode",all.x = T)


#Save merged dataset 
write.csv(all, file = paste0(lifewatch_example_path,"GHG_plus_metadata_per_plot.csv"),row.names = F)




#Integrate day-night ####

#Integrate day-night per strata

#Calculate daylight hours using latitude and date (hours of light)
# install.packages("suncalc")
library(dplyr)
library(suncalc)

# get median coordinates per subsite (including all seasons)
samplings <- ghg %>% 
  select(subsite, latitude, longitude) %>% 
  group_by(subsite) %>% 
  summarise(subsite_latitude=median(latitude,na.rm=T),
            subsite_longitude=median(longitude,na.rm=T)) %>% 
  ungroup()


#Use dplyr and suncalc to calculate daylight hours for each subsite-date combination.
sampling_daylight <- samplings %>%
  merge.data.frame(y=ghg %>% select(subsite, date), by="subsite", all=T) %>% 
  distinct() %>% 
  mutate(date=as.Date(date)) %>% 
  rowwise() %>%
  mutate(
    sunlight_times = list(getSunlightTimes(date = date, lat = subsite_latitude, lon = subsite_longitude, keep = c("sunrise","sunset"))),
    daylight_duration = as.numeric(difftime(sunlight_times$sunset, sunlight_times$sunrise, units = "hours"))
  ) %>%
  ungroup() %>%
  select(subsite, date, subsite_latitude, subsite_longitude, daylight_duration)



#INTEGRATE to net_GHG exchange per plot (integrate dark and transparent with daylight_duration)
net_ghg_per_plot<- ghg %>% 
  mutate(date=as.Date(date)) %>% 
  merge.data.frame(sampling_daylight, by=c("subsite","date"), all = T) %>% 
  select(season, pilot_site, subsite, sampling, subsite_latitude, subsite_longitude, strata, date, daylight_duration, plotcode, latitude, longitude,  
         CO2_best.flux_dark, CO2_best.flux_transparent,
         CH4_diffusive_flux_dark, CH4_diffusive_flux_transparent,
         CH4_ebullitive_flux_dark, CH4_ebullitive_flux_transparent) %>% 
  #Get CH4 total flux (diffusion + ebullition)
  mutate(CH4_total_dark=CH4_diffusive_flux_dark+CH4_ebullitive_flux_dark,
         CH4_total_transparent=CH4_diffusive_flux_transparent+CH4_ebullitive_flux_transparent) %>% 
  #Integrate transparent and dark with daylight_duration (fluxlight*lighthours+fluxdark*darkhours)/24: same units for net flux (umol/m2/s), but integrated for the whole day. For open water and bare plots, use directly the dark flux
  #CO2:
  mutate(net_CO2=case_when(grepl("vegetated",strata)~((CO2_best.flux_transparent*daylight_duration)+(CO2_best.flux_dark*(24-daylight_duration)))/24,
                           grepl("open|bare",strata)~CO2_best.flux_dark),
  #CH4:    
         net_CH4=case_when(grepl("vegetated",strata)~((CH4_total_transparent*daylight_duration)+(CH4_total_dark*(24-daylight_duration)))/24,
                           grepl("open|bare",strata)~CH4_total_dark+CH4_total_dark)) %>% 
  #Discard _transparent/_dark ghg fluxes
  select(-c(CO2_best.flux_dark,CO2_best.flux_dark, CO2_best.flux_transparent,
            CH4_diffusive_flux_dark, CH4_diffusive_flux_transparent,
            CH4_ebullitive_flux_dark, CH4_ebullitive_flux_transparent,CH4_total_dark,CH4_total_transparent)) 


#Save net_ghg per plot
write.csv(net_ghg_per_plot, file=paste0(lifewatch_example_path,"net_ghg_per_plot.csv"),row.names = F)

#Average per subsite and strata
net_ghg_summary<- net_ghg_per_plot %>% 
  select(-c(plotcode,latitude, longitude)) %>% 
  group_by(season, pilot_site, subsite, ,subsite_latitude,subsite_longitude, sampling, date, strata,daylight_duration) %>% 
  summarise(avg_netCO2=mean(net_CO2, na.rm=T),
            sd_netCO2=sd(net_CO2, na.rm=T),
            n_netCO2=sum(!is.na(net_CO2)),
            avg_netCH4=mean(net_CH4, na.rm=T),
            sd_netCH4=sd(net_CH4, na.rm=T),
            n_netCH4=sum(!is.na(net_CH4)),
            .groups = "drop"
  )

#Save net_ghg per subsite-strata
write.csv(net_ghg_summary,file = paste0(lifewatch_example_path, "net_ghg_summary.csv"),row.names = F)

