# ---
# Authors: Miguel Cabrera
# Project: "RESTORE4Cs"
# date: "March 2025"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
# This script loads all the sediment fieldsheets, compiles the cores measurements and stores the resulting compilation in a single csv with the date of compilation.


rm(list = ls()) # clear workspace
cat("/014") # clear console


# ---- packages ----
library(tidyverse)
library(readxl)

# library(lubridate)
# library(zoo)
# library(ggplot2)
# library(grid)
# library(egg)
# library(sp)
# library(sf)
# library(openxlsx)

# source(paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/get_unix_times.R"))

# ---- Settings ----
#dropbox_root <- "C:/Users/Camille Minaudo/Dropbox/RESTORE4Cs - Fieldwork/Data/" #Camille Minaudo
#dropbox_root <- "D:/Dropbox/RESTORE4Cs - Fieldwork/Data/" #Benjamin Misteli PC
dropbox_root <- "C:/Users/Miguel/Dropbox/RESTORE4Cs - Fieldwork/Data/" #Miguel Cabrera 

fieldsheetsed_path <- paste0(dropbox_root,"Sediment/Fieldsheets/")
cores_folder <- paste0(dropbox_root, "Cores/")


# list sediment fieldsheets
f <- list.files(fieldsheetsed_path, pattern = "Fieldsheet", all.files = T, full.names = T, recursive = T)
i <- grep(pattern = ".xlsx", x = f)
myfiles <- f[i]
rm(i,f)

fieldsheets_list <- list()

#Check number of files per subsite sampling
fieldsheets_map<- data.frame(filepath=myfiles) %>% 
  mutate(filename=gsub(pattern = fieldsheetsed_path,x = filepath,replacement = ""),
         sampling=substr(filename, 1,6),
         filename=gsub(pattern = paste(sampling, collapse = "|"), "",filename),
         sampling_code=gsub("/","",sampling),
         subsite_code=substr(filename, 1,8)) %>% 
  select(filepath, filename, sampling_code,subsite_code)

#Any sampling missing fieldsheet? (i.e. less than 6 fieldsheets per sampling?)
fieldsheets_map %>% group_by(sampling_code) %>% mutate(nfiles=n()) %>% filter(nfiles<6) %>% 
  select(filename)
#S2-DA-R2 does not have sediment fieldsheet!
#OK, no sediments or cores were collected


#Load sediment fieldsheets:
fieldsheet_temp <- readxl::read_xlsx(myfiles[1], col_names = TRUE)
colnames_target <- tolower(names(fieldsheet_temp))

for (f in myfiles){
  df <- readxl::read_xlsx(f, skip = 2, col_names = colnames_target, n_max = 30)
  df$date <- as.Date(df$date, tryFormats = c("%d.%m.%Y", "%d/%m/%Y"))
  fieldsheets_list[[f]] <- df
}


complete_list <- do.call(rbind, fieldsheets_list) %>% 
  rownames_to_column(var = "filepath") %>% 
  mutate(filepath= sub("xlsx.*", "xlsx", filepath)) %>%  #remove anything after the end of the filepath (.xslx) 
  merge.data.frame(fieldsheets_map, by = "filepath", all = T) #add identifiers from fieldsheet filenames 

# rm(fieldsheet,fieldsheet_temp, fieldsheets_list)

#Compile all core sampling details
core_list<- complete_list %>% 
  filter(type=="core") %>% 
  mutate(sampling_time=format(sampling_time, format = "%H:%M:%S", tz="utc"),#format sampling time
         )

#Subsites without core details: 
complete_list %>% 
  filter(!is.na(pilot_site))  %>% 
  group_by(subsite_code) %>% 
  summarise(number_ofcores=sum(type=="core")) %>% 
  filter(number_ofcores==0) %>% 
  filter(!subsite_code%in%c("S3-DA-R2","S4-DA-R2"))#NO cores were obtained in DA-R2


#Samples with missing values (other than comments)
incomplete_details<- core_list[!complete.cases(core_list[,!names(core_list) %in% "comments"]),] %>% filter(!is.na(pilot_site)) 
message(paste("The following subsites have incomplete information:",unique( incomplete_details$subsite_code)))


#Check consistency in sample_ID: sample_ID subsite vs filename subsite
inconsistencies_core_id<- core_list %>% 
  separate(sample_id, into = c("core_season","core_site","core_subsite","coreid"),remove = F,sep = "-") %>% 
  filter(!is.na(core_site)) %>% 
  mutate(core_subsite_code=paste(core_season,core_site,core_subsite, sep = "-")) %>% 
  filter(core_subsite_code!=subsite_code)
if(length(unique(inconsistencies_core_id$sample_id))>0){message(paste("The following cores have inconsistent codes:", unique(inconsistencies_core_id$sample_id)))}


#Suspect units for water-depth (suposed to be in cm)

checked_waterdepth<- c("S1-VA-R2","S2-CA-P2","S2-CU-R2","S2-RI-P2","S3-RI-P1","S4-DU-R1","S4-VA-A2")

complete_list %>% 
  filter(water_depth<=2) %>% 
  filter(water_depth!=0) %>% 
  select(subsite_code, water_depth, comments) %>%
  filter(!subsite_code%in%checked_waterdepth) %>% 
  distinct()

#S1-VA-R2. water depth of 2cm are ok, checked
#S2-CA-P2. Water depth of 0.5cm are ok, checked
#S2-CU-R2. water depth of 2cm is ok, checked
#S2-RI-P2. Water depths of 0.5 and 2cm are ok, checked
#S3-RI-P1. Water depths of <2cm are ok, checked
#S4-DU-R1 water depths of <2cm are ok, checked
#S4-VA-A2 water depths of 2 cm are ok, checked

#S2-DU-R1 sediment sample with 1cm in plot 8 (no water in GHG plot 8), not important. 

#Format sampleID to be consistent:
core_list_formated<- core_list %>% 
  filter(!is.na(sample_id)) %>%
  separate(sample_id, into = c("core_season","core_site","core_subsite","coreid"),remove = F,sep = "-") %>% 
  mutate(corenum=parse_number(coreid),
         sample_id=paste0(core_season,"-",core_site,"-", core_subsite,"-C",corenum),
         subsite_code=paste0(core_season,"-",core_site,"-", core_subsite)) %>%
  rename(sampling_date=date, latitude=`gps latitude`,longitude=`gps longitude`, sampling_time_utc=sampling_time) %>% 
  select(pilot_site,subsite_code,sampling_date, sampling_time_utc, latitude, longitude, person_sampling, strata, water_depth, sample_id, comments,filename) %>% 
  #Harmonize naming of strata (small inconsistencies found in fieldsheets)
  mutate(strata=case_when(strata=="vegetation"~"vegetated",
                          strata%in%c("open","Open water")~"open water",
                          strata=="soil"~"bare",
                          TRUE~strata))


#Save core-fieldsheet compilation:
write.csv(core_list_formated, file = paste0(cores_folder, "Core_fieldsheet_compilation.csv"), row.names = F)


#Subset variables for Repository:
all_core_sampling_details<- core_list_formated %>% 
  rename(core_id=sample_id, water_depth_cm=water_depth) %>% 
  select(core_id, pilot_site, subsite_code, sampling_date, latitude, longitude, water_depth_cm)



#Save repository sampling details
write.csv(all_core_sampling_details, file = paste0(cores_folder,"/Cores_flux/All_core_sampling_details.csv"),row.names = F)
