
# ---
# Authors: Camille Minaudo
# Project: "RESTORE4Cs"
# date: "Oct 2023"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
# This script loads all raw measurement files (data level L0a) from one of the gas
# analyzers used in the project, transform it into a unified harmonized csv
# file, and saves a separate file + a pdf plot out of each incubation performed.


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

# ---- List GHG chamber fieldsheets in Dropbox ----
myfieldsheets_list <- list.files(fieldsheetpath, pattern = "Fieldsheet-GHG.xlsx", all.files = T, full.names = T, recursive = T)


# ---- Read all fieldsheets and put them in a single dataframe ----


# Read the first one to get the headers
fieldsheet_temp <- readxl::read_xlsx(myfieldsheets_list[1],
                                     col_names = T)
my_headers <- names(fieldsheet_temp)

# Go through them all and keep the info
isF <- T
for (f in myfieldsheets_list){
  fieldsheet_temp <- readxl::read_xlsx(f, col_names = F, range = "A3:V30")
  names(fieldsheet_temp) <- my_headers
  fieldsheet_temp <- fieldsheet_temp[!is.na(fieldsheet_temp$plot_id),]
  fieldsheet_temp$date <- as.Date( fieldsheet_temp$date, tryFormats = c("%d.%m.%Y", "%d/%m/%Y"))

  if(isF){
    isF <- F
    fieldsheet <- fieldsheet_temp
  } else {
    fieldsheet <- rbind(fieldsheet, fieldsheet_temp)
  }
}


fieldsheet <- fieldsheet[!is.na(fieldsheet$longitude),]




# ---- Import and store measurements to RData ----
data_folders <- list.dirs(datapath, full.names = T, recursive = F)[-2]
message("Here is the list of data folders in here:")
print(data_folders)


i_licor <- grep(pattern = "Licor", x = data_folders)
i_losgatos <- grep(pattern = "Los Gatos", x = data_folders)
i_picarro <- grep(pattern = "Picarro", x = data_folders)
map_analysers <- data.frame(instrument = c("LI-7810","LGR","G2508"),
                            i = c(i_licor, i_losgatos, i_picarro),
                            date.format = c("ymd","mdy","ymd"))
k=0
for (data_folder in data_folders){
  k=k+1
  setwd(data_folder)
  import2RData(path = data_folder, instrument = map_analysers$instrument[k],
               date.format = map_analysers$date.format[k], timezone = 'UTC')

  # load all these R.Data into a single dataframe
  file_list <- list.files(path = paste(data_folder,"/RData",sep=""), full.names = T)
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

  # get read of possible duplicated data
  mydata_imp_clean <- mydata_imp %>%
    group_by(DATE,TIME) %>% # the complete group of interest
    mutate(duplicate = n()) %>% # count number in each group
    filter(duplicate == 1) %>% # select only unique records
    select(-duplicate) # remove group count column


  # save this dataframe as a new RData file
  setwd(data_folder)
  save(mydata_imp_clean, file = "data_all_clean.RData")

}
















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

