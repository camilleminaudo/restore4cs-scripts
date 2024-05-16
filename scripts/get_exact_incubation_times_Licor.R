# ---
# Authors: Camille Minaudo
# Project: "RESTORE4Cs"
# date: "Nov 2023"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
# This script extracts the exact datetime for each chamber measurements performed
# with the LiCor. It loads the corresponding fieldworksheet and
# produce a new fieldwork sheet with corrected time stamps

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

# ---- functions ----


repo_root <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
files.sources = list.files(path = paste0(repo_root,"/functions"), full.names = T)
for (f in files.sources){source(f)}


# ---- Directories ----
dropbox_root <- "C:/Users/Camille Minaudo/Dropbox/RESTORE4Cs - Fieldwork/Data"
datapathRAW <- paste0(dropbox_root,"/GHG/RAW data/RAW Data Licor-7810")
datapath <- paste0(dropbox_root,"/GHG/RAW data/RAW Data Licor-7810/RData")
fieldsheetpath <- paste0(dropbox_root,"/GHG/Fieldsheets")
path2processed <- paste0(dropbox_root,"/GHG/Processed data")

# ---- List folders with Licor data in Dropbox ----
data_folders <- list.dirs(datapath, full.names = T, recursive = F)
subsites <- basename(data_folders)
print(subsites)

# go through the RAW data
fs <- list.files(path = datapathRAW, pattern = c(".txt", ".data"), full.names = T, recursive = T)
r <- grep(pattern = ".RData",x=fs)
fs <- fs[-r]


isF = T
for (f in fs){
  message(paste0("extracting info from ",basename(f)))
  data <- read_Licor(file = f)
  labels <- unique((data$label))
  labels <- labels[!is.na(labels)]
  # labels<- labels[-which(labels=="")]

  # get start and stop times for each label
  for (label in labels){
    ind_corresp <- which(data$label == label)
    
    start <- min(data$unixtime[ind_corresp])
    stop <- max(data$unixtime[ind_corresp])
    

    map_incubations_temp <- data.frame(date = first(data$date[ind_corresp]),
                                       label = label,
                                       start = start,
                                       stop = stop,
                                       time_start = strftime(start, format="%H:%M:%S", tz = 'utc'),
                                       time_stop = strftime(stop, format="%H:%M:%S", tz = 'utc'),
                                       folder = basename(dirname(f)),
                                       file = basename(f),
                                       path = f)
    if(isF){
      isF = F
      map_incubations <- map_incubations_temp
    } else {
      map_incubations <- rbind(map_incubations,map_incubations_temp)
    }
  }
}


# uniqueID <- paste0(map_incubations$date, map_incubations$label)
map_incubations <- map_incubations[!duplicated(map_incubations$start),]
map_incubations <- map_incubations[order(map_incubations$start),]

# Save incubtion map to Dropbox
write.csv(x = map_incubations, file = paste0(datapathRAW,"/map_incubations.csv"), row.names = F)



ggplot(map_incubations, aes(start, stop-start))+geom_point()+theme_article()+scale_y_log10()+geom_hline(yintercept = 900)


map_incubations <- map_incubations[which(map_incubations$stop-map_incubations$start < 15*60),]


# ---- Load fieldsheets ----
# list filenames
myfieldsheets_list <- list.files(fieldsheetpath, pattern = "Fieldsheet-GHG.xlsx", all.files = T, full.names = T, recursive = T)
# i <- grep(pattern = sampling, x = myfieldsheets_list) # selecting the files corresponding to the selected sampling campaign
# myfieldsheets_list <- myfieldsheets_list[i]
# Read all fieldsheets and put them in a single dataframe
fieldsheet <- read_GHG_fieldsheets(myfieldsheets_list)

fieldsheet$uniqID <- tolower(paste(fieldsheet$subsite,fieldsheet$plot_id,substr(fieldsheet$transparent_dark, 1, 1), sep = "-"))



# check the closest incubation in map_incubations for each row in fieldsheet.
# if more than 4 minutes apart, we consider the row in map_incubations out of sampling
fieldsheet_Licor <- fieldsheet[fieldsheet$gas_analyzer=="LI-COR",]
fieldsheet_Licor$corresponding_row <- NA
fieldsheet_Licor$label <- "none"
fieldsheet_Licor$unix_start_corr <- fieldsheet_Licor$unix_stop_corr <- NA


for (i in seq_along(fieldsheet_Licor$plot_id)){
  ind <- which.min(abs(fieldsheet_Licor$unix_start[i] - map_incubations$start))
  if(length(ind)>0){
    if(abs(fieldsheet_Licor$unix_start[i] - map_incubations$start[ind])<3*60){
      fieldsheet_Licor$corresponding_row[i] <- ind
      fieldsheet_Licor$label[i] <- map_incubations$label[ind]
      
      # replacing unix_startand unix_stop with new values
      fieldsheet_Licor$unix_start_corr[i] <- map_incubations$start[ind]
      fieldsheet_Licor$unix_stop_corr[i] <- map_incubations$stop[ind]
    }
  }
}


# # finding corresponding row in case of NA in corresponding_row
# if(sum(is.na(corresponding_row))>0){
#   ind_NAs <- which(is.na(corresponding_row))
#   # it is the most probable that the missing info is in-between two rows with all the data we need
#   interp <- approx(x = seq_along(fieldsheet_Licor$plot_id), y = corresponding_row, xout = ind_NAs)$y
#   is_integer <- (interp - floor(interp)) == 0
#   corresponding_row[ind_NAs][is_integer] <- interp[is_integer]
# }
if(sum(is.na(fieldsheet_Licor$corresponding_row))>0){
  ind_NAs <- which(is.na(fieldsheet_Licor$corresponding_row))
  message("Could not find a corresponding incubation map for the following incubations:")
  fieldsheet_Licor$uniqID[ind_NAs]
  # removing rows from fieldsheet_Licor with missing correspondance
  # fieldsheet_Licor <- fieldsheet_Licor[-ind_NAs,]
  # corresponding_row <- corresponding_row[-ind_NAs]
}


fieldsheet_Licor$error_time_start <- fieldsheet_Licor$unix_start_corr - fieldsheet_Licor$unix_start
fieldsheet_Licor$error_time_stop <- fieldsheet_Licor$unix_stop_corr - fieldsheet_Licor$unix_stop


p_start<- ggplot(fieldsheet_Licor, aes(error_time_start))+geom_density(fill="grey", alpha=0.2)+theme_article()+
  xlab("Start in REMARK - Start in fieldsheet [secs]")+
  geom_vline(xintercept = 0)+
  ggtitle("Difference in start time - LI-COR")

p_stop<- ggplot(fieldsheet_Licor, aes(error_time_stop))+geom_density(fill="grey", alpha=0.2)+theme_article()+
  xlab("Stop in REMARK - Stop in fieldsheet [secs]")+
  geom_vline(xintercept = 0)+
  ggtitle("Difference in stop time - LI-COR")


ggarrange(p_start, p_stop)


as.data.frame(fieldsheet_Licor[which(abs(fieldsheet_Licor$error_time_stop)>500),])



ggplot(fieldsheet_Licor,aes(unix_start, unix_start_corr))+geom_point()+theme_article()

p_start<- ggplot(fieldsheet_Licor,aes(date, unix_start_corr-unix_start))+geom_point()+theme_article()+
  ylab("Start in REMARK - Start in fieldsheet [secs]")+
  geom_hline(yintercept = 0)+
  ggtitle("Difference in start time - LI-COR")
p_start

p_stop<- ggplot(fieldsheet_Licor,aes(date, unix_stop_corr-unix_stop))+geom_point()+theme_article()+
  ylab("Corrected stop - Start in fieldsheet [secs]")+geom_hline(yintercept = 0)+
  ggtitle("stop")


ggarrange(p_start, p_stop)


fieldsheet_Licor$duration_fieldsheet <- fieldsheet_Licor$unix_stop-fieldsheet_Licor$unix_start
fieldsheet_Licor$duration_corrected <- fieldsheet_Licor$unix_stop_corr-fieldsheet_Licor$unix_start_corr

ggplot(fieldsheet_Licor,aes(duration_fieldsheet, duration_corrected))+geom_point()+theme_article()+geom_abline(slope = 1, intercept = 0)+
  xlab("Duration in fieldsheet [secs]")+
  ylab("Duration in REMARK [secs]")

as.data.frame(fieldsheet_Licor[which(abs(fieldsheet_Licor$duration_corrected)<200),])



# 
# # k=0
# for (subsite in subsites){
#   message(paste0("processing data for ",subsite))
#   # k=k+1
#   # data_folder <- data_folders[k]
# 
#   # load corresponding fieldsheet
#   pilotsite <- substr(subsite, 1, 5)
#   f_name <- paste0(fieldsheetpath,"/",pilotsite,"/",subsite,"-Fieldsheet-GHG.xlsx")
#   message(f_name)
# 
#   fieldsheet <- readxl::read_xlsx(f_name, col_names = T)# Read it first to get the headers
#   my_headers <- names(fieldsheet)
#   fieldsheet <- readxl::read_xlsx(f_name, col_names = F, range = "A3:V50")# re-read and affect correctly the headers
#   names(fieldsheet) <- my_headers
#   fieldsheet <- fieldsheet[!is.na(fieldsheet$plot_id),]
#   fieldsheet$date <- as.Date( fieldsheet$date, tryFormats = c("%d.%m.%Y", "%d/%m/%Y"))
# 
#   fieldsheet$unix_start <- get_unix_times(mydate = fieldsheet$date, mytime = fieldsheet$start_time)
#   fieldsheet$unix_stop <- get_unix_times(mydate = fieldsheet$date, mytime = fieldsheet$end_time)
# 
# 
#   # check the closest incubation in map_incubations for each row in fieldsheet.
#   # if more than 3 minutes apart, we consider the row in map_incubations out of sampling
# 
#   corresponding_row <- NA*fieldsheet$plot_id
# 
#   for (i in seq_along(fieldsheet$plot_id)){
#     ind <- which.min(abs(fieldsheet$unix_start[i] - map_incubations$start))
#     # corresponding_row[i] <- ind
#     if(abs(fieldsheet$unix_start[i] - map_incubations$start[ind])<3*60){corresponding_row[i] <- ind}
#   }
#   isNA <- !is.na(corresponding_row)
#   corresponding_row <- corresponding_row[!is.na(corresponding_row)]
#   my_map_incubations <- map_incubations[corresponding_row, ]
# 
#   if(dim(my_map_incubations)[1] >0){
#     message("....found corresponding labels for this one")
#     fieldsheet <- fieldsheet[isNA,]
# 
#     fieldsheet$start_time <- strftime(fieldsheet$start_time, format="%H:%M:%S", tz = 'utc')
#     fieldsheet$end_time <- strftime(fieldsheet$end_time, format="%H:%M:%S", tz = 'utc')
# 
#     fieldsheet$unix_start_corrected <- my_map_incubations$start
#     fieldsheet$unix_stop_corrected <- my_map_incubations$stop
# 
#     fieldsheet$timestamp_start <- as.POSIXct(fieldsheet$unix_start_corrected, tz = "UTC", origin = "1970-01-01")
#     fieldsheet$timestamp_stop <- as.POSIXct(fieldsheet$unix_stop_corrected, tz = "UTC", origin = "1970-01-01")
# 
#     fieldsheet$start_time_corrected <- strftime(fieldsheet$timestamp_start, format="%H:%M:%S", tz = 'utc')
#     fieldsheet$end_time_corrected <- strftime(fieldsheet$timestamp_stop, format="%H:%M:%S", tz = 'utc')
# 
#     setwd(paste0(path2processed,"/corrected_fieldsheets"))
#     f_name <- paste0(subsite,"-Fieldsheet-GHG_corrected.csv")
#     write.csv(file = f_name, x = fieldsheet, row.names = F)
# 
#   }
# 
# 
# }





