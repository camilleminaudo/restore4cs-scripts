
# ---
# Authors: Miguel Cabrera
# Project: "RESTORE4Cs"
# date: "Oct 2024"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
##Working script for decission making into how to process the peaks obtain with the Valencia Li-COR 
#We will use one of the last rawfiles created as an example to build the processing pipeline

#Steps: 

#1 Detect peaks (not always same number of peaks with a single remark, 2-5)
#2 Define width of peaks (generally well constrained, but no unique width)
#3 Integrate area of peak
#4 Identify baseline
#5 Substract baseline from area of peaks


#Generally, much more clear signal for CH4 than for CO2 (much bigger difference between baseline and peaks). Peaks of CH4 and CO2 are simultaneous (make sure!!), so we will define the peaks based on the CH4 signal alone, then pull the CO2 values corresponding to CH4 peak. For the baseline we cannot use the same approach, as CO2 baseline is much more noisy, with some Humps in between and throughout the peaks. 

rm(list = ls()) # clear workspace
cat("/014") # clear console


# ---- packages ----
library(tidyverse)
library(readxl)
library(lubridate)
# library(zoo)
# library(ggplot2)
# library(grid)
# library(egg)
# library(goFlux)
# require(dplyr)
# require(purrr)
# require(msm)
# require(data.table)
# require(tools)
# 
# require(pbapply)

#Browse restore4cs-scripts functions:
repo_root <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
files.sources = list.files(path = paste0(repo_root,"/functions"), full.names = T)

files.sources

#Load functions needed: 

source(files.sources[grep("read_Licor",files.sources)])



# ---- Directories ----
dropbox_root <- "C:/Users/Miguel/Dropbox/RESTORE4Cs - Fieldwork/Data" # You have to make sure this is pointing to the write folder on your local machine
datapath_licor <- paste0(dropbox_root,"/GHG/RAW data/RAW Data Licor-7810")


# load incubation map (run get_exact_incubation_times_Licor.R to update it)
map_incubations<- read.csv(file = paste0(datapath_licor, "/map_incubations.csv"))

# We will eventually need a complete and exclusive list of all cores injection remarks
#For the moment, use an example (day with only core incubations)

example_file<- paste0(dropbox_root, "/GHG/RAW data/RAW Data Licor-7810/S3-CA/TG10-01275-2024-04-26T060000.data.txt")


#Import using camille's function
a<- read_Licor(example_file)


a %>% filter(label!="") %>% 
  group_by(label) %>% 
  summarise(nobs=n()) %>% 
  ggplot(aes(x=1,y=nobs))+
  geom_boxplot()+
  geom_point()

#Typically, 3-4 peaks during a ~60s remark duration (40-120s)
#max intensity of peaks should be reliably among the top 0.8 percentile

a %>% 
  mutate(datetime=ymd_hms(paste(date,UTCtime))) %>% 
  filter(label%in%unique(a$label)[2:12])%>% 
  ggplot(aes(x=datetime))+
  geom_point(aes(y=CH4*1000,col="ch4"))+
  geom_line(aes(y=CH4*1000,col="ch4"))+
  geom_point(aes(y=CO2*4,col="co2"))+
  geom_line(aes(y=CO2*4,col="co2"))


#We can use the findpeaks within the pracma package to find peaks and filter setting the minimum peak height as greater than the 0.8 percentile.
#We will have to manually explore peak identification performance within our dataset


library(pracma)

# Create a sample time series
time_series_data <- a %>%   filter(label==unique(a$label)[3]) %>% pull(CH4)

# Find peaks
peaks <- findpeaks(time_series_data,minpeakheight = quantile(time_series_data, 0.85))

# Print peaks
print(peaks)

# Visualize the data and detected peaks
plot(time_series_data, type = 'b', main = "Peak Detection Example", xlab = "Index", ylab = "Value")
points(peaks[,2], peaks[,1], col = "red", pch = 19)  # Mark peaks in red





#Other options: 

# The zoo package provides methods for dealing with irregular time series data. While it doesn’t directly detect peaks, you can apply smoothing functions and then use which.max() or custom logic to identify peaks.
library(zoo)

data <- c(1, 3, 7, 1, 2, 5, 8, 6, 4, 9, 3, 2)
smoothed <- rollmean(data, k = 3, fill = NA, align = "center")
peaks <- which(diff(sign(diff(smoothed))) == -2) + 1  # Local maxima



# 5. Signal Package
# The signal package provides tools for signal processing, including peak detection through convolution or filtering techniques.
# Example 
library(signal)

data <- c(1, 3, 7, 1, 2, 5, 8, 6, 4, 9, 3, 2)
smoothed <- sgolayfilt(data, p = 3, n = 11)  # Savitzky-Golay filter
peaks <- findpeaks(smoothed)
