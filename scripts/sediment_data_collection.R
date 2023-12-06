# ---
# Authors: Benjamin Misteli
# Project: "RESTORE4Cs"
# date: "Oct 2023"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
# This script loads all the fieldsheets contained in the Dropbox folder, and
# make a shapefile out of it.


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
library(sp)
library(sf)
library(openxlsx)

source(paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/get_unix_times.R"))

# ---- Settings ----
#dropbox_root <- "C:/Users/Camille Minaudo/Dropbox/RESTORE4Cs - Fieldwork/Data/" #Camille Minaudo
#dropbox_root <- "D:/Dropbox/RESTORE4Cs - Fieldwork/Data/Sediment/" #Benjamin Misteli PC
dropbox_root <- "C:/Users/misteli/Dropbox/RESTORE4Cs - Fieldwork/Data/Sediment" #Benjamin Misteli Laptop
sampling <- "S1"
setwd(dropbox_root)

# list files in Dropbox
f <- list.files(dropbox_root, pattern = "Fieldsheet", all.files = T, full.names = T, recursive = T)
i <- grep(pattern = ".xlsx", x = f)
f <- f[i]

myfiles <- f
fieldsheets_list <- list()

for (f in myfiles){
  
  # load file
  fieldsheet_temp <- readxl::read_xlsx(f,col_names = T)
  fieldsheet <- readxl::read_xlsx(f,skip = 2, col_names = F, n_max = 100)
  names(fieldsheet) <- tolower(names(fieldsheet_temp))
  fieldsheet$date <- as.Date(fieldsheet$date, tryFormats = c("%d.%m.%Y", "%d/%m/%Y"))
  fieldsheets_list[[f]] <- fieldsheet
}
complete_list <- do.call(rbind, fieldsheets_list)
row.names(complete_list) <- NULL
complete_list_sediments <- complete_list %>%
  filter(type == "sediment sample")

datasheet_wcl_names <- readxl::read_xlsx("Lab data/Datasheet-Sediment.xlsx",col_names = T)
datasheet_wcl <- readxl::read_xlsx("Lab data/Datasheet-Sediment.xlsx",skip = 2,col_names = F)
names(datasheet_wcl) <- tolower(names(datasheet_wcl_names))

datasheet_wcl <- datasheet_wcl[, -c(2:12)]
datasheet_wcl <- datasheet_wcl[!is.na(datasheet_wcl[, 2]), ]

merged_data <- merge(complete_list_sediments, datasheet_wcl, by = "sample_id", all = TRUE)
merged_data <- merged_data %>%
  arrange(desc(package_arrived))

write.xlsx(merged_data, file = "Lab data/Datasheet-Sediment.xlsx")


# Define the columns for which you want to create boxplots
columns_of_interest <- c("oc_content", "dry_weight_percentage", "ash_free_dry_mass_percentage")  # Add columns as needed

# Extract unique values of 'pilot_site'
unique_pilot_site <- unique(merged_data$pilot_site)

for (column in columns_of_interest) {
  plot_list <- list()
  
  for (CP in unique_pilot_site) {
    # Filter data for the current pilot_site and column
    country_data <- filter(merged_data, pilot_site == CP)
    
    # Create a plot for the current column and pilot_site
    boxplot <- ggplot(data = country_data, aes(x = subsite, y = .data[[column]])) + 
      geom_boxplot() + 
      theme(text = element_text(size = 16)) + 
      scale_color_viridis_d() +
      ggtitle(paste("Case Pilot:", CP))
    
    plot_list[[as.character(CP)]] <- boxplot
  }
  
  # Arrange the combined plots for the current column
  combined_plots <- ggpubr::ggarrange(plotlist = plot_list, ncol = 2, nrow = 3, common.legend = TRUE, legend = "right")
  
  # Save the combined plots for the current column
  png(filename = paste0(column, "_Boxplots.png"), width = 12, height = 10, units = "in", res = 300)
  print(combined_plots)  # Use print() to display and save the plot
  dev.off()
}

