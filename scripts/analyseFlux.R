
# ---
# Authors: Camille Minaudo
# Project: "RESTORE4Cs"
# date: "May 2024"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
# This script loads the fluxes estimates from script raw2flux.R. 
# it is used to explore and analyze the dataset.

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
require(dplyr)
require(purrr)
require(data.table)
require(tools)

repo_root <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
files.sources = list.files(path = paste0(repo_root,"/functions"), full.names = T)
for (f in files.sources){source(f)}


#################################
sampling <- "S2"
# USER, please specify if you want plots to be saved
harmonize2RData <- F
doPlot <- T
#################################



# ---- Directories and data loading ----
dropbox_root <- "C:/Users/Camille Minaudo/Dropbox/RESTORE4Cs - Fieldwork/Data" # You have to make sure this is pointing to the write folder on your local machine
results_path <- paste0(dropbox_root,"/GHG/Processed data/computed_flux/")
setwd(results_path)
listf <- list.files(path = results_path, pattern = ".csv", all.files = T, full.names = T, recursive = F)
table_results_all <- NULL
for (f in listf){
  table_results_all <- rbind(table_results_all, 
                             read.csv(file = f, header = T))
}

table_results_all$sampling <- str_sub(table_results_all$subsite, start = 1, 2)


# ---- quality of model fits ----

ggplot(table_results_all, aes())



# plot only water strata


table_results_all










# ---- Some plots ----


table_results_all$campaign_site <- substr(table_results_all$subsite,start = 1, stop = 5)
table_results_all$subsite_short <- substr(table_results_all$subsite,start = 7, stop = 8)


plt_CO2 <- ggplot(table_results_all, aes(subsite_short, CO2_LM.flux,
                                         fill = lightCondition))+
  geom_hline(yintercept = 0)+
  geom_boxplot(alpha=0.5)+
  # geom_jitter(width = 0.2, aes(colour = strata), size=2)+
  theme_article()+
  xlab("subsite")+
  ylab("CO2 flux mmol/m2/s")+
  ggtitle("CO2 flux")+
  scale_fill_viridis_d(begin = 0.2, end = 0.9)+
  scale_colour_viridis_d(begin = 0.2, end = 0.9, option = "C")+
  facet_grid(.~campaign_site)


plt_CH4diff <- ggplot(table_results_all, aes(subsite_short, CH4_LM.flux, 
                                             fill = lightCondition))+
  geom_hline(yintercept = 0)+
  geom_boxplot(alpha=0.5)+
  # geom_jitter(width = 0.2, aes(colour = strata), size=2)+
  theme_article()+
  xlab("subsite")+
  ylab("CH4 flux nmol/m2/s")+
  ggtitle("CH4 flux")+
  scale_fill_viridis_d(begin = 0.2, end = 0.9)+
  scale_colour_viridis_d(begin = 0.2, end = 0.9, option = "C")+
  facet_grid(.~campaign_site)

plt_all <- ggarrange(plt_CO2, plt_CH4diff, ncol = 1)

ggsave(plot = plt_all, filename = paste0(myfilename,".jpg"), path = results_path, 
       width = 10, height = 8, dpi = 300, units = 'in')



for (cs in unique(table_results_all$campaign_site)){
  
  table_results_cs <- table_results_all[table_results_all$campaign_site == cs,]
  
  plt_CO2 <- ggplot(table_results_cs, aes(subsite_short, CO2_LM.flux,
                                          fill = lightCondition))+
    geom_hline(yintercept = 0)+
    geom_boxplot(alpha=0.5, width=0.5)+
    geom_jitter(width = 0.1, size=2, alpha=0.5, aes(shape = strata))+
    theme_article()+
    xlab("subsite")+
    ylab("CO2 flux mmol/m2/s")+
    ggtitle(paste0(cs, ""))+
    scale_fill_viridis_d(begin = 0.2, end = 0.9)+
    scale_colour_viridis_d(begin = 0.2, end = 0.9, option = "C")+
    # facet_grid(lightCondition~.)+
    theme(legend.position ='none')
  
  
  plt_CH4diff <- ggplot(table_results_cs, aes(subsite_short, CH4_LM.flux, 
                                              fill = lightCondition))+
    geom_hline(yintercept = 0)+
    geom_boxplot(alpha=0.5, width=0.5)+
    geom_jitter(width = 0.1, size=2, alpha=0.5, aes(shape = strata))+
    theme_article()+
    xlab("subsite")+
    ylab("CH4 flux nmol/m2/s")+
    # ggtitle("CH4 flux")+
    scale_fill_viridis_d(begin = 0.2, end = 0.9)+
    scale_colour_viridis_d(begin = 0.2, end = 0.9, option = "C")
  
  plt_cs <- ggarrange(plt_CO2, plt_CH4diff, ncol = 2)
  
  
  setwd(results_path)
  myfilename <- paste(cs,"fluxes",min(as.Date(table_results_all$start.time)),"to",
                      max(as.Date(table_results_all$start.time)), sep = "_")
  
  ggsave(plot = plt_cs, filename = paste0(myfilename,".jpg"), path = results_path, 
         width = 10, height = 5, dpi = 300, units = 'in')
}
