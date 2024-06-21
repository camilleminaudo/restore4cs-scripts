
# ---
# Authors: Camille Minaudo
# Project: "RESTORE4Cs"
# date: "Oct 2023"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
# This script 

rm(list = ls()) # clear workspace
cat("/014") # clear console


# ---- Directories ----

# You have to make sure this is pointing to the write folder on your local machine
dropbox_root <- "C:/Users/Camille Minaudo/Dropbox/RESTORE4Cs - Fieldwork/Data" 


datapath <- paste0(dropbox_root,"/GHG/RAW data")
loggerspath <- paste0(datapath,"/RAW Data Logger")
RData_path <- paste0(dropbox_root,"/GHG/Processed data/RData/")
plots_path <- paste0(dropbox_root,"/GHG/GHG_expert_vs_automated/plots/")
results_path <- paste0(dropbox_root,"/GHG/GHG_expert_vs_automated/results/")

# ---- packages ----
library(tidyverse)
library(lubridate)
library(zoo)
library(ggplot2)
library(egg)
library(goFlux)
require(dplyr)
require(purrr)
require(msm)
require(data.table)
require(tools)
require(pbapply)


list_f <- list.files(path = results_path, pattern = "BLIND_vs_EXPERT_co2_ch4_fluxes_", full.names = T)
list_f <- list_f[grep(pattern = "csv", x = list_f)]

isF <- T
for (f in list_f){
  table_results.tmp <- read.csv(f)
  if(isF){
    isF <- F
    table_results_all <- table_results.tmp
  } else {
    table_results_all <- rbind(table_results_all,table_results.tmp)
  }
}

p_auto_vs_manual <- ggplot(data = table_results_all)+
  geom_abline(slope = 0,intercept = 0, color = 'black')+
  # geom_segment(data = data.frame(UniqueID = CO2_flux_res_auto$UniqueID,
  #                                meth1 = CO2_flux_res_auto$best.flux,
  #                                meth2 = CO2_flux_res_manID$best.flux), aes(x=UniqueID, xend=UniqueID, y = meth1, yend = meth2), linewidth=1, alpha = 0.5)+
  # aes(x=reorder(X_Variable, -Y_Variable, FUN=mean), y=Y_Variable)
  geom_point(aes(reorder(UniqueID, -best.flux, FUN=mean), best.flux, colour = flux_method), size=4, alpha = 0.5)+
  ylab("flux [(mmolCO2 or nmolCH4)/m2/s]")+
  theme_bw()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_colour_viridis_d(begin = 0.1, end = 0.9, option = "F")+facet_grid(.~variable, scales = 'free')+coord_flip()

p_auto_vs_manual
# 
# ggsave(filename = "BLIND_vs_EXPERT_co2_ch4_all.jpeg", plot = p_auto_vs_manual, path = plots_path, width = 10, height = 10, units = 'in', dpi = 300)
# 
# 
# 
# 
# table_results_sprd <- table_results_all[,c("variable","UniqueID","flux_method","best.flux","timestamp_processing")] %>%
#   pivot_wider(names_from = flux_method, values_from = c(best.flux))
# 
# ggplot(data = table_results_sprd)+
#   # geom_abline(slope = 0,intercept = 0, color = 'lightgrey')+
#   # geom_segment(data = data.frame(UniqueID = CO2_flux_res_auto$UniqueID,
#   #                                meth1 = CO2_flux_res_auto$best.flux,
#   #                                meth2 = CO2_flux_res_manID$best.flux), aes(x=UniqueID, xend=UniqueID, y = meth1, yend = meth2), linewidth=1, alpha = 0.5)+
#   geom_violin(aes(variable,abs((Blind-Expert)/Expert)*100), alpha = 0.5)+
#   geom_jitter(aes(variable,abs((Blind-Expert)/Expert)*100), alpha = 0.5, width=0.2)+
#   ylab("CO2 flux relative difference [%]")+
#   theme_bw()+
#   scale_y_log10()+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 
# 
# 
# list_f <- list.files(path = results_path, pattern = "BLIND_vs_EXPERT_ch4_ebullition_", full.names = T)
# list_f <- list_f[grep(pattern = "csv", x = list_f)]
# 
# isF <- T
# for (f in list_f){
#   table_results.tmp <- read.csv(f)
#   if(isF){
#     isF <- F
#     table_results_all <- table_results.tmp
#   } else {
#     table_results_all <- rbind(table_results_all,table_results.tmp)
#   }
# }
# 
# 
# 
# 
# p_auto_vs_manual_ch4_ebullition <- ggplot(data = table_results_all)+
#   geom_abline(slope = 0,intercept = 0, color = 'lightgrey')+
#   geom_segment(data = data.frame(UniqueID = CH4_res_meth1$UniqueID,
#                                  meth1 = CH4_res_meth1$ebullition,
#                                  meth2 = CH4_res_meth2$ebullition), aes(x=UniqueID, xend=UniqueID, y = meth1, yend = meth2), linewidth=1, alpha = 0.5)+
#   geom_point(aes(reorder(UniqueID, -ebullition, FUN=mean), ebullition, colour = flux_method), size=4, alpha = 0.5)+
#   
#   ylab("ebullition component [nmol/m2/s]")+
#   theme_bw()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_colour_viridis_d(begin = 0.1, end = 0.9, option = "F")+coord_flip()
# 
# p_auto_vs_manual_ch4_ebullition
# ggsave(filename = "BLIND_vs_EXPERT_ch4_ebullition_all.jpeg", plot = p_auto_vs_manual_ch4_ebullition, path = plots_path, width = 10, height = 6, units = 'in', dpi = 300)




