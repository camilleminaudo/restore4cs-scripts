# ---
# Authors: Camille Minaudo
# Project: "RESTORE4Cs"
# date: "Oct 2023"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
# This script lists all incubations corresponding to a suite of criterias and create a single RData file out of it
# This data will later be used for automatic versus expert flux assessment.

rm(list = ls()) # clear workspace
cat("/014") # clear console


# --- install goFlux if needed ---
#library(devtools)
#install_github("Qepanna/goFlux")

# ---- packages ----
library(tidyverse)
library(readxl)
library(lubridate)
library(zoo)
library(ggplot2)
library(grid)
library(egg)
library(goFlux)
require(dplyr)
require(purrr)
require(msm)
require(data.table)
require(tools)

require(pbapply)

repo_root <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
files.sources = list.files(path = paste0(repo_root,"/functions"), full.names = T)
for (f in files.sources){source(f)}


#################################
# SPECIFY HERE YOUR NAME
username <- "Camille"

nb_draw <- 10

# You have to make sure this is pointing to the write folder on your local machine
dropbox_root <- "C:/Users/Camille Minaudo/Dropbox/RESTORE4Cs - Fieldwork/Data" 

############################


# ---- Directories ----
datapath <- paste0(dropbox_root,"/GHG/RAW data")
loggerspath <- paste0(datapath,"/RAW Data Logger")
RData_path <- paste0(dropbox_root,"/GHG/Processed data/RData/")
plots_path <- paste0(dropbox_root,"/GHG/GHG_expert_vs_automated/plots/")
results_path <- paste0(dropbox_root,"/GHG/GHG_expert_vs_automated/results/")



# loading auxfile

setwd(dirname(results_path))

auxfile <- read.csv(file = "auxfile.csv")

auxfile$start.time <- as.POSIXct(auxfile$start.time, tz = 'UTC') 



# ----- Random draw of incubations -----

# draw a few incubations randomly
draw <- sample(seq_along(auxfile$subsite), nb_draw)


table_draw <- data.frame(username = username,
                         draw = draw,
                         subsite = auxfile$subsite[draw],
                         UniqueID = auxfile$UniqueID[draw])

myauxfile <- auxfile[draw,]
myauxfile$username <- username

# ---- Load incubation timeseries ----

mydata_all <- NULL
for(subsite in unique(myauxfile$subsite)){
  message("Loading data for ",subsite)
  
  gas <- unique(myauxfile$gas_analiser[which(myauxfile$subsite==subsite)])
  
  for (ga in gas){
    setwd(RData_path)
    if(ga== "LI-COR"){
      gs_suffix <- "LI-7810"
    } else {
      gs_suffix <- ga
    }
    
    load(file = paste0(subsite,"_",gs_suffix,".RData"))
    mydata <- mydata[,c("POSIX.time", 
                        "CO2dry_ppm", "CH4dry_ppb", "H2O_ppm",
                        "CO2_prec",   "CH4_prec",   "H2O_prec"  )]
    
    mydata_all <- rbind(mydata_all, mydata)
    rm(mydata)
  }
}





# ---- Process incubations and compute fluxes ----

# Define the measurements' window of observation

mydata_ow <- obs.win(inputfile = mydata_all, auxfile = myauxfile,
                     obs.length = myauxfile$duration, shoulder = 2)



# ----------- Compute fluxes after manual selection of CO2 data

# Manually identify measurements by clicking on the start and end points
mydata_manID <- lapply(seq_along(mydata_ow), click.peak.loop,
                       flux.unique = mydata_ow,
                       gastype = "CO2dry_ppm",
                       plot.lim = c(200,1000)) %>%
  map_df(., ~as.data.frame(.x))


# Additional auxiliary data required for flux calculation.
mydata_manID <- mydata_manID %>%
  left_join(myauxfile %>% select(UniqueID, Area, Vtot, Tcham, Pcham))


# start <- myauxfile$start.time
# start_corr <- mydata_manID$start.time_corr[match(x = myauxfile$UniqueID, mydata_manID$UniqueID)]
# cbind(start, start_corr)

table_draw$start.time_auto <- as.numeric(myauxfile$start.time)
table_draw$start.time_expert_co2 <- as.numeric(mydata_manID$start.time_corr[match(x = table_draw$UniqueID, mydata_manID$UniqueID)])
table_draw$end.time_auto <- as.numeric(myauxfile$start.time+myauxfile$duration)
table_draw$end.time_expert_co2 <- as.numeric(mydata_manID$end.time[match(x = table_draw$UniqueID, mydata_manID$UniqueID)])




# Add instrument precision for each gas
mydata_manID <- mydata_manID %>%
  mutate(CO2_prec = first(mydata_all$CO2_prec), CH4_prec = first(mydata_all$CH4_prec), 
         N2O_prec = first(mydata_all$N2O_prec), H2O_prec = first(mydata_all$H2O_prec))


# Calculate fluxes
CO2_results_manID <- goFlux(mydata_manID, "CO2dry_ppm")
H2O_results_manID <- goFlux(mydata_manID, "H2O_ppm")
CH4_results_manID <- goFlux(mydata_manID, "CH4dry_ppb")

# Use best.flux to select the best flux estimates (LM or HM)
# based on a list of criteria
criteria <- c("g.factor", "kappa", "MDF", "R2", "SE.rel")

CO2_flux_res_manID <- best.flux(CO2_results_manID, criteria)
H2O_flux_res_manID <- best.flux(H2O_results_manID, criteria)
CH4_flux_res_manID <- best.flux(CH4_results_manID, criteria)

# ----------- Compute fluxes blindly without any manual selection

mydata_ow <- obs.win(inputfile = mydata_all, auxfile = myauxfile,
                     obs.length = myauxfile$duration, shoulder = 2)

# Join mydata_ow with info on start end incubation
mydata_auto <- lapply(seq_along(mydata_ow), join_auxfile_with_data.loop, flux.unique = mydata_ow) %>%
  map_df(., ~as.data.frame(.x))

# Additional auxiliary data required for flux calculation.
mydata_auto <- mydata_auto %>%
  left_join(myauxfile %>% select(username, UniqueID, Area, Vtot, Tcham, Pcham))

# Add instrument precision for each gas
mydata_auto <- mydata_auto %>%
  mutate(CO2_prec = first(mydata_all$CO2_prec), CH4_prec = first(mydata_all$CH4_prec), 
         N2O_prec = first(mydata_all$N2O_prec), H2O_prec = first(mydata_all$H2O_prec))


# Calculate fluxes
CO2_results_auto <- goFlux(dataframe = mydata_auto, gastype = "CO2dry_ppm")
H2O_results_auto <- goFlux(mydata_auto, "H2O_ppm")
CH4_results_auto <- goFlux(mydata_auto, "CH4dry_ppb")

# Use best.flux to select the best flux estimates (LM or HM)
# based on a list of criteria
criteria <- c("g.factor", "kappa", "MDF", "R2", "SE.rel")

CO2_flux_res_auto <- best.flux(CO2_results_auto, criteria)
H2O_flux_res_auto <- best.flux(H2O_results_auto, criteria)
CH4_flux_res_auto <- best.flux(CH4_results_auto, criteria)


# Use best.flux to select the best flux estimates (LM or HM)
# based on a list of criteria
criteria <- c("g.factor", "kappa", "MDF", "R2", "SE.rel")

CO2_flux_res_auto <- best.flux(CO2_results_auto, criteria)
H2O_flux_res_auto <- best.flux(H2O_results_auto, criteria)
CH4_flux_res_auto <- best.flux(CH4_results_auto, criteria)


CO2_flux_res_auto$flux_method <- "Blind"
CO2_flux_res_manID$flux_method <- "Expert"

CH4_flux_res_auto$flux_method <- "Blind"
CH4_flux_res_manID$flux_method <- "Expert"




# ----------- Estimate CH4 diffusion/ebullition contributions
# method 1: density of prob. of first derivative
# method 2: manual selection of linear pattern in CH4

# estimate ch4 diffusion and ebullition components---------| METHOD 1 |-----------
CH4_res_meth1 <- CH4_flux_res_auto
CH4_res_meth1$total_estimated <- NA
CH4_res_meth1$ebullition <- NA
CH4_res_meth1$diffusion <- NA


for (i in seq_along(table_draw$UniqueID)){
  my_myauxfile <- myauxfile[myauxfile$UniqueID==table_draw$UniqueID[i],]
  
  my_incub <- mydata_all[as.numeric(mydata_all$POSIX.time)> table_draw$start.time_expert_co2[i] &
                           as.numeric(mydata_all$POSIX.time)< table_draw$end.time_expert_co2[i],]
  my_incub <- my_incub[!is.na(my_incub$CO2dry_ppm),]
  # calling dedicated function
  df_ebull <- separate_ebullition_from_diffusion(my_incub = my_incub, UniqueID = table_draw$UniqueID[i], doPlot = T)
  # computing fluxes
  H2O_mol = my_incub$H2O_ppm / (1000*1000)
  myfluxterm <- flux.term(my_myauxfile$Vtot, my_myauxfile$Pcham, my_myauxfile$Area,
                          my_myauxfile$Tcham, first(H2O_mol))
  CH4_flux_total <- df_ebull$delta_ch4/df_ebull$duration*myfluxterm # nmol/m2/s
  CH4_flux_diff <- df_ebull$avg_diff_slope*myfluxterm # nmol/m2/s
  CH4_flux_ebull <- CH4_flux_total - CH4_flux_diff
  
  CH4_res_meth1$total_estimated[which(CH4_res_meth1$UniqueID==table_draw$UniqueID[i])] <- CH4_flux_total
  CH4_res_meth1$ebullition[which(CH4_res_meth1$UniqueID==table_draw$UniqueID[i])] <- CH4_flux_ebull
  CH4_res_meth1$diffusion[which(CH4_res_meth1$UniqueID==table_draw$UniqueID[i])] <- CH4_flux_diff
}






# Estimate ch4 diffusion and ebullition components ---------| METHOD 2 |-----------


# Manually identify diffusive (more or less linear) CH4 behaviors by clicking on the start and end points
myCH4_diffusion <- lapply(seq_along(mydata_ow), click.peak.loop,
                          flux.unique = mydata_ow,
                          gastype = "CH4dry_ppb",
                          plot.lim = c(1900,max(mydata_all$CH4dry_ppb, na.rm = T)*1000)) %>%
  map_df(., ~as.data.frame(.x))


table_draw$start.time_expert_ch4 <- as.numeric(myCH4_diffusion$start.time_corr[match(x = table_draw$UniqueID, myCH4_diffusion$UniqueID)])
table_draw$end.time_expert_ch4 <- as.numeric(myCH4_diffusion$end.time[match(x = table_draw$UniqueID, myCH4_diffusion$UniqueID)])


# Add instrument precision for each gas
myCH4_diffusion <- myCH4_diffusion %>%
  mutate(CO2_prec = first(mydata_all$CO2_prec), CH4_prec = first(mydata_all$CH4_prec), 
         N2O_prec = first(mydata_all$N2O_prec), H2O_prec = first(mydata_all$H2O_prec))




# Calculate fluxes for CH4
CH4_results_diffusion <- goFlux(myCH4_diffusion, "CH4dry_ppb")

# Use best.flux to select the best flux estimates (LM or HM)
# based on a list of criteria
criteria <- c("g.factor", "kappa", "MDF", "R2", "SE.rel")

CH4_res_diff <- best.flux(CH4_results_diffusion, criteria)


CH4_res_meth2 <- CH4_flux_res_manID
CH4_res_meth2$total_estimated <- NA
CH4_res_meth2$ebullition <- NA
CH4_res_meth2$diffusion <- CH4_res_diff$best.flux


# Estimating ebullition component
for (id in unique(CH4_res_meth2$UniqueID)){
  i <- which(CH4_res_meth2$UniqueID == id)
  
  CH4_final <- CH4_flux_res_manID$Ct[i]
  CH4_initial <-  CH4_flux_res_manID$C0[i]
  incubation_time <- myauxfile$duration[which(myauxfile$UniqueID == id)]
  CH4_res_meth2$total_estimated[i] <- (CH4_final-CH4_initial)/incubation_time*CH4_flux_res_manID$flux.term[i] # nmol/m2/s
  CH4_res_meth2$ebullition[i] <- CH4_res_meth2$total_estimated[i]  - CH4_res_meth2$diffusion[i] # total flux - diffusive term
}


# compare method 1 and 2 to estimate ebullitive contribution
CH4_res_meth1$flux_method <- "dydt"
CH4_res_meth2$flux_method <- "Expert"


CO2_flux_res_auto$variable <- "CO2"
CO2_flux_res_manID$variable <- "CO2"
CH4_flux_res_auto$variable <- "CH4"
CH4_flux_res_manID$variable <- "CH4"


table_results <- rbind(CO2_flux_res_auto, CO2_flux_res_manID, CH4_flux_res_auto, CH4_flux_res_manID)
table_results$username <- username

table_results_ebull <- rbind(CH4_res_meth1, CH4_res_meth2)
table_results_ebull$username <- username


mytimestamp <- Sys.time() 
table_results$timestamp_processing <- mytimestamp # to keep a track of when the processing was done and ease data analysis in case an incubation is processed multiple times
table_results_ebull$timestamp_processing <- mytimestamp 

# saving fluxes estimates
setwd(results_path)

append_if_exists <- function(filename, data){
  if(file.exists(filename)){
    data <- rbind(read.csv(filename), data)
    write.csv(x = data, file = filename, 
              row.names = F)
  } else {
    write.csv(x = data, file = filename, 
              row.names = F)
  }
  return(data)
}


filename <- paste0("BLIND_vs_EXPERT_co2_ch4_fluxes_",Sys.Date(),".csv")
table_results_all <- append_if_exists(filename, data = table_results)


filename <- paste0("BLIND_vs_EXPERT_ch4_ebullition_",Sys.Date(),".csv")
table_results_ebull_all <- append_if_exists(filename, data = table_results_ebull)


filename <- paste0("BLIND_vs_EXPERT_table_draw_",Sys.Date(),".csv")
table_draw_all <- append_if_exists(filename, data = table_draw)



#----- some plots for the incubations processed in the current session -----



p_auto_vs_manual <- ggplot(data = table_results)+
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
ggsave(filename = paste0("BLIND_vs_EXPERT_co2_ch4_fluxes_",Sys.Date(),".jpeg"), 
       plot = p_auto_vs_manual, path = plots_path, width = 10, height = 10, units = 'in', dpi = 300)


p_auto_vs_manual_ch4_ebullition <- ggplot(data = table_results_ebull)+
  geom_abline(slope = 0,intercept = 0, color = 'lightgrey')+
  geom_segment(data = data.frame(UniqueID = CH4_res_meth1$UniqueID,
                                 meth1 = CH4_res_meth1$ebullition,
                                 meth2 = CH4_res_meth2$ebullition), aes(x=UniqueID, xend=UniqueID, y = meth1, yend = meth2), linewidth=1, alpha = 0.5)+
  geom_point(aes(reorder(UniqueID, -ebullition, FUN=mean), ebullition, colour = flux_method), size=4, alpha = 0.5)+
  
  ylab("ebullition component [nmol/m2/s]")+
  theme_bw()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_colour_viridis_d(begin = 0.1, end = 0.9, option = "F")+coord_flip()

p_auto_vs_manual_ch4_ebullition
ggsave(filename = paste0("BLIND_vs_EXPERT_ch4_ebullition_",Sys.Date(),".jpeg"), 
       plot = p_auto_vs_manual_ch4_ebullition, path = plots_path, width = 10, height = 10, units = 'in', dpi = 300)




table_results_sprd_CO2  <- table_results[table_results$variable=="CO2",c("UniqueID","flux_method","best.flux","timestamp_processing")] %>%
  pivot_wider(names_from = flux_method, values_from = c(best.flux))

median((table_results_sprd_CO2$Blind-table_results_sprd_CO2$Expert)/table_results_sprd_CO2$Expert*100, na.rm = T)

ggplot(data = table_results_sprd_CO2)+
  geom_abline(slope = 0,intercept = 0, color = 'lightgrey')+
  # geom_segment(data = data.frame(UniqueID = CO2_flux_res_auto$UniqueID,
  #                                meth1 = CO2_flux_res_auto$best.flux,
  #                                meth2 = CO2_flux_res_manID$best.flux), aes(x=UniqueID, xend=UniqueID, y = meth1, yend = meth2), linewidth=1, alpha = 0.5)+
  geom_point(aes(UniqueID, (Blind-Expert)/Expert*100), size=4, alpha = 0.5)+
  ylab("CO2 flux relative difference [%]")+
  theme_bw()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(table_results_sprd_CO2, aes((Blind-Expert)/Expert*100))+geom_density()


table_results_sprd_CH4  <- table_results[table_results$variable=="CH4",c("UniqueID","flux_method","best.flux","timestamp_processing")] %>%
  pivot_wider(names_from = flux_method, values_from = c(best.flux))

median((table_results_sprd_CH4$Blind-table_results_sprd_CH4$Expert)/table_results_sprd_CH4$Expert*100, na.rm = T)

ggplot(data = table_results_sprd_CH4)+
  geom_abline(slope = 0,intercept = 0, color = 'lightgrey')+
  # geom_segment(data = data.frame(UniqueID = CH4_flux_res_auto$UniqueID,
  #                                meth1 = CH4_flux_res_auto$best.flux,
  #                                meth2 = CH4_flux_res_manID$best.flux), aes(x=UniqueID, xend=UniqueID, y = meth1, yend = meth2), linewidth=1, alpha = 0.5)+
  geom_point(aes(UniqueID, (Blind-Expert)/Expert*100), size=4, alpha = 0.5)+
  ylab("CH4 flux relative difference [%]")+
  theme_bw()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))




#----- some plots for the all incubations processed so far -----


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

ggsave(filename = "BLIND_vs_EXPERT_co2_ch4_all.jpeg", plot = p_auto_vs_manual, path = plots_path, width = 10, height = 10, units = 'in', dpi = 300)




table_results_sprd <- table_results_all[,c("variable","UniqueID","flux_method","best.flux","timestamp_processing")] %>%
  pivot_wider(names_from = flux_method, values_from = c(best.flux))

ggplot(data = table_results_sprd)+
  # geom_abline(slope = 0,intercept = 0, color = 'lightgrey')+
  # geom_segment(data = data.frame(UniqueID = CO2_flux_res_auto$UniqueID,
  #                                meth1 = CO2_flux_res_auto$best.flux,
  #                                meth2 = CO2_flux_res_manID$best.flux), aes(x=UniqueID, xend=UniqueID, y = meth1, yend = meth2), linewidth=1, alpha = 0.5)+
  geom_violin(aes(variable,abs((Blind-Expert)/Expert)*100), alpha = 0.5)+
  geom_jitter(aes(variable,abs((Blind-Expert)/Expert)*100), alpha = 0.5, width=0.2)+
  ylab("CO2 flux relative difference [%]")+
  theme_bw()+
  scale_y_log10()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



list_f <- list.files(path = results_path, pattern = "BLIND_vs_EXPERT_ch4_ebullition_", full.names = T)
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




p_auto_vs_manual_ch4_ebullition <- ggplot(data = table_results_all)+
  geom_abline(slope = 0,intercept = 0, color = 'lightgrey')+
  geom_segment(data = data.frame(UniqueID = CH4_res_meth1$UniqueID,
                                 meth1 = CH4_res_meth1$ebullition,
                                 meth2 = CH4_res_meth2$ebullition), aes(x=UniqueID, xend=UniqueID, y = meth1, yend = meth2), linewidth=1, alpha = 0.5)+
  geom_point(aes(reorder(UniqueID, -ebullition, FUN=mean), ebullition, colour = flux_method), size=4, alpha = 0.5)+
  
  ylab("ebullition component [nmol/m2/s]")+
  theme_bw()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_colour_viridis_d(begin = 0.1, end = 0.9, option = "F")+coord_flip()

p_auto_vs_manual_ch4_ebullition
ggsave(filename = "BLIND_vs_EXPERT_ch4_ebullition_all.jpeg", plot = p_auto_vs_manual_ch4_ebullition, path = plots_path, width = 10, height = 6, units = 'in', dpi = 300)










