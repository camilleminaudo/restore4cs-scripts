

# ---
# Authors: MIGUEL CABRERA
# Project: "RESTORE4Cs"
# date: "Feb 2025"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
# This script computes fluxes for all core incubations analyzed in the UB with N2O and CO2+CH4 Licors.
# Table with GHG concentrations (ppm) calculated by Script Summary_CORES.R (in https://github.com/MCabreraBrufau/Measuring_discrete_GHG_samples_with_Licor.git)

#Fluxes are computed in ppm as the difference of concentration for each core at Tf minus the average concentration found at t0 in the initial cores analyzed (usually C1, C3 and C5). Significance of flux is added as T-test p-value (assuming unequal variances)


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
# library(goFlux)
# require(dplyr)
# require(purrr)
require(msm)
# require(data.table)
require(tools)
# require(pbapply)

repo_root <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
files.sources = list.files(path = paste0(repo_root,"/functions"), full.names = T)
for (f in files.sources){source(f)}


# ---- Directories ----
dropbox_root <- "C:/Users/Miguel/Dropbox/RESTORE4Cs - Fieldwork/Data" # You have to make sure this is pointing to the write folder on your local machine
results_path<- paste0(dropbox_root,"/Cores/UB_concentrations/")
fieldsheetpath <- paste0(dropbox_root,"/Sediment/Fieldsheets")



#_________________####


#FIELDSHEET PREP --------

#ADAPT for core fieldsheet retrieval + incubation details UVEG (duration, volume air/water/sed, incubation conditions, obs,... )


#OLD, for chambers: 
# {
# #Data pre-processing and harmonization: 
# # ---- List GHG chamber fieldsheets in Dropbox and read them ---
# # list filenames
# myfieldsheets_list <- list.files(fieldsheetpath, pattern = "Fieldsheet-GHG.xlsx", all.files = T, full.names = T, recursive = T)
# i <- grep(pattern = sampling, x = myfieldsheets_list) # selecting the files corresponding to the selected sampling campaign
# myfieldsheets_list <- myfieldsheets_list[i]
# # Read all fieldsheets and put them in a single dataframe
# fieldsheet <- read_GHG_fieldsheets(myfieldsheets_list)
# 
# 
# #Create UniqueID_notime column to join
# fieldsheet$UniqueID_notime<- fieldsheet %>% 
#   mutate(UniqueID_notime=tolower(paste(subsite, plot_id,substr(strata,1,1) ,substr(transparent_dark,1,1), sep = "-"))) %>% pull(UniqueID_notime)
#   
#   
# #Check: UniqueID_notime is unique?
# length(unique(fieldsheet$UniqueID_notime))==dim(fieldsheet)[1]
#   
# # recalculating start and stop in proper formats
# fieldsheet$timestamp_start <- as.POSIXct(fieldsheet$unix_start, tz = "UTC", origin = "1970-01-01")
# fieldsheet$timestamp_stop <- as.POSIXct(fieldsheet$unix_stop, tz = "UTC", origin = "1970-01-01")
# 
# fieldsheet$start_time <- strftime(fieldsheet$timestamp_start, format="%H:%M:%S", tz = 'utc')
# fieldsheet$end_time <- strftime(fieldsheet$timestamp_stop, format="%H:%M:%S", tz = 'utc')
# 
# fieldsheet <- fieldsheet[order(fieldsheet$subsite),]
# }


#_________________####
#FLUX CALCULATION####
# ----- 1. Import & summary -----

injections_ppm<- read.csv(paste0(results_path,"/N2O_CO2_CH4_ppm_all_injections.csv"))
#IMPORTANT: t0 concentrations are very different for many individual cores for CH4 and CO2 (it is not appropriate to average them to calculate the deltappm). 

##---TO DECIDE----
#Initial exploration reveals CO2 and CH4 t0 concentrations very different depending on individual core. 
#CHECK for which subsites it is appropriate to average all initial cores 1,3,5. For those that it is not appropriate, decide if we go back and analyse the t0 of core-replicate 2,4,6 for the available samplings (all S2, all S3, S4-CA, in total 234 additional samples) or if we use the initial concentration from the notes of Valencianos. 

#DECISSION: We trust more the initial t0 measured on-site with gas analyzer than the exetainer t0 values, for CO2 and CH4 we will use the t0 concentrations annotated in the lab by the valencianos. 

#For N2O, we do need initial measurement with exetainers. We will analyze only samples for which the average CV of the measured t0 exetainers very high: S2-CA-A2 CV=0.63, one sample with very very high N2O. The following highest CV is 0.06, which is not that high and arises from general variability, not from one sample being too high. 


#For the moment proceed with average of 3 initial cores analysed for the deltaGHG of cores for which we did not analyse the initial concentration. Use core-specific t0 and tf for the cores with available data.

#Calculate average sd n of initial cores and put in different columns side by side of tf summaries.
#For each core_ID, create 2 sets of columns: t0 and tf with _avg, _sd and _n 
initialfinal_ppm_summary<- injections_ppm %>% 
  select(sample, gas, ppm) %>% 
  separate(sample, into = c("season","site","subsite","coretime"),remove = F) %>% 
  mutate(corenum=parse_number(coretime),
         time=str_extract(coretime, "i|f"),
         sampling=paste(season, site,subsite, sep="-"),
         coreID=paste0(sampling, "-C",corenum),
         t0_ppm=case_when(time=="i"~ppm,
                          TRUE~NA_real_)) %>% 
  #Create allt0 _avg, _sd, _n per sampling and gas
  group_by(sampling,gas) %>% 
  mutate(allt0_avg=mean(t0_ppm, na.rm=T),
         allt0_sd=sd(t0_ppm, na.rm = T),
         allt0_n=sum(!is.na(t0_ppm))) %>% 
  #Create t0 avg, sd and n (if available)
  group_by(gas, sampling, coreID) %>% 
  mutate(t0_avg=mean(t0_ppm, na.rm=T),
         t0_sd=sd(t0_ppm, na.rm=T),
         t0_n=sum(!is.na(t0_ppm))) %>% 
  #Remove t0 rows
  filter(time=="f") %>% 
  #Create tf avg, sd and n
    summarise(allt0_avg=mean(allt0_avg),allt0_sd=mean(allt0_sd),allt0_n=mean(allt0_n), #keep allt0 columns
              t0_avg=mean(t0_avg, na.rm=T), t0_sd=mean(t0_sd, na.rm=T),t0_n=mean(t0_n, na.rm=T),#keep t0 columns
            tf_avg=mean(ppm, na.rm=T),
            tf_sd=sd(ppm, na.rm=T),
            tf_n=sum(!is.na(ppm))) %>% 
  #Create source_t0 column and populate t0 columns with data from allt0 columns for cores for which we did not analyse the actual t0 values 
  mutate(source_t0=case_when(is.na(t0_avg)~"averaget0s",
                             TRUE~"measured"),
         t0_avg=case_when(source_t0=="averaget0s"~allt0_avg,
                          source_t0=="measured"~t0_avg),
         
         t0_sd=case_when(source_t0=="averaget0s"~allt0_sd,
                          source_t0=="measured"~t0_sd),
         
         t0_n=case_when(source_t0=="averaget0s"~allt0_n,
                          source_t0=="measured"~t0_n)) %>% 
  #Remove allt0 columns
  select(-c(allt0_avg,allt0_sd,allt0_n))
         
         
         
  
# ---- 2. DeltaGHG & pvalue  -------
#TO ADAPT: create t-test for valencianos t0 CO2 and CH4 (1 single value)

#Calculate deltaGHG in ppm
#Add pvalue for significance of flux

#Calculate deltaGHG in ppm and the significance of the deltaghg (i.e. test if the concentration at tf is significantly different from that used as t0)
deltaGHG_ppm_pval<- initialfinal_ppm_summary %>% 
  mutate(
    #Calculate deltaghg_ppm
    deltaghg_ppm=tf_avg-t0_avg,
    # Calculate standard error for the difference
    se_diff = sqrt((tf_sd^2 / tf_n) + (t0_sd^2 / t0_n)),
    
    # Calculate the t-statistic for the difference
    t_stat = (tf_avg - t0_avg) / se_diff,
    
    # Calculate the degrees of freedom using the formula for unequal variances
    df_ = (((tf_sd^2 / tf_n) + (t0_sd^2 / t0_n))^2) /
      (((tf_sd^2 / tf_n)^2 / (tf_n - 1)) + ((t0_sd^2 / t0_n)^2 / (t0_n - 1))),
    
    # Calculate the p-value using the t-distribution
    p_value = 2 * pt(-abs(t_stat), df_)
  ) %>% 
  #remove unnecesary variables
  select(-c(df_,t_stat))
  

##---_Save deltaGHGppm-----
#Save deltaGHG (in ppm units), we will calculate the flux in propper units after exchange with carlos and antonio. 

write.csv(x = deltaGHG_ppm_pval, file = paste0(results_path, "deltaGHG_ppm_pval_usingt0average.csv"),row.names = F)




#N2O flux check: 
deltaGHG_ppm_pval %>% 
  filter(gas=="n2o") %>% 
  separate(sampling, into = c("season","site","subsite")) %>% 
  filter(deltaghg_ppm<1) %>% 
  ggplot(aes(y=deltaghg_ppm,  fill=p_value<0.05))+
  geom_histogram(position = "identity", alpha=0.3,bins = 50)

n2o_flux<- deltaGHG_ppm_pval %>% filter(gas=="n2o") %>% filter(!is.na(deltaghg_ppm))

message(paste("N2O Fluxes: Out of ", dim(n2o_flux)[1]), " N2O fluxes, we have ", sum(n2o_flux$p_value<0.05), " significant fluxes (",  round(sum(n2o_flux$p_value<0.05)/dim(n2o_flux)[1]*100,2), "%), at the p<0.05 level")




#CO2 flux check: 
deltaGHG_ppm_pval %>% 
  filter(gas=="co2") %>% 
  separate(sampling, into = c("season","site","subsite")) %>% 
  # filter(deltaghg_ppm<1) %>% 
  ggplot(aes(y=deltaghg_ppm,  fill=p_value<0.05))+
  geom_histogram(position = "identity", alpha=0.3,bins = 50)

co2_flux<- deltaGHG_ppm_pval %>% filter(gas=="co2") %>% filter(!is.na(deltaghg_ppm))

message(paste("CO2 Fluxes: Out of ", dim(co2_flux)[1]), " CO2 fluxes, we have ", sum(co2_flux$p_value<0.05), " significant fluxes (",  round(sum(co2_flux$p_value<0.05)/dim(co2_flux)[1]*100,2), "%), at the p<0.05 level")



#CH4 flux check: 
deltaGHG_ppm_pval %>% 
  filter(gas=="ch4") %>% 
  separate(sampling, into = c("season","site","subsite")) %>% 
  # filter(deltaghg_ppm<1) %>% 
  ggplot(aes(y=deltaghg_ppm,  fill=p_value<0.05))+
  geom_histogram(position = "identity", alpha=0.3,bins = 50)

ch4_flux<- deltaGHG_ppm_pval %>% filter(gas=="ch4") %>% filter(!is.na(deltaghg_ppm))

message(paste("CH4 Fluxes: Out of ", dim(ch4_flux)[1]), " CH4 fluxes, we have ", sum(ch4_flux$p_value<0.05), " significant fluxes (",  round(sum(ch4_flux$p_value<0.05)/dim(ch4_flux)[1]*100,2), "%), at the p<0.05 level")



##---Decide N2O to analyze----
#N2O exploration: 
#Check CV of allt0 
deltaGHG_ppm_pval %>% 
  ungroup() %>% 
  filter(gas=="n2o") %>% 
  filter(source_t0=="averaget0s") %>% 
  mutate(t0_cv=t0_sd/t0_avg) %>% 
  ggplot(aes(x=t0_cv))+
  geom_histogram(bins=50)+
  geom_vline(xintercept=0.05)+
  facet_wrap(.~gas, scales="free")

#Cores without measured N2O at t0 and for which subsite the measured t0 cores have a grouped CV larger than 5%
deltaGHG_ppm_pval %>% 
        ungroup() %>% 
  filter(gas=="n2o") %>% 
  filter(source_t0=="averaget0s") %>% 
  mutate(t0_cv=t0_sd/t0_avg) %>% 
  filter(t0_cv>0.05) 


#N2O general exploration of initial values
testinitial %>% 
  filter(gas=="n2o") %>% 
  ggplot(aes(x=corenum, y=ppm, col=season))+
  geom_point()+
  facet_grid(site~subsite,scales="free")

injections_ppm %>% 
  filter(gas=="n2o") %>% 
  separate(sample, into = c("season","site","subsite","coretime")) %>% 
  mutate(corenum=parse_number(coretime),
         time=str_extract(coretime, pattern = "i|f")) %>% 
  ggplot(aes(x=corenum, y=ppm, col=time))+
  geom_point()+
  facet_grid(site~subsite,scales="free")



head(deltaGHG_ppm_pval %>% 
       filter(gas=="n2o") %>% 
       arrange(desc(tf_avg)))



# ---- 3. Convert to mol-------

#Once we have the necessary details (i.e. P, Vol, T, for cores and incubation time), we can convert the deltaGHG from  ppm to mols per time and per surface 

#Leftover from N2o chamber calcualtion to be used as a guide.
#Overall approach:
#1st calculate deltaN2O in nmols (using chamber volume and Ideal Gas Law)
#2nd calculate change in nmols per second 
#3rd calculate areal flux in nmols s^-1 m^-2  

#Using Ideal Gas Law: PV = nRT, where:
#P: Pressure (atm)
#V: Volume (L)
#n: number of mols
#T: Temperature (ºK)
#R: Gas Constant,  R = 0.08206 (L·atm·K-1·mol-1)
R_constant<- 0.08206

#Re-arranging: n = PV/RT

#to calculate mols of gas in the chamber at t0 and tf, we will 1st calculate the Volume of n2o (Vn2o, in L) using the ppm and total volume of chamber. Then apply the ideal gas law to get the mols at the given Tcham and Pressure.

#the goflux package uses a slightly different function, with P in KPa, and incorporating a water-vapor correction. Unclear why, we will stick to the more simple Ideal Gas Law. 

# 4.1. Import N2O_ppm and transform units.

#Import atm and tf N2O concentration in ppm, identified with "UniqueID_notime", atm_N2Oppm, tf_N2Oppm 
n2o<- read.csv(paste0(results_path,"/S4_restore4cs_N2Oppm_atmchambers.csv"))


#FLUX CAlCULATION: #Auxfile has Pcham in kPa, Area in cm2 and Tcham in ºcelsius, Vtot is already in L
n2o_flux<- n2o %>% 
  #Join with auxfile
  merge.data.frame(auxfile, by="UniqueID_notime", all.x = T) %>% 
  #Transform Units:
  mutate(Vtot, #Vtot is already in L
         Pcham_atm= 1, # set Pcham to 1 atm
         Tcham_k= Tcham+273.15, #Transform Tcham to Kelvin
         duration,#Duration of incubation (in seconds) calculated from fieldsheet times
         Area_m2=Area*1e-4) %>% # divide by 10,000 to transform cm2 to m2  
  #calculate Vn2o, at t0 and tf
  mutate(V_n2o_t0 = (atm_avgN2Oppm/1e6)*Vtot, 
         V_n2o_tf = (tf_avgN2Oppm/1e6)*Vtot) %>% 
  #calculate nmol_n2o at t0 and tf (apply mol = PV/RT and then multiply by 1e9 to get nmol)
  mutate(nmol_n2o_t0 = 1e9* ( (Pcham_atm*V_n2o_t0)/(R_constant*Tcham_k) ),
         nmol_n2o_tf = 1e9* ( (Pcham_atm*V_n2o_tf)/(R_constant*Tcham_k) )) %>% 
  #calculate delta_nmol_n2o (as nmol difference tf-t0)
  mutate(delta_nmol_n2o = nmol_n2o_tf-nmol_n2o_t0) %>% 
  #calculate absolute flux: nmol_n2o_per_second
  mutate(nmol_n2o_per_second=delta_nmol_n2o/duration) %>% 
  #calculate areal flux: nmol_n2o_per_second_m2
  mutate(nmol_n2o_per_second_m2=nmol_n2o_per_second/Area_m2) %>% 


#CHUNK to calculate the significance of the fluxes (i.e. if the concentration at tf is significantly different from atmospheric.)
mutate(
  # Calculate standard deviations from cv (cv = SD / Mean, so SD = Mean * cv)
  tf_sd = tf_avgN2Oppm * tf_cv,
  atm_sd = atm_avgN2Oppm * atm_cv,
  
  # Calculate standard error for the difference
  se_diff = sqrt((tf_sd^2 / tf_n) + (atm_sd^2 / atm_n)),
  
  # Calculate the t-statistic for the difference
  t_stat = (tf_avgN2Oppm - atm_avgN2Oppm) / se_diff,
  
  # Calculate the degrees of freedom using the formula for unequal variances
  df_ = (((tf_sd^2 / tf_n) + (atm_sd^2 / atm_n))^2) /
    (((tf_sd^2 / tf_n)^2 / (tf_n - 1)) + ((atm_sd^2 / atm_n)^2 / (atm_n - 1))),
  
  # Calculate the p-value using the t-distribution
  p_value = 2 * pt(-abs(t_stat), df_)
)





# ---- 3. N2O flux filtering-------


#Add significance and clean up table (remove unncecessary variables.)
#FILTER also fluxes of which we cannot be sure (plotincubation==2, unclear or incomplete ventilation between transparent and dark incubation)
n2o_flux_simple<- n2o_flux %>% 
  separate(UniqueID_notime, into = c("sampling","pilotsite","subsite","plot","d1","d2"),remove = F) %>% 
  mutate(siteID=paste(pilotsite,subsite, sep = "-")) %>% 
  filter(plotincubation==1) %>% 
  select(UniqueID_notime, sampling, pilotsite, subsite, siteID, start.time, duration, water_depth, strata, chamberType,lightCondition, 
         nmol_n2o_per_second_m2, p_value)
  

message(paste("All Fluxes: Out of ", dim(n2o_flux)[1]), " N2O fluxes, we have ", sum(n2o_flux$p_value<0.05), " significant fluxes (",  round(sum(n2o_flux$p_value<0.05)/dim(n2o_flux)[1]*100,2), "%), at the p<0.05 level")


message(paste("Only reliable fluxes (clear ventilation): Out of ", dim(n2o_flux_simple)[1]), " N2O fluxes, we have ", sum(n2o_flux$p_value<0.05), " significant fluxes (",  round(sum(n2o_flux$p_value<0.05)/dim(n2o_flux_simple)[1]*100,2), "%), at the p<0.05 level")


#__________________####



#SAVE RESULTS####
write.csv(n2o_flux, file = paste0(datapath,"/S4_restore4cs_N2O_EXTRADATA_n2ofluxes.csv"),row.names = F)

write.csv(n2o_flux_simple, file=paste0(datapath,"/S4_restore4cs_N2O_arealflux_nmol_s-1_m-2.csv"),row.names = F)



#Miscelanea: 

ggplot(n2o_flux_simple, aes(y=nmol_n2o_per_second_m2,  fill=p_value<0.05))+
  geom_histogram(position = "identity", alpha=0.3)

