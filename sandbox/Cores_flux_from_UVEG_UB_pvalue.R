

# ---
# Authors: MIGUEL CABRERA
# Project: "RESTORE4Cs"
# date: "April 2025"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
# This script computes per-area fluxes for all core incubations. 
#N2O initial concentrations are derived from samples analyzed in the UB. T0 N2O concentrations are either measured (initial exetainer) or averaged from same-subsite cores (after checking the variability of initial exetainers, minimum analyzed: 3 out of 6, when CV of the 3 initial exetainers exceeds 6%, the remaining 3 initial exetainers for that subsite are also analyzed).
#N2O final concentrations are measured from final exetainers in the UB with zero-air baseline. 

#all CO2 and CH4 initial concentrations are obtained from the in-situ measurements with Licor (after estabilization of measure inside the open core, the concentration is written down and the core is closed to start the incubation).

#CO2 and CH4 final concentrations are derived from UB measurements with low baseline for Seasons 2-4. For Season 1, final concentrations of CO2 and CH4 are derived from in-situ injections (with ambient air baseline)

#Incubation details (duration, headspace, water) where measured in-situ and written down, for all incubated cores. Incubation conditions are not homogeneous, depending on case-pilot site a temperature controlled chamber is used or cores are incubated outside the lab at ambient temperature and light conditions.  


# Concentrations derived from injection of samples using the open-loop method were analyzed with the scripts available at: https://github.com/MCabreraBrufau/Measuring_discrete_GHG_samples_with_Licor.git


#ADD calculation description------



rm(list = ls()) # clear workspace
cat("/014") # clear console

# ---- packages ----
library(tidyverse)
library(readxl)
# library(lubridate)
# library(zoo)
library(ggplot2)
# library(grid)
# library(egg)
# require(msm)

# repo_root <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
# files.sources = list.files(path = paste0(repo_root,"/functions"), full.names = T)
# for (f in files.sources){source(f)}


# ---- Directories ----
dropbox_root <- "C:/Users/Miguel/Dropbox/RESTORE4Cs - Fieldwork/Data" # You have to make sure this is pointing to the write folder on your local machine
results_UB_path<- paste0(dropbox_root,"/Cores/UB_concentrations/")
results_uveg_path<- paste0(dropbox_root,"/Cores/UVEG_concentrations/")
# fieldsheetpath <- paste0(dropbox_root,"/Sediment/Fieldsheets")
core_path<- paste0(dropbox_root,"/Cores/")


#_________________####


#incubation-data --------

incub<- read_xlsx(path = paste0(core_path,"Incubation_details_all_coresS1S2S3S4_para rellenar.xlsx")) %>% 
  select(core_id, headspace_cm, water_cm, approximatesediment_cm, coreincubation_hours,comment, t0_co2_ppm, t0_ch4_ppb)


#core-fieldsheet-----
#Load core-fieldsheet compilation (created with Compile_Core_sampling_details.R)

field<- read.csv(file = paste0(core_path, "Core_fieldsheet_compilation_2025-03-17.csv"))


#List subsites with missing sampling info in fieldsheets: 
incub[!incub$core_id%in%field$sample_id,] %>%
  select(core_id) %>% 
  mutate(subsite=substr(core_id, 1,8)) %>%
  group_by(subsite) %>% 
  summarise(cores_missing=n())
#4 subsites without any info for core sampling (core_id not in fieldsheet)




#method CV -----
#Retrieve open-loop method CV from one-point calibration based on repeated injections of standards. 
open_loop_cal<- read.csv(paste0(results_UB_path,"One-point_calibration_factor.csv")) %>% 
  mutate(cv_method=sd_factor/factor)

co2_cv_method<- open_loop_cal[open_loop_cal$gas=="co2",]$cv_method
ch4_cv_method<- open_loop_cal[open_loop_cal$gas=="ch4",]$cv_method
n2o_cv_method<- open_loop_cal[open_loop_cal$gas=="n2o",]$cv_method


#UB data----
#Load exetainer ppm concentrations: 
samplesUB_ppm<- read.csv(file = paste0(results_UB_path,"N2O_CO2_CH4_ppm_exetainer_avg_sd_n.csv"))

#Split and format

#Extract all tf concentrations from UB
##gas, core_id, coretime, avg_ppm, sd_ppm, n_injections:
tfub<- samplesUB_ppm %>% 
  filter(grepl("f", sample)) %>% # get only tf
  select(-dayofanalysis) %>% 
  #obtain core_id and coretime
  separate(sample, into = c("season","site","subsite","coretime"),remove = F) %>% 
  mutate(corenum=parse_number(coretime),
         coretime="final",
         sampling=paste(season, site,subsite, sep="-"),
         core_id=paste0(sampling, "-C",corenum)) %>% 
  #calculate SE for each core
  mutate(se_sample=sd_ppm/sqrt(n_injections),
         estimate_origin="Measured open-loop") %>% 
  #Propagate method precision onto SE
  mutate(se_method=case_when(gas=="co2"~(avg_ppm*co2_cv_method)/sqrt(n_injections),
                             gas=="ch4"~(avg_ppm*ch4_cv_method)/sqrt(n_injections),
                             gas=="n2o"~(avg_ppm*n2o_cv_method)/sqrt(n_injections)),
         se_ppm=sqrt(se_sample^2 + se_method^2)) %>% 
  rename(tf_ppm=avg_ppm, tf_se_ppm=se_ppm, tf_estimate_origin=estimate_origin) %>% 
  select(gas, core_id, tf_ppm, tf_se_ppm, tf_estimate_origin)


#Extract measured t0 N2O from UB and fill gaps
#N2O t0 will be extracted and missing values calculated from the concentration of each measured core (from avg of measured cores within each subsite sampling and propagating the error). Final format: core_id, avg_ppm, se_ppm, n_injections, estimate_origin (peak-integration or averaged)

t0ub_n2o<- samplesUB_ppm %>% 
  filter(gas=="n2o") %>% #get only n2o 
  select(-dayofanalysis) %>% 
  #obtain core_id and coretime
  separate(sample, into = c("season","site","subsite","coretime"),remove = F) %>% 
  mutate(corenum=parse_number(coretime),
         coretime=str_extract(coretime, "i|f"),
         sampling=paste(season, site,subsite, sep="-"),
         core_id=paste0(sampling, "-C",corenum)) %>% 
  select(gas, core_id, coretime, avg_ppm, sd_ppm, n_injections) %>% 
  #separate initial measures for each core_id (leaving NA for final samples)
  mutate(t0_avg_ppm=case_when(coretime=="i"~avg_ppm, TRUE~NA_real_),
         t0_sd_ppm=case_when(coretime=="i"~sd_ppm, TRUE~NA_real_),
         t0_n_injections=case_when(coretime=="i"~n_injections)) %>% 
  #Select only gas, core_id, and t0 values, agregate to leave core_id without measured t0
  select(gas, core_id, t0_avg_ppm,t0_sd_ppm,t0_n_injections) %>% 
  group_by(gas,core_id) %>% 
  summarise(avg_ppm=mean(t0_avg_ppm, na.rm=T),
            sd_ppm=mean(t0_sd_ppm, na.rm=T),
            n_injections=mean(t0_n_injections, na.rm=T)
  ) %>% 
  #Calculate SE for each exetainer
  mutate(se_sample=sd_ppm/sqrt(n_injections),
         sampling=substr(core_id, 1, 8)) %>% 
  #Propagate method precission onto SE
  mutate(se_method=case_when(gas=="co2"~(avg_ppm*co2_cv_method)/sqrt(n_injections),
                             gas=="ch4"~(avg_ppm*ch4_cv_method)/sqrt(n_injections),
                             gas=="n2o"~(avg_ppm*n2o_cv_method)/sqrt(n_injections)),
         se_ppm=sqrt(se_sample^2 + se_method^2)) %>% 
  #Calculate subsite-level average
  ungroup() %>% 
  group_by(sampling) %>%
  #Calculate for each subsite: subsite-level average, within-sample variance and between-sample variance, and finally, combine variances to get subsite-level SE. 
  mutate(subsite_avg_ppm=mean(avg_ppm,na.rm=T),
         k=sum(!is.na(avg_ppm)),
         var_between=sum((avg_ppm - subsite_avg_ppm)^2,na.rm = T) / (k - 1),
         var_within = sum(se_ppm^2,na.rm=T) / k^2, 
         total_variance = var_between / k + var_within / k,
         subsite_se_ppm = sqrt(total_variance),
         coretime="initial"# re-introduce coretime for all core_id
         ) %>% 
  ungroup() %>% 
  select(gas, core_id, coretime, avg_ppm, sd_ppm, n_injections, se_ppm, subsite_avg_ppm,subsite_se_ppm) %>% 
  #Fill gaps for is.na(n_injections) and indicate source of t0 estimate
mutate(estimate_origin= case_when(is.na(n_injections)~"Subsite average, propagated SE",
                           TRUE~"Measured open-loop"),
       avg_ppm= case_when(is.na(n_injections)~subsite_avg_ppm, 
                          TRUE~avg_ppm),
       se_ppm= case_when(is.na(n_injections)~subsite_se_ppm,
                         TRUE~se_ppm),
       n_injections=case_when(!is.na(n_injections)~n_injections,
                              TRUE~0)) %>% 
  rename(t0_ppm=avg_ppm, t0_se_ppm=se_ppm, t0_estimate_origin=estimate_origin) %>% 
  select(gas, core_id, t0_ppm, t0_se_ppm, t0_estimate_origin)


# UVEG data----


#Load S1 final concentrations measured in-situ with UVEG licor

samplesUVEG_ppm <- read.csv(file = paste0(results_uveg_path,"CO2_CH4_ppm_core_avg_sd_n.csv"))

#SE will be calculate for each sample, for estimates from average baseline, we will use the n of baseline. For injections with just one injection, we will impute 0 as the SE of sample. 
#Open-loop method precision (CV) will be propagated onto the final SE
#Instrument precission will be used in the cases where the average of the whole remark was used (due to peak-detection failure)


tfuveg<- samplesUVEG_ppm %>% 
  #Calculate SE for each type of estimate
  mutate(se_sample=case_when(estimate=="peak-integration"~sd_ppm/sqrt(n_injections),
                          TRUE~sd_ppm/sqrt(n_base)),
         estimate_origin=case_when(estimate=="peak-integration"~"Measured open-loop",
                                   TRUE~"Average of remark"),
         n_injections=case_when(is.na(n_injections)~0,
                                TRUE~n_injections),
         se_sample=case_when(is.na(se_sample)~0, 
                                TRUE~se_sample)) %>% 
  #Propagate method precission onto SE for samples injected, for average of remarks use instrument precission (3.5 ppm for CO2 and 0.6ppb for CH4)
  mutate(se_method=case_when(estimate_origin=="Measured open-loop"&gas=="co2"~(avg_ppm*co2_cv_method)/sqrt(n_injections),
                             estimate_origin=="Measured open-loop"&gas=="ch4"~(avg_ppm*ch4_cv_method)/sqrt(n_injections),
                             estimate_origin=="Average of remark"&gas=="co2"~3.5/sqrt(n_base),
                             estimate_origin=="Average of remark"&gas=="ch4"~(0.6/1000)/sqrt(n_base)),
         se_ppm=sqrt(se_sample^2 + se_method^2)) %>% 
  separate(sample, into = c("season","site","subsite","coretime"),remove = F) %>% 
  mutate(corenum=parse_number(coretime),
         sampling=paste(season, site,subsite, sep="-"),
         core_id=paste0(sampling, "-C",corenum)) %>% 
  rename(tf_ppm=avg_ppm, tf_se_ppm=se_ppm, tf_estimate_origin=estimate_origin) %>% 
  select(gas, core_id, tf_ppm, tf_se_ppm, tf_estimate_origin)




mixesUVEG_ppm <- read.csv(file = paste0(results_uveg_path,"CH4_ppm_coremix_avg_sd_cv_n.csv"))


tfmix<- mixesUVEG_ppm %>% 
  #obtain core_id and coretime
  separate(sample, into = c("season","site","subsite","coretime"),remove = F) %>% 
  mutate(corenum=parse_number(coretime),
         coretime="final",
         sampling=paste(season, site,subsite, sep="-"),
         core_id=paste0(sampling, "-C",corenum)) %>% 
  #calculate SE for each core
  mutate(se_sample=sd_ppm/sqrt(n_injections),
         estimate_origin="Measured open-loop") %>% 
  #Impute SE=0 in case of only 1 valid injection
  mutate(se_sample=if_else(is.na(se_sample),true = 0,false = se_sample)) %>% 
  #Propagate method precision onto SE
  mutate(se_method=case_when(gas=="co2"~(avg_ppm*co2_cv_method)/sqrt(n_injections),
                             gas=="ch4"~(avg_ppm*ch4_cv_method)/sqrt(n_injections),
                             gas=="n2o"~(avg_ppm*n2o_cv_method)/sqrt(n_injections)),
         se_ppm=sqrt(se_sample^2 + se_method^2)) %>% 
  rename(tfmix_ppm=avg_ppm, tfmix_se_ppm=se_ppm, tfmix_estimate_origin=estimate_origin) %>% 
  select(gas, core_id, tfmix_ppm, tfmix_se_ppm, tfmix_estimate_origin)





#Calculate SE of t0 measured in situ using instrument precission for each gas and assuming 30s average for each t0 reading. 
t0uveg_co2_ch4<- incub %>% 
  select(core_id, t0_co2_ppm,t0_ch4_ppb) %>% 
  mutate(ch4=t0_ch4_ppb/1000,
         co2=t0_co2_ppm) %>% 
  select(core_id, ch4, co2) %>% 
  pivot_longer(cols = c(co2,ch4), names_to = "gas",values_to = "t0_ppm") %>% 
  mutate(t0_se_ppm=case_when(gas=="co2"~3.5/sqrt(30),
                          gas=="ch4"~(0.6/1000)/sqrt(30)),
         t0_estimate_origin="Average of remark")%>% 
  select(gas, core_id, t0_ppm, t0_se_ppm, t0_estimate_origin)
  


#JOIN concentrations----

t0<- rbind(t0uveg_co2_ch4,t0ub_n2o)
tf<- rbind(tfub,tfuveg)

t0tf<- merge.data.frame(t0, tf, by=c("gas","core_id"), all = T) %>% merge.data.frame(tfmix, by=c("gas","core_id"),all = T)

#In t0tf NAs for gas are implicit, pivot wider will reveal them

#Check missing:
t0tf %>% 
  group_by(gas) %>% 
  summarise(missing=sum(is.na(tf_ppm)), total=n(), percent_missing=(missing/total)*100)

#A few mixes lower than TF concentrations.... what to do?

#Join with incub details for final flux calculation. 

#core diameter is XXX ask Carlos (me suena 7cm, check)










#Calculate average sd n of initial cores and put in different columns side by side of tf summaries.
#For each core_ID, create 2 sets of columns: t0 and tf with _avg, _sd and _n 
initialfinal_ppm_summary<- injectionsUB_ppm %>% 
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














#_________________####
#FLUX CALCULATION####
# ----- 1. Import & summary -----

injections_ppm<- read.csv(paste0(results_UB_path,"/N2O_CO2_CH4_ppm_all_injections.csv"))
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
  # filter(source_t0=="averaget0s") %>% 
  mutate(t0_cv=t0_sd/t0_avg) %>% 
  ggplot(aes(x=t0_cv, fill=source_t0))+
  geom_histogram(position = "identity",bins=50,alpha=0.3)+
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



##---Decide CH4&CO2 t0 to analyze----
#Check CV of allt0 
deltaGHG_ppm_pval %>% 
  ungroup() %>% 
  filter(gas=="ch4") %>% 
  # filter(source_t0=="averaget0s") %>% 
  mutate(t0_cv=t0_sd/t0_avg) %>% 
  ggplot(aes(x=t0_cv, fill=source_t0))+
  geom_histogram(position = "identity",bins=50,alpha=0.3)+
  geom_vline(xintercept=0.05)+
  facet_wrap(.~gas, scales="free")



#Check CV of allt0 CO2
deltaGHG_ppm_pval %>% 
  ungroup() %>% 
  filter(gas=="co2") %>% 
  # filter(source_t0=="averaget0s") %>% 
  mutate(t0_cv=t0_sd/t0_avg) %>% 
  ggplot(aes(x=t0_cv, fill=source_t0))+
  geom_histogram(position = "identity",bins=50,alpha=0.3)+
  geom_vline(xintercept=0.05)+
  facet_wrap(.~gas, scales="free")


unreliable_t0_avg_ch4<- deltaGHG_ppm_pval %>% 
  ungroup() %>% 
  filter(gas=="ch4") %>% 
  filter(source_t0=="averaget0s") %>% 
  mutate(t0_cv=t0_sd/t0_avg) %>% 
  filter((t0_cv>0.05&p_value>0.01)|t0_cv>0.1) %>% #cores with slightly variable t0 average and unclear flux OR cores with very variable t0 average 
  mutate(available_exetainer=grepl("S2|S3|S4-CA", coreID))



unreliable_t0_avg_co2<- deltaGHG_ppm_pval %>% 
  ungroup() %>% 
  filter(gas=="co2") %>% 
  filter(source_t0=="averaget0s") %>% 
  mutate(t0_cv=t0_sd/t0_avg) %>% 
  filter((t0_cv>0.05&p_value>0.01)|t0_cv>0.1) %>% #cores with slightly variable t0 average and unclear flux OR cores with very variable t0 average 
  mutate(available_exetainer=grepl("S2|S3|S4-CA", coreID))

#Potential Cores to repeat if we use t0 CO2 and CH4 from exetainers
cores_to_analyze_co2_ch4<- rbind(unreliable_t0_avg_ch4,unreliable_t0_avg_co2) %>% 
  select(coreID, available_exetainer) %>% 
  distinct() 

#42 cores t0 needed to analyse IF we use t0 concentrations of CO2&CH4 from exetainers
sum(cores_to_analyze_co2_ch4$available_exetainer)

#9 cores for which we cannot recover accurate t0 initial concentrations
sum(!cores_to_analyze_co2_ch4$available_exetainer)

#Out of 9 cores without accurate t0, 2 ch4 fluxes would remain unclear: from sampling S4-DA-A1
unreliable_t0_avg_ch4 %>% filter(available_exetainer==F) %>% mutate(clearflux=p_value<0.01) %>% select(-sampling)

#Out of 9 cores without accurate t0, 2 co2 fluxes would remain unclear: from sampling S4-VA-R2
unreliable_t0_avg_co2 %>% filter(available_exetainer==F) %>% mutate(clearflux=p_value<0.01) %>% select(-sampling)



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

