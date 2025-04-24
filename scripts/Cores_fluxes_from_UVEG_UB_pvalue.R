
# ---
# Authors: MIGUEL CABRERA
# Project: "RESTORE4Cs"
# date: "April 2025"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of script ----
# This script computes per-area daily fluxes of CO2, CH4 and N2O for all core incubations. As well as the CH4 stock in the core. 

#ORIGIN of CONCENTRATIONS (ppm)
#All source concentrations used for flux calcualtion can be found (properly formatted) in the csv All_core_GHG_concentrations_ppm_initial_final_mix.csv in the dropbox Core folder. 

#N2O initial concentrations are derived from samples analyzed in the UB. T0 N2O concentrations are either measured (initial exetainer) or averaged from same-subsite cores (after checking the variability of initial exetainers, minimum analyzed: 3 out of 6, when CV of the 3 initial exetainers exceeds 6%, the remaining 3 initial exetainers for that subsite were also analyzed, only 1 case).

#All CO2 and CH4 initial concentrations are obtained from the in-situ measurements with portable gas analyzer Li-cor/Los Gatos/Picarro (inlet inside core, measure taken after stabilization of measure inside the open core arising from air flushing by analyzer) the concentration of CO2 and CH4 is written down, an exetainer is taken for N2O and the core is closed to start the incubation.

#N2O final concentrations are measured from final exetainers in the UB with zero-air baseline. 
#CO2 and CH4 final concentrations are derived from UB measurements with low baseline for Seasons 2-4. For Season 1, final concentrations of CO2 and CH4 are derived from in-situ injections (with ambient air baseline and sampling directly from core).

#Stocks are calculated as the difference in CH4 concentration between the initial time and the homogenized core (manual shaking until visible mix of sediment, water and air). Only cores with water are homogenized for CH4 stocks. 

#Incubation details (duration, headspace, water) are measured in-situ and written down, for all incubated cores. 

# Incubation conditions are not homogeneous, depending on case-pilot site a temperature controlled chamber is used or cores are incubated outside the lab at ambient temperature and light conditions.  

#Concentrations derived from injection of samples using the open-loop method were analyzed with the scripts available at: https://github.com/MCabreraBrufau/Measuring_discrete_GHG_samples_with_Licor.git
#For injection sequences of final samples where no peak is detected, the average value of the remark is taken. 


#Fluxes are calculated based on core-dimensions, incubation time and GHG concentration difference between initial and final samples. See section FLUX CALCULATION for details. 



rm(list = ls()) # clear workspace
cat("/014") # clear console

# ---- packages ----
library(tidyverse)
library(readxl)
library(ggplot2)



# ---- Directories ----
dropbox_root <- "C:/Users/Miguel/Dropbox/RESTORE4Cs - Fieldwork/Data" # You have to make sure this is pointing to the write folder on your local machine
results_UB_path<- paste0(dropbox_root,"/Cores/UB_concentrations/")
results_uveg_path<- paste0(dropbox_root,"/Cores/UVEG_concentrations/")
core_path<- paste0(dropbox_root,"/Cores/")





#TO-DO--------
#Update fieldsheet missing details (when Carlos fills them)
#Format meta-data (sampling data & incubation details, ask Carlos to decide what to include). 


#_________________####


#Incubation-data --------

incub<- read_xlsx(path = paste0(core_path,"Incubation_details_all_coresS1S2S3S4_para rellenar.xlsx")) %>% 
  select(core_id, headspace_cm, water_cm, approximatesediment_cm, coreincubation_hours,comment, t0_co2_ppm, t0_ch4_ppb)



#Sampling-data-----
#Load core-fieldsheet compilation (created with Compile_Core_sampling_details.R)

field<- read.csv(file = paste0(core_path, "Core_fieldsheet_compilation.csv"))

#List subsites with missing sampling info in fieldsheets: 
incub[!incub$core_id%in%field$sample_id,] %>%
  select(core_id) %>% 
  mutate(subsite=substr(core_id, 1,8)) %>%
  group_by(subsite) %>% 
  summarise(cores_missing=n())
#4 subsites without any info for core sampling (core_id not in fieldsheet)

#Subsites with missing details (other than comments)
incomplete_details<- field[!complete.cases(field[,!names(field) %in% "comments"]),] %>% filter(!is.na(pilot_site)) 
message(paste("The following subsites have incomplete information for cores:", unique(incomplete_details$subsite_code)))
#1 subsite with missin details.




#Method CV -----
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





#Calculate SE of t0 measured in situ using instrument precission for each gas and assuming 10s average for each t0 stabilized reading. 
t0uveg_co2_ch4<- incub %>% 
  select(core_id, t0_co2_ppm,t0_ch4_ppb) %>% 
  mutate(ch4=t0_ch4_ppb/1000,
         co2=t0_co2_ppm) %>% 
  select(core_id, ch4, co2) %>% 
  pivot_longer(cols = c(co2,ch4), names_to = "gas",values_to = "t0_ppm") %>% 
  mutate(t0_se_ppm=case_when(gas=="co2"~3.5/sqrt(10),
                          gas=="ch4"~(0.6/1000)/sqrt(10)),
         t0_estimate_origin="Average of remark")%>% 
  select(gas, core_id, t0_ppm, t0_se_ppm, t0_estimate_origin)
  


#Save concentrations----

t0<- rbind(t0uveg_co2_ch4,t0ub_n2o)
tf<- rbind(tfub,tfuveg)

cores_ppm<- merge.data.frame(t0, tf, by=c("gas","core_id"), all = T) %>% merge.data.frame(tfmix, by=c("gas","core_id"),all = T)

#Save all measured concentrations with proper format (Error propagated for each concentration estimate_origin)
write.csv(cores_ppm,file = paste0(core_path, "Cores_flux/All_core_GHG_concentrations_ppm_initial_final_mix.csv"),row.names = F)

#Remove intermediate data.frames
rm(t0, tf,tfmix, t0ub_n2o, t0uveg_co2_ch4,tfub, tfuveg, mixesUVEG_ppm, open_loop_cal, samplesUB_ppm, samplesUVEG_ppm, n2o_cv_method, ch4_cv_method, co2_cv_method)



#FLUX CALCULATION####

#Calculate fluxes of the 3 GHG and stocks of CH4, as well as the significance (test if they are different than cero using the propagated error)

#Core dimensions: core inner radius is 1.75cm (confirmed with Carlos), headspace is given in cm

#Convert to molar units using Ideal Gas Law: PV = nRT, where:
#P: Pressure (atm), assumed to be 1 atm
#V: Volume (L) determined from core dimensions
#n: number of mols
#T: Temperature (ºK), assumed to be 20ºC (293.15ºK)
#R: Gas Constant,  R = 0.08206 (L·atm·K-1·mol-1)


core_fluxes<- incub %>% 
  select(core_id, headspace_cm, water_cm, approximatesediment_cm, coreincubation_hours) %>% 
  #Get corevol L: h*pi*r^2 (in cm^3), divided by 1000 to get L
  mutate(corevol_litres= headspace_cm*pi*(1.75^2)*0.001) %>% 
  merge.data.frame(cores_ppm, by="core_id") %>% 
  select(-c(t0_estimate_origin, tf_estimate_origin, tfmix_estimate_origin)) %>% 
  #Format into long and fix naming
  pivot_longer(cols = c(t0_ppm,t0_se_ppm, tf_ppm,tf_se_ppm,tfmix_ppm,tfmix_se_ppm),names_to = "measuretype",values_to = "GHG_ppm") %>% 
  mutate(measuretype=gsub("_ppm","", measuretype)) %>% 
  separate(measuretype, into = c("coretime", "measuretype"),fill = "right") %>% 
  mutate(measuretype=if_else(is.na(measuretype),true = "avg", false = measuretype)) %>% 
  #Calculate mols of GHG (1st obtain Litres of gas, then convert to mols assuming 20ºC (293,15ºK) and 1 atm of pressure. PV =nRT, so n=(PV)/RT. Using R= 0.08206 with units: litres*atm/ºK*mol 
  mutate(GHG_litres = (GHG_ppm/1e6)*corevol_litres,
         GHG_mols = (1*GHG_litres)/(293.15*0.08206)) %>% 
  #Select only relevant variables and pivot wider for mol differences
  select(-c(corevol_litres, GHG_ppm, GHG_litres)) %>% 
  pivot_wider(names_from = c(coretime, measuretype),values_from = GHG_mols) %>% 
  #calculate mol difference between t0 and tf, and between t0 and tfmix (propagating the errors)
  mutate(deltamol_t0tf_avg=tf_avg-t0_avg,
         deltamol_t0tf_se= sqrt(tf_se^2+t0_se^2),
         deltamol_t0tfmix_avg=tfmix_avg-t0_avg,
         deltamol_t0tfmix_se= sqrt(tfmix_se^2+t0_se^2)) %>% 
  #Calculate fluxes (mol m^-2 d^-1) and stocks (mol m^-2). Using core-inner area in m^2 and core-incubationtime in days
  mutate(flux_mol_per_m2_per_d=(deltamol_t0tf_avg/(pi*(1.75/100)^2))/(coreincubation_hours/24),
         flux_SE=(deltamol_t0tf_se/(pi*(1.75/100)^2))/(coreincubation_hours/24),
         stock_mol_per_m2=deltamol_t0tfmix_avg/(pi*(1.75/100)^2),
         stock_SE=deltamol_t0tfmix_se/(pi*(1.75/100)^2)) %>% 
  #remove intermediate variables:
  select(-c(deltamol_t0tf_avg,deltamol_t0tf_se,deltamol_t0tfmix_avg,deltamol_t0tfmix_se, t0_avg,t0_se, tf_avg, tf_se, tfmix_avg,tfmix_se)) %>% 
  #Calculate significance of estimates (are they different than cero?)
  #For fluxes we take the absolute value and calcaulate we performed a one-sided z-test using the absolute flux (calculated from the difference in final and initial concentration) divided by the propagated standard error derived from analytical precision of the method and the replicates (
  
  #For GHG fluxes, we test if the absolute flux is significantly greater than zero, we perform a one-sided z-test using the absolute value of the difference divided by the propagated standard error. The propagated standard error accounts for both the analytical precision of the method and the standard deviation among analytical replicates. Significance was evaluated using the upper tail of the standard normal distribution.”
  #For CH4 stocks, we follow a similar procedure to test if the stock is significantly greater than zero. 
  mutate(flux_pvalue= 1-pnorm(abs(flux_mol_per_m2_per_d)/flux_SE),
         stock_pvalue = 1-pnorm(stock_mol_per_m2/stock_SE)) %>% 
  #Format into wide form for each GHG 
  pivot_wider(names_from = gas, values_from = c(flux_mol_per_m2_per_d,flux_SE,flux_pvalue, stock_mol_per_m2,stock_SE,stock_pvalue),names_vary = "slowest") %>% 
  #Rename final variables
  rename(CH4_flux_mol_per_m2_per_d=flux_mol_per_m2_per_d_ch4,
         CH4_flux_SE=flux_SE_ch4,
         CH4_flux_pvalue=flux_pvalue_ch4,
         CO2_flux_mol_per_m2_per_d=flux_mol_per_m2_per_d_co2,
         CO2_flux_SE=flux_SE_co2,
         CO2_flux_pvalue=flux_pvalue_co2,
         N2O_flux_mol_per_m2_per_d=flux_mol_per_m2_per_d_n2o,
         N2O_flux_SE=flux_SE_n2o,
         N2O_flux_pvalue=flux_pvalue_n2o,
         CH4_stock_mol_per_m2=stock_mol_per_m2_ch4,
         CH4_stock_SE=stock_SE_ch4,
         CH4_stock_pvalue=stock_pvalue_ch4) %>% 
  #Transform units to get reasonable values: CH4 and N2O in micromol/m2/d, CO2 in mmol /m2/d
  #CH4 stock and in mmol /m2
  mutate(CH4_flux_micromol_per_m2_per_d= CH4_flux_mol_per_m2_per_d*1e6,
         CH4_flux_SE= CH4_flux_SE*1e6,
         N2O_flux_micromol_per_m2_per_d= N2O_flux_mol_per_m2_per_d*1e6,
         N2O_flux_SE=N2O_flux_SE*1e6,
         CO2_flux_mmol_per_m2_per_d= CO2_flux_mol_per_m2_per_d*1e3,
         CO2_flux_SE=CO2_flux_SE*1e3,
         CH4_stock_mmol_per_m2=CH4_stock_mol_per_m2*1e3,
         CH4_stock_SE=CH4_stock_SE*1e3) %>% 
  #Leave only fluxes and stocks 
  select(core_id,
         CO2_flux_mmol_per_m2_per_d,CO2_flux_SE,CO2_flux_pvalue,
         CH4_flux_micromol_per_m2_per_d,CH4_flux_SE,CH4_flux_pvalue,
         N2O_flux_micromol_per_m2_per_d,N2O_flux_SE,N2O_flux_pvalue,
         CH4_stock_mmol_per_m2,CH4_stock_SE,CH4_stock_pvalue)
  

#Save core fluxes
write.csv(core_fluxes,file = paste0(core_path, "Cores_flux/All_core_GHGfluxes_and_CH4stocks.csv"),row.names = F)

#_________________####



core_fluxes %>% 
  select(CO2_flux_pvalue, CH4_flux_pvalue, N2O_flux_pvalue) %>% 
  pivot_longer(cols = names(.), values_to = "value", names_to = "variable") %>% 
  group_by(variable) %>% 
  summarise(cores=n(), measures=sum(!is.na(value)),percent_significant=sum(value<0.05, na.rm=T)/measures*100, percent_NAs=100-(measures/cores*100))





