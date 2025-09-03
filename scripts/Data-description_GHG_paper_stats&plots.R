#DATA-DESCRIPTION paper: Stats and plots

#Author: Miguel Cabrera-Brufau
#Date: July 2025
#Project: Restore4cs


#Description----
#This scrip is used to explore and plot the main features of the GHG dataset, including:
    #In-situ instantaneous fluxes (CO2, CH4, N2O), Centered on Best-flux, comparing over/underestimation of 3 Models 
    #Ex-situ core incubation fluxes (CO2, CH4, N2O)
    #Water GHG concentrations from Headspace method (CO2, CH4, N2O)
#Conservation status will not be used to explain data


#Ideas for Main figures:

#In-situ Chambers:
  #1. Density plots (by GHG), potentially broken down by lightcondition
  #1b (Doubt): Density plots by GHG and strata (and lightcondition?)
  
#Ex-situ Cores
  #1. Density plots by GHG

#water [GHG]
  #1. Density plots by GHG


#Inter-method comparisons:
  
  #1. Model choice impact chambers: calculate and present over/underestimation for each model, by comparing the relative bias of the two alternative models. Group by Gas, group by best-model. Use general formula:

            #Rel_bias = abs(nonbest - best)/ abs(best) 
    #potential issues with: near-zero best and sign-differences between models (any occurrence?)

  #2.(DOUBT) Chamber vs cores (correlation and biass) by GHG. Issue: unclear 1:1 sample correspondence. Doubt lightconditions (any veg, bare-trans that corresponds to cores?)  

  #3. Chambers (OW) vs water[GHG] correlation 



#DATA PREPARATION: 

#GHG: need all estimates and choice of best (CO2, CH4, N2O)


#TO adapt: -------
#1. Import all fluxes and combine in single table (UniqueID, ghg_best, ghg_se, ghg_model) 
#2. Import harmonized field observations. 


#Inputs: 
#Per-UniqueID: co2, ch4 and n2o best fluxes (and se flux)
#Per-UniqueID and per-plotcode harmonized field observations (created/Updated by  script Harmonize_Fieldsheet-veg_perPLOT_and_perUniqueID). 

rm(list = ls()) # clear workspace
cat("/014") # clear console


# ---- Packages ----
#For data-handling & plotting
library(tidyverse)
library(readxl)
library(ggplot2)
# library(ggpubr)
# library(ggResidpanel)
# library(ggExtra)

#NOT sure if needed: Clean as I progress
# library(lubridate)
# library(zoo)
# library(grid)
# library(egg)
# require(dplyr)
# require(purrr)
# require(data.table)
# require(tools)
# library(hms)
# library(suncalc)
# #For modelling: 
# library(lmerTest) #similar to lme4 but allows to estimate significance of predictors (and interactions)
# library(bestNormalize)
# library(outliers) #Outlier exploration
# library(rstatix) #Outlier exploration
# library(emmeans)
# library(glmmTMB)
# library(DHARMa)
# library(performance)
# library(tibble)
# library(partR2)#for breaking down variance by fixed effect
# library(car)  # for leveneTest
# library(caret)  # for cross-validation
# library(MuMIn)  # for AICc if needed
# library(multcomp)# for emmeans post-hoc comparisons


#Repo-functions
repo_root <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
files.sources = list.files(path = paste0(repo_root,"/functions"), full.names = T)
for (f in files.sources){source(f)}


# ---- Directories ----
dropbox_root <- "C:/Users/Miguel/Dropbox/RESTORE4Cs - Fieldwork/Data" # You have to make sure this is pointing to the write folder on your local machine

#Path to aquaGHG best.flux results:
results_path <- paste0(dropbox_root,"/GHG/Processed data/computed_flux/")
#Path to N2O flux results:
n2o_results_path <- paste0(dropbox_root,"/GHG/N2O_fluxes/")

#Path to save plots from exploration: 
plots_path <- paste0(dropbox_root,"/GHG/Processed data/computed_flux/Stats and plots SEFS14")

#__________________------
#DATA PREP------

#1. Import inst. fluxes-----
#1. Import bestflux results

#All datasets must be identified at the UniqueID level
co2_all<- read.csv(paste0(results_path,"co2_bestflux.csv")) #umol/m2/s
ch4_all<- read.csv(paste0(results_path,"ch4_bestflux.csv")) #nmol/m2/s

#Check N2O identification scheme: UniqueID_notime
n2o_all<- read.csv(paste0(n2o_results_path,"S4_restore4cs_N2O_arealflux_nmol_s-1_m-2.csv"))


#Format to join with others: UniqueID, ghgspecies, ghgmodel(LM/HM/total.flux) bestmodel(T/F), flux
#CO2
co2_long<- co2_all %>% 
  filter(best.model!="None appropriate") %>% #REMOVE wrong incubations
  select(UniqueID, LM.flux, HM.flux, total.flux, best.model) %>% 
  rename(LM=LM.flux, HM=HM.flux, total=total.flux) %>%
  mutate(best.model=if_else(best.model=="total.flux","total",best.model)) %>% 
  pivot_longer(cols = c(LM, HM, total), names_to = "ghgmodel", values_to = "flux") %>% 
  mutate(bestmodel=ghgmodel==best.model, #Create logic bestmodel
         ghgspecies="co2") %>% #Add ghgspecies
  #Get UniqueID_notime
  separate(UniqueID, into = c("c1","c2","c3","c4","c5","c6","c7"), sep = "-", remove = F) %>% 
  mutate(UniqueID_notime=paste(c1,c2,c3,c4,c5,c6,sep = "-")) %>% 
  dplyr::select(-c(c1,c2,c3,c4,c5,c6,c7)) %>% 
  select(UniqueID_notime,ghgspecies,ghgmodel,bestmodel, flux)
  

#CH4
ch4_long<- ch4_all %>% 
  filter(best.model!="None appropriate") %>% #REMOVE wrong incubations
  select(UniqueID, LM.flux, HM.flux, total.flux, best.model) %>% 
  rename(LM=LM.flux, HM=HM.flux, total=total.flux) %>%
  mutate(best.model=if_else(best.model=="total.flux","total",best.model)) %>% 
  pivot_longer(cols = c(LM, HM, total), names_to = "ghgmodel", values_to = "flux") %>% 
  mutate(bestmodel=ghgmodel==best.model, #Create logic bestmodel
         ghgspecies="ch4") %>% #Add ghgspecies
  #Get UniqueID_notime
  separate(UniqueID, into = c("c1","c2","c3","c4","c5","c6","c7"), sep = "-", remove = F) %>% 
  mutate(UniqueID_notime=paste(c1,c2,c3,c4,c5,c6,sep = "-")) %>% 
  dplyr::select(-c(c1,c2,c3,c4,c5,c6,c7)) %>% 
  select(UniqueID_notime,ghgspecies,ghgmodel,bestmodel, flux)

#N2O
n2o_long<- n2o_all %>% 
  mutate(ghgspecies="n2o",
         ghgmodel="total",
         bestmodel=T,
         flux=N2Oflux_nmol_per_second_per_m2
          ) %>% 
  select(UniqueID_notime ,ghgspecies,ghgmodel,bestmodel, flux)


#JOIN all ghg fluxes per UniqueID (useful for model comparison)
ghg_long<- co2_long %>% 
  rbind(ch4_long) %>% 
  rbind(n2o_long) %>% 
  mutate(light=substr(UniqueID_notime, nchar(UniqueID_notime), nchar(UniqueID_notime))) %>%
  mutate(strata=case_when(grepl("v",UniqueID_notime)~"vegetated",
                          grepl("b",UniqueID_notime)~"bare",
                          grepl("o",UniqueID_notime)~"open water"))

#Filter for bestmodel only
ghg_bestlong<- ghg_long %>% 
  filter(bestmodel==T)

#rm(co2_all, co2_long, ch4_all, ch4_long, n2o_all, n2o_long)




#__________________-----

#PLOTS -------

#__________________-----



#A. Chambers------
##1. Density plots by gas (best)-----

#Need log10scales (positive and negative), Choose units to best show distribution when log-transformed

#TO-DO:------
#Create various log10trans options (using different flux units) to inspect and select best way of presenting the data. Create functions and use functions within ggplot call and breaks and labels. 


#CO2_untransformed
ghg_bestlong %>% 
  filter(ghgspecies=="co2") %>% 
  ggplot(aes(x=flux, fill=light))+
  geom_histogram( color="#e9ecef", alpha=0.5, position = 'identity') +
  scale_fill_manual(values=c("#404080", "orange"))+
  theme_bw()


ghg_bestlong %>% 
  filter(ghgspecies=="co2") %>% 
  ggplot(aes(x=sign(flux)*log10(abs(flux*10)+1), fill=light))+
  geom_histogram( color="#e9ecef", alpha=0.5, position = 'identity') +
  scale_fill_manual(values=c("#404080", "orange"))+
  theme_bw()

ghg_bestlong %>% 
  filter(ghgspecies=="co2") %>% 
  ggplot(aes(x=sign(flux)*log10(abs(flux*10)+1), fill=light))+
  geom_histogram( color="#e9ecef", alpha=0.5, position = 'identity') +
  scale_fill_manual(values=c("#404080", "orange"))+
  theme_bw()+
  facet_wrap(~strata)
  

#Alternative: geom_density with after_stat(count) to show "smoothed" histograms
#UNCLEAR WHY COUNT of geom_density does not match that of geom_histogram (possible effect of binwidth + logscale)
ghg_bestlong %>% 
  filter(ghgspecies=="co2") %>% 
  ggplot(aes(x=sign(flux)*log10(abs(flux*10)+1),y=after_stat(count), fill=light))+
  geom_density( color="#e9ecef", alpha=0.5, position = 'identity') +
  geom_histogram(position = 'identity')+
  scale_fill_manual(values=c("#404080", "orange"))+
  theme_bw()+
  facet_wrap(~strata)





#CH4_untransformed
ghg_bestlong %>% 
  filter(ghgspecies=="ch4") %>% 
  ggplot(aes(x=flux, fill=light))+
  geom_histogram( color="#e9ecef", alpha=0.5, position = 'identity') +
  scale_fill_manual(values=c("#404080", "orange"))+
  theme_bw()


ghg_bestlong %>% 
  filter(ghgspecies=="ch4") %>% 
  ggplot(aes(x=sign(flux)*log10(abs(flux*10)+1), fill=light))+
  geom_histogram( color="#e9ecef", alpha=0.5, position = 'identity') +
  scale_fill_manual(values=c("#404080", "orange"))+
  theme_bw()

ghg_bestlong %>% 
  filter(ghgspecies=="ch4") %>% 
  ggplot(aes(x=sign(flux)*log10(abs(flux*10)+1), fill=light))+
  geom_histogram( color="#e9ecef", alpha=0.5, position = 'identity') +
  scale_fill_manual(values=c("#404080", "orange"))+
  theme_bw()+
  facet_wrap(~strata)


#N2O_untransformed
ghg_bestlong %>% 
  filter(ghgspecies=="n2o") %>% 
  ggplot(aes(x=flux, fill=light))+
  geom_histogram( color="#e9ecef", alpha=0.5, position = 'identity') +
  scale_fill_manual(values=c("#404080", "orange"))+
  theme_bw()


ghg_bestlong %>% 
  filter(ghgspecies=="n2o") %>% 
  ggplot(aes(x=sign(flux)*log10(abs(flux*10)+1), fill=light))+
  geom_histogram( color="#e9ecef", alpha=0.5, position = 'identity') +
  scale_fill_manual(values=c("#404080", "orange"))+
  theme_bw()

ghg_bestlong %>% 
  filter(ghgspecies=="n2o") %>% 
  ggplot(aes(x=sign(flux)*log10(abs(flux*10)+1), fill=light))+
  geom_histogram( color="#e9ecef", alpha=0.5, position = 'identity') +
  scale_fill_manual(values=c("#404080", "orange"))+
  theme_bw()+
  facet_wrap(~strata)





##2. Model choice effects/bias -----

#Caculate bias caused by model choice (with respect to bestmodel identified)



#CO2: only LM vs HM, easy in 1 plot (overestimation HM (when LM is best), underestimation LM (when HM is best)), issues for near-zero denominator flux. POtentially biplot with 1:1 line (but how to convey which is the best model?, Two panels? Different color?). THis would allow for logscales using always absolute flux magnitudes (do not care about sign for this model comparison, effect is symetric)

#CH4: need for 3 options, check how to present. 




#B. CORES -------

#Density plots per gas




#C. WATER concentrations------

#Density plots per gas
#Under/oversaturated as colour? Can we calculate it and show it in the same graph?

