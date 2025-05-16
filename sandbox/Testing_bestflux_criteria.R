
# ---
# Authors: Miguel Cabrera
# Project: "RESTORE4Cs"
# date: "April 2025"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script

# This script analyses the quality of fit of all fluxes to check the aproppriateness of selection criteria of best.flux function. This is based on the per-GHG cropped fluxes calculated.

#Biggest risk is overestimating fluxes by choosing HM over LM when is the HM curvature is exaggerated due to artefacts. The Exaggeration of curvature can be assessed based on the K.ratio and based on g-factor. 

#CO2 criteria-----
#CO2 flux selection: after inspection and cropping of artefacts 

#Removal of wrong incubations: incubations with unreliable timeseries (due to artefacts or errors in manipulations) are logged in the same excell for cropping decissions. No CO2 flux will be included for these incubations. 

#Hard thresholds: 
#1. Exaggerated curvature: K.ratio > 1 --> default to LM
#2. Exaggerated curvature: gfactor > 4 --> default to LM. All fluxes with g.fact > 4 have been inspected for apropriateness of LM over HM.
#3. Minimum detectable flux (MDF): LM flux < MDF.lim --> default to LM. All fluxes with MDF have been inspected for aproppriateness of LM over HM. Rationale: the choice of model should be the defining factor of a significant/insignificant flux. Additionally, if the flux is below detection, the premise of HM model (i.e. concentration gradient reduction during incubation leading to progresively lower flux over time) is extremely unlikely.


#Scoring system: we will use AICc improvement (via AICc weights) to select HM over LM when doing so improves fit by more than 5%. 

#We will only assess the scoring system with the fluxes that do not breach any of the hard thresholds. 








# ---- packages ----
rm(list = ls()) # clear workspace
cat("/014") # clear console

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
library(ggpubr)
library(ggExtra)


repo_root <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
files.sources = list.files(path = paste0(repo_root,"/functions"), full.names = T)
for (f in files.sources){source(f)}


# ---- Directories and data loading ----
dropbox_root <- "C:/Users/Miguel/Dropbox/RESTORE4Cs - Fieldwork/Data" # You have to make sure this is pointing to the write folder on your local machine
fieldsheetpath <- paste0(dropbox_root,"/GHG/Fieldsheets")
quality_path<- paste0(dropbox_root,"/GHG/Working data/Incubation_quality/")
# plots_path <- paste0(dropbox_root,"/GHG/Processed data/computed_flux/plots")
results_path<- "C:/Users/Miguel/Dropbox/TEST_quality_raw2flux/Per_GHG_cropped_fluxes/"


#Load per-GHG cropped fluxes
listf <- list.files(path = results_path, pattern = "^S_fluxes", all.files = T, full.names = T, recursive = F)
table_results_all <- NULL

for (f in listf){
  table_results_all <- rbind(table_results_all, 
                             read.csv(file = f, header = T))}

#Add
table_results_all$sampling <- str_sub(table_results_all$subsite, start = 1, 2)
table_results_all$pilotsite <- str_sub(table_results_all$subsite, start = 4, 5)
table_results_all$subsite <- str_sub(table_results_all$subsite, start = 7, 8)
table_results_all$siteID <- str_sub(table_results_all$subsite, start = 4, 8)


table_results_all <- table_results_all[!is.na(table_results_all$lightCondition),]

#Load subsites level_incubation CSV files with details of flux calulation for different models (produced by raw2flux.R script)
get_full_detailed_table <- function(table_results_all, variable){
  
  listf <- list.files(path = paste0(results_path,"level_incubation"), pattern = ".csv", all.files = T, full.names = T, recursive = F)
  mytable <- NULL
  for (f in listf[grep(pattern = variable, x = listf)]){
    mytable.tmp <- read.csv(file = f, header = T)
    mytable <- rbind(mytable, mytable.tmp)
  }
  mytable$sampling <- str_sub(mytable$UniqueID, start = 1, 2)
  mytable$pilotsite <- str_sub(mytable$UniqueID, start = 4, 5)
  mytable$subsite <- str_sub(mytable$UniqueID, start = 7, 8)
  mytable$siteID <- str_sub(mytable$UniqueID, start = 4, 8)
  
  ind_match <- match(mytable$UniqueID, table_results_all$UniqueID)
  
  
  mytable <- cbind(mytable, table_results_all[ind_match,c("gas_analiser","start.time","duration","water_depth","strata","chamberType","lightCondition")])
  mytable <- mytable[!is.na(ind_match),]
  return(mytable)
}


#Get detailed table for each gas
table_co2 <- get_full_detailed_table(table_results_all, variable = "co2")
# table_ch4 <- get_full_detailed_table(table_results_all, variable = "ch4")



#Get Wrong CO2------

#Get wrong incubations: those with too much influence of artefacts or wrong manipulations based on visual inspection of time-series.
inspectiontable<- read_xlsx(path = paste0(quality_path, "Inspection_table_allincubations_tocrop_perGHG.xlsx")) %>% 
  select(UniqueID,co2_decission,co2_obs_inspection, ch4_decission,ch4_obs_inspection)

co2_discard<- inspectiontable %>% filter(co2_decission=="discard") %>% pull(UniqueID)

#these incubations will be flagged in the final dataset and will not influence the criteria decisions. 



#Criteria Hard.threshold checks

#Inspect mdf----
#We want to make sure MDF fluxes do not arise from artefacts

co2_mdf<- table_co2 %>%  
  filter(!UniqueID%in%co2_discard) %>%
  filter(abs(LM.flux)<=MDF.lim)

#MDFs CO2 inspected: true below detection flat-ish white noise incubations
mdf_inspected<- c("s1-ca-a2-1-v-d-10:44","s1-du-p2-12-b-d-11:01","s1-ri-a2-1-b-d-10:20","s1-ri-p1-1-v-d-10:19","s1-ri-p1-3-v-d-10:51","s1-ri-r2-4-b-d-10:01","s2-du-a1-15-v-t-13:29","s2-du-a2-1-o-d-08:54","s3-du-p2-7-b-d-10:26","s3-va-p1-13-v-t-11:48","s1-ca-a2-13-o-d-13:17","s1-ca-a2-4-v-d-11:25","s1-ca-a2-5-o-d-11:33","s1-ca-a2-6-v-d-11:49","s1-ca-a2-8-o-d-12:12","s1-ca-r1-1-o-d-08:41","s1-ca-r1-4-o-d-09:45","s1-ca-r1-6-o-d-10:17","s1-cu-a1-14-v-d-11:14","s1-cu-a2-14-o-d-11:32","s1-cu-a2-16-o-d-11:58","s1-cu-a2-4-o-d-07:53","s1-cu-a2-6-o-d-08:12","s1-cu-p2-9-o-d-10:56","s1-cu-r2-18-v-d-10:36","s1-cu-r2-21-b-d-11:20","s1-da-p1-3-o-d-09:11","s1-da-p1-5-o-d-09:33","s1-da-p1-8-o-d-10:05","s1-da-p1-9-o-d-10:17","s1-du-a1-2-o-d-09:26","s1-du-a1-3-o-d-09:38","s1-du-a1-4-o-d-09:50","s1-du-p2-5-b-t-08:39","s1-du-p2-6-v-t-08:52","s1-du-r2-8-v-t-09:44","s1-ri-a2-10-b-d-11:50","s1-ri-a2-12-b-d-12:05","s1-ri-a2-5-b-d-10:56","s1-ri-a2-6-b-d-11:01","s1-ri-p1-10-v-t-12:59","s1-ri-p1-4-v-d-11:11","s1-ri-r2-6-v-d-10:39","s1-ri-r2-7-b-d-10:50","s1-va-a2-1-o-d-09:02","s1-va-a2-4-v-d-10:01","s1-va-a2-6-o-d-10:35","s1-va-p1-13-v-t-12:41","s1-va-p2-4-b-t-09:56","s1-va-r2-13-o-d-12:34","s1-va-r2-14-o-d-12:53","s1-va-r2-15-o-d-13:13","s1-va-r2-8-b-d-11:08","s2-ca-a1-3-v-t-09:12","s2-ca-a1-8-b-d-11:08","s2-ca-a2-1-v-d-08:15","s2-ca-a2-10-v-d-10:58","s2-ca-a2-7-v-d-09:47","s2-ca-a2-9-v-d-10:31","s2-ca-r1-14-v-d-12:21","s2-cu-a1-4-v-t-08:59","s2-cu-a2-11-o-d-10:30","s2-cu-a2-14-o-d-11:03","s2-cu-a2-15-o-d-11:11","s2-cu-a2-3-v-d-08:54","s2-cu-a2-6-o-d-09:42","s2-cu-a2-7-o-d-09:50","s2-cu-a2-8-o-d-10:04","s2-cu-a2-9-o-d-10:13","s2-cu-p1-2-v-t-07:51","s2-cu-p2-1-v-d-08:26","s2-cu-p2-1-v-t-08:18","s2-cu-p2-10-o-d-10:52","s2-cu-p2-11-o-d-11:08","s2-cu-p2-13-o-d-11:23","s2-cu-p2-14-o-d-11:31","s2-cu-p2-15-o-d-11:39","s2-cu-p2-2-v-d-08:46","s2-cu-p2-2-v-t-08:38","s2-cu-p2-3-v-d-09:08","s2-cu-p2-3-v-t-09:00","s2-cu-p2-9-o-d-10:45","s2-cu-r1-4-v-t-08:40","s2-cu-r2-1-b-d-08:08","s2-cu-r2-1-b-t-08:01","s2-cu-r2-10-o-d-11:10","s2-cu-r2-11-o-d-11:23","s2-cu-r2-12-o-d-11:37","s2-cu-r2-3-b-d-08:41","s2-cu-r2-4-v-d-09:02","s2-cu-r2-4-v-t-08:55","s2-cu-r2-5-v-d-09:23","s2-cu-r2-5-v-t-09:14","s2-cu-r2-6-v-t-09:40","s2-cu-r2-9-o-d-11:00","s2-da-a1-1-o-d-08:17","s2-da-a1-3-o-d-08:43","s2-da-a2-1-o-d-07:32","s2-da-a2-10-o-d-09:31","s2-da-a2-11-o-d-09:57","s2-da-a2-12-o-d-10:06","s2-da-a2-13-o-d-10:17",
                  "s2-da-a2-15-o-d-10:46","s2-da-a2-2-o-d-07:52","s2-da-a2-3-o-d-08:04","s2-da-a2-4-o-d-08:16","s2-da-a2-5-o-d-08:26","s2-da-a2-6-o-d-08:36","s2-da-a2-7-o-d-08:47","s2-da-a2-8-o-d-09:02","s2-da-p1-1-o-d-08:10","s2-da-p2-2-v-t-09:58","s2-da-r2-2-v-d-07:46","s2-du-a1-1-o-d-09:21","s2-du-a1-2-o-d-09:36","s2-du-a1-3-o-d-09:48","s2-du-a2-12-v-d-12:36","s2-du-a2-14-b-d-13:04","s2-du-p2-1-o-d-09:17","s2-du-p2-3-o-d-10:11","s2-du-r1-9-v-d-11:22","s2-ri-a2-1-b-d-08:35","s2-ri-a2-16-o-d-12:08","s2-ri-a2-2-b-d-08:51","s2-ri-a2-7-b-d-10:00","s2-ri-a2-8-b-d-10:10","s2-ri-p1-13-o-d-12:52","s2-ri-p1-14-o-d-13:19","s2-ri-p1-15-o-d-13:41","s2-ri-p1-2-v-d-09:11","s2-ri-p1-3-v-d-09:31","s2-ri-p1-4-v-d-09:54","s2-ri-p1-5-v-d-10:08","s2-ri-p1-7-o-d-10:36","s2-ri-r2-12-v-d-12:40","s2-ri-r2-2-o-d-12:47","s2-ri-r2-3-o-d-13:26","s2-va-a1-2-o-d-10:21","s2-va-a1-4-v-d-11:06","s2-va-a1-4-v-t-10:57","s2-va-a1-6-b-d-11:34","s2-va-a1-7-v-t-11:57","s2-va-p2-12-v-t-11:44","s2-va-p2-13-v-d-12:14","s2-va-p2-15-v-t-12:49","s2-va-p2-5-o-d-10:05","s2-va-p2-6-o-d-10:18","s2-va-p2-7-o-d-10:35","s2-va-p2-9-b-d-11:07","s2-va-r1-1-o-d-09:18","s2-va-r1-11-v-d-12:46","s2-va-r1-11-v-t-12:36","s2-va-r1-12-v-d-13:13","s2-va-r1-12-v-t-13:03","s2-va-r1-5-o-d-10:22","s2-va-r1-8-o-d-11:15","s3-ca-a2-10-o-d-09:29","s3-ca-a2-12-o-d-09:56","s3-ca-a2-14-o-d-10:25","s3-ca-a2-15-v-d-10:39","s3-ca-a2-16-v-d-10:57","s3-ca-a2-2-o-d-07:43","s3-ca-a2-3-v-d-08:05","s3-ca-a2-4-o-d-08:13","s3-ca-a2-5-v-d-08:27","s3-ca-a2-5-v-t-08:20","s3-ca-a2-6-o-d-08:36","s3-ca-a2-7-v-d-08:50","s3-da-a2-1-o-d-06:06","s3-da-a2-4-o-d-07:00","s3-du-a1-6-b-d-09:19","s3-du-a2-11-b-d-10:55","s3-du-a2-3-v-t-07:43","s3-du-p2-5-o-d-09:38","s3-du-r1-9-v-t-09:18","s3-ri-a2-11-b-d-10:58","s3-ri-a2-13-b-d-11:20","s3-ri-a2-14-b-d-11:30","s3-ri-a2-16-o-d-12:03","s3-ri-a2-17-o-d-12:23","s3-ri-a2-6-b-d-09:51","s3-ri-a2-8-b-d-10:21","s3-ri-a2-9-b-d-10:29","s3-ri-r1-13-b-d-14:04","s3-va-r2-1-o-d-08:29","s3-va-r2-3-v-d-09:10","s3-va-r2-3-v-t-09:04","s4-cu-a2-11-o-d-10:31","s4-cu-r2-12-o-d-10:10","s4-da-a2-12-o-d-08:51","s4-da-p1-12-o-d-09:55","s4-ri-a2-6-b-d-07:41","s4-ri-p1-5-o-d-07:54","s2-ca-a2-5-v-t-09:15","s2-ca-a2-2-v-d-08:33")


#any MDF flux not inspected?
co2_mdf %>% filter(!UniqueID%in%mdf_inspected) %>% pull(UniqueID)



#check g.fact ------

#Here we check that the hard-threshold for g.fact (i.e. if g.fact > 4 --> LM) does not cause any truly non-linear pattern to be lost (i.e. that there is no good HM of flux >4*LM)

#These are All the incubations with unreasonably large g-fact (>4). They have all been visually inspected and LM is more representative for them. 
co2_gfact_inspected<- c("s3-va-r1-1-v-d-08:47","s4-ri-p2-6-o-d-08:28","s2-va-a2-14-v-d-13:06","s4-da-p2-5-o-d-09:34","s4-ca-a1-8-v-t-10:52","s1-ca-p2-3-v-t-10:53","s4-da-r2-7-o-d-08:01","s4-du-a1-14-b-d-12:34","s1-va-r1-13-v-d-12:55","s1-ca-p1-14-v-d-14:49","s3-ri-r1-13-b-d-14:04","s2-da-p2-9-o-d-11:15","s4-cu-a1-10-o-d-09:59","s2-cu-p1-13-o-d-10:11","s4-cu-a1-4-v-t-08:22","s2-cu-p1-11-v-t-09:50","s4-ri-p1-14-o-d-11:04","s1-cu-p1-17-v-t-14:01","s1-da-p2-11-o-d-12:38","s2-du-r2-4-b-d-10:07","s4-ca-a2-13-o-d-10:27","s1-da-r2-15-b-t-12:35","s3-da-p2-3-v-d-09:45","s4-da-p2-1-v-t-08:45","s2-cu-a1-15-v-d-11:02","s2-da-p2-1-o-d-09:44","s4-ca-p2-1-v-d-08:01","s2-du-p1-15-v-d-13:01","s4-du-r2-11-v-t-10:49","s4-ri-p1-7-o-d-08:26","s3-cu-r1-2-v-t-07:06","s3-cu-r1-2-v-t-07:06","s1-cu-r1-13-v-t-08:53","s3-cu-r2-12-v-t-11:39","s1-va-a2-12-v-t-12:32","s3-cu-a2-5-v-t-10:37","s4-va-r1-12-v-d-11:39","s4-da-a1-8-v-t-09:40","s4-da-p1-4-v-t-08:06","s1-cu-a1-6-o-d-08:23","s1-du-r2-15-v-d-12:56","s4-ca-p1-5-v-t-08:24","s1-du-r2-11-v-t-11:16","s4-du-r2-10-v-t-10:26","s2-du-r2-14-v-t-13:09","s1-da-p2-8-o-d-12:04","s2-da-p2-12-v-d-11:57")


#Any g.fact >4 to inspect?
table_co2 %>% filter(!UniqueID%in%co2_discard) %>%
  filter(abs(LM.flux)>MDF.lim) %>% 
  filter(g.fact>=4) %>% 
  filter(!UniqueID%in%co2_gfact_inspected)  %>% 
  select(UniqueID, g.fact, LM.flux, HM.flux)






#Table TO decide----
#After checking that g.fact and mdf criteria are appropriate, subset dataset that needs a decision.


co2_todecide<- table_co2 %>% 
  #Remove incubations marked to discard
  filter(!UniqueID%in%co2_discard) %>%
  #Remove HM NA (when HM fails, there is nothing to decide, default to LM)
  filter(!is.na(HM.flux)) %>% 
  #Remove K.ratio >= 1 (HM curvature cannot exceed the max in any case, if so default to LM)
  filter(HM.k<k.max) %>% 
  #Remove g.fact >4 (after our inspection, these will default to LM)
  filter(g.fact<4) %>% 
  #Remove MDF fluxes (after our inspection, these will default to LM)
  filter(abs(LM.flux)>MDF.lim)


#Inspect AICc for selection criteria (using AICc weight)
#Test criteria: using AICc weight (based on the relative difference in AICc, calculate the AICc weight metric, a probabilistic estimate of how likely it is that one model is better than the other). 

co2_todecide <- co2_todecide %>%
  mutate(
    delta_AICc_LM = LM.AICc - pmin(LM.AICc, HM.AICc),  # Compare to the best AICc (min AICc)
    delta_AICc_HM = HM.AICc - pmin(LM.AICc, HM.AICc),
    
    # Calculate the Akaike weight for HM: AICc_weight_HM (0-1), how likely is it that HM is better than LM?
    AICc_weight_HM = exp(-0.5 * delta_AICc_HM) / (exp(-0.5 * delta_AICc_LM) + exp(-0.5 * delta_AICc_HM))
  ) %>% 
  mutate(k.ratio=HM.k/k.max)

#Inspect AICc_weight_HM
co2_todecide %>%
  ggplot(aes(x=AICc_weight_HM, fill=AICc_weight_HM>0.5))+
  geom_histogram()+
  scale_x_continuous(name="AICc weight for HM, probability of HM being best")

#Most cases have very clear separation. 


#Inspect against g.fact 
co2_todecide %>%
  ggplot(aes(x=AICc_weight_HM,y=g.fact, col=AICc_weight_HM>0.5))+
  geom_point()+
  geom_hline(yintercept=1.25)+
  scale_x_continuous(name="AICc weight for HM, probability of HM being best")

#Inspect cases with non-clear separation (those within 0.25-0.9 AICc_weight_HM) and with relatively big differences in flux estimate (g.fact>1.25)

co2_inspect<- co2_todecide %>% 
  filter(between(AICc_weight_HM, 0.25,0.9)) %>% 
  filter(g.fact>1.25)

#The following cases have been inspected and the model choice resulting from the above criteria (based on AICc weight) has been deemed aproppriate
co2_inspected_aicc<- c("s1-da-p1-11-o-d-10:46","s1-du-r2-5-b-d-08:32","s1-du-r2-8-v-d-09:51","s1-ri-a2-1-b-d-10:20","s2-ca-r2-3-o-d-09:09","s2-cu-a1-14-o-d-10:48","s2-cu-a1-4-v-d-09:05","s2-cu-a1-7-o-d-09:37","s2-cu-a2-2-v-d-08:32","s2-cu-a2-3-v-t-08:47","s2-da-r1-11-o-d-11:15","s2-da-r1-9-v-d-10:56","s2-du-a2-10-b-d-12:03","s2-ri-a2-10-b-d-10:48","s2-ri-p1-6-v-t-10:17","s2-ri-p2-11-v-d-12:24","s2-va-a2-14-v-t-12:57","s2-va-r1-9-v-d-11:42","s2-va-r2-14-o-d-12:59","s2-va-r2-7-v-t-11:08","s3-ca-p2-1-v-t-07:40","s3-ca-p2-5-v-t-09:33","s3-da-a1-2-b-d-07:34","s3-da-a2-10-o-d-08:57","s3-ri-a2-2-b-d-09:06","s4-da-a1-8-v-t-09:40","s4-da-r1-2-v-d-09:03","s4-da-r2-8-o-d-08:12","s4-ri-a2-12-b-d-09:02","s4-ri-a2-9-b-d-08:27","s2-va-r2-3-v-t-10:19")

co2_inspect %>% select(UniqueID, LM.flux, HM.flux, g.fact, k.ratio,AICc_weight_HM) %>% filter(!UniqueID%in%co2_inspected_aicc)


#Final step CO2: 
#Check that AICc weight and AICc raw comparison output the same choice (to be able to simplify the criteria), then limit the use of AICc weight to inspection of doubts (already done)
co2_todecide %>% 
  mutate(aicc_weight_choice=if_else(AICc_weight_HM>0.5,"HM","LM"),
         rawaicc_choice=if_else(LM.AICc>HM.AICc, "HM","LM")) %>% 
  mutate(choice_agree=aicc_weight_choice==rawaicc_choice) %>% 
  select(choice_agree) %>% distinct()



#FINAL CO2 table----
co2_final<-table_co2 %>% 
#Set prefered_model variable with current criteria (thresholds for MDF, Kappa and g.fact, plus best AICc as only criteria)
  mutate(prefered_model=case_when(UniqueID%in%co2_discard~"discard",
                                  is.na(HM.flux)~"LM",
                                  abs(LM.flux)<=MDF.lim~"LM",
                                  HM.k>=k.max~"LM",
                                  g.fact>=4~"LM",
                                  LM.AICc<HM.AICc~"LM",
                                  LM.AICc>HM.AICc~"HM",
                                  TRUE~"no selection"))

#Any flux without selection
unique(co2_final$prefered_model)

#TO DO-------
#Add quality_check variable to overwrite go.flux ones: "MDF", "HM kappa exceeds max", "G.factor exceeds 4", "Unreliable flux due to artefacts". Potentially multiple messages. 




#CH4 criteria--------

#before deciding on criteria, adapt new method CAmille and check results. 










#Weighted Score of models-----


## ---- Define your weights ----
weights <- c(RMSE = 0.15, MAE = 0.15, CV = 0.1, AICc_weight = 0.6)  # Adjust as needed
epsilon <- 1e-6  # Small constant to avoid division by zero

## ---- Compute Akaike weights ----
# Calculate the delta AICc for each model
#It provides the probability of being the best model (0-1) 
# Akaike weights give you a probabilistic interpretation of how much better one model is than the other

co2_todecide <- co2_todecide %>%
  mutate(
    delta_AICc_LM = LM.AICc - pmin(LM.AICc, HM.AICc),  # Compare to the best AICc (min AICc)
    delta_AICc_HM = HM.AICc - pmin(LM.AICc, HM.AICc),
    
    # Calculate the Akaike weights
    AICc_weight_LM = exp(-0.5 * delta_AICc_LM) / (exp(-0.5 * delta_AICc_LM) + exp(-0.5 * delta_AICc_HM)),
    AICc_weight_HM = exp(-0.5 * delta_AICc_HM) / (exp(-0.5 * delta_AICc_LM) + exp(-0.5 * delta_AICc_HM))
  )

## ---- Rel. improve and Score ----
co2_todecide <- co2_todecide %>%
  mutate(
    # Compute CVs
    LM.CV = LM.SE / (abs(LM.slope) + epsilon),
    HM.CV = HM.SE / (abs(HM.slope) + epsilon),
    
    # Compute relative improvements (positive = HM is better, relative reduction in error)
    Improvement.RMSE = (LM.RMSE - HM.RMSE) / LM.RMSE, 
    Improvement.MAE  = (LM.MAE  - HM.MAE)  / LM.MAE,  
    Improvement.CV   = (LM.CV   - (HM.CV/4))   / LM.CV,   
    Improvement.AICc = (LM.AICc - HM.AICc) / LM.AICc,
    # Rescale AICc weights to [-1 LM is 100% likely to 1 HM is 100% likely]
    Rescaled_AICc_weight_HM = 2 * AICc_weight_HM - 1,   # Rescale to [-1, 1]
    
    # Composite score using Akaike weights for AICc
    Composite.Score = Rescaled_AICc_weight_HM * weights["AICc_weight"] +
      Improvement.RMSE * weights["RMSE"] +
      Improvement.MAE  * weights["MAE"] +
      Improvement.CV   * weights["CV"],
    
    # Model decision
    Preferred.Model = if_else(Composite.Score > 0, "HM", "LM")
  )


#G.fact inspection----
#Not adapted

#lets see what is the highest g.fact of the best-performing HM models (best 75%)
good_hm<- co2_sig %>% 
  filter(HM.CV<=quantile(HM.CV, 0.8, na.rm=T))

ggplot(good_hm, aes(x=HM.CV, y=g.fact))+
  geom_point()+
  geom_hline(yintercept = 3)

good_hm%>% filter(g.fact>3) %>%arrange(HM.CV) %>%  select(UniqueID, HM.CV, HM.k, k.max, g.fact) %>% mutate(kratio=HM.k/k.max)

ggplot(good_hm, aes(x=HM.CV, y=HM.k/k.max))+
  geom_point()+
  geom_hline(yintercept = 1)+
  geom_label(data=. %>% filter(HM.k/k.max>0.5), aes(label=UniqueID))

co2_sig %>% 
  ggplot(aes(x=HM.CV, y=g.fact))+
  geom_point()+
  geom_hline(yintercept = 1)+
  geom_label(data=. %>% filter(HM.k/k.max>0.5), aes(label=UniqueID))





#SE.rel-----
#We need a metric to distinguish how good the estimates are when the model performance is very bad (even when the estimate is higher than the MDF), this will flag very noisy measurements, we could base this on se.rel (expressed as % of SE/flux). Hard limit could be se.rel< 5%, check number of incubations. 


co2_todecide %>% 
  ggplot(aes(x=abs(LM.se.rel), y=abs(HM.se.rel), col=Preferred.Model))+
  geom_point()+
  geom_abline(slope = 4)


#Proportion of incubations with less than 5% error for LM
n<- dim(co2_todecide)[1]
n_less_5perc <- length(which(abs(co2_todecide$LM.se.rel)<5))
ggplot(co2_todecide[order(abs(co2_todecide$LM.se.rel)),], aes(seq(1,n)/n*100, abs(LM.se.rel)))+
  geom_vline(xintercept = seq(0,100,50), color = "grey70")+
  geom_hline(yintercept = c(0,5,100), color = "grey70")+
  geom_point()+
  xlab("Proportion of timeseries")+
  ylab("Relative error [% of LM flux]")+
  theme_article()+
  ggtitle(paste0("LM SE.rel is below 5% for ",round(n_less_5perc/n*100,digits = 2),"% of the measurements"))


noisyLM<- co2_todecide[which(abs(co2_todecide$LM.se.rel)>=5),]

noisyLM %>% arrange(desc(abs(LM.flux))) %>% select(UniqueID, LM.flux, LM.se.rel, HM.se.rel, g.fact)

noisyLM %>% 
  ggplot(aes(x=abs(LM.flux), y=abs(LM.se.rel), col=Preferred.Model))+
  geom_point()


#Proportion of incubations with less than 5% error for HM
n_less_5perc_hm <- length(which(abs(co2_todecide$HM.se.rel)<5))
ggplot(co2_todecide[order(abs(co2_todecide$HM.se.rel)),], aes(seq(1,n)/n*100, abs(HM.se.rel)))+
  geom_vline(xintercept = seq(0,100,50), color = "grey70")+
  geom_hline(yintercept = c(0,5,100), color = "grey70")+
  geom_point()+
  xlab("Proportion of timeseries")+
  ylab("Relative error [% of HM flux]")+
  theme_article()+
  ggtitle(paste0("HM SE.rel is below 5% for ",round(n_less_5perc_hm/n*100,digits = 2),"% of the measurements"))

noisyHM<- co2_todecide[which(abs(co2_todecide$HM.se.rel)>=5),]
noisyHM %>% arrange(desc(abs(HM.flux))) %>% select(UniqueID, HM.flux, HM.se.rel, LM.se.rel, g.fact)

noisyHM %>% 
  ggplot(aes(x=abs(HM.flux), y=abs(HM.se.rel), col=Preferred.Model))+
  geom_point()




#Inspect composite score------
co2_todecide %>%
  ggplot(aes(x = Composite.Score, fill = Preferred.Model)) +
  geom_histogram(bins = 60, alpha = 0.3) +
  theme_minimal() +
  labs(title = "Composite Model Preference by Timeseries",
       x = "Composite Score (Positive = HM Better)",
       y = "Count")


co2_todecide %>%
  ggplot(aes(x=Rescaled_AICc_weight_HM, fill=Preferred.Model))+
  geom_histogram()+
  scale_x_continuous(name="AICc weight for HM, probability of HM being best")


co2_todecide %>%
  ggplot(aes(x=AICc_weight_HM, fill=Preferred.Model))+
  geom_histogram()+
  scale_x_continuous(name="AICc weight for HM, probability of HM being best")


co2_todecide %>%
  ggplot(aes(x=AICc_weight_HM, y=g.fact, col=Rescaled_AICc_weight_HM>0.05))+
  geom_point()+
  scale_x_continuous(name="AICc weight for HM, probability of HM being best")




co2_todecide %>%
  filter(Improvement.AICc>-0.2) %>% 
  filter(Improvement.AICc<0.5) %>% 
  ggplot(aes(x=Improvement.AICc, fill=Preferred.Model))+
  geom_histogram()+
  geom_vline(xintercept = 0.05)+
  scale_x_continuous(name="AICc relative improvement with HM")


co2_todecide %>% 
  ggplot(aes(x=Improvement.MAE, fill=Preferred.Model))+
  geom_histogram()+
  scale_x_continuous(name="Relative improvement of MAE with HM")

co2_todecide %>% 
  ggplot(aes(x=Improvement.RMSE, fill=Preferred.Model))+
  geom_histogram()+
  scale_x_continuous(name="Relative improvement of RMSE with HM")

co2_todecide %>% 
  ggplot(aes(x=Improvement.CV, fill=Preferred.Model))+
  geom_histogram(bins = 60)+
  scale_x_continuous(name="Relative improvement of SE.rel with HM")

co2_todecide %>% filter(Improvement.CV<(-1))

#DIrect comparison of model quality LM vs HM
#In current run, criteria are "SE", "RMSE", "AICc", 


co2_todecide %>% 
  ggplot(aes(x=LM.AICc, y=HM.AICc, col=Preferred.Model))+
  geom_point()+
  geom_abline(slope = 1)

co2_todecide %>% 
  ggplot(aes(x=LM.RMSE, y=HM.RMSE, col=Preferred.Model))+
  geom_point()+
  geom_abline(slope = 1)

co2_todecide %>% 
  ggplot(aes(x=LM.SE, y=HM.SE, col=Preferred.Model))+
  geom_point()+
  geom_abline(slope = 4)

#IN incubations where LM behaves good (90% of incubations with LM relative error <5%) and same flux estimate for LM and HM (0.95<g.fact<1.05), the relative error of the HM model is exactly 4 times as high as that of the LM. We should take this "penalization of the HM model into account in our weighted Improvement criteria. 
ggMarginal(p=
             co2_todecide %>% 
  filter(abs(LM.se.rel)<5) %>% 
  ggplot(aes(x=abs(LM.CV), y=abs(HM.CV), col=between(g.fact,0.95,1.05)))+
  geom_point()+
  geom_abline(slope = 4)
)

#Our best.flux criteria result in almost always choosing HM due to RMSE and AICc despite SE (HM wins 2:1)


#RMSE and MAE give the exact same info (rmse should be a bit more sensitive to extreme residuals)
co2_sig %>% 
  ggplot(aes(x=LM.MAE, y=HM.MAE, col=Preferred.Model))+
  geom_point()+
  geom_abline(slope = 1)

co2_sig %>% 
  ggplot(aes(x=LM.RMSE, y=HM.RMSE, col=Preferred.Model))+
  geom_point()+
  geom_abline(slope = 1)


co2_sig %>% 
  ggplot(aes(x=LM.CV, y=HM.CV, col=Preferred.Model))+
  geom_point()+
  geom_abline(slope = 1)



#Check Relative improvement of fit vs risk of oversestimation (g.fact)

co2_todecide %>% 
  ggplot(aes(x=(LM.RMSE-HM.RMSE)/LM.RMSE, y=g.fact, col=Preferred.Model))+
  geom_point()+
  scale_x_continuous(name="Relative improvement of RMSE with HM")+
  scale_y_log10()

co2_todecide %>% 
  ggplot(aes(x=(LM.MAE-HM.MAE)/LM.MAE, y=g.fact, col=Preferred.Model))+
  geom_point()+
  scale_x_continuous(name="Relative improvement of MAE with HM")

co2_todecide %>% 
  ggplot(aes(x=(LM.CV-HM.CV)/LM.CV, y=g.fact, col=Preferred.Model))+
  geom_point()+
  scale_x_continuous(name="Relative improvement of CV with HM")

co2_todecide %>% 
  ggplot(aes(x=AICc_weight_HM, y=g.fact, col=Preferred.Model))+
  geom_point()+
  scale_x_continuous(name="AICc weight for HM, probability of HM being best")+
  scale_y_log10()




co2_sig %>% 
  filter(LM.r2>0.15) %>% 
  filter(g.fact<=7.9695) %>% 
  ggplot(aes(x=(HM.r2-LM.r2)/LM.r2, y=g.fact,col=Preferred.Model))+
  geom_point()+
  geom_label(data=. %>% filter(g.fact>4.3)
               ,aes(label=UniqueID))+
  scale_x_continuous(name="Relative Improvement of R2 with HM")





co2_sig %>% 
  filter(LM.r2>0.15) %>% 
  filter(g.fact<=7.9695) %>% 
  arrange(desc(g.fact)) %>% select(UniqueID, g.fact) %>% slice(seq(1:20)) %>% pull(UniqueID)


#s4-du-r2-11-v-t-10:49 BAD FIT, weird artefact smooth DOUBT
#s1-ca-p1-14-v-d-14:49 BAD fit, super noisy, should be linear (arbitrary HM)
# s1-da-r2-12-v-t-11:31  cropped 2nd half, take linear estimate
#S4-cu-r1-9-v-t Noisy, maybe advance start to get better results### DOUBT, test
#"s2-cu-p1-13-o-d-10:11" Should be linear, smooth artefact creates bump for HM to win




#comparison between RMSE and MAE to dectect extreme outliers
co2_sig %>% 
  ggplot(aes(x=LM.RMSE, y=LM.MAE, col=((LM.RMSE-1)>1.2*LM.MAE)))+
  geom_point()+
  geom_abline(intercept=-0,slope = 0.85)
  
#Samples with high-leverage outliers: 
co2_sig  %>% filter((LM.RMSE-1)>1.2*LM.MAE) %>% arrange(desc(LM.RMSE/LM.MAE)) %>% pull(UniqueID)

#"s2-cu-a1-6-o-d-09:23": spikes and point-artefacts, keep like this
# "s2-cu-a1-7-o-d-09:37": keep like this
#"s2-cu-p1-5-o-d-08:30": point artefacts, keep like this
#"s2-va-a2-2-o-d-09:29": co2 spike mid-incubation, cannot crop


#Same exercise with HM:
co2_sig %>% 
  ggplot(aes(x=HM.RMSE, y=HM.MAE, col=(HM.RMSE-1)>1.2*HM.MAE))+
  geom_point()+
  geom_abline(intercept=-0,slope = 0.85)

co2_sig %>% 
  filter(!(LM.RMSE>1.5*LM.MAE+1.5|LM.RMSE>2*LM.MAE)) %>% 
  ggplot(aes(x=HM.RMSE,y=HM.RMSE/HM.MAE))+
  geom_point()+
  geom_label(data = . %>% filter(HM.RMSE>2*HM.MAE), aes(label=UniqueID))

#Fixed what could be fixed for cropping decissions based on HM mae vs rmse






co2_sig %>% 
  ggplot(aes(x=LM.SE, y=HM.SE, col=HM.SE>LM.SE))+
  geom_point()+
  geom_abline(slope = 1)

ggMarginal(p = co2_sig %>% 
             ggplot()+
  geom_point(aes(x=abs(LM.se.rel), y=abs(HM.se.rel), col=abs(HM.se.rel)>abs(LM.se.rel)))+
  geom_abline(slope = 1, intercept=0.6)+
    scale_x_log10()+
    scale_y_log10()
)

ggMarginal(p = co2_sig %>% 
  ggplot(aes(x=LM.r2, y=HM.r2, col=abs(LM.se.rel)>1))+
  geom_point()+
  geom_abline(slope = 1)+
  geom_hline(yintercept = 0.9)

)



ggMarginal(p = co2_sig %>% 
             ggplot(aes(x=LM.r2, y=abs(LM.se.rel), col=abs(LM.se.rel)>1))+
             geom_point()+
             scale_y_log10()

)



co2_sig %>% 
  ggplot(aes(x=abs(LM.flux), y=abs(HM.flux), col=LM.r2>HM.r2))+
  geom_point()+
  geom_abline(slope = 1)+
  scale_x_log10()+
  scale_y_log10()

co2_sig %>% 
  ggplot(aes(x=HM.RMSE/LM.RMSE, y=HM.RMSE, col=HM.r2>0.95))+
  geom_point()+
  facet_wrap(~model)


co2_sig %>% 
  filter(g.fact<20) %>% 
  ggplot()+
  geom_point(aes(x=HM.r2,y=g.fact, col=model))+
  geom_density(aes(x=HM.r2))+
  geom_label(data=. %>% filter(g.fact>5),aes(x=HM.r2,y=g.fact,label=UniqueID))


co2_sig %>% 
  ggplot(aes(y=g.fact, x=HM.RMSE/LM.RMSE))+
  geom_point()+
  geom_label(data=. %>% filter(g.fact>20),aes(label=UniqueID))


co2_sig %>% 
  filter(g.fact>1.40) %>% 
  ggplot(aes(x=HM.r2, y=LM.r2, col=g.fact))+
  geom_point()


#What is the maximum HM.k for very good models?

co2_sig %>% 
  filter(HM.r2>0.95) %>% 
  ggplot(aes(x=g.fact, y=HM.k/k.max))+
  geom_point()+
  geom_label(data=. %>% filter(HM.k/k.max>0.35), aes(label=UniqueID))


co2_todecide %>% 
  mutate(below_detectionLM=abs(LM.flux)<MDF.lim,
         below_detectionHM=abs(HM.flux)<MDF.lim) %>% 
  ggplot(aes(x=(HM.flux-LM.flux)/HM.flux, y=LM.r2, col=Preferred.Model))+
  geom_point()



table_co2 %>% 
  mutate(below_detectionLM=abs(LM.flux)<MDF.lim,
         below_detectionHM=abs(HM.flux)<MDF.lim) %>% 
  summarise(prop_BDL_LM=sum(below_detectionLM)/n(),
            prop_BDL_HM=sum(below_detectionHM, na.rm = T)/sum(!is.na(HM.flux)))



table_co2 %>% 
  ggplot(aes(x=(HM.flux-LM.flux)/HM.flux, y=(HM.r2-LM.r2)))+
  geom_point()


# ---- quality of model fits ----

n <- dim(table_co2)[1]
ind_lm <- which(table_co2$model=="LM")
in_hm <- which(table_co2$model=="HM")

message("Fluxes calculated for ",n, " different incubations.")
message(" --- ",dim(table_results_all[table_results_all$strata=="bare",])[1] ," in bare ")
message(" --- ",dim(table_results_all[table_results_all$strata=="open water",])[1] ," in open water ")
message(" --- ",dim(table_results_all[table_results_all$strata=="vegetated",])[1] ," in vegetated")


message("Best model is non-linear for ",round(length(in_hm)/n*100*100)/100,"% of the measurements")


ind_flagged <- which(table_co2$quality.check!="")
# table_co2$quality.check[ind_flagged]
message(round(length(ind_flagged)/n*100*100)/100,"% of the measurements are flagged")


in_hm_notflagged <- which(table_co2$model=="HM" & table_co2$quality.check=="")
message("Best model is non-linear for ",round(length(in_hm_notflagged)/(n-length(ind_flagged))*100*100)/100,"% of the measurements not flagged")

table.flags <- NULL
for(flag in unique(table_co2$quality.check[ind_flagged])){
  ind_mdf <- which(table_co2$quality.check==flag)
  message(round(length(ind_mdf)/n*100*100)/100,"% of the measurements are flagged because of ",flag)
  
  table.flags <- rbind(table.flags,
                       data.frame(flag = flag,
                                  n = length(ind_mdf)))
}
table.flags$perc <- round(table.flags$n/n*100*100)/100
table.flags <- table.flags[order(table.flags$n, decreasing = T),]
table.flags



table_co2[which(is.na(table_co2$quality.check)),]

ggplot(co2_sig[which(!is.na(co2_sig$quality.check)),], 
       aes(quality.check, (HM.flux-LM.flux)/HM.flux, fill = model))+
  geom_hline(yintercept = c(-1,0,1), color = "grey70")+
  geom_jitter(alpha=0.1, size=2, width = 0.1)+
  geom_violin(alpha=0.2, scale = "width", draw_quantiles = c(0.5), aes(colour = model, fill = model))+
  theme_article()+
  # scale_x_log10()+
  # scale_y_log10()+
  ylim(c(-1,2))+
  scale_fill_viridis_d(option = "A", end = 0.8, direction = -1)+
  scale_colour_viridis_d(option = "A", end = 0.8, direction = -1)+
  coord_flip()

table_co2$diff_model <- (table_co2$HM.flux-table_co2$LM.flux)/table_co2$HM.flux

n_less_10perc <- length(which(abs(table_co2$diff_model)<0.1))

message("Linear and non-linear models have less than 10% difference for ",round(n_less_10perc/n*100*100)/100,"% of the measurements")

ggplot(table_co2[order((table_co2$diff_model)),], aes(seq(1,n)/n*100, (diff_model)*100))+
  geom_vline(xintercept = seq(0,100,50), color = "grey70")+
  geom_hline(yintercept = c(-100,-10,0,10,100), color = "grey70")+
  geom_point()+
  # scale_y_log10()+
  ylim(c(-50,120))+
  xlab("Proportion of timeseries")+
  ylab("Relative difference [% of HM flux]")+
  theme_article()+ggtitle(paste0("HM-LM models difference is below 10% for ",round(n_less_10perc/n*100*10)/10,"% of the measurements"))


ggplot(table_co2, aes(abs(best.flux), abs(diff_model)*100))+
  # geom_vline(xintercept = seq(0,100,50), color = "grey70")+
  geom_hline(yintercept = c(1,10,100), color = "grey70")+
  geom_point(aes(colour = quality.check))+
  scale_x_log10()+
  # xlim(c(-50,50))+
  # ylim(c(-50,120))+
  # xlab("Proportion of timeseries")+
  ylab("Relative difference [% of HM flux]")+
  theme_article()+ggtitle(paste0("HM-LM models difference is below 10% for ",round(n_less_1perc/n*100*10)/10,"% of the measurements"))




ggplot(co2_sig)+
  geom_point(aes(abs(best.flux), HM.MAE, colour = "HM"), alpha=0.2, size=3)+
  geom_point(aes(abs(best.flux), LM.MAE, colour = "LM"), alpha=0.2, size=3)+
  theme_article()+
  # xlim(c(0,50))+
  scale_x_log10()+
  scale_y_log10()+
  ylab("MAE")+
  scale_colour_viridis_d("model", option = "A", end = 0.8, direction = -1)


ggplot(co2_sig, aes(LM.MAE, HM.MAE, shape=model,colour = model))+
  geom_point(alpha=0.4, size=3)+
  theme_article()+
  scale_x_log10()+
  scale_y_log10()+
  scale_colour_viridis_d("best model", option = "A", end = 0.8, direction = -1)



#Check which criteria is the "limiting factor" in the selection for each incubation. 


# ---- quality of model fits ----
# c("MAE", "RMSE", "AICc", "SE") allowed in best.flux

#In current run, criteria are "SE", "RMSE", "AICc", "kappa" ("kappa" included but non-working, due to k.HM already being limited to the max allowed in criteria)



#1. For fluxes below detection (MDF), the choice of model is not relevant (we can leave the best.flux result). 
table_co2 %>% 
  mutate(below_detection=abs(best.flux)<MDF.lim) %>%  
  ggplot(aes(x=LM.r2, fill=below_detection))+
  geom_histogram()+
  facet_wrap(~below_detection, scales="free")



#2. For fluxes with poor fit (R2< threshold), LM should be assumed to avoid noise creating artificially high fluxes.

table_co2 %>% 
  mutate(below_detection=abs(best.flux)<MDF.lim) %>%  
  filter(!below_detection) %>% 
  filter(LM.r2<0.5) %>% 
  ggplot(aes(x=HM.r2, y=LM.r2, col=model))+
  geom_point()


#There are negative r2!! 
table_co2 %>% 
  filter(LM.r2<0|HM.r2<0) %>% 
  select(UniqueID, quality.check, model, LM.r2, HM.r2, MDF.lim, best.flux) %>% 
  mutate(below_detection=abs(best.flux)<MDF.lim) %>%  
  arrange(below_detection)


#General considerations: HM has the risk of overestimating fluxes, there should be a stricter criteria for chosing HM over LM when the difference (g.fact) is very large. Additionally, when fluxes are very low, the fit of the models is going to be worse, and the risk of overestimation is larger. 








#SElect the "worst-HM" and inspect visually, 
#inspect highest difference (g-fact) to see the adequacy. We should only go with HM when it fits really-really well. When the fit is not very good, it will be more prudent to trust the LM. 

table_co2 %>% 
  filter(!grepl("MDF",quality.check)) %>% 
  filter(!is.na(g.fact)) %>%
  ggplot(aes(y=g.fact, fill=model))+
  geom_histogram(bins = 100)+
  geom_boxplot()+
  facet_wrap(~model)




table_co2 %>% 
  ggplot(aes(x=g.fact, y=abs(HM.se.rel), col=model))+
  geom_point()+
  scale_x_log10()

table_co2 %>% 
  ggplot(aes(x=abs(g.fact), y=HM.r2, col=model))+
  geom_point()+
  scale_x_log10()

table_co2 %>% 
  select(UniqueID,g.fact, HM.se.rel, HM.SE,HM.r2, LM.se.rel,LM.SE,LM.r2, model) %>% 
  arrange(desc(g.fact))










table_co2 %>% 
  filter(UniqueID!="s1-ri-a2-8-v-t-11:24") %>% 
  ggplot(aes(x=LM.r2, y=HM.r2, col=model))+
  geom_point()

table_co2 %>% 
  filter(UniqueID!="s1-ri-a2-8-v-t-11:24") %>% 
  ggplot(aes(x=abs(LM.se.rel), y=abs(HM.se.rel), col=model))+
  geom_point()


table_co2 %>% 
  filter(HM.r2>0.2&LM.r2>0.2) %>% 
  filter(UniqueID!="s1-ri-a2-8-v-t-11:24") %>% 
  ggplot(aes(x=(HM.r2/LM.r2), y=HM.k, col=model))+
  geom_point()




# ---- g.factor ------

table_co2 %>% 
  group_by(model) %>% 
  filter(!is.na(g.fact)) %>% 
  # filter(UniqueID!="s1-ri-a2-8-v-t-11:24") %>% #this sample has HM= 0 unclear why, but method fails
  summarise(maxg=max(g.fact),
            ming=min(g.fact),
            avg.g=mean(g.fact))

table_co2 %>% 
  filter(g.fact==min(g.fact,na.rm = T)) %>% pull(UniqueID)

table_ch4 %>% 
  group_by(model) %>% 
  filter(!is.na(g.fact)) %>% 
  # filter(UniqueID!="s1-ri-a2-8-v-t-11:24") %>% #this sample has HM= 0 unclear why, but method fails
  summarise(maxg=max(g.fact),
            ming=min(g.fact),
            avg.g=mean(g.fact))

table_co2 %>% 
  filter(g.fact==min(g.fact,na.rm = T)) %>% pull(UniqueID)





table_co2 %>% 
  filter(!is.na(g.fact)) %>% 
  filter(UniqueID!="s1-ri-a2-8-v-t-11:24") %>%
  ggplot(aes(y=g.fact, col=model))+
  geom_boxplot()+
  scale_y_log10()


# ---- Kappa ------
#K ratio is limited to 1 by our goflux call, and in best.flux we set the limit to 1 so, the kappa criteria is never used to select the LM over the HM 


#Check the K.ratio between model fit and theoretical maximum:
table_co2 %>% 
  ggplot(aes(x=HM.k/k.max, fill=model))+
  geom_histogram(position = "identity")+
  facet_wrap(~model)

table_co2 %>% 
  ggplot(aes(x=HM.k/k.max, y=g.fact, col=model))+
geom_point()


#k-ratio of 1 can represent: noisy measurements in which HM is not adequate or true "exagerated curvatures" when g-factor is very high (i.e. when the HM flux is very very high and the curvature has been limited to its maximum)
table_co2 %>% 
  mutate(kratio=HM.k/k.max) %>% 
  filter(kratio==max(kratio, na.rm=T)) %>% 
  select(UniqueID, HM.RMSE, LM.RMSE, HM.k, g.fact) %>% 
  filter(g.fact>4) %>% 
  arrange(g.fact)
  
# the samples with kratio==1, 

#incubations with exagerated curvature (kappa==max.kappa)
table_co2 %>% 
  mutate(kratio=HM.k/k.max) %>% 
  filter(kratio==1) %>% 
  ggplot(aes(x=g.fact, y=HM.AICc/LM.AICc, col=model))+
  geom_point()



# ---- quality of model fits ----




n <- dim(table_co2)[1]
ind_lm <- which(table_co2$model=="LM")
in_hm <- which(table_co2$model=="HM")

message("Fluxes calculated for ",n, " different incubations.")
message(" --- ",dim(table_results_all[table_results_all$strata=="bare",])[1] ," in bare ")
message(" --- ",dim(table_results_all[table_results_all$strata=="open water",])[1] ," in open water ")
message(" --- ",dim(table_results_all[table_results_all$strata=="vegetated",])[1] ," in vegetated")


message("Best model is non-linear for ",round(length(in_hm)/n*100*100)/100,"% of the measurements")


ind_flagged <- which(table_co2$quality.check!="")
# table_co2$quality.check[ind_flagged]
message(round(length(ind_flagged)/n*100*100)/100,"% of the measurements are flagged")


in_hm_notflagged <- which(table_co2$model=="HM" & table_co2$quality.check=="")
message("Best model is non-linear for ",round(length(in_hm_notflagged)/(n-length(ind_flagged))*100*100)/100,"% of the measurements not flagged")

table.flags_co2 <- NULL
for(flag in unique(table_co2$quality.check[ind_flagged])){
  ind_mdf <- which(table_co2$quality.check==flag)
  message(round(length(ind_mdf)/n*100*100)/100,"% of the measurements are flagged because of ",flag)
  
  table.flags_co2 <- rbind(table.flags_co2,
                           data.frame(flag = flag,
                                      n = length(ind_mdf)))
}
table.flags_co2$perc <- round(table.flags_co2$n/n*100*100)/100
table.flags_co2 <- table.flags_co2[order(table.flags_co2$n, decreasing = T),]
table.flags_co2


#co2 fluxes without any flag:
table_co2[which(is.na(table_co2$quality.check)|table_co2$quality.check==""),]

#Goflux meaning of different flags: 
#SE: "Standard Error" it tells us that the noise in our incubation is larger than the instrument accuracy (i.e. we have more noise than we would expect just from "instrumental noise").
#MDF: "minimum detectable flux". I.e. flux is Below detection limit. 
#gfact > 2: the ratio bettewwn the non-linear(HM) flux  estimate and the linear flux estimate is larger than 2: i.e. HM.flux > 2*LM.flux . IF this happens, we assume that the HM.flux is overestimated, so we take the LM.flux instead (more or less arbitrary decission, different thresholds can be selected (relaxed gfact>4, medium gfact>2, strict gfact>1.25))
#gfact < 0.5: same as above but for the case when LM.flux is higher than HM.flux
#HM.flux is NA: non-linear model failed to return a flux.



ind_flagged <- which(table_ch4$quality.check!="")
# table_co2$quality.check[ind_flagged]
message(round(length(ind_flagged)/n*100*100)/100,"% of the measurements are flagged")

table.flags_ch4 <- NULL
for(flag in unique(table_ch4$quality.check[ind_flagged])){
  ind_mdf <- which(table_ch4$quality.check==flag)
  message(round(length(ind_mdf)/n*100*100)/100,"% of the measurements are flagged because of ",flag)
  
  table.flags_ch4 <- rbind(table.flags_ch4,
                           data.frame(flag = flag,
                                      n = length(ind_mdf)))
}
table.flags_ch4$perc <- round(table.flags_ch4$n/n*100*100)/100
table.flags_ch4 <- table.flags_ch4[order(table.flags_ch4$n, decreasing = T),]
table.flags_ch4


#ch4 fluxes without any goflux-flag:
table_ch4[which(is.na(table_ch4$quality.check)|table_ch4$quality.check==""),]



#Examine R2 for CO2
ggplot(table_co2)+
  geom_histogram(aes(x=LM.r2, fill="linear"),bins = 40,alpha = 0.5)+
  geom_histogram(aes(x=HM.r2, fill="non-linear"),bins = 40,alpha=0.5)+
  labs(fill="R2")+
  geom_vline(xintercept = 0.9)+
  facet_wrap(~model)


#Examine R2 for CH4
ggplot(table_ch4)+
  geom_histogram(aes(x=LM.r2, fill="linear"),bins = 40,alpha = 0.5)+
  geom_histogram(aes(x=HM.r2, fill="non-linear"),bins = 40,alpha=0.5)+
  labs(fill="R2")+
  geom_vline(xintercept = 0.9)



#Many different cases, most likely "manual" cropping for incubations with less than R2=0.8 is recomended (this represents ~1000 inspections, doable semi-automatically). 
table_co2 %>% filter(HM.r2<0.8|LM.r2<0.8) %>% summarise(n=n()) %>% pull(n)
table_ch4 %>% filter(HM.r2<0.8|LM.r2<0.8) %>% summarise(n=n()) %>% pull(n)



#List of fit-Flagged####

#Get full table with fit qualities: 
#UniqueID, start.time_gofluxfit, duration, 
#CO2_best.model, CO2_quality.check, CO2_contains.artefact, CO2_best.model.R2,
#CH4_best.model, CH4_quality.check, CH4_contains.artefact, CH4_best.model.R2,

co2_quality<- table_co2 %>% 
  select(UniqueID, start.time,duration, 
         model, quality.check, contains.artefact, LM.r2, HM.r2) %>% 
  mutate(best.model.r2=case_when(model=="LM"~LM.r2,
                                 model=="HM"~HM.r2)) %>% 
  rename(start.time_gofluxfit=start.time, co2_best.model=model, co2_contains.artefact=contains.artefact, co2_best.model.r2=best.model.r2, co2_quality.check=quality.check) %>% select(-c(LM.r2,HM.r2))

ch4_quality<- table_ch4 %>% 
  select(UniqueID, start.time,duration, 
         model, quality.check, contains.artefact, LM.r2, HM.r2) %>% 
  mutate(best.model.r2=case_when(model=="LM"~LM.r2,
                                 model=="HM"~HM.r2)) %>% 
  rename(start.time_gofluxfit=start.time, ch4_best.model=model, ch4_contains.artefact=contains.artefact, ch4_best.model.r2=best.model.r2,ch4_quality.check=quality.check) %>% select(-c(LM.r2,HM.r2))


table_quality<- co2_quality %>% 
  merge.data.frame(ch4_quality, by=c("UniqueID", "start.time_gofluxfit","duration")) %>% 
  mutate(co2_to_inspect= (co2_contains.artefact!="no artefact") | (co2_best.model.r2<0.8),
         ch4_to_inspect= (ch4_contains.artefact!="no artefact") | (ch4_best.model.r2<0.8),
         any_to_inspect= co2_to_inspect|ch4_to_inspect)


#save table of incubations to inspect.
write.csv(table_quality, file = paste0(quality_path,"fit_models_flags_5smargin.csv"),row.names = F)


table_quality %>% 
  summarise(co2_flags=sum(co2_to_inspect),
            ch4_flags=sum(ch4_to_inspect),
            all_flags=sum(any_to_inspect))

# Fieldsheet notes####

#Compile all fieldsheets with comments
fieldsheets<- list.files(path= fieldsheetpath, recursive = T, pattern = "GHG.xlsx",full.names = T)
field<- read_GHG_fieldsheets(fieldsheets)

fnotes<- field %>% 
  filter(!is.na(comments))

write.csv(fnotes, file=paste0(quality_path,"fieldsheets_with_comments.csv"),row.names = F)

notes_unique <-fnotes %>% group_by(comments) %>% summarise(n=n())

#Manual step: inspect comments and add decissions
#DONE. 

#Import and Join inspectionlist-----

field_withcoments<- read.csv(paste0(quality_path,"decisions_from_fieldsheets_with_comments.csv"))

field_decissions<- field_withcoments %>% 
  filter(decission_for_incubation!="") %>%
  select(uniqID, start_time, end_time, decission_for_incubation, crop_start_s, crop_end_s, comments) %>% rename(UniqueID=uniqID, start_timefield=start_time, end_timefield=end_time, decissionfield=decission_for_incubation, commentfield=comments)


#Join and re-arrange order to check (we will do 1 pdf per season)
table_quality_field<- table_quality %>% 
  merge.data.frame(field_decissions, by="UniqueID",all = T) %>% 
  mutate(season=substr(UniqueID,1,2)) %>% 
  select(season,UniqueID,start.time_gofluxfit,duration,
         co2_best.model,co2_quality.check,co2_contains.artefact,co2_best.model.r2,
         ch4_best.model,ch4_quality.check,ch4_contains.artefact,ch4_best.model.r2,
         co2_to_inspect,ch4_to_inspect,any_to_inspect,
         start_timefield, end_timefield, decissionfield, crop_start_s, crop_end_s, commentfield) %>% 
  arrange(UniqueID, decissionfield, !any_to_inspect)

write.csv(table_quality_field, file=paste0(quality_path,"Inspection_table_allincubations_withfielddecissions_toFILL.csv"),row.names = F)



#MANUAL STEP: 

#Inspect all suspect incubations and decide cropping or discarding fluxes.

#All incubations inspected: crop and discard decissions noted in excel:
# Inspection_table_allincubations_tocrop.xlsx
#Inspection based on plots of fluxes produced by start-stop times of corrected fieldsheets plus a 5 second margin to start_time and stop_time: all start-stop times for final calculations have to include the 5 second margins plus (if aplicable) the cropping decission logged in excell.


#Common patterns during inspection: 
#LM CO2 underestimating photosynthetic flux, curving of HM too limited (Kappa). Inpect examples and set new criteria for best.flux g-factor limit (the maximum disparity allowed for HMflux/LMflux)


#Inspection plots-----

#Use the pdfs created by raw2flux_miguel_edits_5smargin to inspect and note weird things in the inspection_table_TOFILL

#What time exactly is in plots?
#Second 0 in plots is exactly start.time_gofluxfit (wich is 5s after the start-time from fieldsheet-corrections)

#Common patterns: 
#LM CO2 underestimating photosynthetic flux, curving of HM too limited (Kappa). Inpect examples and set new criteria for best.flux g-factor limit (the maximum disparity allowed for HMflux/LMflux)


#Miscelaneous plots####

#extract a few candidates for inspection based on R2:
#CO2 good fits r2>0.95: all looks good
table_co2 %>% 
  filter(between(x=HM.r2,
                 lower= quantile(.$HM.r2,probs = 0.5, na.rm = T),
                 upper =quantile(.$HM.r2,probs = 0.57, na.rm = T))) %>% 
  filter(model=="HM") %>% 
  select(UniqueID, best.flux, HM.r2)%>% 
  filter(best.flux%in%c(max(best.flux),min(best.flux),quantile(best.flux, 0.5)))

#CO2 good fits r2>0.90: A bit more noise, dubious minor outliers
table_co2 %>% 
  filter(between(x=HM.r2,
                 lower= quantile(.$HM.r2,probs = 0.4, na.rm = T),
                 upper =quantile(.$HM.r2,probs = 0.43, na.rm = T))) %>% 
  filter(model=="HM") %>% 
  select(UniqueID, best.flux, HM.r2)%>% 
  filter(best.flux%in%c(max(best.flux),min(best.flux),quantile(best.flux, 0.5)))

#CO2 acceptable fits r2>0.85: A bit more noise, overall good
table_co2 %>% 
  filter(between(x=HM.r2,
                 lower= quantile(.$HM.r2,probs = 0.35, na.rm = T),
                 upper =quantile(.$HM.r2,probs = 0.36, na.rm = T))) %>% 
  filter(model=="HM") %>% 
  select(UniqueID, best.flux, HM.r2)%>% 
  filter(best.flux%in%c(max(best.flux),min(best.flux),quantile(best.flux, 0.5)))

#CO2 acceptable fits r2>0.80: A bit more noise, overall good (reliable fluxes)
table_co2 %>% 
  filter(between(x=HM.r2,
                 lower= quantile(.$HM.r2,probs = 0.30, na.rm = T),
                 upper =quantile(.$HM.r2,probs = 0.32, na.rm = T))) %>% 
  filter(model=="HM") %>% 
  select(UniqueID, best.flux, HM.r2)%>% 
  filter(best.flux%in%c(max(best.flux),min(best.flux),quantile(best.flux, 0.5)))

#CO2 dubious fits r2>0.75: noisy but flux discernible, manipulation artefacts seen.
table_co2 %>% 
  filter(between(x=HM.r2,
                 lower= quantile(.$HM.r2,probs = 0.25, na.rm = T),
                 upper =quantile(.$HM.r2,probs = 0.26, na.rm = T))) %>% 
  filter(model=="HM") %>% 
  select(UniqueID, best.flux, HM.r2) %>% 
  filter(best.flux%in%c(max(best.flux),min(best.flux),quantile(best.flux, 0.5)))

#CO2 fits: r2 0.65: noisy but clear
table_co2 %>% 
  filter(between(x=HM.r2,
                 lower= quantile(.$HM.r2,probs = 0.2, na.rm = T),
                 upper =quantile(.$HM.r2,probs = 0.21, na.rm = T))) %>% 
  filter(model=="HM") %>% 
  select(UniqueID, best.flux, HM.r2) %>% 
  filter(best.flux%in%c(max(best.flux),min(best.flux),quantile(best.flux, 0.5)))
#Co2 fits r2 ~0.5: noise & a few outliers but clear trend
table_co2 %>% 
  filter(between(x=HM.r2,
                 lower= quantile(.$HM.r2,probs = 0.12, na.rm = T),
                 upper =quantile(.$HM.r2,probs = 0.18, na.rm = T))) %>% 
  filter(model=="HM") %>% 
  select(UniqueID, best.flux, HM.r2) %>% 
  filter(best.flux%in%c(max(best.flux),min(best.flux),quantile(best.flux, 0.5)))
#CO2 fits r2 ~0.15-0.3, Artefact instrument, nonconsistent trends, 
table_co2 %>% 
  filter(between(x=HM.r2,
                 lower= quantile(.$HM.r2,probs = 0.05, na.rm = T),
                 upper =quantile(.$HM.r2,probs = 0.11, na.rm = T))) %>% 
  filter(model=="HM") %>% 
  select(UniqueID, best.flux, HM.r2) %>% 
  filter(best.flux%in%c(max(best.flux),min(best.flux),quantile(best.flux, 0.5)))







ggplot(table_co2[which(!is.na(table_co2$quality.check)),], 
       aes(quality.check, (HM.flux-LM.flux)/HM.flux, fill = model))+
  geom_hline(yintercept = c(-1,0,1), color = "grey70")+
  geom_jitter(alpha=0.1, size=2, width = 0.1)+
  geom_violin(alpha=0.2, scale = "width", draw_quantiles = c(0.5), aes(colour = model, fill = model))+
  theme_article()+
  # scale_x_log10()+
  # scale_y_log10()+
  ylim(c(-1,2))+
  scale_fill_viridis_d(option = "A", end = 0.8, direction = -1)+
  scale_colour_viridis_d(option = "A", end = 0.8, direction = -1)+
  coord_flip()

table_co2$diff_model <- (table_co2$HM.flux-table_co2$LM.flux)/table_co2$HM.flux

n_less_1perc <- length(which(abs(table_co2$diff_model)<0.01))

message("Linear and non-linear models have less than 1% difference for ",round(n_less_1perc/n*100*100)/100,"% of the measurements")

ggplot(table_co2[order((table_co2$diff_model)),], aes(seq(1,n)/n*100, (diff_model)*100))+
  geom_vline(xintercept = seq(0,100,50), color = "grey70")+
  geom_hline(yintercept = c(-100,-10,0,10,100), color = "grey70")+
  geom_point()+
  # scale_y_log10()+
  ylim(c(-50,120))+
  xlab("Proportion of timeseries")+
  ylab("Relative difference [% of HM flux]")+
  theme_article()+ggtitle(paste0("HM-LM models difference is below 10% for ",round(n_less_1perc/n*100*10)/10,"% of the measurements"))


ggplot(table_co2, aes(abs(best.flux), abs(diff_model)*100))+
  # geom_vline(xintercept = seq(0,100,50), color = "grey70")+
  geom_hline(yintercept = c(1,10,100), color = "grey70")+
  geom_point(aes(colour = quality.check))+
  scale_x_log10()+
  scale_y_log10()+
  # xlim(c(-50,50))+
  # ylim(c(-50,120))+
  # xlab("Proportion of timeseries")+
  ylab("Relative difference [% of HM flux]")+
  theme_article()+ggtitle(paste0("HM-LM models difference is below 1% for ",round(n_less_1perc/n*100*10)/10,"% of the measurements"))


df_exceedance_co2 <- NULL
for(thresh in c(0.1,0.5,seq(1,100))){
  n_less <- length(which(abs(table_co2$diff_model)<thresh/100))
  df_exceedance_co2 <- rbind(df_exceedance_co2,
                             data.frame(t = thresh,
                                        n = n_less,
                                        p = n_less/dim(table_co2)[1]*100))
}

p_co2 <- ggplot(df_exceedance_co2, aes(t, p))+geom_path()+geom_point()+theme_article()+
  geom_segment(aes(x=-0,xend=10,
                   y=df_exceedance_co2$p[df_exceedance_co2$t==10], yend=df_exceedance_co2$p[df_exceedance_co2$t==10]))+
  geom_segment(aes(x=10,xend=10,
                   y=-Inf, yend=df_exceedance_co2$p[df_exceedance_co2$t==10]))+
  # scale_x_log10()+
  ylab("Proportion of timeseries [%]")+
  xlab("Relative difference [% of HM flux]")+
  ggtitle(paste0("For CO2, HM-LM models difference is below 10% for ",round(df_exceedance_co2$p[df_exceedance_co2$t==10]),"% of the measurements"))


table_ch4$diff_model <- (table_ch4$HM.flux-table_ch4$LM.flux)/table_ch4$HM.flux

df_exceedance_ch4 <- NULL
for(thresh in c(0.1,0.5,seq(1,100))){
  n_less <- length(which(abs(table_ch4$diff_model)<thresh/100))
  df_exceedance_ch4 <- rbind(df_exceedance_ch4,
                             data.frame(t = thresh,
                                        n = n_less,
                                        p = n_less/dim(table_ch4)[1]*100))
}

p_ch4 <- ggplot(df_exceedance_ch4, aes(t, p))+geom_path()+geom_point()+theme_article()+
  geom_segment(aes(x=-0,xend=10,
                   y=df_exceedance_ch4$p[df_exceedance_ch4$t==10], yend=df_exceedance_ch4$p[df_exceedance_ch4$t==10]))+
  geom_segment(aes(x=10,xend=10,
                   y=-Inf, yend=df_exceedance_ch4$p[df_exceedance_ch4$t==10]))+
  # scale_x_log10()+
  ylab("Proportion of timeseries [%]")+
  xlab("Relative difference [% of HM flux]")+
  ggtitle(paste0("For CH4, HM-LM models difference is below 10% for ",round(df_exceedance_ch4$p[df_exceedance_ch4$t==10]),"% of the measurements"))


p_diff_models <- ggarrange(p_co2, p_ch4, ncol = 1)
ggsave(plot = p_diff_models, filename = "HM-LM models difference.jpeg", path = plots_path, 
       width = 7, height = 5, dpi = 300, units = 'in', scale = 1.1)


ggplot(table_co2, aes(HM.MAE, diff_model))+geom_point(aes(colour = model))+
  scale_colour_viridis_d("best model", option = "A", end = 0.8, direction = -1)+
  theme_article()+ 
  scale_x_log10()+  scale_y_log10()


ggplot(table_ch4, aes(HM.MAE, diff_model))+geom_point(aes(colour = model))+
  scale_colour_viridis_d("best model", option = "A", end = 0.8, direction = -1)+
  theme_article()+ 
  scale_x_log10()+  scale_y_log10()




ggplot(table_co2)+
  geom_point(aes(abs(best.flux), HM.MAE, colour = "HM"), alpha=0.2, size=3)+
  geom_point(aes(abs(best.flux), LM.MAE, colour = "LM"), alpha=0.2, size=3)+
  theme_article()+
  # xlim(c(0,50))+
  scale_x_log10()+
  scale_y_log10()+
  ylab("MAE")+
  scale_colour_viridis_d("model", option = "A", end = 0.8, direction = -1)


ggplot(table_ch4, aes(LM.MAE, HM.MAE, shape=model,colour = model))+
  geom_point(alpha=0.4, size=3)+
  theme_article()+
  scale_x_log10()+
  scale_y_log10()+
  scale_colour_viridis_d("best model", option = "A", end = 0.8, direction = -1)





ggplot(table_ch4)+
  geom_density(aes(LM.MAE, fill = "LM"), alpha=0.5)+
  geom_density(aes(HM.MAE, fill = "HM"), alpha=0.5)+
  theme_article()+
  # facet_grid(.~strata)+
  scale_x_log10()


ggplot(table_ch4, aes(LM.MAE, best.flux, colour = HM.MAE))+
  geom_point(alpha=0.9, size=3)+
  theme_article()+
  scale_x_log10()+
  scale_y_log10()+
  scale_colour_viridis_c(option = "A", end = 0.8, direction = -1)





# ---- Select data with reasonable fits ----

ind_sel <- which(table_co2$HM.MAE<10)
table_co2_sel <- table_co2[ind_sel,]
table_co2_sel <- table_co2_sel[]

ind_sel <- which(table_ch4$HM.MAE<100)
table_ch4_sel <- table_ch4[ind_sel,]


# save selection in a single table

table_co2_out <- table_co2_sel[,c(which(names(table_co2_sel)=="UniqueID"),
                                  which(names(table_co2_sel)=="sampling"):which(names(table_co2_sel)=="lightCondition"),
                                  which(names(table_co2_sel)=="best.flux"))]

myfilename <- paste("all_CO2_fluxes",min(as.Date(table_co2_out$start.time)),"to",
                    max(as.Date(table_co2_out$start.time)), sep = "_")
write.csv(x = table_co2_out, file = paste0(myfilename,".csv"), 
          row.names = F)

table_ch4_out <- table_ch4_sel[,c(which(names(table_ch4_sel)=="UniqueID"),
                                  which(names(table_ch4_sel)=="sampling"):which(names(table_ch4_sel)=="lightCondition"),
                                  which(names(table_ch4_sel)=="best.flux"),
                                  which(names(table_ch4_sel)=="ebullition"),
                                  which(names(table_ch4_sel)=="diffusion"))]

myfilename <- paste("all_CH4_fluxes",min(as.Date(table_ch4_out$start.time)),"to",
                    max(as.Date(table_ch4_out$start.time)), sep = "_")
write.csv(x = table_ch4_out, file = paste0(myfilename,".csv"), 
          row.names = F)





