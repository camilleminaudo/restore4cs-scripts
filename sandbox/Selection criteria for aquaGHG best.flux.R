#Selection criteria for CH4 water (aquaGHG)

#Description: this script is used to flag errors in flux separation of CH4 water and to chose the best.flux in each case. 

#IMPORTANT: all goflux-outputs of aquaHG with flux-separation==T refer to cropped incubations based on diffusion detection algorithm. 



#Chosing best flux will follow this logic: 

#0. If incubation is wrong (disard due to artefacts or wrong manipulation), no flux estimate is provided as best NA and quality flag is added to discard. 

#1. If there is water, flux separator method will be calculated

#2. IF No bubles are found (i.e. length difusion = total length-1)--> GOflux criteria (to determine, similar logic to CO2 criteria, after inspection for hard thresholds for g.fact)

#3. IF Bubles are found (difusion length < total length - 1): 

  #IF ebullition is significant (based on flux +- SD of total flux, difusive flux, ebullitive flux): Use total flux as best flux 

  #IF ebullition is non-significant, or method fails (difusion > total flux):

      #Check Goflux LM quality, if above threshold (to determine), best flux is taken from LM
          #IF LM quality is below threshold, best flux is total flux 


#STEPS: 
  #1. First of all, total flux must be calculated for the whole period of the incubation (now it seems it is only calculated from first point of detected difusive chunk)

  #2. Inpsect scripts aquaGHG to determine significance of ebullition

  #3. Determine g.fact threshold for goflux when no bubles are detected

  #4. Determine LM quality thresholds for LM vs total.flux for non-significant ebullition. 

#Things to check------
#Find out or ask Camille

#Why the nb.obs and obs.length_diffusion are sometimes different (>1 diference) in CH4w_sep_flux?? 

#How is significance of ebullition calculated for aquaGHG warnings?:

# warnings in aquaGHG:
# if( mybest.flux$ebullition.flux < abs(mybest.flux$diffusion.flux+SD_diffusion.flux)){
#   warning(paste0("for ",id, ", ebullition term is within range of uncertainty of diffusion."))
#   mybest.flux$quality.check <- "ebullition too low to be trusted"
# }
# 
# if( mybest.flux$ebullition.flux < 0){
#   warning(paste0("for ",id, ", negative ebullition term. It was forced to 0."))
#   mybest.flux$ebullition.flux <- 0
# }
# 
# if( mybest.flux$diffusion.flux > mybest.flux$total.flux){
#   warning(paste0("for ",id, ", diffusion term is larger than total flux estimated."))
#   mybest.flux$quality.check <- "diffusion > total flux"
# }


rm(list=ls())
#Directories---------
#dropbox root path: You have to make sure this is pointing to the write folder on your local machine:
dropbox_root <- "C:/Users/Miguel/Dropbox/RESTORE4Cs - Fieldwork/Data"

#Path to Co2 and Ch4 auxfiles with corrected start.time and duration and with discard decissions:
auxfile_path<- paste0(dropbox_root,"/GHG/Working data/Auxfiles_4aquaGHG/") 

#Path to results from aquaGHG:
results_path<- "C:/Users/Miguel/Dropbox/testing_aquaGHG/"




#Packages and functions----------
library(tidyverse)
library(ggExtra)





#1. Load results------

#Load discard decissions from auxfile ch4
full_auxfilech4<- read.csv(file = paste0(auxfile_path,"ch4_auxfile.csv"))

#Load results for CH4 incubations with flux separation: 
load(paste0(results_path, "CH4all_fluxsep_aquaGHG.Rdata"))

#Load results for CH4 incubations withOUT flux separation: 
load(paste0(results_path, "CH4all_nosep_aquaGHG.Rdata"))


#2. Fluxsep inspection table----

#Create a csv table with the results from fluxsep
fluxsep_4inspection<- full_auxfilech4 %>% 
  #add duration of nosep incubation
  merge.data.frame(CH4all_nosep_flux.auto %>%  
                     select(UniqueID, nb.obs, LM.r2) %>%
                     rename(nosep_nb.obs=nb.obs), by="UniqueID") %>% 
  #add basic results from fluxsep 
  merge.data.frame(CH4all_fluxsep_flux.auto %>% 
                     select(UniqueID, nb.obs, full.total.flux, full.ebullition.flux, diffusion.flux, diffusion.flux.SD) %>% 
                     rename(dif_fluxsep_nb.obs=nb.obs), by="UniqueID") %>% 
  #0. Is incubation correct? (i.e. is not marked for discard)
  mutate(good_incubation=ch4_decission=="ok") %>% 
  #1. Is ebullition detected? (use difference in nb.obs between nosep and sep results)
  mutate(buble_detected=!(nosep_nb.obs==dif_fluxsep_nb.obs)) %>% 
  #2. Is diffusion below total flux? (diffusion.flux<full.total.flux)
  mutate(difusion_below_total=diffusion.flux<full.total.flux) %>% 
  #3. Is separation appropriate? (combination of detected and positive bubles)
  mutate(good_separation=buble_detected&difusion_below_total) %>% 
  #4. Is separation significant? (full.ebullition.flux > diffusion.flux +diffusion.flux.SD)
  mutate(sig_separation=good_separation&(full.ebullition.flux> diffusion.flux+diffusion.flux.SD)) %>% 
  #5. Is incubation very linear? (full incubation LM.r2>0.95)
  mutate(good_fulllinearfit=LM.r2>0.95) %>% 
  #Leave only logical flags and auxfile data
  select(-c(nosep_nb.obs,dif_fluxsep_nb.obs,full.total.flux,full.ebullition.flux,diffusion.flux,diffusion.flux.SD,LM.r2)) %>% 
  #Order
  arrange(subsite, start.time)


#Save csv for inspection
write.csv(fluxsep_4inspection, file = paste0(results_path, "fluxsep_inspectiontable.csv"), row.names = F)


#Load visual inspection flags for CH4 and add them to the data-driven ones: 
fluxsep_inspected<- fluxsep_4inspection %>% 
  merge.data.frame(read.csv(paste0(results_path, "fluxsep_inspectiontable_filled.csv")) %>% select(UniqueID, visual_ebullition,failure_description,common_obs,extra_obs,recheck), by="UniqueID")%>% 
  #Order
  arrange(subsite, start.time)



#2. Criteria for non-ebullition -----

#Filter results to decide goflux criteria:
#To decide the goflux hard thresholds (i.e. g.fact) we discacrd incubations that match any of the following:
  #Wrong incubations (those marked due to strong artefact, backflush after ebullition or wrong manipulation)
  #Incubations with detected ebullition (either by automatic algorithm or by visual inspection)

incub_discard<- fluxsep_inspected %>% filter(!good_incubation) %>% pull(UniqueID)
incub_ebu_auto<- fluxsep_inspected %>% filter(buble_detected) %>% pull(UniqueID)
incub_ebu_visual<- fluxsep_inspected %>% filter(visual_ebullition) %>% pull(UniqueID)

ch4_4gofluxcriteria<- CH4all_nosep_flux.auto %>% 
  filter(!UniqueID%in%c(incub_discard,incub_ebu_auto,incub_ebu_visual))

  

##2.1. Inspect mdf----
#We want to make sure MDF fluxes do not arise from artefacts, all MDF incubations will default to LM 
ch4_mdf_LM<- ch4_4gofluxcriteria %>%  
  filter(abs(LM.flux)<=MDF.lim)

#MDFs CH4 visually inspected: true below detection flat-ish white noise incubations
mdf_inspected<- c("s1-ca-p2-15-v-t-14:54","s1-ca-p2-4-v-d-11:24","s1-ca-p2-4-v-t-11:17","s1-ca-p2-6-b-d-11:53","s1-ca-r1-9-v-t-10:57","s1-ca-r1-17-b-d-12:49","s1-ca-r2-1-b-d-08:47","s1-da-a1-1-b-d-08:22","s1-da-a1-2-v-d-08:46","s1-da-a1-2-v-t-08:38","s1-da-a1-5-v-d-10:01","s1-da-a1-7-b-d-11:51","s1-da-a1-9-b-d-12:23","s1-du-a1-11-b-t-13:12","s1-du-a1-12-v-t-13:27","s1-du-a1-13-v-t-14:07","s1-du-a1-14-v-d-14:32","s1-du-a1-15-v-d-15:01","s1-du-a1-15-v-t-14:52","s1-du-a2-10-b-t-11:21","s1-du-a2-11-v-d-11:49","s1-du-a2-12-v-t-12:06","s1-du-a2-14-b-t-12:50","s1-du-a2-6-v-d-09:18","s1-du-a2-6-v-t-09:09","s1-du-p1-14-v-d-13:13","s1-du-p1-9-v-d-11:37","s1-du-p2-2-v-t-07:39","s1-du-p2-3-o-d-07:56","s1-du-p2-5-b-t-08:39","s1-du-p2-6-v-d-09:00","s1-du-p2-7-o-d-09:16","s1-du-p2-14-v-d-11:51","s1-du-p2-11-v-d-10:43","s1-du-p2-11-v-t-10:36","s1-du-p2-13-v-d-11:31","s1-du-p2-13-v-t-11:24","s1-du-p2-14-v-t-11:43","s1-du-p2-15-v-d-12:17","s1-du-p2-15-v-t-12:09","s1-du-p2-2-v-d-07:48","s1-du-p2-4-v-d-08:23","s1-du-p2-4-v-t-08:12","s1-du-p2-6-v-t-08:52","s1-du-p2-8-v-d-09:36","s1-du-p2-9-b-t-09:44","s1-du-r1-3-v-t-08:24","s1-du-r2-1-v-d-07:29","s1-du-r2-10-v-d-11:03","s1-du-r2-12-v-t-11:34","s1-du-r2-13-v-d-12:21","s1-du-r2-14-v-t-12:33","s1-du-r2-15-v-d-12:56","s1-du-r2-4-b-d-08:20","s1-va-p1-11-v-d-11:52","s1-va-p1-11-v-t-11:45","s1-va-p1-13-v-d-12:48","s1-va-p1-13-v-t-12:41","s1-va-p1-3-b-d-09:20","s1-va-p1-4-b-d-09:30","s1-va-p1-9-v-d-11:22","s1-va-r2-2-b-d-09:33","s1-va-r2-3-v-t-09:41","s1-va-r2-5-v-d-10:22",
                  "s2-ca-a1-4-v-t-09:38","s2-ca-a1-5-v-d-10:15","s2-ca-a1-5-v-t-10:02","s2-ca-a1-6-v-d-10:37","s2-ca-a1-6-v-t-10:30","s2-ca-a1-7-v-d-10:56","s2-ca-a1-7-v-t-10:47","s2-ca-r1-14-v-d-12:21","s2-cu-r1-1-b-d-08:05","s2-cu-r1-2-b-d-08:15","s2-cu-r1-3-b-d-08:28","s2-da-a1-7-b-d-10:36","s2-da-a1-7-b-t-10:41","s2-da-a1-7-v-d-10:28","s2-da-a1-7-v-t-10:20","s2-da-a1-9-b-d-11:27","s2-da-a1-9-b-t-11:21","s2-du-a1-10-v-d-12:11","s2-da-r1-2-v-d-09:22","s2-du-a1-10-v-t-12:04","s2-du-a1-11-v-d-12:28","s2-du-a1-11-v-t-12:22","s2-du-a1-13-v-d-12:58","s2-du-a1-13-v-t-12:53","s2-du-a1-14-v-d-13:16","s2-du-a1-14-v-t-13:11","s2-du-a1-15-v-d-13:35","s2-du-a1-15-v-t-13:29","s2-du-a1-16-v-d-14:05","s2-du-a1-16-v-t-13:57","s2-du-a1-4-v-d-10:13","s2-du-a1-5-b-d-10:30","s2-du-a2-3-v-t-09:34","s2-du-a2-10-b-d-12:03","s2-du-a2-10-b-t-11:50","s2-du-a2-4-v-d-09:57","s2-du-a2-8-v-d-11:19","s2-du-a2-8-v-t-11:13","s2-du-p1-3-o-d-09:31","s2-du-p1-10-b-d-11:29","s2-du-p1-11-v-t-11:44","s2-du-p1-14-v-d-12:42","s2-du-p1-2-v-t-09:15","s2-du-p1-4-v-t-09:49","s2-du-p1-8-b-d-10:56","s2-du-p2-10-v-d-12:14","s2-du-p2-10-v-t-12:09","s2-du-p2-11-v-d-12:34","s2-du-p2-11-v-t-12:27","s2-du-p2-12-v-d-12:46","s2-du-p2-12-v-t-12:41","s2-du-p2-13-b-d-13:01","s2-du-p2-13-b-t-12:54","s2-du-p2-14-v-d-13:16","s2-du-p2-2-v-d-09:41","s2-du-p2-2-v-t-09:33","s2-du-p2-3-o-d-10:11","s2-du-p2-4-b-d-10:32","s2-du-p2-4-b-t-10:25","s2-du-p2-5-o-d-10:43","s2-du-p2-7-v-d-11:31","s2-du-p2-7-v-t-11:21","s2-du-p2-8-v-d-11:47","s2-du-p2-8-v-t-11:41","s2-du-p2-9-b-d-11:59","s2-du-p2-9-b-t-11:55","s2-du-r1-15-v-t-12:55","s2-du-r1-2-v-t-09:35","s2-du-r1-7-v-d-10:54","s2-du-r1-7-v-t-10:41","s2-du-r2-8-o-d-11:11","s2-du-r2-1-v-d-09:29","s2-du-r2-10-v-d-11:48","s2-du-r2-13-v-d-12:59","s2-du-r2-14-v-d-13:16","s2-du-r2-14-v-t-13:09","s2-du-r2-15-v-t-13:25","s2-du-r2-6-v-d-10:40","s2-ri-a2-3-b-d-08:59","s2-ri-a2-9-v-t-10:19","s2-ri-r1-9-v-t-12:33",
                  "s3-ca-a1-9-v-t-10:01","s3-ca-r1-10-v-d-09:53","s3-ca-r1-10-v-t-09:44","s3-ca-r1-15-b-d-11:00","s3-ca-r1-5-v-d-08:42","s3-ca-r1-5-v-t-08:35","s3-cu-r2-2-b-d-07:07","s3-cu-r2-3-b-d-07:18","s3-da-a1-1-b-d-07:07","s3-da-a1-3-b-d-07:51","s3-da-a1-4-v-t-08:05","s3-da-a1-6-v-t-09:10","s3-da-a1-6-b-d-09:23","s3-da-a1-11-b-d-11:29","s3-da-a1-11-v-d-11:21","s3-da-a1-11-v-t-11:12","s3-da-a1-12-b-d-12:09","s3-da-a1-12-v-d-12:00","s3-da-a1-12-v-t-11:51","s3-da-a1-2-v-t-07:17","s3-da-a1-4-b-d-08:22","s3-da-a1-5-b-d-08:55","s3-da-a1-6-v-d-09:17","s3-da-a1-7-b-d-09:40","s3-da-a1-8-b-d-09:59","s3-da-a1-9-b-d-10:19","s3-da-p2-2-v-t-09:16","s3-du-a1-10-v-d-10:42","s3-du-a1-10-v-t-10:33","s3-du-a1-12-v-d-11:30","s3-du-a1-13-v-d-11:47","s3-du-a1-13-v-t-11:41","s3-du-a1-14-v-d-12:05","s3-du-a1-14-v-t-11:58","s3-du-a1-7-b-t-09:28","s3-du-a1-9-v-d-10:19","s3-du-a1-9-v-t-10:09","s3-du-a1-12-v-t-11:23","s3-du-a2-11-b-d-10:55","s3-du-a2-11-b-t-10:49","s3-du-a2-12-v-t-11:09","s3-du-a2-7-v-d-09:23","s3-du-a2-7-v-t-09:14","s3-du-a2-8-v-t-09:37","s3-du-a2-6-v-t-08:53","s3-du-a2-6-v-d-09:01","s3-du-a2-9-b-d-10:09","s3-du-a2-9-b-t-10:02","s3-du-p1-2-o-d-07:47","s3-du-p1-15-v-d-11:31",  "s3-du-p2-11-v-d-11:59","s3-du-p2-11-v-t-11:51","s3-du-p2-12-v-d-12:24","s3-du-p2-12-v-t-12:16","s3-du-p2-15-v-d-13:16","s3-du-p2-15-v-t-13:09","s3-du-p2-2-v-d-08:33","s3-du-p2-2-v-t-08:22","s3-du-p2-3-v-d-08:57","s3-du-p2-4-v-d-09:18","s3-du-p2-4-v-t-09:10","s3-du-p2-8-v-d-11:01","s3-du-p2-8-v-t-10:52","s3-du-r1-2-v-d-07:36","s3-du-r1-2-v-t-07:30","s3-va-p1-13-v-d-11:58","s3-va-p1-5-v-d-09:24","s3-va-p1-5-v-t-09:19",
                  "s4-ca-r1-2-v-t-07:55","s4-ca-r1-3-b-d-08:12","s4-ca-r1-4-v-t-08:19","s4-cu-r1-4-b-d-07:34","s4-cu-r2-2-b-d-06:55","s4-da-a1-3-v-t-08:02","s4-da-a1-5-b-d-08:32","s4-da-a1-7-v-t-09:22","s4-da-a1-10-b-d-10:04","s4-da-a1-2-b-d-07:48","s4-da-a1-8-v-t-09:40","s4-du-p1-12-v-t-10:40","s4-du-p1-14-v-d-11:26","s4-du-p1-4-v-t-08:37","s4-du-p1-6-v-d-09:17","s4-du-p1-2-v-t-08:05","s4-du-p1-13-v-t-10:50","s4-du-r2-1-v-d-08:20","s4-du-r2-12-v-t-11:06","s4-du-r2-2-v-t-08:29","s4-du-r2-3-b-d-08:47","s4-du-r2-9-v-d-10:17","s4-ri-a2-1-b-d-06:52","s4-ri-a2-6-b-d-07:41","s4-ri-a2-9-b-d-08:27","s4-va-p1-12-v-t-10:11","s4-va-p1-7-v-t-08:56","s4-va-p1-9-v-t-09:20","s4-va-r2-1-b-d-08:01","s4-va-r2-11-v-t-10:08","s4-va-r2-14-v-d-10:46","s4-va-r2-15-v-d-11:01","s4-va-r2-15-v-t-10:54","s4-va-r2-4-v-t-08:35"
)

#any MDF flux not inspected?
ch4_mdf_LM %>% filter(!UniqueID%in%mdf_inspected) %>% pull(UniqueID)

#any mdf_inspected that shouldnt?
mdf_inspected[!mdf_inspected%in%ch4_mdf_LM$UniqueID]


##2.2. Inspect Kmax-----
#Inspect HM.k>=k.max () to ensure that the LM flux is not influenced by artefacts
ch4_kmax<- ch4_4gofluxcriteria %>%
  filter(!UniqueID%in%mdf_inspected) %>% #Already inspected (due to MDF)
  filter(HM.k>=k.max)

#In kmax_inspected, incubations with k=kmax with good LM estimates:
kmax_inspected<-c("s1-ca-p1-3-v-d-10:31","s1-ca-p2-3-v-d-10:59","s1-da-a1-4-v-t-09:15","s1-da-a1-4-v-d-09:22","s1-da-a1-7-v-d-11:35","s1-du-a1-10-v-t-12:29","s1-du-a2-4-v-t-07:53","s1-du-a2-4-v-d-08:01","s1-du-p1-1-o-d-08:11","s1-du-p1-6-v-d-10:30","s1-du-p1-9-v-t-11:30","s1-du-p1-13-v-t-12:47","s1-du-p1-10-v-d-12:04","s1-du-p1-15-v-t-13:27","s1-du-p1-15-v-d-13:33","s1-du-r1-9-v-d-10:00","s1-du-r2-1-v-t-07:23","s1-du-r2-6-v-t-08:51","s1-du-r2-10-v-t-10:57","s1-du-r2-11-v-d-11:22","s1-du-r2-12-v-d-11:40","s1-du-r2-15-v-t-12:50","s1-va-a1-8-b-d-12:46",
                  "s2-ca-a1-4-v-d-09:50","s2-ca-a2-7-v-t-09:40","s2-ca-p1-2-v-d-08:54","s2-ca-r1-8-v-d-10:45","s2-da-r1-2-v-t-09:15","s2-du-p1-2-v-d-09:22","s2-du-p1-9-v-t-11:10","s2-du-p1-6-o-d-10:23","s2-du-p1-7-o-d-10:40","s2-du-p1-12-v-t-12:00","s2-du-p1-14-v-t-12:34","s2-du-r2-2-v-t-09:37","s2-du-r2-2-v-d-09:45","s2-du-r2-7-v-d-10:55","s2-du-r2-11-v-t-11:58","s2-du-r2-13-v-t-12:52","s2-ri-a2-7-b-d-10:00","s2-ri-p1-6-v-t-10:17","s2-ri-r2-12-v-t-12:35","s2-va-p1-2-v-d-09:21","s2-va-r2-1-b-t-09:44",
                  "s3-ca-a1-7-v-d-09:35","s3-ca-p2-1-v-d-07:49","s3-da-a1-2-v-d-07:25","s3-da-a1-4-v-d-08:12","s3-du-a2-3-v-t-07:43","s3-du-p1-6-v-d-08:51","s3-du-p1-15-v-t-11:24","s3-du-r1-1-v-t-07:15","s3-du-r1-9-v-t-09:18","s3-va-p1-12-b-d-11:38",
                  "s4-ca-p1-14-v-d-11:52","s4-ca-r2-2-v-d-07:27","s4-da-a1-4-v-t-08:16","s4-du-p1-4-v-d-08:45","s4-du-p1-11-b-d-10:30","s4-du-p2-8-v-d-10:04","s4-du-p2-14-v-d-12:08","s4-du-p2-15-v-d-12:25","s4-du-r2-2-v-d-08:36","s4-du-r2-10-v-t-10:26","s4-du-r2-11-v-d-10:56","s4-ri-a1-5-b-d-07:53","s4-ri-a2-7-v-d-07:58","s4-ri-a2-11-b-d-08:54","s4-ri-p1-1-v-d-06:59","s4-ri-p1-4-v-t-07:38","s4-ri-p1-5-o-d-07:54","s4-ri-p1-11-v-t-09:30","s4-va-p1-14-v-d-10:47","s4-va-p2-10-v-t-08:42","s4-va-p2-14-v-t-09:27","s4-va-r2-3-v-d-08:27","s4-va-r2-11-v-d-10:15")


tocropand_reinspect<- c(
                        )

#Any kmax not inspected?
ch4_kmax %>% filter(!UniqueID%in%c(kmax_inspected,tocropand_reinspect)) %>% pull(UniqueID)

#Any kmax_inspected that shouldnt?
kmax_inspected[!kmax_inspected%in%ch4_kmax$UniqueID]
 


##2.3. Inspect g.fact ------

#Filter the dataset to remove already inspected incubations (mdf, kappamax, tocrop)
ch4_gfact<-ch4_4gofluxcriteria %>%
  filter(!UniqueID%in%c(mdf_inspected,kmax_inspected,tocropand_reinspect)) 

#lets see what is the highest g.fact of the best-performing HM models (best 50%)
ch4_gfact %>% 
  filter(HM.r2>quantile(HM.r2, 0.5, na.rm=T)) %>%
  ggplot( aes(x=HM.r2, y=g.fact))+
  geom_point()+
  geom_hline(yintercept = 3)


#Here we check that the hard-threshold chosen for g.fact (i.e. if g.fact > threshold --> LM as best flux) does not cause any truly non-linear pattern to be lost.


#Based on the g.fact distribution above, we will start the inspection with g.fact=3 as threshold and progressively lower it until clear and clean non-linear patterns appear

#These are All the incubations with unreasonably large g-fact. They have all been visually inspected and LM is more representative for them. 
gfact_inspected<- c("s1-ca-p1-8-v-t-12:12","s1-ca-p1-11-v-d-13:35","s1-ca-p1-14-v-d-14:49","s1-ca-p1-15-v-t-15:02","s1-ca-r1-3-v-d-09:33","s1-ca-r2-6-v-d-09:56","s1-da-r1-5-o-d-09:05","s1-du-p1-8-v-d-11:13","s1-du-r1-9-v-t-09:54","s1-ri-p2-1-v-d-09:37","s1-ri-r1-9-b-t-11:39","s1-ri-r1-10-v-t-11:56","s1-ri-r1-13-v-t-12:48","s1-va-a1-7-b-d-12:37","s1-va-p2-8-b-t-11:12","s2-ca-p1-6-v-t-10:06","s2-ca-p2-8-b-t-10:50","s2-ri-a1-1-b-t-08:50","s2-ri-a1-13-b-t-10:57","s2-ri-p2-3-v-d-09:15","s2-ri-p2-8-v-d-10:58","s3-ca-a1-10-v-d-10:23","s3-cu-p1-2-v-d-06:50","s3-du-p1-9-v-d-09:38","s3-ri-r1-8-v-d-11:36","s3-ri-r1-14-v-t-14:12","s4-ca-a1-12-v-d-12:55","s4-ca-p1-8-v-d-09:50","s4-ca-p2-10-v-t-12:50","s4-ca-r2-7-v-d-08:46","s4-du-a1-7-v-d-10:09","s4-du-a1-11-v-t-11:40","s4-du-a1-12-v-d-12:04","s4-du-p2-2-v-d-08:16","s4-ri-a1-10-b-d-08:40","s4-ri-r2-4-v-t-07:32","s4-va-p1-6-b-d-08:47","s4-va-p2-3-v-t-07:19","s4-va-p2-13-v-t-09:11","s4-du-a1-15-v-t-12:53","s4-va-a1-3-b-d-09:10","s1-ca-r2-15-v-d-13:41","s1-du-p1-5-v-t-09:56","s1-du-p1-14-v-t-13:07","s1-ri-p2-11-v-d-12:44","s1-ri-r1-2-v-d-09:58","s1-ri-r1-3-v-t-10:06","s2-ca-p1-1-v-t-08:27","s2-ri-r1-7-v-t-11:59","s3-du-r2-11-v-d-11:30","s4-ca-p1-15-v-t-12:13","s4-ca-p2-11-v-t-13:08","s4-du-p1-6-v-t-09:11","s4-ri-a1-2-b-d-07:22","s4-ri-r1-5-v-t-08:11","s4-ri-r1-6-v-d-08:35","s4-va-p1-10-v-t-09:34","s1-ca-r2-15-v-t-13:34","s1-ca-p1-6-v-d-11:24","s1-cu-a1-13-v-d-10:47","s1-du-p1-7-b-d-10:45","s1-du-p1-10-v-t-11:57","s1-ri-a1-8-b-t-12:05","s1-ri-r1-15-v-d-13:35","s1-va-a1-1-o-d-09:49","s2-cu-p1-11-v-d-09:57","s2-du-a1-1-o-d-09:21","s2-ri-p2-1-o-d-08:32","s2-va-p1-6-o-d-10:34","s2-va-r2-1-b-d-09:52","s3-da-p2-9-v-t-10:51","s3-du-r1-11-v-t-09:53","s3-du-r2-7-o-d-10:12","s3-va-a1-2-v-d-09:13","s4-du-a2-7-v-d-10:19","s4-du-a2-11-v-t-11:34","s4-du-p1-13-v-d-11:03","s4-du-r2-6-o-d-09:19","s1-ca-p1-1-o-d-09:25","s3-da-p1-14-v-d-11:58")#gfact>3


#To fix and re-inspect: too long incubations, aproaching asyntotic concentration
toreinspect_gfact<- c(
                      )


#Any g.fact to inspect?
ch4_gfact %>% 
  filter(!UniqueID%in%c(gfact_inspected,toreinspect_gfact)) %>% #already inspected
  filter(g.fact>=3) %>% 
  select(UniqueID, g.fact)

#any gfact_inspected that shouldnt?
gfact_inspected[!gfact_inspected%in%ch4_gfact$UniqueID]


#g.fact threshold of 3 is the optimal solution, many true non-linear with g.fact >2.5


##2.4. AICc & MAE improvement-----

#To decide for the remaining incubations (those non-ebullitive, non-mdf, non-kmax, non-gfact>3), we will use a combination of AICc weight (probability 0-1 of HM model being best) and relative improvement on MAE with HM.
#Compare distributions for each criteria, inspect critical decissions. 


#After checking that g.fact and mdf criteria are appropriate, subset dataset that needs a decision between LM and HM:
ch4_todecide<- ch4_4gofluxcriteria %>% #Non-wrong, non-ebullitive (auto or visual)
  filter(!is.na(HM.flux)) %>% #Non-NA HM (nothing to decide in those cases, default LM)
  filter(!UniqueID%in%c(mdf_inspected,kmax_inspected,tocropand_reinspect, gfact_inspected, toreinspect_gfact))#Non-mdf, non-kmax, non-gfact>3


#Inspect AICc for selection criteria (using AICc weight)
#Test criteria: using AICc weight (based on the relative difference in AICc, calculate the AICc weight metric, a probabilistic estimate of how likely it is that one model is better than the other). 
#ADD relative improvement in MAE as additional threshold criteria (for HM to be selected, it has to reduce model error by more than 5%).

ch4_todecide <- ch4_todecide %>%
  mutate(
    delta_AICc_LM = LM.AICc - pmin(LM.AICc, HM.AICc),  # Compare to the best AICc (min AICc)
    delta_AICc_HM = HM.AICc - pmin(LM.AICc, HM.AICc),
    # Calculate the Akaike weight for HM: AICc_weight_HM, how likely is that HM is better (0-1)?
    AICc_weight_HM = exp(-0.5 * delta_AICc_HM) / (exp(-0.5 * delta_AICc_LM) + exp(-0.5 * delta_AICc_HM)),
    # Calculate relative improvement in MAE using HM
    rel_reduct_MAE_HM  = (LM.MAE  - HM.MAE)  / LM.MAE,
    k.ratio=HM.k/k.max
  )

#Inspect AICc_weight_HM
ch4_todecide %>%
  ggplot(aes(x=AICc_weight_HM, fill=AICc_weight_HM>0.5))+
  geom_histogram(bins=100)+
  scale_x_continuous(name="AICc weight for HM, probability of HM being best")
#Most cases have very clear separation. 

#Inspect against LM r2 (rel_reduction in MAE as a guide)
ch4_todecide %>%
  ggplot(aes(x=AICc_weight_HM, y=LM.r2,col=rel_reduct_MAE_HM>0.05))+
  geom_point()+
  scale_x_continuous(name="AICc weight for HM, probability of HM being best")
#Many incubations that win with AICc weight, have poor LM r2 (indicating noise) and do not improve MAE by more than 5%


#Inspect rel_reduct_MAE_HM
ch4_todecide %>%
  ggplot(aes(x=rel_reduct_MAE_HM, fill=AICc_weight_HM>0.5))+
  geom_histogram(bins=100)+
  scale_x_continuous(name="Relative reduction in MAE with HM")

ch4_todecide %>%
  ggplot(aes(x=rel_reduct_MAE_HM, y=LM.r2,col=AICc_weight_HM>0.5))+
  geom_point()+
  geom_vline(xintercept=0.025)+
  scale_x_continuous(name="Relative reduction in MAE with HM")


#Inspect against k.ratio 
ggMarginal(
  ch4_todecide %>%
    ggplot(aes(x=AICc_weight_HM,y=k.ratio, col=rel_reduct_MAE_HM>0.05))+
    geom_point()+
    # geom_hline(yintercept=1.25)+
    labs(col=NULL)+
    scale_x_continuous(name="AICc weight for HM, probability of HM being best")
)



#Inpsect cases with better AICc for HM, but low reduction in MAE (<5%) and with relatively big differences in flux (gfact>1.5)

#Inspect to decide 
ch4_lowMAEimprovement<-ch4_todecide %>% 
  filter(AICc_weight_HM>0.5) %>% 
  filter(rel_reduct_MAE_HM<0.05) %>% 
  filter(g.fact>1.5)

#This cases have been inspected: tipical behaviour is linear with a few artefacts forcing a curve in the model. 
bestLM<- c("s1-ca-p1-4-v-t-10:43","s1-ca-p1-9-v-t-12:32","s1-ca-p1-15-v-d-15:09","s1-ca-r2-8-v-t-10:25","s1-da-p2-3-v-d-11:07","s1-da-r2-14-b-t-12:26","s1-du-p1-11-v-t-12:17","s2-ca-p1-4-v-t-09:26","s2-ca-p1-10-o-d-12:22","s2-ca-p1-11-v-d-12:44","s2-ca-p2-12-v-t-12:36","s2-ca-r1-8-v-t-10:39","s2-ca-r2-10-v-t-11:07","s2-cu-r2-4-v-d-09:02","s2-du-r1-3-v-t-09:52","s2-du-r1-15-v-d-13:01","s2-du-r2-7-v-t-10:49","s2-du-r2-9-o-d-11:25","s2-ri-a1-2-b-d-09:15","s2-ri-a1-10-b-t-10:30","s2-ri-p2-4-v-t-09:27","s2-ri-p2-10-v-d-12:08","s2-ri-p2-11-v-t-12:17","s2-ri-r1-2-v-t-10:21","s2-ri-r1-3-v-d-10:51","s2-ri-r1-4-v-d-11:12","s2-ri-r1-5-v-t-11:20","s2-ri-r1-5-v-d-11:30","s2-ri-r1-6-v-d-11:47","s2-va-r2-4-o-d-10:34","s3-ca-p1-3-v-t-08:00","s3-da-a1-5-v-t-08:40","s3-du-r1-12-v-d-10:16","s3-du-r2-2-v-t-09:04","s3-du-r2-9-v-t-10:46","s3-ri-p2-2-v-t-09:18","s3-ri-r1-6-v-t-10:41","s3-ri-r1-6-b-d-11:04","s3-va-a1-2-v-t-09:01","s4-ca-a1-13-v-d-13:31","s4-ca-p1-15-v-d-12:21","s4-ca-p2-8-b-d-12:16","s4-ca-p2-9-v-t-12:29","s4-du-a1-6-b-d-09:38","s4-du-p1-1-b-d-07:55","s4-du-p1-8-b-d-09:43","s4-du-p2-12-b-d-11:28","s4-du-p2-13-v-t-11:40","s4-du-r1-3-v-t-08:52","s4-ri-a2-4-v-t-07:19","s4-ri-a2-4-v-d-07:25","s4-ri-a2-10-v-t-08:40","s4-ri-p1-7-o-d-08:26","s4-ri-p1-12-v-t-09:44","s4-ri-p1-12-v-d-09:49","s4-ri-p2-1-v-d-07:09","s4-ri-p2-8-v-t-09:04","s4-ri-p2-10-v-t-09:38","s4-ri-r1-1-v-d-07:11","s4-ri-r1-9-v-d-10:13","s4-ri-r2-12-v-t-10:12","s4-va-a2-11-b-d-10:18","s4-va-p2-5-v-t-07:42","s4-va-p2-14-v-d-09:34","s4-va-r2-3-v-t-08:19","s4-va-r2-4-v-d-08:42","s1-ri-p2-4-v-t-10:24","s1-ri-p2-12-o-d-12:51","s1-ri-r1-8-v-t-11:24","s2-du-p1-4-v-d-09:57","s2-ri-r1-8-v-t-12:16","s2-ri-r1-10-v-t-12:51","s3-da-a1-5-v-d-08:47","s3-du-r2-9-v-d-10:53","s3-va-r2-3-v-d-09:10","s4-ca-p1-10-v-t-10:44","s4-cu-r1-1-b-d-07:06","s4-du-a2-10-v-t-11:12","s4-du-p2-4-o-d-08:47","s4-ri-a1-3-b-d-07:32","s4-ri-a1-4-b-d-07:42","s4-ri-p2-4-o-d-07:55","s4-ri-p2-11-v-t-10:06","s4-ri-r1-5-v-d-08:19","s4-ri-r1-8-v-t-09:22","s4-ri-r2-10-v-t-09:40")

noclear<- c("s1-cu-r1-15-v-d-09:44","s1-du-a1-6-b-t-10:30","s1-ri-r1-14-v-d-13:20","s2-ca-p1-9-o-d-12:12","s2-cu-p1-11-v-t-09:50","s2-va-p1-8-o-d-11:07","s3-ri-a1-1-b-t-09:30","s3-ri-p2-7-v-t-10:58","s3-ri-r1-13-b-d-14:04","s3-ri-r1-14-v-d-14:19","s3-va-r2-4-b-d-09:18","s3-va-r2-5-v-t-09:27","s4-ri-p1-10-v-t-09:10","s4-ri-p2-2-o-d-07:17","s4-ri-r2-11-v-t-09:53","s4-ri-r2-11-v-d-10:01","s1-ca-p1-13-o-d-14:13","s1-cu-r1-13-v-d-09:02","s1-ri-r1-14-v-t-13:12","s2-ca-r2-7-v-d-10:19","s2-va-r2-9-v-d-11:47")

bestHM<-c()

discard_artefacts<-c()

crop_reinspect<- c()


#Any to inspect?
ch4_lowMAEimprovement %>% 
  filter(!UniqueID%in%c(bestLM,bestHM,noclear,discard_artefacts,crop_reinspect)) %>%
  select(UniqueID, g.fact,k.ratio, rel_reduct_MAE_HM) %>% head()

#Any inspected that shouldnt?
c(bestLM,noclear)[!c(bestLM,noclear)%in%ch4_lowMAEimprovement$UniqueID]


#These criteria work well for all, a few-low g.fact HM unclear, but overall very good. 

#quality check: inspect k.ratio of selected HM and inspect visually.
k_final_good<-c("s1-da-a1-3-b-d-08:56","s3-va-p1-11-b-d-11:28","s2-ca-a1-8-b-d-11:08")

ch4_todecide %>% 
  filter(AICc_weight_HM>0.5) %>% 
  filter(rel_reduct_MAE_HM>=0.05) %>% 
  filter(k.ratio>0.5) %>% 
  filter(!UniqueID%in%c(crop_reinspect, k_final_good)) %>% 
  select(UniqueID,g.fact,k.ratio)


#quality check: inspect MAE distribution of selected HM and LM. Visually review worse incubations. 

#HM chosen: check high MAE, worse-scenario is ebullitive dynamics and HM chosen.
ch4_todecide %>% 
  filter(AICc_weight_HM>0.5) %>% 
  filter(rel_reduct_MAE_HM>=0.05) %>% 
  filter(HM.MAE>15) %>% 
  arrange(desc(HM.MAE)) %>% 
  filter(g.fact>1.1) %>% 
  filter(LM.flux>5) %>%
  select(UniqueID, HM.MAE, g.fact, k.ratio)%>% 
  filter(!UniqueID%in%c(hm_apropriate,lm_apropriate,non_detected_ebullition)) %>% head()

hm_apropriate<- c("s2-ca-p2-15-b-t-13:23","s4-cu-a2-3-o-d-08:08","s2-ri-r2-2-o-d-12:47","s4-cu-a2-2-v-d-07:54","s3-va-r2-7-b-d-10:02","s4-va-a2-1-b-d-07:13","s4-va-a2-2-o-d-07:21","s2-cu-r1-9-o-d-10:03","s4-cu-a2-2-v-t-07:45","s2-da-r1-1-o-d-09:00")
non_detected_ebullition<-c()


lm_apropriate<-c("s1-ri-r2-3-o-d-14:05","s3-ca-a1-2-v-d-07:56","s1-ca-a2-12-v-d-13:07")


ch4_todecide %>% 
  filter(AICc_weight_HM>0.5) %>% 
  filter(rel_reduct_MAE_HM>=0.05) %>% 
  filter(!UniqueID%in%c(hm_apropriate,lm_apropriate,non_detected_ebullition)) %>% 
  ggplot(aes(x=HM.MAE))+
  geom_histogram()


#quality check: inspect MAE of selected LM (not critical, it is what it is)
ch4_todecide %>% 
  filter(rel_reduct_MAE_HM<0.05) %>% 
  filter(!UniqueID%in%c(non_detected_ebullition,lm_apropriate)) %>% 
  ggplot(aes(x=LM.MAE))+
  geom_histogram()

ch4_todecide %>% 
  filter(rel_reduct_MAE_HM<0.05) %>% 
  filter(LM.flux>5) %>% 
  filter(!UniqueID%in%c(non_detected_ebullition,lm_apropriate)) %>%
  ggplot(aes(x=LM.MAE,y=LM.flux/full.total.flux))+
  geom_point()+
  geom_hline(yintercept = 1)

non_detected_ebullition<-c()

lm_apropriate<- c("s1-ri-a1-16-o-d-14:00","s1-va-r1-7-o-d-10:38")


ch4_todecide %>% 
  filter(rel_reduct_MAE_HM<0.05) %>% 
  filter(LM.MAE>10) %>% 
  filter(LM.flux>5) %>% 
  arrange(desc(LM.MAE)) %>% 
  filter(!UniqueID%in%c(non_detected_ebullition,lm_apropriate)) %>% 
  mutate(lmtotalratio=LM.flux/full.total.flux) %>% 
  filter(!between(lmtotalratio,0.9,1.1)) %>% 
  select(UniqueID,LM.MAE, g.fact,LM.r2,lmtotalratio) %>% head()



#3. Criteria for ebullition-----


##3.1. LM.r2 threshold-----
#If ebullition is detected, either by automatic algorithm or by visual inspection: 

  #IF LM.r2 > threshold ---> select LM
  #IF LM.r2 < threshold ---> select totalflux 

#Check threshold to chose LM over totalflux, when ebullition is present, total.flux should be more reliable than linear depending on the shape. check difference Linear/total to decide LM.r2 threshold 

ch4_ebullitive<- CH4all_nosep_flux.auto %>% 
  filter(!UniqueID%in%incub_discard) %>% #Not wrong incubations
  filter(UniqueID%in%c(incub_ebu_visual,incub_ebu_auto)) %>% #Contains ebullition (auto or visual)
  mutate(linear_total_ratio=LM.flux/full.total.flux)

ggMarginal(
  ch4_ebullitive %>% 
    filter(LM.r2>0.95) %>% 
    ggplot(aes(x=LM.r2, y=linear_total_ratio))+
    geom_point()
)

#Threshold to chose LM over total.flux will be set at LM.r2 0.99, as below this threshold there are a lot of incubations whose linear flux deviates more than 10% from the total.flux


#POR AQUI bien_____________-----

#APPLY combined criteria (goflux+ebullition) for best.flux

#IMPLEMENT flags and final format for CH4 bestflux

#IMPLEMENT same goflux criteria (including MAE improvement for CO2), copy-paste and re-inspect thresholds (for cases not already done in testingbesflux criteria)

#IMPLEMENT flags for ch4 flux separation.

#QUality check: compare look which method is assigned to biggest fluxes (HM should not be present among the highest fluxes, as these most likely arise from ebullitive patterns)


#4. CH4 flux selection-------

#We apply the criteria above to select the best.flux in every case (ebullition or not).






#5. Flux separation -----

#Join and flag results from Flux separation process







##3.1. Inspect LM vs total.flux-----

linear_better<-c()
totalflux_better<-c("s1-cu-a1-3-o-d-07:40")

ch4_ebullitive %>% 
  filter(prefered_model=="LM") %>% 
  filter(!UniqueID%in%c(linear_better,totalflux_better)) %>% 
  filter(!between(linear_total_ratio, 0.9,1.1)) %>% 
  select(UniqueID,LM.r2, g.fact, linear_total_ratio) %>% 
  arrange(desc(LM.r2))




ch4_ebullitive %>% 
  filter(prefered_model=="LM") %>% 
  filter(!UniqueID%in%c(linear_better,totalflux_better)) %>% 
  ggplot(aes(x=prefered_model, y=linear_total_ratio))+
  geom_violin()+
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               dotsize = 0.5)


##3.1. Inpsect Non-detected----
# Add a flag ("Visual ebullition, HM not appropriate") and default to ebullition criteria (LM.r2>0.95-->LM, otherwise totalflux)

ch4_non_detected<-CH4all_nosep_flux.auto %>% 
  filter(!UniqueID%in%incub_discard) %>% 
  filter(UniqueID%in%incub_ebu_visual) %>% #Contains visual ebullition
  filter(!UniqueID%in%incub_ebu_auto) %>% #Was not detected as ebullition
  mutate(prefered_model=if_else(LM.r2>0.95,"LM","total.flux"),
         linear_total_ratio=LM.flux/full.total.flux)

#check linear to total.flux of non-detected incubations with
ch4_non_detected %>%
  filter(prefered_model=="LM") %>% 
  ggplot(aes(x=prefered_model, y=linear_total_ratio, col=prefered_model))+
  geom_violin()+
  geom_point()

lm_correct<- c("s1-ca-a2-4-v-d-11:25","s1-ca-a2-9-o-d-12:27","s1-cu-r1-4-o-d-06:55","s1-da-p1-14-v-t-11:32","s2-cu-r2-14-o-d-11:55","s2-cu-r2-13-o-d-11:45","s2-da-a2-1-o-d-07:32","s2-da-a2-2-o-d-07:52","s2-da-a2-3-o-d-08:04","s2-da-a2-4-o-d-08:16","s2-da-a2-5-o-d-08:26","s2-da-a2-6-o-d-08:36","s2-da-a2-8-o-d-09:02","s2-da-a2-9-o-d-09:17","s2-da-a2-12-o-d-10:06","s2-da-a2-13-o-d-10:17","s2-da-a2-14-o-d-10:29","s2-da-a2-15-o-d-10:46","s2-da-p1-10-o-d-09:56","s2-da-p1-11-o-d-10:08","s3-ca-a2-17-o-d-11:08","s4-va-a2-9-v-t-09:24","s4-va-a2-9-v-d-09:32","s4-va-p1-5-b-d-08:19")

totalflux_correct<- c("s1-ca-a2-4-v-t-11:19","s1-ca-a2-10-v-d-12:40","s1-ca-a2-13-o-d-13:17", "s1-ca-a2-15-o-d-13:43","s1-ca-r1-12-o-d-11:45", "s1-cu-p1-19-v-d-14:57", "s1-cu-r2-13-o-d-08:55", "s1-ri-r2-4-b-d-10:01", "s1-va-a2-4-v-t-09:54", "s1-va-a2-4-v-d-10:01", "s1-va-a2-8-v-t-11:07","s1-va-a2-12-v-d-12:40","s1-va-r2-15-o-d-13:13","s2-ca-a2-6-b-d-09:31","s2-cu-a2-8-o-d-10:04","s2-cu-p2-5-o-d-09:36","s2-da-a2-10-o-d-09:31","s2-da-a2-11-o-d-09:57","s2-va-r1-1-o-d-09:18","s3-da-a2-11-b-d-09:20","s3-da-r1-8-v-d-10:06","s3-da-r1-8-v-d-10:06","s4-ca-a2-15-v-d-10:54","s4-va-a1-8-v-t-10:58")

reclasify<- c()

#check
ch4_non_detected %>% 
  filter(prefered_model=="LM") %>%
  filter(between(linear_total_ratio,0.95,1.05)) %>% 
  select(UniqueID, g.fact, prefered_model) %>% 
  filter(!UniqueID%in%c(lm_correct, totalflux_correct,reclasify)) %>% head()




##3.2. Inspect False-positives------
#Some incubations were marked as ebullitive by the algorithm due to artefacts and noise, but they do not show bubbles (visual_ebullition==F), check what flux is assigned to them (LM vs total.flux).

ch4_falsepositives<- CH4all_nosep_flux.auto %>% 
  filter(!UniqueID%in%incub_discard) %>% #Not wrong incubations
  filter(UniqueID%in%c(incub_ebu_auto)) %>% #Detected bubbles
  filter(!UniqueID%in%c(incub_ebu_visual)) %>% #Not visual ebullition
  mutate(prefered_model=if_else(LM.r2>0.95,"LM","total.flux"),
         linear_total_ratio=LM.flux/full.total.flux)

ch4_falsepositives %>% group_by(prefered_model) %>% 
  summarise(numcases=n(),
            avg_LMtotalratio=mean(linear_total_ratio))

#Most flase-positives are assigned linear flux estimate, with good correspondence LM-to-totalflux
ch4_falsepositives %>%
  ggplot(aes(x=prefered_model, y=linear_total_ratio, col=prefered_model))+
  geom_boxplot()+
  geom_point()

#Check the cases with total.flux assigned: total.flux is more appropriate than linear fit, flat-ish with a lot of noise in the middle
falsepos_totalflux_inspected<- c("s2-cu-p1-5-o-d-08:30","s2-cu-p1-6-v-t-08:39","s2-ri-a2-13-o-d-11:24","s2-va-r2-3-v-d-10:26","s3-cu-p2-3-v-d-09:06","s3-du-r2-5-b-d-09:44","s4-ca-a1-8-v-t-10:52","s4-ca-a1-8-v-d-11:00","s4-ca-a1-12-v-t-12:48","s4-cu-p2-9-o-d-08:43","s4-cu-r1-2-b-d-07:13","s4-du-a2-13-b-d-12:30","s4-du-a2-14-v-t-12:41","s4-du-a2-14-v-d-12:47","s4-du-a2-15-v-d-13:04","s4-ri-p1-9-o-d-08:56","s4-ri-p1-10-v-d-09:15","s4-ri-p1-13-o-d-10:07","s4-va-p1-12-v-d-10:18","s4-va-r1-13-v-d-12:03")

ch4_falsepositives %>%
  filter(prefered_model=="total.flux") %>% 
  filter(!UniqueID%in%falsepos_totalflux_inspected) %>% 
  select(UniqueID)








ch4_ebullitive<- CH4all_nosep_flux.auto %>% 
  filter(!UniqueID%in%incub_discard) %>% #Not wrong incubations
  filter(UniqueID%in%c(incub_ebu_visual,incub_ebu_auto)) %>% #Contains ebullition (auto or visual)
  mutate(prefered_model=if_else(LM.r2>0.95,"LM","total.flux"),
         linear_total_ratio=LM.flux/full.total.flux)


ch4_ebullitive %>% 
  filter(prefered_model=="LM") %>% 
  ggplot(aes(x=LM.r2, y=linear_total_ratio, col=between(linear_total_ratio,0.9,1.1)))+
  geom_point()

#-----falta de adaptar----








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







#2. Classify CH4w incubations----

  #make sure the order between the two datasets is the same and match UniqueIDs
dataflags<- CH4w_nosep_flux.auto %>% 
  select(UniqueID, nb.obs) %>%
  rename(total_nb.obs=nb.obs) %>% 
  merge.data.frame(CH4w_sep_flux.auto %>% 
                     select(UniqueID, nb.obs, full.total.flux, full.ebullition.flux, diffusion.flux, diffusion.flux.SD) %>% 
                     rename(difusion_nb.obs=nb.obs), by="UniqueID") %>% 
  merge.data.frame(auxfile_inspected %>% select(UniqueID, ch4_decission), by="UniqueID") %>% 
  #0. Is incubation correct? (i.e. is not marked for discard)
  mutate(good_incubation=ch4_decission=="ok") %>% 
  #1. Is ebullition detected? (use difference in nb.obs between nosep and sep results)
  mutate(buble_detected=!(total_nb.obs==difusion_nb.obs)) %>% 
  #2. Is ebullition positive? i.e. diffusion below total flux? (diffusion.flux<full.total.flux)
  mutate(positive_ebullition=diffusion.flux<full.total.flux) %>% 
  #3. Is separation_appropriate? (combination of detected and positive bubles)
  mutate(good_separation=buble_detected&positive_ebullition) %>% 
  #4. Is separation significant? (full.ebullition.flux > diffusion.flux +diffusion.flux.SD)
  mutate(sig_separation=full.ebullition.flux> diffusion.flux+diffusion.flux.SD)


#Separate UniqueIDs based on dataflags
#Get wrong incubations (best.flux <- NA, quality.check<- "Wrong incubation, no reliable flux")
wrong_incubations<- dataflags %>% filter(good_incubation==F) %>% pull(UniqueID) 

#Get incubations to set goflux criteria for difusion only
difusive_incubations <- dataflags %>% filter(good_incubation==T&buble_detected==F) %>% pull(UniqueID) 

#Get incubations for which LM fit threshold is needed (those with non-significant separation and those with negative ebullition)
lm_or_totalflux_incubations <- dataflags %>% filter(good_incubation==T&buble_detected==T & (sig_separation==F|positive_ebullition==F)) %>% pull(UniqueID)

#Get incubations for which flux separation works ok (buble detected, buble positive, significant separation)
ebullitive_incubations<- dataflags %>% filter(good_incubation==T&good_separation==T&sig_separation==T) %>% pull(UniqueID)


#CHECKS: 
#Check, any incubation not classified?
dataflags %>% 
  filter(!(UniqueID%in%c(wrong_incubations,difusive_incubations,lm_or_totalflux_incubations,ebullitive_incubations)))

#Check any incubation in two classifications?
anyDuplicated(c(wrong_incubations,difusive_incubations,lm_or_totalflux_incubations,ebullitive_incubations))


#3. LM.R2 limit for buble-containing-----
#Set the LM fit threshold for incubations contatined in lm_or_totalflux_incubations (those with wrong or non-significant separation)

#Use visual inspection descriptions in auxfile_inspected

ch4_lmortotal<- CH4w_nosep_flux.auto %>% 
  filter(UniqueID%in%lm_or_totalflux_incubations) %>% 
  merge.data.frame(auxfile_inspected %>% select(UniqueID, visual_ebullition, method_sucess, Common_obs), by="UniqueID") %>% 
  #add full.total.flux estimate (to compare missmatch LM.flux vs full.total.flux)
  merge.data.frame(CH4w_sep_flux.auto %>% select(UniqueID, full.total.flux, SD_full.total.flux), by="UniqueID") %>% 
  mutate(visual_recommendation=case_when(grepl("linear fit", Common_obs)~"LM",
                                         grepl("totalflux", Common_obs)~"total_flux",
                                         TRUE~NA)) %>% 
  #ratio_linear_total
  mutate(ratio_linear_total=LM.flux/full.total.flux)
  

ch4_lmortotal %>% filter(!is.na(visual_recommendation)) %>% 
  ggplot(aes(y=LM.r2, fill=visual_recommendation))+
  geom_histogram()+
  facet_wrap(~visual_recommendation)

ch4_lmortotal %>% filter(!is.na(visual_recommendation)) %>% 
  ggplot(aes(y=LM.AICc, fill=visual_recommendation))+
  geom_histogram()+
  facet_wrap(~visual_recommendation)

ch4_lmortotal %>% filter(!is.na(visual_recommendation)) %>% 
  ggplot(aes(x=LM.r2, fill=visual_recommendation))+
  geom_histogram()+
  geom_vline(xintercept = 0.95)+
  facet_wrap(~visual_recommendation)

ch4_lmortotal %>% filter(!is.na(visual_recommendation)) %>% 
  ggplot(aes(x=LM.r2, y=ratio_linear_total, col=visual_recommendation))+
  geom_point()


#Any linear visual decissions with LMR2 < 0.95 and more than 10% difference in estimate (LM or total_flux)(would totalflux be apropriate?)
ch4_lmortotal %>% filter(visual_recommendation=="LM"&LM.r2<0.95) %>%
  select(UniqueID, LM.r2, ratio_linear_total, LM.flux, full.total.flux) %>% 
  filter(abs(ratio_linear_total-1)>0.1)


#Any ebullitive visual decissions with LMR2>=0.95 and more than 10% difference in estimate to inspect
ch4_lmortotal %>% filter(visual_recommendation=="total_flux"&LM.r2>=0.95) %>%
  select(UniqueID, LM.r2, ratio_linear_total, LM.flux, full.total.flux) %>% 
  filter(abs(ratio_linear_total-1)>0.1)
#s3-cu-p2-5-o-d-10:02 linear ok-ish
#s4-cu-a1-9-o-d-09:48 linear ok-ish
#s4-cu-p2-3-o-d-07:44 linear ok-ish
#s4-cu-p2-7-o-d-08:17 linear ok-ish
#all incubations could be represented well by a linear model


#Re-check very high difference in estimates
ch4_lmortotal %>% 
  mutate(selection=if_else(LM.r2>=0.95, "LM","total_flux"),
         agree_model=selection==visual_recommendation) %>% 
  select(UniqueID, LM.r2, ratio_linear_total, LM.flux, full.total.flux,visual_recommendation, agree_model) %>%
  filter(abs(ratio_linear_total-1)>0.5)
#Re-Checked, visual inspection ok: all total flux, differences caused by leverage of bubles on LM, affecting slope.  


#Decission: hard threshold to select LM or total_flux as best.flux is set to LM.r2 0.95: best.flux=if_else(LM.r2>=0.95, "LM","total_flux")









# 4.Criteria for diffusive water incubations: -----
#Follow the same logic as CO2: 
#Must have the exact same thresholds of dryCH4 incubations
#1. Copy-paste approach from CO2 criteria, 
#2. join dataset with ch4dry
#3. Explore threshold options and implications, with special attention to not-detected ebullition cases: this should never be fitted with HM model (to avoid overestimation of fluxes, fitting a strong curve to a buble-flat type of incubation). Potentially, same approach than for dry-ebullitive patterns (High-curvature gfact > threshold --> LM or total-flux, depending on which provides less uncertainty. 

difusivech4<- CH4w_nosep_flux.auto %>% 
  filter(UniqueID%in%difusive_incubations) %>% 
  #Add comments from visual inspection: 
  merge.data.frame(auxfile_inspected %>% select(UniqueID, visual_ebullition, method_sucess, Common_obs), by="UniqueID") %>% 
  #add full.total.flux estimate (to compare missmatch LM.flux vs full.total.flux)
  merge.data.frame(CH4w_sep_flux.auto %>% select(UniqueID, full.total.flux, SD_full.total.flux), by="UniqueID") %>% 
  mutate(visual_recommendation=case_when(grepl("linear fit", Common_obs)~"LM",
                                         grepl("totalflux", Common_obs)~"total_flux",
                                         TRUE~NA)) %>% 
  #ratio_linear_total
  mutate(ratio_linear_total=LM.flux/full.total.flux)


#Check general parameters of method_sucess==F (non-detected ebullition) cases
difusivech4 %>% 
  filter(method_sucess==F) %>% 
  ggplot(aes(x=LM.r2, y=HM.r2, col=Common_obs))+
  geom_point()

difusivech4 %>% 
  filter(method_sucess==F) %>% 
  ggplot(aes(x=g.fact, y=HM.r2, col=Common_obs))+
  geom_point()+
  geom_vline(xintercept = 4)



unique(difusivech4$Common_obs)

difusivech4 %>%
  ggplot(aes(y=LM.r2, fill=Common_obs))+
  geom_histogram()+
  facet_wrap(~Common_obs)


difusivech4 %>%
  ggplot(aes(x=g.fact, y=ratio_linear_total, col=Common_obs))+
  geom_point()



#THINGS to decide:------
#WHAT TO DO WITH NON-detected ebullition when LM and HM are not appropriate? 
#Potentially, calculate total.flux for all CH4 incubations, and set thresholds for goflux fit under which we use define best flux as total flux

#Check plots and define when the above logic fails:

#1. Cases of visual ebullition not detected algorithm (buble_detected==F) for which LM/HM are not appropriate. 

#2. Cases with good_separation==T but non-reliable diffusion (difusion contains buble), including those with significant separation. 

#What to do with negative ebullition cases: either LM or full.total.flux  

#Save to inspect missmatches and define errors.  
write.csv(inspection_table, paste0(results_path,"Inspection_table_criteria.csv"),row.names = F)


