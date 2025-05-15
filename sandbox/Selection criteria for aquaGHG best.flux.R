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






#1. Load results------

#Load results for water incubations with flux separation: 
load(paste0(results_path, "CH4w_sep_results_aquaGHG.Rdata"))

#Load results for water incubations withOUT flux separation: 
load(paste0(results_path, "CH4w_nosep_results_aquaGHG.Rdata"))

#Load visual inspection flags for CH4 water: 
auxfile_inspected<- read.csv(paste0(results_path, "ch4_water_auxfile4inspection_filled.csv"))


#Load discard decissions from auxfile ch4
full_auxfilech4<- read.csv(file = paste0(auxfile_path,"ch4_auxfile.csv"))
#Load results for ch4 dry incubations: 
load(paste0(results_path, "CH4dry_results_aquaGHG.Rdata"))


#2. Criteria for dryCH4 -----



##Discard ch4dry------

#Get wrong incubations: those with too much influence of artefacts or wrong manipulations based on visual inspection of time-series (recorded in auxfile).

ch4dry_discard<- full_auxfilech4 %>% filter(water_depth>0&ch4_decission=="discard") %>% pull(UniqueID)

#these incubations will be flagged as unreliable in the final dataset and will not influence the criteria decisions. 



#Criteria Hard.threshold checks

##Inspect mdf----
#We want to make sure MDF fluxes do not arise from artefacts, all MDF incubations will default to LM 
ch4dry_mdf<- CH4dry_flux.auto %>%  
  filter(!UniqueID%in%ch4dry_discard) %>%
  filter(abs(LM.flux)<=MDF.lim)

#POr aqui------
#Need to re-run goflux per-ghg-cropped to produce plots to be able to inspect final incubations (those with final cropping decissions)
#MDFs CH4 inspected: true below detection flat-ish white noise incubations
mdf_inspected<- c("s1-ca-p2-15-v-t-14:54","s1-ca-p2-4-v-d-11:24","s1-ca-p2-4-v-t-11:17","s1-ca-p2-6-b-d-11:53","s1-ca-r1-17-b-d-12:49","s1-ca-r2-1-b-d-08:47","s1-da-a1-1-b-d-08:22","s1-da-a1-2-v-d-08:46","s1-da-a1-2-v-t-08:38","s1-da-a1-5-v-d-10:01","s1-da-a1-9-b-d-12:23","s1-du-a1-11-b-t-13:12")

#"s1-ca-r1-16-v-d-12:41" cropstart 60s
#"s1-ca-r1-9-v-t-10:57" cropstart 20s
#"s1-da-a1-7-b-d-11:51" cropstart 40s (co2&ch4)


#any MDF flux not inspected?
ch4dry_mdf %>% filter(!UniqueID%in%mdf_inspected) %>% pull(UniqueID)


#-----falta de adaptar----
#falta de adaptar desde aqui hasta 2.classify ch4w incubations.

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
co2_inspected_aicc<- c("s1-da-p1-11-o-d-10:46","s1-du-r2-5-b-d-08:32","s1-du-r2-8-v-d-09:51","s2-ca-r2-3-o-d-09:09","s2-cu-a1-14-o-d-10:48","s2-cu-a1-4-v-d-09:05","s2-cu-a1-7-o-d-09:37","s2-cu-a2-2-v-d-08:32","s2-cu-a2-3-v-t-08:47","s2-da-r1-11-o-d-11:15","s2-da-r1-9-v-d-10:56","s2-du-a2-10-b-d-12:03","s2-ri-a2-10-b-d-10:48","s2-ri-p1-6-v-t-10:17","s2-ri-p2-11-v-d-12:24","s2-va-a2-14-v-t-12:57","s2-va-r1-9-v-d-11:42","s2-va-r2-14-o-d-12:59","s2-va-r2-7-v-t-11:08","s3-ca-p2-1-v-t-07:40","s3-ca-p2-5-v-t-09:33","s3-da-a1-2-b-d-07:34","s3-da-a2-10-o-d-08:57","s3-ri-a2-2-b-d-09:06","s4-da-a1-8-v-t-09:40","s4-da-r1-2-v-d-09:03","s4-da-r2-8-o-d-08:12","s4-ri-a2-12-b-d-09:02","s4-ri-a2-9-b-d-08:27","s2-va-r2-3-v-t-10:19")

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
#3. Explore threshold options and implications, with special attention to not-detected ebullition cases: this should never be fitted with HM model (to avoid overestimation of fluxes, fitting a hard curve to a buble-flat type of incubation)


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


