
# ---
# Authors: Miguel Cabrera
# Project: "RESTORE4Cs"
# date: "Feb 2025"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script

# This script analyses the quality of fit of all fluxes. First by ranking and flagging incubations based on regression fit. Then by flagging incubations that include artefacts (examining deltaGHG and flagging incubations with more than 5% of deltaGHG==0)
#It saves all flags and quality issues based on the 5s margin incubation run to inspect every incubation and save cropping decissions for final flux calculation

#LAST update on inspection files (including cropping decisions): 8/5/2025

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
library(ggpubr)

repo_root <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
files.sources = list.files(path = paste0(repo_root,"/functions"), full.names = T)
for (f in files.sources){source(f)}


# ---- Directories and data loading ----
dropbox_root <- "C:/Users/Miguel/Dropbox/RESTORE4Cs - Fieldwork/Data" # You have to make sure this is pointing to the write folder on your local machine

# results_path <- paste0(dropbox_root,"/GHG/Processed data/computed_flux/")

# plots_path <- paste0(dropbox_root,"/GHG/Processed data/computed_flux/plots")
results_path<- "C:/Users/Miguel/Dropbox/TEST_quality_raw2flux/5smargin fluxes/"

fieldsheetpath <- paste0(dropbox_root,"/GHG/Fieldsheets")

quality_path<- paste0(dropbox_root, "/GHG/Working data/Incubation_quality/") #to save quality assesment and flags.



#Load csv files produced by raw2flux_miguel_edits_5smargin.R. script
#Last run with 5second crop of start-times and end-times (croped after all map_incubations corrections). 

# setwd(results_path)

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
table_ch4 <- get_full_detailed_table(table_results_all, variable = "ch4")



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
  geom_vline(xintercept = 0.9)

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





