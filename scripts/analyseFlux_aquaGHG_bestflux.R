
# ---
# Authors: Miguel Cabrera, Camille Minaudo
# Project: "RESTORE4Cs"
# date: "June 2025"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script----
# This script loads the fluxes estimates from Selection criteria for aquaGHG best.flux.R script
# it is used to explore and analyze the dataset.

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

repo_root <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
files.sources = list.files(path = paste0(repo_root,"/functions"), full.names = T)
for (f in files.sources){source(f)}


# ---- Directories ----
dropbox_root <- "C:/Users/Miguel/Dropbox/RESTORE4Cs - Fieldwork/Data" # You have to make sure this is pointing to the write folder on your local machine
#Path to Co2 and Ch4 auxfiles:
auxfile_path<- paste0(dropbox_root,"/GHG/Working data/Auxfiles_4aquaGHG/") 
#Path to aquaGHG best.flux results:
results_path <- paste0(dropbox_root,"/GHG/Processed data/computed_flux/")
#Path to save plots from exploration: 
plots_path <- paste0(dropbox_root,"/GHG/Processed data/computed_flux/plots")
#Path with fieldsheets:
fieldsheet_path <- paste0(dropbox_root,"/GHG/Fieldsheets")
#Path with vegetation data: 
vegetation_path<- paste0(dropbox_root, "/Vegetation/")


#----Import fluxes-----
#1. Import bestflux results
co2_best<- read.csv(paste0(results_path,"co2_bestflux.csv"))
ch4_best<- read.csv(paste0(results_path,"ch4_bestflux.csv"))



#Import meta-data------
#Read all GHGfieldsheets (using function from Camille)
fieldsheets<- list.files(path= fieldsheet_path, recursive = T, pattern = "GHG.xlsx",full.names = T)
field<- read_GHG_fieldsheets(fieldsheets)

#Import harmnonized vegetation
veg_data<- read.csv(paste0(vegetation_path,"RESTORE4Cs_finalvegetation_ABG_biomass_description.csv"))

#Join and keep only relevant info:
field_veg<- field %>% 
  rename(UniqueID=uniqID) %>% 
  separate(UniqueID, into = c("c1","c2","c3","c4","c5","c6","c7"), sep = "-", remove = F) %>% 
  mutate(plotcode=toupper(paste(c1,c2,c3,c4,sep = "-"))) %>% 
  merge.data.frame(veg_data %>% select(plotcode, vegetation_description), by="plotcode", all=T) %>% 
  select(UniqueID,pilot_site,subsite,plot_id,plotcode,strata,transparent_dark,water_depth,comments, vegetation_description) %>% 
  rename(sampling=subsite, og_strata=strata) %>% 
  mutate(season=substr(sampling, 1,2), 
         subsite=substr(sampling,4,8),
         strata=case_when(og_strata=="vegetated"&water_depth>0~"vegetated water",
                          og_strata=="vegetated"&water_depth<=0~"vegetated dry",
                          TRUE~og_strata)) %>% 
  select(UniqueID, pilot_site, subsite, sampling, plot_id, plotcode, og_strata, strata, water_depth, transparent_dark,vegetation_description, comments)




#Join data-----

table_co2<- merge.data.frame(co2_best,field_veg, by="UniqueID", all.x = T) %>% 
  #Remove Unreliable fluxes 
  filter(!grepl("Discard", best.model.flags)) %>% 
  #Add column best.flux
  mutate(best.flux=case_when(best.model=="LM"~LM.flux,
                             best.model=="HM"~HM.flux),
         best.flux.se=case_when(best.model=="LM"~LM.flux.se,
                                best.model=="HM"~HM.flux.se),
         best.model.MAE=case_when(best.model=="LM"~LM.MAE,
                                  best.model=="HM"~HM.MAE)
         )

table_ch4<- merge.data.frame(ch4_best,field_veg, by="UniqueID", all.x = T) %>% 
  #Remove Unreliable fluxes 
  filter(!grepl("Discard", best.model.flags)) %>% 
  #Add column best.flux
  mutate(best.flux=case_when(best.model=="LM"~LM.flux,
                             best.model=="HM"~HM.flux,
                             best.model=="total.flux"~total.flux),
         best.flux.se=case_when(best.model=="LM"~LM.flux.se,
                                best.model=="HM"~HM.flux.se,
                                best.model=="total.flux"~total.flux.se))


#Separate ch4 from difusive incubations
table_ch4_difusive<- table_ch4 %>% filter(!grepl("Ebullitive", best.model.flags)) %>% 
  mutate( best.model.MAE=case_when(best.model=="LM"~LM.MAE,
                                   best.model=="HM"~HM.MAE))


#Inspect C0.10s-----

#Inpect and flag incubations that start with high CH4 concentrations (unreliable fluxes)
#Fluxes are expected to be under-estimated when incubation starts with very high concentrations (due to lack of ventilation) this pattern should presumably only affect HM/LM fluxes in a relevant way. Fluxes driven by ebullition, should not suffer this under-estimation too much. 

table_ch4_difusive %>% 
  ggplot(aes(x=best.flux, y=C0.10s, col=best.model))+
  geom_point()+
  geom_errorbar(aes(ymax = C0.10s+C0.10s.sd, ymin=C0.10s-C0.10s.sd))


table_ch4_difusive %>% 
  filter(C0.10s<3000) %>% 
  ggplot(aes(x=C0.10s))+
  geom_histogram(bins = 100)+
  scale_x_continuous(breaks = seq(1500,12000,by=100))

#A few fluxes with starting CH4 concentration above 2500: make sure they do not come from cropped incubations: 

cropdecissions<- read_xlsx(paste0(dropbox_root,"/GHG/Working data/Incubation_quality/Inspection_table_allincubations_tocrop_perGHG.xlsx")) %>% select(UniqueID,co2_start_crop_s,ch4_start_crop_s, ch4_end_crop_s)
####FINAL CHECKS and RE-calculations-----

#ALL re-calculated
uncrop<- c()
           

#artefacts cropped, start high ok
ok<- c("s1-va-a2-15-v-d-13:46", "s1-va-a2-8-v-d-11:13", "s2-ca-a2-2-v-d-08:33",  "s2-ca-a2-2-v-t-08:25","s2-ca-p2-15-b-t-13:23", "s2-ca-p2-14-b-t-13:13", "s2-ri-r2-12-v-d-12:40", "s2-va-a2-8-v-d-11:13", "s3-va-a2-12-o-d-11:05","s3-va-a2-5-v-d-09:17", "s3-va-p1-9-v-t-10:42", "s4-da-a2-1-o-d-06:31","s4-va-a1-1-o-d-08:30","s4-va-a2-2-o-d-07:21", "s4-va-p2-2-b-d-07:12","s1-ca-a2-7-v-t-11:58","s1-ri-p1-1-v-d-10:19","s2-du-a2-3-v-t-09:34","s2-du-a2-4-v-d-09:57","s2-du-a2-5-v-d-10:24","s2-du-a2-5-v-t-10:16","s2-ri-a2-2-b-d-08:51","s2-ri-r2-4-v-t-09:15","s3-va-a2-5-v-t-09:09","s4-ca-a2-6-v-t-08:55","s4-va-a2-1-b-d-07:13","s1-du-a2-15-b-t-13:05","s2-ca-a1-1-o-d-08:47","s2-ca-a1-2-o-d-08:59","s2-ca-a2-10-v-d-10:58","s2-ca-a2-4-b-d-09:01","s2-ca-r2-5-o-d-09:42","s2-du-a2-6-v-d-10:40","s2-du-r1-3-v-d-09:59","s2-du-r1-4-o-d-10:07","s2-du-r1-5-o-d-10:23","s2-du-r1-7-v-t-10:41","s2-du-r1-9-v-t-11:15","s2-ri-a2-14-o-d-11:39","s3-du-p2-1-o-d-08:05","s3-va-a2-13-v-t-11:22","s4-ca-a2-4-v-t-08:25","s4-da-r1-3-v-t-09:25","s1-cu-p2-13-o-d-12:37","s1-du-p1-3-o-d-08:43","s2-du-a2-9-v-d-11:46","s2-du-r1-14-b-d-12:46","s2-du-r1-15-v-d-13:01","s2-ri-p2-3-v-d-09:15","s2-ri-p2-7-v-d-10:29"
       )

table_ch4_difusive %>% 
  filter(C0.10s>2250) %>% 
  merge.data.frame(cropdecissions, by="UniqueID", all.x = T) %>% 
  mutate(croped=if_else(!is.na(ch4_start_crop_s),T,F)) %>% 
  filter(croped) %>% 
  filter(!UniqueID%in%c(uncrop, ok)) %>% 
  filter(co2_start_crop_s!=ch4_start_crop_s) %>% 
  select(UniqueID, C0.10s, best.model, g.fact,ch4_start_crop_s,ch4_end_crop_s) %>% head()




table_ch4 %>% 
  ggplot(aes(x=C0.10s, y=Cf.10s, col=best.model))+
  geom_point()+
  geom_abline(slope=1, intercept = 0)+
  facet_wrap(~C0.10s>5000,scales="free")


  ggplot(aes(x=C0.10s, y=total.flux.se/total.flux, col=best.model))+
  geom_point()+
  geom_abline(slope=1, intercept = 0)+
  facet_wrap(~C0.10s>5000,scales="free")


table_ch4 %>% 
  ggplot(aes(x=C0.10s, y=C0.10s.sd/C0.10s, col=best.model))+
  geom_point()+
  facet_wrap(~C0.10s>3000,scales="free")


#How many plots will be left without any CH4flux at different thresholds of max C0.10s?

# Define thresholds
thresholds <- seq(2500, 10000, by = 500)

# Apply logic for each threshold, keeping full data
annotated <- lapply(thresholds, function(thresh) {
  table_ch4 %>%
    group_by(plotcode) %>%
    mutate(all_above = all(C0.10s > thresh)) %>%
    ungroup() %>%
    mutate(threshold = thresh)
}) %>%
  bind_rows()

annotated %>% 
  group_by(threshold) %>% 
  summarise(count=sum(all_above)) %>% 
  ggplot(aes(x=threshold, y=count))+
  geom_point()


  



#CO2 overview----

n <- dim(table_co2)[1]
ind_lm <- which(table_co2$best.model=="LM")
in_hm <- which(table_co2$best.model=="HM")

message("CO2 Fluxes calculated for ",n, " different incubations.")
message(" --- ",dim(table_co2[table_co2$og_strata=="bare",])[1] ," in bare ")
message(" --- ",dim(table_co2[table_co2$og_strata=="open water",])[1] ," in open water ")
message(" --- ",dim(table_co2[table_co2$og_strata=="vegetated",])[1] ," in vegetated")


message("CO2 Best model is non-linear for ",round(length(in_hm)/n*100*100)/100,"% of the measurements")



#Flags CO2 ----

table.flags_co2 <- table_co2 %>%
  filter(best.model.flags != "") %>%
  pull(best.model.flags) %>%
  strsplit(" \\| ") %>%
  unlist() %>%
  trimws() %>%
  table() %>%
  as.data.frame()

names(table.flags_co2) <- c("flag", "n")

table.flags_co2 <- table.flags_co2 %>%
  mutate(perc = round(n / nrow(table_co2) * 100, 2)) %>%
  arrange(desc(n))

table.flags_co2


table_co2_long <- table_co2 %>%
  filter(!is.na(best.model.flags), best.model.flags != "") %>%
  mutate(row_id = row_number()) %>%  # Unique ID to rejoin later
  separate_rows(best.model.flags, sep = " \\| ") %>%
  mutate(best.model.flags = trimws(best.model.flags))

# Step 2: Plot using separated flags
ggplot(table_co2_long, 
       aes(x = best.model.flags, 
           y = (HM.flux - LM.flux)/HM.flux, 
           fill = best.model)) +
  geom_hline(yintercept = c(-1, 0, 1), color = "grey70") +
  geom_jitter(alpha = 0.1, size = 2, width = 0.1) +
  geom_violin(alpha = 0.2, scale = "width", draw_quantiles = c(0.5), 
              aes(colour = best.model, fill = best.model)) +
  theme_article() +  # Assuming you've defined or loaded this theme
  ylim(c(-1, 2)) +
  scale_fill_viridis_d(option = "A", end = 0.8, direction = -1) +
  scale_colour_viridis_d(option = "A", end = 0.8, direction = -1) +
  coord_flip()


#Model flux diff----
table_co2$diff_model <- (table_co2$HM.flux-table_co2$LM.flux)/table_co2$HM.flux

n_less_1perc <- length(which(abs(table_co2$diff_model)<0.01))

message("Linear and non-linear models have less than 1% difference for ",round(n_less_1perc/n*100*100)/100,"% of the measurements")

n_less_10perc<- length(which(abs(table_co2$diff_model)<0.1))
ggplot(table_co2[order((table_co2$diff_model)),], aes(seq(1,n)/n*100, (diff_model)*100))+
  geom_vline(xintercept = seq(0,100,50), color = "grey70")+
  geom_hline(yintercept = c(-100,-10,0,10,100), color = "grey70")+
  geom_point()+
  # scale_y_log10()+
  ylim(c(-50,120))+
  xlab("Proportion of timeseries")+
  ylab("Relative difference [% of HM flux]")+
  theme_article()+ggtitle(paste0("CO2 HM-LM models difference is below 10% for ",round(n_less_10perc/n*100*10)/10,"% of the measurements"))


ggplot(table_co2, aes(abs(best.flux), abs(diff_model)*100))+
  # geom_vline(xintercept = seq(0,100,50), color = "grey70")+
  geom_hline(yintercept = c(1,10,100), color = "grey70")+
  geom_point(aes(col=best.model))+
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


#adapt Ch4 model dif----
table_ch4_difusive$diff_model <- (table_ch4_difusive$HM.flux-table_ch4_difusive$LM.flux)/table_ch4_difusive$HM.flux

df_exceedance_ch4 <- NULL
for(thresh in c(0.1,0.5,seq(1,100))){
  n_less <- length(which(abs(table_ch4_difusive$diff_model)<thresh/100))
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
  ggtitle(paste0("For difusive CH4, HM-LM models difference is below 10% for ",round(df_exceedance_ch4$p[df_exceedance_ch4$t==10]),"% of the measurements"))




p_diff_models <- ggarrange(p_co2, p_ch4, ncol = 1)
ggsave(plot = p_diff_models, filename = "HM-LM models difference.jpeg", path = plots_path, 
       width = 7, height = 5, dpi = 300, units = 'in', scale = 1.1)




ggplot(table_co2, aes(HM.MAE, diff_model))+geom_point(aes(colour = best.model))+
  scale_colour_viridis_d("best model", option = "A", end = 0.8, direction = -1)+
  theme_article()+ 
  scale_x_log10()+  scale_y_log10()


ggplot(table_co2)+
  geom_point(aes(x=abs(best.flux), y= best.model.MAE, colour = best.model), alpha=0.2, size=3)+
  theme_article()+
  # xlim(c(0,50))+
  scale_x_log10()+
  scale_y_log10()+
  ylab("best.model.MAE")+
  scale_colour_viridis_d("model", option = "A", end = 0.8, direction = -1)


ggplot(table_co2, aes(LM.MAE, HM.MAE, shape=best.model,colour = best.model))+
  geom_point(alpha=0.4, size=3)+
  theme_article()+
  scale_x_log10()+
  scale_y_log10()+
  scale_colour_viridis_d("best model", option = "A", end = 0.8, direction = -1)


ggplot(table_co2)+
  geom_density(aes(best.model.MAE, fill = best.model), alpha=0.5)+
  theme_article()+
  scale_x_log10()



ggplot(table_ch4_difusive, aes(HM.MAE, diff_model))+geom_point(aes(colour = best.model))+
  scale_colour_viridis_d("best model", option = "A", end = 0.8, direction = -1)+
  theme_article()+ 
  scale_x_log10()+  scale_y_log10()


ggplot(table_ch4_difusive)+
  geom_point(aes(x=abs(best.flux), y= best.model.MAE, colour = best.model), alpha=0.2, size=3)+
  theme_article()+
  # xlim(c(0,50))+
  scale_x_log10()+
  scale_y_log10()+
  ylab("best.model.MAE")+
  scale_colour_viridis_d("model", option = "A", end = 0.8, direction = -1)


ggplot(table_ch4_difusive, aes(LM.MAE, HM.MAE, shape=best.model,colour = best.model))+
  geom_point(alpha=0.4, size=3)+
  theme_article()+
  scale_x_log10()+
  scale_y_log10()+
  scale_colour_viridis_d("best model", option = "A", end = 0.8, direction = -1)


ggplot(table_ch4_difusive)+
  geom_density(aes(best.model.MAE, fill = best.model), alpha=0.5)+
  theme_article()+
  scale_x_log10()

#There are two clearly defined groups of incubations according to their best.model.MAE. This is related to instrument precision: 
ggplot(table_ch4_difusive)+
  geom_density(aes(best.model.MAE, fill = best.model), alpha=0.5)+
  theme_article()+
  scale_x_log10()+
  facet_grid(gas_analiser~.)



# ---- plot CO2 and CH4 flux across sites and seasons ----


p_overview_co2 <- ggplot(table_co2, aes(lightCondition, best.flux, colour = lightCondition))+
  geom_hline(yintercept = 0)+
  geom_jitter(alpha=0.5)+
  # geom_boxplot(alpha=0.8)+
  theme_article()+
  ylab("CO2 flux umol/m2/s")+
  scale_fill_viridis_d(begin = 0.2, end = 0.9)+
  scale_colour_viridis_d(begin = 0.2, end = 0.9, option = "C")+
  # facet_wrap(strata~.)+
  theme(legend.position = 'none')+xlab("")+
  ggtitle(paste0("CO2 fluxes, ",length(unique(table_co2$UniqueID)), " valid incubations"))


p_overview_ch4 <- ggplot(table_ch4, aes(lightCondition, best.flux, colour = lightCondition))+
  geom_hline(yintercept = 0)+
  geom_jitter(alpha=0.5)+
  # geom_boxplot(alpha=0.8)+
  theme_article()+
  ylab("CH4 flux nmol/m2/s")+
  scale_fill_viridis_d(begin = 0.2, end = 0.9)+
  scale_colour_viridis_d(begin = 0.2, end = 0.9, option = "C")+
  # facet_wrap(strata~.)+
  theme(legend.position = 'right')+xlab("")+
  ggtitle(paste0("CH4 fluxes, ",length(unique(table_ch4$UniqueID)), " valid incubations"))

p_overview_co2_ch4 <- ggarrange(p_overview_co2, p_overview_ch4, nrow = 1)

myfilename <- paste("000_overview_CO2_CH4_lightCond",paste(unique(table_ch4$sampling), collapse = "_"), sep = "_")
ggsave(plot = p_overview_co2_ch4, filename = paste0(myfilename,".jpg"), path = plots_path,
       width = 10, height = 4, dpi = 300, units = 'in', scale = 1.1)






p_overview_co2 <- ggplot(table_co2, aes(strata, best.flux, colour = lightCondition))+
  geom_hline(yintercept = 0)+
  geom_jitter(alpha=0.5)+
  # geom_boxplot(alpha=0.8)+
  theme_article()+
  ylab("CO2 flux umol/m2/s")+
  scale_fill_viridis_d(begin = 0.2, end = 0.9)+
  scale_colour_viridis_d(begin = 0.2, end = 0.9, option = "C")+
  # facet_wrap(sampling~.)+
  theme(legend.position = 'none')+xlab("")

p_overview_ch4 <- ggplot(table_ch4, aes(strata, best.flux, colour = lightCondition))+
  geom_hline(yintercept = 0)+
  geom_jitter(alpha=0.5)+
  # geom_boxplot(alpha=0.8)+
  theme_article()+
  ylab("CH4 flux nmol/m2/s")+
  scale_fill_viridis_d(begin = 0.2, end = 0.9)+
  scale_colour_viridis_d(begin = 0.2, end = 0.9, option = "C")+
  # facet_wrap(strata~.)+
  theme(legend.position = 'right')+xlab("")
p_overview_co2_ch4 <- ggarrange(p_overview_co2, p_overview_ch4, nrow = 1)

myfilename <- paste("000_overview_CO2_CH4_strata",paste(unique(table_ch4$sampling), collapse = "_"), sep = "_")
ggsave(plot = p_overview_co2_ch4, filename = paste0(myfilename,".jpg"), path = plots_path,
       width = 10, height = 4, dpi = 300, units = 'in', scale = 1.1)



p_overview_co2 <- ggplot(table_co2, aes(sampling, best.flux))+
  geom_hline(yintercept = 0)+
  geom_jitter(alpha=0.4, aes(colour = lightCondition))+
  geom_boxplot(alpha=0.5, width=0.3)+
  theme_article()+
  ylab("CO2 flux umol/m2/s")+
  scale_fill_viridis_d(begin = 0.2, end = 0.9)+
  scale_colour_viridis_d(begin = 0.2, end = 0.9, option = "C")+
  # facet_wrap(sampling~.)+
  theme(legend.position = 'none')+xlab("")

p_overview_ch4 <- ggplot(table_ch4, aes(sampling, best.flux))+
  geom_hline(yintercept = 0)+
  geom_jitter(alpha=0.4)+
  geom_boxplot(alpha=0.5, width=0.3)+
  theme_article()+
  ylab("CH4 flux nmol/m2/s")+
  scale_fill_viridis_d(begin = 0.2, end = 0.9)+
  scale_colour_viridis_d(begin = 0.2, end = 0.9, option = "C")+
  # facet_wrap(strata~.)+
  theme(legend.position = 'right')+xlab("")
p_overview_co2_ch4 <- ggarrange(p_overview_co2, p_overview_ch4, nrow = 1)

myfilename <- paste("000_overview_CO2_CH4_season",paste(unique(table_ch4$sampling), collapse = "_"), sep = "_")
ggsave(plot = p_overview_co2_ch4, filename = paste0(myfilename,".jpg"), path = plots_path,
       width = 10, height = 4, dpi = 300, units = 'in', scale = 1.1)




p_overview_co2 <- ggplot(table_co2, aes(sampling, best.flux))+
  geom_hline(yintercept = 0)+
  geom_jitter(alpha=0.5, aes(colour = lightCondition))+
  geom_boxplot(alpha=0.5, width=0.3)+
  theme_article()+
  ylab("CO2 flux umol/m2/s")+
  scale_fill_viridis_d(begin = 0.2, end = 0.9)+
  scale_colour_viridis_d(begin = 0.2, end = 0.9, option = "C")+
  facet_grid(strata~., scales = "free_y")+
  theme(legend.position = 'none')+xlab("")

p_overview_ch4 <- ggplot(table_ch4, aes(sampling, best.flux))+
  geom_hline(yintercept = 0)+
  geom_jitter(alpha=0.5, aes(colour = lightCondition))+
  geom_boxplot(alpha=0.5, width=0.3)+
  theme_article()+
  ylab("CH4 flux nmol/m2/s")+
  scale_fill_viridis_d(begin = 0.2, end = 0.9)+
  scale_colour_viridis_d(begin = 0.2, end = 0.9, option = "C")+
  facet_grid(strata~., scales = "free_y")+
  theme(legend.position = 'right')+xlab("")
p_overview_co2_ch4 <- ggarrange(p_overview_co2, p_overview_ch4, nrow = 1)

myfilename <- paste("000_overview_CO2_CH4_strata_season",paste(unique(table_ch4$sampling), collapse = "_"), sep = "_")
ggsave(plot = p_overview_co2_ch4, filename = paste0(myfilename,".jpg"), path = plots_path,
       width = 10, height = 4, dpi = 300, units = 'in', scale = 1.1)





p_overview_co2 <- ggplot(table_co2, aes(pilotsite , best.flux, colour = lightCondition))+
  geom_hline(yintercept = 0)+
  # geom_boxplot(alpha=0.99)+
  geom_jitter(alpha=0.5, width = 0.2)+
  theme_article()+
  ylab("CO2 flux umol/m2/s")+
  scale_fill_viridis_d(begin = 0.2, end = 0.9)+
  scale_colour_viridis_d(begin = 0.2, end = 0.9, option = "C")+
  # facet_grid(. ~strata)+
  theme(legend.position = 'none')+xlab("")


p_overview_ch4 <- ggplot(table_ch4, aes(pilotsite , best.flux, colour = lightCondition))+
  geom_hline(yintercept = 0)+
  # geom_boxplot(alpha=0.99)+
  geom_jitter(alpha=0.5, width = 0.2)+
  theme_article()+
  ylab("CH4 flux nmol/m2/s")+
  scale_fill_viridis_d(begin = 0.2, end = 0.9)+
  scale_colour_viridis_d(begin = 0.2, end = 0.9, option = "C")+
  # facet_grid(. ~strata)+
  theme(legend.position = 'right')+xlab("")

p_overview_co2_ch4 <- ggarrange(p_overview_co2, p_overview_ch4, nrow = 1)

myfilename <- paste("000_overview_CO2_CH4_sites",paste(unique(table_ch4$sampling), collapse = "_"), sep = "_")
ggsave(plot = p_overview_co2_ch4, filename = paste0(myfilename,".jpg"), path = plots_path, 
       width = 10, height = 4, dpi = 300, units = 'in', scale = 1.1)





p_overview_co2 <- ggplot(table_co2, aes(subsite , best.flux, colour = lightCondition))+
  geom_hline(yintercept = 0)+
  geom_boxplot(alpha=0.99)+
  # geom_jitter(alpha=0.5, width = 0.2)+
  theme_article()+
  ylab("CO2 flux umol/m2/s")+
  scale_fill_viridis_d(begin = 0.2, end = 0.9)+
  scale_colour_viridis_d(begin = 0.2, end = 0.9, option = "C")+
  facet_grid(lightCondition~pilotsite)+
  theme(legend.position = 'top')+xlab("")+
  # ggtitle(paste0("CO2 fluxes, ",length(unique(table_co2$UniqueID)), " valid incubations"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))




p_overview_ch4 <- ggplot(table_ch4, aes(subsite , best.flux))+
  geom_hline(yintercept = 0)+
  geom_boxplot(alpha=0.99)+
  # geom_jitter(alpha=0.5, width = 0.2)+
  theme_article()+
  ylab("CH4 flux nmol/m2/s")+
  scale_fill_viridis_d(begin = 0.2, end = 0.9)+
  scale_colour_viridis_d(begin = 0.2, end = 0.9, option = "C")+
  facet_wrap(. ~pilotsite, scales = "free_y", nrow = 1)+
  theme(legend.position = 'none')+xlab("")+
  # ggtitle(paste0("CH4 fluxes, ",length(unique(table_ch4$UniqueID)), " valid incubations"))
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


p_overview_co2_ch4 <- ggarrange(p_overview_co2, p_overview_ch4)

myfilename <- paste("000_overview_CO2_CH4_subsites",paste(unique(table_ch4$sampling), collapse = "_"), sep = "_")
ggsave(plot = p_overview_co2_ch4, filename = paste0(myfilename,".jpg"), path = plots_path, 
       width = 8, height = 5, dpi = 300, units = 'in', scale = 1.1)







p_overview_co2 <- ggplot(table_co2, aes(sampling , best.flux, colour = lightCondition))+
  geom_hline(yintercept = 0)+
  # geom_boxplot(alpha=0.99)+
  geom_jitter(alpha=0.5, width = 0.2)+
  theme_article()+
  ylab("CO2 flux umol/m2/s")+
  scale_fill_viridis_d(begin = 0.2, end = 0.9)+
  scale_colour_viridis_d(begin = 0.2, end = 0.9, option = "C")+
  # facet_grid(. ~strata)+
  theme(legend.position = 'right')+xlab("")+
  ggtitle(paste0("CO2 fluxes, ",length(unique(table_co2$UniqueID)), " valid incubations"))


p_overview_ch4 <- ggplot(table_ch4, aes(sampling , best.flux, colour = lightCondition))+
  geom_hline(yintercept = 0)+
  # geom_boxplot(alpha=0.99)+
  geom_jitter(alpha=0.5, width = 0.2)+
  theme_article()+
  ylab("CH4 flux nmol/m2/s")+
  scale_fill_viridis_d(begin = 0.2, end = 0.9)+
  scale_colour_viridis_d(begin = 0.2, end = 0.9, option = "C")+
  # facet_grid(. ~strata)+
  theme(legend.position = 'none')+xlab("")+
  ggtitle(paste0("CH4 fluxes, ",length(unique(table_ch4$UniqueID)), " valid incubations"))


p_overview_co2_ch4 <- ggarrange(p_overview_co2, p_overview_ch4)

myfilename <- paste("000_overview_CO2_CH4_season",paste(unique(table_ch4$sampling), collapse = "_"), sep = "_")
ggsave(plot = p_overview_co2_ch4, filename = paste0(myfilename,".jpg"), path = plots_path, 
       width = 4, height = 5, dpi = 300, units = 'in', scale = 1.1)








plot_overview <- function(mytable, label, title, variable){
  
  plt <- ggplot(mytable, aes(sampling, best.flux,
                             fill = lightCondition))+
    geom_hline(yintercept = 0)+
    geom_boxplot(alpha=0.5)+
    # geom_jitter(width = 0.2, aes(colour = strata), size=2)+
    theme_article()+
    xlab("subsite")+
    ylab(label)+
    ggtitle(title)+
    scale_fill_viridis_d(begin = 0.2, end = 0.9)+
    scale_colour_viridis_d(begin = 0.2, end = 0.9, option = "C")+
    facet_grid(pilotsite~subsite, scale="free_y")+theme(legend.position = "top")
  
  setwd(results_path)
  myfilename <- paste(variable,paste(unique(mytable$sampling), collapse = "_"), sep = "_")
  ggsave(plot = plt, filename = paste0(myfilename,".jpg"), path = plots_path, 
         width = 6, height = 6, dpi = 300, units = 'in', scale = 1.1)
  
  
  for (ps in unique(mytable$pilotsite)){
    plt <- ggplot(mytable[mytable$pilotsite==ps,], 
                  aes(sampling, best.flux, colour = strata, fill=strata))+
      geom_hline(yintercept = 0)+
      geom_jitter(width = 0.2, size=2)+
      geom_boxplot(alpha=0.5)+
      theme_article()+
      xlab("subsite")+
      ylab(label)+
      ggtitle(paste0(title, " in ",ps))+
      scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "C")+
      scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "C")+
      facet_grid(lightCondition~subsite, scales = "free_y")+theme(legend.position = "top")
    
    setwd(plots_path)
    myfilename <- paste(variable,ps,paste(unique(mytable$sampling), collapse = "_"), sep = "_")
    ggsave(plot = plt, filename = paste0(myfilename,".jpg"), path = plots_path, 
           width = 8, height = 5, dpi = 300, units = 'in', scale = 1.1)
  }
}


plot_overview(mytable = table_co2, label = "CO2 flux umol/m2/s", title = "CO2 flux", variable = "CO2")
plot_overview(mytable = table_ch4, label = "CH4 flux nmol/m2/s", title = "CH4 flux", variable = "CH4")



plot_overview(mytable = table_co2[which(table_co2$strata=="open water" & table_co2$lightCondition=='dark'),],
              label = "CO2 flux umol/m2/s", title = "CO2 flux", variable = "CO2_open_water_dark")
plot_overview(mytable = table_ch4[which(table_ch4$strata=="open water" & table_ch4$lightCondition=='dark'),],
              label = "CH4 flux nmol/m2/s", title = "CH4 flux", variable = "CH4_open_water_dark")



plot_overview(mytable = table_co2[which(table_co2$strata=="vegetated"),],
              label = "CO2 flux umol/m2/s", title = "CO2 flux", variable = "CO2_vegetation")
plot_overview(mytable = table_ch4[which(table_ch4$strata=="vegetated"),],
              label = "CH4 flux nmol/m2/s", title = "CH4 flux", variable = "CH4_vegetation")






# ---- tentative NET ECOSYSTEM EXCHANGE ----

isFs <- T
for (s in unique(paste0(table_co2$sampling, table_co2$siteID))){
  
  sel <- table_co2[paste0(table_co2$sampling, table_co2$siteID) == s,]
  
  sel_light <- sel[sel$lightCondition == "transparent",]
  sel_dark <- sel[sel$lightCondition == "dark",]
  
  if(dim(sel_light)[1]>=3){
    co2_light = mean(sel_light$best.flux, na.rm = T) 
  } else {
    co2_light <- 0
  }
  
  if(dim(sel_dark)[1]>=3){
    co2_dark = mean(sel_dark$best.flux, na.rm = T) 
  } else {
    co2_dark <- 0
  }
  
  
  df_NEE.tmp <- data.frame(pilotsite = unique(sel$pilotsite),
                           subsite = unique(sel$subsite),
                           siteID = unique(sel$siteID),
                           sampling = unique (sel$sampling),
                           co2_light = co2_light, # umol/m2/s
                           co2_dark = co2_dark, # umol/m2/s
                           ch4 = 25 * mean(table_ch4$best.flux[paste0(table_ch4$sampling, table_ch4$siteID) == s], na.rm = T)/1000) # umol/m2/s
  
  df_NEE.tmp$NEE <- (12*df_NEE.tmp$co2_light+12*df_NEE.tmp$co2_dark+24*df_NEE.tmp$ch4)*3600*1e-6 * 12 #  # umol/m2/s to gC m-2 day-1
  
  if(isFs){
    isFs <- F
    df_NEE <- df_NEE.tmp
  } else {
    df_NEE <- rbind(df_NEE, df_NEE.tmp)
  }
}

plt_NEE <- ggplot(df_NEE, aes(subsite, NEE))+
  geom_hline(yintercept = 0)+
  geom_boxplot(alpha=0.99)+
  # geom_jitter(alpha=0.5, width = 0.2)+
  theme_article()+
  ylab("Global Warming Potential gC m-2 day-1")+
  scale_fill_viridis_d(begin = 0.2, end = 0.9)+
  scale_colour_viridis_d(begin = 0.2, end = 0.9, option = "C")+
  facet_wrap(. ~pilotsite, scales = "free_y", nrow = 1)+
  theme(legend.position = 'right')+xlab("")

myfilename <- paste("TestNet_Ecosystem_Exchange",paste(unique(table_ch4$sampling), collapse = "_"), sep = "_")
ggsave(plot = plt_NEE, filename = paste0(myfilename,".jpg"), path = plots_path, 
       width = 10, height = 3.5, dpi = 300, units = 'in', scale = 1.1)



# ---- Effect of bubbling (not adapted)----

table_ch4$p_ebullition <- table_ch4$ebullition/table_ch4$total_estimated *100

ggplot(table_ch4[table_ch4$p_ebullition >=0 & table_ch4$p_ebullition <=1,], 
       aes(LM.p.val, ebullition/total_estimated))+
  geom_point()+theme_article()

ggplot(table_ch4[which(table_ch4$p_ebullition >=0 & table_ch4$p_ebullition <=100),], 
       aes(subsite, ebullition))+
  # geom_hline(yintercept = 100)+
  geom_jitter(alpha=0.5, width = 0.2)+
  geom_boxplot(alpha=0.8, width=0.5)+
  theme_article()+
  scale_y_log10()+
  facet_wrap(pilotsite~.)+
  ylab("CH4 ebullition flux [nmol/m2/s]")+
  xlab("")



p_ebull <- ggplot(table_ch4[which(table_ch4$p_ebullition >=0 & table_ch4$p_ebullition <=100),], 
                  aes(pilotsite, ebullition))+
  # geom_hline(yintercept = 100)+
  geom_jitter(alpha=0.5, width = 0.2)+
  geom_boxplot(alpha=0.8, width=0.5)+
  theme_article()+
  scale_y_log10()+
  ylab("CH4 ebullition flux [nmol/m2/s]")+
  xlab("")


p_perc_ebull <- ggplot(table_ch4[which(table_ch4$p_ebullition >=0 & table_ch4$p_ebullition <=100),], 
                       aes(pilotsite, p_ebullition))+
  # geom_hline(yintercept = 100)+
  geom_jitter(alpha=0.5, width = 0.2)+
  geom_boxplot(alpha=0.8, width=0.5)+
  theme_article()+
  ylab("Contribution of bubbling to total CH4 flux [%]")+
  xlab("")

p_ebull <- ggarrange(p_ebull, p_perc_ebull)
myfilename <- paste("000_overview_bubbling",paste(unique(table_ch4$sampling), collapse = "_"), sep = "_")
ggsave(plot = p_ebull, filename = paste0(myfilename,".jpg"), path = plots_path, 
       width = 8, height = 5, dpi = 300, units = 'in', scale = 1.1)


ggplot(table_ch4[which(table_ch4$p_ebullition >=0 & table_ch4$p_ebullition <=100),], 
       aes(subsite, ebullition, fill=subsite))+
  # geom_hline(yintercept = 100)+
  geom_jitter(alpha=0.5)+
  # geom_boxplot(alpha=0.8, width=0.3)+
  theme_article()+
  # ylab("CH4 flux nmol/m2/s")+
  scale_fill_viridis_d(begin = 0.2, end = 0.9)+
  scale_colour_viridis_d(begin = 0.2, end = 0.9, option = "C")+
  # ylim(c(0,100))+
  # scale_y_log10()+
  facet_grid(pilotsite~., scales = "free_y")+
  # theme(legend.position = 'none')+
  xlab("")

table_ch4_ebull <- table_ch4[!is.na(table_ch4$p_ebullition),]
n_50 <- length(which(table_ch4_ebull$p_ebullition > 50))

p_prop <- ggplot(table_ch4_ebull[order(table_ch4_ebull$p_ebullition),], 
                 aes(seq_along(table_ch4_ebull$p_ebullition)/dim(table_ch4_ebull)[1]*100,p_ebullition))+
  geom_hline(yintercept = 100)+
  geom_hline(yintercept = 50, alpha=0.2)+
  # geom_jitter(alpha=0.5)+
  geom_point(alpha=0.5, width=0.3, aes(colour = table_ch4_ebull$p_ebullition[order(table_ch4_ebull$p_ebullition)] > 50))+
  theme_article()+
  # ylab("CH4 flux nmol/m2/s")+
  scale_fill_viridis_d(begin = 0.2, end = 0.9)+
  scale_colour_viridis_d("over 50%",begin = 0.2, end = 0.9, option = "A", direction = -1)+
  ylim(c(0,120))+
  # scale_y_log10()+
  # facet_wrap(strata~.)+
  theme(legend.position = 'none')+
  xlab("Proportion of timeseries [%]")+
  ylab("Contribution of bubbling to total CH4 flux [%]")+
  ggtitle(paste0("Bubbling > 50% of F_total for ",round(n_50/dim(table_ch4_ebull)[1]*100),"% of the measurements"))

p_ebull <- ggplot(table_ch4_ebull, 
                  aes(pilotsite, ebullition))+
  # geom_hline(yintercept = 100)+
  geom_jitter(alpha=0.5, width = 0.2, aes(colour = table_ch4_ebull$p_ebullition > 50))+
  geom_boxplot(alpha=0.8, width=0.5)+
  scale_colour_viridis_d("over 50%",begin = 0.2, end = 0.9, option = "A", direction = -1)+
  theme_article()+
  scale_y_log10()+
  ylab("CH4 ebullition flux [nmol/m2/s]")+
  theme(legend.position = 'bottom')+
  xlab("")



ggarrange(p_ebull, p_prop, nrow = 1)

