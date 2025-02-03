
# ---
# Authors: Miguel Cabrera
# Project: "RESTORE4Cs"
# date: "Feb 2025"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
# This script analyses the quality of fit of all fluxes. First by ranking and flagging incubations based on regression fit. Then by flagging incubations that include artefacts (examining deltaGHG and flagging incubations with more than 5% of deltaGHG==0)




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
results_path <- paste0(dropbox_root,"/GHG/Processed data/computed_flux/")

plots_path <- paste0(dropbox_root,"/GHG/Processed data/computed_flux/plots")

fieldsheetpath <- paste0(dropbox_root,"/GHG/Fieldsheets")

#Load csv files produced by raw2flux.R. script
setwd(results_path)

listf <- list.files(path = results_path, pattern = "^S[0-9]{1}_fluxes", all.files = T, full.names = T, recursive = F)
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

# table_results_all <- table_results_all[table_results_all$CO2_best.flux<1000,]

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


table_co2 <- get_full_detailed_table(table_results_all, variable = "co2")
table_ch4 <- get_full_detailed_table(table_results_all, "ch4")



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


#Only a few co2 fluxes without any flag:
table_co2[which(is.na(table_co2$quality.check)),]

#Goflux meaning of different flags: 
#SE: "Standard Error" it tells us that the noise in our incubation is larger than the instrument accuracy (i.e. we have more noise than we would expect just from "instrumental noise").
#MDF: "minimum detectable flux". I.e. flux is Below detection limit. 
#gfact > 2: the ratio bettewwn the non-linear(HM) flux  estimate and the linear flux estimate is larger than 2: i.e. HM.flux > 2*LM.flux . IF this happens, we assume that the HM.flux is overestimated, so we take the LM.flux instead (more or less arbitrary decission, different thresholds can be selected (relaxed gfact>4, medium gfact>2, strict gfact>1.25))
#gfact < 0.5: same as above but for the case when LM.flux is higher than HM.flux
#HM.flux is NA: non-linear model failed to return a flux.



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


#No ch4 fluxes without any flag:
table_ch4[which(is.na(table_ch4$quality.check)),]


####poraqui####

#Examine R2 for CO2

#LM.flux R2

ggplot(table_co2)+
  geom_histogram(aes(x=LM.r2, fill="linear"),bins = 40)+
  geom_histogram(aes(x=HM.r2, fill="non-linear"),bins = 40)+
  labs(fill="R2")

#Examine R2 for CH4
ggplot(table_ch4)+
  geom_histogram(aes(x=LM.r2, fill="linear"),bins = 40)+
  geom_histogram(aes(x=HM.r2, fill="non-linear"),bins = 40)+
  labs(fill="R2")


#extract a few candidates for inspection based on R2:
#CO2 good fits r2>0.95: all looks good
table_co2 %>% 
  filter(between(x=HM.r2,
                 lower= quantile(.$HM.r2,probs = 0.5, na.rm = T),
                 upper =quantile(.$HM.r2,probs = 0.57, na.rm = T))) %>% 
  filter(model=="HM") %>% 
  select(UniqueID, best.flux, HM.r2)%>% 
  filter(best.flux%in%c(max(best.flux),min(best.flux),quantile(best.flux, 0.5)))

#Some more noise, dubious minor outliers
table_co2 %>% 
  filter(between(x=HM.r2,
                 lower= quantile(.$HM.r2,probs = 0.4, na.rm = T),
                 upper =quantile(.$HM.r2,probs = 0.43, na.rm = T))) %>% 
  filter(model=="HM") %>% 
  select(UniqueID, best.flux, HM.r2)%>% 
  filter(best.flux%in%c(max(best.flux),min(best.flux),quantile(best.flux, 0.5)))

#Poraqui
table_co2 %>% 
  filter(between(x=HM.r2,
                 lower= quantile(.$HM.r2,probs = 0.35, na.rm = T),
                 upper =quantile(.$HM.r2,probs = 0.36, na.rm = T))) %>% 
  filter(model=="HM") %>% 
  select(UniqueID, best.flux, HM.r2)%>% 
  filter(best.flux%in%c(max(best.flux),min(best.flux),quantile(best.flux, 0.5)))

table_co2 %>% 
  filter(between(x=HM.r2,
                 lower= quantile(.$HM.r2,probs = 0.30, na.rm = T),
                 upper =quantile(.$HM.r2,probs = 0.32, na.rm = T))) %>% 
  filter(model=="HM") %>% 
  select(UniqueID, best.flux, HM.r2)%>% 
  filter(best.flux%in%c(max(best.flux),min(best.flux),quantile(best.flux, 0.5)))

table_co2 %>% 
  filter(between(x=HM.r2,
                 lower= quantile(.$HM.r2,probs = 0.25, na.rm = T),
                 upper =quantile(.$HM.r2,probs = 0.26, na.rm = T))) %>% 
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




# ---- plot CO2 and CH4 flux across sites and seasons ----





p_overview_co2 <- ggplot(table_co2_sel, aes(lightCondition, best.flux, colour = lightCondition))+
  geom_hline(yintercept = 0)+
  geom_jitter(alpha=0.5)+
  # geom_boxplot(alpha=0.8)+
  theme_article()+
  ylab("CO2 flux umol/m2/s")+
  scale_fill_viridis_d(begin = 0.2, end = 0.9)+
  scale_colour_viridis_d(begin = 0.2, end = 0.9, option = "C")+
  # facet_wrap(strata~.)+
  theme(legend.position = 'none')+xlab("")+
  ggtitle(paste0("CO2 fluxes, ",length(unique(table_co2_sel$UniqueID)), " valid incubations"))


p_overview_ch4 <- ggplot(table_ch4_sel, aes(lightCondition, best.flux, colour = lightCondition))+
  geom_hline(yintercept = 0)+
  geom_jitter(alpha=0.5)+
  # geom_boxplot(alpha=0.8)+
  theme_article()+
  ylab("CH4 flux nmol/m2/s")+
  scale_fill_viridis_d(begin = 0.2, end = 0.9)+
  scale_colour_viridis_d(begin = 0.2, end = 0.9, option = "C")+
  # facet_wrap(strata~.)+
  theme(legend.position = 'right')+xlab("")+
  ggtitle(paste0("CH4 fluxes, ",length(unique(table_ch4_sel$UniqueID)), " valid incubations"))

p_overview_co2_ch4 <- ggarrange(p_overview_co2, p_overview_ch4, nrow = 1)

myfilename <- paste("000_overview_CO2_CH4_lightCond",paste(unique(table_ch4_sel$sampling), collapse = "_"), sep = "_")
ggsave(plot = p_overview_co2_ch4, filename = paste0(myfilename,".jpg"), path = plots_path,
       width = 10, height = 4, dpi = 300, units = 'in', scale = 1.1)






p_overview_co2 <- ggplot(table_co2_sel, aes(strata, best.flux, colour = lightCondition))+
  geom_hline(yintercept = 0)+
  geom_jitter(alpha=0.5)+
  # geom_boxplot(alpha=0.8)+
  theme_article()+
  ylab("CO2 flux umol/m2/s")+
  scale_fill_viridis_d(begin = 0.2, end = 0.9)+
  scale_colour_viridis_d(begin = 0.2, end = 0.9, option = "C")+
  # facet_wrap(sampling~.)+
  theme(legend.position = 'none')+xlab("")

p_overview_ch4 <- ggplot(table_ch4_sel, aes(strata, best.flux, colour = lightCondition))+
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

myfilename <- paste("000_overview_CO2_CH4_strata",paste(unique(table_ch4_sel$sampling), collapse = "_"), sep = "_")
ggsave(plot = p_overview_co2_ch4, filename = paste0(myfilename,".jpg"), path = plots_path,
       width = 10, height = 4, dpi = 300, units = 'in', scale = 1.1)



p_overview_co2 <- ggplot(table_co2_sel, aes(sampling, best.flux))+
  geom_hline(yintercept = 0)+
  geom_jitter(alpha=0.4, aes(colour = lightCondition))+
  geom_boxplot(alpha=0.5, width=0.3)+
  theme_article()+
  ylab("CO2 flux umol/m2/s")+
  scale_fill_viridis_d(begin = 0.2, end = 0.9)+
  scale_colour_viridis_d(begin = 0.2, end = 0.9, option = "C")+
  # facet_wrap(sampling~.)+
  theme(legend.position = 'none')+xlab("")

p_overview_ch4 <- ggplot(table_ch4_sel, aes(sampling, best.flux))+
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

myfilename <- paste("000_overview_CO2_CH4_season",paste(unique(table_ch4_sel$sampling), collapse = "_"), sep = "_")
ggsave(plot = p_overview_co2_ch4, filename = paste0(myfilename,".jpg"), path = plots_path,
       width = 10, height = 4, dpi = 300, units = 'in', scale = 1.1)




p_overview_co2 <- ggplot(table_co2_sel, aes(sampling, best.flux))+
  geom_hline(yintercept = 0)+
  geom_jitter(alpha=0.5, aes(colour = lightCondition))+
  geom_boxplot(alpha=0.5, width=0.3)+
  theme_article()+
  ylab("CO2 flux umol/m2/s")+
  scale_fill_viridis_d(begin = 0.2, end = 0.9)+
  scale_colour_viridis_d(begin = 0.2, end = 0.9, option = "C")+
  facet_grid(strata~., scales = "free_y")+
  theme(legend.position = 'none')+xlab("")

p_overview_ch4 <- ggplot(table_ch4_sel, aes(sampling, best.flux))+
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

myfilename <- paste("000_overview_CO2_CH4_strata_season",paste(unique(table_ch4_sel$sampling), collapse = "_"), sep = "_")
ggsave(plot = p_overview_co2_ch4, filename = paste0(myfilename,".jpg"), path = plots_path,
       width = 10, height = 4, dpi = 300, units = 'in', scale = 1.1)





p_overview_co2 <- ggplot(table_co2_sel, aes(pilotsite , best.flux, colour = lightCondition))+
  geom_hline(yintercept = 0)+
  # geom_boxplot(alpha=0.99)+
  geom_jitter(alpha=0.5, width = 0.2)+
  theme_article()+
  ylab("CO2 flux umol/m2/s")+
  scale_fill_viridis_d(begin = 0.2, end = 0.9)+
  scale_colour_viridis_d(begin = 0.2, end = 0.9, option = "C")+
  # facet_grid(. ~strata)+
  theme(legend.position = 'none')+xlab("")

p_overview_ch4 <- ggplot(table_ch4_sel, aes(pilotsite , best.flux, colour = lightCondition))+
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

myfilename <- paste("000_overview_CO2_CH4_sites",paste(unique(table_ch4_sel$sampling), collapse = "_"), sep = "_")
ggsave(plot = p_overview_co2_ch4, filename = paste0(myfilename,".jpg"), path = plots_path, 
       width = 10, height = 4, dpi = 300, units = 'in', scale = 1.1)





p_overview_co2 <- ggplot(table_co2_sel, aes(subsite , best.flux, colour = lightCondition))+
  geom_hline(yintercept = 0)+
  geom_boxplot(alpha=0.99)+
  # geom_jitter(alpha=0.5, width = 0.2)+
  theme_article()+
  ylab("CO2 flux umol/m2/s")+
  scale_fill_viridis_d(begin = 0.2, end = 0.9)+
  scale_colour_viridis_d(begin = 0.2, end = 0.9, option = "C")+
  facet_grid(lightCondition~pilotsite)+
  theme(legend.position = 'top')+xlab("")+
  # ggtitle(paste0("CO2 fluxes, ",length(unique(table_co2_sel$UniqueID)), " valid incubations"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))




p_overview_ch4 <- ggplot(table_ch4_sel, aes(subsite , best.flux))+
  geom_hline(yintercept = 0)+
  geom_boxplot(alpha=0.99)+
  # geom_jitter(alpha=0.5, width = 0.2)+
  theme_article()+
  ylab("CH4 flux nmol/m2/s")+
  scale_fill_viridis_d(begin = 0.2, end = 0.9)+
  scale_colour_viridis_d(begin = 0.2, end = 0.9, option = "C")+
  facet_wrap(. ~pilotsite, scales = "free_y", nrow = 1)+
  theme(legend.position = 'none')+xlab("")+
  # ggtitle(paste0("CH4 fluxes, ",length(unique(table_ch4_sel$UniqueID)), " valid incubations"))
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


p_overview_co2_ch4 <- ggarrange(p_overview_co2, p_overview_ch4)

myfilename <- paste("000_overview_CO2_CH4_subsites",paste(unique(table_ch4_sel$sampling), collapse = "_"), sep = "_")
ggsave(plot = p_overview_co2_ch4, filename = paste0(myfilename,".jpg"), path = plots_path, 
       width = 8, height = 5, dpi = 300, units = 'in', scale = 1.1)







p_overview_co2 <- ggplot(table_co2_sel, aes(sampling , best.flux, colour = lightCondition))+
  geom_hline(yintercept = 0)+
  # geom_boxplot(alpha=0.99)+
  geom_jitter(alpha=0.5, width = 0.2)+
  theme_article()+
  ylab("CO2 flux umol/m2/s")+
  scale_fill_viridis_d(begin = 0.2, end = 0.9)+
  scale_colour_viridis_d(begin = 0.2, end = 0.9, option = "C")+
  # facet_grid(. ~strata)+
  theme(legend.position = 'right')+xlab("")+
  ggtitle(paste0("CO2 fluxes, ",length(unique(table_co2_sel$UniqueID)), " valid incubations"))


p_overview_ch4 <- ggplot(table_ch4_sel, aes(sampling , best.flux, colour = lightCondition))+
  geom_hline(yintercept = 0)+
  # geom_boxplot(alpha=0.99)+
  geom_jitter(alpha=0.5, width = 0.2)+
  theme_article()+
  ylab("CH4 flux nmol/m2/s")+
  scale_fill_viridis_d(begin = 0.2, end = 0.9)+
  scale_colour_viridis_d(begin = 0.2, end = 0.9, option = "C")+
  # facet_grid(. ~strata)+
  theme(legend.position = 'none')+xlab("")+
  ggtitle(paste0("CH4 fluxes, ",length(unique(table_ch4_sel$UniqueID)), " valid incubations"))


p_overview_co2_ch4 <- ggarrange(p_overview_co2, p_overview_ch4)

myfilename <- paste("000_overview_CO2_CH4_season",paste(unique(table_ch4_sel$sampling), collapse = "_"), sep = "_")
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


plot_overview(mytable = table_co2_sel, label = "CO2 flux umol/m2/s", title = "CO2 flux", variable = "CO2")
plot_overview(mytable = table_ch4_sel, label = "CH4 flux nmol/m2/s", title = "CH4 flux", variable = "CH4")



plot_overview(mytable = table_co2_sel[which(table_co2_sel$strata=="open water" & table_co2_sel$lightCondition=='dark'),],
              label = "CO2 flux umol/m2/s", title = "CO2 flux", variable = "CO2_open_water_dark")
plot_overview(mytable = table_ch4_sel[which(table_ch4_sel$strata=="open water" & table_ch4_sel$lightCondition=='dark'),],
              label = "CH4 flux nmol/m2/s", title = "CH4 flux", variable = "CH4_open_water_dark")



plot_overview(mytable = table_co2_sel[which(table_co2_sel$strata=="vegetated"),],
              label = "CO2 flux umol/m2/s", title = "CO2 flux", variable = "CO2_vegetation")
plot_overview(mytable = table_ch4_sel[which(table_ch4_sel$strata=="vegetated"),],
              label = "CH4 flux nmol/m2/s", title = "CH4 flux", variable = "CH4_vegetation")






# ---- tentative NET ECOSYSTEM EXCHANGE ----

isFs <- T
for (s in unique(paste0(table_co2_sel$sampling, table_co2_sel$siteID))){
  
  sel <- table_co2_sel[paste0(table_co2_sel$sampling, table_co2_sel$siteID) == s,]
  
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

myfilename <- paste("Net_Ecosystem_Exchange",paste(unique(table_ch4_sel$sampling), collapse = "_"), sep = "_")
ggsave(plot = plt_NEE, filename = paste0(myfilename,".jpg"), path = plots_path, 
       width = 10, height = 3.5, dpi = 300, units = 'in', scale = 1.1)



# ---- Effect of bubbling ----

table_ch4_sel$p_ebullition <- table_ch4_sel$ebullition/table_ch4_sel$total_estimated *100

ggplot(table_ch4_sel[table_ch4_sel$p_ebullition >=0 & table_ch4_sel$p_ebullition <=1,], 
       aes(LM.p.val, ebullition/total_estimated))+
  geom_point()+theme_article()

ggplot(table_ch4_sel[which(table_ch4_sel$p_ebullition >=0 & table_ch4_sel$p_ebullition <=100),], 
       aes(subsite, ebullition))+
  # geom_hline(yintercept = 100)+
  geom_jitter(alpha=0.5, width = 0.2)+
  geom_boxplot(alpha=0.8, width=0.5)+
  theme_article()+
  scale_y_log10()+
  facet_wrap(pilotsite~.)+
  ylab("CH4 ebullition flux [nmol/m2/s]")+
  xlab("")



p_ebull <- ggplot(table_ch4_sel[which(table_ch4_sel$p_ebullition >=0 & table_ch4_sel$p_ebullition <=100),], 
                  aes(pilotsite, ebullition))+
  # geom_hline(yintercept = 100)+
  geom_jitter(alpha=0.5, width = 0.2)+
  geom_boxplot(alpha=0.8, width=0.5)+
  theme_article()+
  scale_y_log10()+
  ylab("CH4 ebullition flux [nmol/m2/s]")+
  xlab("")


p_perc_ebull <- ggplot(table_ch4_sel[which(table_ch4_sel$p_ebullition >=0 & table_ch4_sel$p_ebullition <=100),], 
                       aes(pilotsite, p_ebullition))+
  # geom_hline(yintercept = 100)+
  geom_jitter(alpha=0.5, width = 0.2)+
  geom_boxplot(alpha=0.8, width=0.5)+
  theme_article()+
  ylab("Contribution of bubbling to total CH4 flux [%]")+
  xlab("")

p_ebull <- ggarrange(p_ebull, p_perc_ebull)
myfilename <- paste("000_overview_bubbling",paste(unique(table_ch4_sel$sampling), collapse = "_"), sep = "_")
ggsave(plot = p_ebull, filename = paste0(myfilename,".jpg"), path = plots_path, 
       width = 8, height = 5, dpi = 300, units = 'in', scale = 1.1)


ggplot(table_ch4_sel[which(table_ch4_sel$p_ebullition >=0 & table_ch4_sel$p_ebullition <=100),], 
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

table_ch4_ebull <- table_ch4_sel[!is.na(table_ch4_sel$p_ebullition),]
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

# repo_root <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
# files.sources = list.files(path = paste0(repo_root,"/functions"), full.names = T)
# for (f in files.sources){source(f)}


# ---- Directories ----
dropbox_root <- "C:/Users/Miguel/Dropbox/RESTORE4Cs - Fieldwork/Data" # You have to make sure this is pointing to the write folder on your local machine
# datapath <- paste0(dropbox_root,"/GHG/RAW data")
# fieldsheetpath <- paste0(dropbox_root,"/GHG/Fieldsheets")
# loggerspath <- paste0(datapath,"/RAW Data Logger")
# RData_path <- paste0(dropbox_root,"/GHG/Processed data/RData/")
results_path <- paste0(dropbox_root,"/GHG/Processed data/computed_flux/")
summary_path <- paste0(dropbox_root,"/GHG/Processed data/computed_flux/summaries/")
plots_path <- paste0(dropbox_root,"/GHG/Processed data/computed_flux/summaries/sumplots/")


#Read all computed fluxes per season
season_csvs<- list.files(path = results_path, pattern = "^S[0-9]{1}_fluxes")
dat<- data.frame()
for (i in season_csvs){
  s<- read.csv(paste0(results_path,i))
  
  dat<- rbind(dat,s)
}
rm(s, i)


####OPEN wATER DATASET####
ow<- dat %>% 
  separate(subsite, into = c("season", "site", "subsite_only"), remove = F) %>% 
  filter(strata=="open water")

ow_wrongchamber<-ow %>% filter(chamberType!="floating") %>% pull(UniqueID)
#S1-VA-P2-16 open water with tube, transparent and dark, test for biofilm
#S1-VA-R2-(13,14,15) only 2cm of water, used tube chamber instead of floating

ow_wronglight<- ow %>% filter(lightCondition!="dark") %>% pull(UniqueID)
#S1-DA-R2-11 wrong lightcondition, fixed in fieldsheet but not recalculated as of 20241203
#S1-DU-A1-(2,3,4) 3 open water with floating chamber but lightcondition written "transparent"
#S1-VA-P2-16 open water with tube, transparent and dark, test for biofilm
#S4-CA-R1-8 wrong lightcondition, fixed in fieldsheet but not recalculated as of 20241203



#Filter OW dataset to remove anything suspicious:

unique(ow$CO2_quality.check)
  #check decissions for CO2_quality.check and CH4_quality.check:
    # SE indicates noisy timeseries, 
    # g-fact. indicates ratio between nonlinear model(HM) and linear one(LM)  
    # MDF: below minimum detectable flow for instrument accuracy (~below detection limit)
    #"HM.flux is NA" 


ow_tointegrate<- ow %>% 
  filter(!UniqueID%in%c(ow_wrongchamber,ow_wronglight))
  # filter(!grepl("g-fact. > 2", CO2_quality.check))
  
unique(ow$CO2_quality.check)

ow_long<-ow %>% 
  filter(chamberType=="floating"&lightCondition=="dark") %>% 
  select(UniqueID, season,site,subsite_only, water_depth, CO2_best.flux,CO2_best.model,CO2_quality.check, )




#CH4 from Vegetation dark vs light
dat %>% 
  filter(strata=="vegetated") %>%
  filter(!UniqueID%in%c("s1-va-a1-12-v-t-14:58","s1-va-a1-12-v-d-15:06")) %>% #Not well done incubations, there are repetittions
  filter(!grepl("s1-da-r1-8|s1-da-r2-12",UniqueID)) %>% #wrong lightcondition, fixed in fieldsheets but not recalculated yet
  mutate(plotID=str_extract(string = UniqueID, pattern="s[0-9]{1}-[a-z]{2}-[a-z]{1}[0-9]{1}-[0-9]{1,2}")) %>% 
  filter(CH4_quality.check=="") %>%
  select(plotID, lightCondition,CH4_best.flux) %>% 
  pivot_wider(id_cols = plotID,names_from = lightCondition, values_from = CH4_best.flux) %>% 
  ggplot(aes(x=transparent, y=dark))+
  geom_point()+
  geom_abline(slope = 1)
str(a)



####GHG veg-land vs veg-water####


#Load CH4 data after filter implemented table_ch4$HM.MAE<100
ch4<- read.csv(paste0(results_path,"all_good_CH4_fluxes_2023-10-03_to_2024-08-01.csv"))

#Load co2 data after filter implemented table_co2$Hm.MAE<10
co2<- read.csv(paste0(results_path,"all_good_CO2_fluxes_2023-10-03_to_2024-08-01.csv"))

#####1. Strata distribution####

#Check the strata distribution (veg-land vs veg-water) in each subsite
#I.e. for which subsites would it make a difference to separate?
#Use distribution from fieldsheets

fieldsheets<- list.files(path= fieldsheetpath, recursive = T, pattern = "GHG.xlsx",full.names = T)
field<- read_GHG_fieldsheets(fieldsheets)

lifewatch_example_path<- paste0(dropbox_root,"/GHG/Processed data/computed_flux/LifeWatch Italy GHG dataset example/")
incubations_without_flux<- read.csv(file = paste0(lifewatch_example_path,"Incubations_without_flux.csv"))

#How many incubations performed per strata and light-condition
chamber_deployments<- field %>% 
  select(-c(logger_floating_chamber,logger_transparent_chamber,chamber_type,longitude,latitude, start_time,end_time, initial_co2,initial_ch4,final_co2,final_ch4, chamber_height_cm, unix_start,unix_stop, plot_id,comments,date, person_sampling, gas_analyzer)) %>% 
  filter(!uniqID%in%incubations_without_flux$uniqID) %>% 
  mutate(campaign=substr(subsite, 1,2),
         status=substr(subsite, 7,7),
         sampling=subsite,
         subsite=substr(sampling,4,8),
         water_presence=case_when(water_depth>0~T,
                                  water_depth==0~F,
                                  TRUE~NA),
         vegetation_presence=case_when(strata=="vegetated"~T,
                                       TRUE~F),
         transparent_condition=case_when(transparent_dark=="dark"~F,
                                         transparent_dark=="transparent"~T)
         
  ) %>% 
  select(campaign, pilot_site, status, subsite,sampling, uniqID,water_depth, water_presence, vegetation_presence, transparent_condition)


chamber_deployments %>% 
  group_by(pilot_site, campaign, subsite, vegetation_presence) %>% 
  summarise(Yes=sum(water_presence),
            No=sum(!water_presence)) %>% 
  filter(vegetation_presence==T) %>% 
  pivot_longer(cols=c(No, Yes),names_to = "water_presence",values_to = "num_fluxes") %>% 
  rbind(tibble(pilot_site=c("RI","RI"), campaign=c("S1","S3"), subsite=c("RI-A1","RI-A1"), vegetation_presence=c(T,T), water_presence=c("No","Yes"), num_fluxes=c(0,0))) %>% 
  rowwise() %>% 
  mutate(subsite_only=gsub(paste0(pilot_site,"-"), "", subsite)) %>% 
  ggplot(aes(x=campaign, y=num_fluxes, fill=water_presence))+
  scale_fill_manual(values = c("No" = "brown", "Yes" = "blue")) +
  geom_bar(stat = "identity", position = "dodge")+
  scale_y_continuous(name = "Number of vegetated fluxes")+
  facet_grid(pilot_site~subsite_only, scales = "free")+
  theme_bw()+
  ggtitle("Water presence (depth>0cm) in vegetated plots (number of fluxes)")
  

#CASE PILOTS with uniform water presence in each subsite (for relevant fluxes): CU, DA, RI  

#CASE PILOTS with variable water presence: CA, DU, VA. 
#How relevant is the water variability? Is  just for a couple of chambers?

#CA: subsites with very mixed water presence 
#VA: subsites with mixed water presence, strong seasonality in some cases
#DU: most subsites with fairly uniform water absence, some with more mixed presence.


#DECISSION: separating vegetated-land and vegetated-water will be the best option for coherence across case-pilots IF there is any obvious effect of water presence on the fluxes of vegetated plots. 



#####2.GHG Veg-land vs Veg-water####

######2.1. CH4#####
#Overview CA, DU, VA CH4
#CA shows significant differences in CH4 emissions between vegetated land and vegetated water for most subsites, not so clear for DU and VA using the best.flux CH4 dataset (after filtering out non-reliable fluxes based on model fit)
ch4 %>% 
  filter(strata=="vegetated") %>% 
  filter(pilotsite%in%c("ca","du","va")) %>%
  mutate(water_presence=case_when(water_depth>0~"Yes",
                                  water_depth<=0~"No",
                                  T~NA)) %>% 
  ggplot(aes(x=water_presence,y=best.flux, col=water_presence))+
  geom_boxplot()+
  scale_color_manual(values = c("No" = "brown", "Yes" = "blue")) +
  stat_compare_means(label.sep = "\n",label.y.npc = 0.8)+
  scale_y_continuous(name="CH4 best.flux")+
  stat_summary(fun = length, aes(label = paste("n =", ..y..)), geom = "text", vjust = -2.5, size = 3) +
  theme_bw()+
  facet_grid(pilotsite~subsite,scales = "free")

######2.1. CO2#####
#Only a few subsites show significant differences in transparent CO2 flux, even considering that veg-water vegetation is usually much bigger

co2 %>% 
  filter(strata=="vegetated") %>% 
  filter(lightCondition=="transparent") %>% 
  filter(pilotsite%in%c("ca","du","va")) %>%
  mutate(water_presence=case_when(water_depth>0~"Yes",
                                  water_depth<=0~"No",
                                  T~NA)) %>% 
  ggplot(aes(x=water_presence,y=best.flux, col=water_presence))+
  geom_boxplot()+
  scale_color_manual(values = c("No" = "brown", "Yes" = "blue")) +
  stat_compare_means(label.sep = "\n",label.y.npc = 0.1)+
  scale_y_continuous(name="CO2 best.flux transparent", limits = c(-45,45))+
  stat_summary(fun = length, 
               aes(label = paste("n =", ..y..)), 
               geom = "text", 
               size = 3, 
               vjust = -2.5) +
  theme_bw()+
  facet_grid(pilotsite~subsite)

#CO2 dark
#As for transparent, only a few subsites show a significant difference, with veg-water ER being usually lower than veg-land.
co2 %>% 
  filter(strata=="vegetated") %>% 
  filter(lightCondition=="dark") %>% 
  filter(pilotsite%in%c("ca","du","va")) %>%
  mutate(water_presence=case_when(water_depth>0~"Yes",
                                  water_depth<=0~"No",
                                  T~NA)) %>% 
  ggplot(aes(x=water_presence,y=best.flux, col=water_presence))+
  geom_boxplot()+
  scale_color_manual(values = c("No" = "brown", "Yes" = "blue")) +
  stat_compare_means(label.sep = "\n",label.y.npc = 0.1)+
  scale_y_continuous(name="CO2 best.flux dark", limits = c(-30,30))+
  stat_summary(fun = length, 
               aes(label = paste("n =", ..y..)), 
               geom = "text", 
               size = 3, 
               vjust = -2.5) +
  theme_bw()+
  facet_grid(pilotsite~subsite)


###CH4 OW vs Veg-water####
#####1. Strata distribution####

#Check the strata distribution (Open-water vs veg-water) in each subsite
#Use distribution from fieldsheets

chamber_deployments %>% 
  group_by(pilot_site, campaign, subsite, water_presence) %>% 
  summarise(Yes=sum(vegetation_presence),
            No=sum(!vegetation_presence)) %>% 
  filter(water_presence==T) %>% 
  pivot_longer(cols=c(No, Yes),names_to = "vegetation_presence",values_to = "num_fluxes") %>% 
  # rbind(tibble(pilot_site=c("RI","RI"), campaign=c("S1","S3"), subsite=c("RI-A1","RI-A1"), vegetation_presence=c(T,T), water_presence=c("No","Yes"), num_fluxes=c(0,0))) %>% 
  rowwise() %>% 
  mutate(subsite_only=gsub(paste0(pilot_site,"-"), "", subsite)) %>% 
  ggplot(aes(x=campaign, y=num_fluxes, fill=vegetation_presence))+
  scale_fill_manual(values = c("No" = "lightblue", "Yes" = "darkgreen")) +
  geom_bar(stat = "identity", position = "dodge")+
  scale_y_continuous(name = "Number of Water fluxes")+
  facet_grid(pilot_site~subsite_only, scales = "free")+
  theme_bw()+
  ggtitle("Distribution of Vegetation presence in water plots (depth>0cm)")


######2.1. CH4#####
#CH4 from go.flux
#Significant increase due to vegetation_presence for CH4 fluxes in most subsites 
ch4 %>% 
  filter(strata!="bare") %>% 
  filter(water_depth>0) %>% 
  filter(pilotsite%in%c("ca","cu","da","du","va")) %>%
  mutate(vegetation_presence=case_when(strata=="vegetated"~"Yes",
                                       strata=="open water"~"No",
                                  T~NA)) %>% 
  ggplot(aes(x=vegetation_presence,y=best.flux, col=vegetation_presence))+
  geom_boxplot()+
  scale_color_manual(values = c("No" = "lightblue", "Yes" = "darkgreen")) +
  stat_compare_means(label.sep = "\n",label.y.npc = 0.8)+
  scale_y_continuous(name="CH4 best.flux")+
  stat_summary(fun = length, aes(label = paste("n =", ..y..)), geom = "text", vjust = -2.5, size = 3) +
  theme_bw()+
  facet_grid(pilotsite~subsite,scales = "free")


#CH4 ebullition (when detected)
#Significant increase due to vegetation_presence for CH4 ebullitive fluxes in almost all subsites
ch4 %>% 
  filter(strata!="bare") %>% 
  filter(water_depth>0) %>% 
  filter(pilotsite%in%c("ca","cu","da","du","va")) %>%
  mutate(vegetation_presence=case_when(strata=="vegetated"~"Yes",
                                       strata=="open water"~"No",
                                       T~NA)) %>% 
  filter(ebullition>0) %>% 
  ggplot(aes(x=vegetation_presence,y=ebullition, col=vegetation_presence))+
  geom_boxplot()+
  scale_color_manual(values = c("No" = "lightblue", "Yes" = "darkgreen")) +
  stat_compare_means(label.sep = "\n",label.y.npc = 0.8)+
  scale_y_continuous(name="CH4 ebullition")+
  stat_summary(fun = length, aes(label = paste("n =", ..y..)), geom = "text", vjust = -2.5, size = 3) +
  theme_bw()+
  facet_grid(pilotsite~subsite,scales = "free")


#CH4 total flux (ebullition+diffusion with camille aproach), Including cases without detectable ebullition 
#Significant increase due to vegetation_presence for CH4 total fluxes in almost all subsites
ch4 %>% 
  filter(strata!="bare") %>% 
  filter(water_depth>0) %>% 
  filter(pilotsite%in%c("ca","cu","da","du","va")) %>%
  mutate(vegetation_presence=case_when(strata=="vegetated"~"Yes",
                                       strata=="open water"~"No",
                                       T~NA)) %>% 
  mutate(total_ch4_flux=ebullition+diffusion) %>% 
  ggplot(aes(x=vegetation_presence,y=total_ch4_flux, col=vegetation_presence))+
  geom_boxplot()+
  scale_color_manual(values = c("No" = "lightblue", "Yes" = "darkgreen")) +
  stat_compare_means(label.sep = "\n",label.y.npc = 0.8)+
  scale_y_continuous(name="CH4 total flux (ebu+diff)")+
  stat_summary(fun = length, aes(label = paste("n =", ..y..)), geom = "text", vjust = -2.5, size = 3) +
  theme_bw()+
  facet_grid(pilotsite~subsite,scales = "free")



#Sampling dates per subsite####
samplingdates<- field %>% 
  select(subsite, date) %>% 
  separate(subsite, into = c("season", "site", "status_num"), sep = "-") %>% 
  mutate(subsite=paste(site, status_num,sep = "-")) %>% 
  select(season,subsite, date) %>% 
  distinct() %>% 
  filter(!(subsite=="RI-R2"&date%in%as.POSIXct(c("2023-11-15","2024-01-24","2024-04-10","2024-07-03"))))#Remove days for high-tide sampling only (3samples) of Restored site RI-R2, leaving the day when the bulk of sampling was performed.


write.csv(samplingdates, file = paste0(lifewatch_example_path,"Restore4Cs_ghg_sampling_dates_per_subsite.csv"),row.names = F)

