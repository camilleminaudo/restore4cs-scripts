
# ---
# Authors: Camille Minaudo
# Project: "RESTORE4Cs"
# date: "May 2024"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
# This script loads the fluxes estimates from script raw2flux.R. 
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
# files.sources = list.files(path = paste0(repo_root,"/functions"), full.names = T)
# for (f in files.sources){source(f)}


# ---- Directories and data loading ----
dropbox_root <- "C:/Users/Camille Minaudo/Dropbox/RESTORE4Cs - Fieldwork/Data" # You have to make sure this is pointing to the write folder on your local machine
results_path <- paste0(dropbox_root,"/GHG/Processed data/computed_flux/")

plots_path <- paste0(dropbox_root,"/GHG/Processed data/computed_flux/plots")


setwd(results_path)
listf <- list.files(path = results_path, pattern = ".csv", all.files = T, full.names = T, recursive = F)
table_results_all <- NULL
for (f in listf){
  table_results_all <- rbind(table_results_all, 
                             read.csv(file = f, header = T))
}

table_results_all$sampling <- str_sub(table_results_all$subsite, start = 1, 2)
table_results_all$pilotsite <- str_sub(table_results_all$subsite, start = 4, 5)
table_results_all$subsite <- str_sub(table_results_all$subsite, start = 7, 8)
table_results_all$siteID <- str_sub(table_results_all$subsite, start = 4, 8)


# table_results_all <- table_results_all[table_results_all$CO2_best.flux<1000,]


get_full_detailed_table <- function(table_results_all, variable){
  
  listf <- list.files(path = paste0(results_path,"level_incubation"), pattern = ".csv", all.files = T, full.names = T, recursive = F)
  mytable <- NULL
  for (f in listf[grep(pattern = variable, x = listf)]){
    mytable.tmp <- read.csv(file = f, header = T)
    if(dim(mytable.tmp)[2]==44){
      mytable <- rbind(mytable, mytable.tmp)
    }
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
  scale_x_log10()+
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
  scale_x_log10()+
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



table_co2$UniqueID[table_co2$HM.MAE>10]
table_ch4$UniqueID[table_ch4$HM.MAE>100]


# ---- plot CO2 and CH4 flux across sites and seasons ----

ind_sel <- which(table_co2$HM.MAE<10)
table_co2_sel <- table_co2[ind_sel,]

ind_sel <- which(table_ch4$HM.MAE<100)
table_ch4_sel <- table_ch4[ind_sel,]


p_overview_co2 <- ggplot(table_co2_sel, aes(pilotsite, best.flux, colour = lightCondition))+
  geom_hline(yintercept = 0)+
  geom_jitter(alpha=0.5)+
  # geom_boxplot(alpha=0.8)+
  theme_article()+
  ylab("CO2 flux mmol/m2/s")+
  scale_fill_viridis_d(begin = 0.2, end = 0.9)+
  scale_colour_viridis_d(begin = 0.2, end = 0.9, option = "C")+facet_wrap(strata~.)+theme(legend.position = 'top')+xlab("")


p_overview_ch4 <- ggplot(table_ch4_sel, aes(pilotsite, best.flux, colour = lightCondition))+
  geom_hline(yintercept = 0)+
  geom_jitter(alpha=0.5)+
  # geom_boxplot(alpha=0.8)+
  theme_article()+
  ylab("CH4 flux nmol/m2/s")+
  scale_fill_viridis_d(begin = 0.2, end = 0.9)+
  scale_colour_viridis_d(begin = 0.2, end = 0.9, option = "C")+facet_wrap(strata~.)+theme(legend.position = 'none')

ggarrange(p_overview_co2, p_overview_ch4)




ggplot(table_co2_sel, aes(subsite, best.flux, colour = lightCondition))+
  geom_hline(yintercept = 0)+
  geom_jitter(alpha=0.5)+
  geom_boxplot(alpha=0.8)+
  theme_article()+
  ylab("CO2 flux mmol/m2/s")+
  scale_fill_viridis_d(begin = 0.2, end = 0.9)+
  scale_colour_viridis_d(begin = 0.2, end = 0.9, option = "C")+facet_wrap(strata~., scales = "free_y")+theme(legend.position = 'top')+xlab("")



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


plot_overview(mytable = table_co2_sel, label = "CO2 flux mmol/m2/s", title = "CO2 flux", variable = "CO2")
plot_overview(mytable = table_ch4_sel, label = "CH4 flux nmol/m2/s", title = "CH4 flux", variable = "CH4")



plot_overview(mytable = table_co2_sel[which(table_co2_sel$strata=="open water" & table_co2_sel$lightCondition=='dark'),],
              label = "CO2 flux mmol/m2/s", title = "CO2 flux", variable = "CO2_open_water_dark")
plot_overview(mytable = table_ch4_sel[which(table_ch4_sel$strata=="open water" & table_ch4_sel$lightCondition=='dark'),],
              label = "CH4 flux nmol/m2/s", title = "CH4 flux", variable = "CH4_open_water_dark")



plot_overview(mytable = table_co2_sel[which(table_co2_sel$strata=="vegetated"),],
              label = "CO2 flux mmol/m2/s", title = "CO2 flux", variable = "CO2_vegetation")
plot_overview(mytable = table_ch4_sel[which(table_ch4_sel$strata=="vegetated"),],
              label = "CH4 flux nmol/m2/s", title = "CH4 flux", variable = "CH4_vegetation")





# ---- Only open water + dark ----
ind_sel <- which(table_results_all$water_depth>0 & table_results_all$lightCondition=='dark')

table_water_d <- table_results_all[ind_sel,]


ggplot(table_water_d, aes(subsite, CO2_best.flux, colour = sampling))+
  geom_jitter()+geom_boxplot()+
  theme_article()+facet_wrap(pilotsite~.)







# ---- Some plots ----


table_results_all$campaign_site <- substr(table_results_all$subsite,start = 1, stop = 5)
table_results_all$subsite_short <- substr(table_results_all$subsite,start = 7, stop = 8)


plt_CO2 <- ggplot(table_results_all, aes(subsite_short, CO2_LM.flux,
                                         fill = lightCondition))+
  geom_hline(yintercept = 0)+
  geom_boxplot(alpha=0.5)+
  # geom_jitter(width = 0.2, aes(colour = strata), size=2)+
  theme_article()+
  xlab("subsite")+
  ylab("CO2 flux mmol/m2/s")+
  ggtitle("CO2 flux")+
  scale_fill_viridis_d(begin = 0.2, end = 0.9)+
  scale_colour_viridis_d(begin = 0.2, end = 0.9, option = "C")+
  facet_grid(.~campaign_site)


plt_CH4diff <- ggplot(table_results_all, aes(subsite_short, CH4_LM.flux, 
                                             fill = lightCondition))+
  geom_hline(yintercept = 0)+
  geom_boxplot(alpha=0.5)+
  # geom_jitter(width = 0.2, aes(colour = strata), size=2)+
  theme_article()+
  xlab("subsite")+
  ylab("CH4 flux nmol/m2/s")+
  ggtitle("CH4 flux")+
  scale_fill_viridis_d(begin = 0.2, end = 0.9)+
  scale_colour_viridis_d(begin = 0.2, end = 0.9, option = "C")+
  facet_grid(.~campaign_site)

plt_all <- ggarrange(plt_CO2, plt_CH4diff, ncol = 1)

ggsave(plot = plt_all, filename = paste0(myfilename,".jpg"), path = plots_path, 
       width = 10, height = 8, dpi = 300, units = 'in')



for (cs in unique(table_results_all$campaign_site)){
  
  table_results_cs <- table_results_all[table_results_all$campaign_site == cs,]
  
  plt_CO2 <- ggplot(table_results_cs, aes(subsite_short, CO2_LM.flux,
                                          fill = lightCondition))+
    geom_hline(yintercept = 0)+
    geom_boxplot(alpha=0.5, width=0.5)+
    geom_jitter(width = 0.1, size=2, alpha=0.5, aes(shape = strata))+
    theme_article()+
    xlab("subsite")+
    ylab("CO2 flux mmol/m2/s")+
    ggtitle(paste0(cs, ""))+
    scale_fill_viridis_d(begin = 0.2, end = 0.9)+
    scale_colour_viridis_d(begin = 0.2, end = 0.9, option = "C")+
    # facet_grid(lightCondition~.)+
    theme(legend.position ='none')
  
  
  plt_CH4diff <- ggplot(table_results_cs, aes(subsite_short, CH4_LM.flux, 
                                              fill = lightCondition))+
    geom_hline(yintercept = 0)+
    geom_boxplot(alpha=0.5, width=0.5)+
    geom_jitter(width = 0.1, size=2, alpha=0.5, aes(shape = strata))+
    theme_article()+
    xlab("subsite")+
    ylab("CH4 flux nmol/m2/s")+
    # ggtitle("CH4 flux")+
    scale_fill_viridis_d(begin = 0.2, end = 0.9)+
    scale_colour_viridis_d(begin = 0.2, end = 0.9, option = "C")
  
  plt_cs <- ggarrange(plt_CO2, plt_CH4diff, ncol = 2)
  
  
  setwd(plots_path)
  myfilename <- paste(cs,"fluxes",min(as.Date(table_results_all$start.time)),"to",
                      max(as.Date(table_results_all$start.time)), sep = "_")
  
  ggsave(plot = plt_cs, filename = paste0(myfilename,".jpg"), path = plots_path, 
         width = 10, height = 5, dpi = 300, units = 'in')
}
