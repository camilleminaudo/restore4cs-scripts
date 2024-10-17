
# ---
# Authors: Camille Minaudo
# Project: "RESTORE4Cs"
# date: "Oct 2023"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
# This script


# my questions:
# 1- how different are CO2 fluxes estimates between expert timeseries selection and "blind" automation
# 2- how different are CH4 fluxes estimates between expert diffusion/ebullition separation and "blind" automation
# 3- how do experts usually select timeseries?




rm(list = ls()) # clear workspace
cat("/014") # clear console


# ---- Directories ----

# You have to make sure this is pointing to the write folder on your local machine
dropbox_root <- "C:/Users/Camille Minaudo/Dropbox/RESTORE4Cs - Fieldwork/Data" 


datapath <- paste0(dropbox_root,"/GHG/RAW data")
loggerspath <- paste0(datapath,"/RAW Data Logger")
RData_path <- paste0(dropbox_root,"/GHG/Processed data/RData/")
# plots_path <- "C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/meetings_presentations/2024_PPNW_Girona"
plots_path <- paste0(dropbox_root,"/GHG/GHG_expert_vs_automated/plots/")
results_path <- paste0(dropbox_root,"/GHG/GHG_expert_vs_automated/results/")

# ---- packages ----
library(tidyverse)
library(lubridate)
library(zoo)
library(ggplot2)
library(egg)
library(goFlux)
require(dplyr)
require(purrr)
require(msm)
require(data.table)
require(tools)
require(pbapply)


# ---- Loading data ----

load_fs <- function(path, pattern){
  
  list_f <- list.files(path = path, pattern = pattern, full.names = T)
  list_f <- list_f[grep(pattern = "csv", x = list_f)]
  
  isF <- T
  for (f in list_f){
    df_tmp <- read.csv(f)
    if(isF){
      isF <- F
      df_out <- df_tmp
    } else {
      df_out <- rbind(df_out,df_tmp)
    }
  }
  return(df_out)
}

table_all <- load_fs(path = results_path, pattern = "BLIND_vs_EXPERT_co2_ch4_fluxes_")
table_ebull <- load_fs(path = results_path, pattern = "ch4_ebullition")
table_ebull$ebullition[which(table_ebull$ebullition>table_ebull$total_estimated)] <- table_ebull$total_estimated[which(table_ebull$ebullition>table_ebull$total_estimated)]
table_ebull$ebullition[which(table_ebull$ebullition < 0)] <- 0


table_draws <- load_fs(path = results_path, pattern = "table_draw")
table_draws$diff_t_start_co2 <- table_draws$start.time_expert_co2-table_draws$start.time_auto
table_draws$diff_t_start_ch4 <- table_draws$start.time_expert_ch4-table_draws$start.time_auto
table_draws$diff_t_end_co2 <- table_draws$end.time_expert_co2-table_draws$end.time_auto
table_draws$diff_t_end_ch4 <- table_draws$end.time_expert_ch4-table_draws$end.time_auto
table_draws$duration_expert_co2 <- table_draws$end.time_expert_co2 - table_draws$start.time_expert_co2
table_draws$duration_expert_ch4 <- table_draws$end.time_expert_ch4 - table_draws$start.time_expert_ch4



# ---- overview database ----

sprd_CO2 <- table_all[table_all$variable=="CO2",c("variable","UniqueID","flux_method","best.flux","timestamp_processing","username")] %>%
  pivot_wider(names_from = flux_method, values_from = c(best.flux))

n <- length((sprd_CO2$UniqueID))
message(paste0("A total of ",n, " incubations were analysed"))

table_n <- data.frame(id = names(table(sprd_CO2$UniqueID)),
                      n = as.numeric(table(sprd_CO2$UniqueID)))


p_overview <- ggplot(table_n, aes(n))+geom_histogram()+theme_article()+xlab("Repetitions")+ylab("Counts")+
  ggtitle(paste0("n = ",n))

ggsave(plot = p_overview, filename = "overview_expert_vs_blind.jpeg", path = plots_path, 
       width = 5, height = 4, dpi = 300, units = 'in', scale = 0.8)




x <- seq(1,630)
already <- sort(unique(table_draws$draw))
areleft <- x[-already]

nb_draw <- 10
draw <- sample(areleft, nb_draw)



# ---- a function to get exceedance curves ----

get_df_exceedance <- function(relative_diff){
  relative_diff <- relative_diff[!is.na(relative_diff)]
  df_exceedance <- NULL
  for(thresh in sort(unique(relative_diff))){
    n_less <- length(which(abs(relative_diff) < thresh))
    df_exceedance <- rbind(df_exceedance,
                           data.frame(t = thresh,
                                      n = n_less,
                                      p = n_less/length(relative_diff)*100))
  }
  return(df_exceedance)
}



# ---- CO2 fluxes expert vs "blind" ----

sprd_CO2 <- table_all[table_all$variable=="CO2",c("variable","UniqueID","flux_method","best.flux","timestamp_processing","username")] %>%
  pivot_wider(names_from = flux_method, values_from = c(best.flux))
list_rm <- table_draws$UniqueID[which(table_draws$duration_expert_co2<100)]

ind_rmv <- match(list_rm, sprd_CO2$UniqueID)
sprd_CO2 <- sprd_CO2[-ind_rmv,]

sprd_CO2 <- sprd_CO2[sprd_CO2$Expert> -0.4,]

sprd_CO2$relative_diff <- abs((sprd_CO2$Expert-sprd_CO2$Blind)/sprd_CO2$Expert)*100

df_exceedance_CO2 <- get_df_exceedance(sprd_CO2$relative_diff)
ind_closest_10 <- which.min(abs(df_exceedance_CO2$t-10))


p_xy <- ggplot(data = sprd_CO2)+
  geom_abline(slope = 1,intercept = 0, color = 'grey')+
  geom_point(aes(Blind, Expert))+
  scale_x_log10()+
  scale_y_log10()+
  xlab("Blind automated CO2 flux [mmol/m2/s]")+
  ylab("Expert CO2 flux [mmol/m2/s]")+
  theme_bw()+ggtitle(paste0("CO2 flux",", difference < 10% for ",round(df_exceedance_CO2$p[ind_closest_10]),"% of the measurements"))

p_exceed <- ggplot(df_exceedance_CO2, aes(t, p))+geom_path()+geom_point()+
  geom_hline(yintercept = 0)+
  theme_bw()+
  geom_segment(aes(x=-0,xend=10,
                   y=df_exceedance_CO2$p[ind_closest_10], yend=df_exceedance_CO2$p[ind_closest_10]), color ="red")+
  geom_segment(aes(x=10,xend=10,
                   y=-Inf, yend=df_exceedance_CO2$p[ind_closest_10]), color ="red")+
  # scale_x_log10()+
  scale_x_continuous(transform = "log10", breaks = c(1,10,100,1000))+
  ylab("Proportion of timeseries [%]")+
  xlab("Relative difference [% of Expert flux]")

p_exceed_co2 <- ggarrange(p_xy, p_exceed, nrow = 1)


ggsave(plot = p_exceed_co2, filename = "overview_expert_vs_blind_co2.jpeg", path = plots_path, 
       width = 8, height = 4, dpi = 300, units = 'in', scale = 1.2)



sprd_CO2$UniqueID[which(sprd_CO2$Blind>0.8)]

sprd_CO2$UniqueID[which(sprd_CO2$Blind < -.4)]



sprd_CO2$absolute_diff <- abs(sprd_CO2$Expert - sprd_CO2$Blind)

table_all_co2 <- table_all[table_all$variable=="CO2",]
sprd_CO2$LM.RMSE <- table_all_co2$LM.RMSE[match(sprd_CO2$UniqueID, table_all_co2$UniqueID)]


p_diff_LM.RMSE <- ggplot(sprd_CO2, aes(LM.RMSE, absolute_diff))+geom_point()+
  # geom_abline(slope = 1, intercept = 0)+
  scale_x_log10()+
  scale_y_log10()+
  theme_article()+
  xlab("RMSE Linear model [mmol/m2/s]")+
  ylab("abs(Expert - Automated) [mmol/m2/s]")+
  # ggtitle("Poor linear fits correspond to largest absolute differences")+
  theme(legend.position = "top")


ggsave(plot = p_diff_LM.RMSE, filename = "overview_expert_vs_blind_co2_LM_RMSE.jpeg", path = plots_path, 
       width = 8, height = 4, dpi = 300, units = 'in', scale = 1.2)






# ---- CH4 fluxes expert vs "blind" ----



ggplot(table_ebull, aes(best.flux, total_estimated))+
  geom_abline(slope = 1,intercept = 0, color = 'grey')+
    geom_point()+
  scale_x_log10()+
  scale_y_log10()+
    xlab("Blind automated CH4 flux [nmol/m2/s]")+
    ylab("Expert CH4 flux [nmol/m2/s]")+
    theme_bw()+ggtitle("CH4 total flux")




# ---- CH4 diffusion expert vs "dydt" ----

sprd_diffus <- table_ebull[,c("UniqueID","flux_method","diffusion","timestamp_processing","username")] %>%
  pivot_wider(names_from = flux_method, values_from = c(diffusion))
list_rm <- table_draws$UniqueID[which(table_draws$duration_expert_ch4<10)]

ind_rmv <- match(list_rm, sprd_diffus$UniqueID)
sprd_diffus <- sprd_diffus[-ind_rmv,]


sprd_diffus$relative_diff <- abs((sprd_diffus$Expert-sprd_diffus$dydt)/sprd_diffus$Expert)*100

df_exceedance_CH4 <- get_df_exceedance(sprd_diffus$relative_diff)
ind_closest_10 <- which.min(abs(df_exceedance_CH4$t-10))

p_xy <- ggplot(data = sprd_diffus)+
  geom_abline(slope = 1,intercept = 0, color = 'grey')+
  geom_point(aes(dydt, Expert))+
  scale_x_log10()+
  scale_y_log10()+
  xlab("Automated CH4 diffusive flux [nmol/m2/s]")+
  ylab("Expert CH4 diffusive flux [nmol/m2/s]")+
  theme_bw()+ggtitle(paste0("CH4 diffusion",", difference < 10% for ",round(df_exceedance_CH4$p[ind_closest_10]),"% of the measurements"))

p_exceed <- ggplot(df_exceedance_CH4, aes(t, p))+geom_path()+geom_point()+
  geom_hline(yintercept = 0)+
  theme_bw()+
  geom_segment(aes(x=-0,xend=10,
                   y=df_exceedance_CH4$p[ind_closest_10], yend=df_exceedance_CH4$p[ind_closest_10]), color ="red")+
  geom_segment(aes(x=10,xend=10,
                   y=-Inf, yend=df_exceedance_CH4$p[ind_closest_10]), color ="red")+
  # scale_x_log10()+
  scale_x_continuous(transform = "log10", breaks = c(1,10,100,1000))+
  ylab("Proportion of timeseries [%]")+
  xlab("Relative difference [% of Expert flux]")

ggarrange(p_xy, p_exceed, nrow = 1)


p_exceed_ch4_diff <- ggarrange(p_xy, p_exceed, nrow = 1)

ggsave(plot = p_exceed_ch4_diff, filename = "overview_expert_vs_blind_ch4_diff.jpeg", path = plots_path, 
       width = 8, height = 4, dpi = 300, units = 'in', scale = 1.2)


# ---- CH4 ebullition expert vs "dydt" ----

sprd_ebull <- table_ebull[,c("UniqueID","flux_method","ebullition","timestamp_processing","username")] %>%
  pivot_wider(names_from = flux_method, values_from = c(ebullition))
list_rm <- table_draws$UniqueID[which(table_draws$duration_expert_ch4<10)]

ind_rmv <- match(list_rm, sprd_ebull$UniqueID)
sprd_ebull <- sprd_ebull[-ind_rmv,]



sprd_ebull$relative_diff <- abs((sprd_ebull$Expert-sprd_ebull$dydt)/sprd_ebull$Expert)*100

df_exceedance_CH4 <- get_df_exceedance(sprd_ebull$relative_diff)
ind_closest_10 <- which.min(abs(df_exceedance_CH4$t-10))


p_xy <- ggplot(data = sprd_ebull)+
  geom_abline(slope = 1,intercept = 0, color = 'grey')+
  geom_point(aes(dydt, Expert))+
  scale_x_log10()+
  scale_y_log10()+
  xlab("Automated CH4 ebullition [nmol/m2/s]")+
  ylab("Expert CH4 ebullition [nmol/m2/s]")+
  theme_bw()+ggtitle(paste0("CH4 ebullition",", difference < 10% for ",round(df_exceedance_CH4$p[ind_closest_10]),"% of the measurements"))

p_exceed <- ggplot(df_exceedance_CH4, aes(t, p))+geom_path()+geom_point()+
  geom_hline(yintercept = 0)+
  theme_bw()+
  geom_segment(aes(x=-0,xend=10,
                   y=df_exceedance_CH4$p[ind_closest_10], yend=df_exceedance_CH4$p[ind_closest_10]), color ="red")+
  geom_segment(aes(x=10,xend=10,
                   y=-Inf, yend=df_exceedance_CH4$p[ind_closest_10]), color ="red")+
  # scale_x_log10()+
  scale_x_continuous(transform = "log10", breaks = c(1,10,100,1000))+
  ylab("Proportion of timeseries [%]")+
  xlab("Relative difference [% of Expert flux]")

p_exceed_ch4_ebull <- ggarrange(p_xy, p_exceed, nrow = 1)

ggsave(plot = p_exceed_ch4_ebull, filename = "overview_expert_vs_blind_ch4_ebull.jpeg", path = plots_path, 
       width = 8, height = 4, dpi = 300, units = 'in', scale = 1.2)




# ---- CH4 total expert vs "dydt" ----

sprd_total <- table_ebull[,c("UniqueID","flux_method","total_estimated","timestamp_processing","username")] %>%
  pivot_wider(names_from = flux_method, values_from = c(total_estimated))
list_rm <- table_draws$UniqueID[which(table_draws$duration_expert_ch4<10)]

ind_rmv <- match(list_rm, sprd_total$UniqueID)
sprd_total <- sprd_total[-ind_rmv,]
sprd_total <- sprd_total[which(sprd_total$dydt<1500),]


sprd_total$relative_diff <- abs((sprd_total$Expert-sprd_total$dydt)/sprd_total$Expert)*100

df_exceedance_CH4 <- get_df_exceedance(sprd_total$relative_diff)
ind_closest_10 <- which.min(abs(df_exceedance_CH4$t-10))

p_xy <- ggplot(data = sprd_total)+
  geom_abline(slope = 1,intercept = 0, color = 'grey')+
  geom_point(aes(dydt, Expert))+
  xlab("Automated CH4 total flux [nmol/m2/s]")+
  ylab("Expert CH4 total flux [nmol/m2/s]")+
  theme_bw()+ggtitle(paste0("CH4 total flux",", difference < 10% for ",round(df_exceedance_CH4$p[ind_closest_10]),"% of the measurements"))

p_exceed <- ggplot(df_exceedance_CH4, aes(t, p))+geom_path()+geom_point()+
  geom_hline(yintercept = 0)+
  theme_bw()+
  geom_segment(aes(x=-0,xend=10,
                   y=df_exceedance_CH4$p[ind_closest_10], yend=df_exceedance_CH4$p[ind_closest_10]), color ="red")+
  geom_segment(aes(x=10,xend=10,
                   y=-Inf, yend=df_exceedance_CH4$p[ind_closest_10]), color ="red")+
  # scale_x_log10()+
  scale_x_continuous(transform = "log10", breaks = c(1,10,100,1000))+
  ylab("Proportion of timeseries [%]")+
  xlab("Relative difference [% of Expert flux]")

p_exceed_ch4_tot <- ggarrange(p_xy, p_exceed, nrow = 1)

ggsave(plot = p_exceed_ch4_tot, filename = "overview_expert_vs_blind_ch4.jpeg", path = plots_path, 
       width = 8, height = 4, dpi = 300, units = 'in', scale = 1.2)




sprd_total$UniqueID[which(sprd_total$Expert< (-5))]




sprd_total$absolute_diff <- abs(sprd_total$Expert - sprd_total$dydt)

sprd_total$ebullition_expert <- sprd_ebull$Expert[match(sprd_total$UniqueID, sprd_ebull$UniqueID)]
sprd_total$diffusion_expert <- sprd_diffus$Expert[match(sprd_total$UniqueID, sprd_diffus$UniqueID)]
sprd_total$LM.RMSE <- table_ebull$LM.RMSE[match(sprd_total$UniqueID, table_ebull$UniqueID)]
sprd_total$LM.r2 <- table_ebull$LM.r2[match(sprd_total$UniqueID, table_ebull$UniqueID)]



p_diff_LM.RMSE <- ggplot(sprd_total[sprd_total$Expert>0,], aes(LM.RMSE, absolute_diff, size = ebullition_expert))+geom_point()+
  # geom_abline(slope = 1, intercept = 0)+
  scale_x_log10()+
  scale_y_log10()+
  theme_article()+
  xlab("RMSE Linear model [nmol/m2/s]")+
  ylab("abs(Expert - Automated) [nmol/m2/s]")+
  scale_size("Ebullition [nmol/m2/s]")+
  # ggtitle("Poor linear fits correspond to largest absolute differences")+
  theme(legend.position = "top")

ggsave(plot = p_diff_LM.RMSE, filename = "overview_expert_vs_blind_LM_RMSE.jpeg", path = plots_path, 
       width = 8, height = 4, dpi = 300, units = 'in', scale = 1.2)





ggplot(sprd_total[sprd_total$Expert>0,], aes(Expert, diffusion_expert, size = ebullition_expert))+geom_point()+
  # geom_abline(slope = 1, intercept = 0)+
  scale_x_log10()+
  scale_y_log10()+
  theme_article()+
  xlab("Expert total [nmol/m2/s]")+
  ylab("Expert diffusive [nmol/m2/s]")+
  scale_size("Expert ebullition [nmol/m2/s]")+
  # ggtitle("Poor linear fits correspond to largest absolute differences")+
  theme(legend.position = "top")



# 
# 
# p_prop <- ggplot(sprd_ebull[order(sprd_ebull$Expert),], 
#                  aes(seq_along(sprd_ebull$Expert)/dim(sprd_ebull)[1]*100,Expert))+
#   geom_hline(yintercept = 100)+
#   geom_hline(yintercept = 50, alpha=0.2)+
#   # geom_jitter(alpha=0.5)+
#   geom_point(alpha=0.5, width=0.3, aes(colour = sprd_ebull$Expert[order(sprd_ebull$Expert)] > 50))+
#   theme_article()+
#   # ylab("CH4 flux nmol/m2/s")+
#   scale_fill_viridis_d(begin = 0.2, end = 0.9)+
#   scale_colour_viridis_d("over 50%",begin = 0.2, end = 0.9, option = "A", direction = -1)+
#   ylim(c(0,120))+
#   # scale_y_log10()+
#   # facet_wrap(strata~.)+
#   theme(legend.position = 'none')+
#   xlab("Proportion of timeseries [%]")+
#   ylab("Contribution of bubbling to total CH4 flux [%]")+
#   ggtitle(paste0("Bubbling > 50% of F_total for ",round(n_50/dim(sprd_ebull)[1]*100),"% of the measurements"))





table_ch4_ebull <- sprd_total[!is.na(sprd_total$ebullition_expert),]
table_ch4_ebull$p_ebullition <- table_ch4_ebull$ebullition_expert/table_ch4_ebull$Expert *100

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
  theme(legend.position = 'right')+
  xlab("Proportion of timeseries [%]")+
  ylab("Contribution of bubbling to total CH4 flux [%]")


p_ebull <- ggplot(table_ch4_ebull, 
                  aes(" ", Expert))+
  # geom_hline(yintercept = 100)+
  geom_jitter(alpha=0.5, width = 0.2, aes(colour = table_ch4_ebull$p_ebullition > 50))+
  geom_boxplot(alpha=0.8, width=0.2)+
  scale_colour_viridis_d("over 50%",begin = 0.2, end = 0.9, option = "A", direction = -1)+
  theme_article()+
  scale_y_log10()+
  ylab("CH4 total flux [nmol/m2/s]")+
  theme(legend.position = 'none')+
  xlab("")#+ggtitle(paste0("Bubbling > 50% of total flux for ",round(n_50/dim(table_ch4_ebull)[1]*100),"% of the measurements"))


p_contrib_ebull_expert <- ggarrange(p_ebull, p_prop, nrow = 1)

ggsave(plot = p_contrib_ebull_expert, filename = "overview_expert_vs_blind_contrib_ebull.jpeg", path = plots_path, 
       width = 8, height = 4, dpi = 300, units = 'in', scale = 1.2)






sprd_total$ebullition_blind <- sprd_ebull$dydt[match(sprd_total$UniqueID, sprd_ebull$UniqueID)]
sprd_total$diffusion_blind <- sprd_diffus$dydt[match(sprd_total$UniqueID, sprd_diffus$UniqueID)]


ggplot(sprd_total[sprd_total$dydt>0,], aes(dydt, diffusion_blind, size = ebullition_blind))+geom_point()+
  # geom_abline(slope = 1, intercept = 0)+
  scale_x_log10()+
  scale_y_log10()+
  theme_article()+
  xlab("Expert total [nmol/m2/s]")+
  ylab("Expert diffusive [nmol/m2/s]")+
  scale_size("Expert ebullition [nmol/m2/s]")+
  # ggtitle("Poor linear fits correspond to largest absolute differences")+
  theme(legend.position = "top")




table_ch4_ebull <- sprd_total[!is.na(sprd_total$ebullition_blind),]
table_ch4_ebull$p_ebullition <- table_ch4_ebull$ebullition_blind/table_ch4_ebull$dydt *100

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
  theme(legend.position = 'right')+
  xlab("Proportion of timeseries [%]")+
  ylab("Contribution of bubbling to total CH4 flux [%]")


p_ebull <- ggplot(table_ch4_ebull, 
                  aes(" ", Expert))+
  # geom_hline(yintercept = 100)+
  geom_jitter(alpha=0.5, width = 0.2, aes(colour = table_ch4_ebull$p_ebullition > 50))+
  geom_boxplot(alpha=0.8, width=0.2)+
  scale_colour_viridis_d("over 50%",begin = 0.2, end = 0.9, option = "A", direction = -1)+
  theme_article()+
  scale_y_log10()+
  ylab("CH4 total flux [nmol/m2/s]")+
  theme(legend.position = 'none')+
  xlab("")+
  ggtitle(paste0("Bubbling > 50% of total flux for ",round(n_50/dim(table_ch4_ebull)[1]*100),"% of the measurements"))


p_contrib_ebull_blind <- ggarrange(p_ebull, p_prop, nrow = 1)





# ---- Differences between experts ----
 
ind_multiple_assess <- which(duplicated(table_draws$UniqueID))
table_draws$UniqueID[ind_multiple_assess]


list_multiple_experts <- unique(table_draws$UniqueID[ind_multiple_assess])
length(list_multiple_experts)

isF <- T
runID <- 0
for(i in list_multiple_experts){
  runID = runID+1
  
  table__i <- table_all[which(table_all$UniqueID == i & table_all$flux_method=="Expert"),]
  n = dim(table__i)[1]
  
  for(user in unique(table__i$username)){
    for (var in c("CO2","CH4")){
      df.tmp <- data.frame(id = i,
                           runID = as.factor(runID),
                           var = var,
                           n = n,
                           user = user,
                           n_user = dim(table__i[which(table__i$variable==var & table__i$username == user),])[1],
                           flux_mean = mean(table__i$best.flux[which(table__i$variable==var & table__i$username == user)]),
                           flux_sd = sd(table__i$best.flux[which(table__i$variable==var & table__i$username == user)]),
                           flux_mean_overall = mean(table__i$best.flux[which(table__i$variable==var)]),
                           flux_sd_overall = sd(table__i$best.flux[which(table__i$variable==var)]))
      
      if (isF){
        isF = F
        df_multiple_users <- df.tmp
      } else {
        df_multiple_users <- rbind(df_multiple_users, df.tmp)
      }
    }
  }
}



df_multiple_users$flux_CV_overall <- df_multiple_users$flux_sd_overall/df_multiple_users$flux_mean_overall


# plot(df_multiple_users$flux_mean_overall, df_multiple_users$flux_sd_overall)

df_multiple_users_co2 <- df_multiple_users[df_multiple_users$var=="CO2",]
df_multiple_users_co2$runID <- factor(df_multiple_users_co2$runID, levels = unique(df_multiple_users_co2$runID[order(df_multiple_users_co2$flux_CV_overall)]))
p_CO2 <- ggplot(df_multiple_users_co2[which(df_multiple_users_co2$n>3),], aes(runID, flux_mean))+geom_boxplot()+
  # coord_flip()+
  theme_article()+facet_wrap(var~., scales = "free")

ggplot(df_multiple_users_co2[which(df_multiple_users_co2$n>3),], aes(runID, flux_sd_overall))+geom_point()+
  # coord_flip()+
  theme_article()+facet_wrap(var~., scales = "free")


df_multiple_users_ch4 <- df_multiple_users[df_multiple_users$var=="CH4",]
df_multiple_users_ch4$runID <- factor(df_multiple_users_ch4$runID, levels = unique(df_multiple_users_ch4$runID[order(df_multiple_users_ch4$flux_sd_overall)]))
p_CH4 <- ggplot(df_multiple_users_ch4[which(df_multiple_users_ch4$n>3),], aes(runID, flux_mean))+geom_boxplot()+
  # coord_flip()+
  theme_article()+facet_wrap(var~., scales = "free")

ggarrange(p_CO2, p_CH4, ncol = 1)

sprd_multipl <- df_multiple_users[which(df_multiple_users$n>3),c("var","id","n","flux_mean","user", "flux_mean_overall", "flux_sd_overall")] %>%
  pivot_wider(names_from = user, values_from = c(flux_mean))

ggplot(sprd_multipl, aes(abs(flux_mean_overall), flux_sd_overall/abs(flux_mean_overall)*100))+
  geom_point()+theme_article()+facet_wrap(var~., scales="free")+
  # scale_x_log10()+
  xlab("Average flux across all replicates")+
  ylab("CV among replicates [% of average]")




# ---- Stats in timeseries cuts ----


ggplot(table_draws, aes(draw))+geom_density()

table_draws$diff_t_start_co2 <- table_draws$start.time_expert_co2-table_draws$start.time_auto
table_draws$diff_t_start_ch4 <- table_draws$start.time_expert_ch4-table_draws$start.time_auto
table_draws$diff_t_end_co2 <- table_draws$end.time_expert_co2-table_draws$end.time_auto
table_draws$diff_t_end_ch4 <- table_draws$end.time_expert_ch4-table_draws$end.time_auto


ggplot(table_draws, aes(diff_t_start_co2))+geom_density()

table_draws$duration_expert_co2 <- table_draws$end.time_expert_co2 - table_draws$start.time_expert_co2

ggplot(table_draws, aes(reorder(draw, duration_expert_co2, FUN=min), duration_expert_co2))+geom_point()+theme_article()

table_draws$UniqueID[which(table_draws$duration_expert_co2 < 100)]


# ---- Stats in fluxes differences ----

table_results_sprd <- table_all[,c("variable","UniqueID","flux_method","best.flux","timestamp_processing","username")] %>%
  pivot_wider(names_from = flux_method, values_from = c(best.flux))

ggplot(data = table_results_sprd)+
  # geom_abline(slope = 0,intercept = 0, color = 'lightgrey')+
  # geom_segment(data = data.frame(UniqueID = CO2_flux_res_auto$UniqueID,
  #                                meth1 = CO2_flux_res_auto$best.flux,
  #                                meth2 = CO2_flux_res_manID$best.flux), aes(x=UniqueID, xend=UniqueID, y = meth1, yend = meth2), linewidth=1, alpha = 0.5)+
  geom_violin(aes(variable,abs((Blind-Expert)/Expert)*100), alpha = 0.5)+
  geom_jitter(aes(variable,abs((Blind-Expert)/Expert)*100), alpha = 0.5, width=0.2)+
  ylab("CO2 flux relative difference [%]")+
  theme_bw()+
  # scale_y_log10()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))






plot(sort(table(paste0(table_results_sprd$variable,table_results_sprd$UniqueID))))

plot(sort(table(table_results_sprd$username)))


# ---- Ebullition ----

ggplot(data = table_ebull[which(table_ebull$flux_method=="Expert" & table_ebull$diffusion>0 & table_ebull$ebullition>0),])+ # [table_ebull$flux_method=="Expert",]
  # geom_abline(slope = 1,intercept = 0, color = 'grey')+
  # geom_segment(data = data.frame(UniqueID = CH4_res_meth1$UniqueID,
  #                                meth1 = CH4_res_meth1$ebullition,
                                 # meth2 = CH4_res_meth2$ebullition), aes(x=UniqueID, xend=UniqueID, y = meth1, yend = meth2), linewidth=1, alpha = 0.5)+
  geom_point(aes(diffusion, ebullition,
                 colour = log10(total_estimated)), size=4, alpha = 0.5)+
  # ylab("ebullition component [nmol/m2/s]")+
  scale_x_log10()+
  scale_y_log10()+
  theme_bw()+ 
  # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_colour_viridis_c(begin = 0.1, end = 0.9, option = "F", direction = -1)


ggplot(table_ebull[which(table_ebull$flux_method=="Expert" & table_ebull$diffusion>0 & table_ebull$ebullition>0),])+
  geom_density(aes(ebullition/diffusion, fill="ebullition/diffusion"))+
  # geom_boxplot(aes(ebullition, fill="ebullition"))+
  theme_article()+coord_flip()


table_ebull$UniqueID[which(table_ebull$total_estimated>1000)]

















# ---- Co2 vs CH4 ----

# draw <- sample(seq_along(auxfile$subsite), nb_draw)
draw <- seq_along(auxfile$subsite)

# draw <- which(auxfile$UniqueID == "s1-cu-a2-3-o-d-06:59")

table_draw <- data.frame(username = username,
                         draw = draw,
                         subsite = auxfile$subsite[draw],
                         UniqueID = auxfile$UniqueID[draw])
myauxfile <- auxfile[draw,]
myauxfile$username <- username

mydata_all <- NULL
for(k in seq_along(myauxfile$UniqueID)){
  message("Loading data for ",myauxfile$UniqueID[k])
  gas <- unique(myauxfile$gas_analiser[k])
  
  setwd(RData_path)
  if(gas== "LI-COR"){
    gs_suffix <- "LI-7810"
  } else {
    gs_suffix <- gas
  }
  load(file = paste0(myauxfile$subsite[k],"_",gs_suffix,".RData"))
  mydata <- mydata[,c("POSIX.time", 
                      "CO2dry_ppm", "CH4dry_ppb", "H2O_ppm",
                      "CO2_prec",   "CH4_prec",   "H2O_prec"  )]
  
  mydata <- mydata[which(mydata$POSIX.time>=myauxfile$start.time[k] & 
                           mydata$POSIX.time<=myauxfile$start.time[k]+myauxfile$duration[k]),]
  mydata$UniqueID <- myauxfile$UniqueID[k]
  mydata$Etime <- mydata$POSIX.time - min(mydata$POSIX.time)
  
  table_draw$corr_co2_ch4[k] <- cor(mydata$CO2dry_ppm, mydata$CH4dry_ppb)
  mod_lm <- lm(data = mydata, formula = CH4dry_ppb~CO2dry_ppm)
  table_draw$slope_co2_ch4[k] <- coefficients(mod_lm)[2]
  table_draw$r2_co2_ch4[k] <- summary(mod_lm)$adj.r.squared
  
  mydata_all <- rbind(mydata_all, mydata)
  rm(mydata)
}



# ggplot(mydata_all, aes(Etime, CH4dry_ppb))+geom_point()+geom_smooth(method = 'lm')+facet_wrap(UniqueID~.)+theme_article()
# ggplot(mydata_all, aes(CO2dry_ppm, CH4dry_ppb))+geom_point()+geom_smooth(method = 'lm')+facet_wrap(UniqueID~.)+theme_article()



ggplot(table_draw, aes(corr_co2_ch4))+geom_density()

ggplot(table_draw, aes(slope_co2_ch4, corr_co2_ch4))+geom_point()


table_draw$UniqueID[which(table_draw$slope_co2_ch4>5000)]

