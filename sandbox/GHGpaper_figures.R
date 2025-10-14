#GHGpaper_figures



#Author: Miguel Cabrera-Brufau
#Date: September 2025
#Project: Restore4cs

#Description----
#This scrip is used to produce the figures that will be included in the mainGHG paper (main + supplementary), as well a other exploratory figures. 

#Inputs:  In situ data + outputs of modelling
#ChamberData4paper.csv
#Allghg_emmeans_posthoc_chambermodels.csv
# Stratacomposition_in-situ.csv

#Outputs: 
#FIGURES


rm(list = ls()) # clear workspace
cat("/014") # clear console


# ---- Packages ----
#For data-handling & plotting
library(tidyverse)
library(readxl)
library(lubridate)
library(zoo)
# library(ggplot2)
# library(ggpubr)
# library(grid)
# library(egg)
require(purrr)
require(data.table)
require(tools)
library(hms)
library(suncalc)
library(cowplot) #Alligned axis
library(ggforce) #geom_sina


#Repo-functions
repo_root <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
files.sources = list.files(path = paste0(repo_root,"/functions"), full.names = T)
for (f in files.sources){source(f)}




# ---- Directories ----
dropbox_root <- "C:/Users/Miguel/Dropbox/RESTORE4Cs - Fieldwork/Data" # You have to make sure this is pointing to the write folder on your local machine

#Path to save datasets for main paper: 
paper_path<- paste0(dropbox_root,"/GHG/Main_paper_results/")

#Folders to save plots:
  #MAIN
main_figures<- paste0(paper_path,"Main_figures/")
  #Supplementary
sup_figures<- paste0(paper_path,"Sup_figures/")
  #Extra/exploration/miscelanea
extra_figures<- paste0(paper_path,"Exploratory_figures/")


#0. Import and format------

#Import Chamberdata4paper.csv
data4paper<-read.csv(file = paste0(paper_path,"ChamberData4paper.csv"))

#Format main data: 
data4paper<- data4paper %>% 
  #Remove NAs
  filter(!is.na(dailyflux)) %>% 
  #Keep only ghg of interest
  filter(ghgspecies%in%c("co2","ch4","gwp_co2andch4")) %>% 
  #Rename to GWPco2andch4
  mutate(ghgspecies=if_else(ghgspecies=="gwp_co2andch4","GWPco2andch4",ghgspecies)) %>% 
  #Factor grouping variables:
  mutate(season=factor(season, levels = c("S1","S2","S3","S4"), ordered = T),
         casepilot=factor(casepilot, levels = c("DU","RI","CA","VA","DA","CU"), ordered = T),
         status=factor(status, levels = c("Preserved","Altered","Restored"), ordered = T),
         subsite=factor(subsite, ordered = F),
         ghgspecies=factor(ghgspecies, ordered=F),
         strata=factor(strata, levels = c("open water","bare","vegetated"), ordered = F),
         sampling=factor(sampling, ordered = F))


##general models------
#Import results from models (R2 and significances)

summary_models_all<- read.csv(paste0(paper_path, "Allghg_summary_chambermodels.csv"))

significances_alleffects<- summary_models_all %>% 
  dplyr::select(dataset, status_pval, season_pval, interaction_pval) %>% 
  separate(dataset, into = c("casepilot","ghgspecies")) %>% 
  pivot_longer(c(status_pval,season_pval,interaction_pval),names_sep = "_",values_to = "pvalue_effect", names_to = c("effect","drop")) %>% 
  dplyr::select(-drop)
  
head(significances_alleffects)

#Import Emmeans post-hoc (all ghgs Combined into single table)
emmeans_all<-read.csv(paste0(paper_path,"Allghg_emmeans-posthoc_chambermodels.csv"))

#Format Emmeans: remove groupletters from comparisons when effect was not significant. 
  #For status_within_season group-letters, use the significance of "interaction" 

emmeans_all<- emmeans_all %>% 
  mutate(effect=if_else(comparison=="status_within_season","interaction",comparison)) %>% 
  #Join signifcances of model effects
  left_join(significances_alleffects , by=c("casepilot","ghgspecies","effect")) %>% 
  #Replace post-hoc letters with "" to avoid non-significant effect comparisons to be shown
  #We are leaving the pvalue threshold at 0.1, (va status_within_season for ch4 is almmost significant)
  mutate(group_letter=if_else(pvalue_effect>0.1,"",group_letter)) %>%
  mutate(ghgspecies=factor(ghgspecies, ordered=F),
         casepilot=factor(casepilot, levels = c("DU","RI","CA","VA","DA","CU"), ordered = T),
         status=factor(status, levels = c("Preserved","Altered","Restored"), ordered = T),
         season=factor(season, levels = c("S1","S2","S3","S4"), ordered = T)) %>% 
  arrange(casepilot, ghgspecies, comparison, season,status)


#separate emmeans per type of comparison
emmeans_status<- emmeans_all %>% 
  filter(comparison=="status") %>% 
  group_by(ghgspecies,casepilot) %>% 
  #If there are no differences (same letter for all status of a given ghgspecies*casepilot combination), remove the group_letter. 
  mutate(group_letter = if (n_distinct(group_letter) == 1) "" else group_letter) %>%
  ungroup()


emmeans_season<- emmeans_all %>% 
  filter(comparison=="season") %>% 
  group_by(ghgspecies,casepilot) %>% 
  #If there are no differences (same letter for all status of a given ghgspecies*casepilot combination), remove the group_letter. 
  mutate(group_letter = if (n_distinct(group_letter) == 1) "" else group_letter) %>%
  ungroup()


#Subset the emmeans and group-letters of status_within_season comparison. 
emmeans_status_within_season<- emmeans_all %>% 
  filter(comparison=="status_within_season") %>% 
  group_by(ghgspecies, casepilot, season) %>%
  #If there are no differences between status within a given combination of ghgspecies*casepilot*season, remove the group_letter. This avoids printing group_letters when no difference was seen. 
  mutate(group_letter = if (n_distinct(group_letter) == 1) "" else group_letter) %>%
  ungroup()


#Import emmeans contrasts
emmeans_contrasts<- read.csv(paste0(paper_path,"Allghg_pairwise_diffemmeans-posthoc_chambermodels.csv"))

emmeans_contrasts

#Function to get symbols for pvalues
pval_to_symbol <- function(p) {
  ifelse(p < 0.001, "(***)",
         ifelse(p < 0.01, "(**)",
                ifelse(p < 0.05, "(*)",
                       ifelse(p < 0.1, "(.)", "(ns)")
                )
         )
  )
}

emmeans_contrasts<- emmeans_contrasts %>% 
  mutate(pval_symbol=pval_to_symbol(p.value))


# formated_pvalues_all<- significances_alleffects %>% 
#   mutate(comparison=effect) %>% 
#   mutate(pval_symbol=pval_to_symbol(pvalue_effect)) %>% 
#   mutate(p_value_label = paste0("p = ", formatC(pvalue_effect, format = "e", digits = 2))) %>% 
#   mutate(p_value_label=if_else(pvalue_effect>0.05, "n.s.", p_value_label))
  


##cores models-----

#Import Chamberdata4paper.csv
datacores4paper<-read.csv(file = paste0(paper_path,"CoresData4paper.csv"))

#Format main data: 
datacores4paper<- datacores4paper %>% 
  #Remove NAs
  filter(!is.na(dailyflux)) %>% 
  #Factor grouping variables:
  mutate(season=factor(season, levels = c("S1","S2","S3","S4"), ordered = T),
         casepilot=factor(casepilot, levels = c("DU","RI","CA","VA","DA","CU"), ordered = T),
         status=factor(status, levels = c("Preserved","Altered","Restored"), ordered = T),
         subsite=factor(subsite, ordered = F),
         ghgspecies=factor(ghgspecies, ordered=F),
         sampling=factor(sampling, ordered = F))



#Import results from models (R2 and significances)

summary_coresmodels_all<- read.csv(paste0(paper_path, "Allghg_summary_coresmodels.csv"))

significances_cores_alleffects<- summary_coresmodels_all %>% 
  dplyr::select(dataset, status_pval, season_pval, interaction_pval) %>% 
  separate(dataset, into = c("casepilot","ghgspecies")) %>% 
  pivot_longer(c(status_pval,season_pval,interaction_pval),names_sep = "_",values_to = "pvalue_effect", names_to = c("effect","drop")) %>% 
  dplyr::select(-drop)

head(significances_cores_alleffects)

#Import Emmeans post-hoc (all ghgs Combined into single table)
emmeans_cores<-read.csv(paste0(paper_path,"Allghg_emmeans-posthoc_coresmodels.csv"))
head(emmeans_cores)
#Format Emmeans:
emmeans_cores<- emmeans_cores %>% 
  #add effect (status,season, interaction), to join post-hoc comparisons with significances of effects. The significance of interaction is used to filter the status_within_season comparisons
  mutate(effect=if_else(comparison!="status_within_season", comparison, "interaction")) %>% 
  #Join signifcances of model effects
  left_join(significances_cores_alleffects, by=c("casepilot","ghgspecies","effect")) %>% 
  #Replace post-hoc letters with "" to avoid non-significant effect comparisons to be shown
  mutate(group_letter=if_else(pvalue_effect>0.05,"",group_letter)) %>% 
  mutate(ghgspecies=factor(ghgspecies, ordered=F),
         casepilot=factor(casepilot, levels = c("DU","RI","CA","VA","DA","CU"), ordered = T),
         status=factor(status, levels = c("Preserved","Altered","Restored"), ordered = T),
         season=factor(season, levels = c("S1","S2","S3","S4"), ordered = T)) 


#Function to get symbols for pvalues
pval_to_symbol <- function(p) {
  ifelse(p < 0.001, "(***)",
         ifelse(p < 0.01, "(**)",
                ifelse(p < 0.05, "(*)",
                       ifelse(p < 0.1, "(.)", "(ns)")
                )
         )
  )
}


formated_pvalues_all<- significances_cores_alleffects %>% 
  mutate(comparison=effect) %>% 
  mutate(pval_symbol=pval_to_symbol(pvalue_effect)) %>% 
  mutate(p_value_label = paste0("p = ", formatC(pvalue_effect, format = "e", digits = 2))) %>% 
  mutate(p_value_label=if_else(pvalue_effect>0.05, "n.s.", p_value_label))





##Specific models------

#TO ADAPT. 
#Import results from models (R2 and significances)

summary_specificmodels_all<- read.csv(paste0(paper_path, "Allghg_summary_Specificchambermodels.csv"))

significances_specificmodels_alleffects<- summary_specificmodels_all %>% 
  dplyr::select(dataset, status_pval, season_pval, interaction_pval) %>% 
  separate(dataset, into = c("casepilot","ghgspecies")) %>% 
  pivot_longer(c(status_pval,season_pval,interaction_pval),names_sep = "_",values_to = "pvalue_effect", names_to = c("effect","drop")) %>% 
  dplyr::select(-drop)

head(significances_specificmodels_alleffects)

#IMPORTANT-----
#Missing strata significance and 3-way interaction, need adapt function in modelling script to work with specific models (to add strata significance when present in fixed effects)

#Import Emmeans post-hoc (all ghgs Combined into single table)
emmeans_specificmodels_all<-read.csv(paste0(paper_path,"Allghg_emmeans-posthoc_specificchambermodels.csv"))

#Format Emmeans: CANNOT filter letters for effect significance until I fix the summary function above. TO remove groupletters from comparisons when effect was not significant.
unique(emmeans_specificmodels_all$comparison)

emmeans_specificmodels_all<- emmeans_specificmodels_all %>% 
  # mutate(effect=if_else(comparison=="status_within_season","interaction",comparison)) %>% 
  # #Join signifcances of model effects
  # left_join(significances_alleffects , by=c("casepilot","ghgspecies","effect")) %>% 
  # #Replace post-hoc letters with "" to avoid non-significant effect comparisons to be shown
  # #We are leaving the pvalue threshold at 0.1, (va status_within_season for ch4 is almmost significant)
  # mutate(group_letter=if_else(pvalue_effect>0.1,"",group_letter)) %>%
  mutate(ghgspecies=factor(ghgspecies, ordered=F),
         casepilot=factor(casepilot, levels = c("DU","RI","CA","VA","DA","CU"), ordered = T),
         status=factor(status, levels = c("Preserved","Altered","Restored"), ordered = T),
         season=factor(season, levels = c("S1","S2","S3","S4"), ordered = T)) %>% 
  arrange(casepilot, ghgspecies, comparison, season,status)


#separate emmeans per type of comparison
emmeans_status_specificmodels<- emmeans_specificmodels_all %>% 
  filter(comparison=="status") %>% 
  group_by(ghgspecies,casepilot) %>% 
  #If there are no differences (same letter for all status of a given ghgspecies*casepilot combination), remove the group_letter. 
  mutate(group_letter = if (n_distinct(group_letter) == 1) "" else group_letter) %>%
  ungroup()


emmeans_season_specificmodels<- emmeans_specificmodels_all %>% 
  filter(comparison=="season") %>% 
  group_by(ghgspecies,casepilot) %>% 
  #If there are no differences (same letter for all status of a given ghgspecies*casepilot combination), remove the group_letter. 
  mutate(group_letter = if (n_distinct(group_letter) == 1) "" else group_letter) %>%
  ungroup()


#Subset the emmeans and group-letters of status_within_season comparison. 
emmeans_status_within_season_specificmodels<- emmeans_specificmodels_all %>% 
  filter(comparison=="status_within_season") %>% 
  group_by(ghgspecies, casepilot, season) %>%
  #If there are no differences between status within a given combination of ghgspecies*casepilot*season, remove the group_letter. This avoids printing group_letters when no difference was seen. 
  mutate(group_letter = if (n_distinct(group_letter) == 1) "" else group_letter) %>%
  ungroup()


#Subset the emmeans and group-letters of status_within_strata comparison. 
emmeans_status_within_strata_specificmodels<- emmeans_specificmodels_all %>% 
  filter(comparison=="status_within_strata") %>% 
  group_by(ghgspecies, casepilot, strata) %>%
  #If there are no differences between status within a given combination of ghgspecies*casepilot*strata, remove the group_letter. This avoids printing group_letters when no difference was seen. 
  mutate(group_letter = if (n_distinct(group_letter) == 1) "" else group_letter) %>%
  ungroup()









#DONE CREATE:SEASONAL Plots (1per ghg, casepilot per pannel)(with posthoc letter groups), no status.
#DONE status_within_season plots 1per ghg, casepilot per pannel (with symbol per pair status|season): x=season, col=status (FAIL), instead posthoc letter groups



#TO-DO: --------

#CREATE: Plots of specific models with data (1 per casepilot, ghg per pannel).

#Compare gwp estimates for status: via GWP model vs via CO2 and CH4 model and error propagation. 
  #We could compare the obs-vs-pred of the two methods. 


#USE actual strata distribution of valid plotcodes for plots (the same way it is implemented for permanova in modelChamber script). Create: strata-composition plot (all casepilots together, how to add permanova significance?). Export from there to create stratacomp to plot. 

#Import stratacomposition_in-situ 
# stratacomp<- read.csv(file = paste0(paper_path,"Stratacomposition_in-situ.csv"))

#____________________--------
#MAIN PLOTS------


#Fig. 1 (Preserved Baseline)-----
#A 3-panel vertical plot. 1 panel per gas. Boxplot (without outliers). Linear scale. 
#For CH4, inset with DU,RI,CA,VA zoomed in. 

base_co2_ident<-
  data4paper %>% 
  filter(status=="Preserved", ghgspecies%in%c("co2")) %>% 
  ggplot(aes(x=casepilot, y=(dailyflux*1000)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_boxplot(width=0.2,outliers = F, fill="#00BA38",size = 0.7)+
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))+
  guides(color = "none", fill="none")+
  labs(y= expression(CO[2]~NEE~(mmol~d^-1~m^-2)),
       x=NULL,
       fill=paste0("Status"))



# Main plot
base_ch4_ident_main <- data4paper %>%
  filter(status == "Preserved", ghgspecies %in% c("ch4")) %>%
  ggplot(aes(x = casepilot, y = dailyflux)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(width = 0.2, outliers = FALSE, fill = "#00BA38", size = 0.7) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5)
  ) +
  guides(color = "none", fill = "none") +
  labs(
    y = expression(CH[4] ~ NEE ~ (mmol ~ d^-1 ~ m^-2)),
    x = NULL,
    fill = "Status"
  )

# Inset plot
inset_plot <- data4paper %>%
  filter(status == "Preserved", ghgspecies %in% c("ch4"),
         casepilot %in% c("DU", "RI", "CA", "VA")) %>%
  ggplot(aes(x = casepilot, y = dailyflux)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(width = 0.2, outliers = FALSE, fill = "#00BA38", size = 0.7) +
  theme_bw(base_size = 8) +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size = 6, face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)
  )

#Create the grob of your inset plot
inset_grob <- ggplotGrob(inset_plot)

#define position of inset and geom_rect (in main-plot scale)
xmin<- 1
xmax<- 4
ymin<- 5
ymax<- 20

#Extract and store the y-range (without outliers) of the main plot (they are overriden when using annotation_custom
y_range <- range(ggplot_build(base_ch4_ident_main)$layout$panel_params[[1]]$y.range)


#Combine plot with inset (using positions and limits defined above) using annotation_custom
base_ch4_ident_comb <- base_ch4_ident_main +
  annotation_custom(
    grob = inset_grob,
    xmin = xmin, xmax = xmax,
    ymin = ymin, ymax = ymax
  ) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            color = "black", fill = NA) +
  coord_cartesian(ylim = c(y_range[1], y_range[2]))


  
base_gwp_ident<-data4paper %>% 
  filter(status=="Preserved", ghgspecies%in%c("GWPco2andch4")) %>% 
  mutate(casepilot=factor(casepilot, levels = c("DU","RI","CA","VA","DA","CU"))) %>% 
  ggplot(aes(x=casepilot, y=(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_boxplot(width=0.2,outliers = F, fill="#00BA38",size = 0.7)+
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))+
  guides(color = "none", fill="none")+
  labs(y= expression(GWP~NEE~(g*CO[2~eq]~d^-1~m^-2)),
       x=paste0("Casepilot"),
       fill=paste0("Status"))


#COMBINE all 3 panels (with CH4 inset)
baseline_3sp_ident_withinset <- plot_grid(
  base_co2_ident,
  base_ch4_ident_comb,
  base_gwp_ident,
  ncol = 1,
  align = "v",
  axis = "lr"
  # labels = c("A", "B", "C"),
  # label_size = 12,
  # label_x = 0.02,  # move label closer to left
  # label_y = 1      # top-aligned
)

#Save:
ggsave(plot = baseline_3sp_ident_withinset, 
       filename = "PAPERPLOTS_Fig1_baseline_3ghg_boxplot_withinset.png",
       path = main_figures,
       device = "png",
       dpi = 400,
       units = "mm",
       height = 200, 
       width = 90)


#Extra alternative: NO INSET for CH4
#COMBINE all 3 panels (without CH4 inset)
baseline_3sp_ident_noinset <- plot_grid(
  base_co2_ident,
  base_ch4_ident_main,
  base_gwp_ident,
  ncol = 1,
  align = "v",
  axis = "lr"
  # labels = c("A", "B", "C"),
  # label_size = 12,
  # label_x = 0.02,  # move label closer to left
  # label_y = 1      # top-aligned
)

#Save:
ggsave(plot = baseline_3sp_ident_noinset, 
       filename = "PAPERPLOTS_Fig1_baseline_3ghg_boxplot_noinset.png",
       path = extra_figures,
       device = "png",
       dpi = 400,
       units = "mm",
       height = 200, 
       width = 90)




#Template Plot Status------

#Preserved-Altered-Restored fluxes (boxplots, no "outliers") and group differences (EMMeans post-hoc), for each casepilot (vertical pannels). 

#We create a single plot template object (for which we will change the data-origin and yaxis label for each ghgspecies). It takes 4 arguments: GHG ("co2", "ch4","GWPco2andch4"),  y_label (the expression call to use for the y axis), y_multiplier: (multiplier of dailyflux to use, 1000 for CO2 in mmol, 1 for CH4mmol, 1 for GWP in gCO2eq), drawboxplot_outliers (T/F)

plot_ghg_faceted <- function(GHG,
                             y_label, 
                             y_multiplier = 1,
                             drawboxplot_outliers=F) {
  
  #filter data for GHG
  data<- data4paper %>% filter(ghgspecies==GHG)
  emmeans_data<- emmeans_status %>% filter(ghgspecies==GHG) 
  
  #Produce plot: 
  ggplot(data, aes(x = status, y = dailyflux * y_multiplier)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_boxplot(width = 0.2, outliers = drawboxplot_outliers, aes(fill = status), size = 0.7) +
    
    # Add emmeans points
    geom_point(data = emmeans_data, aes(x = status, y = emmean_bt * y_multiplier),
               shape = 23, size = 3, fill = "black") +
    
    # Add group letters
    geom_text(data = emmeans_data, aes(x = status, label = group_letter, y = Inf),
              vjust = 1.2, hjust = 0.5, color = "black", inherit.aes = FALSE,
              size = 5, fontface = "bold") +
    theme_bw() +
    scale_fill_manual(values = c(
      "Preserved" = "#009E73",
      "Altered"   = "#D55E00",
      "Restored"  = "#56B4E9"
    )) +
    scale_color_manual(values = c(
      "Preserved" = "#009E73",
      "Altered"   = "#D55E00",
      "Restored"  = "#56B4E9"
    )) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.2), 0)) +
    theme(
      axis.text = element_text(face = "bold")
    ) +
    guides(color = "none", fill = "none") +
    labs(
      y = y_label,
      x = "Conservation status",
      fill = "Status"
    ) +
    facet_grid(rows = vars(casepilot), scales = "free")
}


#Fig. 2 (CO2 status)-----

#Produce plot using status template and subset data for CO2: 
plot_ghg_faceted(GHG="co2",
                 y_label = expression(CO[2]~NEE~(mmol~d^-1~m^-2)),
                 y_multiplier = 1000)


#Save plot: 
ggsave(filename = "PAPERPLOTS_Fig2_status_CO2_molar_boxplot_allCP.png",
       path = main_figures,
       device = "png",
       dpi = 400,
       units = "mm",
       height = 229,
       width = 90)



#Fig. 3 (CH4 status)-----

#Produce plot using status template and subset data for CH4: 
plot_ghg_faceted(GHG="ch4",
                 y_label = expression(CH[4]~NEE~(mmol~d^-1~m^-2)),
                 y_multiplier = 1)

#Save plot: 
ggsave(filename = "PAPERPLOTS_Fig3_status_CH4_molar_boxplot_allCP.png",
       path = main_figures,
       device = "png",
       dpi = 400,
       units = "mm",
       height = 229,
       width = 90)


#Fig. 4 (GWP status)-----

#Produce plot using status template and subset data for GWP: 
plot_ghg_faceted(GHG="GWPco2andch4",
                 y_label = expression(GWP~NEE~(g*CO[2~eq]~d^-1~m^-2)),
                 y_multiplier = 1)

#Save plot: 
ggsave(filename = "PAPERPLOTS_Fig4_status_GWP_boxplot_allCP.png",
       path = main_figures,
       device = "png",
       dpi = 400,
       units = "mm",
       height = 229,
       width = 90)




#TEST plots cores------

#Models dont behave as expected, I believe subsite-specific variability is masking effect of status in some cases, many other cases core-flux variability is too high to see any effect (also linked to high geographical variability, how representative is a core of the whole subsite?)


#Produce plot using status template and subset data for CO2: 
plot_ghg_faceted(data = datacores4paper %>% filter(ghgspecies=="co2"),
                 emmeans_data = emmeans_cores %>% filter(comparison=="status",ghgspecies=="co2"),
                 y_label = expression(CO[2]~NEE~(mmol~d^-1~m^-2)),
                 y_multiplier = 1,
                 drawboxplot_outliers= F)


#Produce plot using status template and subset data for CH4: 
plot_ghg_faceted(data = datacores4paper %>% filter(ghgspecies=="ch4"),
                 emmeans_data = emmeans_cores %>% filter(comparison=="status",ghgspecies=="ch4"),
                 y_label = "ch4",
                 y_multiplier = 1,
                 drawboxplot_outliers= F)

plot_ghg_faceted(data = datacores4paper %>% filter(ghgspecies=="n2o"),
                 emmeans_data = emmeans_cores %>% filter(comparison=="status",ghgspecies=="n2o"),
                 y_label = "n2o",
                 y_multiplier = 1,
                 drawboxplot_outliers= F)


plot_ghg_faceted(data = datacores4paper %>% filter(ghgspecies=="GWPco2andch4"),
                 emmeans_data = emmeans_cores %>% filter(comparison=="status",ghgspecies=="GWPco2andch4"),
                 y_label = "GWPco2andch4",
                 y_multiplier = 1,
                 drawboxplot_outliers= F)

plot_ghg_faceted(data = datacores4paper %>% filter(ghgspecies=="GWPtotal"),
                 emmeans_data = emmeans_cores %>% filter(comparison=="status",ghgspecies=="GWPtotal"),
                 y_label = "GWPtotal",
                 y_multiplier = 1,
                 drawboxplot_outliers= F)



#____________________--------
#TABLES--------

#Significance of effects------
#Format model summary into table with effect significance, R2c and R2m
#Arrange by co2,ch4,GWP, then casepilot (DU,RI,CA,VA,DA,CU) , then effect (status,season,interaction)


formated_summary_models<- summary_models_all %>% 
  dplyr::select(dataset, status_pval, season_pval, interaction_pval, homoced_R2m, homoced_R2c) %>% 
  pivot_longer(cols=c(status_pval, season_pval, interaction_pval), names_to = c("effect","drop"),names_sep = "_", values_to = "p_value") %>% 
  separate(dataset, into = c("casepilot", "ghgspecies")) %>% 
  mutate(casepilot=factor(casepilot, levels = c("DU","RI","CA","VA","DA","CU"), ordered = T),
         ghgspecies=factor(ghgspecies, levels = c("co2","ch4","GWPco2andch4"), ordered = T),
         effect=factor(effect, levels = c("status","season","interaction"), ordered = T),
         sig_symbol=gsub("\\)","",gsub("\\(","",pval_to_symbol(p_value)))) %>% 
  rename(R2c=homoced_R2c, R2m=homoced_R2m) %>% 
  dplyr::select(ghgspecies, casepilot, effect, sig_symbol, R2m, R2c, p_value) %>% 
  mutate(R2m=round(R2m,digits = 3), R2c=round(R2c, digits = 3)) %>% 
  arrange(ghgspecies, casepilot, effect)

write.csv(formated_summary_models, file=paste0(main_figures, "EffectSignif_and_Rsquared_chambermodels.csv"), row.names = F)



#____________________--------
#SUP. PLOTS ------



#Season------

#Per-season boxplots (no "outliers") and group differences (EMMeans post-hoc), for each casepilot (vertical pannels). 

#We create a single plot template object (for which we will change the data-origin and yaxis label for each ghgspecies). It takes 4 arguments: data (the ghg-specific in-situ dataframe), emmeans_data (the ghg-specific emmeans dataframe, comparison==status), y_label (the expression call to use for the y axis), y_multiplier: (multiplier of dailyflux to use, 1000 for CO2 in mmol, 1 for CH4mmol, 1 for GWP in gCO2eq)

season_plot <- function(GHG, y_label, y_multiplier = 1, drawboxplot_outliers=F) {
  #Filter and format
  data<- data4paper %>% filter(ghgspecies==GHG) %>% 
    mutate(season=case_when(season=="S1"~"Autumn",
                            season=="S2"~"Winter",
                            season=="S3"~"Spring",
                            season=="S4"~"Summer")) %>% 
    mutate(season=factor(season, levels=c("Autumn", "Winter","Spring","Summer"), ordered = T))
  #Filter and format
  emmeans_data<- emmeans_season %>% filter(ghgspecies==GHG)%>% 
    mutate(season=case_when(season=="S1"~"Autumn",
                            season=="S2"~"Winter",
                            season=="S3"~"Spring",
                            season=="S4"~"Summer")) %>% 
    mutate(season=factor(season, levels=c("Autumn", "Winter","Spring","Summer"), ordered = T))
  
  #Plot
  ggplot(data, aes(x = season, y = dailyflux * y_multiplier)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_boxplot(width = 0.2, outliers = drawboxplot_outliers, aes(fill = season), size = 0.7) +
    
    # Add emmeans points
    geom_point(data = emmeans_data, aes(x = season, y = emmean_bt * y_multiplier),
               shape = 23, size = 3, fill = "black") +
    
    # Add group letters
    geom_text(data = emmeans_data, aes(x = season, label = group_letter, y = Inf),
              vjust = 1.2, hjust = 0.5, color = "black", inherit.aes = FALSE,
              size = 5, fontface = "bold") +
    scale_fill_manual(values = c(
      "Autumn" = "#FF6961",# Deep red for Autumn
      "Winter"   = "#AEC6CF",# Cool blue for Winter
      "Spring"  = "#77DD77",# Light green for Spring
      "Summer" ="#FFB347"# Warm orange for Summer
    )) +
    theme_bw() +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.2), 0)) +
    theme(
      axis.text = element_text(face = "bold")
    ) +
    guides(color = "none", fill = "none") +
    labs(
      y = y_label,
      x = "Season",
      fill= "Season"
    ) +
    facet_grid(rows = vars(casepilot), scales = "free")
}


## Save plots-----
#CO2 
season_plot(GHG="co2",
            y_label = expression(CO[2]~NEE~(mmol~d^-1~m^-2)),
            y_multiplier = 1000,
            drawboxplot_outliers = F)


#Save plot: 
ggsave(filename = "Sup_Fig_season_CO2.png",
       path = sup_figures,
       device = "png",
       dpi = 400,
       units = "mm",
       height = 229,
       width = 100)

#CH4
season_plot(GHG="ch4",
            y_label = expression(CH[4]~NEE~(mmol~d^-1~m^-2)),
            y_multiplier = 1,
            drawboxplot_outliers = F)

#Save plot: 
ggsave(filename = "Sup_Fig_season_CH4.png",
       path = sup_figures,
       device = "png",
       dpi = 400,
       units = "mm",
       height = 229,
       width = 100)

#GWPco2andch4
season_plot(GHG="GWPco2andch4",
            y_label = expression(GWP~NEE~(g*CO[2~eq]~d^-1~m^-2)),
            y_multiplier = 1,
            drawboxplot_outliers = F)

#Save plot: 
ggsave(filename = "Sup_Fig_season_GWP.png",
       path = sup_figures,
       device = "png",
       dpi = 400,
       units = "mm",
       height = 229,
       width = 100)

#Status per season ------

#Plot both data and emmeans across status*season for each casepilot. 
#Color for status, x-position for season. 
#Points for emmeans and group-letters for status_within_season comparisons


##with groupletters (ok)-----

#Works finee. Difficult to adjust y-axis of letters dynamically (now set to y=Inf)
plot_status_per_season <- function(GHG, y_label, 
                                        y_multiplier = 1, 
                                       dodge_width = 0.4) {
  #Filter datasets for GHG:
  data<- data4paper %>% filter(ghgspecies==GHG)%>% 
    mutate(season=case_when(season=="S1"~"Autumn",
                            season=="S2"~"Winter",
                            season=="S3"~"Spring",
                            season=="S4"~"Summer")) %>% 
    mutate(season=factor(season, levels=c("Autumn", "Winter","Spring","Summer"), ordered = T))
  
  emmeans_data<- emmeans_status_within_season %>% filter(ghgspecies==GHG)%>% 
    mutate(season=case_when(season=="S1"~"Autumn",
                            season=="S2"~"Winter",
                            season=="S3"~"Spring",
                            season=="S4"~"Summer")) %>% 
    mutate(season=factor(season, levels=c("Autumn", "Winter","Spring","Summer"), ordered = T))
  
  #Produce plots:
  ggplot(data, aes(x = season, y = dailyflux * y_multiplier)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    # Boxplots with dodge
    geom_boxplot(aes(fill = status), width = 0.2, outliers = F,
                 position = position_dodge(width = dodge_width), size = 0.7) +
    # Emmeans points with same dodge
    geom_point(data = emmeans_data,
               aes(x = season, y = emmean_bt * y_multiplier, group = status),
               shape = 23, size = 2, fill = "black",
               position = position_dodge(width = dodge_width)) +
    # Add group letters with same dodge
    geom_text(data = emmeans_data, 
              aes(x = season, label = group_letter, y = Inf, group=status),
              vjust = 1.2, hjust = 0.5, color = "black", inherit.aes = FALSE,
              size = 4.5, fontface = "bold",
              position = position_dodge(width = dodge_width)) +
    #expand y-scale to avoid overlap between group_letters and wiskers
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.2), 0)) +
    theme_bw() +
    scale_fill_manual(values = c(
      "Preserved" = "#009E73",
      "Altered"   = "#D55E00",
      "Restored"  = "#56B4E9"
    )) +
    theme(
      axis.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 0, hjust = 0.5,vjust=0.5)
    ) +
    guides(color = "none") +
    labs(
      y = y_label,
      x = "Season",
      fill = "Status"
    ) +
    facet_grid(rows = vars(casepilot), scales = "free") 
  
}


## Save plots-----

#CO2: 
plot_status_per_season(GHG = "co2",
                       y_label = expression(CO[2]~NEE~(mmol~d^-1~m^-2)),
                       y_multiplier = 1000, #to change from molar to milli-molar units
                       dodge_width = 0.6)

#Save plot: 
ggsave(filename = "Sup_Fig_status_per_season_CO2.png",
       path = sup_figures,
       device = "png",
       dpi = 400,
       units = "mm",
       height = 229,
       width = 180)



#CH4
plot_status_per_season(GHG = "ch4",
                       y_label = expression(CH[4]~NEE~(mmol~d^-1~m^-2)),
                       y_multiplier = 1,
                       dodge_width = 0.6)

#Save plot: 
ggsave(filename = "Sup_Fig_status_per_season_CH4.png",
       path = sup_figures,
       device = "png",
       dpi = 400,
       units = "mm",
       height = 229,
       width = 180)


#GWPco2andch4
plot_status_per_season(GHG = "GWPco2andch4",
                       y_label = expression(GWP~NEE~(g*CO[2~eq]~d^-1~m^-2)),
                       y_multiplier = 1,
                       dodge_width = 0.6)

#Save plot: 
ggsave(filename = "Sup_Fig_status_per_season_GWP.png",
       path = sup_figures,
       device = "png",
       dpi = 400,
       units = "mm",
       height = 229,
       width = 180)



##(omit)Test sig. and brackets----

#Very difficult and complex to do (skip for the moment)
#TO Add significance symbols, we need to subset the emmeans_contrasts and calculate the appropriate y-axis positions. WE will set the preseved vs altered and altered vs restored to X% above the max value of the data (filtered for 1.5*IQR, same as boxplot(outliers=F) ), for preserved vs restored we will set the significance to be above the previous two. 

#Calculate max of wiskers position for each pannel
status_within_season_panel_max<- data4paper %>% 
  group_by(status, season, casepilot, ghgspecies) %>% 
  summarise(
    Q1 = quantile(dailyflux, 0.25, na.rm = TRUE),
    Q3 = quantile(dailyflux, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    #Calculate upper_whisker position of boxplots
    upper_whisker = min(max(dailyflux, na.rm = TRUE), Q3 + 1.5 * IQR),
    .groups = "drop"
  ) %>%
  group_by(casepilot, season,ghgspecies) %>%
  summarise(y_max = max(upper_whisker), .groups = "drop")


#Filter sig.symbols and set y-position for each pannel
status_within_season_signif<-   emmeans_contrasts %>%
  filter(comparison == "status_within_season") %>%
  separate(contrast, into = c("group1", "group2"), sep = " - ", remove = T) %>% 
  dplyr::select(ghgspecies, casepilot, season, group1, group2, p.value,pval_symbol) %>% 
  left_join(status_within_season_panel_max, by = c("casepilot", "season","ghgspecies")) %>%
  mutate(
    y.position = y_max * 1.05,  # 5% above max value
    label = pval_symbol
  ) %>% 
  mutate(casepilot=factor(casepilot, levels = c("DU","RI","CA","VA","DA","CU"), ordered = T))


#Function to plot: 
plot_status_by_season <- function(data, emmeans_data, signif_df, y_label, y_multiplier = 1) {
  ggplot(data, aes(x = status, y = dailyflux * y_multiplier)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    
    # Boxplots
    geom_boxplot(aes(fill = status), width = 0.6, outliers = F, size = 0.7) +
    
    # Emmeans points
    geom_point(data = emmeans_data,
               aes(x = status, y = emmean_bt * y_multiplier),
               shape = 23, size = 3, fill = "black") +
    
    # Significance annotations
    stat_pvalue_manual(
      data = signif_df,
      label = "pval_symbol",
      xmin = "group1",
      xmax = "group2",
      y.position = "y.position",
      tip.length = 0.01,
      bracket.size = 0.5,
      size = 5
    ) +
    
    facet_grid(cols = vars(season), rows = vars(casepilot), scales = "free") +
    
    scale_fill_manual(values = c(
      "Preserved" = "#009E73",
      "Altered"   = "#D55E00",
      "Restored"  = "#56B4E9"
    )) +
    
    theme_bw() +
    theme(
      axis.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)
    ) +
    
    guides(color = "none") +
    labs(
      y = y_label,
      x = "Status",
      fill = "Status"
    )
}


#example: 
plot_status_by_season(data=data4paper %>% filter(ghgspecies=="ch4"),
                      emmeans_data = emmeans_all %>% filter(comparison=="status_within_season",ghgspecies=="ch4"),
                      signif_df = status_within_season_signif %>% filter(ghgspecies=="ch4"),
                      y_label = "ch4",
                      y_multiplier = 1
)

#ISSUES: need to set offset for y-position for each of the 3 comparisons, additionally, pannel limits clip significance symbols. Potentially, set a single y-position for each casepilot (to be reused across seasons) and offset the 3 comparisons proportionally to this y-position. This will result in same-height of sig.symbols across seasons, set by the maximum wisker height across all boxplots of the same casepilot. 


#____________________--------
#EXTRA PLOTS-------

#Specific model plots----
##DA 4-level status-----

#Fitst need to re-level status in data4paper 
da_specific_data <- data4paper %>% 
  filter(casepilot=="DA") %>% 
  mutate(status=case_when(subsite=="DA-A1"~"Altered 1",
                          subsite=="DA-A2"~"Altered 2",
                          status=="Preserved"~"Preserved",
                          status=="Restored"~"Restored")) %>% 
  mutate(status=factor(status, levels = c("Preserved","Altered 1","Altered 2", "Restored"), ordered = T))

head(da_specific_data)
str(da_specific_data)

da_emmeans_status<- emmeans_status_specificmodels %>%
  filter(casepilot=="DA") %>% 
  mutate(status=factor(status, levels = c("Preserved","Altered 1","Altered 2", "Restored"), ordered = T))
  
str(da_emmeans_status)

#CO2 plot: 
da_co2_specific_plot<-
  ggplot(da_specific_data %>% filter(ghgspecies=="co2"), 
         aes(x = status, y = dailyflux * 1000)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(width = 0.2, outliers = F, aes(fill = status), size = 0.7) +
  
  # Add emmeans points
  geom_point(data = da_emmeans_status %>% filter(ghgspecies=="co2"),
             aes(x = status, y = emmean_bt * 1000),
             shape = 23, size = 3, fill = "black") +
  
  # Add group letters
  geom_text(data = da_emmeans_status %>% filter(ghgspecies=="co2"),
            aes(x = status, label = group_letter, y = Inf),
            vjust = 1.2, hjust = 0.5, color = "black", inherit.aes = FALSE,
            size = 5, fontface = "bold") +
  theme_bw() +
  scale_fill_manual(values = c(
    "Preserved" = "#009E73",
    "Altered 1" = "#D55E00",
    "Altered 2" = "orange",
    "Restored"  = "#56B4E9"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2), 0)) +
  theme(
    axis.text = element_text(face = "bold")
  ) +
  guides(color = "none", fill = "none") +
  labs(
    y = expression(CO[2]~NEE~(mmol~d^-1~m^-2)),
    x = NULL,
    fill = "Status"
  )

da_co2_specific_plot



#CH4 plot: 
da_ch4_specific_plot<-
  ggplot(da_specific_data %>% filter(ghgspecies=="ch4"), 
         aes(x = status, y = dailyflux)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(width = 0.2, outliers = F, aes(fill = status), size = 0.7) +
  
  # Add emmeans points
  geom_point(data = da_emmeans_status %>% filter(ghgspecies=="ch4"),
             aes(x = status, y = emmean_bt),
             shape = 23, size = 3, fill = "black") +
  
  # Add group letters
  geom_text(data = da_emmeans_status %>% filter(ghgspecies=="ch4"),
            aes(x = status, label = group_letter, y = Inf),
            vjust = 1.2, hjust = 0.5, color = "black", inherit.aes = FALSE,
            size = 5, fontface = "bold") +
  theme_bw() +
  scale_fill_manual(values = c(
    "Preserved" = "#009E73",
    "Altered 1" = "#D55E00",
    "Altered 2" = "orange",
    "Restored"  = "#56B4E9"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2), 0)) +
  theme(
    axis.text = element_text(face = "bold")
  ) +
  guides(color = "none", fill = "none") +
  labs(
    y = expression(CH[4]~NEE~(mmol~d^-1~m^-2)),
    x = NULL,
    fill = "Status"
  )

da_ch4_specific_plot


#GWP plot: 
da_gwp_specific_plot<-
  ggplot(da_specific_data %>% filter(ghgspecies=="GWPco2andch4"), 
         aes(x = status, y = dailyflux)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(width = 0.2, outliers = F, aes(fill = status), size = 0.7) +
  # Add emmeans points
  geom_point(data = da_emmeans_status %>% filter(ghgspecies=="GWPco2andch4"),
             aes(x = status, y = emmean_bt),
             shape = 23, size = 3, fill = "black") +
  # Add group letters
  geom_text(data = da_emmeans_status %>% filter(ghgspecies=="GWPco2andch4"),
            aes(x = status, label = group_letter, y = Inf),
            vjust = 1.2, hjust = 0.5, color = "black", inherit.aes = FALSE,
            size = 5, fontface = "bold") +
  theme_bw() +
  scale_fill_manual(values = c(
    "Preserved" = "#009E73",
    "Altered 1" = "#D55E00",
    "Altered 2" = "orange",
    "Restored"  = "#56B4E9"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2), 0)) +
  theme(
    axis.text = element_text(face = "bold")
  ) +
  guides(color = "none", fill = "none") +
  labs(
    y = expression(GWP~NEE~(g*CO[2~eq]~d^-1~m^-2)),
    x = "Conservation status",
    fill = "Status"
  )

da_gwp_specific_plot

#COMBINE pannels for each GHG
da_specific_plot_3ghg <- plot_grid(
  da_co2_specific_plot,
  da_ch4_specific_plot,
  da_gwp_specific_plot,
  ncol = 1,
  align = "v",
  axis = "lr"
  # labels = c("A", "B", "C"),
  # label_size = 12,
  # label_x = 0.02,  # move label closer to left
  # label_y = 1      # top-aligned
)

da_specific_plot_3ghg

#Save:
ggsave(plot = da_specific_plot_3ghg, 
       filename = "SpecificModels_DA_3ghg_4levelstatus.png",
       path = sup_figures,
       device = "png",
       dpi = 400,
       units = "mm",
       height = 200, 
       width = 90)






##CU considering strata------



#Fitst need to re-level status in data4paper 
cu_specific_data <- data4paper %>% 
  filter(casepilot=="CU") %>% 
  mutate(status=factor(status, levels = c("Preserved","Altered", "Restored"), ordered = T),
         strata=factor(strata))

head(cu_specific_data)
str(cu_specific_data)

cu_emmeans_status_within_strata<- emmeans_status_within_strata_specificmodels %>%
  filter(casepilot=="CU")%>% 
  mutate(status=factor(status, levels = c("Preserved","Altered", "Restored"), ordered = T),
         strata=factor(strata))

str(cu_emmeans_status_within_strata)


#CO2 plot: 
cu_co2_specific_plot<-
  ggplot(cu_specific_data %>% filter(ghgspecies=="co2"), 
         aes(x = status, y = dailyflux * 1000)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(width = 0.2, outliers = F, aes(fill = status), size = 0.7) +
  
  # Add emmeans points
  geom_point(data = cu_emmeans_status_within_strata %>% filter(ghgspecies=="co2"),
             aes(x = status, y = emmean_bt * 1000),
             shape = 23, size = 3, fill = "black") +
  
  # Add group letters
  geom_text(data = cu_emmeans_status_within_strata %>% filter(ghgspecies=="co2"),
            aes(x = status, label = group_letter, y = Inf),
            vjust = 1.2, hjust = 0.5, color = "black", inherit.aes = FALSE,
            size = 5, fontface = "bold") +
  theme_bw() +
  scale_fill_manual(values = c(
    "Preserved" = "#009E73",
    "Altered" = "#D55E00",
    "Restored"  = "#56B4E9"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2), 0)) +
  theme(
    axis.text = element_text(face = "bold")
  ) +
  guides(color = "none", fill = "none") +
  labs(
    y = expression(CO[2]~NEE~(mmol~d^-1~m^-2)),
    x = NULL,
    fill = "Status"
  )+
  facet_wrap(facets=vars(strata), scales="free")

cu_co2_specific_plot



#CH4 plot: 
cu_ch4_specific_plot<-
  ggplot(cu_specific_data %>% filter(ghgspecies=="ch4"), 
         aes(x = status, y = dailyflux)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(width = 0.2, outliers = F, aes(fill = status), size = 0.7) +
  
  # Add emmeans points
  geom_point(data = cu_emmeans_status_within_strata %>% filter(ghgspecies=="ch4"),
             aes(x = status, y = emmean_bt),
             shape = 23, size = 3, fill = "black") +
  
  # Add group letters
  geom_text(data = cu_emmeans_status_within_strata %>% filter(ghgspecies=="ch4"),
            aes(x = status, label = group_letter, y = Inf),
            vjust = 1.2, hjust = 0.5, color = "black", inherit.aes = FALSE,
            size = 5, fontface = "bold") +
  theme_bw() +
  scale_fill_manual(values = c(
    "Preserved" = "#009E73",
    "Altered" = "#D55E00",
    "Restored"  = "#56B4E9"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2), 0)) +
  theme(
    axis.text = element_text(face = "bold")
  ) +
  guides(color = "none", fill = "none") +
  labs(
    y = expression(CH[4]~NEE~(mmol~d^-1~m^-2)),
    x = NULL,
    fill = "Status"
  )+
  facet_wrap(facets=vars(strata), scales="free")

cu_ch4_specific_plot


#GWP plot: 
cu_gwp_specific_plot<-
  ggplot(cu_specific_data %>% filter(ghgspecies=="GWPco2andch4"), 
         aes(x = status, y = dailyflux)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(width = 0.2, outliers = F, aes(fill = status), size = 0.7) +
  # Add emmeans points
  geom_point(data = cu_emmeans_status_within_strata %>% filter(ghgspecies=="GWPco2andch4"),
             aes(x = status, y = emmean_bt),
             shape = 23, size = 3, fill = "black") +
  # Add group letters
  geom_text(data = cu_emmeans_status_within_strata %>% filter(ghgspecies=="GWPco2andch4"),
            aes(x = status, label = group_letter, y = Inf),
            vjust = 1.2, hjust = 0.5, color = "black", inherit.aes = FALSE,
            size = 5, fontface = "bold") +
  theme_bw() +
  scale_fill_manual(values = c(
    "Preserved" = "#009E73",
    "Altered" = "#D55E00",
    "Restored"  = "#56B4E9"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2), 0)) +
  theme(
    axis.text = element_text(face = "bold")
  ) +
  guides(color = "none", fill = "none") +
  labs(
    y = expression(GWP~NEE~(g*CO[2~eq]~d^-1~m^-2)),
    x = "Conservation status",
    fill = "Status"
  )+
  facet_wrap(facets=vars(strata), scales="free")

cu_gwp_specific_plot

#COMBINE pannels for each GHG
cu_specific_plot_3ghg <- plot_grid(
  cu_co2_specific_plot,
  cu_ch4_specific_plot,
  cu_gwp_specific_plot,
  ncol = 1,
  align = "v",
  axis = "lr"
  # labels = c("A", "B", "C"),
  # label_size = 12,
  # label_x = 0.02,  # move label closer to left
  # label_y = 1      # top-aligned
)

cu_specific_plot_3ghg

#Save:
ggsave(plot = cu_specific_plot_3ghg, 
       filename = "SpecificModels_CU_3ghg_2levelstrata.png",
       path = sup_figures,
       device = "png",
       dpi = 400,
       units = "mm",
       height = 200, 
       width = 150)


#OVERALL status effect for CU with updated model: 

##CU_specific_statusonly----

cu_emmeans_statusonly<- emmeans_status_specificmodels %>% 
  filter(casepilot=="CU") %>% 
  mutate(status=factor(status, levels = c("Preserved","Altered", "Restored"), ordered = T))
  

#CO2 plot: 
cu_co2_status_plot<-
  ggplot(cu_specific_data %>% filter(ghgspecies=="co2"), 
         aes(x = status, y = dailyflux * 1000)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(width = 0.2, outliers = F, aes(fill = status), size = 0.7) +
  
  # Add emmeans points
  geom_point(data = cu_emmeans_statusonly %>% filter(ghgspecies=="co2"),
             aes(x = status, y = emmean_bt * 1000),
             shape = 23, size = 3, fill = "black") +
  
  # Add group letters
  geom_text(data = cu_emmeans_statusonly %>% filter(ghgspecies=="co2"),
            aes(x = status, label = group_letter, y = Inf),
            vjust = 1.2, hjust = 0.5, color = "black", inherit.aes = FALSE,
            size = 5, fontface = "bold") +
  theme_bw() +
  scale_fill_manual(values = c(
    "Preserved" = "#009E73",
    "Altered" = "#D55E00",
    "Restored"  = "#56B4E9"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2), 0)) +
  theme(
    axis.text = element_text(face = "bold")
  ) +
  guides(color = "none", fill = "none") +
  labs(
    y = expression(CO[2]~NEE~(mmol~d^-1~m^-2)),
    x = NULL,
    fill = "Status"
  )

cu_co2_status_plot



#CH4 plot: 
cu_ch4_status_plot<-
  ggplot(cu_specific_data %>% filter(ghgspecies=="ch4"), 
         aes(x = status, y = dailyflux)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(width = 0.2, outliers = F, aes(fill = status), size = 0.7) +
  
  # Add emmeans points
  geom_point(data = cu_emmeans_statusonly %>% filter(ghgspecies=="ch4"),
             aes(x = status, y = emmean_bt),
             shape = 23, size = 3, fill = "black") +
  
  # Add group letters
  geom_text(data = cu_emmeans_statusonly %>% filter(ghgspecies=="ch4"),
            aes(x = status, label = group_letter, y = Inf),
            vjust = 1.2, hjust = 0.5, color = "black", inherit.aes = FALSE,
            size = 5, fontface = "bold") +
  theme_bw() +
  scale_fill_manual(values = c(
    "Preserved" = "#009E73",
    "Altered" = "#D55E00",
    "Restored"  = "#56B4E9"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2), 0)) +
  theme(
    axis.text = element_text(face = "bold")
  ) +
  guides(color = "none", fill = "none") +
  labs(
    y = expression(CH[4]~NEE~(mmol~d^-1~m^-2)),
    x = NULL,
    fill = "Status"
  )

cu_ch4_status_plot


#GWP plot: 
cu_gwp_status_plot<-
  ggplot(cu_specific_data %>% filter(ghgspecies=="GWPco2andch4"), 
         aes(x = status, y = dailyflux)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(width = 0.2, outliers = F, aes(fill = status), size = 0.7) +
  # Add emmeans points
  geom_point(data = cu_emmeans_statusonly %>% filter(ghgspecies=="GWPco2andch4"),
             aes(x = status, y = emmean_bt),
             shape = 23, size = 3, fill = "black") +
  # Add group letters
  geom_text(data = cu_emmeans_statusonly %>% filter(ghgspecies=="GWPco2andch4"),
            aes(x = status, label = group_letter, y = Inf),
            vjust = 1.2, hjust = 0.5, color = "black", inherit.aes = FALSE,
            size = 5, fontface = "bold") +
  theme_bw() +
  scale_fill_manual(values = c(
    "Preserved" = "#009E73",
    "Altered" = "#D55E00",
    "Restored"  = "#56B4E9"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2), 0)) +
  theme(
    axis.text = element_text(face = "bold")
  ) +
  guides(color = "none", fill = "none") +
  labs(
    y = expression(GWP~NEE~(g*CO[2~eq]~d^-1~m^-2)),
    x = "Conservation status",
    fill = "Status"
  )

cu_gwp_status_plot

#COMBINE pannels for each GHG
cu_specific_statusonlyplot_3ghg <- plot_grid(
  cu_co2_status_plot,
  cu_ch4_status_plot,
  cu_gwp_status_plot,
  ncol = 1,
  align = "v",
  axis = "lr"
  # labels = c("A", "B", "C"),
  # label_size = 12,
  # label_x = 0.02,  # move label closer to left
  # label_y = 1      # top-aligned
)

cu_specific_statusonlyplot_3ghg

ggsave(plot = cu_specific_statusonlyplot_3ghg, 
       filename = "SpecificModels_CU_3ghg_generalstatus.png",
       path = sup_figures,
       device = "png",
       dpi = 400,
       units = "mm",
       height = 200, 
       width = 90)




##VA considering strata------

#Fitst need to re-level strata  
va_specific_data <- data4paper %>% 
  filter(casepilot=="VA") %>% 
  mutate(status=factor(status, levels = c("Preserved","Altered", "Restored"), ordered = T),
         strata=as.factor(if_else(strata=="vegetated","vegetated", "non-vegetated")))
         

head(va_specific_data)
str(va_specific_data)

va_emmeans_status_within_strata<- emmeans_status_within_strata_specificmodels %>%
  filter(casepilot=="VA")%>% 
  mutate(status=factor(status, levels = c("Preserved","Altered", "Restored"), ordered = T),
         strata=as.factor(if_else(strata=="vegetated","vegetated", "non-vegetated")))

str(va_emmeans_status_within_strata)


#CO2 plot: 
va_co2_specific_plot<-
  ggplot(va_specific_data %>% filter(ghgspecies=="co2"), 
         aes(x = status, y = dailyflux * 1000)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(width = 0.2, outliers = F, aes(fill = status), size = 0.7) +
  
  # Add emmeans points
  geom_point(data = va_emmeans_status_within_strata %>% filter(ghgspecies=="co2"),
             aes(x = status, y = emmean_bt * 1000),
             shape = 23, size = 3, fill = "black") +
  
  # Add group letters
  geom_text(data = va_emmeans_status_within_strata %>% filter(ghgspecies=="co2"),
            aes(x = status, label = group_letter, y = Inf),
            vjust = 1.2, hjust = 0.5, color = "black", inherit.aes = FALSE,
            size = 5, fontface = "bold") +
  theme_bw() +
  scale_fill_manual(values = c(
    "Preserved" = "#009E73",
    "Altered" = "#D55E00",
    "Restored"  = "#56B4E9"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2), 0)) +
  theme(
    axis.text = element_text(face = "bold")
  ) +
  guides(color = "none", fill = "none") +
  labs(
    y = expression(CO[2]~NEE~(mmol~d^-1~m^-2)),
    x = NULL,
    fill = "Status"
  )+
  facet_wrap(facets=vars(strata), scales="free")

va_co2_specific_plot



#CH4 plot: 
va_ch4_specific_plot<-
  ggplot(va_specific_data %>% filter(ghgspecies=="ch4"), 
         aes(x = status, y = dailyflux)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(width = 0.2, outliers = F, aes(fill = status), size = 0.7) +
  
  # Add emmeans points
  geom_point(data = va_emmeans_status_within_strata %>% filter(ghgspecies=="ch4"),
             aes(x = status, y = emmean_bt),
             shape = 23, size = 3, fill = "black") +
  
  # Add group letters
  geom_text(data = va_emmeans_status_within_strata %>% filter(ghgspecies=="ch4"),
            aes(x = status, label = group_letter, y = Inf),
            vjust = 1.2, hjust = 0.5, color = "black", inherit.aes = FALSE,
            size = 5, fontface = "bold") +
  theme_bw() +
  scale_fill_manual(values = c(
    "Preserved" = "#009E73",
    "Altered" = "#D55E00",
    "Restored"  = "#56B4E9"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2), 0)) +
  theme(
    axis.text = element_text(face = "bold")
  ) +
  guides(color = "none", fill = "none") +
  labs(
    y = expression(CH[4]~NEE~(mmol~d^-1~m^-2)),
    x = NULL,
    fill = "Status"
  )+
  facet_wrap(facets=vars(strata), scales="free")

va_ch4_specific_plot


#GWP plot: 
va_gwp_specific_plot<-
  ggplot(va_specific_data %>% filter(ghgspecies=="GWPco2andch4"), 
         aes(x = status, y = dailyflux)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(width = 0.2, outliers = F, aes(fill = status), size = 0.7) +
  # Add emmeans points
  geom_point(data = va_emmeans_status_within_strata %>% filter(ghgspecies=="GWPco2andch4"),
             aes(x = status, y = emmean_bt),
             shape = 23, size = 3, fill = "black") +
  # Add group letters
  geom_text(data = va_emmeans_status_within_strata %>% filter(ghgspecies=="GWPco2andch4"),
            aes(x = status, label = group_letter, y = Inf),
            vjust = 1.2, hjust = 0.5, color = "black", inherit.aes = FALSE,
            size = 5, fontface = "bold") +
  theme_bw() +
  scale_fill_manual(values = c(
    "Preserved" = "#009E73",
    "Altered" = "#D55E00",
    "Restored"  = "#56B4E9"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2), 0)) +
  theme(
    axis.text = element_text(face = "bold")
  ) +
  guides(color = "none", fill = "none") +
  labs(
    y = expression(GWP~NEE~(g*CO[2~eq]~d^-1~m^-2)),
    x = "Conservation status",
    fill = "Status"
  )+
  facet_wrap(facets=vars(strata), scales="free")

va_gwp_specific_plot

#COMBINE pannels for each GHG
va_specific_plot_3ghg <- plot_grid(
  va_co2_specific_plot,
  va_ch4_specific_plot,
  va_gwp_specific_plot,
  ncol = 1,
  align = "v",
  axis = "lr"
  # labels = c("A", "B", "C"),
  # label_size = 12,
  # label_x = 0.02,  # move label closer to left
  # label_y = 1      # top-aligned
)

va_specific_plot_3ghg

#Save:
ggsave(plot = va_specific_plot_3ghg, 
       filename = "SpecificModels_va_3ghg_2levelstrata.png",
       path = sup_figures,
       device = "png",
       dpi = 400,
       units = "mm",
       height = 200, 
       width = 150)



#OVERALL status effect for VA with updated model: 

##VA_specific_statusonly----

va_emmeans_statusonly<- emmeans_status_specificmodels %>% 
  filter(casepilot=="VA") %>% 
  mutate(status=factor(status, levels = c("Preserved","Altered", "Restored"), ordered = T))


#CO2 plot: 
va_co2_status_plot<-
  ggplot(va_specific_data %>% filter(ghgspecies=="co2"), 
         aes(x = status, y = dailyflux * 1000)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(width = 0.2, outliers = F, aes(fill = status), size = 0.7) +
  
  # Add emmeans points
  geom_point(data = va_emmeans_statusonly %>% filter(ghgspecies=="co2"),
             aes(x = status, y = emmean_bt * 1000),
             shape = 23, size = 3, fill = "black") +
  
  # Add group letters
  geom_text(data = va_emmeans_statusonly %>% filter(ghgspecies=="co2"),
            aes(x = status, label = group_letter, y = Inf),
            vjust = 1.2, hjust = 0.5, color = "black", inherit.aes = FALSE,
            size = 5, fontface = "bold") +
  theme_bw() +
  scale_fill_manual(values = c(
    "Preserved" = "#009E73",
    "Altered" = "#D55E00",
    "Restored"  = "#56B4E9"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2), 0)) +
  theme(
    axis.text = element_text(face = "bold")
  ) +
  guides(color = "none", fill = "none") +
  labs(
    y = expression(CO[2]~NEE~(mmol~d^-1~m^-2)),
    x = NULL,
    fill = "Status"
  )

va_co2_status_plot



#CH4 plot: 
va_ch4_status_plot<-
  ggplot(va_specific_data %>% filter(ghgspecies=="ch4"), 
         aes(x = status, y = dailyflux)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(width = 0.2, outliers = F, aes(fill = status), size = 0.7) +
  
  # Add emmeans points
  geom_point(data = va_emmeans_statusonly %>% filter(ghgspecies=="ch4"),
             aes(x = status, y = emmean_bt),
             shape = 23, size = 3, fill = "black") +
  
  # Add group letters
  geom_text(data = va_emmeans_statusonly %>% filter(ghgspecies=="ch4"),
            aes(x = status, label = group_letter, y = Inf),
            vjust = 1.2, hjust = 0.5, color = "black", inherit.aes = FALSE,
            size = 5, fontface = "bold") +
  theme_bw() +
  scale_fill_manual(values = c(
    "Preserved" = "#009E73",
    "Altered" = "#D55E00",
    "Restored"  = "#56B4E9"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2), 0)) +
  theme(
    axis.text = element_text(face = "bold")
  ) +
  guides(color = "none", fill = "none") +
  labs(
    y = expression(CH[4]~NEE~(mmol~d^-1~m^-2)),
    x = NULL,
    fill = "Status"
  )

va_ch4_status_plot


#GWP plot: 
va_gwp_status_plot<-
  ggplot(va_specific_data %>% filter(ghgspecies=="GWPco2andch4"), 
         aes(x = status, y = dailyflux)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(width = 0.2, outliers = F, aes(fill = status), size = 0.7) +
  # Add emmeans points
  geom_point(data = va_emmeans_statusonly %>% filter(ghgspecies=="GWPco2andch4"),
             aes(x = status, y = emmean_bt),
             shape = 23, size = 3, fill = "black") +
  # Add group letters
  geom_text(data = va_emmeans_statusonly %>% filter(ghgspecies=="GWPco2andch4"),
            aes(x = status, label = group_letter, y = Inf),
            vjust = 1.2, hjust = 0.5, color = "black", inherit.aes = FALSE,
            size = 5, fontface = "bold") +
  theme_bw() +
  scale_fill_manual(values = c(
    "Preserved" = "#009E73",
    "Altered" = "#D55E00",
    "Restored"  = "#56B4E9"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2), 0)) +
  theme(
    axis.text = element_text(face = "bold")
  ) +
  guides(color = "none", fill = "none") +
  labs(
    y = expression(GWP~NEE~(g*CO[2~eq]~d^-1~m^-2)),
    x = "Conservation status",
    fill = "Status"
  )

va_gwp_status_plot

#COMBINE pannels for each GHG
va_specific_statusonlyplot_3ghg <- plot_grid(
  va_co2_status_plot,
  va_ch4_status_plot,
  va_gwp_status_plot,
  ncol = 1,
  align = "v",
  axis = "lr"
  # labels = c("A", "B", "C"),
  # label_size = 12,
  # label_x = 0.02,  # move label closer to left
  # label_y = 1      # top-aligned
)

va_specific_statusonlyplot_3ghg

ggsave(plot = va_specific_statusonlyplot_3ghg, 
       filename = "SpecificModels_VA_3ghg_generalstatus.png",
       path = sup_figures,
       device = "png",
       dpi = 400,
       units = "mm",
       height = 200, 
       width = 90)




##CA considering strata------

#Fitst need to re-level strata  
ca_specific_data <- data4paper %>% 
  filter(casepilot=="CA") %>% 
  mutate(status=factor(status, levels = c("Preserved","Altered", "Restored"), ordered = T),
         strata=as.factor(if_else(strata=="vegetated","vegetated", "non-vegetated")))


head(ca_specific_data)
str(ca_specific_data)

ca_emmeans_status_within_strata<- emmeans_status_within_strata_specificmodels %>%
  filter(casepilot=="CA")%>% 
  mutate(status=factor(status, levels = c("Preserved","Altered", "Restored"), ordered = T),
         strata=as.factor(if_else(strata=="vegetated","vegetated", "non-vegetated")))

str(ca_emmeans_status_within_strata)


#CO2 plot: 
ca_co2_specific_plot<-
  ggplot(ca_specific_data %>% filter(ghgspecies=="co2"), 
         aes(x = status, y = dailyflux * 1000)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(width = 0.2, outliers = F, aes(fill = status), size = 0.7) +
  
  # Add emmeans points
  geom_point(data = ca_emmeans_status_within_strata %>% filter(ghgspecies=="co2"),
             aes(x = status, y = emmean_bt * 1000),
             shape = 23, size = 3, fill = "black") +
  
  # Add group letters
  geom_text(data = ca_emmeans_status_within_strata %>% filter(ghgspecies=="co2"),
            aes(x = status, label = group_letter, y = Inf),
            vjust = 1.2, hjust = 0.5, color = "black", inherit.aes = FALSE,
            size = 5, fontface = "bold") +
  theme_bw() +
  scale_fill_manual(values = c(
    "Preserved" = "#009E73",
    "Altered" = "#D55E00",
    "Restored"  = "#56B4E9"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2), 0)) +
  theme(
    axis.text = element_text(face = "bold")
  ) +
  guides(color = "none", fill = "none") +
  labs(
    y = expression(CO[2]~NEE~(mmol~d^-1~m^-2)),
    x = NULL,
    fill = "Status"
  )+
  facet_wrap(facets=vars(strata), scales="free")

ca_co2_specific_plot



#CH4 plot: 
ca_ch4_specific_plot<-
  ggplot(ca_specific_data %>% filter(ghgspecies=="ch4"), 
         aes(x = status, y = dailyflux)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(width = 0.2, outliers = F, aes(fill = status), size = 0.7) +
  
  # Add emmeans points
  geom_point(data = ca_emmeans_status_within_strata %>% filter(ghgspecies=="ch4"),
             aes(x = status, y = emmean_bt),
             shape = 23, size = 3, fill = "black") +
  
  # Add group letters
  geom_text(data = ca_emmeans_status_within_strata %>% filter(ghgspecies=="ch4"),
            aes(x = status, label = group_letter, y = Inf),
            vjust = 1.2, hjust = 0.5, color = "black", inherit.aes = FALSE,
            size = 5, fontface = "bold") +
  theme_bw() +
  scale_fill_manual(values = c(
    "Preserved" = "#009E73",
    "Altered" = "#D55E00",
    "Restored"  = "#56B4E9"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2), 0)) +
  theme(
    axis.text = element_text(face = "bold")
  ) +
  guides(color = "none", fill = "none") +
  labs(
    y = expression(CH[4]~NEE~(mmol~d^-1~m^-2)),
    x = NULL,
    fill = "Status"
  )+
  facet_wrap(facets=vars(strata), scales="free")

ca_ch4_specific_plot


#GWP plot: 
ca_gwp_specific_plot<-
  ggplot(ca_specific_data %>% filter(ghgspecies=="GWPco2andch4"), 
         aes(x = status, y = dailyflux)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(width = 0.2, outliers = F, aes(fill = status), size = 0.7) +
  # Add emmeans points
  geom_point(data = ca_emmeans_status_within_strata %>% filter(ghgspecies=="GWPco2andch4"),
             aes(x = status, y = emmean_bt),
             shape = 23, size = 3, fill = "black") +
  # Add group letters
  geom_text(data = ca_emmeans_status_within_strata %>% filter(ghgspecies=="GWPco2andch4"),
            aes(x = status, label = group_letter, y = Inf),
            vjust = 1.2, hjust = 0.5, color = "black", inherit.aes = FALSE,
            size = 5, fontface = "bold") +
  theme_bw() +
  scale_fill_manual(values = c(
    "Preserved" = "#009E73",
    "Altered" = "#D55E00",
    "Restored"  = "#56B4E9"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2), 0)) +
  theme(
    axis.text = element_text(face = "bold")
  ) +
  guides(color = "none", fill = "none") +
  labs(
    y = expression(GWP~NEE~(g*CO[2~eq]~d^-1~m^-2)),
    x = "Conservation status",
    fill = "Status"
  )+
  facet_wrap(facets=vars(strata), scales="free")

ca_gwp_specific_plot

#COMBINE pannels for each GHG
ca_specific_plot_3ghg <- plot_grid(
  ca_co2_specific_plot,
  ca_ch4_specific_plot,
  ca_gwp_specific_plot,
  ncol = 1,
  align = "v",
  axis = "lr"
  # labels = c("A", "B", "C"),
  # label_size = 12,
  # label_x = 0.02,  # move label closer to left
  # label_y = 1      # top-aligned
)

ca_specific_plot_3ghg

#Save:
ggsave(plot = ca_specific_plot_3ghg, 
       filename = "SpecificModels_ca_3ghg_2levelstrata.png",
       path = sup_figures,
       device = "png",
       dpi = 400,
       units = "mm",
       height = 200, 
       width = 150)



#OVERALL status effect for CA with updated model: 

##CA_specific_statusonly----

ca_emmeans_statusonly<- emmeans_status_specificmodels %>% 
  filter(casepilot=="CA") %>% 
  mutate(status=factor(status, levels = c("Preserved","Altered", "Restored"), ordered = T))


#CO2 plot: 
ca_co2_status_plot<-
  ggplot(ca_specific_data %>% filter(ghgspecies=="co2"), 
         aes(x = status, y = dailyflux * 1000)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(width = 0.2, outliers = F, aes(fill = status), size = 0.7) +
  
  # Add emmeans points
  geom_point(data = ca_emmeans_statusonly %>% filter(ghgspecies=="co2"),
             aes(x = status, y = emmean_bt * 1000),
             shape = 23, size = 3, fill = "black") +
  
  # Add group letters
  geom_text(data = ca_emmeans_statusonly %>% filter(ghgspecies=="co2"),
            aes(x = status, label = group_letter, y = Inf),
            vjust = 1.2, hjust = 0.5, color = "black", inherit.aes = FALSE,
            size = 5, fontface = "bold") +
  theme_bw() +
  scale_fill_manual(values = c(
    "Preserved" = "#009E73",
    "Altered" = "#D55E00",
    "Restored"  = "#56B4E9"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2), 0)) +
  theme(
    axis.text = element_text(face = "bold")
  ) +
  guides(color = "none", fill = "none") +
  labs(
    y = expression(CO[2]~NEE~(mmol~d^-1~m^-2)),
    x = NULL,
    fill = "Status"
  )

ca_co2_status_plot



#CH4 plot: 
ca_ch4_status_plot<-
  ggplot(ca_specific_data %>% filter(ghgspecies=="ch4"), 
         aes(x = status, y = dailyflux)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(width = 0.2, outliers = F, aes(fill = status), size = 0.7) +
  
  # Add emmeans points
  geom_point(data = ca_emmeans_statusonly %>% filter(ghgspecies=="ch4"),
             aes(x = status, y = emmean_bt),
             shape = 23, size = 3, fill = "black") +
  
  # Add group letters
  geom_text(data = ca_emmeans_statusonly %>% filter(ghgspecies=="ch4"),
            aes(x = status, label = group_letter, y = Inf),
            vjust = 1.2, hjust = 0.5, color = "black", inherit.aes = FALSE,
            size = 5, fontface = "bold") +
  theme_bw() +
  scale_fill_manual(values = c(
    "Preserved" = "#009E73",
    "Altered" = "#D55E00",
    "Restored"  = "#56B4E9"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2), 0)) +
  theme(
    axis.text = element_text(face = "bold")
  ) +
  guides(color = "none", fill = "none") +
  labs(
    y = expression(CH[4]~NEE~(mmol~d^-1~m^-2)),
    x = NULL,
    fill = "Status"
  )

ca_ch4_status_plot


#GWP plot: 
ca_gwp_status_plot<-
  ggplot(ca_specific_data %>% filter(ghgspecies=="GWPco2andch4"), 
         aes(x = status, y = dailyflux)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(width = 0.2, outliers = F, aes(fill = status), size = 0.7) +
  # Add emmeans points
  geom_point(data = ca_emmeans_statusonly %>% filter(ghgspecies=="GWPco2andch4"),
             aes(x = status, y = emmean_bt),
             shape = 23, size = 3, fill = "black") +
  # Add group letters
  geom_text(data = ca_emmeans_statusonly %>% filter(ghgspecies=="GWPco2andch4"),
            aes(x = status, label = group_letter, y = Inf),
            vjust = 1.2, hjust = 0.5, color = "black", inherit.aes = FALSE,
            size = 5, fontface = "bold") +
  theme_bw() +
  scale_fill_manual(values = c(
    "Preserved" = "#009E73",
    "Altered" = "#D55E00",
    "Restored"  = "#56B4E9"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2), 0)) +
  theme(
    axis.text = element_text(face = "bold")
  ) +
  guides(color = "none", fill = "none") +
  labs(
    y = expression(GWP~NEE~(g*CO[2~eq]~d^-1~m^-2)),
    x = "Conservation status",
    fill = "Status"
  )

ca_gwp_status_plot

#COMBINE pannels for each GHG
ca_specific_statusonlyplot_3ghg <- plot_grid(
  ca_co2_status_plot,
  ca_ch4_status_plot,
  ca_gwp_status_plot,
  ncol = 1,
  align = "v",
  axis = "lr"
  # labels = c("A", "B", "C"),
  # label_size = 12,
  # label_x = 0.02,  # move label closer to left
  # label_y = 1      # top-aligned
)

ca_specific_statusonlyplot_3ghg

ggsave(plot = ca_specific_statusonlyplot_3ghg, 
       filename = "SpecificModels_CA_3ghg_generalstatus.png",
       path = sup_figures,
       device = "png",
       dpi = 400,
       units = "mm",
       height = 200, 
       width = 90)




##DU considering strata------

#Fitst need to re-level strata  
du_specific_data <- data4paper %>% 
  filter(casepilot=="DU") %>% 
  mutate(status=factor(status, levels = c("Preserved","Altered", "Restored"), ordered = T),
         strata=as.factor(strata))


head(du_specific_data)
str(du_specific_data)

du_emmeans_status_within_strata<- emmeans_status_within_strata_specificmodels %>%
  filter(casepilot=="DU")%>% 
  mutate(status=factor(status, levels = c("Preserved","Altered", "Restored"), ordered = T),
         strata=as.factor(strata))

str(du_emmeans_status_within_strata)


#CO2 plot: 
du_co2_specific_plot<-
  ggplot(du_specific_data %>% filter(ghgspecies=="co2"), 
         aes(x = status, y = dailyflux * 1000)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(width = 0.2, outliers = F, aes(fill = status), size = 0.7) +
  
  # Add emmeans points
  geom_point(data = du_emmeans_status_within_strata %>% filter(ghgspecies=="co2"),
             aes(x = status, y = emmean_bt * 1000),
             shape = 23, size = 3, fill = "black") +
  
  # Add group letters
  geom_text(data = du_emmeans_status_within_strata %>% filter(ghgspecies=="co2"),
            aes(x = status, label = group_letter, y = Inf),
            vjust = 1.2, hjust = 0.5, color = "black", inherit.aes = FALSE,
            size = 5, fontface = "bold") +
  theme_bw() +
  scale_fill_manual(values = c(
    "Preserved" = "#009E73",
    "Altered" = "#D55E00",
    "Restored"  = "#56B4E9"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2), 0)) +
  theme(
    axis.text = element_text(face = "bold")
  ) +
  guides(color = "none", fill = "none") +
  labs(
    y = expression(CO[2]~NEE~(mmol~d^-1~m^-2)),
    x = NULL,
    fill = "Status"
  )+
  facet_wrap(facets=vars(strata), scales="free")

du_co2_specific_plot



#CH4 plot: 
du_ch4_specific_plot<-
  ggplot(du_specific_data %>% filter(ghgspecies=="ch4"), 
         aes(x = status, y = dailyflux)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(width = 0.2, outliers = F, aes(fill = status), size = 0.7) +
  
  # Add emmeans points
  geom_point(data = du_emmeans_status_within_strata %>% filter(ghgspecies=="ch4"),
             aes(x = status, y = emmean_bt),
             shape = 23, size = 3, fill = "black") +
  
  # Add group letters
  geom_text(data = du_emmeans_status_within_strata %>% filter(ghgspecies=="ch4"),
            aes(x = status, label = group_letter, y = Inf),
            vjust = 1.2, hjust = 0.5, color = "black", inherit.aes = FALSE,
            size = 5, fontface = "bold") +
  theme_bw() +
  scale_fill_manual(values = c(
    "Preserved" = "#009E73",
    "Altered" = "#D55E00",
    "Restored"  = "#56B4E9"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2), 0)) +
  theme(
    axis.text = element_text(face = "bold")
  ) +
  guides(color = "none", fill = "none") +
  labs(
    y = expression(CH[4]~NEE~(mmol~d^-1~m^-2)),
    x = NULL,
    fill = "Status"
  )+
  facet_wrap(facets=vars(strata), scales="free")

du_ch4_specific_plot


#GWP plot: 
du_gwp_specific_plot<-
  ggplot(du_specific_data %>% filter(ghgspecies=="GWPco2andch4"), 
         aes(x = status, y = dailyflux)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(width = 0.2, outliers = F, aes(fill = status), size = 0.7) +
  # Add emmeans points
  geom_point(data = du_emmeans_status_within_strata %>% filter(ghgspecies=="GWPco2andch4"),
             aes(x = status, y = emmean_bt),
             shape = 23, size = 3, fill = "black") +
  # Add group letters
  geom_text(data = du_emmeans_status_within_strata %>% filter(ghgspecies=="GWPco2andch4"),
            aes(x = status, label = group_letter, y = Inf),
            vjust = 1.2, hjust = 0.5, color = "black", inherit.aes = FALSE,
            size = 5, fontface = "bold") +
  theme_bw() +
  scale_fill_manual(values = c(
    "Preserved" = "#009E73",
    "Altered" = "#D55E00",
    "Restored"  = "#56B4E9"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2), 0)) +
  theme(
    axis.text = element_text(face = "bold")
  ) +
  guides(color = "none", fill = "none") +
  labs(
    y = expression(GWP~NEE~(g*CO[2~eq]~d^-1~m^-2)),
    x = "Conservation status",
    fill = "Status"
  )+
  facet_wrap(facets=vars(strata), scales="free")

du_gwp_specific_plot

#COMBINE pannels for each GHG
du_specific_plot_3ghg <- plot_grid(
  du_co2_specific_plot,
  du_ch4_specific_plot,
  du_gwp_specific_plot,
  ncol = 1,
  align = "v",
  axis = "lr"
  # labels = c("A", "B", "C"),
  # label_size = 12,
  # label_x = 0.02,  # move label closer to left
  # label_y = 1      # top-aligned
)

du_specific_plot_3ghg

#Save:
ggsave(plot = du_specific_plot_3ghg, 
       filename = "SpecificModels_du_3ghg_2levelstrata.png",
       path = sup_figures,
       device = "png",
       dpi = 400,
       units = "mm",
       height = 200, 
       width = 150)



#OVERALL status effect for DU with updated model: 

##DU_specific_statusonly----

du_emmeans_statusonly<- emmeans_status_specificmodels %>% 
  filter(casepilot=="DU") %>% 
  mutate(status=factor(status, levels = c("Preserved","Altered", "Restored"), ordered = T))


#CO2 plot: 
du_co2_status_plot<-
  ggplot(du_specific_data %>% filter(ghgspecies=="co2"), 
         aes(x = status, y = dailyflux * 1000)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(width = 0.2, outliers = F, aes(fill = status), size = 0.7) +
  
  # Add emmeans points
  geom_point(data = du_emmeans_statusonly %>% filter(ghgspecies=="co2"),
             aes(x = status, y = emmean_bt * 1000),
             shape = 23, size = 3, fill = "black") +
  
  # Add group letters
  geom_text(data = du_emmeans_statusonly %>% filter(ghgspecies=="co2"),
            aes(x = status, label = group_letter, y = Inf),
            vjust = 1.2, hjust = 0.5, color = "black", inherit.aes = FALSE,
            size = 5, fontface = "bold") +
  theme_bw() +
  scale_fill_manual(values = c(
    "Preserved" = "#009E73",
    "Altered" = "#D55E00",
    "Restored"  = "#56B4E9"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2), 0)) +
  theme(
    axis.text = element_text(face = "bold")
  ) +
  guides(color = "none", fill = "none") +
  labs(
    y = expression(CO[2]~NEE~(mmol~d^-1~m^-2)),
    x = NULL,
    fill = "Status"
  )

du_co2_status_plot



#CH4 plot: 
du_ch4_status_plot<-
  ggplot(du_specific_data %>% filter(ghgspecies=="ch4"), 
         aes(x = status, y = dailyflux)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(width = 0.2, outliers = F, aes(fill = status), size = 0.7) +
  
  # Add emmeans points
  geom_point(data = du_emmeans_statusonly %>% filter(ghgspecies=="ch4"),
             aes(x = status, y = emmean_bt),
             shape = 23, size = 3, fill = "black") +
  
  # Add group letters
  geom_text(data = du_emmeans_statusonly %>% filter(ghgspecies=="ch4"),
            aes(x = status, label = group_letter, y = Inf),
            vjust = 1.2, hjust = 0.5, color = "black", inherit.aes = FALSE,
            size = 5, fontface = "bold") +
  theme_bw() +
  scale_fill_manual(values = c(
    "Preserved" = "#009E73",
    "Altered" = "#D55E00",
    "Restored"  = "#56B4E9"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2), 0)) +
  theme(
    axis.text = element_text(face = "bold")
  ) +
  guides(color = "none", fill = "none") +
  labs(
    y = expression(CH[4]~NEE~(mmol~d^-1~m^-2)),
    x = NULL,
    fill = "Status"
  )

du_ch4_status_plot


#GWP plot: 
du_gwp_status_plot<-
  ggplot(du_specific_data %>% filter(ghgspecies=="GWPco2andch4"), 
         aes(x = status, y = dailyflux)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(width = 0.2, outliers = F, aes(fill = status), size = 0.7) +
  # Add emmeans points
  geom_point(data = du_emmeans_statusonly %>% filter(ghgspecies=="GWPco2andch4"),
             aes(x = status, y = emmean_bt),
             shape = 23, size = 3, fill = "black") +
  # Add group letters
  geom_text(data = du_emmeans_statusonly %>% filter(ghgspecies=="GWPco2andch4"),
            aes(x = status, label = group_letter, y = Inf),
            vjust = 1.2, hjust = 0.5, color = "black", inherit.aes = FALSE,
            size = 5, fontface = "bold") +
  theme_bw() +
  scale_fill_manual(values = c(
    "Preserved" = "#009E73",
    "Altered" = "#D55E00",
    "Restored"  = "#56B4E9"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2), 0)) +
  theme(
    axis.text = element_text(face = "bold")
  ) +
  guides(color = "none", fill = "none") +
  labs(
    y = expression(GWP~NEE~(g*CO[2~eq]~d^-1~m^-2)),
    x = "Conservation status",
    fill = "Status"
  )

du_gwp_status_plot

#COMBINE pannels for each GHG
du_specific_statusonlyplot_3ghg <- plot_grid(
  du_co2_status_plot,
  du_ch4_status_plot,
  du_gwp_status_plot,
  ncol = 1,
  align = "v",
  axis = "lr"
  # labels = c("A", "B", "C"),
  # label_size = 12,
  # label_x = 0.02,  # move label closer to left
  # label_y = 1      # top-aligned
)

du_specific_statusonlyplot_3ghg

ggsave(plot = du_specific_statusonlyplot_3ghg, 
       filename = "SpecificModels_DU_3ghg_generalstatus.png",
       path = sup_figures,
       device = "png",
       dpi = 400,
       units = "mm",
       height = 200, 
       width = 90)









#Strata comp. plots-----
strata_colors <- c(
  "bare" = "#d95f02",      # brown/orange
  "vegetated" = "#228B22", # green
  "open water" = "#1f78b4" # blue
)

##CA strata_prop plot
strataplot_ca<- stratacomp %>% filter(casepilot == "CA") %>% 
  ggplot(aes(x = season, y = strata_prop*100, fill = strata)) +
  geom_bar(stat = "identity", col="black") +
  facet_wrap(~subsite) +
  scale_fill_manual(values = strata_colors, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format(scale = 1),expand=c(0,0)) +
  labs(y = "Strata Chambers (%)", x = "Season", fill = "Strata") +
  theme_minimal(base_line_size = 0.5)

ggsave(plot = strataplot_ca,filename = "Strata_comp_CA.png",
       width = 115,height = 100,units = "mm",device = "png",dpi = 400,
       path = extra_figures)

##CU strata_prop plot
strataplot_cu<- stratacomp %>% filter(casepilot == "CU") %>% 
  ggplot(aes(x = season, y = strata_prop*100, fill = strata)) +
  geom_bar(stat = "identity", col="black") +
  facet_wrap(~subsite) +
  scale_fill_manual(values = strata_colors, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format(scale = 1),expand=c(0,0)) +
  labs(y = "Strata Chambers (%)", x = "Season", fill = "Strata") +
  theme_minimal()

ggsave(plot = strataplot_cu,filename = "Strata_comp_CU.png",
       width = 115,height = 100,units = "mm",device = "png",dpi = 400,
       path = extra_figures)

##DA strata_prop plot
strataplot_da<- stratacomp %>% filter(casepilot == "DA") %>% 
  ggplot(aes(x = season, y = strata_prop*100, fill = strata)) +
  geom_bar(stat = "identity", col="black") +
  facet_wrap(~subsite) +
  scale_fill_manual(values = strata_colors, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format(scale = 1),expand=c(0,0)) +
  labs(y = "Strata Chambers (%)", x = "Season", fill = "Strata") +
  theme_minimal()

ggsave(plot = strataplot_da,filename = "Strata_comp_DA.png",
       width = 115,height = 100,units = "mm",device = "png",dpi = 400,
       path = extra_figures)

##DU strata_prop plot
strataplot_du<- stratacomp %>% filter(casepilot == "DU") %>% 
  ggplot(aes(x = season, y = strata_prop*100, fill = strata)) +
  geom_bar(stat = "identity", col="black") +
  facet_wrap(~subsite) +
  scale_fill_manual(values = strata_colors, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format(scale = 1),expand=c(0,0)) +
  labs(y = "Strata Chambers (%)", x = "Season", fill = "Strata") +
  theme_minimal()

ggsave(plot = strataplot_du,filename = "Strata_comp_DU.png",
       width = 115,height = 100,units = "mm",device = "png",dpi = 400,
       path = extra_figures)


##RI strata_prop plot
strataplot_ri<- stratacomp %>% filter(casepilot == "RI") %>% 
  ggplot(aes(x = season, y = strata_prop*100, fill = strata)) +
  geom_bar(stat = "identity", col="black") +
  facet_wrap(~subsite) +
  scale_fill_manual(values = strata_colors, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format(scale = 1),expand=c(0,0)) +
  labs(y = "Strata Chambers (%)", x = "Season", fill = "Strata") +
  theme_minimal()

ggsave(plot = strataplot_ri,filename = "Strata_comp_RI.png",
       width = 115,height = 100,units = "mm",device = "png",dpi = 400,
       path = extra_figures)


##VA strata_prop plot
strataplot_va<- stratacomp %>% filter(casepilot == "VA") %>% 
  ggplot(aes(x = season, y = strata_prop*100, fill = strata)) +
  geom_bar(stat = "identity", col="black") +
  facet_wrap(~subsite) +
  scale_fill_manual(values = strata_colors, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format(scale = 1),expand=c(0,0)) +
  labs(y = "Strata Chambers (%)", x = "Season", fill = "Strata") +
  theme_minimal()

ggsave(plot = strataplot_va,filename = "Strata_comp_VA.png",
       width = 115,height = 100,units = "mm",device = "png",dpi = 400,
       path = extra_figures)





#Status pseudolog scale------

#Using geom_sina and geom_violin to show all data, With scale transformation of pseudolog 

status_plot_pseudolog <- function(GHG,
                             y_label, 
                             y_multiplier = 1) {
  
  #filter data for GHG
  data<- data4paper %>% filter(ghgspecies==GHG)
  emmeans_data<- emmeans_status %>% filter(ghgspecies==GHG) 
  
  #Produce plot: 
  ggplot(data, aes(x = status, y = dailyflux * y_multiplier)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    
    geom_violin(aes(col=status),position = position_dodge(0.9),linewidth = 1,scale = "width")+
    geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
    # geom_boxplot(width = 0.2, outliers = drawboxplot_outliers, aes(fill = status), size = 0.7) +
    
    # Add emmeans points
    geom_point(data = emmeans_data, aes(x = status, y = emmean_bt * y_multiplier),
               shape = 23, size = 3, fill = "black") +
    
    # Add group letters
    geom_text(data = emmeans_data, aes(x = status, label = group_letter, y = Inf),
              vjust = 1.2, hjust = 0.5, color = "black", inherit.aes = FALSE,
              size = 5, fontface = "bold") +
    theme_bw() +
    scale_fill_manual(values = c(
      "Preserved" = "#009E73",
      "Altered"   = "#D55E00",
      "Restored"  = "#56B4E9"
    )) +
    scale_color_manual(values = c(
      "Preserved" = "#009E73",
      "Altered"   = "#D55E00",
      "Restored"  = "#56B4E9"
    )) +
    scale_y_continuous(trans = pseudo_log_trans(sigma = 0.1,base = 10),
                       guide = guide_axis_logticks(negative.small = 0.1),
                       expand = expansion(mult = c(0.1, 0.2), 0), #IF needed only
                       breaks = c(-1000,-100,-10,-1,-0,1,10,100,1000))+
    
    theme(
      axis.text = element_text(face = "bold")
    ) +
    guides(color = "none", fill = "none") +
    labs(
      y = y_label,
      x = "Conservation status",
      fill = "Status"
    ) +
    facet_grid(rows = vars(casepilot), scales = "free")
}


#

#CO2: (sigma=0.1),hour-glass shapes (not-intuitive visualization)
status_plot_pseudolog(GHG="co2",
                 y_label = expression(CO[2]~NEE~(mmol~d^-1~m^-2)),
                 y_multiplier = 1000)

ggsave(filename = "Alternative_status_CO2_pseudolog.png",
       path = extra_figures,
       device = "png",
       dpi = 400,
       units = "mm",
       height = 229,
       width = 90)

#CH4: (sigma=0.1), OK for CA,VA,DA,CU (good comparison), NOT appropriate for DU, RI
status_plot_pseudolog(GHG="ch4",
                      y_label = expression(CH[4]~NEE~(mmol~d^-1~m^-2)),
                      y_multiplier = 1)

ggsave(filename = "Alternative_status_CH4_pseudolog.png",
       path = extra_figures,
       device = "png",
       dpi = 400,
       units = "mm",
       height = 229,
       width = 90)

#CH4: (sigma=0.1), Hour-glass shape for mots, non-intuitive visualization
status_plot_pseudolog(GHG="GWPco2andch4",
                      y_label = expression(GWP~NEE~(g*CO[2~eq]~d^-1~m^-2)),
                      y_multiplier = 1)

ggsave(filename = "Alternative_status_GWP_pseudolog.png",
       path = extra_figures,
       device = "png",
       dpi = 400,
       units = "mm",
       height = 229,
       width = 90)




#____________________--------


#OLD PLOTS (NOT TO USE)-------

##CO2 --------
#Vertical plot CO2 molar units panels per case-pilot, independent scales, emmeans + letters. 

#Visualization issues: Each casepilot would benefit from its own scale (different compression ranges for logsign, interaction logtransformation & units). However, I think this would be very confusing. I could potentially show only boxplot (no outliers) + emmean, this "hides" weird distribution of data and makes the plots much clearer. 

co2_4paper<- data4paper %>% filter(ghgspecies=="co2")

co2_emmeans_4paper<-  emmeans_all %>% filter(comparison=="status",ghgspecies=="co2")


#Identity scale: boxplot without outliers works ok
ggplot(co2_4paper, aes(x=status, y=dailyflux*1000))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  #ADD EMMEANS and letters
  geom_point(data=co2_emmeans_4paper, aes(x=status,y=emmean_bt*1000),
             shape=23,size=3,fill = "black")+
  geom_text(data=co2_emmeans_4paper, aes(x=status, label = group_letter, y=Inf),
            vjust=1.2,hjust=0.5, color = "black", inherit.aes = F, size=5, fontface="bold") +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0.1,0.2),0))+
  theme(
    axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs(y= expression(CO[2]~NEE~(mmol~d^-1~m^-2)),
       x=paste0("Conservation status"),
       fill=paste0("Status"))+
  facet_grid(rows="casepilot", scales="free")

ggsave(filename = "PAPERPLOTS_CO2_molar_boxplot_allCP.png.png",path = plots_path,device = "png",
       dpi = 400,
       units = "mm",
       height = 229,
       width = 90)









# save.image(paste0(plots_path,"/model_results_last.Rdata"))

# load(paste0(plots_path,"/model_results_last.Rdata"))


#0. Pseudolog effect. 
library(tidyverse)
library(scales)
library(cowplot)

# Dummy data
set.seed(42)


df <- tibble(
  casepilot = rep(rep(c("Low", "Mid", "High"), each = 200),each=2),
  dailyflux = c(heavy_tails <- rlnorm(600, meanlog = 2, sdlog = 1) * sample(c(-1, 1), 600, replace = TRUE), rlnorm(600, meanlog = 0.5, sdlog = 1)),
  type = rep(c("Centered", "LogNormal"), each = 600)
)


df <- tibble(casepilot = rep(rep(c("Low", "Mid", "High"), each = 200), 2), dailyflux = c(unlist(lapply(c(-10, 0, 10), function(m) rlnorm(200, meanlog = 2, sdlog = 1) * if (m == 0) sample(c(-1, 1), 200, replace = TRUE) else sign(m))), unlist(lapply(c(0.1, 1, 10), function(m) rlnorm(200, meanlog = log(m), sdlog = 1)))), type = rep(c("Centered", "LogNormal"), each = 600))

df <- tibble(casepilot = rep(rep(c("Low", "Mid", "High"), each = 200), 2), dailyflux = c(unlist(lapply(c(-10, 0, 10), function(m) rlnorm(200, meanlog = 2, sdlog = 1) * if (m == 0) sample(c(-1, 1), 200, replace = TRUE) else sign(m))), unlist(lapply(c(0.1, 1, 10), function(m) rlnorm(200, meanlog = log(m), sdlog = 1)))), type = rep(c("Centered", "LogNormal"), each = 600))


df <- tibble(casepilot = rep(rep(c("Low", "Mid", "High"), each = 200), 2), dailyflux = c(unlist(lapply(c(-10, 0, 10), function(shift) rlnorm(200, meanlog = 2, sdlog = 1) * sample(c(-1, 1), 200, replace = TRUE) + shift)), unlist(lapply(c(0.1, 1, 10), function(m) rlnorm(200, meanlog = log(m), sdlog = 1)))), type = rep(c("Centered", "LogNormal"), each = 600))

df <- tibble(casepilot = rep(rep(c("Low", "Mid", "High"), each = 200), 2), dailyflux = c(unlist(lapply(c(-10, 0, 10), function(shift) rlnorm(200, meanlog = 2, sdlog = 1) * sample(c(-1, 1), 200, replace = TRUE) + shift)), unlist(lapply(c(0.1, 1, 10), function(m) rlnorm(200, meanlog = log(m), sdlog = 1) - 0.007))), type = rep(c("Centered", "LogNormal"), each = 600))


# Sigma values to test
sigmas <- c(NA,1e-2, 1e-1, 1, 10)

# Function to generate one plot
make_plot <- function(data, sigma = NULL, title_suffix = "", linear = FALSE) {
  p <- ggplot(data, aes(x = casepilot, y = dailyflux)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_violin(position = position_dodge(0.9),linewidth = 1,scale = "width")+
    geom_boxplot(width = 0.2, fill = "#00BA38", size = 0.7) +
    
    theme_bw() +
    theme(
      axis.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5)
    ) +
    labs(
      y = expression(Flux~(units)),
      x = "Casepilot",
      title = if (linear) paste0("Linear scale (", title_suffix, ")") else paste0("Sigma = ", sigma, " (", title_suffix, ")")
    )
  
  if (!is.na(sigma)) {
    
    base_breaks <- c(1e-2, 1e-1, 1, 10, 100, 1000)
    pos_breaks <- base_breaks[base_breaks >= sigma]
    neg_breaks <- -rev(pos_breaks)
    breaks <- c(neg_breaks, 0, pos_breaks)
    
    
    
    p <- p + scale_y_continuous(
      trans = pseudo_log_trans(sigma = sigma, base = 10),
      guide = guide_axis_logticks(negative.small = sigma),
      breaks = breaks
    )
  }
  
  return(p)
}

plots_centered <- c(
  map(sigmas, ~make_plot(filter(df, type == "Centered"), sigma = .x, title_suffix = "Centered"))
)

plots_lognorm <- c(
  map(sigmas, ~make_plot(filter(df, type == "LogNormal"), sigma = .x, title_suffix = "LogNormal"))
)


# Arrange plots using cowplot
row1 <- plot_grid(plotlist = plots_centered, nrow = 1)
row2 <- plot_grid(plotlist = plots_lognorm, nrow = 1)
final_plot <- plot_grid(row1, row2, ncol = 1)

ggsave(plot = final_plot, 
       filename = "PAPERPLOTS_example_pseudolog_effect.png",path = plots_path,device = "png",
       dpi = 400,
       units = "mm",
       height = 200, 
       width = 300)



##With-out_linear-----
#A 3-panel vertical plot. 1 panel per gas. Boxplot (With outliers). Linear-scale

base_co2_ident_out<-data4paper %>% 
  filter(status=="Preserved", ghgspecies%in%c("co2")) %>% 
  mutate(casepilot=factor(casepilot, levels = c("DU","RI","CA","VA","DA","CU"))) %>% 
  ggplot(aes(x=casepilot, y=(dailyflux*1000)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_boxplot(width=0.2,outliers = T, fill="#00BA38",size = 0.7)+
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))+
  guides(color = "none", fill="none")+
  labs(y= expression(CO[2]~NEE~(mmol~d^-1~m^-2)),
       x=NULL,
       fill=paste0("Status"))


base_ch4_ident_out<-data4paper %>% 
  filter(status=="Preserved", ghgspecies%in%c("ch4")) %>% 
  mutate(casepilot=factor(casepilot, levels = c("DU","RI","CA","VA","DA","CU"))) %>% 
  ggplot(aes(x=casepilot, y=(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_boxplot(width=0.2,outliers = T, fill="#00BA38",size = 0.7)+
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))+
  guides(color = "none", fill="none")+
  labs(y= expression(CH[4]~NEE~(mmol~d^-1~m^-2)),
       x=NULL,
       fill=paste0("Status"))

base_gwp_ident_out<-data4paper %>% 
  filter(status=="Preserved", ghgspecies%in%c("gwp_co2andch4")) %>% 
  mutate(casepilot=factor(casepilot, levels = c("DU","RI","CA","VA","DA","CU"))) %>% 
  ggplot(aes(x=casepilot, y=(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_boxplot(width=0.2,outliers = T, fill="#00BA38",size = 0.7)+
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))+
  guides(color = "none", fill="none")+
  labs(y= expression(GWP~NEE~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Casepilot"),
       fill=paste0("Status"))


baseline_3sp_ident_out<-plot_grid(base_co2_ident_out,base_ch4_ident_out,base_gwp_ident_out,align = "v",nrow = 3)

ggsave(plot = baseline_3sp_ident_out, 
       filename = "PAPERPLOTS_baseline_3ghg_boxplot_with-outlier_indentity.png",path = plots_path,device = "png",
       dpi = 400,
       units = "mm",
       height = 200, 
       width = 90)

##With-out_pseudolog-----
#A 3-panel vertical plot. 1 panel per gas. Boxplot (With outliers). Pseudolog-scale (optimized)

base_co2_pseudolog<- data4paper %>% 
  filter(status=="Preserved", ghgspecies%in%c("co2")) %>% 
  mutate(casepilot=factor(casepilot, levels = c("DU","RI","CA","VA","DA","CU"))) %>% 
  ggplot(aes(x=casepilot, y=(dailyflux*1000)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_boxplot(width=0.2,outliers = T, fill="#00BA38",size = 0.7)+
  theme_bw() +
  scale_y_continuous(trans = pseudo_log_trans(sigma = 1,base = 10),
                     guide = guide_axis_logticks(negative.small = 1),
                     breaks = c(-1000,-100,-10,-0,10,100,1000))+
  theme(
    axis.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))+
  guides(color = "none", fill="none")+
  labs(y= expression(CO[2]~NEE~(mmol~d^-1~m^-2)),
       x=NULL,
       fill=paste0("Status"))



base_ch4_pseudolog<- data4paper %>% 
  filter(status=="Preserved", ghgspecies%in%c("ch4")) %>% 
  mutate(casepilot=factor(casepilot, levels = c("DU","RI","CA","VA","DA","CU"))) %>% 
  ggplot(aes(x=casepilot, y=(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_boxplot(width=0.2,outliers = F, fill="#00BA38",size = 0.7)+
  theme_bw() +
  scale_y_continuous(trans = pseudo_log_trans(sigma = 0.01,base = 10),
                     guide = guide_axis_logticks(negative.small = 0.01),
                     breaks = c(-10,-0.1,0,0.1,1,10,100,1000))+
  theme(
    axis.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))+
  guides(color = "none", fill="none")+
  labs(y= expression(CH[4]~NEE~(mmol~d^-1~m^-2)),
       x=NULL,
       fill=paste0("Status"))

base_gwp_pseudolog<- data4paper %>% 
  filter(status=="Preserved", ghgspecies%in%c("gwp_co2andch4")) %>% 
  mutate(casepilot=factor(casepilot, levels = c("DU","RI","CA","VA","DA","CU"))) %>% 
  ggplot(aes(x=casepilot, y=(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_boxplot(width=0.2,outliers = F, fill="#00BA38",size = 0.7)+
  theme_bw() +
  scale_y_continuous(trans = pseudo_log_trans(sigma = 1e-1,base = 10),
                     guide = guide_axis_logticks(negative.small =  1e-1),
                     breaks = c(-100,-10,-1,-0.1,0,0.1,1,10,100,1000))+
  theme(
    axis.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))+
  guides(color = "none", fill="none")+
  labs(y= expression(CH[4]~NEE~(mmol~d^-1~m^-2)),
       x=paste0("Casepilot"),
       fill=paste0("Status"))

baseline_3sp_pseudolog<-plot_grid(base_co2_pseudolog,base_ch4_pseudolog,base_gwp_pseudolog,align = "v",nrow = 3)

ggsave(plot = baseline_3sp_pseudolog, 
       filename = "PAPERPLOTS_baseline_3ghg_boxplot_with-outlier_pseudolog.png",path = plots_path,device = "png",
       dpi = 400,
       units = "mm",
       height = 200, 
       width = 90)






#PLOTS per GHG ------



##CO2 --------
#Vertical plot CO2 molar units panels per case-pilot, independent scales, emmeans + letters. 

#Visualization issues: Each casepilot would benefit from its own scale (different compression ranges for logsign, interaction logtransformation & units). However, I think this would be very confusing. I could potentially show only boxplot (no outliers) + emmean, this "hides" weird distribution of data and makes the plots much clearer. 

co2_4paper<- data4paper %>% filter(ghgspecies=="co2") %>%   mutate(casepilot=factor(casepilot, levels = c("DU","RI","CA","VA","DA","CU")))
co2_emmeans_4paper<-  emmeans_all %>% filter(comparison=="status",ghgspecies=="co2")%>%   
  mutate(casepilot=factor(casepilot, levels = c("DU","RI","CA","VA","DA","CU")))


#Identity scale: boxplot without outliers works ok
ggplot(co2_4paper, aes(x=status, y=dailyflux*1000))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  #ADD EMMEANS and letters
  geom_point(data=co2_emmeans_4paper, aes(x=status,y=emmean_bt*1000),
             shape=23,size=3,fill = "black")+
  geom_text(data=co2_emmeans_4paper, aes(x=status, label = group_letter, y=Inf),
            vjust=1.2,hjust=0.5, color = "black", inherit.aes = F, size=5, fontface="bold") +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0.1,0.2),0))+
  theme(
    axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs(y= expression(CO[2]~NEE~(mmol~d^-1~m^-2)),
       x=paste0("Conservation status"),
       fill=paste0("Status"))+
  facet_grid(rows="casepilot", scales="free")

ggsave(filename = "PAPERPLOTS_CO2_molar_boxplot_allCP.png.png",path = plots_path,device = "png",
       dpi = 400,
       units = "mm",
       height = 229,
       width = 90)


#Identity scale: boxplot(no outliers) + sina (all data), DOES not work, outliers force scale so that group differences are not clear. 
ggplot(co2_4paper, aes(x=status, y=dailyflux*1000))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(aes(col=status),position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  #ADD EMMEANS and letters
  geom_point(data=co2_emmeans_4paper, aes(x=status,y=emmean_bt*1000),
             shape=23,size=3,fill = "black")+
  geom_text(data=co2_emmeans_4paper, aes(x=status, label = group_letter, y=Inf),
            vjust=1.2,hjust=0.5, color = "black", inherit.aes = F, size=5, fontface="bold") +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0.1,0.2),0))+
  theme(
    axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs(y= expression(CO[2]~NEE~(mmol~d^-1~m^-2)),
       x=paste0("Conservation status"),
       fill=paste0("Status"))+
  facet_grid(rows="casepilot", scales="free")


ggsave(filename = "PAPERPLOTS_CO2_molar_boxplot_allCP_with-outliers.png.png",path = plots_path,device = "png",
       dpi = 400,
       units = "mm",
       height = 229,
       width = 90)

#Identity scale: boxplot(no outliers, all data) + sina (85central), FILTEREDcentral distribution of each status*cp: Still not good, tails still preclude from seing group differences in CU
central85_co2<- co2_4paper %>% 
  group_by(casepilot, status) %>% 
  mutate(upper=quantile(dailyflux, 0.925),
         lower=quantile(dailyflux, 0.075)) %>% 
  filter(between(dailyflux, lower, upper))

ggplot(co2_4paper, aes(x=status, y=dailyflux*1000))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(data=central85_co2, aes(col=status),position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(data=central85_co2, position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  #ADD EMMEANS and letters
  geom_point(data=co2_emmeans_4paper, aes(x=status,y=emmean_bt*1000),
             shape=23,size=3,fill = "black")+
  geom_text(data=co2_emmeans_4paper, aes(x=status, label = group_letter, y=Inf),
            vjust=1.2,hjust=0.5, color = "black", inherit.aes = F, size=5, fontface="bold") +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0.1,0.2),0))+
  theme(
    axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs(y= expression(CO[2]~NEE~(mmol~d^-1~m^-2)),
       x=paste0("Conservation status"),
       fill=paste0("Status"))+
  facet_grid(rows="casepilot", scales="free")








#Test using per-CP pseudolog transformation ISSUE (sign is not preserved due to standarize=T when performing the transformation, sligth shift) CANNOT present transformed data in plots. 


##CH4 --------
#Vertical plot CH4 molar units panels per case-pilot, independent scales, emmeans + letters. 

#Visualization issues: Each casepilot would benefit from its own scale (different compression ranges for logsign, interaction logtransformation & units). However, I think this would be very confusing. I could potentially show only boxplot (no outliers) + emmean, this "hides" weird distribution of data and makes the plots much clearer. 

ch4_4paper<- data4paper %>% filter(ghgspecies=="ch4")%>%   
  mutate(casepilot=factor(casepilot, levels = c("DU","RI","CA","VA","DA","CU")))
ch4_emmeans_4paper<-  emmeans_all %>% filter(comparison=="status",ghgspecies=="ch4")%>%   
  mutate(casepilot=factor(casepilot, levels = c("DU","RI","CA","VA","DA","CU")))


#Identity scale: boxplot without outliers works ok
ggplot(ch4_4paper, aes(x=status, y=dailyflux))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  #ADD EMMEANS and letters
  geom_point(data=ch4_emmeans_4paper, aes(x=status,y=emmean_bt),
             shape=23,size=3,fill = "black")+
  geom_text(data=ch4_emmeans_4paper, aes(x=status, label = group_letter, y=Inf),
            vjust=1.2,hjust=0.5, color = "black", inherit.aes = F, size=5, fontface="bold") +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0.1,0.2),0))+
  theme(
    axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs(y= expression(CH[4]~NEE~(mmol~d^-1~m^-2)),
       x=paste0("Conservation status"),
       fill=paste0("Status"))+
  facet_grid(rows="casepilot", scales="free")

ggsave(filename = "PAPERPLOTS_CH4_molar_boxplot_allCP.png.png",path = plots_path,device = "png",
       dpi = 400,
       units = "mm",
       height = 229,
       width = 90)

#Identity scale: boxplot(no outliers) + sina (all data), DOES not work, outliers force scale so that group differences are not clear. 
ggplot(ch4_4paper, aes(x=status, y=dailyflux))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(aes(col=status),position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  #ADD EMMEANS and letters
  geom_point(data=ch4_emmeans_4paper, aes(x=status,y=emmean_bt),
             shape=23,size=3,fill = "black")+
  geom_text(data=ch4_emmeans_4paper, aes(x=status, label = group_letter, y=Inf),
            vjust=1.2,hjust=0.5, color = "black", inherit.aes = F, size=5, fontface="bold") +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0.1,0.2),0))+
  theme(
    axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs(y= expression(CH[4]~NEE~(mmol~d^-1~m^-2)),
       x=paste0("Conservation status"),
       fill=paste0("Status"))+
  facet_grid(rows="casepilot", scales="free")

ggsave(filename = "PAPERPLOTS_CH4_molar_boxplot_allCP_with-outliers.png.png",path = plots_path,device = "png",
       dpi = 400,
       units = "mm",
       height = 229,
       width = 90)

#Identity scale: boxplot(no outliers, all data) + sina (85central), FILTEREDcentral distribution of each status*cp: Still not good, tails still preclude from seing group differences in CU
central85_ch4<- ch4_4paper %>% 
  group_by(casepilot, status) %>% 
  mutate(upper=quantile(dailyflux, 0.925),
         lower=quantile(dailyflux, 0.075)) %>% 
  filter(between(dailyflux, lower, upper))

ggplot(ch4_4paper, aes(x=status, y=dailyflux))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(data=central85_ch4, aes(col=status),position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(data=central85_ch4, position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  #ADD EMMEANS and letters
  geom_point(data=ch4_emmeans_4paper, aes(x=status,y=emmean_bt),
             shape=23,size=3,fill = "black")+
  geom_text(data=ch4_emmeans_4paper, aes(x=status, label = group_letter, y=Inf),
            vjust=1.2,hjust=0.5, color = "black", inherit.aes = F, size=5, fontface="bold") +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0.1,0.2),0))+
  theme(
    axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs(y= expression(CH[4]~NEE~(mmol~d^-1~m^-2)),
       x=paste0("Conservation status"),
       fill=paste0("Status"))+
  facet_grid(rows="casepilot", scales="free")



#Test using per-CP pseudolog transformation ISSUE (sign is not preserved due to standarize=T when performing the transformation, sligth shift) CANNOT present transformed data in plots. 



##GWP (co2+ch4) --------
#Vertical plot CO2 molar units panels per case-pilot, independent scales, emmeans + letters. 

#Visualization issues: Each casepilot would benefit from its own scale (different compression ranges for logsign, interaction logtransformation & units). However, I think this would be very confusing. I could potentially show only boxplot (no outliers) + emmean, this "hides" weird distribution of data and makes the plots much clearer. 

gwp_4paper<- data4paper %>% filter(ghgspecies=="gwp_co2andch4")%>%   
  mutate(casepilot=factor(casepilot, levels = c("DU","RI","CA","VA","DA","CU")))
gwp_emmeans_4paper<-  emmeans_all %>% filter(comparison=="status",ghgspecies=="gwp_co2andch4")%>%   
  mutate(casepilot=factor(casepilot, levels = c("DU","RI","CA","VA","DA","CU")))


#Identity scale: boxplot without outliers works ok
ggplot(gwp_4paper, aes(x=status, y=dailyflux))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  #ADD EMMEANS and letters
  geom_point(data=gwp_emmeans_4paper, aes(x=status,y=emmean_bt),
             shape=23,size=3,fill = "black")+
  geom_text(data=gwp_emmeans_4paper, aes(x=status, label = group_letter, y=Inf),
            vjust=1.2,hjust=0.5, color = "black", inherit.aes = F, size=5, fontface="bold") +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0.1,0.2),0))+
  theme(
    axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs(y= expression(CO[2]+CH[4]~GWP~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Conservation status"),
       fill=paste0("Status"))+
  facet_grid(rows="casepilot", scales="free")

ggsave(filename = "PAPERPLOTS_CO2+CH4_gwp_boxplot_allCP.png.png",path = plots_path,device = "png",
       dpi = 400,
       units = "mm",
       height = 229,
       width = 90)

#Identity scale: boxplot(no outliers) + sina (all data), DOES not work, outliers force scale so that group differences are not clear. 
ggplot(gwp_4paper, aes(x=status, y=dailyflux))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(aes(col=status),position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  #ADD EMMEANS and letters
  geom_point(data=gwp_emmeans_4paper, aes(x=status,y=emmean_bt),
             shape=23,size=3,fill = "black")+
  geom_text(data=gwp_emmeans_4paper, aes(x=status, label = group_letter, y=Inf),
            vjust=1.2,hjust=0.5, color = "black", inherit.aes = F, size=5, fontface="bold") +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0.1,0.2),0))+
  theme(
    axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs(y= expression(CO[2]+CH[4]~GWP~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Conservation status"),
       fill=paste0("Status"))+
  facet_grid(rows="casepilot", scales="free")

ggsave(filename = "PAPERPLOTS_CO2+CH4_gwp_boxplot_allCP_with-outliers.png.png",path = plots_path,device = "png",
       dpi = 400,
       units = "mm",
       height = 229,
       width = 90)

#Identity scale: boxplot(no outliers, all data) + sina (85central), FILTEREDcentral distribution of each status*cp: Still not good, tails still preclude from seing group differences in CU
central85_gwp<- gwp_4paper %>% 
  group_by(casepilot, status) %>% 
  mutate(upper=quantile(dailyflux, 0.925),
         lower=quantile(dailyflux, 0.075)) %>% 
  filter(between(dailyflux, lower, upper))

ggplot(gwp_4paper, aes(x=status, y=dailyflux))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(data=central85_gwp, aes(col=status),position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(data=central85_gwp, position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  #ADD EMMEANS and letters
  geom_point(data=gwp_emmeans_4paper, aes(x=status,y=emmean_bt),
             shape=23,size=3,fill = "black")+
  geom_text(data=gwp_emmeans_4paper, aes(x=status, label = group_letter, y=Inf),
            vjust=1.2,hjust=0.5, color = "black", inherit.aes = F, size=5, fontface="bold") +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0.1,0.2),0))+
  theme(
    axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs(y= expression(CO[2]+CH[4]~GWP~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Conservation status"),
       fill=paste0("Status"))+
  facet_grid(rows="casepilot", scales="free")





#____________________--------
#SEFS PLOTS-----------

#CHECK PLOT TUTORIAL-----------
# https://r-graph-gallery.com/web-violinplot-with-ggstatsplot.html


#3 panel plots, CO2, CH4, CO2+CH4
#Boxplots + violin plot (sina dots)
#Letters from GLMM post-hoc groups
#EMMean from GLMM post-hoc groups
#Log10-scales when needed.

library(ggforce)

log10sign<-  function(x) sign(x) * log10(abs(x) + 1)
#CO2:  needs log10scale for visual clarity. 
log10breaks_real<- c(-300,-100,-30,-10,-3,-1,0,1,3,10,30,100,300)
log10breaks_trans<- log10sign(log10breaks_real)


#1. Camargue-------

##___co2--------
#CO2:  needs log10scale for visual clarity. 


ca_co2<- data4paper %>% filter(casepilot=="CA", ghgspecies=="co2") %>% filter(!is.na(dailyflux))
ca_co2_emmean<- emmeans_all %>% filter(casepilot=="CA",comparison=="status",ghgspecies=="co2")

ca_co2_sefsplot<-
  ggplot(ca_co2, aes(x=status, y=log10sign(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(aes(col=status),position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  # geom_errorbar(data=ca_co2_emmean, aes(x=status,
  #                                       ymin=log10sign(lower.CL_bt),
  #                                       ymax=log10sign(upper.CL_bt)),inherit.aes = F,
  #               width=0.1)+ #ADD EMMEANS CI
  geom_point(data=ca_co2_emmean, aes(x=status,y=log10sign(emmean_bt)),shape=23,size=3,fill = "black")+ #ADD EMMEANS
  geom_text(data=ca_co2_emmean, aes(x=status, label = group_letter, y=log10sign(10)), vjust=0.9, color = "black", inherit.aes = F, size=6, fontface="bold") +
  scale_y_continuous(minor_breaks = NULL,
                     # limits = c(log10sign(-100),log10sign(100)),
                     breaks = log10breaks_trans,
                     labels = log10breaks_real
                     
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs(y= expression(CO[2]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Conservation status (n = ",nrow(ca_co2),")"),
       fill=paste0("Status"))

ca_co2_sefsplot

##___ch4 ------
# needs log10scale for visual clarity. 

ca_ch4<- ch4_formodels %>% filter(casepilot=="CA") %>% 
  # filter(dailyflux>-0.025) %>% #Remove 3 outliers  
  filter(!is.na(dailyflux))
ca_ch4_emmean<- emmeans_all %>% 
  filter(casepilot=="CA",comparison=="status",ghgspecies=="ch4")

ca_ch4_sefsplot<-
  ggplot(ca_ch4, aes(x=status, y=log10sign(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(aes(col=status),position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  # geom_errorbar(data=ca_ch4_emmean, aes(x=status,
  #                                       ymin=log10sign(lower.CL_bt),
  #                                       ymax=log10sign(upper.CL_bt)),inherit.aes = F,
  #               width=0.1)+ #ADD EMMEANS CI
  geom_point(data=ca_ch4_emmean, aes(x=status,y=log10sign(emmean_bt)),shape=23,size=3,fill = "black")+ #ADD EMMEANS
  geom_text(data=ca_ch4_emmean, aes(x=status, label = group_letter, y=log10sign(300)), vjust=0.9, color = "black", inherit.aes = F, size=6, fontface="bold") +
  scale_y_continuous(minor_breaks = NULL,
                     # limits = c(log10sign(-100),log10sign(100)),
                     breaks = log10breaks_trans,
                     labels = log10breaks_real
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs(y= expression(CH[4]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Conservation status (n = ",nrow(ca_ch4),")"),
       fill=paste0("Status"))

ca_ch4_sefsplot

##___co2+ch4-------
#: same trans than CO2. 

ca_co2andch4<- data4models %>% 
  filter(casepilot=="CA", ghgspecies=="gwp_co2andch4") %>% filter(!is.na(dailyflux))
ca_co2andch4_emmean<- emmeans_all %>% filter(casepilot=="CA",comparison=="status",ghgspecies=="gwp_co2andch4")

ca_co2andch4_sefsplot<- 
  ggplot(ca_co2andch4, aes(x=status, y=log10sign(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(aes(col=status),position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  # geom_errorbar(data=ca_co2andch4_emmean, aes(x=status,
  #                                       ymin=log10sign(lower.CL_bt),
  #                                       ymax=log10sign(upper.CL_bt)),inherit.aes = F,
  #               width=0.1)+ #ADD EMMEANS CI
  geom_point(data=ca_co2andch4_emmean, aes(x=status,y=log10sign(emmean_bt)),shape=23,size=3,fill = "black")+ #ADD EMMEANS
  geom_text(data=ca_co2andch4_emmean, aes(x=status, label = group_letter, y=log10sign(280)), vjust=0.9, color = "black", inherit.aes = F, size=6, fontface="bold") +
  scale_y_continuous(minor_breaks = NULL,
                     # limits = c(log10sign(-100),log10sign(100)),
                     breaks = log10breaks_trans,
                     labels = log10breaks_real
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs(y= expression(CO[2]+CH[4]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Conservation status (n = ",nrow(ca_co2andch4),")"),
       fill=paste0("Status"))

ca_co2andch4_sefsplot

ca_COMBO_sefsplot <- grid.arrange(ca_co2_sefsplot,ca_ch4_sefsplot,ca_co2andch4_sefsplot, ncol=3)
ggsave(plot = ca_COMBO_sefsplot, 
       filename = "SEFSPLOTS_CA_3ghg_status_violin.png",path = plots_path,device = "png",
       dpi = 400,
       units = "mm",
       height = 100, width = 350)


#2. Curonian-------

##___co2--------
#CO2:  needs log10scale for visual clarity. 


cu_co2<- data4models %>% filter(casepilot=="CU", ghgspecies=="co2") %>% filter(!is.na(dailyflux))
cu_co2_emmean<- emmeans_all %>% filter(casepilot=="CU",comparison=="status",ghgspecies=="co2")

cu_co2_sefsplot<-
  ggplot(cu_co2, aes(x=status, y=log10sign(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(aes(col=status),position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  # geom_errorbar(data=cu_co2_emmean, aes(x=status,
  #                                       ymin=log10sign(lower.CL_bt),
  #                                       ymax=log10sign(upper.CL_bt)),inherit.aes = F,
  #               width=0.1)+ #ADD EMMEANS CI
  geom_point(data=cu_co2_emmean, aes(x=status,y=log10sign(emmean_bt)),shape=23,size=3,fill = "black")+ #ADD EMMEANS
  geom_text(data=cu_co2_emmean, aes(x=status, label = group_letter, y=log10sign(1)), vjust=0.9, color = "black", inherit.aes = F, size=6, fontface="bold") +
  scale_y_continuous(minor_breaks = NULL,
                     # limits = c(log10sign(-100),log10sign(100)),
                     breaks = log10breaks_trans,
                     labels = log10breaks_real
                     
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs(y= expression(CO[2]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Conservation status (n = ",nrow(cu_co2),")"),
       fill=paste0("Status"))

cu_co2_sefsplot


##___ch4 ------
# needs log10scale for visual clarity. 

cu_ch4<- data4models %>% 
  filter(casepilot=="CU", ghgspecies=="ch4") %>% filter(!is.na(dailyflux))
cu_ch4_emmean<- emmeans_all %>% 
  filter(casepilot=="CU",comparison=="status",ghgspecies=="ch4")

cu_ch4_sefsplot<-
  ggplot(cu_ch4, aes(x=status, y=log10sign(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(aes(col=status),position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  geom_point(data=cu_ch4_emmean, aes(x=status,y=log10sign(emmean_bt)),shape=23,size=3,fill = "black")+ #ADD EMMEANS
  geom_text(data=cu_ch4_emmean, aes(x=status, label = group_letter, y=log10sign(300)), vjust=0.9, color = "black", inherit.aes = F, size=6, fontface="bold") +
  scale_y_continuous(minor_breaks = NULL,
                     # limits = c(log10sign(-100),log10sign(100)),
                     breaks = log10breaks_trans,
                     labels = log10breaks_real
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs(y= expression(CH[4]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Conservation status (n = ",nrow(cu_ch4),")"),
       fill=paste0("Status"))

cu_ch4_sefsplot

##___co2+ch4-------
#: same trans than CO2. 

cu_co2andch4<- data4models %>% 
  filter(casepilot=="CU", ghgspecies=="gwp_co2andch4") %>% filter(!is.na(dailyflux))
cu_co2andch4_emmean<- emmeans_all %>% filter(casepilot=="CU",comparison=="status",ghgspecies=="gwp_co2andch4")

cu_co2andch4_sefsplot<- 
  ggplot(cu_co2andch4, aes(x=status, y=log10sign(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(aes(col=status),position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  # geom_errorbar(data=cu_co2andch4_emmean, aes(x=status,
  #                                       ymin=log10sign(lower.CL_bt),
  #                                       ymax=log10sign(upper.CL_bt)),inherit.aes = F,
  #               width=0.1)+ #ADD EMMEANS CI
  geom_point(data=cu_co2andch4_emmean, aes(x=status,y=log10sign(emmean_bt)),shape=23,size=3,fill = "black")+ #ADD EMMEANS
  geom_text(data=cu_co2andch4_emmean, aes(x=status, label = group_letter, y=log10sign(250)), vjust=0.9, color = "black", inherit.aes = F, size=6, fontface="bold") +
  scale_y_continuous(minor_breaks = NULL,
                     # limits = c(log10sign(-100),log10sign(100)),
                     breaks = log10breaks_trans,
                     labels = log10breaks_real
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs(y= expression(CO[2]+CH[4]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Conservation status (n = ",nrow(cu_co2andch4),")"),
       fill=paste0("Status"))

cu_co2andch4_sefsplot

cu_COMBO_sefsplot <- grid.arrange(cu_co2_sefsplot,cu_ch4_sefsplot,cu_co2andch4_sefsplot, ncol=3)
ggsave(plot = cu_COMBO_sefsplot, 
       filename = "SEFSPLOTS_CU_3ghg_status_violin.png",path = plots_path,device = "png",
       dpi = 400,
       units = "mm",
       height = 100, width = 350)



#3. Danube-------

##___co2--------
#CO2:  needs log10scale for visual clarity. 


da_co2<- data4models %>% filter(casepilot=="DA", ghgspecies=="gwp_co2") %>% filter(!is.na(dailyflux))
da_co2_emmean<- emmeans_all %>% filter(casepilot=="DA",comparison=="status",ghgspecies=="gwp_co2")

da_co2_sefsplot<-
  ggplot(da_co2, aes(x=status, y=log10sign(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(aes(col=status),position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  # geom_errorbar(data=da_co2_emmean, aes(x=status,
  #                                       ymin=log10sign(lower.CL_bt),
  #                                       ymax=log10sign(upper.CL_bt)),inherit.aes = F,
  #               width=0.1)+ #ADD EMMEANS CI
  geom_point(data=da_co2_emmean, aes(x=status,y=log10sign(emmean_bt)),shape=23,size=3,fill = "black")+ #ADD EMMEANS
  geom_text(data=da_co2_emmean, aes(x=status, label = group_letter, y=log10sign(200)), vjust=0.9, color = "black", inherit.aes = F, size=6, fontface="bold") +
  scale_y_continuous(minor_breaks = NULL,
                     # limits = c(log10sign(-100),log10sign(100)),
                     breaks = log10breaks_trans,
                     labels = log10breaks_real
                     
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs(y= expression(CO[2]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Conservation status (n = ",nrow(da_co2),")"),
       fill=paste0("Status"))

da_co2_sefsplot


##___ch4 ------
# needs log10scale for visual clarity. 

da_ch4<- data4models %>% 
  filter(casepilot=="DA", ghgspecies=="gwp_ch4") %>% filter(!is.na(dailyflux))
da_ch4_emmean<- emmeans_all %>% 
  filter(casepilot=="DA",comparison=="status",ghgspecies=="gwp_ch4")

da_ch4_sefsplot<-
  ggplot(da_ch4, aes(x=status, y=log10sign(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(aes(col=status),position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  # geom_errorbar(data=da_ch4_emmean, aes(x=status,
  #                                       ymin=log10sign(lower.CL_bt),
  #                                       ymax=log10sign(upper.CL_bt)),inherit.aes = F,
  #               width=0.1)+ #ADD EMMEANS CI
  geom_point(data=da_ch4_emmean, aes(x=status,y=log10sign(emmean_bt)),shape=23,size=3,fill = "black")+ #ADD EMMEANS
  geom_text(data=da_ch4_emmean, aes(x=status, label = group_letter, y=log10sign(700)), vjust=0.9, color = "black", inherit.aes = F, size=6, fontface="bold") +
  scale_y_continuous(minor_breaks = NULL,
                     # limits = c(log10sign(-100),log10sign(100)),
                     breaks = log10breaks_trans,
                     labels = log10breaks_real
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs(y= expression(CH[4]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Conservation status (n = ",nrow(da_ch4),")"),
       fill=paste0("Status"))

da_ch4_sefsplot

##___co2+ch4-------
#: same trans than CO2. 

da_co2andch4<- data4models %>% 
  filter(casepilot=="DA", ghgspecies=="gwp_co2andch4") %>% filter(!is.na(dailyflux))
da_co2andch4_emmean<- emmeans_all %>% filter(casepilot=="DA",comparison=="status",ghgspecies=="gwp_co2andch4")

da_co2andch4_sefsplot<- 
  ggplot(da_co2andch4, aes(x=status, y=log10sign(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(aes(col=status),position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  # geom_errorbar(data=da_co2andch4_emmean, aes(x=status,
  #                                       ymin=log10sign(lower.CL_bt),
  #                                       ymax=log10sign(upper.CL_bt)),inherit.aes = F,
  #               width=0.1)+ #ADD EMMEANS CI
  geom_point(data=da_co2andch4_emmean, aes(x=status,y=log10sign(emmean_bt)),shape=23,size=3,fill = "black")+ #ADD EMMEANS
  geom_text(data=da_co2andch4_emmean, aes(x=status, label = group_letter, y=log10sign(600)), vjust=0.9, color = "black", inherit.aes = F, size=6, fontface="bold") +
  scale_y_continuous(minor_breaks = NULL,
                     # limits = c(log10sign(-100),log10sign(100)),
                     breaks = log10breaks_trans,
                     labels = log10breaks_real
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs(y= expression(CO[2]+CH[4]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Conservation status (n = ",nrow(da_co2andch4),")"),
       fill=paste0("Status"))

da_co2andch4_sefsplot

da_COMBO_sefsplot <- grid.arrange(da_co2_sefsplot,da_ch4_sefsplot,da_co2andch4_sefsplot, ncol=3)
ggsave(plot = da_COMBO_sefsplot, 
       filename = "SEFSPLOTS_DA_3ghg_status_violin.png",path = plots_path,device = "png",
       dpi = 400,
       units = "mm",
       height = 100, width = 350)



#4. Dutch delta-------
#WE might try non-transformed CH4 removing a couple of outliers. 

##___co2--------
#CO2:  needs log10scale for visual clarity. 


du_co2<- data4models %>% 
  filter(casepilot=="DU", ghgspecies=="gwp_co2") %>% filter(!is.na(dailyflux))
du_co2_emmean<- emmeans_all %>% 
  filter(casepilot=="DU",comparison=="status",ghgspecies=="gwp_co2")

du_co2_sefsplot<-
  ggplot(du_co2, aes(x=status, y=log10sign(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(aes(col=status),position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5,,scale = "width")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  geom_point(data=du_co2_emmean, aes(x=status,y=log10sign(emmean_bt)),shape=23,size=3,fill = "black")+ #ADD EMMEANS
  geom_text(data=du_co2_emmean, aes(x=status, label = group_letter, y=log10sign(100)), vjust=0.9, color = "black", inherit.aes = F,fontface = "bold",size=6 ) +
  scale_y_continuous(minor_breaks = NULL,
                     # limits = c(log10sign(-100),log10sign(100)),
                     breaks = log10breaks_trans,
                     labels = log10breaks_real
                     
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs(y= expression(CO[2]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Conservation status (n = ",nrow(du_co2),")"),
       fill=paste0("Status"))

du_co2_sefsplot


##___ch4 ------
# 

du_ch4<- data4models %>% 
  filter(casepilot=="DU", ghgspecies=="gwp_ch4") %>% filter(!is.na(dailyflux)) %>% 
  filter(dailyflux<0.1)# remove  outliers

du_ch4_emmean<- emmeans_all %>% 
  filter(casepilot=="DU",comparison=="status",ghgspecies=="gwp_ch4")

du_ch4_sefsplot<-
  ggplot(du_ch4, aes(x=status, y=log10sign(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(aes(col=status),position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5,,scale = "width")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  geom_point(data=du_ch4_emmean, aes(x=status,y=log10sign(emmean_bt)),shape=23,size=3,fill = "black")+ #ADD EMMEANS
  geom_text(data=du_ch4_emmean, aes(x=status, label = group_letter, y=log10sign(0.2)), vjust=0.9, color = "black", inherit.aes = F,fontface = "bold",size=6 ) +
  scale_y_continuous(minor_breaks = NULL,
                     limits = c(log10sign(-0.1),log10sign(0.3)),
                     breaks = log10sign(c(-0.1,-0.03,-0.01,0,0.01,0.03,0.1,0.3)),
                     labels =c(-0.1,-0.03,-0.01,0,0.01,0.03,0.1,0.3)
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs(y= expression(CH[4]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Conservation status (n = ",nrow(du_ch4),")"),
       fill=paste0("Status"))

du_ch4_sefsplot

##___co2+ch4-------
#: same trans than CO2. 

du_co2andch4<- data4models %>% 
  filter(casepilot=="DU", ghgspecies=="gwp_co2andch4") %>% filter(!is.na(dailyflux))
du_co2andch4_emmean<- emmeans_all %>%
  filter(casepilot=="DU",comparison=="status",ghgspecies=="gwp_co2andch4")

du_co2andch4_sefsplot<- 
  ggplot(du_co2andch4, aes(x=status, y=log10sign(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(aes(col=status),position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5,scale = "width")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  geom_point(data=du_co2andch4_emmean, aes(x=status,y=log10sign(emmean_bt)),shape=23,size=3,fill = "black")+ #ADD EMMEANS
  geom_text(data=du_co2andch4_emmean, aes(x=status, label = group_letter, y=log10sign(90)), vjust=0.9, color = "black", inherit.aes = F, size=6, fontface="bold") +
  scale_y_continuous(minor_breaks = NULL,
                     # limits = c(log10sign(-100),log10sign(100)),
                     breaks = log10breaks_trans,
                     labels = log10breaks_real
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs(y= expression(CO[2]+CH[4]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Conservation status (n = ",nrow(du_co2andch4),")"),
       fill=paste0("Status"))

du_co2andch4_sefsplot

du_COMBO_sefsplot <- grid.arrange(du_co2_sefsplot,du_ch4_sefsplot,du_co2andch4_sefsplot, ncol=3)
ggsave(plot = du_COMBO_sefsplot, 
       filename = "SEFSPLOTS_DU_3ghg_status_violin.png",path = plots_path,device = "png",
       dpi = 400,
       units = "mm",
       height = 100, width = 350)



#5. Ria de Aveiro-------

##___co2 (log)--------

ri_co2<- data4models %>% 
  filter(casepilot=="RI", ghgspecies=="gwp_co2") %>% filter(!is.na(dailyflux))%>% 
  filter(between(dailyflux, -7.5,2.5 )) #Remove 3 outliers

ri_co2_emmean<- emmeans_all %>% 
  filter(casepilot=="RI",comparison=="status",ghgspecies=="gwp_co2")

ri_co2_sefsplot<-
  ggplot(ri_co2, aes(x=status, y=log10sign(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(aes(col=status),position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  geom_point(data=ri_co2_emmean, aes(x=status,y=log10sign(emmean_bt)),shape=23,size=3,fill = "black")+ #ADD EMMEANS
  geom_text(data=ri_co2_emmean, aes(x=status, label = group_letter, y=log10sign(3)), vjust=0.9, color = "black", inherit.aes = F, fontface = "bold",size = 6) +
  scale_y_continuous(minor_breaks = NULL,
                     limits = c(log10sign(-10),log10sign(10)),
                     breaks = log10breaks_trans,
                     labels = log10breaks_real
                     
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs(y= expression(CO[2]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Conservation status (n = ",nrow(ri_co2),")"),
       fill=paste0("Status"))

ri_co2_sefsplot


##___ch4 (log) ------


ri_ch4<- data4models %>% 
  filter(casepilot=="RI", ghgspecies=="gwp_ch4") %>% filter(!is.na(dailyflux))%>% 
  filter(between(dailyflux, -0.025,0.025)) #CH4: limits = c(-0.025,0.025), (removes 8 outliers)

ri_ch4_emmean<- emmeans_all %>% 
  filter(casepilot=="RI",comparison=="status",ghgspecies=="gwp_ch4")

ri_ch4_sefsplot<-
  ggplot(ri_ch4, aes(x=status, y=log10sign(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(aes(col=status),position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  # geom_errorbar(data=ri_ch4_emmean, aes(x=status,
  #                                       ymin=log10sign(lower.CL_bt),
  #                                       ymax=log10sign(upper.CL_bt)),inherit.aes = F,
  #               width=0.1)+ #ADD EMMEANS CI
  geom_point(data=ri_ch4_emmean, aes(x=status,y=log10sign(emmean_bt)),shape=23,size=3,fill = "black")+ #ADD EMMEANS
  geom_text(data=ri_ch4_emmean, aes(x=status, label = group_letter, y=log10sign(0.03)), vjust=0.9, color = "black", inherit.aes = F, fontface = "bold",size = 6) +
  scale_y_continuous(minor_breaks = NULL,
                     limits = c(log10sign(-0.03),log10sign(0.03)),
                     breaks = log10sign(c(-0.03,-0.01,0,0.01,0.03))
                     ,
                     labels = c(-0.03,-0.01,0,0.01,0.03)
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs( 
    y= expression(CH[4]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
    x=paste0("Conservation status (n = ",nrow(ri_ch4),")"),
    fill=paste0("Status"))
ri_ch4_sefsplot

##___co2+ch4 (log)-------
#: same trans than CO2. 

ri_co2andch4<- data4models %>% 
  filter(casepilot=="RI", ghgspecies=="gwp_co2andch4") %>% filter(!is.na(dailyflux))%>% 
  filter(between(dailyflux, -7.5,2.5 )) #Remove 3 outliers

ri_co2andch4_emmean<- emmeans_all %>%
  filter(casepilot=="RI",comparison=="status",ghgspecies=="gwp_co2andch4")

ri_co2andch4_sefsplot<- 
  ggplot(ri_co2andch4, aes(x=status, y=log10sign(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(aes(col=status),position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  geom_point(data=ri_co2andch4_emmean, aes(x=status,y=log10sign(emmean_bt)),shape=23,size=3,fill = "black")+ #ADD EMMEANS
  geom_text(data=ri_co2andch4_emmean, aes(x=status, label = group_letter, y=log10sign(3)), vjust=0.9, color = "black", inherit.aes = F, size = 6,fontface = "bold") +
  scale_y_continuous(minor_breaks = NULL,
                     limits = c(log10sign(-10),log10sign(10)),
                     breaks = log10breaks_trans,
                     labels = log10breaks_real
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs(
    y= expression(CO[2]+CH[4]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
    x=paste0("Conservation status (n = ",nrow(ri_ch4),")"),
    fill=paste0("Status"))

ri_co2andch4_sefsplot

ri_COMBO_sefsplot <- grid.arrange(ri_co2_sefsplot,ri_ch4_sefsplot,ri_co2andch4_sefsplot, ncol=3)
ggsave(plot = ri_COMBO_sefsplot, 
       filename = "SEFSPLOTS_RI_3ghg_status_violin.png",path = plots_path,device = "png",
       dpi = 400,
       units = "mm",
       height = 100, width = 350)

##___co2 (linear)--------

ri_co2<- data4models %>% 
  filter(casepilot=="RI", ghgspecies=="gwp_co2") %>% filter(!is.na(dailyflux))%>% 
  filter(between(dailyflux, -7.5,2.5 )) #Remove 3 outliers

ri_co2_emmean<- emmeans_all %>% 
  filter(casepilot=="RI",comparison=="status",ghgspecies=="gwp_co2")

ri_co2_sefsplot_linear<-
  ggplot(ri_co2, aes(x=status, y=(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(aes(col=status),position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  geom_point(data=ri_co2_emmean, aes(x=status,y=(emmean_bt)),shape=23,size=3,fill = "black")+ #ADD EMMEANS
  geom_text(data=ri_co2_emmean, aes(x=status, label = group_letter, y=(3)), vjust=0.9, color = "black", inherit.aes = F, fontface = "bold",size = 6) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs(y= expression(CO[2]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Conservation status (n = ",nrow(ri_co2),")"),
       fill=paste0("Status"))

ri_co2_sefsplot_linear


##___ch4 (linear)------


ri_ch4<- data4models %>% 
  filter(casepilot=="RI", ghgspecies=="gwp_ch4") %>% filter(!is.na(dailyflux))%>% 
  filter(between(dailyflux, -0.025,0.025)) #CH4: limits = c(-0.025,0.025), (removes 8 outliers)

ri_ch4_emmean<- emmeans_all %>% 
  filter(casepilot=="RI",comparison=="status",ghgspecies=="gwp_ch4")

ri_ch4_sefsplot_linear<-
  ggplot(ri_ch4, aes(x=status, y=(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(aes(col=status),position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  geom_point(data=ri_ch4_emmean, aes(x=status,y=(emmean_bt)),shape=23,size=3,fill = "black")+ #ADD EMMEANS
  geom_text(data=ri_ch4_emmean, aes(x=status, label = group_letter, y=(0.03)), vjust=0.9, color = "black", inherit.aes = F, fontface = "bold",size = 6) +
  
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs( 
    y= expression(CH[4]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
    x=paste0("Conservation status (n = ",nrow(ri_ch4),")"),
    fill=paste0("Status"))
ri_ch4_sefsplot_linear

##___co2+ch4 (linear)-------
#: same trans than CO2. 

ri_co2andch4<- data4models %>% 
  filter(casepilot=="RI", ghgspecies=="gwp_co2andch4") %>% filter(!is.na(dailyflux))%>% 
  filter(between(dailyflux, -7.5,2.5 )) #Remove 3 outliers

ri_co2andch4_emmean<- emmeans_all %>%
  filter(casepilot=="RI",comparison=="status",ghgspecies=="gwp_co2andch4")

ri_co2andch4_sefsplot_linear<- 
  ggplot(ri_co2andch4, aes(x=status, y=(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(aes(col=status),position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  geom_point(data=ri_co2andch4_emmean, aes(x=status,y=(emmean_bt)),shape=23,size=3,fill = "black")+ #ADD EMMEANS
  geom_text(data=ri_co2andch4_emmean, aes(x=status, label = group_letter, y=(3)), vjust=0.9, color = "black", inherit.aes = F, size = 6,fontface = "bold") +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs(
    y= expression(CO[2]+CH[4]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
    x=paste0("Conservation status (n = ",nrow(ri_ch4),")"),
    fill=paste0("Status"))

ri_co2andch4_sefsplot_linear

ri_COMBO_sefsplot_linear <- grid.arrange(ri_co2_sefsplot_linear,ri_ch4_sefsplot_linear,ri_co2andch4_sefsplot_linear, ncol=3)
ggsave(plot = ri_COMBO_sefsplot_linear, 
       filename = "SEFSPLOTS_RI_linear_3ghg_status_violin.png",path = plots_path,device = "png",
       dpi = 400,
       units = "mm",
       height = 100, width = 350)




#6. Valencia-------

##___co2--------
#CO2:  needs log10scale for visual clarity. 


va_co2<- data4models %>% filter(casepilot=="VA", ghgspecies=="gwp_co2") %>% filter(!is.na(dailyflux))
va_co2_emmean<- emmeans_all %>% filter(casepilot=="VA",comparison=="status",ghgspecies=="gwp_co2")

va_co2_sefsplot<-
  ggplot(va_co2, aes(x=status, y=log10sign(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(aes(col=status),position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  # geom_errorbar(data=va_co2_emmean, aes(x=status,
  #                                       ymin=log10sign(lower.CL_bt),
  #                                       ymax=log10sign(upper.CL_bt)),inherit.aes = F,
  #               width=0.1)+ #ADD EMMEANS CI
  geom_point(data=va_co2_emmean, aes(x=status,y=log10sign(emmean_bt)),shape=23,size=3,fill = "black")+ #ADD EMMEANS
  geom_text(data=va_co2_emmean, aes(x=status, label = group_letter, y=log10sign(120)), vjust=0.9, color = "black", inherit.aes = F, nudge_x = 0.1, size = 6,fontface = "bold") +
  scale_y_continuous(minor_breaks = NULL,
                     # limits = c(log10sign(-100),log10sign(100)),
                     breaks = log10breaks_trans,
                     labels = log10breaks_real
                     
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs(y= expression(CO[2]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Conservation status (n = ",nrow(va_co2),")"),
       fill=paste0("Status"))

va_co2_sefsplot


##___ch4 ------
# needs log10scale for visual clarity. 

va_ch4<- data4models %>% 
  filter(casepilot=="VA", ghgspecies=="gwp_ch4") %>% filter(!is.na(dailyflux))
va_ch4_emmean<- emmeans_all %>% 
  filter(casepilot=="VA",comparison=="status",ghgspecies=="gwp_ch4")

va_ch4_sefsplot<-
  ggplot(va_ch4, aes(x=status, y=log10sign(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(aes(col=status),position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  # geom_errorbar(data=va_ch4_emmean, aes(x=status,
  #                                       ymin=log10sign(lower.CL_bt),
  #                                       ymax=log10sign(upper.CL_bt)),inherit.aes = F,
  #               width=0.1)+ #ADD EMMEANS CI
  geom_point(data=va_ch4_emmean, aes(x=status,y=log10sign(emmean_bt)),shape=23,size=3,fill = "black")+ #ADD EMMEANS
  geom_text(data=va_ch4_emmean, aes(x=status, label = group_letter, y=log10sign(30)), vjust=0.9, color = "black", inherit.aes = F, nudge_x = 0.1, size = 6,fontface = "bold") +
  scale_y_continuous(minor_breaks = NULL,
                     # limits = c(log10sign(-100),log10sign(100)),
                     breaks = log10breaks_trans,
                     labels = log10breaks_real
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs(y= expression(CH[4]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Conservation status (n = ",nrow(va_ch4),")"),
       fill=paste0("Status"))

va_ch4_sefsplot

##___co2+ch4-------
#: same trans than CO2. 

va_co2andch4<- data4models %>% 
  filter(casepilot=="VA", ghgspecies=="gwp_co2andch4") %>% filter(!is.na(dailyflux))
va_co2andch4_emmean<- emmeans_all %>% filter(casepilot=="VA",comparison=="status",ghgspecies=="gwp_co2andch4")

va_co2andch4_sefsplot<- 
  ggplot(va_co2andch4, aes(x=status, y=log10sign(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(aes(col=status),position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, aes(fill=status),size = 0.7)+
  geom_point(data=va_co2andch4_emmean, aes(x=status,y=log10sign(emmean_bt)),shape=23,size=3,fill = "black")+ #ADD EMMEANS
  geom_text(data=va_co2andch4_emmean, aes(x=status, label = group_letter, y=log10sign(120)), vjust=0.9, color = "black", inherit.aes = F, nudge_x = 0.1, size = 6,fontface = "bold") +
  scale_y_continuous(minor_breaks = NULL,
                     # limits = c(log10sign(-100),log10sign(100)),
                     breaks = log10breaks_trans,
                     labels = log10breaks_real
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"))+
  guides(color = "none", fill="none")+
  labs(y= expression(CO[2]+CH[4]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Conservation status (n = ",nrow(va_co2andch4),")"),
       fill=paste0("Status"))

va_co2andch4_sefsplot

va_COMBO_sefsplot <- grid.arrange(va_co2_sefsplot,va_ch4_sefsplot,va_co2andch4_sefsplot, ncol=3)
ggsave(plot = va_COMBO_sefsplot, 
       filename = "SEFSPLOTS_VA_3ghg_status_violin.png",path = plots_path,device = "png",
       dpi = 400,
       units = "mm",
       height = 100, width = 350)



#Preserved across CP--------


#___Boxplot + violin -------
data4models %>% 
  filter(status=="Preserved", ghgspecies%in%c("gwp_co2")) %>% 
  mutate(casepilot=factor(casepilot, levels = c("DU","RI","CA","VA","DA","CU"))) %>% 
  ggplot(aes(x=casepilot, y=log10sign(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(col ="#00BA38",position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, fill="#00BA38",size = 0.7)+
  scale_y_continuous(minor_breaks = NULL,
                     # limits = c(log10sign(-100),log10sign(100)),
                     breaks = log10breaks_trans,
                     labels = log10breaks_real
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
  guides(color = "none", fill="none")+
  labs(title=expression(CO[2]~GWP) ,y= expression(CO[2]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Casepilot"),
       fill=paste0("Status"))


ggsave( 
  filename = "Preserved only_CO2_gwp.png",path = plots_path,device = "png",
  dpi = 400,
  units = "mm",
  height = 100, width = 125)




data4models %>% 
  filter(status=="Preserved", ghgspecies%in%c("gwp_ch4")) %>% 
  mutate(casepilot=factor(casepilot, levels = c("DU","RI","CA","VA","DA","CU"))) %>% 
  ggplot(aes(x=casepilot, y=log10sign(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(col ="#00BA38",position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, fill="#00BA38",size = 0.7)+
  scale_y_continuous(minor_breaks = NULL,
                     # limits = c(log10sign(-100),log10sign(100)),
                     breaks = log10breaks_trans,
                     labels = log10breaks_real
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
  guides(color = "none", fill="none")+
  labs(title=expression(CH[4]~GWP) ,y= expression(CH[4]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Casepilot"),
       fill=paste0("Status"))


ggsave( 
  filename = "Preserved only_CH4_gwp.png",path = plots_path,device = "png",
  dpi = 400,
  units = "mm",
  height = 100, width = 125)


#___Boxplot only-------

data4models %>% 
  filter(status=="Preserved", ghgspecies%in%c("gwp_co2")) %>% 
  mutate(casepilot=factor(casepilot, levels = c("DU","RI","CA","VA","DA","CU"))) %>% 
  ggplot(aes(x=casepilot, y=log10sign(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  # geom_violin(col ="#00BA38",position = position_dodge(0.9),linewidth = 1,scale = "width")+
  # geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, fill="#00BA38",size = 0.7)+
  scale_y_continuous(minor_breaks = NULL,
                     # limits = c(log10sign(-100),log10sign(100)),
                     breaks = log10breaks_trans,
                     labels = log10breaks_real
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
  guides(color = "none", fill="none")+
  labs(title=expression(CO[2]~GWP) ,y= expression(CO[2]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Casepilot"),
       fill=paste0("Status"))


ggsave( 
  filename = "Preserved only_just_boxplot_CO2_gwp.png",path = plots_path,device = "png",
  dpi = 400,
  units = "mm",
  height = 100, width = 125)




data4models %>% 
  filter(status=="Preserved", ghgspecies%in%c("gwp_ch4")) %>% 
  mutate(casepilot=factor(casepilot, levels = c("DU","RI","CA","VA","DA","CU"))) %>% 
  ggplot(aes(x=casepilot, y=log10sign(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  # geom_violin(col ="#00BA38",position = position_dodge(0.9),linewidth = 1,scale = "width")+
  # geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, fill="#00BA38",size = 0.7)+
  scale_y_continuous(minor_breaks = NULL,
                     # limits = c(log10sign(-100),log10sign(100)),
                     breaks = log10breaks_trans,
                     labels = log10breaks_real
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
  guides(color = "none", fill="none")+
  labs(title=expression(CH[4]~GWP) ,y= expression(CH[4]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Casepilot"),
       fill=paste0("Status"))


ggsave( 
  filename = "Preserved only_just boxplot_CH4_gwp.png",path = plots_path,device = "png",
  dpi = 400,
  units = "mm",
  height = 100, width = 125)




#____________________--------


#____________________--------
#3. Plots CO2 Status effect-----
#TO-ADAPT: 
#Emmeans are already back-transformed into dailyflux co2 units (same as ChamberData4paper), fix scale labels. 
#Path to save is not correct (decide if/where to save these plots)

co2_casepilot_comparisons <- read.csv(paste0(paper_path,"CO2_emmeans-posthoc_chambermodels.csv"))


#Plot casepilot-per-casepilot, emmeans +- SE/CL (in molar scale)
#Camargue: 
plot_co2_status_CA <- co2_casepilot_comparisons %>%
  filter(casepilot == "CA", comparison == "status") %>%
  ggplot(., aes(x = status, y = emmean_bt, label = group_letter)) +
  geom_errorbar(aes(ymin = emmean_bt - SE_bt, ymax = emmean_bt + SE_bt),
                linewidth = 0.2, width = 0.2) +
  geom_point(shape = 15, size = 4, aes(col=status)) +
  geom_text(nudge_x = 0.2, color = "black") +
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw() +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  labs(title = expression(Status~effect~on~net~CO[2]~flux),
       subtitle = paste0("Camargue"), 
       y=expression(EMM~(g~CO[2~eq]~d^-1~m^-2)),
       x="Conservation status",
       col=paste0("Status"))

plot_co2_status_CA
ggsave(plot=plot_co2_status_CA, filename="CO2_status_plot_CA.png",
       device = "png",path = plots_path,dpi = 400,
       width = 115, height = 100,units = "mm")


#Curonian Lagoon: 
plot_co2_status_CU<- co2_casepilot_comparisons %>%
  filter(casepilot == "CU", comparison == "status") %>%
  ggplot(., aes(x = status, y = emmean_bt, label = group_letter)) +
  geom_errorbar(aes(ymin = emmean_bt - SE_bt, ymax = emmean_bt + SE_bt),
                linewidth = 0.2, width = 0.2) +
  geom_point(shape = 15, size = 4, aes(col=status)) +
  geom_text(nudge_x = 0.2, color = "black") +
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw() +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  labs(title = expression(Status~effect~on~net~CO[2]~flux),
       subtitle = paste0("Curonian lagoon"), 
       y=expression(EMM~(g~CO[2~eq]~d^-1~m^-2)),
       x="Conservation status",
       col=paste0("Status"))

plot_co2_status_CU
ggsave(plot=plot_co2_status_CU, filename="CO2_status_plot_CU.png",
       device = "png",path = plots_path,dpi = 400,
       width = 115, height = 100,units = "mm")

#Danube delta: 
plot_co2_status_DA<- co2_casepilot_comparisons %>%
  filter(casepilot == "DA", comparison == "status") %>%
  ggplot(., aes(x = status, y = emmean_bt, label = group_letter)) +
  geom_errorbar(aes(ymin = emmean_bt - SE_bt, ymax = emmean_bt + SE_bt),
                linewidth = 0.2, width = 0.2) +
  geom_point(shape = 15, size = 4, aes(col=status)) +
  geom_text(nudge_x = 0.2, color = "black") +
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw() +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  labs(title = expression(Status~effect~on~net~CO[2]~flux),
       subtitle = paste0("Danube delta"),
       y=expression(EMM~(g~CO[2~eq]~d^-1~m^-2)),
       x="Conservation status",
       col=paste0("Status"))

plot_co2_status_DA
ggsave(plot=plot_co2_status_DA, filename="CO2_status_plot_DA.png",
       device = "png",path = plots_path,dpi = 400,
       width = 115, height = 100,units = "mm")

#Dutch delta: 
plot_co2_status_DU<- co2_casepilot_comparisons %>%
  filter(casepilot == "DU", comparison == "status") %>%
  ggplot(., aes(x = status, y = emmean_bt, label = group_letter)) +
  geom_errorbar(aes(ymin = emmean_bt - SE_bt, ymax = emmean_bt + SE_bt),
                linewidth = 0.2, width = 0.2) +
  geom_point(shape = 15, size = 4, aes(col=status)) +
  geom_text(nudge_x = 0.2, color = "black") +
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw() +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  labs(title = expression(Status~effect~on~net~CO[2]~flux),
       subtitle = paste0("Dutch delta"),
       y=expression(EMM~(g~CO[2~eq]~d^-1~m^-2)),
       x="Conservation status",
       col=paste0("Status"))

plot_co2_status_DU
ggsave(plot=plot_co2_status_DU, filename="CO2_status_plot_DU.png",
       device = "png",path = plots_path,dpi = 400,
       width = 115, height = 100,units = "mm")


#Ria d'Aveiro: 
plot_co2_status_RI<- co2_casepilot_comparisons %>%
  filter(casepilot == "RI", comparison == "status") %>%
  ggplot(., aes(x = status, y = emmean_bt, label = group_letter)) +
  geom_errorbar(aes(ymin = emmean_bt - SE_bt, ymax = emmean_bt + SE_bt),
                linewidth = 0.2, width = 0.2) +
  geom_point(shape = 15, size = 4, aes(col=status)) +
  geom_text(nudge_x = 0.2, color = "black") +
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw() +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  labs(title = expression(Status~effect~on~net~CO[2]~flux),
       subtitle = paste0("Ria d'Aveiro"), 
       y=expression(EMM~(g~CO[2~eq]~d^-1~m^-2)),
       x="Conservation status",
       col=paste0("Status"))

plot_co2_status_RI
ggsave(plot=plot_co2_status_RI, filename="CO2_status_plot_RI.png",
       device = "png",path = plots_path,dpi = 400,
       width = 115, height = 100,units = "mm")

#Valencian wetland: 
plot_co2_status_VA<- co2_casepilot_comparisons %>%
  filter(casepilot == "VA", comparison == "status") %>%
  ggplot(., aes(x = status, y = emmean_bt, label = group_letter)) +
  geom_errorbar(aes(ymin = emmean_bt - SE_bt, ymax = emmean_bt + SE_bt),
                linewidth = 0.2, width = 0.2) +
  geom_point(shape = 15, size = 4, aes(col=status)) +
  geom_text(nudge_x = 0.2, color = "black") +
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw() +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  labs(title = expression(Status~effect~on~net~CO[2]~flux),
       subtitle = paste0("Valencian wetland"), 
       y=expression(EMM~(g~CO[2~eq]~d^-1~m^-2)),
       x="Conservation status",
       col=paste0("Status"))

plot_co2_status_VA
ggsave(plot=plot_co2_status_VA, filename="CO2_status_plot_VA.png",
       device = "png",path = plots_path,dpi = 400,
       width = 115, height = 100,units = "mm")


##Extra plots (to-do)------

#Plots Seasonal effect
#To do, little interest for overview ppt (seasonal effect is averaged across the three status)

#Example:
#Camargue: 
plot_co2_season_CA<- co2_casepilot_comparisons %>%
  filter(casepilot == "CA", comparison == "season") %>%
  ggplot(., aes(x = season, y = emmean_bt, label = group_letter)) +
  geom_errorbar(aes(ymin = emmean_bt - SE_bt, ymax = emmean_bt + SE_bt),
                linewidth = 0.2, width = 0.2) +
  geom_point(shape = 15, size = 4, aes(col=season)) +
  geom_text(nudge_x = 0.2, color = "black") +
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw() +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  labs(title = expression(Seasonal~effect~on~net~CO[2]~flux),
       subtitle = paste0("Camargue"), 
       caption = paste0("Seasonal effect is averaged across the three conservation status.\nBoxes indicate the EM mean, error bars indicate SE,\n",              "Letters indicate different groups"),
       y=expression(EMM~(g~CO[2~eq]~d^-1~m^-2)),
       x="Sampling campaign",
       col=paste0("Sampling"))

plot_co2_season_CA


##Plots Interaction effect
#To do, unclear interest for overview ppt (maybe interesting for some casepilots)

#Example:
#Camargue: interaction plot
plot_co2_interaction_CA<- co2_casepilot_comparisons %>%
  filter(casepilot == "CA", comparison == "interaction") %>%
  #ggplot inside brackets to be able to access casepilot name with ggtitle or define dodge for all
  { 
    #Define dodge width for separating status.
    dodge <- position_dodge(width = 0.6)
    
    ggplot(., aes(x = season, y = emmean_bt, label = group_letter, group=status)) +
      geom_errorbar(aes(ymin = emmean_bt - SE_bt, ymax = emmean_bt + SE_bt),
                    linewidth = 0.2, width = 0.2, position = dodge) +
      geom_point(shape = 15, size = 4, aes(col=status), position = dodge) +
      geom_text(vjust = -2, color = "black", position = dodge) +
      geom_hline(yintercept=0, linetype="dashed")+
      theme_bw() +
      theme(axis.title = element_text(face = "bold"),
            axis.text = element_text(face = "bold"),
            plot.caption = element_text(hjust = 0)) +
      labs(title = expression(Seasonal~effect~on~net~CO[2]~flux),
           subtitle = paste0("Camargue"), 
           caption = paste0("Boxes indicate the EM mean, error bars indicate SE,\n",              "Letters indicate different groups"),
           y=expression(EMM~(mol~CO[2]~d^-1~m^-2)),
           x="Sampling campaign",
           col=paste0("Sampling"))
  }

plot_co2_interaction_CA


#Plot casepilot-per-casepilot, emmeans +- SE (in gwp scale)

#3. Status CH4 effect plots-----

#Camargue: 
plot_ch4_status_CA<-ch4_casepilot_comparisons %>%
  filter(casepilot == "CA", comparison == "status") %>%
  ggplot(., aes(x = status, y = emmean_bt, label = group_letter)) +
  geom_errorbar(aes(ymin = emmean_bt - SE_bt, ymax = emmean_bt + SE_bt),
                linewidth = 0.2, width = 0.2) +
  geom_point(shape = 15, size = 4, aes(col=status)) +
  geom_text(nudge_x = 0.2, color = "black") +
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw() +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  labs(title = expression(Status~effect~on~net~CH[4]~flux),
       subtitle = paste0("Camargue"), 
       y=expression(EMM~(mmol~CH[4]~d^-1~m^-2)),
       x="Conservation status",
       col=paste0("Status"))

plot_ch4_status_CA
ggsave(plot=plot_ch4_status_CA, filename="CH4_status_plot_CA.png",
       device = "png",path = plots_path,dpi = 400,
       width = 115, height = 100,units = "mm")


#Curonian Lagoon: 
plot_ch4_status_CU<- ch4_casepilot_comparisons %>%
  filter(casepilot == "CU", comparison == "status") %>%
  ggplot(., aes(x = status, y = emmean_bt, label = group_letter)) +
  geom_errorbar(aes(ymin = emmean_bt - SE_bt, ymax = emmean_bt + SE_bt),
                linewidth = 0.2, width = 0.2) +
  geom_point(shape = 15, size = 4, aes(col=status)) +
  geom_text(nudge_x = 0.2, color = "black") +
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw() +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  labs(title = expression(Status~effect~on~net~CH[4]~flux),
       subtitle = paste0("Curonian lagoon"), 
       y=expression(EMM~(mmol~CH[4]~d^-1~m^-2)),
       x="Conservation status",
       col=paste0("Status"))

plot_ch4_status_CU
ggsave(plot=plot_ch4_status_CU, filename="CH4_status_plot_CU.png",
       device = "png",path = plots_path,dpi = 400,
       width = 115, height = 100,units = "mm")

#Danube delta: 
plot_ch4_status_DA<- ch4_casepilot_comparisons %>%
  filter(casepilot == "DA", comparison == "status") %>%
  ggplot(., aes(x = status, y = emmean_bt, label = group_letter)) +
  geom_errorbar(aes(ymin = emmean_bt - SE_bt, ymax = emmean_bt + SE_bt),
                linewidth = 0.2, width = 0.2) +
  geom_point(shape = 15, size = 4, aes(col=status)) +
  geom_text(nudge_x = 0.2, color = "black") +
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw() +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  labs(title = expression(Status~effect~on~net~CH[4]~flux),
       subtitle = paste0("Danube delta"), 
       y=expression(EMM~(mmol~CH[4]~d^-1~m^-2)),
       x="Conservation status",
       col=paste0("Status"))

plot_ch4_status_DA

ggsave(plot=plot_ch4_status_DA, filename="CH4_status_plot_DA.png",
       device = "png",path = plots_path,dpi = 400,
       width = 115, height = 100,units = "mm")

#Dutch delta: 
plot_ch4_status_DU<- ch4_casepilot_comparisons %>%
  filter(casepilot == "DU", comparison == "status") %>%
  ggplot(., aes(x = status, y = emmean_bt, label = group_letter)) +
  geom_errorbar(aes(ymin = emmean_bt - SE_bt, ymax = emmean_bt + SE_bt),
                linewidth = 0.2, width = 0.2) +
  geom_point(shape = 15, size = 4, aes(col=status)) +
  geom_text(nudge_x = 0.2, color = "black") +
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw() +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  labs(title = expression(Status~effect~on~net~CH[4]~flux),
       subtitle = paste0("Dutch delta"),
       y=expression(EMM~(mmol~CH[4]~d^-1~m^-2)),
       x="Conservation status",
       col=paste0("Status"))

plot_ch4_status_DU

ggsave(plot=plot_ch4_status_DU, filename="CH4_status_plot_DU.png",
       device = "png",path = plots_path,dpi = 400,
       width = 115, height = 100,units = "mm")


#Ria d'Aveiro: 
plot_ch4_status_RI<-ch4_casepilot_comparisons %>%
  filter(casepilot == "RI", comparison == "status") %>%
  ggplot(., aes(x = status, y = emmean_bt, label = group_letter)) +
  geom_errorbar(aes(ymin = emmean_bt - SE_bt, ymax = emmean_bt + SE_bt),
                linewidth = 0.2, width = 0.2) +
  geom_point(shape = 15, size = 4, aes(col=status)) +
  geom_text(nudge_x = 0.2, color = "black") +
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw() +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  labs(title = expression(Status~effect~on~net~CH[4]~flux),
       subtitle = paste0("Ria d'Aveiro"), 
       y=expression(EMM~(mmol~CH[4]~d^-1~m^-2)),
       x="Conservation status",
       col=paste0("Status"))

plot_ch4_status_RI
ggsave(plot=plot_ch4_status_RI, filename="CH4_status_plot_RI.png",
       device = "png",path = plots_path,dpi = 400,
       width = 115, height = 100,units = "mm")

#Valencian wetland: 
plot_ch4_status_VA<-ch4_casepilot_comparisons %>%
  filter(casepilot == "VA", comparison == "status") %>%
  ggplot(., aes(x = status, y = emmean_bt, label = group_letter)) +
  geom_errorbar(aes(ymin = emmean_bt - SE_bt, ymax = emmean_bt + SE_bt),
                linewidth = 0.2, width = 0.2) +
  geom_point(shape = 15, size = 4, aes(col=status)) +
  geom_text(nudge_x = 0.2, color = "black") +
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw() +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  labs(title = expression(Status~effect~on~net~CH[4]~flux),
       subtitle = paste0("Valencian wetland"), 
       y=expression(EMM~(mmol~CH[4]~d^-1~m^-2)),
       x="Conservation status",
       col=paste0("Status"))

plot_ch4_status_VA
ggsave(plot=plot_ch4_status_VA, filename="CH4_status_plot_VA.png",
       device = "png",path = plots_path,dpi = 400,
       width = 115, height = 100,units = "mm")


# grid.arrange(plot_ch4_status_CA,plot_ch4_status_CU,plot_ch4_status_DA, plot_ch4_status_DU, plot_ch4_status_RI, plot_ch4_status_VA, ncol=3)
# grid.arrange(plot_ch4_status_CA, plot_ch4_status_CU,p3, ncol = 3)



##Extra plots (to-do)-----

#Seasonal plots
#To do, not too much interest for overview ppt (seasonal effect is averaged across the three status)

#Example: CA
#Camargue: seasonal plot
plot_ch4_season_CA<-ch4_casepilot_comparisons %>%
  filter(casepilot == "CA", comparison == "season") %>%
  #ggplot inside brackets to be able to access casepilot name with ggtitle
  { 
    ggplot(., aes(x = season, y = emmean_gwp, label = group_letter)) +
      geom_errorbar(aes(ymin = emmean_gwp - SE_gwp, ymax = emmean_gwp + SE_gwp),
                    linewidth = 0.2, width = 0.2) +
      geom_point(shape = 15, size = 4, aes(col=season)) +
      geom_text(nudge_y = 0.4, color = "black") +
      geom_hline(yintercept=0, linetype="dashed")+
      theme_bw() +
      theme(axis.title = element_text(face = "bold"),
            axis.text = element_text(face = "bold"),
            plot.caption = element_text(hjust = 0)) +
      labs(title = expression(Seasonal~effect~on~net~CH[4]~flux),
           subtitle = paste0("Camargue"), 
           caption = paste0("Seasonal effect is averaged across the three conservation status.\nBoxes indicate the EM mean, error bars indicate SE,\n",              "Letters indicate different groups"),
           y=expression(EMM~(g~CO[2~eq]~d^-1~m^-2)),
           x="Sampling campaign",
           col=paste0("Sampling"))
  }

plot_ch4_season_CA


##Interaction effect plots
#To do, unclear interest for overview ppt (maybe interesting for some casepilots)

#Example: CA
#Camargue: interaction plot
plot_ch4_interaction_CA<-ch4_casepilot_comparisons %>%
  filter(casepilot == "CA", comparison == "interaction") %>%
  #ggplot inside brackets to be able to access casepilot name with ggtitle
  { 
    #Define dodge width for separating status.
    dodge <- position_dodge(width = 0.6)
    
    ggplot(., aes(x = season, y = emmean_gwp, label = group_letter, group=status)) +
      geom_errorbar(aes(ymin = emmean_gwp - SE_gwp, ymax = emmean_gwp + SE_gwp),
                    linewidth = 0.2, width = 0.2, position = dodge) +
      geom_point(shape = 15, size = 4, aes(col=status), position = dodge) +
      geom_text(vjust = -2, color = "black", position = dodge) +
      geom_hline(yintercept=0, linetype="dashed")+
      theme_bw() +
      theme(axis.title = element_text(face = "bold"),
            axis.text = element_text(face = "bold"),
            plot.caption = element_text(hjust = 0)) +
      labs(title = expression(Seasonal~effect~on~net~CO[2]~flux),
           subtitle = paste0("Camargue"), 
           caption = paste0("Boxes indicate the EM mean, error bars indicate SE,\n",              "Letters indicate different groups"),
           y=expression(EMM~(g~CO[2~eq]~d^-1~m^-2)),
           x="Sampling campaign",
           col=paste0("Sampling"))
  }

plot_ch4_interaction_CA


#3. Plots GWP Status effect-----
#Plot casepilot-per-casepilot, emmeans +- SE/CL (in gwp scale)
#Camargue: 
plot_co2andch4_status_CA <- co2andch4_casepilot_comparisons %>%
  filter(casepilot == "CA", comparison == "status") %>%
  ggplot(., aes(x = status, y = emmean_bt, label = group_letter)) +
  geom_errorbar(aes(ymin = emmean_bt - SE_bt, ymax = emmean_bt + SE_bt),
                linewidth = 0.2, width = 0.2) +
  geom_point(shape = 15, size = 4, aes(col=status)) +
  geom_text(nudge_x = 0.2, color = "black") +
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw() +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  labs(title = expression(Status~effect~on~NEE~(CO[2]+CH[4])),
       subtitle = paste0("Camargue"), 
       y=expression(EMM~(g~CO[2~eq]~d^-1~m^-2)),
       x="Conservation status",
       col=paste0("Status"))

plot_co2andch4_status_CA
ggsave(plot=plot_co2andch4_status_CA, filename="CO2plusCH4_status_plot_CA.png",
       device = "png",path = plots_path,dpi = 400,
       width = 115, height = 100,units = "mm")


#Curonian Lagoon: 
plot_co2andch4_status_CU<- co2andch4_casepilot_comparisons %>%
  filter(casepilot == "CU", comparison == "status") %>%
  ggplot(., aes(x = status, y = emmean_bt, label = group_letter)) +
  geom_errorbar(aes(ymin = emmean_bt - SE_bt, ymax = emmean_bt + SE_bt),
                linewidth = 0.2, width = 0.2) +
  geom_point(shape = 15, size = 4, aes(col=status)) +
  geom_text(nudge_x = 0.2, color = "black") +
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw() +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  labs(title = expression(Status~effect~on~NEE~(CO[2]+CH[4])),
       subtitle = paste0("Curonian lagoon"), 
       y=expression(EMM~(g~CO[2~eq]~d^-1~m^-2)),
       x="Conservation status",
       col=paste0("Status"))

plot_co2andch4_status_CU
ggsave(plot=plot_co2andch4_status_CU, filename="CO2plusCH4_status_plot_CU.png",
       device = "png",path = plots_path,dpi = 400,
       width = 115, height = 100,units = "mm")

#Danube delta: 
plot_co2andch4_status_DA<- co2andch4_casepilot_comparisons %>%
  filter(casepilot == "DA", comparison == "status") %>%
  ggplot(., aes(x = status, y = emmean_bt, label = group_letter)) +
  geom_errorbar(aes(ymin = emmean_bt - SE_bt, ymax = emmean_bt + SE_bt),
                linewidth = 0.2, width = 0.2) +
  geom_point(shape = 15, size = 4, aes(col=status)) +
  geom_text(nudge_x = 0.2, color = "black") +
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw() +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  labs(title = expression(Status~effect~on~NEE~(CO[2]+CH[4])),
       subtitle = paste0("Danube delta"),
       y=expression(EMM~(g~CO[2~eq]~d^-1~m^-2)),
       x="Conservation status",
       col=paste0("Status"))

plot_co2andch4_status_DA
ggsave(plot=plot_co2andch4_status_DA, filename="CO2plusCH4_status_plot_DA.png",
       device = "png",path = plots_path,dpi = 400,
       width = 115, height = 100,units = "mm")

#Dutch delta: 
plot_co2andch4_status_DU<- co2andch4_casepilot_comparisons %>%
  filter(casepilot == "DU", comparison == "status") %>%
  ggplot(., aes(x = status, y = emmean_bt, label = group_letter)) +
  geom_errorbar(aes(ymin = emmean_bt - SE_bt, ymax = emmean_bt + SE_bt),
                linewidth = 0.2, width = 0.2) +
  geom_point(shape = 15, size = 4, aes(col=status)) +
  geom_text(nudge_x = 0.2, color = "black") +
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw() +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  labs(title = expression(Status~effect~on~NEE~(CO[2]+CH[4])),
       subtitle = paste0("Dutch delta"),
       y=expression(EMM~(g~CO[2~eq]~d^-1~m^-2)),
       x="Conservation status",
       col=paste0("Status"))

plot_co2andch4_status_DU
ggsave(plot=plot_co2andch4_status_DU, filename="CO2plusCH4_status_plot_DU.png",
       device = "png",path = plots_path,dpi = 400,
       width = 115, height = 100,units = "mm")


#Ria d'Aveiro: 
plot_co2andch4_status_RI<- co2andch4_casepilot_comparisons %>%
  filter(casepilot == "RI", comparison == "status") %>%
  ggplot(., aes(x = status, y = emmean_bt, label = group_letter)) +
  geom_errorbar(aes(ymin = emmean_bt - SE_bt, ymax = emmean_bt + SE_bt),
                linewidth = 0.2, width = 0.2) +
  geom_point(shape = 15, size = 4, aes(col=status)) +
  geom_text(nudge_x = 0.2, color = "black") +
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw() +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  labs(title = expression(Status~effect~on~NEE~(CO[2]+CH[4])),
       subtitle = paste0("Ria d'Aveiro"), 
       y=expression(EMM~(g~CO[2~eq]~d^-1~m^-2)),
       x="Conservation status",
       col=paste0("Status"))

plot_co2andch4_status_RI
ggsave(plot=plot_co2andch4_status_RI, filename="CO2plusCH4_status_plot_RI.png",
       device = "png",path = plots_path,dpi = 400,
       width = 115, height = 100,units = "mm")

#Valencian wetland: 
plot_co2andch4_status_VA<- co2andch4_casepilot_comparisons %>%
  filter(casepilot == "VA", comparison == "status") %>%
  ggplot(., aes(x = status, y = emmean_bt, label = group_letter)) +
  geom_errorbar(aes(ymin = emmean_bt - SE_bt, ymax = emmean_bt + SE_bt),
                linewidth = 0.2, width = 0.2) +
  geom_point(shape = 15, size = 4, aes(col=status)) +
  geom_text(nudge_x = 0.2, color = "black") +
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw() +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  labs(title = expression(Status~effect~on~NEE~(CO[2]+CH[4])),
       subtitle = paste0("Valencian wetland"), 
       y=expression(EMM~(g~CO[2~eq]~d^-1~m^-2)),
       x="Conservation status",
       col=paste0("Status"))

plot_co2andch4_status_VA
ggsave(plot=plot_co2andch4_status_VA, filename="CO2plusCH4_status_plot_VA.png",
       device = "png",path = plots_path,dpi = 400,
       width = 115, height = 100,units = "mm")


##Extra plots (to-do)-----
#Plots Seasonal effect 
#To do, little interest for overview ppt (seasonal effect is averaged across the three status)

#Example:
#Camargue: 
plot_co2_season_CA<- co2_casepilot_comparisons %>%
  filter(casepilot == "CA", comparison == "season") %>%
  ggplot(., aes(x = season, y = emmean_bt, label = group_letter)) +
  geom_errorbar(aes(ymin = emmean_bt - SE_bt, ymax = emmean_bt + SE_bt),
                linewidth = 0.2, width = 0.2) +
  geom_point(shape = 15, size = 4, aes(col=season)) +
  geom_text(nudge_x = 0.2, color = "black") +
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw() +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  labs(title = expression(Seasonal~effect~on~NEE~(CO[2]+CH[4])),
       subtitle = paste0("Camargue"), 
       caption = paste0("Seasonal effect is averaged across the three conservation status.\nBoxes indicate the EM mean, error bars indicate SE,\n",              "Letters indicate different groups"),
       y=expression(EMM~(g~CO[2~eq]~d^-1~m^-2)),
       x="Sampling campaign",
       col=paste0("Sampling"))

plot_co2_season_CA


##Plots Interaction effect
#To do, unclear interest for overview ppt (maybe interesting for some casepilots)

#Example:
#Camargue: interaction plot
plot_co2_interaction_CA<- co2_casepilot_comparisons %>%
  filter(casepilot == "CA", comparison == "interaction") %>%
  #ggplot inside brackets to be able to access casepilot name with ggtitle or define dodge for all
  { 
    #Define dodge width for separating status.
    dodge <- position_dodge(width = 0.6)
    
    ggplot(., aes(x = season, y = emmean_gwp, label = group_letter, group=status)) +
      geom_errorbar(aes(ymin = emmean_gwp - SE_gwp, ymax = emmean_gwp + SE_gwp),
                    linewidth = 0.2, width = 0.2, position = dodge) +
      geom_point(shape = 15, size = 4, aes(col=status), position = dodge) +
      geom_text(vjust = -2, color = "black", position = dodge) +
      geom_hline(yintercept=0, linetype="dashed")+
      theme_bw() +
      theme(axis.title = element_text(face = "bold"),
            axis.text = element_text(face = "bold"),
            plot.caption = element_text(hjust = 0)) +
      labs(title = expression(Seasonal~effect~on~NEE~(CO[2]+CH[4])),
           subtitle = paste0("Camargue"), 
           caption = paste0("Boxes indicate the EM mean, error bars indicate SE,\n",              "Letters indicate different groups"),
           y=expression(EMM~(g~CO[2~eq]~d^-1~m^-2)),
           x="Sampling campaign",
           col=paste0("Sampling"))
  }

plot_co2_interaction_CA







