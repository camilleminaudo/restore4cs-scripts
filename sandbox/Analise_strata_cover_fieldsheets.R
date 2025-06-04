
#Chamber strata distribution.

# ---
# Authors: Miguel Cabrera
# Project: "RESTORE4Cs"
# date: "June 2025"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---


# --- Description of this script
# This script is used to analyse the camber distribution across subsites and seasons. Results will help in identifying dominant strata per subsite, seasonal variability in dominance and relevant strata for comparisons. 

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
library(tools)


repo_root <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
files.sources = list.files(path = paste0(repo_root,"/functions"), full.names = T)
for (f in files.sources){source(f)}


# ---- Directories ----
dropbox_root <- "C:/Users/Miguel/Dropbox/RESTORE4Cs - Fieldwork/Data" # You have to make sure this is pointing to the write folder on your local machine
results_path <- paste0(dropbox_root,"/GHG/Processed data/computed_flux/")
export_path<- paste0(dropbox_root,"/GHG/Processed data/computed_flux/Summaries for MA/")
fieldsheet_path <- paste0(dropbox_root,"/GHG/Fieldsheets")

vegetation_path<- paste0(dropbox_root, "/Vegetation/")

#----Import fieldsheets----

#Read all GHGfieldsheets (using function from Camille)
fieldsheets<- list.files(path= fieldsheet_path, recursive = T, pattern = "GHG.xlsx",full.names = T)
field<- read_GHG_fieldsheets(fieldsheets)


#----Import vegetation data----

veg_data<- read.csv(paste0(vegetation_path,"RESTORE4Cs_finalvegetation_ABG_biomass_description.csv"))


#----Join&format----
chambers<- field %>% 
  mutate(plotcode=paste(subsite, plot_id,sep = "-")) %>% 
  select(pilot_site,subsite,date,plot_id,plotcode,chamber_type,strata,water_depth,comments,transparent_dark) %>% 
  rename(sampling=subsite, og_strata=strata) %>% 
  mutate(comments=if_else(is.na(comments),"",comments)) %>% 
  pivot_wider(names_from = transparent_dark, values_from = comments) %>% 
  merge.data.frame(veg_data %>% select(plotcode,vegetation_description), by="plotcode", all = T) %>% 
  mutate(season=substr(sampling, 1,2),
         subsite=substr(sampling,4,8),
         strata=case_when(og_strata=="vegetated"&water_depth>0~"vegetated water",
                          og_strata=="vegetated"&water_depth<=0~"vegetated dry",
                          TRUE~og_strata)) 

#---Seasonal_coverage----
#Obtain plot-level summary (all subsites, RI includes high-tide)
seasonal_coverage<- chambers %>%
  group_by(pilot_site, season, subsite, sampling, strata) %>%
  summarise(count = n(), .groups = "drop") %>% 
  group_by(pilot_site, season, subsite, sampling) %>% 
  mutate(total_plots=sum(count),
         percent=count/total_plots*100) %>% 
  mutate(strata=factor(strata, levels=c("bare", "vegetated dry", "vegetated water", "open water")))


ri_lowtide_seasonal_coverage<- chambers %>% 
  filter(pilot_site=="RI") %>% 
  #Remove high-tide chambers (based on harmonized comments:High-tide, Rising-tide, Receding-tide) leaving low-tide: Tidal-pool
  filter(!grepl("High-tide|Rising-tide|Receding-tide",dark)) %>% 
  group_by(pilot_site, season, subsite, sampling, strata) %>%
  summarise(count = n(), .groups = "drop") %>% 
  group_by(pilot_site, season, subsite, sampling) %>% 
  mutate(total_plots=sum(count),
         percent=count/total_plots*100) %>% 
  mutate(strata=factor(strata, levels=c("bare", "vegetated dry", "vegetated water", "open water")))
  



strata_colors <- c(
  "bare" = "#d95f02",              # brown/orange
  "vegetated dry" = "#228B22",     # green
  "vegetated water" = "#66c2a5",   # green-blue mix
  "open water" = "#1f78b4"         # blue
)


#Pilot-site plots of seasonal strata distribution in each subsite

#CA chambers-----
seasonal_coverage %>% filter(pilot_site == "CA") %>% 
ggplot(aes(x = season, y = percent, fill = strata)) +
  geom_bar(stat = "identity") +
  facet_wrap(~subsite) +
  geom_text(aes(label = count), 
            position = position_stack(vjust = 0.5), 
            color = "white", size = 3) +
  geom_text(aes(x = season, y = 100, label = total_plots),  # y = 100 puts it just above stacked bar
            inherit.aes = FALSE,
            vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = strata_colors, drop = FALSE) +  # keep empty strata in legend

  scale_y_continuous(labels = scales::percent_format(scale = 1), expand = expansion(mult = c(0, 0.1))) +
  labs(y = "Strata Chambers (%)", x = "Season", fill = "Strata") +
  theme_minimal()
  

#CU chambers-----
seasonal_coverage %>% filter(pilot_site == "CU") %>% 
  ggplot(aes(x = season, y = percent, fill = strata)) +
  geom_bar(stat = "identity") +
  facet_wrap(~subsite) +
  geom_text(aes(label = count), 
            position = position_stack(vjust = 0.5), 
            color = "white", size = 3) +
  geom_text(aes(x = season, y = 100, label = total_plots),  # y = 100 puts it just above stacked bar
            inherit.aes = FALSE,
            vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = strata_colors, drop = FALSE) +  # keep empty strata in legend
  
  scale_y_continuous(labels = scales::percent_format(scale = 1), expand = expansion(mult = c(0, 0.1))) +
  labs(y = "Strata Chambers (%)", x = "Season", fill = "Strata") +
  theme_minimal()

#Only 5 vegetated dry, phragmites, likely similar GHG profile to vegetation water 
chambers %>% filter(pilot_site=="CU") %>% filter(strata=="vegetated dry") %>%
  select(sampling, water_depth, vegetation_description,dark, transparent)

#Vegetation dominated by reeds, phragmites mostly
chambers %>% filter(pilot_site=="CU") %>% filter(og_strata=="vegetated") %>% 
  group_by(vegetation_description) %>% count()


#DA chambers-----
seasonal_coverage %>% filter(pilot_site == "DA") %>% 
  ggplot(aes(x = season, y = percent, fill = strata)) +
  geom_bar(stat = "identity") +
  facet_wrap(~subsite) +
  geom_text(aes(label = count), 
            position = position_stack(vjust = 0.5), 
            color = "white", size = 3) +
  geom_text(aes(x = season, y = 100, label = total_plots),  # y = 100 puts it just above stacked bar
            inherit.aes = FALSE,
            vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = strata_colors, drop = FALSE) +  # keep empty strata in legend
  
  scale_y_continuous(labels = scales::percent_format(scale = 1), expand = expansion(mult = c(0, 0.1))) +
  labs(y = "Strata Chambers (%)", x = "Season", fill = "Strata") +
  theme_minimal()

#DA-A1 bare includes incubations after plant removal
chambers %>% filter(subsite=="DA-A1") %>% filter(strata=="bare") %>% 
  arrange(sampling, plot_id) %>% 
  select(plotcode, transparent,dark)



#DU chambers-----
seasonal_coverage %>% filter(pilot_site == "DU") %>% 
  ggplot(aes(x = season, y = percent, fill = strata)) +
  geom_bar(stat = "identity") +
  facet_wrap(~subsite) +
  geom_text(aes(label = count), 
            position = position_stack(vjust = 0.5), 
            color = "white", size = 3) +
  geom_text(aes(x = season, y = 100, label = total_plots),  # y = 100 puts it just above stacked bar
            inherit.aes = FALSE,
            vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = strata_colors, drop = FALSE) +  # keep empty strata in legend
  
  scale_y_continuous(labels = scales::percent_format(scale = 1), expand = expansion(mult = c(0, 0.1))) +
  labs(y = "Strata Chambers (%)", x = "Season", fill = "Strata") +
  theme_minimal()

#Very few vegetated water: re-check the need for distinction with updated GHG fluxes.
chambers %>% filter(pilot_site=="DU") %>% filter(strata=="vegetated water") %>% 
 select(plotcode, water_depth, transparent,dark,vegetation_description)


#RI chambers (high & lowtide)-----
#INCLUDING High tide
seasonal_coverage %>% filter(pilot_site == "RI") %>% 
  ggplot(aes(x = season, y = percent, fill = strata)) +
  geom_bar(stat = "identity") +
  facet_wrap(~subsite) +
  geom_text(aes(label = count), 
            position = position_stack(vjust = 0.5), 
            color = "white", size = 3) +
  geom_text(aes(x = season, y = 98, label = total_plots),  # y = 100 puts it just above stacked bar
            inherit.aes = FALSE,
            vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = strata_colors, drop = FALSE) +  # keep empty strata in legend
  scale_y_continuous(labels = scales::percent_format(scale = 1), expand = expansion(mult = c(0, 0.1))) +
  labs(y = "Strata Chambers (%)", x = "Season", fill = "Strata") +
  theme_minimal()+
  ggtitle("High-tide included")

#Lowtide ONLY
ri_lowtide_seasonal_coverage %>% filter(pilot_site == "RI") %>% 
  ggplot(aes(x = season, y = percent, fill = strata)) +
  geom_bar(stat = "identity") +
  facet_wrap(~subsite) +
  geom_text(aes(label = count), 
            position = position_stack(vjust = 0.5), 
            color = "white", size = 3) +
  geom_text(aes(x = season, y = 98, label = total_plots),  # y = 100 puts it just above stacked bar
            inherit.aes = FALSE,
            vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = strata_colors, drop = FALSE) +  # keep empty strata in legend
  scale_y_continuous(labels = scales::percent_format(scale = 1), expand = expansion(mult = c(0, 0.1))) +
  labs(y = "Strata Chambers (%)", x = "Season", fill = "Strata") +
  theme_minimal()+
  ggtitle("High-tide excluded")

#Check vegetated water: same vegetation, no differenciation needed. 
chambers %>% filter(pilot_site=="RI") %>% 
  filter(strata=="vegetated water") %>% 
  select(plotcode, water_depth, transparent, dark)

#Check bare in Preserved and restored
chambers %>% filter(pilot_site=="RI") %>% 
  filter(strata=="bare") %>%
  filter(!subsite%in%c("RI-A1","RI-A2")) %>% 
  arrange(subsite, sampling, plotcode) %>% 
  select(sampling, dark)

#RI-A1: check incubations with water_depth
chambers %>% filter(pilot_site=="RI") %>% filter(subsite=="RI-A1") %>% 
  filter(!grepl("High-tide|Rising-tide|Receding-tide",dark)) %>% 
  filter(water_depth!=0) %>% select(plotcode, water_depth, transparent,dark)



#VA chambers-----

seasonal_coverage %>% filter(pilot_site == "VA") %>% 
  ggplot(aes(x = season, y = percent, fill = strata)) +
  geom_bar(stat = "identity") +
  facet_wrap(~subsite) +
  geom_text(aes(label = count), 
            position = position_stack(vjust = 0.5), 
            color = "white", size = 3) +
  geom_text(aes(x = season, y = 100, label = total_plots),  # y = 100 puts it just above stacked bar
            inherit.aes = FALSE,
            vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = strata_colors, drop = FALSE) +  # keep empty strata in legend
  scale_y_continuous(labels = scales::percent_format(scale = 1), expand = expansion(mult = c(0, 0.1))) +
  labs(y = "Strata Chambers (%)", x = "Season", fill = "Strata") +
  theme_minimal()


chambers %>% filter(subsite=="VA-P1"&og_strata=="vegetated") %>% select(sampling, vegetation_description)
