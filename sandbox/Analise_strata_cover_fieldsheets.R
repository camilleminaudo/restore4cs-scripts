
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
fieldsheetpath <- paste0(dropbox_root,"/GHG/Fieldsheets")



#----Import fieldsheets----

#Read all GHGfieldsheets (using function from Camille)
fieldsheets<- list.files(path= fieldsheetpath, recursive = T, pattern = "GHG.xlsx",full.names = T)
field<- read_GHG_fieldsheets(fieldsheets)

names(field)
#Obtain plot-level summary
strata2<- field %>% 
  select(pilot_site,subsite,date,plot_id,chamber_type,strata,water_depth,comments,transparent_dark) %>% 
  rename(sampling=subsite, og_strata=strata) %>% 
  pivot_wider(names_from = transparent_dark, values_from = comments) %>% 
  mutate(season=substr(sampling, 1,2),
         subsite=substr(sampling,4,8),
         strata=case_when(og_strata=="vegetated"&water_depth>1~"vegetated water",
                          og_strata=="vegetated"&water_depth<=1~"vegetated dry",
                          TRUE~og_strata)) 

seasonal_coverage<- strata2 %>%
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

#RI chambers-----
seasonal_coverage %>% filter(pilot_site == "RI") %>% 
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



