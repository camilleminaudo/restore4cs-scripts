
# ---
# Authors: Miguel Cabrera
# Project: "RESTORE4Cs"
# date: "March 2025"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
# This script is used to gather the vegetated chambers information, matching them with the vegetation data (dry weight), and the GAIA dataset. A manual step is done to fill-in missing vegetation description and to harmonize the way it is presented. A final csv is created with plotcode (season-site-subsite-plotnum), Abovegroundbiomass per square meter, vegetation description and origin of vegetation description (for internal use).


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
fieldsheetpath <- paste0(dropbox_root,"/GHG/Fieldsheets")
vegetation_path<- paste0(dropbox_root,"/Vegetation/")



#1. Import GHG fieldsheets----

#Read all GHGfieldsheets 
fieldsheets<- list.files(path= fieldsheetpath, recursive = T, pattern = "GHG.xlsx",full.names = T)
field<- read_GHG_fieldsheets(fieldsheets)

#Drop unwanted variables and rename key UniqueID:
field<- field %>% 
  select(-c(person_sampling,gas_analyzer,logger_floating_chamber,logger_transparent_chamber,start_time,initial_co2,initial_ch4,end_time,final_ch4,final_co2,chamber_height_cm, unix_start,unix_stop)) %>% 
  rename(UniqueID=uniqID)



# 2. Vegetation ID & DW ----

#List all vegetation files 
#Import pattern: RESTORE4Cs_vegetation.........xlsx (rbind all). 
veg_files<- list.files(path = vegetation_path, pattern = "^RESTORE4Cs_vegetation_",full.names = T,recursive = T)

#1. Read all vegetation files and combine them
veg<- read_vegetation_fieldsheets(veg_files)

#Check Duplicated plant plotcode
veg %>% 
  dplyr::summarise(n = dplyr::n(), .by = c(plotcode)) |>dplyr::filter(n > 1L) %>% 
  merge.data.frame(veg, by=c("plotcode"), all = F)


#2. Extract from fieldsheets the comments of vegetated plots. Format into comments_dark and comments_transparent
vegetated_plots<- field %>% 
  filter(strata=="vegetated") %>%#Select only vegetated plots
  mutate(plotcode=paste(subsite, plot_id,sep="-")) %>% 
  select(plotcode, strata,transparent_dark, comments) %>% 
  pivot_wider(values_from = comments, names_from = transparent_dark) %>% 
  distinct()


#2.1 Add GAIA notes and picture availability to vegetated_plots

gaiadata_veg<- read.csv(paste0(dropbox_root,"/GHG/Gaia_GPS/full_gaia_table_allusers.csv")) %>% 
  select(plotcode, notes, fullsize_url) %>% 
  mutate(plotcode=gsub("plot","",plotcode),
         picture_available=!is.na(fullsize_url)) %>% 
  rename(comment_gaia=notes)%>%
  select(plotcode, comment_gaia, picture_available) %>% 
  group_by(plotcode) %>% 
  summarise(pictures_available=sum(picture_available, na.rm = T),
            comments_gaia = paste(unique(comment_gaia[comment_gaia != "" & !is.na(comment_gaia)]), collapse = ", "),#Combined notes of all observations with same plotcode (those with multiple pictures) 
            .groups="drop") %>% 
  filter(plotcode%in%vegetated_plots$plotcode)



#Join fieldsheet and vegetation info
vegetated_plots_veg<- vegetated_plots %>% 
  merge.data.frame(veg, by=c("plotcode"),all = T) %>% 
  merge.data.frame(gaiadata_veg, by=c("plotcode"),all=T)

#Check vegetation DW data without vegetated plotGHG: 7 unmatched vegetation samples
veg[!veg$plotcode%in%vegetated_plots$plotcode,]


#S2-RI, S1-CA, S2-CA, S2-CU, S3-VA, S4-CA: unmatched vegetation data was corrected based on sample distribution of chambers and vegetation samples. Modifications are commented in vegetation excels.  
#S2-RI had some correction, but remaining unmatched vegetation data CANNOT be assigned to vegetated chambers (no obvious correspondence or typos)


#Vegetated plotsGHG without vegetationDW 90 chambers without vegetation data
print(vegetated_plots[!vegetated_plots$plotcode%in%veg$plotcode,],n=150)



#3. Join and combine comments to get vegetation description (sp/genus/family/), calculate ABG_biomass using tube area (pi*(12.1cm radius)^2)
#Missidentified vegetation data (veg data without chamber) are excluded
vegetation_dw_final<- vegetated_plots %>% 
  merge.data.frame(x=.,y = veg, by=c("plotcode"),all.x = T) %>% 
  merge.data.frame(x=.,y= gaiadata_veg, by=c("plotcode"),all.x = T) %>% 
  mutate(ABG_biomass_gpersquaremeter=dry_weight/(pi*0.121^2)) %>%# g/square m (radius in meters) 
  select(plotcode,
         strata, ABG_biomass_gpersquaremeter, transparent, dark, comments, comments_gaia, pictures_available) %>% 
  rename(comment1=transparent, comment2=dark, comment3=comments, gaiaNotes=comments_gaia, gaiaPictures=pictures_available) %>% 
  mutate(comment1=case_when(is.na(comment1)~"",TRUE~comment1),
         comment2=case_when(is.na(comment2)~"",TRUE~comment2),
         comment3=case_when(is.na(comment3)~"",TRUE~comment3),
         gaiaNotes=case_when(is.na(gaiaNotes)~"",TRUE~gaiaNotes))


#Save in vegetation folder
write.csv(vegetation_dw_final,file = paste0(vegetation_path, "RESTORE4Cs_ABGbiomass_vegetated_plots.csv"), row.names = F)


rm(veg, vegetated_plots,vegetated_plots_veg,vegetation_dw_final, gaiadata_veg)


#3. Manual step (done)----
#Harmonize plot vegetation description from all fieldsheet comments (GHG and Veg_files), additionally incorporating gaiaNotes (notes saved in GAIA app with the location) and inspecting gaiaPictures to identify missing vegetation.

#Done, the harmonization is saved as an excel in the vegetation folder.

#4. Final Format-----
#Import the harmonized dataset
vegetation<- read_xlsx(path = paste0(vegetation_path,"RESTORE4Cs_ABGbiomass_vegetated_plots_picturedescription.xlsx"),na = "NA")

#Perform the final formatting (combine as text the 3 vegetation description columns)
vegetation_final_format<- vegetation %>% 
  select(plotcode, strata, ABG_biomass_gpersquaremeter, originVegID, 
         Vegetation_description_1,
         Vegetation_description_2,
         Vegetation_description_3) %>% 
  mutate(vegetation_description=paste0(Vegetation_description_1,"; ",Vegetation_description_2,"; ",Vegetation_description_3)) %>% 
  mutate(vegetation_description=gsub("; NA","", vegetation_description)) %>% 
  mutate(vegetation_description=if_else(vegetation_description=="NA",NA,vegetation_description),
         ABG_biomass_gpersquaremeter=as.numeric(ABG_biomass_gpersquaremeter)) %>% 
  select(plotcode, ABG_biomass_gpersquaremeter, originVegID, vegetation_description)



#Save final formatted vegetation
write.csv(vegetation_final_format, file = paste0(vegetation_path,"RESTORE4Cs_finalvegetation_ABG_biomass_description.csv"),row.names = F)



#Check final dataset:
vegetation_final_format %>% 
  select(-originVegID) %>% 
  summarise(total_samples=sum(!is.na(plotcode)),
            n_biomasses=sum(!is.na(ABG_biomass_gpersquaremeter)),
            n_vegetation_description=sum(!is.na(vegetation_description)))
