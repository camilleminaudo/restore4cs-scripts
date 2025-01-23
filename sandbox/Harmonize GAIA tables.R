#Harmonize waypoint names of different users
#Clear WP
rm(list=ls())

library(tidyverse)


gaia_path<- "C:/Users/Miguel/Dropbox/RESTORE4Cs - Fieldwork/Data/GHG/Gaia_GPS/"

#Jorge

rawjorge<- read.csv(paste0(gaia_path,"GAIA_table_Jorge.csv"))

jorge_harmonized<- rawjorge %>% 
  mutate(subsite=toupper(gsub(".geojson", "", gaia_file)),
         plotnum=gsub(" ", "", tolower(str_extract(string = title, "Plot  [0-9]{1,2}|Plot [0-9]{1,2}|plot [0-9]{1,2}"))),
         subsite=case_when(subsite=="S1-CA-A1"~"S1-CA-R2",
                           TRUE~subsite),#Subsite was missidentified
         plotcode=case_when(!is.na(plotnum)~paste(subsite, plotnum, sep="-"),
                            TRUE~NA),
         user="Jorge") %>% 
  select(waypoint_id, title, time_created, latitude, longitude, notes, fullsize_url, gaia_file, gaia_user, subsite, plotcode)

#Miguel

rawmiguel<- read.csv(paste0(gaia_path, "GAIA_table_Miguel.csv"))

miguel_harmonized<- rawmiguel %>% 
  # filter(!is.na(fullsize_url)) %>% 
  mutate(plotnum=tolower(str_extract(string=title, pattern="plot [0-9]{1,2}|plot[0-9]{1,2}")),
         subsite=toupper(str_extract(string=title, pattern="S[0-9]{1}-[A-Z]{2}-[A-Z]{1}[0-9]{1}")),
         plotcode=case_when(!is.na(plotnum)~paste(subsite, plotnum, sep="-"),
                            TRUE~NA),
         user="Miguel") %>% 
  select(waypoint_id, title, time_created, latitude, longitude, notes, fullsize_url, gaia_file, gaia_user, subsite, plotcode)

#Benj

rawbenj<- read.csv(paste0(gaia_path, "GAIA_table_Benj.csv"))

benj_harmonized<- rawbenj %>%
  mutate(title=gsub("–","-",title)) %>% 
  mutate(title=gsub("S44","S4",title)) %>% 
  mutate(title=gsub("S1-R1-P1-","S1-RI-P1-", title)) %>% 
  mutate(title=gsub("S4-CU—R1-","S4-CU-R1-", title)) %>% 
  mutate(title=gsub("S4-VaR2-","S4-VA-R2-",title)) %>%  
  mutate(title=gsub("S4-Va-P1-","S4-VA-P1-",title)) %>% 
  mutate(title=gsub("SE-VA-A2-","S3-VA-A2-",title)) %>% 
  mutate(title=case_when(waypoint_id=="8f22ec6d9d9063b4d1c2f76f2816adb8"~"S2-CU-R1-GHG11",
                         TRUE~title)) %>% 
  mutate(title=case_when(waypoint_id=="16457d617158dff9e286773613003430"~"S3-RI-P1-GHG10",
                         TRUE~title)) %>% # waypoint_id 16457d617158dff9e286773613003430 has picture and is floating chamber, title: My Photo - 09.04.24, 12:35:42, from time it should be: S3-RI-P1-GHG10 

  mutate(plotnum=gsub("GHG","", str_extract(string=title, pattern="GHG[0-9]{1,2}")),
         subsite=toupper(str_extract(string=title, pattern="S[0-9]{1}-[A-Z]{2}-[A-Z]{1}[0-9]{1}")),
         plotcode=case_when(!is.na(plotnum)~paste0(subsite, "-plot",plotnum),
                            TRUE~NA),
         user="Benj") %>% 
  select(waypoint_id, title, time_created, latitude, longitude, notes, fullsize_url, gaia_file, gaia_user, subsite, plotcode)


#Dani Morant: 
#note from Dani: Descarta el S1 CA P1, porque las coordenadas estaban mal, y tendrás de esa campaña lo de Jorge (si no lo tubieras, puedes coger los datos que te envío, pero sin coordenadas correctas). En esa campaña, no tenía los guiones entre los códigos, si me dices, los cambio, aquí o en el Teams cuando lo subas, en un segundo. Si ya están inlcuidos los de Jorge, nada.

rawdanimorant<- read.csv(paste0(gaia_path, "GAIA_table_DaniMorant.csv"))

danimorant_harmonized<- rawdanimorant%>%
  filter(!grepl("rising tide", title)) %>% #remove water samples
  mutate(sampling=str_extract(title, pattern="S[0-9]{1}"),
         case_pilot=str_extract(title, pattern="CA|CU|DA|DU|RI|VA"),
         subsitestatus=str_extract(title, pattern="A1|A2|P1|P2|R1|R2"),
         subsite=paste(sampling,case_pilot, subsitestatus,sep="-"),
         plot_number=parse_number(str_extract(title, pattern="plot[0-9]{1,2}")),
         space_plotnumber_space=parse_number(str_extract(title, pattern=" [0-9]{1,2} | [0-9]{1,2}")),
         plotnum=case_when(!is.na(plot_number)~plot_number, 
                           !is.na(space_plotnumber_space)~space_plotnumber_space,
                           TRUE~NA),
         plotcode=paste0(subsite, "-plot",plotnum)) %>% 
  #Manual modifications based on inspection of data:
  mutate(plotcode=case_when(waypoint_id=="c863ac1b-355e-4f10-a152-933ffc59598b"~"S1-CA-P1-plot3",
                            waypoint_id=="54ea1e19-c56c-4449-a98a-1244d4fbb849"~"S1-CA-P1-plot4",
                            waypoint_id=="4ac4019b-da1c-41bb-9a3a-676164e7867f"~"S1-CA-P2-plot3",
                            waypoint_id=="186602ed-4349-4bf1-9801-2e14aa6e72dd"~"S1-CA-R2-plot4",
                            waypoint_id=="b21abc29-dde4-48dc-9dab-f35acfd01b18"~"S1-VA-P1-plot2",
                            waypoint_id=="f600757c-c372-482a-a7f5-77ba721a9500"~"S1-VA-P1-plot3",
                            waypoint_id=="390c7a67-dbbc-45f6-92cd-746d5783e28b"~"S1-VA-P1-plot4",
                            waypoint_id=="8844f7d1-3dfb-438a-81d1-679922426121"~"S1-VA-P1-plot5",
                            waypoint_id=="1a615086-c767-413c-b3c1-7a256925d937"~"S1-CA-P2-plot8",
                            TRUE~plotcode),
         subsite=substr(plotcode,1,8)) %>%
  #Discard lat/long from S1-CA-P1 as per Dani's comments
  mutate(latitude=case_when(subsite=="S1-CA-P1"~NA,
                            TRUE~latitude),
         longitude=case_when(subsite=="S1-CA-P1"~NA,
                             TRUE~longitude)) %>% 
  select(waypoint_id, title, time_created, latitude, longitude, notes, fullsize_url, gaia_file, gaia_user, subsite, plotcode)


#Camille
rawcamille<- read.csv(paste0(gaia_path, "GAIA_table_Camille.csv"))

camille_harmonized<- rawcamille%>%
  filter(!grepl("sed|core|zoo|possible|crazy|invertebrates|R6", title)) %>% #remove not GHG plots 
  filter(!waypoint_id%in%c("86bbeab3-9fba-4df1-8ad0-de79f44c1b08","1fa0c680-60d4-4084-b228-c73d320c5cbb")) %>% #Remove not GHG plots
  mutate(subsite=str_extract(string=title, pattern="S[0-9]{1}-[A-Z]{2}-[A-Z]{1}[0-9]{1}"),
         folder=toupper(gsub(".geojson","",gaia_file)),
         subsite=case_when(is.na(subsite)~folder,
                           TRUE~subsite),
         plot_space_number=parse_number(str_extract(title, "plot [0-9]{1,2}")),
         plot_number=parse_number(str_extract(title, "plot[0-9]{1,2}")),
         plotnum=case_when(!is.na(plot_space_number)~plot_space_number, 
                           !is.na(plot_number)~plot_number,
                           TRUE~NA),
         plotcode=paste0(subsite,"-plot",plotnum)) %>% 
  select(waypoint_id, title, time_created, latitude, longitude, notes, fullsize_url, gaia_file, gaia_user, subsite, plotcode)

#OBS: duplicated plotcodes with different waypoint_id and coordinates for subsite S1-DA-P1 (not sure why but comming directly from gaia files)




#JOIN ALL DATASETS
harmonized_all<- rbind(benj_harmonized,miguel_harmonized,jorge_harmonized, danimorant_harmonized, camille_harmonized)

#Save
write.csv(harmonized_all, file = paste0(gaia_path,"full_gaia_table_allusers.csv"),row.names = F)


pictures_only<- harmonized_all %>% 
  filter(!is.na(fullsize_url))

write.csv(pictures_only, file = paste0(gaia_path,"Pictures_urls.csv"),row.names = F)
 

