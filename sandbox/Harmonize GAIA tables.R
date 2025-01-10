#Harmonize waypoint names of different users


gaia_path<- "C:/Users/Miguel/Dropbox/RESTORE4Cs - Fieldwork/Data/GHG/Gaia_GPS/"

#Jorge

rawjorge<- read.csv(paste0(gaia_path,"GAIA_table_Jorge.csv"))

harmonized_jorge<- rawjorge %>% 
  mutate(subsite=toupper(gsub(".geojson", "", gaia_file)),
         plotnum=gsub(" ", "", tolower(str_extract(string = title, "Plot  [0-9]{1,2}|Plot [0-9]{1,2}|plot [0-9]{1,2}"))),
         subsite=case_when(subsite=="S1-CA-A1"~"S1-CA-R2",
                           TRUE~subsite),#Subsite was missidentified
         plotcode=case_when(!is.na(plotnum)~paste(subsite, plotnum, sep="-"),
                            TRUE~NA),
         user="Jorge")

#Miguel

rawmiguel<- read.csv(paste0(gaia_path, "GAIA_table_Miguel.csv"))

harmonized_miguel<- rawmiguel %>% 
  # filter(!is.na(fullsize_url)) %>% 
  mutate(plotnum=tolower(str_extract(string=title, pattern="plot [0-9]{1,2}|plot[0-9]{1,2}")),
         subsite=toupper(str_extract(string=title, pattern="S[0-9]{1}-[A-Z]{2}-[A-Z]{1}[0-9]{1}")),
         plotcode=case_when(!is.na(plotnum)~paste(subsite, plotnum, sep="-"),
                            TRUE~NA),
         user="Miguel")

#Benj

rawbenj<- read.csv(paste0(gaia_path, "GAIA_table_Benj.csv"))

harmonized_benj<- rawbenj %>%
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
         user="Benj")


#Dani Morant: STILL TO WRITE HARMONIZATION 

#note from Dani: Descarta el S1 CA P1, porque las coordenadas estaban mal, y tendrás de esa campaña lo de Jorge (si no lo tubieras, puedes coger los datos que te envío, pero sin coordenadas correctas). En esa campaña, no tenía los guiones entre los códigos, si me dices, los cambio, aquí o en el Teams cuando lo subas, en un segundo. Si ya están inlcuidos los de Jorge, nada.

rawdanimorant<- read.csv(paste0(gaia_path, "GAIA_table_DaniMorant.csv"))

harmonized_danimorant<- rawdanimorant



harmonized_all<- rbind(harmonized_jorge,harmonized_miguel,harmonized_benj)
                       
harmonized_all %>% group_by(subsite) %>% 
  summarise(n=sum(n())) %>%  
  print(n=105)


harmonized_all %>% filter(is.na(subsite))


