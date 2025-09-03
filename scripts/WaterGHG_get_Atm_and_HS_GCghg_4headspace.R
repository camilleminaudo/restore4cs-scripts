# Get Atm and HS CO2 & CH4 from GCdata



#Description----
# ---
# Authors: Miguel Cabrera Brufau
# Project: "RESTORE4Cs"
# date: "August 2025"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
#MOTIVATION: We need a daily calcurve for the GC data (CO2 and CH4). WE will only use the GC data for headspaces (and corresponding atmospheres). Cores and chambers already have more reliable fluxes from Licor injections and in-situ GHG timeseries. We will not use any N2O data from GC (not reliable).
#Based on appropriate tests, the waterGHG concentration is linearly dependent on the calibration factor of the source atm and headspace concentrations used to calculate it. 


#FUNCTIONING: This script uses the main GC excel to calculate daily calibration curves using the initial calcurve of each day. And 
#We will check any potential biass by comparing with internal standards (inter-batch) and external standards (atmospheres with concentration measured in-situ). Depending on the identified biasses, we will correct the daily calfactors using the internal and external standards.


#Requirements: properly formated GC excell, Fieldsheet exetainer correspondence for atmospheres, and in-situ atmosphere concenttrations. 

rm(list = ls())

# ---- Directories ----
dropbox_root <- "C:/Users/Miguel/Dropbox/RESTORE4Cs - Fieldwork/Data/"

fieldsheet_path <- paste0(dropbox_root,"GHG/Fieldsheets")
headspace_path<- paste0(dropbox_root,"GHG/Headspaces/")
gc_path<- paste0(dropbox_root,"GHG/GC_data/")


#---Packages and functions -----

library(tidyverse)
library(readxl)
library(broom)
library(ggpmisc)
library(egg)

repo_root <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
files.sources = list.files(path = paste0(repo_root,"/functions"), full.names = T)
for (f in files.sources){source(f)}


# ---- function to save a list of plots into pdf file
gg_save_pdf = function(list, filename) {
  pdf(filename,width = 14,height = 7)
  for (p in list) {
    print(p)
  }
  dev.off()
  invisible(NULL)
}




#Steps-----

#DONE: produce dailycalcurve only with calibration vials (batch==0).
#DONE: inspect&remove suspicious calibration vials. 
#DONE: check calibration bias with P5 vials and atmosphere exetainers
#DONE: decide on corrective action for calibration factors (using internal and external standards)

#TODO: save definitive calibration factors (area-to-ppm for CO2 and CH4)


#1. Import data ----

##1.1. GCdata-----
#Data import and formatting
gcog<- read_xlsx(path = paste0(gc_path, "DATA_R4Cs_with_datetimes_rcolnames_daybatch_v20241119.xlsx")) %>% 
  mutate(minuteson=difftime(datetime,startup, units = "mins"),
         dayofanalysis=as.Date(as.character(yearmonthday), format="%Y%m%d"),
         julianday=as.numeric(julian(dayofanalysis,origin=as.Date("2024-01-01")))) %>% 
  filter(dayofanalysis>as.Date("2024-01-01")) %>% #Data before this date is not reliable
  mutate(uniqueinj_id=paste(yearmonthday,daybatch,seq_position,sep="_")) %>% #add unique id with explicit info
  rename(vol=volstd) %>% 
  select(-c(n2o,base_n2o_rt6.7))#remove n2o data

#Identify days when samples were analyzed (the rest we do not care)
daysofsamples<- gcog %>% 
  group_by(dayofanalysis) %>% 
  summarise(num_vials=sum(grepl("^S", sample_id))) %>% 
  filter(num_vials>0) %>% 
  pull(dayofanalysis)

#Filter for days with samples
gcsamples<- gcog %>% 
  filter(dayofanalysis%in%daysofsamples)


##1.2. Exetainer fieldsheet-----
#Import fieldsheets exetainers field samples
fieldsheets_exe<-list.files(fieldsheet_path, pattern = "exetainers.xlsx", recursive = T, full.names = T)

for(i in fieldsheets_exe){
  
  a<- readxl::read_xlsx(i, trim_ws = T)
  a<- a[!is.na(a[,1]),]
  
  if(i==fieldsheets_exe[1]){A<- a}else {A<- rbind(A,a)}
}

all_exetainer_field<- A

rm(a,i,A)


##1.3. Atm ppm-----

subsite_atm<- read.csv(paste0(headspace_path, "in_situ_atm_co2&ch4.csv")) %>% 
  rename(subsite_ID=subsite) %>% 
  select(subsite_ID, ghgspecies, avg) %>% 
  pivot_wider(names_from = ghgspecies, values_from = avg) %>% 
  rename(co2=co2_ppm, ch4=ch4_ppb) %>% 
  mutate(ch4=ch4/1000)#get ch4 in ppm

head(subsite_atm)


#OPTION A: no blank correct-----
#One calibration curve per day (not forced to intercept), including blanks in calibration curves.



#2. Prepare caldata-----

##Get calvials------
calvials<- gcsamples %>% 
  filter(grepl("^P",sample_id)) %>% #get only standards
  filter(batch==0) %>%  #select only standards from cal-curves
  filter(!is.na(vol)) %>%  #remove Pmix standards
  # filter(vol!=11.2) %>% #remove P20 vials, not reliable
  pivot_longer(cols=c("co2","ch4"), names_to = "gas", values_to = "area") %>% #pivot longer by gas
  mutate(uniqueareaid=paste(uniqueinj_id, gas,sep="_")) %>% #create separate injid per gas
  mutate(ppmstd=case_when(gas=="co2"~3000*(vol/11.2),
                          gas=="ch4"~15.29*(vol/11.2)))


## Get stdvials (internal std)-----

#Inter-batch internal standards (P20, P5, some P2)
stdvials<- gcsamples %>% 
  filter(grepl("^P",sample_id)) %>% #get only standards
  filter(batch!=0) %>%  #Remove standards from calcurves
  filter(!is.na(vol)) %>%  #remove Pmix standards
  # filter(vol!=11.2) %>% #remove P20 vials, not reliable
  pivot_longer(cols=c("co2","ch4"), names_to = "gas", values_to = "area") %>% #pivot longer by gas
  mutate(uniqueareaid=paste(uniqueinj_id, gas,sep="_")) %>% #create separate injid per gas
  mutate(ppmstd=case_when(gas=="co2"~3000*(vol/11.2),
                          gas=="ch4"~15.29*(vol/11.2)))


##Get atmvials (external std)----

#Atmosphere exetainers with ppm imputed from in-situ data (between incubations)
atm_exetainers<- all_exetainer_field %>% 
  filter(`headspace/trapped/atm`=="atmosphere") %>% 
  select(subsite_ID, exetainer_ID) %>% 
  left_join(subsite_atm, by="subsite_ID") %>% 
  pivot_longer(cols = c(co2, ch4), names_to = "gas", values_to = "ppmstd") %>% 
  select(-subsite_ID) %>% 
  rename(sample_id=exetainer_ID)

head(atm_exetainers)

#Joined with GC areas based on sample_id
atmvials<- gcsamples %>% 
  filter(sample_id%in%atm_exetainers$sample_id) %>% 
  pivot_longer(cols=c("co2","ch4"), names_to = "gas", values_to = "area") %>% #pivot longer by gas
  mutate(uniqueareaid=paste(uniqueinj_id, gas,sep="_")) %>% #create separate injid per gas
  left_join(atm_exetainers, by=c("sample_id", "gas"))

head(atmvials)

rm(atm_exetainers)


#3. Inspect calibrations------

#Blank evolution (area co2 and ch4)

calvials %>% rbind(stdvials) %>% 
  filter(vol==0) %>% 
  mutate(caldate=as.Date(as.character(yearmonthday),format="%Y%m%d")) %>% 
  ggplot(aes(x=caldate, y=area,  col = uniqueareaid%in%calvials$uniqueareaid))+
  geom_boxplot(outliers = F)+
  geom_point()+
  facet_wrap(facets=.~gas, scales="free")


ch4_mean_blankarea<- calvials %>% rbind(stdvials) %>% 
  filter(vol==0) %>% 
  filter(gas=="ch4") %>% 
  filter(area<quantile(area, 0.85)) %>% 
  summarise(ch4_blank=mean(area, na.rm=T)) %>% pull(ch4_blank)



co2_mean_blankarea<- calvials %>% rbind(stdvials) %>% 
  filter(vol==0) %>% 
  filter(gas=="co2") %>% 
  filter(area<quantile(area, 0.6)) %>% 
  summarise(co2_blank=mean(area, na.rm=T)) %>% pull(co2_blank)

calvials %>% rbind(stdvials) %>% 
  filter(vol==0) %>% 
  filter(gas=="co2") %>% 
  filter(area<quantile(area, 0.6)) %>%
  ggplot(aes(x=1, y=area))+
  geom_hline(aes(yintercept = co2_mean_blankarea))+
  geom_violin()+
  geom_boxplot()




##Outlier areas-----
#Mannually write down unreliable areas for each gas from calvials and stdvials

co2_out<- c("20240109_1_6_co2","20240109_1_7_co2","20240110_1_7_co2","20240221_1_7_co2","20240222_1_7_co2","20240311_1_8_co2","20240402_1_8_co2","20240321_1_8_co2","20240402_1_2_co2","20240402_1_3_co2","20240404_1_3_co2","20240411_1_1_co2","20240411_1_2_co2","20240412_1_5_co2","20240415_1_3_co2","20240415_1_5_co2","20240424_1_8_co2","20240426_1_8_co2","20240502_1_8_co2","20240503_1_8_co2","20240508_1_8_co2","20240510_1_8_co2","20240521_1_1_co2","20240521_1_8_co2","20240522_1_1_co2","20240522_1_8_co2","20240523_1_1_co2","20240530_1_1_co2","20240610_1_8_co2","20240611_1_2_co2","20240703_1_4_co2","20240729_1_1_co2","20240729_1_3_co2","20240807_1_3_co2","20240809_1_3_co2","20240813_1_3_co2","20240814_1_2_co2","20240814_1_9_co2")
ch4_out<- c("20240109_1_6_ch4","20240109_1_7_ch4","20240110_1_7_ch4","20240221_1_7_ch4","20240222_1_7_ch4","20240311_1_8_ch4","20240402_1_8_ch4","20240402_1_2_ch4","20240402_1_3_ch4","20240321_1_8_ch4","20240404_1_3_ch4","20240405_1_7_ch4","20240410_1_8_ch4","20240411_1_1_ch4","20240411_1_2_ch4","20240415_1_3_ch4","20240415_1_5_ch4","20240424_1_8_ch4","20240426_1_8_ch4","20240502_1_8_ch4","20240503_1_8_ch4","20240506_1_8_ch4","20240507_1_8_ch4","20240508_1_8_ch4","20240507_3_12_ch4","20240510_1_8_ch4","20240513_1_3_ch4","20240521_1_1_ch4","20240521_1_8_ch4","20240531_4_12_ch4","20240522_1_1_ch4","20240522_1_8_ch4","20240523_1_1_ch4","20240530_1_1_ch4","20240531_1_8_ch4","20240610_1_8_ch4","20240611_1_2_ch4","20240612_1_8_ch4","20240613_1_8_ch4","20240703_1_4_ch4","20240729_1_1_ch4","20240729_1_3_ch4","20240807_1_3_ch4","20240809_1_3_ch4","20240813_1_3_ch4","20240814_1_2_ch4","20240814_1_9_ch4")

##Get cals-----

#Calculate all the cal curves parameters (one per each dayofanalysis)
cal<- calvials %>% 
  filter(!(uniqueareaid)%in%c(co2_out,ch4_out)) %>% 
  select(uniqueinj_id,uniqueareaid,yearmonthday,ppmstd,gas,area)%>%
  nest_by(yearmonthday,gas) %>% 
  mutate(mod = list(lm(area~ppmstd, data=data))) %>%
  reframe(tidy(mod),glance(mod))

##Plot cals-----


#Plot factors 
cal%>%
  filter(term=="ppmstd")%>%
  mutate(caldate=as.Date(as.character(yearmonthday),format="%Y%m%d")) %>% 
  # filter(r.squared>0.95) %>% 
  ggplot(aes(x=caldate, y=estimate, col=gas))+
  geom_point()+
  geom_point(data=data.frame(caldate=as.Date(as.character(c(20240502)),format="%Y%m%d"),
                             estimate=c(0,0),
                             gas=c("ch4","co2")))+
  facet_wrap(facets=.~gas, scales = "free")

#PLot intercepts

cal%>%
  filter(term=="(Intercept)")%>%
  mutate(caldate=as.Date(as.character(yearmonthday),format="%Y%m%d")) %>% 
  # filter(r.squared>0.95) %>% 
  ggplot(aes(x=caldate, y=estimate, col=gas))+
  geom_point()+
  geom_hline(aes(yintercept = co2_mean_blankarea))+
  facet_wrap(facets=.~gas, scales = "free")

#CH4cal evolution
cal %>% 
  filter(gas=="ch4") %>% 
  mutate(caldate=as.Date(as.character(yearmonthday),format="%Y%m%d")) %>% 
  # filter(r.squared>0.95) %>% 
  ggplot(aes(x=caldate, y=estimate, col=gas))+
  geom_point()+
  facet_wrap(facets=.~term, scales = "free")

#CO2cal evolution
cal %>% 
  filter(gas=="co2") %>% 
  mutate(caldate=as.Date(as.character(yearmonthday),format="%Y%m%d")) %>% 
  # filter(r.squared>0.95) %>% 
  ggplot(aes(x=caldate, y=estimate, col=gas))+
  geom_point()+
  facet_wrap(facets=.~term, scales = "free")


##Plot R2
cal%>%
  filter(term=="ppmstd")%>%
  mutate(caldate=as.Date(as.character(yearmonthday),format="%Y%m%d")) %>% 
  ggplot(aes(x=caldate, y=r.squared, col=gas))+
  geom_point()



#We need to inspect all the calibration curves one by one! 

#PLot every calcurve (after excluding outliers explictly, specified in "gas_out" vectors)

#Remove objects created in loop (for repeated executions of loop) 
rm(calcurve,cal_out,int_std,ext_std,dummydata,plt_CH4,plt_CO2,plt_N2O,datecalcurve)


max_co2<- calvials %>% 
  filter(!(uniqueareaid)%in%c(co2_out,ch4_out)) %>%
  filter(gas=="co2")%>% 
  summarise(max_area=max(area, na.rm = T)) %>% pull(max_area)

max_ch4<- calvials %>% 
  filter(!(uniqueareaid)%in%c(co2_out,ch4_out)) %>%
  filter(gas=="ch4")%>% 
  summarise(max_area=max(area, na.rm = T)) %>% pull(max_area)


#create vector to store calcurve plots
plt_list <- vector('list', length(unique(calvials$yearmonthday)))

#Loop over each date and plot the calcurve for all 2 ghgs
for(datecalcurve in unique(calvials$yearmonthday)){
  
  calcurve<- calvials %>% 
    filter(yearmonthday==datecalcurve) %>% 
    filter(!(uniqueareaid)%in%c(co2_out,ch4_out))
  
  cal_out<- calvials %>% 
    filter(yearmonthday==datecalcurve) %>% 
    filter(!(uniqueareaid%in%calcurve$uniqueareaid))
  
  int_std<- stdvials %>% 
    filter(yearmonthday==datecalcurve) %>% 
    filter(!uniqueareaid%in%c(co2_out, ch4_out))
  
  ext_std<- atmvials %>% 
    filter(yearmonthday==datecalcurve) %>% 
    filter(!uniqueareaid%in%c(co2_out, ch4_out))
  
  dummydata<- data.frame(ppmstd=c(0), area=as.numeric(NA_real_))
  
  plt_CO2 <- ggplot(subset(calcurve, gas=="co2"), aes(ppmstd, area))+
    geom_point()+
    # geom_text(aes(label=uniqueareaid))+
    stat_poly_line() +
    stat_poly_eq(use_label(c("eq", "adj.R2", "n")))+
    geom_point(data = subset(cal_out, gas=="co2"), aes(ppmstd,area, col="cal_out"))+
    geom_point(data = subset(int_std, gas=="co2"), aes(ppmstd,area, col="int_std"))+
    geom_point(data = subset(ext_std, gas=="co2"), aes(ppmstd,area, col="ext_std"))+
    # geom_text(data=subset(qual_inj, gas=="co2"), aes(vol, area,label=uniqueareaid))+
    geom_point(data=dummydata,aes(ppmstd,area, col = "cal_out"))+
    geom_point(data=dummydata,aes(ppmstd,area, col = "int_std"))+
    geom_point(data=dummydata,aes(ppmstd,area, col = "ext_std"))+
    scale_y_continuous(limits = c(0,max_co2))+
    theme_article()+
    labs(col = "") +
    theme(legend.position = "inside", legend.position.inside = c(0.8,0.1))+
    xlab("Std gas conc (ppm)")+
    ylab("CO2 Area (FID)")+
    ggtitle(paste0(datecalcurve," cal. plots"))
  
  plt_CH4 <- ggplot(subset(calcurve, gas=="ch4"), aes(ppmstd, area))+
    geom_point()+
    # geom_text(aes(label=uniqueareaid))+
    stat_poly_line() +
    stat_poly_eq(use_label(c("eq", "adj.R2", "n")))+
    geom_point(data = subset(cal_out, gas=="ch4"), aes(ppmstd,area, col="cal_out"))+
    geom_point(data = subset(int_std, gas=="ch4"), aes(ppmstd,area, col="int_std"))+
    geom_point(data = subset(ext_std, gas=="ch4"), aes(ppmstd,area, col="ext_std"))+
    # geom_text(data=subset(int_std, gas=="ch4"), aes(ppmstd, area,label=uniqueareaid))+
    geom_point(data=dummydata,aes(ppmstd,area, col = "cal_out"))+
    geom_point(data=dummydata,aes(ppmstd,area, col = "int_std"))+
    geom_point(data=dummydata,aes(ppmstd,area, col = "ext_std"))+
    scale_y_continuous(limits = c(0,max_ch4))+
    theme_article()+
    labs(col = "") +
    theme(legend.position = "inside", legend.position.inside = c(0.8,0.1))+
    xlab("Std gas conc (ppm)")+
    ylab("CH4 Area (FID)")
  
  
  plt_list[[which(unique(calvials$yearmonthday)==datecalcurve)]] <- ggarrange(plt_CO2, plt_CH4, ncol = 2,nrow = 1)
  
}

# Print pdf
# setwd(plots_path)
gg_save_pdf(list = plt_list, filename = paste0(gc_path,"calcurvesGC_co2&ch4_R4Cs.pdf"))


rm(calcurve,cal_out,int_std,ext_std,dummydata,plt_CH4,plt_CO2,datecalcurve)




#4. Calculate ppm -----
#for all samples

calforall<- cal %>% 
  select(yearmonthday,gas,term,estimate) %>% 
  pivot_wider(names_from = term, values_from = estimate) %>% 
  rename(intercept=`(Intercept)`,
         slope=ppmstd)


conc_all<- gcsamples %>% 
  pivot_longer(cols=c("co2","ch4"), names_to = "gas", values_to = "area") %>% #pivot longer by gas
  mutate(uniqueareaid=paste(uniqueinj_id, gas,sep="_")) %>% #create separate injid per gas
  left_join(calforall, by = c("gas", "yearmonthday")) %>% 
  filter(!(uniqueareaid)%in%c(co2_out,ch4_out)) %>% 
  mutate(ppm=(area-intercept)/slope) %>% 
  mutate(ppmstd=case_when(gas=="co2"~3000*(vol/11.2),
                          gas=="ch4"~15.29*(vol/11.2)))


#Check intercept vs blanks: 
conc_all %>% 
  filter(vol==0) %>% 
  mutate(dateofanalysis=as.Date(as.character(yearmonthday), format = "%Y%m%d")) %>% 
  ggplot(aes(x=dateofanalysis, y=ppm))+
  geom_point()+
  facet_wrap(facets=.~gas, scales="free")


conc_all %>% 
  filter(vol==0) %>% 
  mutate(dateofanalysis=as.Date(as.character(yearmonthday), format = "%Y%m%d")) %>% 
  ggplot(aes(x=intercept, y=ppm))+
  geom_point()+
  facet_wrap(facets=.~gas, scales="free")


#Check agreement dailycal vs all stdvials
conc_all %>% 
  filter(gas=="co2") %>% 
  ggplot(aes(x=ppmstd, y=ppm))+
           geom_point()+
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "n")))+
  scale_y_continuous(limits = c(-1000,4000))+
  scale_x_continuous(limits = c(-1000,4000))
  

conc_all %>% 
  filter(gas=="ch4") %>% 
  ggplot(aes(x=ppmstd, y=ppm))+
  geom_point()+
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "n")))+
  scale_y_continuous(limits = c(-2,20))+
  scale_x_continuous(limits = c(-2,20))

#Check atm 1:1 ppm (in situ vs GC-derived)

atmvials %>%
  left_join(calforall, by = c("gas", "yearmonthday")) %>% 
  filter(!(uniqueareaid)%in%c(co2_out,ch4_out)) %>% 
  mutate(ppm=(area-intercept)/slope) %>% 
  filter(gas=="co2") %>% 
  ggplot(aes(x=ppmstd, y=ppm))+
  geom_point()+
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "n")))


atmvials %>%
  left_join(calforall, by = c("gas", "yearmonthday")) %>% 
  filter(!(uniqueareaid)%in%c(co2_out,ch4_out)) %>% 
  mutate(ppm=(area-intercept)/slope) %>% 
  filter(gas=="ch4") %>% 
  ggplot(aes(x=ppmstd, y=ppm))+
  geom_point()+
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "n")))


atmvials %>%
  left_join(calforall, by = c("gas", "yearmonthday")) %>% 
  filter(!(uniqueareaid)%in%c(co2_out,ch4_out)) %>% 
  mutate(ppm=(area-intercept)/slope) %>% 
  select(-c(uniqueareaid,area,slope,intercept,ppmstd)) %>% 
  pivot_wider(names_from = gas, values_from = ppm,names_prefix = "ppm_") %>% 
  ggplot(aes(x=ppm_co2, y=ppm_ch4))+
  geom_point()+
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "n")))

#Check discrepancy ppmstd vs ppm
bias_cal_atm<- atmvials %>% 
  left_join(calforall, by = c("gas", "yearmonthday")) %>% 
  filter(!(uniqueareaid)%in%c(co2_out,ch4_out)) %>% 
  mutate(ppm=(area-intercept)/slope) %>% 
  mutate(bias_cal=ppm-ppmstd) %>% 
  select(sample_id,bias_cal, gas) %>% 
  pivot_wider(names_from = gas, values_from = bias_cal) %>% 
  #Calculate central 90% distribution of co2 and ch4 bias
  mutate(co2_bias_cal90central=case_when(co2<quantile(co2,0.05,na.rm = T)~NA_real_,
                                         co2>quantile(co2,0.95,na.rm = T)~NA_real_,
                                         TRUE~co2)) %>% 
  
  mutate(ch4_bias_cal90central=case_when(ch4<quantile(ch4,0.05,na.rm = T)~NA_real_,
                                         ch4>quantile(ch4,0.95,na.rm = T)~NA_real_,
                                         TRUE~ch4)) %>% 
  mutate(co2_mean90calbilas=mean(co2_bias_cal90central, na.rm=T),
         ch4_mean90calbilas=mean(ch4_bias_cal90central, na.rm=T),
         co2_rmse_calbias= sqrt(mean((co2_bias_cal90central)^2, na.rm=T)),
         ch4_rmse_calbias= sqrt(mean((ch4_bias_cal90central)^2, na.rm=T))
  )
  

#There are some atmospheric vials with likely erroneous identification, check values of neighbouring vials to try and fix it. 


all_with_identity<- conc_all %>% 
  arrange(gas, order) %>% 
  mutate(exetainer_ID=sample_id) %>% 
  select(-c(uniqueareaid, intercept, slope)) %>% 
  pivot_wider(names_from = gas, values_from = c(area, ppm, ppmstd)) %>% 
  full_join(all_exetainer_field, by = "exetainer_ID") %>% 
  rename(exe_type= `headspace/trapped/atm`) %>% 
  filter(!grepl("^S4", exetainer_ID)) %>% #remove season 4 (No gc_data)
  filter(!grepl("t", exetainer_ID)) %>%  #remove cores
  filter(!exetainer_ID%in%c("S1-CU-A1-1","S1-CU-A1-2","S1-CU-A1-3","S1-CU-A1-4","S1-CU-A1-5","S1-CU-A1-6","S1-CU-A1-7","S1-CU-A1-8","S1-CU-A1-9","S1-CU-A1-10",
                            "S1-RI-A2-18","S2-RI-A1-11a","S2-RI-P2-22b","S2-RI-P2-24b","S3-DA-A1-4b","S1-DU-A1-3","S1-DU-A1-3discard")) #Remove erroneous data (already checked, this should not be included)


all_samples_with_identity<-  all_with_identity %>% 
  filter(!grepl("^P", sample_id))


#Here samples and gcdata that do not match (duplicate labeling)
samples_without_gcdata<-all_samples_with_identity %>% 
  filter(is.na(order))

gcdata_without_sample<- all_samples_with_identity %>%
  filter(is.na(exe_type))


#Save ppm for all exetainers with their correspondece (exetainer_ID vs identity sample)

allsamples_ppm<- all_samples_with_identity %>% 
  select(exetainer_ID, pilot_site, subsite, campaign, subsite_ID, `date (dd.mm.yyyy)`, person_sampling, plot_ID, Strata, `transparent or dark`, exe_type, comment, ppm_co2, ppm_ch4, obs)

names(all_samples_with_identity)

write.csv(allsamples_ppm, file = paste0(gc_path, "allsamples_ppm_noblankcorrection.csv"),row.names = F)



#WE CAN OMIT EVERYTHING AFTER THIS: 
# use allsamples_ppm, this is the best-calibration possible. 







#OPTION B: blank correct-----
#WE CAN OMIT THIS: already tested, not improved bias. 

#Calculate ppm using the average (trimmed) blank area to correct the data before aplying a cero-intercept daily calibration.

#Use already calculated mean-blanks
co2_mean_blankarea
ch4_mean_blankarea


#2. Prepare caldata-----
gcsamples_bc<- gcsamples %>% 
  mutate(co2=co2-co2_mean_blankarea,
         ch4=ch4-ch4_mean_blankarea)



##Get calvials------
calvials_bc<- gcsamples_bc %>% 
  filter(grepl("^P",sample_id)) %>% #get only standards
  filter(batch==0) %>%  #select only standards from cal-curves
  filter(!is.na(vol)) %>%  #remove Pmix standards
  # filter(vol!=11.2) %>% #remove P20 vials, not reliable
  pivot_longer(cols=c("co2","ch4"), names_to = "gas", values_to = "area") %>% #pivot longer by gas
  mutate(uniqueareaid=paste(uniqueinj_id, gas,sep="_")) %>% #create separate injid per gas
  mutate(ppmstd=case_when(gas=="co2"~3000*(vol/11.2),
                          gas=="ch4"~15.29*(vol/11.2)))


## Get stdvials (internal std)-----

#Inter-batch internal standards (P20, P5, some P2)
stdvials_bc<- gcsamples_bc %>% 
  filter(grepl("^P",sample_id)) %>% #get only standards
  filter(batch!=0) %>%  #Remove standards from calcurves
  filter(!is.na(vol)) %>%  #remove Pmix standards
  # filter(vol!=11.2) %>% #remove P20 vials, not reliable
  pivot_longer(cols=c("co2","ch4"), names_to = "gas", values_to = "area") %>% #pivot longer by gas
  mutate(uniqueareaid=paste(uniqueinj_id, gas,sep="_")) %>% #create separate injid per gas
  mutate(ppmstd=case_when(gas=="co2"~3000*(vol/11.2),
                          gas=="ch4"~15.29*(vol/11.2)))


##Get atmvials (external std)----

#Atmosphere exetainers with ppm imputed from in-situ data (between incubations)
atm_exetainers<- all_exetainer_field %>% 
  filter(`headspace/trapped/atm`=="atmosphere") %>% 
  select(subsite_ID, exetainer_ID) %>% 
  left_join(subsite_atm, by="subsite_ID") %>% 
  pivot_longer(cols = c(co2, ch4), names_to = "gas", values_to = "ppmstd") %>% 
  select(-subsite_ID) %>% 
  rename(sample_id=exetainer_ID)

head(atm_exetainers)

#Joined with GC areas based on sample_id
atmvials_bc<- gcsamples_bc %>% 
  filter(sample_id%in%atm_exetainers$sample_id) %>% 
  pivot_longer(cols=c("co2","ch4"), names_to = "gas", values_to = "area") %>% #pivot longer by gas
  mutate(uniqueareaid=paste(uniqueinj_id, gas,sep="_")) %>% #create separate injid per gas
  left_join(atm_exetainers, by=c("sample_id", "gas"))

head(atmvials_bc)

rm(atm_exetainers)


#3. Inspect calibrations------

#Blank evolution (area co2 and ch4): already blank-corrected

calvials_bc %>% rbind(stdvials_bc) %>% 
  filter(vol==0) %>% 
  mutate(caldate=as.Date(as.character(yearmonthday),format="%Y%m%d")) %>% 
  ggplot(aes(x=caldate, y=area,  col = uniqueareaid%in%calvials$uniqueareaid))+
  geom_boxplot(outliers = F)+
  geom_point()+
  facet_wrap(facets=.~gas, scales="free")


##Outlier areas-----
#Mannually write down unreliable areas for each gas from calvials and stdvials

co2_out<- c("20240109_1_6_co2","20240109_1_7_co2","20240110_1_7_co2","20240221_1_7_co2","20240222_1_7_co2","20240311_1_8_co2","20240402_1_8_co2","20240321_1_8_co2","20240402_1_2_co2","20240402_1_3_co2","20240404_1_3_co2","20240411_1_1_co2","20240411_1_2_co2","20240412_1_5_co2","20240415_1_3_co2","20240415_1_5_co2","20240424_1_8_co2","20240426_1_8_co2","20240502_1_8_co2","20240503_1_8_co2","20240508_1_8_co2","20240510_1_8_co2","20240521_1_1_co2","20240521_1_8_co2","20240522_1_1_co2","20240522_1_8_co2","20240523_1_1_co2","20240530_1_1_co2","20240610_1_8_co2","20240611_1_2_co2","20240703_1_4_co2","20240729_1_1_co2","20240729_1_3_co2","20240807_1_3_co2","20240809_1_3_co2","20240813_1_3_co2","20240814_1_2_co2","20240814_1_9_co2")
ch4_out<- c("20240109_1_6_ch4","20240109_1_7_ch4","20240110_1_7_ch4","20240221_1_7_ch4","20240222_1_7_ch4","20240311_1_8_ch4","20240402_1_8_ch4","20240402_1_2_ch4","20240402_1_3_ch4","20240321_1_8_ch4","20240404_1_3_ch4","20240405_1_7_ch4","20240410_1_8_ch4","20240411_1_1_ch4","20240411_1_2_ch4","20240415_1_3_ch4","20240415_1_5_ch4","20240424_1_8_ch4","20240426_1_8_ch4","20240502_1_8_ch4","20240503_1_8_ch4","20240506_1_8_ch4","20240507_1_8_ch4","20240508_1_8_ch4","20240507_3_12_ch4","20240510_1_8_ch4","20240513_1_3_ch4","20240521_1_1_ch4","20240521_1_8_ch4","20240531_4_12_ch4","20240522_1_1_ch4","20240522_1_8_ch4","20240523_1_1_ch4","20240530_1_1_ch4","20240531_1_8_ch4","20240610_1_8_ch4","20240611_1_2_ch4","20240612_1_8_ch4","20240613_1_8_ch4","20240703_1_4_ch4","20240729_1_1_ch4","20240729_1_3_ch4","20240807_1_3_ch4","20240809_1_3_ch4","20240813_1_3_ch4","20240814_1_2_ch4","20240814_1_9_ch4")

##Get cals-----

#Calculate all the cal curves parameters (one per each dayofanalysis), forced through intercept, without P0s
cal_bc<- calvials_bc %>% 
  filter(!(uniqueareaid)%in%c(co2_out,ch4_out)) %>% 
  filter(vol!=0) %>% 
  select(uniqueinj_id,uniqueareaid,yearmonthday,ppmstd,gas,area)%>%
  nest_by(yearmonthday,gas) %>% 
  mutate(mod = list(lm(area~ppmstd+0, data=data))) %>%
  reframe(tidy(mod),glance(mod))


##Plot cals-----


#Plot factors 
cal_bc%>%
  filter(term=="ppmstd")%>%
  mutate(caldate=as.Date(as.character(yearmonthday),format="%Y%m%d")) %>% 
  # filter(r.squared>0.95) %>% 
  ggplot(aes(x=caldate, y=estimate, col=gas))+
  geom_point()+
  geom_point(data=data.frame(caldate=as.Date(as.character(c(20240502)),format="%Y%m%d"),
                             estimate=c(0,0),
                             gas=c("ch4","co2")))+
  facet_wrap(facets=.~gas, scales = "free")



##Plot R2
cal_bc%>%
  filter(term=="ppmstd")%>%
  mutate(caldate=as.Date(as.character(yearmonthday),format="%Y%m%d")) %>% 
  ggplot(aes(x=caldate, y=r.squared, col=gas))+
  geom_point()


#PLot every calcurve (after excluding outliers explictly, specified in "gas_out" vectors)

#Remove objects created in loop (for repeated executions of loop) 
rm(calcurve,cal_out,int_std,ext_std,dummydata,plt_CH4,plt_CO2,plt_N2O,datecalcurve)


max_co2_bc<- calvials_bc %>% 
  filter(!(uniqueareaid)%in%c(co2_out,ch4_out)) %>%
  filter(gas=="co2")%>% 
  summarise(max_area=max(area, na.rm = T)) %>% pull(max_area)

max_ch4_bc<- calvials_bc %>% 
  filter(!(uniqueareaid)%in%c(co2_out,ch4_out)) %>%
  filter(gas=="ch4")%>% 
  summarise(max_area=max(area, na.rm = T)) %>% pull(max_area)


#create vector to store calcurve plots
plt_list <- vector('list', length(unique(calvials_bc$yearmonthday)))

#Loop over each date and plot the calcurve for all 2 ghgs
for(datecalcurve in unique(calvials_bc$yearmonthday)){
  
  calcurve<- calvials_bc %>% 
    filter(yearmonthday==datecalcurve) %>% 
    filter(!(uniqueareaid)%in%c(co2_out,ch4_out))
  
  cal_out<- calvials_bc %>% 
    filter(yearmonthday==datecalcurve) %>% 
    filter(!(uniqueareaid%in%calcurve$uniqueareaid))
  
  int_std<- stdvials_bc %>% 
    filter(yearmonthday==datecalcurve) %>% 
    filter(!uniqueareaid%in%c(co2_out, ch4_out))
  
  ext_std<- atmvials_bc %>% 
    filter(yearmonthday==datecalcurve) %>% 
    filter(!uniqueareaid%in%c(co2_out, ch4_out))
  
  dummydata<- data.frame(ppmstd=c(0), area=as.numeric(NA_real_))
  
  plt_CO2 <- ggplot(subset(calcurve, gas=="co2"), aes(ppmstd, area))+
    geom_point()+
    # geom_text(aes(label=uniqueareaid))+
    stat_poly_line(formula = y~x+0) +
    stat_poly_eq(formula = y~x+0, use_label(c("eq", "adj.R2", "n")))+
    geom_point(data = subset(cal_out, gas=="co2"), aes(ppmstd,area, col="cal_out"))+
    geom_point(data = subset(int_std, gas=="co2"), aes(ppmstd,area, col="int_std"))+
    geom_point(data = subset(ext_std, gas=="co2"), aes(ppmstd,area, col="ext_std"))+
    # geom_text(data=subset(qual_inj, gas=="co2"), aes(vol, area,label=uniqueareaid))+
    geom_point(data=dummydata,aes(ppmstd,area, col = "cal_out"))+
    geom_point(data=dummydata,aes(ppmstd,area, col = "int_std"))+
    geom_point(data=dummydata,aes(ppmstd,area, col = "ext_std"))+
    scale_y_continuous(limits = c(0,max_co2))+
    theme_article()+
    labs(col = "") +
    theme(legend.position = "inside", legend.position.inside = c(0.8,0.1))+
    xlab("Std gas conc (ppm)")+
    ylab("CO2 Area (FID), blank-corrected")+
    ggtitle(paste0(datecalcurve," cal. plots"))
  
  plt_CH4 <- ggplot(subset(calcurve, gas=="ch4"), aes(ppmstd, area))+
    geom_point()+
    # geom_text(aes(label=uniqueareaid))+
    stat_poly_line(formula = y~x+0) +
    stat_poly_eq(formula = y~x+0, use_label(c("eq", "adj.R2", "n")))+
    geom_point(data = subset(cal_out, gas=="ch4"), aes(ppmstd,area, col="cal_out"))+
    geom_point(data = subset(int_std, gas=="ch4"), aes(ppmstd,area, col="int_std"))+
    geom_point(data = subset(ext_std, gas=="ch4"), aes(ppmstd,area, col="ext_std"))+
    # geom_text(data=subset(int_std, gas=="ch4"), aes(ppmstd, area,label=uniqueareaid))+
    geom_point(data=dummydata,aes(ppmstd,area, col = "cal_out"))+
    geom_point(data=dummydata,aes(ppmstd,area, col = "int_std"))+
    geom_point(data=dummydata,aes(ppmstd,area, col = "ext_std"))+
    scale_y_continuous(limits = c(0,max_ch4))+
    theme_article()+
    labs(col = "") +
    theme(legend.position = "inside", legend.position.inside = c(0.8,0.1))+
    xlab("Std gas conc (ppm)")+
    ylab("CH4 Area (FID), blank-corrected")
  
  
  plt_list[[which(unique(calvials_bc$yearmonthday)==datecalcurve)]] <- ggarrange(plt_CO2, plt_CH4, ncol = 2,nrow = 1)
  
}

# Print pdf: DO NOT PRINT PDF, not useful
# gg_save_pdf(list = plt_list, filename = paste0(gc_path,"calcurvesGC_co2&ch4_R4Cs_blank-corrected_origin-forzed.pdf"))


rm(calcurve,cal_out,int_std,ext_std,dummydata,plt_CH4,plt_CO2,datecalcurve)




#4. Calculate ppm -----
#for all samples

calforall_bc<- cal_bc %>% 
  select(yearmonthday,gas,term,estimate) %>% 
  pivot_wider(names_from = term, values_from = estimate) %>% 
  rename(slope=ppmstd)


conc_all_bc<- gcsamples_bc %>% 
  pivot_longer(cols=c("co2","ch4"), names_to = "gas", values_to = "area") %>% #pivot longer by gas
  mutate(uniqueareaid=paste(uniqueinj_id, gas,sep="_")) %>% #create separate injid per gas
  left_join(calforall_bc, by = c("gas", "yearmonthday")) %>% 
  filter(!(uniqueareaid)%in%c(co2_out,ch4_out)) %>% 
  mutate(ppm=(area)/slope) %>% 
  mutate(ppmstd=case_when(gas=="co2"~3000*(vol/11.2),
                          gas=="ch4"~15.29*(vol/11.2)))


#Check intercept vs blanks: 
conc_all_bc %>% 
  filter(vol==0) %>% 
  mutate(dateofanalysis=as.Date(as.character(yearmonthday), format = "%Y%m%d")) %>% 
  ggplot(aes(x=dateofanalysis, y=ppm))+
  geom_point()+
  facet_wrap(facets=.~gas, scales="free")




conc_all_bc %>% 
  filter(gas=="co2") %>% 
  ggplot(aes(x=ppmstd, y=ppm))+
  geom_point()+
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "n")))+
  scale_y_continuous(limits = c(-1000,4000))+
  scale_x_continuous(limits = c(-1000,4000))


conc_all_bc %>% 
  filter(gas=="ch4") %>% 
  ggplot(aes(x=ppmstd, y=ppm))+
  geom_point()+
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "n")))+
  scale_y_continuous(limits = c(-2,20))+
  scale_x_continuous(limits = c(-2,20))

#Check atm 1:1 ppm (in situ vs GC-derived)

atmvials_bc %>%
  left_join(calforall_bc, by = c("gas", "yearmonthday")) %>% 
  filter(!(uniqueareaid)%in%c(co2_out,ch4_out)) %>% 
  mutate(ppm=(area)/slope) %>% 
  filter(gas=="co2") %>% 
  ggplot(aes(x=ppmstd, y=ppm))+
  geom_point()+
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "n")))


atmvials_bc %>%
  left_join(calforall_bc, by = c("gas", "yearmonthday")) %>% 
  filter(!(uniqueareaid)%in%c(co2_out,ch4_out)) %>% 
  mutate(ppm=(area)/slope) %>% 
  filter(gas=="ch4") %>% 
  ggplot(aes(x=ppmstd, y=ppm))+
  geom_point()+
  scale_y_continuous(limits = c(1,7))+
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "n")))


atmvials_bc %>%
  left_join(calforall_bc, by = c("gas", "yearmonthday")) %>% 
  filter(!(uniqueareaid)%in%c(co2_out,ch4_out)) %>% 
  mutate(ppm=(area)/slope) %>% 
  select(-c(uniqueareaid,area,slope,ppmstd)) %>% 
  pivot_wider(names_from = gas, values_from = ppm,names_prefix = "ppm_") %>% 
  ggplot(aes(x=ppm_co2, y=ppm_ch4))+
  geom_point()+
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "n")))


#Check discrepancy ppmstd vs ppm
bias_cal_atm_bc<- atmvials_bc %>% 
  left_join(calforall_bc, by = c("gas", "yearmonthday")) %>% 
  filter(!(uniqueareaid)%in%c(co2_out,ch4_out)) %>% 
  mutate(ppm=(area)/slope) %>% 
  mutate(bias_cal=ppm-ppmstd) %>% 
  select(sample_id,bias_cal, gas) %>% 
  pivot_wider(names_from = gas, values_from = bias_cal) %>% 
  #Calculate central 90% distribution of co2 and ch4 bias
  mutate(co2_bias_cal90central=case_when(co2<quantile(co2,0.05,na.rm = T)~NA_real_,
                                         co2>quantile(co2,0.95,na.rm = T)~NA_real_,
                                         TRUE~co2)) %>% 
  
  mutate(ch4_bias_cal90central=case_when(ch4<quantile(ch4,0.05,na.rm = T)~NA_real_,
                                         ch4>quantile(ch4,0.95,na.rm = T)~NA_real_,
                                         TRUE~ch4)) %>% 
  mutate(co2_mean90calbilas=mean(co2_bias_cal90central, na.rm=T),
         ch4_mean90calbilas=mean(ch4_bias_cal90central, na.rm=T),
         co2_rmse_calbias= sqrt(mean((co2_bias_cal90central)^2, na.rm=T)),
         ch4_rmse_calbias= sqrt(mean((ch4_bias_cal90central)^2, na.rm=T))
  )


#There are some atmospheric vials with likely erroneous identification, check values of neighbouring vials to try and fix it. 


all_with_identity_bc<- conc_all_bc %>% 
  arrange(gas, order) %>% 
  mutate(exetainer_ID=sample_id) %>% 
  select(-c(uniqueareaid, slope)) %>% 
  pivot_wider(names_from = gas, values_from = c(area, ppm, ppmstd)) %>% 
  full_join(all_exetainer_field, by = "exetainer_ID") %>% 
  rename(exe_type= `headspace/trapped/atm`) %>% 
  filter(!grepl("^S4", exetainer_ID)) %>% #remove season 4 (No gc_data)
  filter(!grepl("t", exetainer_ID)) %>%  #remove cores
  filter(!exetainer_ID%in%c("S1-CU-A1-1","S1-CU-A1-2","S1-CU-A1-3","S1-CU-A1-4","S1-CU-A1-5","S1-CU-A1-6","S1-CU-A1-7","S1-CU-A1-8","S1-CU-A1-9","S1-CU-A1-10",
                            "S1-RI-A2-18","S2-RI-A1-11a","S2-RI-P2-22b","S2-RI-P2-24b","S3-DA-A1-4b","S1-DU-A1-3")) #Remove erroneous data (already checked, this should not be included)


all_samples_with_identity_bc<-  all_with_identity_bc %>% 
  filter(!grepl("^P", sample_id))


#Here samples and gcdata that do not match (duplicate labeling)
samples_without_gcdata_bc<-all_samples_with_identity_bc %>% 
  filter(is.na(order))

gcdata_without_sample_bc<- all_samples_with_identity_bc %>%
  filter(is.na(exe_type))



#Save ppm for all exetainers with their correspondece (exetainer_ID vs identity sample)

allsamples_ppm_bc<- all_samples_with_identity_bc %>% 
  select(exetainer_ID, pilot_site, subsite, campaign, subsite_ID, `date (dd.mm.yyyy)`, person_sampling, plot_ID, Strata, `transparent or dark`, exe_type, comment, ppm_co2, ppm_ch4, obs)

#DO NOT SAVE, not useful. 
# write.csv(allsamples_ppm_bc, file = paste0(gc_path, "allsamples_ppm_blank-correction.csv"),row.names = F)


#EVALUATE A vs B-------

#Decide calibration approach based on bias (mean bias) and accuracy (rmse bias) 

#Bias with dailycal no blank-correction:
bias_cal_atm %>% select(co2_mean90calbilas, co2_rmse_calbias, ch4_mean90calbilas, ch4_rmse_calbias) %>% distinct()

#Bias with blank-corrected data, and cero-intercept calfactors:
bias_cal_atm_bc %>% select(co2_mean90calbilas, co2_rmse_calbias, ch4_mean90calbilas, ch4_rmse_calbias) %>% distinct()

#DECISSION: 
#Using average blank-correction and forcing through origin does not improve calibration bias. Therefore, we will use the dailycal of option A. 


#Inspect and fix----

#Check that all values have a calibration and identify atm outliers. 

#ANY yearmonthday with samples but with no calibration?
unique(gcsamples$yearmonthday)[which(!unique(gcsamples$yearmonthday)%in%unique(cal$yearmonthday))]


atm_head<- allsamples_ppm %>% 
  filter(exe_type!="air trapped")



atm_head %>% 
  ggplot(aes(x=subsite_ID, y=ppm_co2, col=exe_type))+
  geom_point()

atm_head %>% 
  ggplot(aes(x=subsite_ID, y=ppm_ch4, col=exe_type))+
  geom_point()




#Filter and save individually Atm and HS: exetainer ID, sample identity, co2, ch4, n2o=NA, method:"GC-derived"



#Filter and save individually atm and HS from Li-COR: exetainerID, sample identity, co2, ch4, n2o, method: Li-COR
#Calculate in a different script water GHG concentrations using jorge's functions (incorporating also the . 



