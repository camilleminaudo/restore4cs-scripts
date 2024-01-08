# ---
# Authors: Camille Minaudo
# Project: "RESTORE4Cs"
# date: "Oct 2023"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
# This script


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
library(GoFluxYourself)
require(dplyr)
require(purrr)
require(msm)

source(paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/get_unix_times.R"))
source(paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/read_GHG_fieldsheets.R"))
source(paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/flux.term.R"))
source(paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/LM.flux.R"))
source(paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/obs.win.R"))
source(paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/separate_ebullition_from_diffusion.R"))



# ---- Directories ----
dropbox_root <- "C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/data/RESTORE4Cs/ghg_flux_R_scripts_example"
datapath <- paste0(dropbox_root,"/processed_data")
fieldsheetpath <- paste0(dropbox_root,"/fieldsheets")
corrfieldsheetpath <- fieldsheetpath
loggerspath <- paste0(dropbox_root,"/raw_data/temperature_logger")
plots_path <- paste0(dropbox_root,"/fluxes_estimates")
results_path <- paste0(dropbox_root,"/fluxes_estimates")


doPlot <- F

# ---- List GHG chamber fieldsheets in Dropbox and read them ---
# list filenames
myfieldsheets_list <- list.files(fieldsheetpath, pattern = "Fieldsheet-GHG.xlsx", all.files = T, full.names = T, recursive = T)
# Read all fieldsheets and put them in a single dataframe
fieldsheet <- read_GHG_fieldsheets(myfieldsheets_list)
# selecting only the rows with Picarro measurements
# fieldsheet <- fieldsheet[fieldsheet$gas_analyzer == "Picarro",]

# ---- List GHG chamber corrected fieldsheets in Dropbox ---
# list filenames
myfieldsheets_corrected_list <- list.files(corrfieldsheetpath, pattern = "Fieldsheet-GHG_corrected.csv", all.files = T, full.names = T, recursive = T)
subsites_corrected_fs <- gsub(pattern = "-Fieldsheet-GHG_corrected.csv",replacement = "",x = basename(myfieldsheets_corrected_list))


# ---- function to save a list of plots into pdf file ----

gg_save_pdf = function(list, filename) {
  pdf(filename)
  for (p in list) {
    print(p)
  }
  dev.off()
  invisible(NULL)
}

# ---- Go through each incubation in fieldsheet and compute linear model for co2 and ch4 ----
subsites <- unique(fieldsheet$subsite)
isF_incub <- T
for (subsite in subsites){
  message("Now processing ",subsite)

  corresp_fs <- fieldsheet[fieldsheet$subsite == subsite,]
  gs <- first(corresp_fs$gas_analyzer)

  # check if there is a corresponding corrected fieldsheet for this subsite
  existsCorr_fs <- which(subsites_corrected_fs == subsite)
  if (length(existsCorr_fs) > 0){
    message("using exact incubation start and stop times")
    fs_corr <- read.csv(file = myfieldsheets_corrected_list[existsCorr_fs])
    corresp_fs$unix_start <- fs_corr$unix_start_corrected[match(corresp_fs$start_time, fs_corr$start_time)]
    corresp_fs$unix_stop <- fs_corr$unix_stop_corrected[match(corresp_fs$start_time, fs_corr$start_time)]
    corresp_fs <- corresp_fs[!is.na(corresp_fs$unix_start),]
  }

  gs_folder <- "Picarro"
  # if(gs == "LI-COR"){
  #   gs_folder <- "RAW Data Licor-7810"
  # } else if (gs == "Los Gatos"){
  #   gs_folder <- "RAW Data Los Gatos"
  # } else if (gs == "Picarro"){
  #   gs_folder <- "Picarro"
  # } else{
  #   warning("------> gas analyser not properly detected!")
  # }


  # read corresponding temperature logger file and keep initial temperature

  # --- Read corresponding Loggers data ----
  SN_logger_float <- first(corresp_fs$logger_floating_chamber)
  SN_logger_tube <- first(corresp_fs$logger_transparent_chamber)
  site_ID <- str_sub(subsite, start = 1, end = 5)

  # finding out if corresponding file exists, and its extension
  require(tools)
  dir_exists_loggdata <- dir.exists(paste0(loggerspath,"/",site_ID,"/"))
  if(dir_exists_loggdata){
    f <- list.files(paste0(loggerspath,"/",site_ID,"/"), full.names = T)
    r <- grep(pattern = ".hobo", x = f)
    if(length(r)>0){f <- f[-r]}
    r <- grep(pattern = ".txt", x = f)
    if(length(r)>0){f <- f[-r]}
    f_ext <- file_ext(f)
    i_f_float <- grep(pattern = SN_logger_float, x = f)[1]
    i_f_tube <- grep(pattern = SN_logger_tube, x = f)[1]
  }

  if(!is.na(SN_logger_float) & !is.na(i_f_float)){
    is_data_logger_float = T
    message("...reading corresponding temperature logger file for floating chamber")
    if(f_ext[i_f_float]=="xlsx"){
      data_logger_float <- readxl::read_xlsx(f[i_f_float],col_names = T)
    } else if(f_ext[i_f_float]=="csv"){
      data_logger_float <- read.csv(f[i_f_float], header = T)
    }
    data_logger_float <- data_logger_float[,seq(1,3)]
    names(data_logger_float) <- c("sn","datetime","temperature")
    data_logger_float$sn <- SN_logger_float
    if(is.character(data_logger_float$datetime)){
      data_logger_float$datetime <- as.POSIXct(data_logger_float$datetime, tz = 'utc', tryFormats = c("%m/%d/%y %r", "%d/%m/%Y %H:%M"))
    }
    data_logger_float$unixtime <- as.numeric(data_logger_float$datetime)
  } else {
    message("===> no data logger linked to the floating chamber!")
    is_data_logger_float = F}

  if(!is.na(SN_logger_tube) & !is.na(i_f_tube)){
    is_data_logger_tube = T
    message("...reading corresponding temperature logger file for tube chamber")
    if(f_ext[i_f_tube]=="xlsx"){
      data_logger_tube <- readxl::read_xlsx(f[i_f_tube],col_names = T)
    } else if(f_ext[i_f_tube]=="csv"){
      data_logger_tube <- read.csv(f[i_f_tube], header = T, fill = T)
    }
    data_logger_tube <- data_logger_tube[,seq(1,4)]
    names(data_logger_tube) <- c("sn","datetime","temperature","light")
    data_logger_tube$sn <- SN_logger_tube
    if(is.character(data_logger_tube$datetime)){
      if(length(grep(pattern = "AM", x = first(data_logger_tube$datetime)))>0 | length(grep(pattern = "PM", x = first(data_logger_tube$datetime)))>0){
        data_logger_tube$datetime <- as.POSIXct(data_logger_tube$datetime, tz = 'UTC', format = c("%m/%d/%y %r"))
      } else {
        data_logger_tube$datetime <- as.POSIXct(data_logger_tube$datetime, tz = 'UTC', format = "%m/%d/%Y %H:%M")
      }
    }
    data_logger_tube$unixtime <- as.numeric(data_logger_tube$datetime)
  } else {
    is_data_logger_tube = F
    message("===> no data logger linked to the tube chamber!")
  }

  path2data <- paste0(datapath,"/",gs_folder,"/RData/",subsite)
  if(dir.exists(path2data)){

    setwd(path2data)
    load(file = paste0("data_",subsite,".RData"))

    if(doPlot){plt_list <- vector('list', length(corresp_fs$plot_id))}

    auxfile <- NULL

    for (incub in seq_along(corresp_fs$plot_id)){
      my_incub <- mydata[as.numeric(mydata$POSIX.time)> corresp_fs$unix_start[incub] &
                           as.numeric(mydata$POSIX.time)< corresp_fs$unix_stop[incub],]
      my_incub <- my_incub[!is.na(my_incub$CO2dry_ppm),]

      if (dim(my_incub)[1]>0){

        # Compute Temperature, Area and Volume from fieldsheet info
        myTemp <- 15 # a default temperature to run the flux calculation...
        if (corresp_fs$chamber_type[incub] == "floating"){
          if(is_data_logger_float){
            myTemp <- median(approx(data_logger_float$unixtime, data_logger_float$temperature, xout = as.numeric(my_incub$POSIX.time))$y )
          }
          myArea = 14365.4439 # cm2
          myVtot = 115/2 # L
        } else if (corresp_fs$chamber_type[incub] == "tube"){
          if(is_data_logger_tube){
            myTemp <- median(approx(data_logger_tube$unixtime, data_logger_tube$temperature, xout = as.numeric(my_incub$POSIX.time))$y )
          }
          myArea = pi*12.1**2 # cm2
          myVtot = myArea*as.numeric(corresp_fs$chamber_height_cm[incub])*1e-3 # L
        } else {
          warning("chamber type not correct!")
        }
        myPcham <- 100.1 #kPa

        # --- Create auxfile table ----
        # An auxfile table, made of fieldsheet, adding important variables. The auxfile
        # requires start.time and UniqueID.
        # start.time must be in the format "%Y-%m-%d %H:%M:%S"

        UniqueID = paste(seq_along(corresp_fs$pilot_site),corresp_fs$strata,corresp_fs$transparent_dark,sep = "-")[incub]

        auxfile_tmp <- data.frame(subsite = subsite,
                                  UniqueID = gsub(" ", "", UniqueID, fixed = TRUE),
                                  gas_analiser = gs,
                                  start.time = as.POSIXct((corresp_fs$unix_start[incub]), tz = "UTC"),
                                  duration = (corresp_fs$unix_stop[incub]) - (corresp_fs$unix_start[incub]),
                                  water_depth = corresp_fs$water_depth[incub],
                                  Area = myArea,
                                  Vtot = myVtot,
                                  Tcham = myTemp,
                                  Pcham = myPcham,
                                  strata = corresp_fs$strata[incub],
                                  chamberType = corresp_fs$chamber_type[incub],
                                  lightCondition = corresp_fs$transparent_dark[incub])
        if(is.null(auxfile)){
          auxfile <- auxfile_tmp
        } else {
          auxfile <- rbind(auxfile, auxfile_tmp)
        }

        if(doPlot){
          plt_CO2 <- ggplot(my_incub, aes(POSIX.time, CO2dry_ppm))+geom_line()+
            theme_article()+
            xlab("time UTC")+
            ylab("CO2dry [ppm]")+
            ggtitle(paste0(subsite," plot ",
                           corresp_fs$plot_id[incub]," ",corresp_fs$strata[incub]," ",
                           corresp_fs$transparent_dark[incub], ", depth = ",corresp_fs$water_depth[incub], " cm"))
          plt_CH4 <- ggplot(my_incub, aes(POSIX.time, CH4dry_ppb))+geom_line()+
            theme_article()+
            xlab("time UTC")+
            ylab("CH4dry [ppm]")
          plt_H2O <- ggplot(my_incub, aes(POSIX.time, H2O_ppm))+geom_line()+
            theme_article()+
            xlab("time UTC")+
            ylab("H2O [ppm]")

          # plt <- ggarrange(plt_CO2, plt_CH4, plt_H2O, ncol = 1)
          plt_list[[incub]] <- ggarrange(plt_CO2, plt_CH4, plt_H2O, ncol = 1)
        }

      }
    }

    # draw a few incubations randomly
    draw <- sample(seq_along(auxfile$subsite), 10)

    # Define the measurements' window of observation
    myauxfile <- auxfile[draw,]
    mydata_ow <- obs.win(inputfile = mydata, auxfile = myauxfile,
                         obs.length = myauxfile$duration, shoulder = 2)




    # ----------- Compute fluxes after manual selection of CO2 data

    # Manually identify measurements by clicking on the start and end points
    mydata_manID <- lapply(seq_along(mydata_ow), click.peak.loop,
                           flux.unique = mydata_ow,
                           gastype = "CO2dry_ppm",
                           plot.lim = c(200,1000)) %>%
      map_df(., ~as.data.frame(.x))


    # Additional auxiliary data required for flux calculation.
    mydata_manID <- mydata_manID %>%
      left_join(auxfile %>% select(UniqueID, Area, Vtot, Tcham, Pcham))

    # Add instrument precision for each gas
    prec = c(3.5, 0.6, 0.4, 45, 45)
    mydata_manID <- mydata_manID %>%
      mutate(CO2_prec = prec[1], CH4_prec = prec[2], N2O_prec = prec[3],
             H2O_prec = prec[4])


    # Calculate fluxes
    CO2_results_manID <- goFlux(mydata_manID, "CO2dry_ppm")
    H2O_results_manID <- goFlux(mydata_manID, "H2O_ppm")
    CH4_results_manID <- goFlux(mydata_manID, "CH4dry_ppb")

    # Use best.flux to select the best flux estimates (LM or HM)
    # based on a list of criteria
    criteria <- c("g.factor", "kappa", "MDF", "R2", "SE.rel")

    CO2_flux_res_manID <- best.flux(CO2_results_manID, criteria)
    H2O_flux_res_manID <- best.flux(H2O_results_manID, criteria)
    CH4_flux_res_manID <- best.flux(CH4_results_manID, criteria)

    # CO2_flux_res_manID <- CO2_flux_res_manID %>%
    #   left_join(myauxfile %>% select(UniqueID, strata, chamberType, lightCondition))


    # Plots results
    # Make a list of plots of all measurements, for each gastype
    CO2_flux_plots <- flux.plot(CO2_flux_res_manID, mydata_manID, "CO2dry_ppm")
    H2O_flux_plots <- flux.plot(H2O_flux_res_manID, mydata_manID, "H2O_ppm")
    CH4_flux_plots <- flux.plot(CH4_flux_res_manID, mydata_manID, "CH4dry_ppb")

    # Print pdf
    setwd(plots_path)
    flux_plot.ls <- c(CO2_flux_plots, CH4_flux_plots, H2O_flux_plots)
    flux2pdf(flux_plot.ls, outfile = paste0(subsite,".pdf"))


    # ----------- Compute fluxes blindly without any manual selection

    mydata_auto <- mydata_manID

    for(id in unique(mydata_auto$UniqueID)){
      ind_id <- which(mydata_auto$UniqueID == id)
      mydata_auto$start.time_corr[ind_id] <- first(mydata_auto$start.time[ind_id])
      mydata_auto$end.time_corr[ind_id] <- first(mydata_auto$start.time[ind_id]) + first(mydata_auto$obs.length[ind_id])
      mydata_auto$Etime[ind_id] <- as.numeric(mydata_auto$POSIX.time[ind_id]) - as.numeric(mydata_auto$start.time_corr)[ind_id]
    }

    mydata_auto$obs.length_corr <- as.numeric(mydata_auto$end.time_corr)-as.numeric(mydata_auto$start.time_corr)

    # Calculate fluxes
    CO2_results_auto <- goFlux(dataframe = mydata_auto, gastype = "CO2dry_ppm")
    H2O_results_auto <- goFlux(mydata_auto, "H2O_ppm")
    CH4_results_auto <- goFlux(mydata_auto, "CH4dry_ppb")

    # Use best.flux to select the best flux estimates (LM or HM)
    # based on a list of criteria
    criteria <- c("g.factor", "kappa", "MDF", "R2", "SE.rel")

    CO2_flux_res_auto <- best.flux(CO2_results_auto, criteria)
    H2O_flux_res_auto <- best.flux(H2O_results_auto, criteria)
    CH4_flux_res_auto <- best.flux(CH4_results_auto, criteria)


    CO2_flux_res_auto$flux_method <- "Blind"
    CO2_flux_res_manID$flux_method <- "Expert"

    ggplot(data = rbind(CO2_flux_res_auto, CO2_flux_res_manID))+
      geom_abline(slope = 0,intercept = 0, color = 'lightgrey')+
      geom_segment(data = data.frame(UniqueID = CO2_flux_res_auto$UniqueID,
                                     meth1 = CO2_flux_res_auto$best.flux,
                                     meth2 = CO2_flux_res_manID$best.flux), aes(x=UniqueID, xend=UniqueID, y = meth1, yend = meth2), size=1, alpha = 0.5)+
      geom_point(aes(UniqueID, best.flux, colour = flux_method), size=4, alpha = 0.5)+
      ylab("CO2 flux [mmol/m2/s]")+
      theme_article()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_colour_viridis_d(begin = 0.1, end = 0.9, option = "F")



    CH4_flux_res_auto$flux_method <- "Blind"
    CH4_flux_res_manID$flux_method <- "Expert"

    ggplot(data = rbind(CH4_flux_res_auto, CH4_flux_res_manID))+
      geom_abline(slope = 0,intercept = 0, color = 'lightgrey')+
      geom_segment(data = data.frame(UniqueID = CH4_flux_res_auto$UniqueID,
                                     meth1 = CH4_flux_res_auto$best.flux,
                                     meth2 = CH4_flux_res_manID$best.flux), aes(x=UniqueID, xend=UniqueID, y = meth1, yend = meth2), size=1, alpha = 0.5)+
      geom_point(aes(UniqueID, best.flux, colour = flux_method), size=4, alpha = 0.5)+
      ylab("CH4 flux [nmol/m2/s]")+
      theme_article()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_colour_viridis_d(begin = 0.1, end = 0.9, option = "F")



    # ----------- Estimate CH4 diffusion/ebullition contributions
    # method 1: density of prob. of first derivative
    # method 2: manual selection of linear pattern in CH4

    # estimate ch4 diffusion and ebullition components---------| METHOD 1 |-----------
    CH4_res_meth1 <- CH4_flux_res_manID
    CH4_res_meth1$total_estimated <- NA
    CH4_res_meth1$ebullition <- NA
    CH4_res_meth1$diffusion <- NA


    for (i in draw){
      my_incub <- mydata[as.numeric(mydata$POSIX.time)> auxfile$start.time[i] &
                           as.numeric(mydata$POSIX.time)< auxfile$start.time[i]+auxfile$duration[i],]
      my_incub <- my_incub[!is.na(my_incub$CO2dry_ppm),]
      # calling dedicated function
      df_ebull <- separate_ebullition_from_diffusion(my_incub = my_incub)
      # computing fluxes
      H2O_mol = my_incub$H2O_ppm / (1000*1000)
      myfluxterm <- flux.term(auxfile$Vtot[i], auxfile$Pcham[i], auxfile$Area[i],
                              auxfile$Tcham[i], first(H2O_mol))
      CH4_flux_total <- df_ebull$delta_ch4/df_ebull$duration*myfluxterm # nmol/m2/s
      CH4_flux_diff <- df_ebull$avg_diff_slope*myfluxterm # nmol/m2/s
      CH4_flux_ebull <- CH4_flux_total - CH4_flux_diff

      CH4_res_meth1$total_estimated[which(CH4_res_meth1$UniqueID==auxfile$UniqueID[i])] <- CH4_flux_total
      CH4_res_meth1$ebullition[which(CH4_res_meth1$UniqueID==auxfile$UniqueID[i])] <- CH4_flux_ebull
      CH4_res_meth1$diffusion[which(CH4_res_meth1$UniqueID==auxfile$UniqueID[i])] <- CH4_flux_diff
    }




    # Estimate ch4 diffusion and ebullition components ---------| METHOD 2 |-----------


    # Manually identify diffusive (more or less linear) CH4 behaviors by clicking on the start and end points
    myCH4_diffusion <- lapply(seq_along(mydata_ow), click.peak.loop,
                              flux.unique = mydata_ow,
                              gastype = "CH4dry_ppb",
                              plot.lim = c(1900,max(fieldsheet$final_ch4)*1000)) %>%
      map_df(., ~as.data.frame(.x))


    myCH4_diffusion <- myCH4_diffusion %>%
      mutate(CO2_prec = prec[1], CH4_prec = prec[2], N2O_prec = prec[3],
             H2O_prec = prec[4])


    # Calculate fluxes for CH4
    CH4_results_diffusion <- goFlux(myCH4_diffusion, "CH4dry_ppb")

    # Use best.flux to select the best flux estimates (LM or HM)
    # based on a list of criteria
    criteria <- c("g.factor", "kappa", "MDF", "R2", "SE.rel")

    CH4_res_diff <- best.flux(CH4_results_diffusion, criteria)


    CH4_res_meth2 <- CH4_flux_res_manID
    CH4_res_meth2$total_estimated <- NA
    CH4_res_meth2$ebullition <- NA
    CH4_res_meth2$diffusion <- CH4_res_diff$best.flux


    # Estimating ebullition component
    for (id in unique(CH4_res_meth2$UniqueID)){
      i <- which(CH4_res_meth2$UniqueID == id)

      CH4_final <- CH4_flux_res_manID$Ct[i]
      CH4_initial <-  CH4_flux_res_manID$C0[i]
      incubation_time <- myauxfile$duration[which(myauxfile$UniqueID == id)]
      CH4_res_meth2$total_estimated[i] <- (CH4_final-CH4_initial)/incubation_time*CH4_flux_res_manID$flux.term[i] # nmol/m2/s
      CH4_res_meth2$ebullition[i] <- CH4_res_meth2$total_estimated[i]  - CH4_res_meth2$diffusion[i] # total flux - diffusive term
    }


    # compare method 1 and 2 to estimate ebullitive contribution
    CH4_res_meth1$flux_method <- "dydt"
    CH4_res_meth2$flux_method <- "Expert"

    ggplot(data = rbind(CH4_res_meth1, CH4_res_meth2))+
      geom_abline(slope = 0,intercept = 0, color = 'lightgrey')+
      geom_segment(data = data.frame(UniqueID = CH4_res_meth1$UniqueID,
                                     meth1 = CH4_res_meth1$ebullition,
                                     meth2 = CH4_res_meth2$ebullition), aes(x=UniqueID, xend=UniqueID, y = meth1, yend = meth2), size=1, alpha = 0.5)+
      geom_point(aes(UniqueID, ebullition, colour = flux_method), size=4, alpha=0.5)+
      ylab("ebullition component [nmol/m2/s]")+
      theme_article()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_colour_viridis_d(begin = 0.1, end = 0.9, option = "F")





    # saving fluxes estimates
    sv_co2 <- rbind(CO2_flux_res_auto[,c('UniqueID','best.flux','flux_method')], CO2_flux_res_manID[,c('UniqueID','best.flux','flux_method')])
    sv_co2$variable <- "CO2"
    sv_ch4 <- rbind(CH4_flux_res_auto[,c('UniqueID','best.flux','flux_method')], CH4_flux_res_manID[,c('UniqueID','best.flux','flux_method')])
    sv_ch4$variable <- "CH4"
    flux_co2_ch4 <- spread(rbind(sv_co2, sv_ch4), variable, best.flux)

    setwd(results_path)
    write.csv(x = flux_co2_ch4, file = paste0(subsite,"_co2_ch4_fluxes.csv"), row.names = F, col.names = T)





  } # closing if(dir.exists(path2data)){
}


#----- joining CO2 and CH4 fluxes estimates into a single table -----


table_results <- auxfile %>%
  left_join(CO2_flux_res %>% select(UniqueID, best.flux, model, quality.check)) %>%
  rename(CO2_flux = best.flux, CO2_model = model, CO2_quality.check = quality.check) %>%
  left_join(CH4_flux_res %>% select(UniqueID, best.flux, model, quality.check)) %>%
  rename(CH4_diffusive_flux = best.flux, CH4_model = model, CH4_quality.check = quality.check)


plt_CO2 <- ggplot(table_results, aes(lightCondition, CO2_flux, fill = lightCondition))+
  geom_hline(yintercept = 0)+
  geom_boxplot(alpha=0.2)+geom_jitter(width = 0.2)+
  theme_article()+facet_wrap(.~strata, scales = "free")+
  ylab("CO2 flux mmol/m2/s")+
  ggtitle(paste0(subsite,", CO2 flux"))+
  scale_fill_viridis_d(begin = 0.2, end = 0.9)


plt_CH4diff <- ggplot(table_results, aes(lightCondition, CH4_diffusive_flux, fill = lightCondition))+
  geom_hline(yintercept = 0)+
  geom_boxplot(alpha=0.2)+geom_jitter(width = 0.2)+
  theme_article()+facet_wrap(.~strata, scales = "free")+
  ylab("CH4 flux nmol/m2/s")+
  ggtitle(paste0(subsite,", CH4 diffusive flux"))+
  scale_fill_viridis_d(begin = 0.2, end = 0.9)


# plt_CH4ebull <- ggplot(table_results, aes(lightCondition, CH4_ebullition_flux, fill = lightCondition))+
#   geom_hline(yintercept = 0)+
#   geom_boxplot(alpha=0.2)+geom_jitter(width = 0.2)+
#   theme_article()+facet_wrap(.~strata, scales = "free")+
#   ylab("CH4 flux nmol/m2/s")+
#   ggtitle(paste0(subsite,", CH4 ebullition"))+
#   scale_fill_viridis_d(begin = 0.2, end = 0.9)

ggarrange(plt_CO2, plt_CH4diff, ncol = 1)













setwd(results_path)






