# ---
# Authors: Camille Minaudo
# Project: "RESTORE4Cs"
# date: "Oct 2023"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---



read_GHG_fieldsheets <- function(myfieldsheets_list){

  # Read the first one to get the headers
  fieldsheet_temp <- readxl::read_xlsx(myfieldsheets_list[1],
                                       col_names = T)
  my_headers <- names(fieldsheet_temp)

  # Go through them all and keep the info
  isF <- T
  for (f in myfieldsheets_list){
    message(paste0("reading fieldsheet ",f))
    fieldsheet_temp <- readxl::read_xlsx(f, col_names = F, range = "A3:V50", n_max = 100, .name_repair = "unique_quiet")
    names(fieldsheet_temp) <- my_headers
    fieldsheet_temp <- fieldsheet_temp[!is.na(fieldsheet_temp$plot_id),]
    fieldsheet_temp$date <- as.Date( fieldsheet_temp$date, tryFormats = c("%d.%m.%Y", "%d/%m/%Y"))
    fieldsheet_temp$subsite <- gsub(pattern = "-Fieldsheet-GHG.xlsx",replacement = "",x = basename(f))

    if(isF){
      isF <- F
      fieldsheet <- fieldsheet_temp
    } else {
      fieldsheet <- rbind(fieldsheet, fieldsheet_temp)
    }
  }

  fieldsheet <- fieldsheet[!is.na(fieldsheet$longitude),]

  fieldsheet$water_depth <- as.numeric(fieldsheet$water_depth)

  fieldsheet$unix_start <- get_unix_times(mydate = fieldsheet$date, mytime = fieldsheet$start_time)
  fieldsheet$unix_stop <- get_unix_times(mydate = fieldsheet$date, mytime = fieldsheet$end_time)

  fieldsheet$start_time <- strftime(fieldsheet$start_time, format="%H:%M:%S", tz = 'utc')
  fieldsheet$end_time <- strftime(fieldsheet$end_time, format="%H:%M:%S", tz = 'utc')

  
  fieldsheet$initial_co2 <- as.numeric(fieldsheet$initial_co2)
  fieldsheet$final_co2 <- as.numeric(fieldsheet$final_co2)
  fieldsheet$initial_ch4 <- as.numeric(fieldsheet$initial_ch4)
  fieldsheet$final_ch4 <- as.numeric(fieldsheet$final_ch4)
  
  fieldsheet$uniqID <- tolower(paste(fieldsheet$subsite,
                                     fieldsheet$plot_id,
                                     substr(fieldsheet$strata, 1, 1),
                                     substr(fieldsheet$transparent_dark, 1, 1), 
                                     substr(fieldsheet$start_time, 1, 5),
                                     sep = "-"))
  
  return(fieldsheet)
}
