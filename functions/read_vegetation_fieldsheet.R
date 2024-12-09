# ---
# Authors: Miguel Cabrera
# Project: "RESTORE4Cs"
# date: "Dec 2024"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---


read_vegetation_fieldsheets <- function(myfieldsheets_list){
  # Read the first one to get the headers
  fieldsheet_temp <- readxl::read_xlsx(myfieldsheets_list[1],
                                       col_names = T)
  my_headers <- names(fieldsheet_temp)
  
  # Go through them all and keep the info
  isF <- T
  for (f in myfieldsheets_list){
    message(paste0("reading vegetation ",f))
    fieldsheet_temp <- readxl::read_xlsx(f, col_names = F, range = "A3:G150", .name_repair = "unique_quiet")
    names(fieldsheet_temp) <- my_headers
    fieldsheet_temp <- fieldsheet_temp[!is.na(fieldsheet_temp$season),]#Drop rows with NAs in season column
    
    if(isF){
      isF <- F
      fieldsheet <- fieldsheet_temp
    } else {
      fieldsheet <- rbind(fieldsheet, fieldsheet_temp)
    }
  }
  
  fieldsheet$dry_weight<- as.numeric(fieldsheet$dry_weight)
  fieldsheet$plot_id<- as.numeric(fieldsheet$plot_id)
  fieldsheet$subsite<- paste(fieldsheet$season, 
                             fieldsheet$pilot_site,
                             paste0(substr(fieldsheet$subsite,1,1),parse_number(fieldsheet$subsite)),sep = "-")
  fieldsheet$plotcode<- paste(fieldsheet$subsite, fieldsheet$plot_id, sep = "-")
  fieldsheet<-fieldsheet %>% 
    select(season,pilot_site,subsite, plot_id, plotcode,dry_weight, comments, person_measuring)
  
  return(fieldsheet)
}