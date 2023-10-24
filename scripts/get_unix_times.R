



get_unix_times <- function(mydate, mytime){
  timestamp <- paste(hour(mytime),minute(mytime),0,sep = ":")
  unix_time <- as.numeric(as.POSIXct(paste(as.Date(mydate),timestamp,
                                           sep = " "), tz = 'UTC'))
  return(unix_time)
}
