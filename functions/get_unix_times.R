



get_unix_times <- function(mydate, mytime){
  if(sum(is.na(mydate)) > 0 | sum(is.na(mytime))>0){
    unix_time <- NA*seq_along(mydate)
    for(i in seq_along(mydate)){
      if(!is.na(mydate[i]) & !is.na(mytime[i])){
        timestamp <- paste(hour(mytime[i]),minute(mytime[i]),0,sep = ":")
        unix_time[i] <- as.numeric(as.POSIXct(paste(as.Date(mydate[i]),timestamp,
                                                 sep = " "), tz = 'UTC'))
      }
    }

  } else {
    timestamp <- paste(hour(mytime),minute(mytime),0,sep = ":")
    unix_time <- as.numeric(as.POSIXct(paste(as.Date(mydate),timestamp,
                                             sep = " "), tz = 'UTC'))
  }
  return(unix_time)
}
