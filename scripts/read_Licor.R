



read_Licor <- function(file){
  message(paste0("reading ",file_to_read," with script for ",analyser," gas analyser"))
  filename <- paste(datapath,directory_analyser,file_to_read, sep = "/")
  data_raw <- read_lines(filename)
  prefix <- substr(data_raw, start = 1,stop = 5) # isolte first 5 characters for each line

  # find line corresponding to headers
  headers <- unlist(strsplit(data_raw[which(prefix == "DATAH")], "\t"))
  units <- unlist(strsplit(data_raw[which(prefix == "DATAU")], "\t"))

  data <- read.delim(filename, sep = "\t", header = F, skip = which(prefix == "DATAU"), na.strings = "nan")
  names(data) <- headers

  my_data <- data.frame(date = data$DATE,
                        UTCtime = data$TIME,
                        unixtime = data$SECONDS,
                        H2O = data$H2O,
                        CO2 = data$CO2,
                        CH4 = data$CH4/1000, #ppm
                        Press = data$CAVITY_P,
                        label = data$REMARK)

  return(my_data)
}
