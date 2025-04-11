
# Authors: Miguel Cabrera
# Project: "RESTORE4Cs"
# date: "Feb 2025"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---


# Simplified function to detect constant slope for a specified duration in selected columns: 
#This function essentially provides a check for artefacts: constant concentration or constant slope for a given duration for any ghg, also checks for negative ghg. 
is_constant_slope_or_neg_any_ghg <- function(my_incub, POSIX.time, check_cols=c("CO2dry_ppm","CH4dry_ppb"), duration = 5) {
  
  # Ensure POSIX.time is present in the data
  if (!(POSIX.time %in% colnames(my_incub))) {
    stop(paste("Column", POSIX.time, "not found in the data"))
  }
  
  # Ensure check_cols are in the data
  if (!all(check_cols %in% colnames(my_incub))) {
    stop("One or more specified columns not found in the data")
  }
  
  # Extract the time column
  time <- as.numeric(my_incub$POSIX.time-my_incub$POSIX.time[1])
  
  # Iterate over each selected column to check
  for (col_name in check_cols) {
    # Calculate the slope (difference between consecutive values) for the current column
    slopes <- diff(my_incub[[col_name]])
    #Calculate minimum value for current column
    minimum<- min(my_incub[[col_name]],na.rm = T)
    
    #Return TRUE if negative GHG found
    if(minimum<0){return("negative values")} else{
    
    #IF all positive,  Check for constant slope in the current column for the specified duration
    for (i in 1:(length(slopes) - duration)) {
      if (all(slopes[i] == slopes[i:(i + duration - 1)])) {
        return("constant slope")  # Return TRUE as soon as we find a constant slope
      }
    }
  }
  }
  return("no artefact")  # Return FALSE if no constant slope is found and all values are positive
}


