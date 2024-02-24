





join_auxfile_with_data <- function(flux.unique) { # flux.unique <- mydata_ow

  # Assign NULL to variables without binding
  flag <- . <- NULL

  gastype <- "CO2dry_ppm"
  # Extract data from data.frame
  flux.meas <- Reduce("c", flux.unique[, gastype])
  time.meas <- Reduce("c", flux.unique[, "POSIX.time"])



  flux.start <- min(time.meas)
  flux.start_corr <- flux.start

  c0 <- flux.meas[which.min(abs(time.meas-flux.start_corr))]
  flux.end <- max(time.meas)
  cf <- flux.meas[which.min(abs(time.meas-flux.end))]
  flux.flag <- which(between(time.meas, flux.start_corr, flux.end))
  flux.diff <- as.numeric(flux.start - flux.start_corr, units = "secs") -1

  # Based on these identifications, the flagging and Etime columns are added
  flux.corr <- flux.unique %>%
    # 0 == no measurement, 1 == measurement point to be used for flux calculation
    mutate(flag = 1) %>%
    # Set to 0 at start of measurement and count seconds to end of measurement
    mutate(Etime = seq(flux.diff, length(time.meas) + flux.diff - 1)) %>%
    # Add start.time_corr, end.time, initial and final concentrations
    mutate(start.time_corr = flux.start_corr,
           end.time = flux.end,
           c0 = c0,
           cf = cf)
  # Return results
  return(flux.corr)

}

join_auxfile_with_data.loop <- function(x, flux.unique) {

  # Function to apply in the loop. Adapt parameters to your needs.
  flux.corr <- join_auxfile_with_data(flux.unique[[x]])

  return(flux.corr)
}

