
# Authors: Camille Minaudo
# Project: "RESTORE4Cs"
# date: "Dec 2023"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
# This script identify the most probable slopes observed in a given chamber incubation
# allowing to distinguish the diffusive flux from the ebullitive flux


source('C:/Projects/myGit/my_r_tools/scripts/get_dxdy.R')


separate_ebullition_from_diffusion <- function(my_incub){
  mych4 <- data.frame(time = as.numeric(my_incub$POSIX.time-my_incub$POSIX.time[1]),
                      ch4 = my_incub$CH4dry_ppb)
  mych4$ch4smooth <- smooth.spline(x = mych4$time, y = mych4$ch4, nknots = round(max(mych4$time)/3), spar = 0.6)$y

  # compute first derivative
  mych4$dydt <- get_dxdy(mych4$time, mych4$ch4smooth)

  # compute density probability of first derivative
  d <- density(mych4$dydt)
  half_dens_max <- max(d$y)/2
  ind_over_half_dens_max <- which(d$y>half_dens_max)

  # we define lower and upper boundaries of the "main slope" as the slopes with
  # a density of probability higher than half of the maximum density
  lower_bound <- d$x[first(ind_over_half_dens_max)]
  upper_bound <- d$x[last(ind_over_half_dens_max)]

  # plot(d)
  # abline(h=half_dens_max, col='grey')
  # abline(v=c(lower_bound, upper_bound), col='lightblue')

  avg_slope <- mean(d$x[ind_over_half_dens_max])
  sd_slope <- sd(d$x[ind_over_half_dens_max])

  plot(mych4$time, mych4$ch4, type='l')
  lines(mych4$time, mych4$ch4smooth, col='black')
  points(mych4$time[mych4$dydt>lower_bound & mych4$dydt<upper_bound], mych4$ch4smooth[mych4$dydt>lower_bound & mych4$dydt<upper_bound], col = 'red')
  abline(a = first(mych4$ch4smooth), b = avg_slope, col= 'red')
  abline(a = first(mych4$ch4smooth), b = avg_slope+sd_slope/2, col= 'grey')
  abline(a = first(mych4$ch4smooth), b = avg_slope-sd_slope/2, col= 'grey')

  delta_C <- last(mych4$ch4smooth) - min(mych4$ch4smooth)
  duration <- last(mych4$time)


  df_out <- data.frame(duration = last(mych4$time),
                       delta_ch4 = last(mych4$ch4smooth) - min(mych4$ch4smooth),
                       avg_diff_slope = avg_slope,
                       sd_diff_slope = sd_slope)

  return(df_out)
}







