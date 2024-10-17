
# Authors: Camille Minaudo
# Project: "RESTORE4Cs"
# date: "Dec 2023"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
# This script identify the most probable slopes observed in a given chamber incubation
# allowing to distinguish the diffusive flux from the ebullitive flux


separate_ebullition_from_diffusion <- function(my_incub, UniqueID, doPlot){
  mych4 <- data.frame(time = as.numeric(my_incub$POSIX.time-my_incub$POSIX.time[1]),
                      ch4 = my_incub$CH4dry_ppb)
  
  ch4smooth <- smooth.spline(x = mych4$time, y = mych4$ch4, nknots = round(length(mych4$time)/3), spar = 0.6)
  
  mych4$ch4smooth <- approx(ch4smooth$x, ch4smooth$y, xout = mych4$time, rule = 2)$y 
  
  # compute first derivative
  mych4$dydt <- get_dxdy(mych4$time, mych4$ch4smooth)
  
  # compute density probability of first derivative
  d <- density(mych4$dydt)
  half_dens_max <- 0.5*max(d$y)
  ind_over_half_dens_max <- which(d$y>half_dens_max)
  
  # we define lower and upper boundaries of the "main slope" as the slopes with
  # a density of probability higher than half of the maximum density
  lower_bound <- d$x[first(ind_over_half_dens_max)]
  upper_bound <- d$x[last(ind_over_half_dens_max)]
  
  # most_probable_slope <- d$x[which.max(d$y)]
  
  avg_slope <- mean(d$x[ind_over_half_dens_max])
  sd_slope <- sd(d$x[ind_over_half_dens_max])
  
  
  mych4_sel <- mych4[mych4$dydt>lower_bound & mych4$dydt<upper_bound,]
  # my_lm <- lm(data = mych4_sel, formula = ch4smooth~time)
  # my_c0 <- coefficients(my_lm)[1]
  my_c0 <- min(mych4$ch4smooth[seq(1,100)])
  
  if(doPlot){
    
    p_density <- ggplot(mych4, aes(dydt))+
      geom_rect(aes(xmin = lower_bound, ymin = -Inf, xmax = upper_bound, ymax = Inf), fill = "#FF7F50", alpha=0.1)+
      geom_hline(yintercept = half_dens_max, alpha=0.5)+
      geom_density()+
      xlab("Rate of change of CH4 [ppb/secs]")+
      theme_article()
    
    p_fit <- ggplot(mych4)+
      geom_path(aes(time, ch4, colour = "raw"), linewidth=0.8)+
      geom_path(aes(time, ch4smooth, colour = "smooth"), linewidth = 1)+
      geom_point(data = mych4_sel, aes(time, ch4smooth, colour = "selected"), size=2, alpha=0.5)+
      geom_abline(slope = avg_slope, intercept = my_c0, color = "red")+
      geom_abline(slope = avg_slope+sd_slope/2, intercept = my_c0, color = "grey70")+
      geom_abline(slope = avg_slope-sd_slope/2, intercept = my_c0, color = "grey70")+
      xlab("Elapsed time [secs]")+
      ylab("CH4 [ppb]")+
      scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
      # scale_color_manual(values=c("#999999", "#56B4E9"))+
      ggtitle(UniqueID)+theme_article()+theme(legend.title = element_blank())+theme(legend.position = c(0.8,0.2))
    
    ggarrange(p_density, p_fit, widths = c(0.3,0.7))
    
  }
  
  delta_C <- last(mych4$ch4smooth) - min(mych4$ch4smooth)
  duration <- last(mych4$time)
  
  
  df_out <- data.frame(duration = last(mych4$time),
                       delta_ch4 = last(mych4$ch4smooth) - min(mych4$ch4smooth),
                       avg_diff_slope = avg_slope,
                       sd_diff_slope = sd_slope)
  
  return(df_out)
}







