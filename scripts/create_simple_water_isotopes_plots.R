
# ---
# Authors: Benjamin Misteli
# Project: "RESTORE4Cs"
# date: "Nov 2023"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
# This script loads the Stable Isotope data and creates an overview plot


rm(list = ls()) # clear workspace
cat("/014") # clear console


# ---- packages ----
library(SIBER)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggpubr)

# ---- Directories ----
dropbox_root <- "C:/Users/misteli/Dropbox/RESTORE4Cs - Fieldwork/Data/"

datapath <- paste0(dropbox_root,"/Water/Stable Isotopes")
setwd(datapath)
water_isotopes <- readxl::read_xlsx("Water_Isotopes_S1_LT.xlsx",
                                     col_names = T)

water_isotopes <- water_isotopes %>% mutate(conservation = factor(Conservation), 
                                        number = factor(Number),
                                        d2H = `δ2H, in ‰`, 
                                        d18O = `δ18O, in ‰`,
                                        .keep = "unused") 

unique_case_pilot <- unique(water_isotopes$`Case Pilot`)
plot_list <- list()

# Loop through each country code
for (case_pilot in unique_case_pilot) {

  
  # Filter data for the current country code
  country_data <- filter(water_isotopes, `Case Pilot` == case_pilot)
  
  # Create a plot for the current country
  plot <- ggplot(data = country_data, 
                 aes(x = d18O, 
                     y = d2H)) + 
    geom_point(aes(color = conservation, shape = number), size = 5) +
    ylab(expression(paste(delta^{2}, "H (permille)"))) +
    xlab(expression(paste(delta^{18}, "O (permille)"))) + 
    theme(text = element_text(size = 16)) + 
    scale_color_viridis_d() +
    ggtitle(paste("Case Pilot:", case_pilot))
  
  # Calculate summary statistics for the current country
  summary_data <- country_data %>% 
    group_by(conservation, number) %>% 
    summarise(count = n(),
              mC = mean(d18O), 
              sdC = sd(d18O),
              mN = mean(d2H), 
              sdN = sd(d2H))
  
  # Create plot for the current country
  second_plot <- plot +
    geom_errorbar(data = summary_data, 
                  mapping = aes(x = mC, y = mN,
                                ymin = mN - 1.96 * sdN, 
                                ymax = mN + 1.96 * sdN), 
                  width = 0) +
    geom_errorbarh(data = summary_data, 
                   mapping = aes(x = mC, y = mN,
                                 xmin = mC - 1.96 * sdC,
                                 xmax = mC + 1.96 * sdC),
                   height = 0) + 
    geom_point(data = summary_data, aes(x = mC, 
                                        y = mN,
                                        fill = conservation), 
               color = "black", shape = 22, size = 5,
               alpha = 0.7, show.legend = FALSE) + 
    scale_fill_viridis_d()+  
    geom_abline(intercept = 10.8, slope =  8.2, color = "red", linetype = "dashed")
  
  plot_list[[case_pilot]] <- second_plot
  
  # Save the second_plot as a PNG file
  filename <- paste("S1_stable_isotopes_", case_pilot, ".png", sep = "")
  ggsave(filename, second_plot, width = 10, height = 8, units = "in", dpi = 300)
}

combined_plots <- ggarrange(plotlist = plot_list, ncol = 2, nrow=3, common.legend = TRUE, legend = "right")
png(filename = "S1_combined_plots.png",
    width = 12, height = 10, units = "in", res=300)
combined_plots
dev.off()

water_isotopes_1 <- readxl::read_xlsx("Water_Isotopes_S1_LT.xlsx",
                                     col_names = T)

water_isotopes_1 <- water_isotopes_1 %>% mutate(CasePilot = factor(`Case Pilot`), 
                                            conservation = factor(Conservation),
                                            d2H = `δ2H, in ‰`, 
                                            d18O = `δ18O, in ‰`,
                                            .keep = "unused") 




plot <- ggplot(data = water_isotopes_1, 
               aes(x = d18O, 
                   y = d2H)) + 
  geom_point(aes(color = CasePilot, shape = conservation), size = 4) +
  ylab(expression(paste(delta^{2}, "H (permille)"))) +
  xlab(expression(paste(delta^{18}, "O (permille)"))) + 
  theme(text = element_text(size = 16)) + 
  scale_color_viridis_d() 

summary_data <- water_isotopes_1 %>% 
  group_by(CasePilot, conservation) %>% 
  summarise(count = n(),
            mC = mean(d18O), 
            sdC = sd(d18O), 
            mN = mean(d2H), 
            sdN = sd(d2H))

second_plot <- plot +
  geom_errorbar(data = summary_data, 
                mapping = aes(x = mC, y = mN,
                              ymin = mN - 1.96 * sdN, 
                              ymax = mN + 1.96 * sdN), 
                width = 0) +
  geom_errorbarh(data = summary_data, 
                 mapping = aes(x = mC, y = mN,
                               xmin = mC - 1.96 * sdC,
                               xmax = mC + 1.96 * sdC),
                 height = 0) + 
  geom_point(data = summary_data, aes(x = mC, 
                                      y = mN,
                                      fill = conservation), 
             color = "black", shape = 22, size = 5,
             alpha = 0.7, show.legend = FALSE) + 
  scale_fill_viridis_d()+  
  geom_abline(intercept = 10.8, slope =  8.2, color = "red", linetype = "dashed")

ggsave("S1_all_in_one.png", second_plot, width = 10, height = 8, units = "in", dpi = 300)
