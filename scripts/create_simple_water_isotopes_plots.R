
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
library(dplyr)
library(ggplot2)
library(ggpubr)
library(lme4)

# ---- Directories ----
dropbox_root <- "C:/Users/misteli/Dropbox/RESTORE4Cs - Fieldwork/Data/"

datapath <- paste0(dropbox_root,"/Water/Stable Isotopes")
setwd(datapath)
water_isotopes <- readxl::read_xlsx("Water_Isotopes_Dataset.xlsx",
                                     col_names = T)

# Remove rows with high tide
water_isotopes <- water_isotopes %>%
  filter(Tide != "HT")

# Convert isotope and environmental data to numeric
water_isotopes <- water_isotopes %>%
  mutate(`δ2H, in ‰` = as.numeric(`δ2H, in ‰`),
         `δ18O, in ‰` = as.numeric(`δ18O, in ‰`),
         Conductivity = as.numeric(Conductivity),
         Temperature = as.numeric(Temperature),
         pH = as.numeric(pH))

# Convert categorical variables to factors
water_isotopes <- water_isotopes %>%
  mutate(Season = factor(Season),
         Subsite = factor(Subsite),
         Conservation = factor(Conservation),
         Replicate = factor(Replicate))

# Calculate d-excess and add as a new column
water_isotopes <- water_isotopes %>%
  mutate(d_excess = `δ2H, in ‰` - 8 * `δ18O, in ‰`)

# Summary statistics of isotopes per Subsite and Season
isotope_summary <- water_isotopes %>%
  group_by(Subsite, Season, Conservation) %>%
  summarise(mean_δ2H = mean(`δ2H, in ‰`, na.rm = TRUE),
            sd_δ2H = sd(`δ2H, in ‰`, na.rm = TRUE),
            mean_δ18O = mean(`δ18O, in ‰`, na.rm = TRUE),
            sd_δ18O = sd(`δ18O, in ‰`, na.rm = TRUE))

# ---- Boxplot for δ²H with statistical comparison ----
p1 <- ggplot(water_isotopes, aes(x = Conservation, y = `δ2H, in ‰`, fill = Conservation)) +
  geom_boxplot() +
  facet_wrap(~ `Case Pilot`) +  # Create a separate plot for each site
  stat_compare_means(method = "anova", label = "p.format") +  # Add ANOVA p-value
  stat_compare_means(comparisons = list(c("P", "R"), c("P", "A"), c("R", "A")), 
                     method = "t.test", label = "p.signif") +  # Add pairwise comparisons
  labs(title = "δ²H by Conservation Type with Statistical Comparison", 
       x = "Conservation Type", y = "δ²H (‰)")
ggsave("Boxplot_δ2H_by_Conservation.png", plot = p1, width = 10, height = 6, dpi = 300)

# ---- Boxplot for δ¹⁸O with statistical comparison ----
p2 <- ggplot(water_isotopes, aes(x = Conservation, y = `δ18O, in ‰`, fill = Conservation)) +
  geom_boxplot() +
  facet_wrap(~ `Case Pilot`) +  # Create a separate plot for each site
  stat_compare_means(method = "anova", label = "p.format") +  # Add ANOVA p-value
  stat_compare_means(comparisons = list(c("P", "R"), c("P", "A"), c("R", "A")), 
                     method = "t.test", label = "p.signif") +  # Add pairwise comparisons
  labs(title = "δ¹⁸O by Conservation Type with Statistical Comparison", 
       x = "Conservation Type", y = "δ¹⁸O (‰)")
ggsave("Boxplot_δ18O_by_Conservation.png", plot = p2, width = 10, height = 6, dpi = 300)

# ---- Boxplot of d-excess by Conservation Type ----
p3 <- ggplot(water_isotopes, aes(x = Conservation, y = d_excess, fill = Conservation)) +
  geom_boxplot() +
  facet_wrap(~ `Case Pilot`) +  # Separate plots for each site
  stat_compare_means(method = "anova", label = "p.format") +  # Add ANOVA p-value for overall comparison
  stat_compare_means(comparisons = list(c("P", "R"), c("P", "A"), c("R", "A")), 
                     method = "t.test", label = "p.signif") +  # Add pairwise comparisons
  labs(title = "d-excess by Conservation Type with Statistical Comparison", 
       x = "Conservation Type", y = "d-excess (‰)")
ggsave("Boxplot_d-excess_by_Conservation.png", plot = p3, width = 10, height = 6, dpi = 300)

# ---- Scatter plot of δ²H vs. δ¹⁸O by Conservation Type ----
p4 <- ggplot(water_isotopes, aes(x = `δ18O, in ‰`, y = `δ2H, in ‰`, color = Conservation)) +
  geom_point() +
  facet_wrap(~ `Case Pilot`, scales = "free") +  # Allow each site to have its own axis scales
  labs(title = "δ²H vs. δ¹⁸O by Conservation Type within Each Site", 
       x = "δ¹⁸O (‰)", 
       y = "δ²H (‰)") +
  geom_smooth(method = "lm", se = FALSE)  # Add linear regression lines without confidence intervals
ggsave("Scatter_δ2H_vs_δ18O_by_Conservation.png", plot = p4, width = 10, height = 6, dpi = 300)

# ---- Local Meteoric Water Line (LMWL) vs GMWL ----
p5 <- ggplot(water_isotopes, aes(x = `δ18O, in ‰`, y = `δ2H, in ‰`, color = Conservation)) +
  geom_point() +
  geom_abline(intercept = 10, slope = 8, color = "blue", linetype = "dashed") +  # GMWL
  facet_wrap(~ `Case Pilot`, scales = "free") +  
  labs(title = "Local Meteoric Water Line (LMWL) vs GMWL", 
       x = "δ¹⁸O (‰)", 
       y = "δ²H (‰)") +
  geom_smooth(method = "lm", se = FALSE)  # LMWL (local line)
ggsave("LMWL_vs_GMWL.png", plot = p5, width = 10, height = 6, dpi = 300)
