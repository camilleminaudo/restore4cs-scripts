#CORES status comparisons and stats




# ---
# Authors: MIGUEL CABRERA
# Project: "RESTORE4Cs"
# date: "April 2025"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of script ----
# This script computes per-area daily fluxes of CO2, CH4 and N2O for all core incubations. As well as the CH4 stock in the core. 


#----Decissions to make: ------

#1. How to treat non-significant fluxes? 
#There are fluxes that are non-significant due to high uncertainty in the initial and final concentration estimates. Although non-significant (p>0.05), they might be of non-trivial magnitude, potentially obscuring status-comparison results (CHECK). Including them will help with N harmonization... Options: A_exclude non-significant fluxes, B_include as-is, C_impute flux? Preferred option: Include non-significant fluxes as-is. Even though they are not significant, excluding them will bias the N. 

#2. Statistical treatment
#What is the most-appropriate statistical treatment to answer the question: does restoration have an effect on GHG fluxes?
#CV of flux as precission weight? give more statistical importance to high-precission estimates? 


#3. GHG-method relevance
#Cores integrate fluxes during a longer time but with ex-situ conditions and only from small areas of wetlands (with biass towards bare sediment and small-depth open water). Given these limitations, CO2 fluxes might not be very representative of real in-situ conditions for the entire wetland. N2O and CH4 fluxes might be more representative when derived from cores (N2O: only viable data, CH4: longer integration is more trustworthy for ebullitive dynamics). HOW DO WE RECONCILE/COMBINE GHGfluxes from in-situ vs ex-situ?


rm(list=ls())

# ---- Directories ----
dropbox_root <- "C:/Users/Miguel/Dropbox/RESTORE4Cs - Fieldwork/Data" # You have to make sure this is pointing to the write folder on your local machine

#Path to CORES fluxes results:
cores_path<- paste0(dropbox_root,"/Cores/Cores_flux/")

#Path to save plots from exploration: 
plots_path <- paste0(dropbox_root,"/Cores/Stats and plots/")


# ---- Packages ----
library(tidyverse)
library(readxl)
library(ggplot2)
library(ggpubr)
library(ggResidpanel)
library(ggforce)

#Data transformation:
library(bestNormalize)

#ALL PACKAGES LOADED FOR in-situ status comparisons:  un-comment as needed
# #For data-handling & plotting
# library(lubridate)
# library(zoo)
# library(grid)
# library(egg)
# require(purrr)
# require(data.table)
# require(tools)
# library(hms)
# library(suncalc)
# #For modelling: 
# library(lmerTest) #similar to lme4 but allows to estimate significance of predictors (and interactions)
# library(bestNormalize)
# library(outliers) #Outlier exploration
# library(rstatix) #Outlier exploration
# library(emmeans)
# library(glmmTMB)
# library(DHARMa)
# library(performance)
# library(tibble)
# library(partR2)#for breaking down variance by fixed effect
# library(car)  # for leveneTest
# library(caret)  # for cross-validation
# library(MuMIn)  # for AICc if needed
# library(multcomp)# for emmeans post-hoc comparisons





#Repo-functions
repo_root <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
files.sources = list.files(path = paste0(repo_root,"/functions"), full.names = T)
for (f in files.sources){source(f)}



#1. Import & format-----

coreflux<- read.csv(paste0(cores_path,"All_core_GHGfluxes_and_CH4stocks.csv"))

head(coreflux)

#Add season, casepilot, subsite, status 
wideflux<- coreflux %>% 
  select(core_id,CO2_flux_mmol_per_m2_per_d,CH4_flux_micromol_per_m2_per_d,N2O_flux_micromol_per_m2_per_d) %>% 
  #ADD relevant qualitative variables for stats: 
  separate(core_id, into = c("season", "casepilot","status_num","corenum"),remove = F,sep = "-") %>% 
  mutate(subsite=paste(casepilot, status_num, sep = "-"),
         status=substr(status_num,1,1),
         #Add status_difalter (to be able to differentiate DA alterations)
         status_difalter=case_when(subsite%in%c("DA-A1","DA-A2")~status_num,
                                   TRUE~status)) %>% 
  #HOMOGENIZE units to gwp100 and calculate gwp_co2ch4 and gwp_all
  mutate(
    #Obtain mol per m2 per day: 
    co2_mol=CO2_flux_mmol_per_m2_per_d*1e-3,
    ch4_mol=CH4_flux_micromol_per_m2_per_d*1e-6,
    n2o_mol=N2O_flux_micromol_per_m2_per_d*1e-6,
    #Calculate Global warming potential (GWP100):  g of CO2eq per m2 per day
    gwp_co2=co2_mol*44.0095, #molar flux * molar mass* GWP100factor(1)
    gwp_ch4=ch4_mol*16.04246*27,#molar flux * molar mass* GWP100factor(27)
    gwp_n2o=n2o_mol*44.0128*273, #molar flux * molar mass * GWP100factor(273)
    gwp_co2ch4= gwp_co2+gwp_ch4, #Sum of gwp of co2 and ch4 only
    gwp_all=gwp_co2+gwp_ch4+gwp_n2o) %>%  #Sum of gwp of all GHGs (co2,ch4,n2o)
  #Organize variables:
  select(core_id, season, casepilot, subsite, status, status_difalter, gwp_co2, gwp_ch4, gwp_n2o, gwp_co2ch4, gwp_all)


#Pivot longer
longflux<- wideflux %>% 
  pivot_longer(cols = c(gwp_co2,gwp_ch4,gwp_n2o,gwp_co2ch4,gwp_all),names_to = "ghgspecies", values_to = "flux") %>% 
  mutate(ghgspecies=factor(ghgspecies, levels = c("gwp_co2", "gwp_ch4","gwp_n2o", "gwp_co2ch4","gwp_all")))



#2. Exploratory plots ----

##___status_per site-----
#DU.
longflux %>% 
  filter(casepilot=="DU") %>% 
  ggplot(aes(x=status, y=flux))+
  geom_boxplot()+
  facet_wrap(.~ghgspecies, scales="free")

#RI.
longflux %>% 
  filter(casepilot=="RI") %>% 
  ggplot(aes(x=status, y=flux))+
  geom_boxplot()+
  facet_wrap(.~ghgspecies, scales="free")

#CA.
longflux %>% 
  filter(casepilot=="CA") %>% 
  ggplot(aes(x=status, y=flux))+
  geom_boxplot()+
  facet_wrap(.~ghgspecies, scales="free")

#VA.
longflux %>% 
  filter(casepilot=="VA") %>% 
  ggplot(aes(x=status, y=flux))+
  geom_boxplot()+
  facet_wrap(.~ghgspecies, scales="free")

#DA.
longflux %>% 
  filter(casepilot=="DA") %>% 
  ggplot(aes(x=status_difalter, y=flux))+
  geom_boxplot()+
  facet_wrap(.~ghgspecies, scales="free")

#CU.
longflux %>% 
  filter(casepilot=="CU") %>% 
  ggplot(aes(x=status, y=flux))+
  geom_boxplot()+
  facet_wrap(.~ghgspecies, scales="free")



##___Preserved across sites-------
#CO2:
longflux %>% 
  filter(status=="P", ghgspecies%in%c("gwp_co2")) %>% 
  mutate(casepilot=factor(casepilot, levels = c("DU","RI","CA","VA","DA","CU"))) %>% 
  ggplot(aes(x=casepilot, y=flux))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(col ="#00BA38",position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, fill="#00BA38",size = 0.7)+
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
  guides(color = "none", fill="none")+
  labs(title=expression(CO[2]~GWP) ,y= expression(CO[2]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Casepilot"),
       fill=paste0("Status"))


#CH4:
longflux %>% 
  filter(status=="P", ghgspecies%in%c("gwp_ch4")) %>% 
  mutate(casepilot=factor(casepilot, levels = c("DU","RI","CA","VA","DA","CU"))) %>% 
  ggplot(aes(x=casepilot, y=flux))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(col ="#00BA38",position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, fill="#00BA38",size = 0.7)+
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
  guides(color = "none", fill="none")+
  labs(title=expression(CH[4]~GWP) ,y= expression(CH[4]~NEE~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Casepilot"),
       fill=paste0("Status"))


#N2O:
longflux %>% 
  filter(status=="P", ghgspecies%in%c("gwp_n2o")) %>% 
  mutate(casepilot=factor(casepilot, levels = c("DU","RI","CA","VA","DA","CU"))) %>% 
  ggplot(aes(x=casepilot, y=flux))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(col ="#00BA38",position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, fill="#00BA38",size = 0.7)+
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
  guides(color = "none", fill="none")+
  labs(title=expression(N2O~GWP) ,y= expression(N2O~NEE~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Casepilot"),
       fill=paste0("Status"))

#CO2+CH4:
longflux %>% 
  filter(status=="P", ghgspecies%in%c("gwp_co2ch4")) %>% 
  mutate(casepilot=factor(casepilot, levels = c("DU","RI","CA","VA","DA","CU"))) %>% 
  ggplot(aes(x=casepilot, y=flux))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(col ="#00BA38",position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, fill="#00BA38",size = 0.7)+
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
  guides(color = "none", fill="none")+
  labs(title=expression(Combined~CO2+CH4~GWP) ,y= expression(NEE~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Casepilot"),
       fill=paste0("Status"))

#CO2+CH4+N2O:
longflux %>% 
  filter(status=="P", ghgspecies%in%c("gwp_all")) %>% 
  mutate(casepilot=factor(casepilot, levels = c("DU","RI","CA","VA","DA","CU"))) %>% 
  ggplot(aes(x=casepilot, y=flux))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_violin(col ="#00BA38",position = position_dodge(0.9),linewidth = 1,scale = "width")+
  geom_sina(position = position_dodge(0.9),size = 0.5, scale="width")+
  geom_boxplot(width=0.2,outliers = F, fill="#00BA38",size = 0.7)+
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
  guides(color = "none", fill="none")+
  labs(title=expression(Combined~CO2+CH4+N2O~GWP) ,y= expression(NEE~(g~CO[2~eq]~d^-1~m^-2)),
       x=paste0("Casepilot"),
       fill=paste0("Status"))



#3. MODELS-----
#We will follow the same approach as for the in-situ data GLMM for all casepilots:

#flux~status*season+(1|subsite)


#3.0. Data trans options------

#IMPORTANT CONSIDERATION--------  
#DO WE REALLY NEED sign-preserving transformations for stats? 

#true that we have positive and negative fluxes but the transformed data is only used to compare between status, not to interpret, and, in any case, the interpretation will be based on the inverse-transformed estimates. 



#Options for data-transformation includes:
  
  #GOOD OPTIONS: Signed log, signed root, signed arcsin, 
  #POTENTIAL OPTION: Yeo-Johnson (maybe, sometimes swaps sign of very low-magnitude values,)
  #NON-APPROPRIATE: Box-Cox (non-negatives only). 

#TO IMPLEMENT: pseudolog. equivalent to logsing but with a tunning parameter for each dataset that optimizes normality and avoids effects of units (potential bimodality)

##Trans functions----
#Transformation functions (with inv_transformations for regaining og scale to interpret)
#Signed log transformation:
logsign <- function(x) {sign(x) * log(abs(x) + 1)}
inv_logsign <- function(x) {sign(x) * (exp(abs(x)) - 1)}
#Signed root transformation:
rootsign <- function(x) {sign(x) * sqrt(abs(x))}
inv_rootsign <- function(x) {sign(x) * (x^2)}
# #Signed arcsinh transformation:
arcsinhsign <- function(x) {sign(x) * asinh(abs(x))}
inv_arcsinhsign <- function(x) {sign(x) * sinh(abs(x))}



#Yeo-Johnson transformation (sometimes swaps signs of low-magnitude values)

# Yeo-Johnson transformation
yeojohnson_transform <- function(x) {
  # Estimate lambda using bestNormalize
  lambda <- yeojohnson(x)$lambda
  yeojohnson(x, lambda = lambda)
}

# Inverse Yeo-Johnson transformation
inv_yeojohnson <- function(y, lambda) {
  # Inverse based on the transformation formula
  ifelse(y >= 0,
         if (lambda == 0) exp(y) - 1 else (y * lambda + 1)^(1 / lambda) - 1,
         if (lambda == 2) -exp(-y) + 1 else -((-y * (2 - lambda) + 1)^(1 / (2 - lambda)) - 1))
}

# Wrapper to store lambda and transformed values
yeojohnson_wrapper <- function(x) {
  result <- yeojohnson(x)
  attr(result$x.t, "lambda") <- result$lambda
  result$x.t
}


#SELECT AND INSPECT BEST TRANSFORMATION: 

# Lista de transformaciones
transformations <- list(
  "Original" = identity,
  "Signed Log" = logsign,
  "Signed Sqrt" = rootsign,
  "Signed arcsin" = arcsinhsign,
  "Yeo-Johnson" = yeojohnson_wrapper
)

# Evaluate normality with Shapiro-Wilk
evaluate_transformations <- function(x) {
  results <- lapply(names(transformations), function(name) {
    x_t <- transformations[[name]](x)
    shapiro <- tryCatch(shapiro.test(x_t), error = function(e) NULL)
    data.frame(
      Transformation = name,
      Shapiro_W = if (!is.null(shapiro)) shapiro$statistic else NA,
      Shapiro_p = if (!is.null(shapiro)) shapiro$p.value else NA
    )
  })
  bind_rows(results)
}

# Visualize transformations: histograms + density
plot_histograms <- function(x) {
  df <- lapply(names(transformations), function(name) {
    data.frame(
      value = transformations[[name]](x),
      transformation = name
    )
  }) %>% bind_rows()
  
  ggplot(df, aes(x = value)) +
    geom_histogram(aes(y = ..density..), bins = 30, fill = "skyblue", color = "black") +
    geom_density(color = "red", linewidth = 1) +
    facet_wrap(~transformation, scales = "free") +
    theme_minimal() +
    labs(title = "Histograms")
}

# Visualize transformations: Q-Q plots
plot_qq <- function(x) {
  df <- lapply(names(transformations), function(name) {
    data.frame(
      value = transformations[[name]](x),
      transformation = name
    )
  }) %>% bind_rows()
  
  ggplot(df, aes(sample = value)) +
    stat_qq() +
    stat_qq_line() +
    facet_wrap(~transformation, scales = "free") +
    theme_minimal() +
    labs(title = "Q-Q Plots")
}



##Evaluate transformations----


#Check normality on CO2: datasets:
co2_table_trans<- tibble()
for (cp in unique(longflux$casepilot)) {
  for (ghg in c("gwp_ch4")) { #"gwp_ch4","gwp_n2o","gwp_co2andch4"
    x<- longflux %>% filter(casepilot==cp&ghgspecies==ghg) %>% filter(!is.na(flux)) %>% pull(flux)
    table_trans<- evaluate_transformations(x)
    table_trans$casepilot<- cp
    table_trans$ghgspecies<- ghg
    co2_table_trans<- bind_rows(co2_table_trans, table_trans)
    print(plot_histograms(x)+ggtitle(paste(cp, ghg)))
    print(plot_qq(x)+ggtitle(paste(cp, ghg)))
  }
}

#get most-normal transformation for each case
co2_best_trans<- co2_table_trans %>% 
  group_by(ghgspecies,casepilot) %>% 
  filter(Shapiro_p==max(Shapiro_p)) %>% 
  arrange(ghgspecies,casepilot)
co2_best_trans


#Function bestNormalize, it evaluates multiple transformations, providing the most-normal one (similar to our manual approach, but does not include sign-preserving transformations).
bestNormalize(longflux %>% filter(casepilot=="CA"&ghgspecies=="gwp_ch4") %>% pull(flux),allow_orderNorm=F )
  
longflux

a<-longflux %>% filter(casepilot=="VA"&ghgspecies=="gwp_ch4") %>% pull(flux)

a_norm<-bestNormalize(a, allow_orderNorm = F)

#bestNormalize can store the normalization so that afterwards we can use it to back-transform to original scale (unclear how it works)



#Pseudolog tests---------

#Pseudolog transformation can be implemented via the scales package: it contains two different functions:
    #pseudo_log_trans() it returns a trans object, in which we can "store" the parameters and use directly within ggplot: For use with ggplot2 axes (trans =)
    #transform_pseudo_log() it returns a vector of transformed values directly

#We could include the pseudolog into the BestNormalize approach. By creating a custom function:




#Notes on pseudolog (~signed log) transformation: https://stratosida.github.io/regression-regrets/Pseudo_log_explainer.html
#finding-a-parameter-that-best-achieves-normality
#Our type of CO2 distribution, is the result of substracting two log-normal distributions (emission-uptake), which leads to a mostly simetrical distribution with heavy tails and mostly centered around 0. 
pseudo_log <-
  function(x, sigma = 1, base = 10)
    asinh(x / (2 * sigma)) / log(base)

inv_pseudo_log <-
  function(x, sigma = 1, base = 10)
    2 * sigma * sinh(x * log(base))
#THis requires to optimice sigma for each dataset, which can make the distribution slimmer, fatter or even introduce bimodality. Only for modelling, not to use for visual representation.
#CA-CO2:
#example using CA.  DOES NOT WORK APPROPRIATELY, likely due to negatives (it should, adjust code.)
x<- sort(data4models %>% filter(casepilot=="DA", ghgspecies=="gwp_co2") %>% filter(!is.na(dailyflux)) %>% pull(dailyflux)*10)
#Creates pool of sigmas to chose best
sigmas <- 2 ** seq(-10, 10, 0.2)

#Calculates correlation between normal distribution and my data untransformed.
origcor <- cor(qnorm((1:length(x) - 0.5) / length(x)), x)

#Recalculates correlation with normal distribution for all options of sigma.
ncorx <-
  sapply(sigmas, function(X)
    cor(qnorm((1:length(
      x
    ) - 0.5) / length(x)), pseudo_log(x, X)))

cat("Optimal sigma: ")
(optsigma <- sigmas[ncorx == max(ncorx)])

plot(log2(sigmas), ncorx, ylab = "Correlation with normal deviates", ylim =
       c(0, 1))
points(log2(sigmas)[ncorx == max(ncorx)], max(ncorx), pch = 19)
abline(h = origcor)
legend(
  "bottomright",
  lty = c(1, NA, NA),
  pch = c(NA, 1, 19),
  legend = c("Original", "Pseudo-log", "Pseudo-log with optimal sigma")
)
box()
hist(pseudo_log(x, optsigma))




#Integration of pseudolog to BestNormalize:-------
rm(list=ls())
# 1. Constructor: fits pseudo-log transformation with optimal sigma
best_pseudo_log <- function(x, base = 10, sigma_range = seq(0.01, 5, 0.1), standardize = TRUE) {
  tryCatch({
    # Evaluate normality across candidate sigma values
    stats <- sapply(sigma_range, function(sigma) {
      x_t <- tryCatch(
        pseudo_log_trans(base = base, sigma = sigma)$transform(x),
        error = function(e) rep(NA, length(x))
      )
      if (length(unique(x_t[!is.na(x_t)])) < 3) return(NA)
      tryCatch(stats::shapiro.test(x_t[!is.na(x_t)])$statistic, error = function(e) NA)
    })
    
    # Select best sigma based on Shapiro-Wilk statistic
    if (all(is.na(stats))) stop("No valid transformation found.")
    best_sigma <- sigma_range[which.max(stats)]
    
    # Apply transformation with best sigma
    trans <- pseudo_log_trans(base = base, sigma = best_sigma)
    x_t <- trans$transform(x)
    
    # Optionally standardize
    mu <- mean(x_t, na.rm = TRUE)
    sd <- sd(x_t, na.rm = TRUE)
    if (standardize) x_t <- (x_t - mu) / sd
    
    # Compute normality statistic (Pearson P/df)
    ptest <- nortest::pearson.test(x_t[!is.na(x_t)])
    norm_stat <- unname(ptest$statistic / ptest$df)
    
    # Return S3 object
    val <- list(
      x = x,
      x_t = x_t,
      base = base,
      sigma = best_sigma,
      trans = trans,
      mean = mu,
      sd = sd,
      standardize = standardize,
      n = length(x) - sum(is.na(x)),
      norm_stat = norm_stat
    )
    class(val) <- c("best_pseudo_log", class(val))
    return(val)
  }, error = function(e) {
    structure(list(norm_stat = NA), class = "best_pseudo_log")
  })
}

# 2. Predict method: applies transformation or inverse using fitted parameters


#Predict method with NA handling:
predict.best_pseudo_log <- function(object, newdata = NULL, inverse = FALSE, ...) {
  if (is.null(newdata) && !inverse) newdata <- object$x
  if (is.null(newdata) && inverse) newdata <- object$x_t
  
  # Preserve NA positions
  na_idx <- is.na(newdata)
  
  if (inverse) {
    newdata[!na_idx] <- object$trans$inverse(newdata[!na_idx])
    if (object$standardize) {
      newdata[!na_idx] <- newdata[!na_idx] * object$sd + object$mean
    }
  } else {
    newdata[!na_idx] <- object$trans$transform(newdata[!na_idx])
    if (object$standardize) {
      newdata[!na_idx] <- (newdata[!na_idx] - object$mean) / object$sd
    }
  }
  
  return(unname(newdata))
}



# 3. Print method (optional): displays transformation details
print.best_pseudo_log <- function(x, ...) {
  cat(ifelse(x$standardize, "Standardized", "Non-Standardized"),
      "pseudo_log(x) Transformation with", x$n, "nonmissing obs.:\n",
      "Relevant statistics:\n",
      "- base =", x$base, "\n",
      "- sigma =", x$sigma, "\n",
      "- mean (before standardization) =", x$mean, "\n",
      "- sd (before standardization) =", x$sd, "\n")
}

# 4. Register custom transformation structure to use within bestNormalize
pseudolog_transform <- list(
  best_pseudo_log = best_pseudo_log,
  predict.best_pseudo_log = predict.best_pseudo_log,
  print.best_pseudo_log = print.best_pseudo_log
)


#5. Example usage and real-world workflow

# Simulate data with NA values
set.seed(123)
x <- c(abs(rt(100, df = 2)) * sample(c(-1, 1), 100, replace = TRUE), NA_real_, NA_real_)

# Apply transformation
obj <- best_pseudo_log(x)

# Print transformation details
print(obj)

# Transform the data
x_transformed <- predict(obj, x)
cat("\nTransformed values:\n")
print(head(x_transformed, 10))

# Backtransform the data
x_backtransformed <- predict(obj, x_transformed, inverse = TRUE)
cat("\nBacktransformed values:\n")
print(head(x_backtransformed, 10))

# Check that original and backtransformed values match (excluding NAs) & NAs are preserved
comparison <- data.frame(
  original = x,
  backtransformed = x_backtransformed,
  match = (round(x, 6) == round(x_backtransformed, 6)) | (is.na(x) & is.na(x_backtransformed))
)

cat("\nComparison of original and backtransformed values:\n")
print(head(comparison, 10))



# Apply bestNormalize with custom transformation
bn_result <- bestNormalize(x, new_transforms = pseudolog_transform, standardize = FALSE, warn = TRUE)
print(bn_result)

# --- Use Case 1: Apply transformation to training data ---
x_trans <- predict(bn_result$chosen_transform)
# Back-transform to original scale
x_back <- predict(bn_result$chosen_transform, inverse = TRUE)

# --- Use Case 2: Apply transformation to new data ---
x_new <- x + rnorm(length(x), sd = 0.05)  # simulate similar new data
x_new_trans <- predict(bn_result$chosen_transform, newdata = x_new)
x_new_back <- predict(bn_result$chosen_transform, newdata = x_new_trans, inverse = TRUE)

# --- Use Case 3: Validate consistency ---
max_diff <- max(abs(x - x_back), na.rm = TRUE)
max_diff_new <- max(abs(x_new - x_new_back), na.rm = TRUE)

cat("Max diff (original vs back-transformed):", max_diff, "\n")
cat("Max diff (new vs back-transformed):", max_diff_new, "\n")
