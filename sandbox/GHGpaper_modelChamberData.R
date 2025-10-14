#GHGpaper_modelchamberdata.R

#Author: Miguel Cabrera-Brufau
#Date: September 2025
#Project: Restore4cs


#Description----
#This scrip is used to model the effect of restoration in each casepilot. Using net daily GHG exchange rates of appropriate incubations, data preparation is in GHGpaper_prepChamberData.R script. 
#We will use only CO2 and CH4 data, the N2O fluxes from chambers are very few (only S4) and not very reliable (many non-significant fluxes).

#DECISSIONS: 
  #Modelling approach: GLMMtmb with general formula dailyflux~status*season + (1|subsite), subsite as random to account for repeated samplings and subsite-specific patterns

  #Contrast: Set contrast options to "contr. sum", so that we are not using any one level of a factor as the reference level, but rather the tests will asses whether there is an overall average effect of one factor across all levels of the other factors. For example, is there an effect of status across all levels of season?


#STEPS: 
#1. Import and format data: filter-out N2O, remove Nas, as.factor levels, 
#2. Transform data: Using BestNormalize package, find and apply for each casepilot, the transformation that maximizes normality. 
#3. Optimize and select model for each casepilot*GHGspecies combo (co2, ch4, GWPco2andch4).
  # 3.1. Decide final structure: family distribution, heterocedasticity
  # 3.2. Detect and treat outliers (if present).

#4. Export model results:
  # 4.1. Summary_chambermodels.csv (normalization, structure, distributionfamily, dispformula, significances, Homoced_R2c, Homoced_R2m, PredObs_cor, nobs, extra(AIC?, BIC?, check common practices for reporting).   
  # 4.2. Residualtests_chambermodels.csv: pvalues of residual tests for models (model assumptions)
  # 4.3. Emmeans_chambermodels.csv: estimated marginal means for each casepilot*GHGspecies and post-hoc groups. 




rm(list = ls()) # clear workspace
cat("/014") # clear console


# ---- Packages ----
#For data-handling & plotting
library(tidyverse)
library(readxl)
library(lubridate)
library(zoo)
# library(ggplot2)
library(ggpubr)
library(rlang)
# library(grid)
# library(egg)
library(ggtext)#for plot captions
require(purrr)
require(data.table)
require(tools)
library(hms)
library(suncalc)
# library(cowplot) #Alligned axis
# library(ggforce)


#For modelling: 
library(vegan) #to test for strata-compositional differences
library(bestNormalize) #For data-transformations
library(scales)#For pseudolog transformation
library(glmmTMB) #For glmm modelling
library(DHARMa) #For residual evaluation
library(performance)#Checking singluar models
library(rstatix) #Outlier exploration & Anova of models
library(emmeans) #To estimate EMMs
library(multcomp)# for emmeans post-hoc comparisons
library(multcompView) #for clc (letters groups)


#Repo-functions
repo_root <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
files.sources = list.files(path = paste0(repo_root,"/functions"), full.names = T)
for (f in files.sources){source(f)}



# ---- Directories ----
dropbox_root <- "C:/Users/Miguel/Dropbox/RESTORE4Cs - Fieldwork/Data" # You have to make sure this is pointing to the write folder on your local machine

#Path to datasets for and to save model outputs: 
paper_path<- paste0(dropbox_root,"/GHG/Main_paper_results/")

extraplots_path<- paste0(paper_path,"Exploratory_figures/")

#0.Contrast options-----
#NOTES on contrasts: in R by default, contrast is "contr.treatment" wich uses the first level of each factor as the "reference" level and all subsequent as "treatment" levels. This implies that, with default contrast type, when looking at the main effects of my model, what is shown for each main effect is whether there is a significant effect of 1 factor (eg. status) at the reference level of all other factors (i.e. season). What we want is to asses whether there is an overall average effect of status across all levels of season. For this purpose we need to set contrasts to "contr. sum", and always call Anova (model, type="III").

#Set contrasts to contr.sum (for the whole R session)
options(contrasts = c("contr.sum", "contr.poly"))


#0. Import and format--------

data4models<- read.csv(paste0(paper_path,"ChamberData4paper.csv"))

#Format and filter:
data4models<- data4models %>% 
  #Remove N2O-data
  filter(ghgspecies%in%c("co2","ch4","gwp_co2andch4")) %>% 
  #Rename to GWPco2andch4
  mutate(ghgspecies=if_else(ghgspecies=="gwp_co2andch4","GWPco2andch4",ghgspecies)) %>% 
  #Remove NAs
  filter(!is.na(dailyflux)) %>% 
  mutate(season=factor(season, ordered = F),
         status=factor(status, ordered = F),
         strata=factor(strata, ordered = F)) %>% 
  mutate(vegpresence=if_else(strata=="vegetated","Vegetated","Non-vegetated")) %>% 
  mutate(vegpresence=factor(vegpresence, ordered=F))
  



#MODELS------


#IMPORTANT! -----
#Evaluate model structure and assumptions (residuals), not significance of effects nor even explained variability (we expect very little variance explained by status in many cases). 

#1st step is to decide the approach based on our data structure and distribution. 

#Fixed decissions: 
#Each case-pilot is modelled independently (different data-distributions and transformation needs)
#Need to include subsite as random effect (explicitly account for 1.repeated samplings)

#Approach: use Generalized Linear Mixed Models (glmmTMB)
#They can account for non-normal data (even after transformation) by using different distribution families
#They can incorporate heterocedasticity across predictors 



#checks for data (in all combinations for models)
#1. Normality of response variable (not essential but will suggest need for transformation)
#2. Normality of residuals (with the various potential transformations). Moderately important, severe non-normality can bias p-values and confidence intervals. 
#3. Homocedasticity of residuals across factors. Important for valid standard errors and p-values. Can be omited by including dispformula 
#4. Linearity. Residuals vs fitted should not have patterns. Important 
#5. Independence of observations, accounted for by random effect subsite


#STEPS in model optimization: 
#1. Decide transformation: maximize normality of data
#2. Decide modelling family and structure (Gaussian vs T-family, homoscedastic vs dispformula)
#3. Detect and evaluate impact of outliers--> Decide outlier exclusion (for model)


#Note on Outliers -
#All data has been screened for wrong incubations, however, by chance, we might encounter more extreme outliers than predicted based on data-distribution. 


#for GLMMtmb, the most appropriate outlier detection test is that of the DHARMA package, the test first identifies if the number of “outliers” that you get is significantly more than one would expect by chance, then identifies the suspects, if any. Common procedure is to re-fit the model without the suspects (if the test was significant), and check whether the AIC improves substantially without these: 
#example: 

# outlier_logical <- testOutliers(sim_res) #outlier detection, produces plot and returns logical vector with flagged outliers
# model_clean <- update(model, data = mydata[!outlier_logical, ]) #recalculate model without outliers
# AIC(model, model_clean)# Compare AIC
# summary(model)$coefficients # Compare fixed effects of original model
# summary(model_clean)$coefficients #Compare fixed effects of model without outliers

#Criteria for outlier removal: ΔAIC ≥ 2 suggests meaningful model improvement, ΔAIC ≥ 4–7 suggests strong evidence

# OUTLIERS should be considered for each model, so they cannot be removed in advance but after model structure has been optimized. 


#0. Custom Functions -------
#Here functions to access relevant results from models. Will be used to summarise information.  


#Formula to get marginal and conditional R2s from homocedastic model list:
#We  can only partition the explanatory power between fixed and random effects (marginal vs conditional R2) for models that assume homocedasticity (i.e. no dispformula).
get_homoc_R2s <- function(model_list) {
  library(glmmTMB)
  library(performance)
  library(dplyr)
  library(purrr)
  
  map_df(names(model_list), function(cpghg) {
    model <- model_list[[cpghg]]
    transformation<- table_trans %>% filter(dataset==cpghg) %>% pull(trans_name)
    
    r2_vals<- NULL
    # Try to extract R2 values
    r2_vals <- tryCatch({
      MuMIn::r.squaredGLMM(model) #Using MuMin to be able to get marginal R2 even with singularities
    }, error = function(e) {
      warning(paste("R2 extraction failed for", cp))
      return(c(NA_real_, NA_real_))
      
    })
    
    tibble(
      dataset = cpghg,
      transformation=transformation,
      formula = deparse(formula(model)),
      dispformula = deparse(model$modelInfo$allForm$dispformula),
      family = model$modelInfo$family$family,
      homoced_R2m = r2_vals[1,1],
      homoced_R2c = r2_vals[1,2]
    )
  })
}




#Function to get various pseudoR2 estimates:
#We Cannot calculate R2 (conditional & marginal) for heterocedastic models (those with dispformula), (or gamma family, unclear for t_family distributions), instead we calculate various pseudo-R2 metrics (log-Likelyhood based) to get the "improvement in model fit relative to a null model (intercept + random effects), capturing the combined explanatory contribution of fixed effects".
#We have 4 estimates: 
#Efron pseudoR2: squared correlation between observed and predicted values. 
#The other 3 pseudoR2 come from Log-likelihood comparisons between actual and null model. Nagelkerke pseudoR2 is scaled to 0-1 values, we will use this as it is the easiest to interpret.

get_pseudoR2s<- function(model_list){
  
  map_df(names(model_list), function(cpghg){
    model<- model_list[[cpghg]]
    transformation<- table_trans %>% filter(dataset==cpghg) %>% pull(trans_name)
    
    
    obs<- model$frame[,1]
    pred<- predict(model, type="response")
    #Efron_pseudoR2 is equal to the squared correlation between predicted and observed values
    Efron_pseudoR2<- 1 - sum((obs-pred)^2)/sum((obs-mean(obs))^2)
    
    #rcompanion calcualtes 3 different pseudoR2 based on comparision of actual vs null model (mantaining the structure, but removing the predictors from the model formula)
    res<-rcompanion::nagelkerke(model)  
    
    tibble(
      dataset = cpghg,
      transformation = transformation,
      formula = deparse(formula(model)),
      dispformula = deparse(model$modelInfo$allForm$dispformula),
      family = model$modelInfo$family$family,
      pseudoR2_Efron = Efron_pseudoR2,
      pseudoR2_McFadden = res$Pseudo.R.squared.for.model.vs.null[1],
      pseudoR2_Cox_and_snell=res$Pseudo.R.squared.for.model.vs.null[2],
      pseudoR2_Nagelkerke = res$Pseudo.R.squared.for.model.vs.null[3]
    )
  })
}


#Function to extract model structure and all fixed effects significance (regardless of their naming, can differ in order or presence across the models in model_list)
test_get_all_anova_results <- function(model_list) {
  library(car)
  library(glmmTMB)
  library(dplyr)
  library(purrr)
  
  map(names(model_list), function(cpghg) {
    model <- model_list[[cpghg]]
    # transformation <- table_trans %>% filter(dataset == cpghg) %>% pull(trans_name)
    
    # Run Type III ANOVA
    anova_res <- tryCatch(car::Anova(model, type = "III"), error = function(e) NULL)
    
    # Extract p-values and rename terms
    pval_tibble <- if (!is.null(anova_res)) {
      terms <- rownames(anova_res)
      pvals <- anova_res[, "Pr(>Chisq)"]
      names(pvals) <- paste0(gsub(":", ".", terms), "_pval")
      as_tibble(as.list(pvals))
    } else {
      tibble()
    }
    
    # Combine with metadata
    tibble(
      dataset = cpghg,
      # transformation = transformation,
      formula = deparse(formula(model)),
      family = model$modelInfo$family$family,
      dispformula = deparse(model$modelInfo$allForm$dispformula),
      nobs = dim(model$frame)[1]
    ) %>%
      bind_cols(pval_tibble)
  }) %>%
    bind_rows()
}



#Function to summarise ANOVA results: get main effects significance for list of best_models
get_anova_results <- function(model_list) {
  library(car)
  library(glmmTMB)
  library(dplyr)
  library(purrr)
  
  map_df(names(model_list), function(cpghg) {
    model <- model_list[[cpghg]]
    transformation<- table_trans %>% filter(dataset==cpghg) %>% pull(trans_name)
    # Run Type III ANOVA
    anova_res <- tryCatch({
      car::Anova(model, type = "III")
    }, error = function(e) {
      warning(paste("ANOVA failed for", cp))
      return(NULL)
    })
    
    # Extract p-values safely
    get_pval <- function(term) {
      if (is.null(anova_res)) return(NA_real_)
      if (term %in% rownames(anova_res)) return(anova_res[term, "Pr(>Chisq)"])
      return(NA_real_)
    }
    
    tibble(
      dataset = cpghg,
      transformation=transformation, 
      formula = deparse(formula(model)),
      family = model$modelInfo$family$family,
      dispformula = deparse(model$modelInfo$allForm$dispformula),
      nobs = dim(model$frame)[1],
      intercept_pval = get_pval("(Intercept)"),
      status_pval = get_pval("status"),
      season_pval = get_pval("season"),
      interaction_pval = get_pval("status:season")
    )
  })
}


#Function to summarise ANOVA results: get main effects significance for list of best_models
#Effects for significance are named, if they do not exist within the model, the p-value returned is NA. 
get_anova_results_3way <- function(model_list) {
  library(car)
  library(glmmTMB)
  library(dplyr)
  library(purrr)
  
  map_df(names(model_list), function(cpghg) {
    model <- model_list[[cpghg]]
    transformation<- table_trans %>% filter(dataset==cpghg) %>% pull(trans_name)
    # Run Type III ANOVA
    anova_res <- tryCatch({
      car::Anova(model, type = "III")
    }, error = function(e) {
      warning(paste("ANOVA failed for", cp))
      return(NULL)
    })
    
    # Extract p-values safely
    get_pval <- function(term) {
      if (is.null(anova_res)) return(NA_real_)
      if (term %in% rownames(anova_res)) return(anova_res[term, "Pr(>Chisq)"])
      return(NA_real_)
    }
    
    tibble(
      dataset = cpghg,
      transformation=transformation, 
      formula = deparse(formula(model)),
      family = model$modelInfo$family$family,
      dispformula = deparse(model$modelInfo$allForm$dispformula),
      nobs = dim(model$frame)[1],
      intercept_pval = get_pval("(Intercept)"),
      status_pval = get_pval("status"),
      season_pval = get_pval("season"),
      status.season_pval = get_pval("status:season"),
      status.strata_pval =get_pval("status:strata"),
      status.season.strata_pval = get_pval("status:season:strata")
    )
  })
}




#Summary of Model assumptions: for best_model_list, Requires named lists: best_model_list.
#We want to structure the results in a table: model structure,  distribution family, dispformula, 
#add pvalue for: normality of residuals, heteroscedasticity (levene test).

summarize_dharma_diagnostics <- function(model_list) {
  library(DHARMa)
  library(dplyr)
  library(purrr)
  
  map_df(names(model_list), function(cpghg) {
    model <- model_list[[cpghg]]
    data  <- model$frame
    transformation<- table_trans %>% filter(dataset==cpghg) %>% pull(trans_name)
    # Simulate residuals
    sim_res <- tryCatch({
      simulateResiduals(fittedModel = model, plot = FALSE)
    }, error = function(e) {
      warning(paste("DHARMa simulation failed for", cp))
      return(NULL)
    })
    
    if (is.null(sim_res)) {
      return(tibble(
        dataset = cpghg,
        transformation=transformation, 
        uniformity_pval = NA_real_,
        dispersion_pval = NA_real_,
        hetero_status_pval = NA_real_,
        hetero_season_pval = NA_real_
      ))
    }
    
    # Run tests
    tibble(
      dataset = cpghg,
      transformation=transformation, 
      formula = deparse(formula(model)),
      dispformula = deparse(model$modelInfo$allForm$dispformula),
      family = model$modelInfo$family$family,
      convergence=check_convergence(model),
      singular=check_singularity(model,),
      uniformity_pval = testUniformity(sim_res)$p.value,
      dispersion_pval = testDispersion(sim_res)$p.value,
      hetero_status_pval = testCategorical(sim_res,catPred =data$status,  plot = F)$homogeneity$`Pr(>F)`[1],
      hetero_season_pval = testCategorical(sim_res,catPred =data$season,  plot = F)$homogeneity$`Pr(>F)`[1]
    )
  })
}


#Function to Produce comparative histograms and qqplots of untransformed and best-Normalize transformed data. plot_bestNormalize <- function(bn_obj, n_bins = 30, title = NULL) {
plot_bestNormalize <- function(bn_obj, n_bins = 30, title = NULL) {
  library(ggplot2)
  library(gridExtra)
  
  # Extract original and transformed data
  x_orig <- bn_obj$x
  x_trans <- predict(bn_obj)
  
  # Get transformation name from class attribute
  trans_name <- class(bn_obj$chosen_transform)[1]
  
  # Make combined title (if title is provided)
  plot_title <- paste0("Transformation: ", trans_name, if (!is.null(title)) paste0(" - ", title) else "")
  
  # Perform KS tests
  ks_orig <- ks.test(x_orig, "pnorm", mean(x_orig), sd(x_orig))
  ks_trans <- ks.test(x_trans, "pnorm", mean(x_trans), sd(x_trans))
  
  # Create data frame for plotting
  df <- data.frame(
    value = c(x_orig, x_trans),
    type = factor(rep(c("Original", "Transformed"), each = length(x_orig)), levels = c("Original", "Transformed"))
  )
  
  # Histogram
  p_hist <- ggplot(df, aes(x = value, fill = type)) +
    geom_histogram(bins = n_bins, alpha = 0.6, position = "identity") +
    facet_wrap(~type, scales = "free") +
    labs(title = paste("Histogram:", plot_title), x = "Value", y = "Count") +
    theme_minimal() +
    scale_fill_manual(values = c("Original" = "steelblue", "Transformed" = "darkgreen")) +
    theme(legend.position = "none")
  
  # QQ plot data
  qq_data <- data.frame(
    sample = c(x_orig, x_trans),
    type = factor(rep(c("Original", "Transformed"), each = length(x_orig)), levels = c("Original", "Transformed"))
  )
  
  # KS p-values for annotation
  ks_labels <- data.frame(
    type = c("Original", "Transformed"),
    label = c(
      paste0("KS p-value: ", signif(ks_orig$p.value, 4)),
      paste0("KS p-value: ", signif(ks_trans$p.value, 4))
    ),
    x = c(Inf, Inf),
    y = c(-Inf, -Inf)
  )
  
  # QQ plot with annotation
  p_qq <- ggplot(qq_data, aes(sample = sample)) +
    stat_qq(color = "black") +
    stat_qq_line(color = "red", linetype = "dashed") +
    facet_wrap(~type, scales = "free") +
    geom_text(data = ks_labels, aes(x = x, y = y, label = label), hjust = 1.1, vjust = -0.5, inherit.aes = FALSE) +
    labs(title = paste("QQ Plot:", plot_title), x = "Theoretical Quantiles", y = "Sample Quantiles") +
    theme_minimal()
  
  # Combine plots
  gridExtra::grid.arrange(p_hist, p_qq, nrow = 2)
}


#Function that plots observed vs predicted values for a model (or pair of models fitted under the same data), coloring points via color_var
plot_obs_vs_pred_models <- function(model1, data, model2 = NULL, color_var = "vegpresence") {
  # Get name of the data object
  data_name <- deparse(substitute(data))
  
  # Extract and clean fixed effect formula (RHS only)
  model1_name <- gsub("^dailyflux_trans\\s*~\\s*", "", deparse(formula(model1)))
  
  df1 <- data.frame(
    obs = model1$frame$dailyflux_trans,
    pred = predict(model1, type = "response"),
    color = data[[color_var]],
    model = model1_name
  )
  
  if (!is.null(model2)) {
    model2_name <- gsub("^dailyflux_trans\\s*~\\s*", "", deparse(formula(model2)))
    
    df2 <- data.frame(
      obs = model2$frame$dailyflux_trans,
      pred = predict(model2, type = "response"),
      color = data[[color_var]],
      model = model2_name
    )
    
    df_combined <- bind_rows(df1, df2) %>%
      mutate(model = factor(model, levels = c(model1_name, model2_name)))
  } else {
    df_combined <- df1 %>%
      mutate(model = factor(model, levels = model1_name))
  }
  
  # Compute R² and label positions
  r2_labels <- df_combined %>%
    group_by(model) %>%
    summarise(
      r2 = summary(lm(obs ~ pred))$r.squared,
      label = paste0("R² = ", round(r2, 2)),
      .groups = "drop"
    )
  
  # Get range of observed values
  obs_range <- range(df_combined$obs, na.rm = TRUE)
  
  # Plot
  ggplot(df_combined, aes(x = pred, y = obs)) +
    geom_point(aes(col = color), size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    geom_smooth(method = "lm", col = "black", se = FALSE) +
    geom_smooth(method = "lm", aes(col = color), se = FALSE) +
    facet_wrap(~model) +
    geom_text(
      data = r2_labels,
      aes(x = min(obs_range), y = max(obs_range), label = label),
      inherit.aes = FALSE,
      hjust = 0, vjust = 1,
      size = 4
    ) +
    scale_x_continuous(limits = obs_range) +
    scale_y_continuous(limits = obs_range) +
    labs(
      title = paste(data_name, "Observed vs Predicted"),
      subtitle = paste("Best-fit lines and R² from linear regression in model scale\nColored by", color_var),
      x = "Predicted",
      y = "Observed",
      color = color_var
    ) +
    
    theme_bw() +
    theme(
      plot.caption = element_textbox_simple(),
      legend.position = "bottom",
      legend.box.margin = margin(t = -5, b = -5),  # reduce top/bottom margin
      legend.spacing = unit(0.2, "lines")          # reduce spacing between legend items
    )
}

# #Example usage:
# plot_obs_vs_pred_models(model1=ca_simplemodel_co2,
#                         data=ca_co2, 
#                         model2=ca_best_co2,
#                         color_var = "status")


##Pseudolog-------

#Notes on pseudolog (~signed log) transformation: https://stratosida.github.io/regression-regrets/Pseudo_log_explainer.html#finding-a-parameter-that-best-achieves-normality
#Pseudolog tranformation needs tunning of its parameter to each data-set (similar to yeo-johnson). To incorporate this transformation into the BestNormalize, we need to create several functions. 

#Pseudolog transformation allows and maintains NAs 
#Range of potential sigmas are data-driven: based on 10th-90th percentile ranges. 

# 1. Constructor: fits pseudo-log transformation with optimal sigma
#best_pseudo_log with data-driven optimization of sigma using 10th and 90th percentile range
best_pseudo_log <- function(x, base = 10, sigma_range = NULL, standardize = TRUE) {
  tryCatch({
    # Clean input
    x_clean <- x[!is.na(x) & is.finite(x)]
    if (length(x_clean) < 3) stop("Not enough valid data points.")
    
    # Compute data-driven sigma_range if not provided (using 10th-90th percentile range)
    if (is.null(sigma_range)) {
      q <- quantile(abs(x_clean), probs = c(0.1, 0.9), na.rm = TRUE)
      lower <- max(q[1] / 10, .Machine$double.eps)  # Avoid zero or near-zero
      upper <- max(q[2], lower * 10)
      sigma_range <- 10^seq(log10(lower), log10(upper), length.out = 100)#100 potential sigmas
    }
    
    # Evaluate normality across candidate sigma values
    stats <- sapply(sigma_range, function(sigma) {
      x.t <- tryCatch(
        pseudo_log_trans(base = base, sigma = sigma)$transform(x_clean),
        error = function(e) rep(NA, length(x_clean))
      )
      if (length(unique(x.t[!is.na(x.t)])) < 3) return(NA)
      tryCatch(stats::shapiro.test(x.t[!is.na(x.t)])$statistic, error = function(e) NA)
    })
    
    # Select best sigma
    if (all(is.na(stats))) stop("No valid transformation found.")
    best_sigma <- sigma_range[which.max(stats)]
    
    # Apply transformation
    trans <- pseudo_log_trans(base = base, sigma = best_sigma)
    x.t <- trans$transform(x)
    
    # Optionally standardize
    mu <- mean(x.t, na.rm = TRUE)
    sd <- sd(x.t, na.rm = TRUE)
    if (standardize) x.t <- (x.t - mu) / sd
    
    # Pearson normality statistic
    ptest <- nortest::pearson.test(x.t[!is.na(x.t)])
    norm_stat <- unname(ptest$statistic / ptest$df)
    
    # Return result
    val <- list(
      x = x,
      x.t = x.t,
      base = base,
      sigma = best_sigma,
      trans = trans,
      mean = mu,
      sd = sd,
      standardize = standardize,
      n = length(x_clean),
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
  if (is.null(newdata) && inverse) newdata <- object$x.t
  
  na_idx <- is.na(newdata)
  
  if (inverse) {
    if (object$standardize) {
      newdata[!na_idx] <- newdata[!na_idx] * object$sd + object$mean
    }
    newdata[!na_idx] <- object$trans$inverse(newdata[!na_idx])
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



#0.Transformations------
#The data-transformations do not necessarily have to mantain the flux sign. Even though it has essential ecological meaning (emission vs capture), the models do not need to consider it, we will back-transform the output of the models for plotting and interpretation purposes. 

#We will use a completely objective decision for transformation: using BestNormalize function (with pseudolog tranformation as a potential option), we will chose the transformation that results in the most-normal data distribution. 

# Apply bestNormalize with custom transformation pseudolog option

#USE bestNormalize (allowing for pseudo-log trans) to obtain for each dataset (ghg*casepilot combo) the most-Normal transformation: 
rm(table_trans, table_trans_i)
for (cp in unique(data4models$casepilot)) {
  for (ghg in  unique(data4models$ghgspecies)) {
    
    x<- data4models %>% filter(casepilot==cp&ghgspecies==ghg) %>%
      pull(dailyflux)
    bn_result<- bestNormalize(x = x, new_transforms = pseudolog_transform, standardize = TRUE, warn = TRUE, k = 5,allow_orderNorm = F)
    #calculate normality test (shapiro)
    x_trans <- bn_result$chosen_transform$x.t
    pval <- shapiro.test(x_trans)$p.value
    
    table_trans_i <- tibble(
      dataset = paste(cp, ghg,sep = "_"),
      casepilot = cp,
      ghgspecies = ghg,
      nobs= length(x),
      trans_name = attr(x = bn_result$chosen_transform,which = "class")[1],
      trans_Pearson_p=bn_result$chosen_transform$norm_stat, #Pearson's P (the lower, the more-normal data, only used for ranking our transformation options)
      shapiro_pval=pval, #P value of Normality
      bn_result = list(bn_result)  # wrap BestNormalize object in list to store as list-column
    )
    
    
    #Create table_trans for first round of loop, then apend to table
    if(cp==unique(data4models$casepilot)[1]&
       ghg==unique(data4models$ghgspecies)[1]){
      table_trans<- table_trans_i
    }else{ table_trans<- bind_rows(table_trans, table_trans_i)}
  }
}

#View all bestNormalizing transformations:
table_trans

rm(cp, bn_result, table_trans_i, ghg, x, x_trans, pval)



# table_trans[14,]
# 
# #See effect of transformation: 
# table_trans %>%
#   filter(casepilot=="RI"&
#            ghgspecies=="ch4") %>%
#   pull(bn_result) %>%
#   pluck(1) %>% 
#   plot_bestNormalize()


#Automatically apply the best transformation to each dataset within data4models using the bestNormalize objects stored in table_trans 

# Create a named list of BestNormalize objects for easy lookup
bn_list <- table_trans %>%
  dplyr::select(dataset, bn_result) %>%
  deframe()

# Add a key column to data4models for matching
data4models <- data4models %>%
  mutate(dataset = paste(casepilot, ghgspecies, sep = "_"))

# Apply the corresponding BestNormalize transformation
data4models <- data4models %>%
  rowwise() %>%
  mutate(
    trans_name = table_trans$trans_name[match(dataset, paste(table_trans$casepilot, table_trans$ghgspecies, sep = "_"))],
    dailyflux_trans = if (!is.null(bn_list[[dataset]])) predict(bn_list[[dataset]], dailyflux) else NA_real_
  ) %>%
  ungroup()


#__________------

#1. IS strata random?-----

#The type of strata of a particular flux measurement (vegetated, open water, bare) is very likely to have an influence over the flux magnitude. Thus, including this information will likely improve the explanatory capacity of the two fixed effects (status and season). HOWEVER: if strata composition is also dependent on status, we would be adding a random effect that is highly correlated to the main effect that we want to test and this would considerably hinder our capacity to identify the effect of status (variance could be "stolen" by strata, reducing the expanatory power of status).

#In order to decide whether we can incorporate strata as a random effect, we have to first check that the strata composition is actually random (i.e. not related to the principal effect that we want to test: status). 

#APROACH: for each casepilot, we will test whether the strata composition of the sites (as derived from chamber distribution during our seasonal  samplings) is related to their status. We will use PERMANOVA for each casepilot, via the vegan::adonis2 function, using bray-curtis dissilarity as compositional distance, status as the predictor, and including the site to control for repeated measures across seasons. 

#Insead of stratacomp (not directly sampled but calculated and "tuned" to representativeness, we have eliminated seasonal variability where it was due to  arbitrary sampling decission) we need to use actual strata composition of the data is going to be modeled. I.e. re-calculate strata composition of data4models.


#OLD: using stratacomp previously calculated & manipulated to account for sampling inconsistencies 
# stratacomp<- read.csv(paste0(paper_path, "Stratacomposition_in-situ.csv"))
# head(stratacomp)
# 
# #Prepare stratacomp:
# stratacomp_formated<- stratacomp %>% 
#   dplyr::select(sampling, strata, strata_prop) %>% 
#   separate(sampling, into = c("season","casepilot","subsite")) %>% 
#   mutate(status=case_when(grepl("A", subsite)~"Altered",
#                           grepl("P", subsite)~"Preserved",
#                           grepl("R", subsite)~"Restored")) %>% 
#   mutate(strata=gsub("open ","", strata)) %>% 
#   pivot_wider(names_from = strata, values_from = strata_prop) %>% 
#   dplyr::select(casepilot, season, subsite,status, bare, vegetated, water) %>% 
#   mutate(status=as.factor(status), subsite=as.factor(subsite), casepilot=as.factor(casepilot))


#NEW: using actual composition of data4models (incorrect incubations removed, not adjusted to constant composition where no seasonality was expected). vegpresence used instead of strata 

vegpresenceprop_formated<- data4models %>% 
  dplyr::select(plotcode, casepilot, season, status, vegpresence) %>% 
  #remove duplicates (same plotcode, different ghgspecies)
  distinct() %>% 
  group_by(casepilot, season, status, vegpresence) %>%
  #Calculate vegpresence deployment counts and pivot_wider
  summarise(n_vegpresence=sum(!is.na(vegpresence))) %>% 
  pivot_wider(names_from = vegpresence, values_from = n_vegpresence,values_fill = 0) %>% 
  #Calculate vegpresence proportion
  mutate(sum_all=`Non-vegetated`+Vegetated) %>% 
  mutate(Vegetated=Vegetated/sum_all,
         `Non-vegetated`=`Non-vegetated`/sum_all) %>% 
  dplyr::select(-sum_all) %>% 
  #Override CU composition: constant actual composition, differences in chamber deployment are due to systematic sampling biass in Restored sites
  mutate(Vegetated=if_else(casepilot=="CU", 0.2,Vegetated),
         `Non-vegetated`=if_else(casepilot=="CU",0.8, `Non-vegetated`))

head(vegpresenceprop_formated)

#Create function to perform test for each casepilot: 
# Function to run adonis2 for a subset of data
run_adonis_for_casepilot <- function(casepilot_id, df) {
  subset_df <- df %>% filter(casepilot == casepilot_id)
  
  # Extract composition matrix
  composition_vars <- c("Vegetated", "Non-vegetated")
  composition_df <- subset_df[, composition_vars]
  
  if (length(unique(subset_df$status)) < 2) {
    return(data.frame(
      casepilot = casepilot_id,
      R2 = NA,
      F = NA,
      p_value = NA,
      Note = "Not enough Status levels"
    ))
  }
  
  # Check for constant composition
  if (all(apply(composition_df, 2, function(x) length(unique(x)) == 1))) {
    return(data.frame(
      casepilot = casepilot_id,
      R2 = NA,
      F = NA,
      p_value = NA,
      Note = "Constant composition - test skipped"
    ))
  }
  
  # Run adonis2
  adonis_result <- adonis2(
    composition_df ~ status + season,
    data = subset_df,
    method = "bray",
    permutations = 9999
  )
  
  # Extract results
  result_row <- adonis_result[1, ]  # First row is the model term (status)
  
  data.frame(
    casepilot = casepilot_id,
    R2 = result_row$R2,
    F = result_row$F,
    p_value = result_row$`Pr(>F)`,
    Note = NA
  )
}

# Apply to each unique Casepilot
casepilot_ids <- unique(vegpresenceprop_formated$casepilot)
adonis_results <- lapply(casepilot_ids, run_adonis_for_casepilot, df = vegpresenceprop_formated)

# Combine into a single data frame
adonis_summary <- do.call(rbind, adonis_results) 
adonis_summary

#Status is only significantly linked to the vegpresence composition for:
  #RI, clear and strong link,


#GENERAL MODEL potential OPTIONS: 

#For RI and DA: cannot include strata in the general model (linked to status). USE ONLY BASIC FORMULA: dailyflux_trans~season*status + (1|subsite)

#For CA, CU, DU, VA: allow for strata as random in formula, decide final model based on residuals and fit. random effects: (1|subsite) or (1|subsite/strata)


#Always use transformed data (most parsimonious option).
#IF residuals are not normal, try using t_family instead of gaussian. 
#Do not include dispformula (do not include heterocedasticity across predictors)
#Force R2 calculation for singular models (especially relevant if strata is added as nested random effect). 


#2. NEW General models-------

#CO2______-------
##2.1. CA_co2 (ok) -------

#Is vegpresence random?
adonis_summary %>% filter(casepilot=="CA")
#YES, vegpresence is not significantly linked to status, we can allow vegpresnce as a fixed effect

#Subset data and check transformation: 
ca_co2<- data4models %>% filter(casepilot=="CA"&ghgspecies=="co2")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="CA"&
           ghgspecies=="co2") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = ca_co2,
                        family = gaussian(),
                        dispformula = ~1)
#Evaluate basic model:
check_convergence(m1_gaus_nostrata) #OK
check_singularity(m1_gaus_nostrata,tolerance = 1e-8) #OK
r2(m1_gaus_nostrata, tolerance = 1e-10) 
Anova(m1_gaus_nostrata)
summary(m1_gaus_nostrata)
res<- simulateResiduals(m1_gaus_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_nostrata)
em<- emmeans(m1_gaus_nostrata, ~ status,weights = "equal") #using equal weights as we are only averaging across season and all should have the same impact.
summary(em)
pairs(em)
#DECISSION: Good residuals, keep



#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                   data = ca_co2,
                   family = gaussian(),
                   dispformula = ~1)
#Evaluate model:
check_convergence(m2_gaus_vegpresence) #OK
check_singularity(m2_gaus_vegpresence,tolerance = 1e-8) #Singular
r2(m2_gaus_vegpresence, tolerance = 1e-10)
Anova(m2_gaus_vegpresence)
summary(m2_gaus_vegpresence)
res<- simulateResiduals(m2_gaus_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m2_gaus_vegpresence)
em<- emmeans(m2_gaus_vegpresence, ~ status, weights = "proportional")
summary(em)
pairs(em)
#DECISSION: BAD residuals



#T_family:
#Most basic: t_family, only status and season 
m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = ca_co2,
                        family = t_family(),
                        dispformula = ~1)
#Evaluate basic model:
check_convergence(m3_t_nostrata)
check_singularity(m3_t_nostrata,tolerance = 1e-8) 
r2(m3_t_nostrata, tolerance=1e-10) 
Anova(m3_t_nostrata)
summary(m3_t_nostrata)
res<- simulateResiduals(m3_t_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m3_t_nostrata)
em<- emmeans(m3_t_nostrata, ~ status)
summary(em)
pairs(em)
#DECISSION: does not converge. 


#Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                      data = ca_co2,
                      family = t_family(),
                      dispformula = ~1)
#Evaluate basic model:
check_convergence(m4_t_vegpresence) #OK
check_singularity(m4_t_vegpresence,tolerance = 1e-8)
r2(m4_t_vegpresence, tolerance = 1e-9)
Anova(m4_t_vegpresence)
summary(m4_t_vegpresence)
res<- simulateResiduals(m4_t_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m4_t_vegpresence)
check_overdispersion(m4_t_vegpresence)
em<- emmeans(m4_t_vegpresence, ~ status, weights = "proportional")
summary(em)
pairs(em)
#GOOD residuals, keep


#Compare models 
anova(m1_gaus_nostrata,m4_t_vegpresence)
#Adding vegpresence is significantly better (better fit worth complexity increase)

#Save best models (simple and complex best)
ca_simplemodel_co2<- m1_gaus_nostrata
ca_complexmodel_co2<- m4_t_vegpresence


#PLot observed vs predicted of two options: 
plot_obs_vs_pred_models(model1 = ca_simplemodel_co2, 
                        model2 = ca_complexmodel_co2,
                        data=ca_co2,
                        color_var = "vegpresence")

ggsave(filename = "CA_co2_Observed_VS_predicted.png", 
       width = 160,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)




##2.2. CU_co2 (ok) -------

#Is vegpresence random?
adonis_summary %>% filter(casepilot=="CU")
#YES, by definition, vegpresence is not linked to status, we can allow vegpresnce as a fixed effect

#Subset data and check transformation: 
cu_co2<- data4models %>% filter(casepilot=="CU"&ghgspecies=="co2")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="CU"&
           ghgspecies=="co2") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = cu_co2,
                           family = gaussian(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m1_gaus_nostrata) #OK
check_singularity(m1_gaus_nostrata,tolerance = 1e-8)
r2(m1_gaus_nostrata, tolerance = 1e-10) 
Anova(m1_gaus_nostrata)
summary(m1_gaus_nostrata)
res<- simulateResiduals(m1_gaus_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_nostrata)
em<- emmeans(m1_gaus_nostrata, ~ status,weights = "equal") #using equal weights as we are only averaging across season and all should have the same impact.
summary(em)
pairs(em)
#DECISSION: Bad residuals


#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                              data = cu_co2,
                              family = gaussian(),
                              dispformula = ~1)
#Evaluate model:
check_convergence(m2_gaus_vegpresence) #OK
check_singularity(m2_gaus_vegpresence,tolerance = 1e-8) #Singular
r2(m2_gaus_vegpresence, tolerance = 1e-10)
Anova(m2_gaus_vegpresence)
summary(m2_gaus_vegpresence)
res<- simulateResiduals(m2_gaus_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m2_gaus_vegpresence)
em<- emmeans(m2_gaus_vegpresence, ~ status, weights = "proportional")
summary(em)
pairs(em)
#DECISSION: BAD residuals



#T_family:
#Most basic: t_family, only status and season 
m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = cu_co2,
                        family = t_family(),
                        dispformula = ~1)
#Evaluate basic model:
check_convergence(m3_t_nostrata)
check_singularity(m3_t_nostrata,tolerance = 1e-8) 
r2(m3_t_nostrata, tolerance=1e-10) 
Anova(m3_t_nostrata)
summary(m3_t_nostrata)
res<- simulateResiduals(m3_t_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m3_t_nostrata)
em<- emmeans(m3_t_nostrata, ~ status)
summary(em)
pairs(em)
#DECISSION: Better residuals, although they still fail


#Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                           data = cu_co2,
                           family = t_family(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m4_t_vegpresence) #OK
check_singularity(m4_t_vegpresence,tolerance = 1e-8)
r2(m4_t_vegpresence, tolerance = 1e-9)
Anova(m4_t_vegpresence)
summary(m4_t_vegpresence)
res<- simulateResiduals(m4_t_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m4_t_vegpresence)
em<- emmeans(m4_t_vegpresence, ~ status, weights = "proportional")
summary(em)
pairs(em)
#GOOD residuals, keep


#Compare models 
anova(m3_t_nostrata,m4_t_vegpresence)
#Adding vegpresence is significantly better (better fit worth complexity increase)

#Save best models (simple and complex best)
cu_simplemodel_co2<- m3_t_nostrata
cu_complexmodel_co2<- m4_t_vegpresence


#PLot observed vs predicted of two options: 
plot_obs_vs_pred_models(model1 = cu_simplemodel_co2, 
                        model2 = cu_complexmodel_co2,
                        data=cu_co2,
                        color_var = "vegpresence")

ggsave(filename = "CU_co2_Observed_VS_predicted.png", 
       width = 160,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)

##2.3. DA_co2 (ok-ish) -------

#Is vegpresence random?
adonis_summary %>% filter(casepilot=="DA")
#YES, vegpresence is not linked to status, we can allow vegpresnce as a fixed effect

#Subset data and check transformation: 
da_co2<- data4models %>% filter(casepilot=="DA"&ghgspecies=="co2")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="DA"&
           ghgspecies=="co2") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = da_co2,
                           family = gaussian(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m1_gaus_nostrata) #OK
check_singularity(m1_gaus_nostrata,tolerance = 1e-8)
r2(m1_gaus_nostrata, tolerance = 1e-10) 
Anova(m1_gaus_nostrata)
summary(m1_gaus_nostrata)
res<- simulateResiduals(m1_gaus_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_nostrata)
em<- emmeans(m1_gaus_nostrata, ~ status,weights = "equal") #using equal weights as we are only averaging across season and all should have the same impact.
summary(em)
pairs(em)
#DECISSION: Bad residuals


#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                              data = da_co2,
                              family = gaussian(),
                              dispformula = ~1)
#Evaluate model:
check_convergence(m2_gaus_vegpresence) #OK
check_singularity(m2_gaus_vegpresence,tolerance = 1e-8) #Singular
r2(m2_gaus_vegpresence, tolerance = 1e-10)
Anova(m2_gaus_vegpresence)
summary(m2_gaus_vegpresence)
res<- simulateResiduals(m2_gaus_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m2_gaus_vegpresence)
em<- emmeans(m2_gaus_vegpresence, ~ status, weights = "proportional")
summary(em)
pairs(em)
#DECISSION: BAD residuals



#T_family:
#Most basic: t_family, only status and season 
m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = da_co2,
                        family = t_family(),
                        dispformula = ~1)
#Evaluate basic model:
check_convergence(m3_t_nostrata)
check_singularity(m3_t_nostrata,tolerance = 1e-8) 
r2(m3_t_nostrata, tolerance=1e-10) 
Anova(m3_t_nostrata)
summary(m3_t_nostrata)
res<- simulateResiduals(m3_t_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m3_t_nostrata)
em<- emmeans(m3_t_nostrata, ~ status, weights = "equal")
summary(em)
pairs(em)
#DECISSION: Better residuals, although they still fail


#Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                           data = da_co2,
                           family = t_family(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m4_t_vegpresence) #OK
check_singularity(m4_t_vegpresence,tolerance = 1e-8)
r2(m4_t_vegpresence, tolerance = 1e-10)
Anova(m4_t_vegpresence)
summary(m4_t_vegpresence)
res<- simulateResiduals(m4_t_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m4_t_vegpresence)
em<- emmeans(m4_t_vegpresence, ~ status, weights = "proportional")
summary(em)
pairs(em)
#Better residuals, although they present issues

#Compare models 
anova(m3_t_nostrata,m4_t_vegpresence)
#Adding vegpresence is significantly better (better fit worth complexity increase)

#Save best models (simple and complex best)
da_simplemodel_co2<- m3_t_nostrata
da_complexmodel_co2<- m4_t_vegpresence


#PLot observed vs predicted of two options: 
plot_obs_vs_pred_models(model1 = da_simplemodel_co2, 
                        model2 = da_complexmodel_co2,
                        data=da_co2,
                        color_var = "vegpresence")

ggsave(filename = "DA_co2_Observed_VS_predicted.png", 
       width = 160,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)


##2.1. DU_co2 (ok) -------
#Is vegpresence random?
adonis_summary %>% filter(casepilot=="DU")
#YES, vegpresence is not linked to status, we can allow vegpresnce as a fixed effect

#Subset data and check transformation: 
du_co2<- data4models %>% filter(casepilot=="DU"&ghgspecies=="co2")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="DU"&
           ghgspecies=="co2") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = du_co2,
                           family = gaussian(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m1_gaus_nostrata) #OK
check_singularity(m1_gaus_nostrata,tolerance = 1e-8)
r2(m1_gaus_nostrata, tolerance = 1e-10) 
Anova(m1_gaus_nostrata)
summary(m1_gaus_nostrata)
res<- simulateResiduals(m1_gaus_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_nostrata)
em<- emmeans(m1_gaus_nostrata, ~ status,weights = "equal") #using equal weights as we are only averaging across season and all should have the same impact.
summary(em)
pairs(em)
#DECISSION: GOOD residuals


#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                              data = du_co2,
                              family = gaussian(),
                              dispformula = ~1)
#Evaluate model:
check_convergence(m2_gaus_vegpresence) #OK
check_singularity(m2_gaus_vegpresence,tolerance = 1e-8) #Singular
r2(m2_gaus_vegpresence, tolerance = 1e-10)
Anova(m2_gaus_vegpresence)
summary(m2_gaus_vegpresence)
res<- simulateResiduals(m2_gaus_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m2_gaus_vegpresence)
em<- emmeans(m2_gaus_vegpresence, ~ status, weights = "proportional")
summary(em)
pairs(em)
#DECISSION: BAD residuals



#T_family:
#Most basic: t_family, only status and season 
m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = du_co2,
                        family = t_family(),
                        dispformula = ~1)
#Evaluate basic model:
check_convergence(m3_t_nostrata)
check_singularity(m3_t_nostrata,tolerance = 1e-8) 
r2(m3_t_nostrata, tolerance=1e-10) 
Anova(m3_t_nostrata)
summary(m3_t_nostrata)
res<- simulateResiduals(m3_t_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m3_t_nostrata)
em<- emmeans(m3_t_nostrata, ~ status, weights = "equal")
summary(em)
pairs(em)
#DECISSION: Model does not converge


#Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                           data = du_co2,
                           family = t_family(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m4_t_vegpresence) #OK
check_singularity(m4_t_vegpresence,tolerance = 1e-8)
r2(m4_t_vegpresence, tolerance = 1e-10)
Anova(m4_t_vegpresence)
summary(m4_t_vegpresence)
res<- simulateResiduals(m4_t_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m4_t_vegpresence)
em<- emmeans(m4_t_vegpresence, ~ status, weights = "proportional")
summary(em)
pairs(em)
#DECISSION: good residuals

#Compare models 

anova(m1_gaus_nostrata,m4_t_vegpresence)
#Adding vegpresence is significantly better (better fit worth complexity increase)

#Save best models (simple and complex best)
du_simplemodel_co2<- m1_gaus_nostrata
du_complexmodel_co2<- m4_t_vegpresence


#PLot observed vs predicted of two options: 
plot_obs_vs_pred_models(model1 = du_simplemodel_co2, 
                        model2 = du_complexmodel_co2,
                        data=du_co2,
                        color_var = "vegpresence")

ggsave(filename = "DU_co2_Observed_VS_predicted.png", 
       width = 160,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)

##2.1. RI_co2 (ok) -------
#Is vegpresence random?
adonis_summary %>% filter(casepilot=="RI")
#NO! vegpresence is SIGNIFICANTLY linked to status, we cannot allow vegpresnce as a fixed effect

#Subset data and check transformation: 
ri_co2<- data4models %>% filter(casepilot=="RI"&ghgspecies=="co2")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="RI"&
           ghgspecies=="co2") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = ri_co2,
                           family = gaussian(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m1_gaus_nostrata) #OK
check_singularity(m1_gaus_nostrata,tolerance = 1e-8)
r2(m1_gaus_nostrata, tolerance = 1e-10) 
Anova(m1_gaus_nostrata)
summary(m1_gaus_nostrata)
res<- simulateResiduals(m1_gaus_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_nostrata)
em<- emmeans(m1_gaus_nostrata, ~ status,weights = "equal") #using equal weights as we are only averaging across season and all should have the same impact.
summary(em)
pairs(em)
#DECISSION: GOOD residuals


#SKIP! cannot include vegpresence as fixed with status
if(F){
#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                              data = ri_co2,
                              family = gaussian(),
                              dispformula = ~1)
#Evaluate model:
check_convergence(m2_gaus_vegpresence) #OK
check_singularity(m2_gaus_vegpresence,tolerance = 1e-8) #
r2(m2_gaus_vegpresence, tolerance = 1e-10)
Anova(m2_gaus_vegpresence)
summary(m2_gaus_vegpresence)
res<- simulateResiduals(m2_gaus_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m2_gaus_vegpresence)
em<- emmeans(m2_gaus_vegpresence, ~ status, weights = "proportional")
summary(em)
pairs(em)
#DECISSION: 
}


#SKIP t_options: 
if (F){
#T_family:
#Most basic: t_family, only status and season 
m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = ri_co2,
                        family = t_family(),
                        dispformula = ~1)
#Evaluate basic model:
check_convergence(m3_t_nostrata)
check_singularity(m3_t_nostrata,tolerance = 1e-8) 
r2(m3_t_nostrata, tolerance=1e-10) 
Anova(m3_t_nostrata)
summary(m3_t_nostrata)
res<- simulateResiduals(m3_t_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m3_t_nostrata)
em<- emmeans(m3_t_nostrata, ~ status, weights = "equal")
summary(em)
pairs(em)
#DECISSION: 


#Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                           data = ri_co2,
                           family = t_family(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m4_t_vegpresence)
check_singularity(m4_t_vegpresence,tolerance = 1e-8)
r2(m4_t_vegpresence, tolerance = 1e-10)
Anova(m4_t_vegpresence)
summary(m4_t_vegpresence)
res<- simulateResiduals(m4_t_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m4_t_vegpresence)
em<- emmeans(m4_t_vegpresence, ~ status, weights = "proportional")
summary(em)
pairs(em)
#DECISSION: 
}


#Compare models 
anova(m1_gaus_nostrata,m3_t_nostrata)

#CANNOT HAVE VEGPRESENCE AS EFFECT FOR RIA DE AVEIRO, most basic model is best (and only allowed)

#Save best models (simple and complex best)
ri_simplemodel_co2<- m1_gaus_nostrata
# ri_complexmodel_co2<- NA  NO complex model possible


#PLot observed vs predicted ONLY option: 
plot_obs_vs_pred_models(model1 = ri_simplemodel_co2, 
                        # model2 = ri_complexmodel_co2,
                        data=ri_co2,
                        color_var = "vegpresence")

ggsave(filename = "RI_co2_Observed_VS_predicted.png", 
       width = 160,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)

##2.1. VA_co2 (ok) -------

#Is vegpresence random?
adonis_summary %>% filter(casepilot=="VA")
#YES, vegpresence is not linked to status, we can allow vegpresnce as a fixed effect

#Subset data and check transformation: 
va_co2<- data4models %>% filter(casepilot=="VA"&ghgspecies=="co2")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="VA"&
           ghgspecies=="co2") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = va_co2,
                           family = gaussian(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m1_gaus_nostrata) #OK
check_singularity(m1_gaus_nostrata,tolerance = 1e-8)
r2(m1_gaus_nostrata, tolerance = 1e-10) 
Anova(m1_gaus_nostrata)
summary(m1_gaus_nostrata)
res<- simulateResiduals(m1_gaus_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_nostrata)
em<- emmeans(m1_gaus_nostrata, ~ status,weights = "equal") #using equal weights as we are only averaging across season and all should have the same impact.
summary(em)
pairs(em)
#DECISSION: BAD residuals (bad model)


#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                              data = va_co2,
                              family = gaussian(),
                              dispformula = ~1)
#Evaluate model:
check_convergence(m2_gaus_vegpresence) #OK
check_singularity(m2_gaus_vegpresence,tolerance = 1e-8) #Singular
r2(m2_gaus_vegpresence, tolerance = 1e-10)
Anova(m2_gaus_vegpresence)
summary(m2_gaus_vegpresence)
res<- simulateResiduals(m2_gaus_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m2_gaus_vegpresence)
em<- emmeans(m2_gaus_vegpresence, ~ status, weights = "proportional")
summary(em)
pairs(em)
#DECISSION: BAD residuals



#T_family:
#Most basic: t_family, only status and season 
m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = va_co2,
                        family = t_family(),
                        dispformula = ~1)
#Evaluate basic model:
check_convergence(m3_t_nostrata)
check_singularity(m3_t_nostrata,tolerance = 1e-8) 
r2(m3_t_nostrata, tolerance=1e-10) 
Anova(m3_t_nostrata)
summary(m3_t_nostrata)
res<- simulateResiduals(m3_t_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m3_t_nostrata)
em<- emmeans(m3_t_nostrata, ~ status, weights = "equal")
summary(em)
pairs(em)
#DECISSION: Model does not converge


#Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                           data = va_co2,
                           family = t_family(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m4_t_vegpresence) #OK
check_singularity(m4_t_vegpresence,tolerance = 1e-8)
r2(m4_t_vegpresence, tolerance = 1e-10)
Anova(m4_t_vegpresence)
summary(m4_t_vegpresence)
res<- simulateResiduals(m4_t_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m4_t_vegpresence)
em<- emmeans(m4_t_vegpresence, ~ status, weights = "proportional")
summary(em)
pairs(em)
#DECISSION: good residuals

#Compare models 

anova(m1_gaus_nostrata,m4_t_vegpresence)
#Adding vegpresence is significantly better (better fit worth complexity increase)

#Save best models (simple and complex best)
va_simplemodel_co2<- m1_gaus_nostrata
va_complexmodel_co2<- m4_t_vegpresence


#PLot observed vs predicted of two options: 
plot_obs_vs_pred_models(model1 = va_simplemodel_co2, 
                        model2 = va_complexmodel_co2,
                        data=va_co2,
                        color_var = "vegpresence")

ggsave(filename = "VA_co2_Observed_VS_predicted.png", 
       width = 160,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)


#CH4______-------

##2.1. CA_ch4 (ok) -------

#Is vegpresence random?
adonis_summary %>% filter(casepilot=="CA")
#YES, vegpresence is not significantly linked to status, we can allow vegpresnce as a fixed effect

#Subset data and check transformation: 
ca_ch4<- data4models %>% filter(casepilot=="CA"&ghgspecies=="ch4")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="CA"&
           ghgspecies=="ch4") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = ca_ch4,
                           family = gaussian(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m1_gaus_nostrata) #OK
check_singularity(m1_gaus_nostrata,tolerance = 1e-8) #OK
r2(m1_gaus_nostrata, tolerance = 1e-10) 
Anova(m1_gaus_nostrata)
res<- simulateResiduals(m1_gaus_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_nostrata)
em<- emmeans(m1_gaus_nostrata, ~ status,weights = "equal") #using equal weights as we are only averaging across season and all should have the same impact.
summary(em)
pairs(em)
#DECISSION: Good residuals, keep



#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                              data = ca_ch4,
                              family = gaussian(),
                              dispformula = ~1)
#Evaluate model:
check_convergence(m2_gaus_vegpresence) #OK
check_singularity(m2_gaus_vegpresence,tolerance = 1e-8) #Singular
r2(m2_gaus_vegpresence, tolerance = 1e-10)
Anova(m2_gaus_vegpresence)
res<- simulateResiduals(m2_gaus_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m2_gaus_vegpresence)
em<- emmeans(m2_gaus_vegpresence, ~ status, weights = "proportional")
summary(em)
pairs(em)
#DECISSION: GOOD residuals, keep


#SKIP, gaussians already good
if(F){
#T_family:
#Most basic: t_family, only status and season 
m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = ca_ch4,
                        family = t_family(),
                        dispformula = ~1)
#Evaluate basic model:
check_convergence(m3_t_nostrata)
check_singularity(m3_t_nostrata,tolerance = 1e-8) 
r2(m3_t_nostrata, tolerance=1e-10) 
Anova(m3_t_nostrata)
res<- simulateResiduals(m3_t_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m3_t_nostrata)
em<- emmeans(m3_t_nostrata, ~ status)
summary(em)
pairs(em)
#DECISSION: 


#Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                           data = ca_ch4,
                           family = t_family(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m4_t_vegpresence) #OK
check_singularity(m4_t_vegpresence,tolerance = 1e-8)
r2(m4_t_vegpresence, tolerance = 1e-9)
Anova(m4_t_vegpresence)
res<- simulateResiduals(m4_t_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m4_t_vegpresence)
check_overdispersion(m4_t_vegpresence)
em<- emmeans(m4_t_vegpresence, ~ status, weights = "proportional")
summary(em)
pairs(em)
#DECISSION: 
}

#Compare models 
anova(m1_gaus_nostrata,m2_gaus_vegpresence)
#Adding vegpresence does NOT significantly improve results (too complex for the increased fit)

#Save best models (simple and complex best)
ca_simplemodel_ch4<- m1_gaus_nostrata
ca_complexmodel_ch4<- m2_gaus_vegpresence


#PLot observed vs predicted of two options: 
plot_obs_vs_pred_models(model1 = ca_simplemodel_ch4, 
                        model2 = ca_complexmodel_ch4,
                        data=ca_ch4,
                        color_var = "vegpresence")

ggsave(filename = "CA_ch4_Observed_VS_predicted.png", 
       width = 160,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)




##2.2. CU_ch4 (ok) -------

#Is vegpresence random?
adonis_summary %>% filter(casepilot=="CU")
#YES, by definition, vegpresence is not linked to status, we can allow vegpresnce as a fixed effect

#Subset data and check transformation: 
cu_ch4<- data4models %>% filter(casepilot=="CU"&ghgspecies=="ch4")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="CU"&
           ghgspecies=="ch4") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = cu_ch4,
                           family = gaussian(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m1_gaus_nostrata) #OK
check_singularity(m1_gaus_nostrata,tolerance = 1e-8)
r2(m1_gaus_nostrata, tolerance = 1e-10) 
Anova(m1_gaus_nostrata)
summary(m1_gaus_nostrata)
res<- simulateResiduals(m1_gaus_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_nostrata)
em<- emmeans(m1_gaus_nostrata, ~ status,weights = "equal") #using equal weights as we are only averaging across season and all should have the same impact.
summary(em)
pairs(em)
#DECISSION: Good-enough residuals


#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                              data = cu_ch4,
                              family = gaussian(),
                              dispformula = ~1)
#Evaluate model:
check_convergence(m2_gaus_vegpresence) #OK
check_singularity(m2_gaus_vegpresence,tolerance = 1e-8) #Singular
r2(m2_gaus_vegpresence, tolerance = 1e-10)
Anova(m2_gaus_vegpresence)
summary(m2_gaus_vegpresence)
res<- simulateResiduals(m2_gaus_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m2_gaus_vegpresence)
em<- emmeans(m2_gaus_vegpresence, ~ status, weights = "proportional")
summary(em)
pairs(em)
#DECISSION: Good residuals


#SKIP, GAUSSIAN ALREADY GOOD
if(F){
#T_family:
#Most basic: t_family, only status and season 
m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = cu_ch4,
                        family = t_family(),
                        dispformula = ~1)
#Evaluate basic model:
check_convergence(m3_t_nostrata)
check_singularity(m3_t_nostrata,tolerance = 1e-8) 
r2(m3_t_nostrata, tolerance=1e-10) 
Anova(m3_t_nostrata)
summary(m3_t_nostrata)
res<- simulateResiduals(m3_t_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m3_t_nostrata)
em<- emmeans(m3_t_nostrata, ~ status)
summary(em)
pairs(em)
#DECISSION:


#Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                           data = cu_ch4,
                           family = t_family(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m4_t_vegpresence) #OK
check_singularity(m4_t_vegpresence,tolerance = 1e-8)
r2(m4_t_vegpresence, tolerance = 1e-9)
Anova(m4_t_vegpresence)
summary(m4_t_vegpresence)
res<- simulateResiduals(m4_t_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m4_t_vegpresence)
em<- emmeans(m4_t_vegpresence, ~ status, weights = "proportional")
summary(em)
pairs(em)
#DECISSION: 
}

#Compare models 
anova(m1_gaus_nostrata,m2_gaus_vegpresence)
#Adding vegpresence is significantly better (better fit worth complexity increase)

#Save best models (simple and complex best)
cu_simplemodel_ch4<- m1_gaus_nostrata
cu_complexmodel_ch4<- m2_gaus_vegpresence


#PLot observed vs predicted of two options: 
plot_obs_vs_pred_models(model1 = cu_simplemodel_ch4, 
                        model2 = cu_complexmodel_ch4,
                        data=cu_ch4,
                        color_var = "vegpresence")

ggsave(filename = "CU_ch4_Observed_VS_predicted.png", 
       width = 160,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)

##2.3. DA_ch4 (ok-ish) -------

#Is vegpresence random?
adonis_summary %>% filter(casepilot=="DA")
#YES, vegpresence is not linked to status, we can allow vegpresnce as a fixed effect

#Subset data and check transformation: 
da_ch4<- data4models %>% filter(casepilot=="DA"&ghgspecies=="ch4")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="DA"&
           ghgspecies=="ch4") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = da_ch4,
                           family = gaussian(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m1_gaus_nostrata) #OK
check_singularity(m1_gaus_nostrata,tolerance = 1e-8)
r2(m1_gaus_nostrata, tolerance = 1e-10) 
Anova(m1_gaus_nostrata)
summary(m1_gaus_nostrata)
res<- simulateResiduals(m1_gaus_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_nostrata)
em<- emmeans(m1_gaus_nostrata, ~ status,weights = "equal") #using equal weights as we are only averaging across season and all should have the same impact.
summary(em)
pairs(em)
#DECISSION: Bad residuals


#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                              data = da_ch4,
                              family = gaussian(),
                              dispformula = ~1)
#Evaluate model:
check_convergence(m2_gaus_vegpresence) #OK
check_singularity(m2_gaus_vegpresence,tolerance = 1e-8) #Singular
r2(m2_gaus_vegpresence, tolerance = 1e-10)
Anova(m2_gaus_vegpresence)
summary(m2_gaus_vegpresence)
res<- simulateResiduals(m2_gaus_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m2_gaus_vegpresence)
em<- emmeans(m2_gaus_vegpresence, ~ status, weights = "proportional")
summary(em)
pairs(em)
#DECISSION: BAD residuals



#T_family:
#Most basic: t_family, only status and season 
m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = da_ch4,
                        family = t_family(),
                        dispformula = ~1)
#Evaluate basic model:
check_convergence(m3_t_nostrata)
check_singularity(m3_t_nostrata,tolerance = 1e-8) 
r2(m3_t_nostrata, tolerance=1e-10) 
Anova(m3_t_nostrata)
summary(m3_t_nostrata)
res<- simulateResiduals(m3_t_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m3_t_nostrata)
em<- emmeans(m3_t_nostrata, ~ status, weights = "equal")
summary(em)
pairs(em)
#DECISSION: Better residuals, although they still fail


#Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                           data = da_ch4,
                           family = t_family(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m4_t_vegpresence) #OK
check_singularity(m4_t_vegpresence,tolerance = 1e-8)
r2(m4_t_vegpresence, tolerance = 1e-10)
Anova(m4_t_vegpresence)
summary(m4_t_vegpresence)
res<- simulateResiduals(m4_t_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m4_t_vegpresence)
em<- emmeans(m4_t_vegpresence, ~ status, weights = "proportional")
summary(em)
pairs(em)
#Better residuals, although they present issues

#Compare models 
anova(m3_t_nostrata,m4_t_vegpresence)
#Adding vegpresence is significantly better (better fit worth complexity increase)

#Save best models (simple and complex best)
da_simplemodel_ch4<- m3_t_nostrata
da_complexmodel_ch4<- m4_t_vegpresence


#PLot observed vs predicted of two options: 
plot_obs_vs_pred_models(model1 = da_simplemodel_ch4, 
                        model2 = da_complexmodel_ch4,
                        data=da_ch4,
                        color_var = "vegpresence")

ggsave(filename = "DA_ch4_Observed_VS_predicted.png", 
       width = 160,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)


##2.1. DU_ch4 (ok) -------
#Is vegpresence random?
adonis_summary %>% filter(casepilot=="DU")
#YES, vegpresence is not linked to status, we can allow vegpresnce as a fixed effect

#Subset data and check transformation: 
du_ch4<- data4models %>% filter(casepilot=="DU"&ghgspecies=="ch4")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="DU"&
           ghgspecies=="ch4") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = du_ch4,
                           family = gaussian(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m1_gaus_nostrata) #OK
check_singularity(m1_gaus_nostrata,tolerance = 1e-8)
r2(m1_gaus_nostrata, tolerance = 1e-10) 
Anova(m1_gaus_nostrata)
res<- simulateResiduals(m1_gaus_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_nostrata)
em<- emmeans(m1_gaus_nostrata, ~ status,weights = "equal") #using equal weights as we are only averaging across season and all should have the same impact.
summary(em)
pairs(em)
#DECISSION: GOOD residuals


#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                              data = du_ch4,
                              family = gaussian(),
                              dispformula = ~1)
#Evaluate model:
check_convergence(m2_gaus_vegpresence) #OK
check_singularity(m2_gaus_vegpresence,tolerance = 1e-8) #Singular
r2(m2_gaus_vegpresence, tolerance = 1e-10)
Anova(m2_gaus_vegpresence)
summary(m2_gaus_vegpresence)
res<- simulateResiduals(m2_gaus_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m2_gaus_vegpresence)
em<- emmeans(m2_gaus_vegpresence, ~ status, weights = "proportional")
summary(em)
pairs(em)
#DECISSION: GOOD residuals


#SKIP, gaussian already good: 
if(F){
#T_family:
#Most basic: t_family, only status and season 
m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = du_ch4,
                        family = t_family(),
                        dispformula = ~1)
#Evaluate basic model:
check_convergence(m3_t_nostrata)
check_singularity(m3_t_nostrata,tolerance = 1e-8) 
r2(m3_t_nostrata, tolerance=1e-10) 
Anova(m3_t_nostrata)
res<- simulateResiduals(m3_t_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m3_t_nostrata)
em<- emmeans(m3_t_nostrata, ~ status, weights = "equal")
summary(em)
pairs(em)
#DECISSION:


#Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                           data = du_ch4,
                           family = t_family(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m4_t_vegpresence) #OK
check_singularity(m4_t_vegpresence,tolerance = 1e-8)
r2(m4_t_vegpresence, tolerance = 1e-10)
Anova(m4_t_vegpresence)
res<- simulateResiduals(m4_t_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m4_t_vegpresence)
em<- emmeans(m4_t_vegpresence, ~ status, weights = "proportional")
summary(em)
pairs(em)
#DECISSION:

}

#Compare models 
anova(m1_gaus_nostrata,m2_gaus_vegpresence)
#Adding vegpresence does not significantly improve fit (not worth the extra complexity)

#Save best models (simple and complex best)
du_simplemodel_ch4<- m1_gaus_nostrata
du_complexmodel_ch4<- m2_gaus_vegpresence


#PLot observed vs predicted of two options: 
plot_obs_vs_pred_models(model1 = du_simplemodel_ch4, 
                        model2 = du_complexmodel_ch4,
                        data=du_ch4,
                        color_var = "vegpresence")

ggsave(filename = "DU_ch4_Observed_VS_predicted.png", 
       width = 160,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)

##2.1. RI_ch4 (ok) -------
#Is vegpresence random?
adonis_summary %>% filter(casepilot=="RI")
#NO! vegpresence is SIGNIFICANTLY linked to status, we cannot allow vegpresnce as a fixed effect

#Subset data and check transformation: 
ri_ch4<- data4models %>% filter(casepilot=="RI"&ghgspecies=="ch4")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="RI"&
           ghgspecies=="ch4") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = ri_ch4,
                           family = gaussian(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m1_gaus_nostrata) #OK
check_singularity(m1_gaus_nostrata,tolerance = 1e-8)
r2(m1_gaus_nostrata, tolerance = 1e-10) 
Anova(m1_gaus_nostrata)
summary(m1_gaus_nostrata)
res<- simulateResiduals(m1_gaus_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_nostrata)
em<- emmeans(m1_gaus_nostrata, ~ status,weights = "equal") #using equal weights as we are only averaging across season and all should have the same impact.
summary(em)
pairs(em)
#DECISSION: BAD residuals


#SKIP! cannot include vegpresence as fixed with status
if(F){
  #Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
  m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                                data = ri_ch4,
                                family = gaussian(),
                                dispformula = ~1)
  #Evaluate model:
  check_convergence(m2_gaus_vegpresence) #OK
  check_singularity(m2_gaus_vegpresence,tolerance = 1e-8) #
  r2(m2_gaus_vegpresence, tolerance = 1e-10)
  Anova(m2_gaus_vegpresence)
  summary(m2_gaus_vegpresence)
  res<- simulateResiduals(m2_gaus_vegpresence)
  plotQQunif(res)
  plotResiduals(res)
  check_residuals(m2_gaus_vegpresence)
  em<- emmeans(m2_gaus_vegpresence, ~ status, weights = "proportional")
  summary(em)
  pairs(em)
  #DECISSION: 
}


#SKIP t_options: 

  #T_family:
  #Most basic: t_family, only status and season 
  m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                          data = ri_ch4,
                          family = t_family(),
                          dispformula = ~1)
  #Evaluate basic model:
  check_convergence(m3_t_nostrata)
  check_singularity(m3_t_nostrata,tolerance = 1e-8) 
  r2(m3_t_nostrata, tolerance=1e-10) 
  Anova(m3_t_nostrata)
  res<- simulateResiduals(m3_t_nostrata)
  plotQQunif(res)
  plotResiduals(res)
  check_residuals(m3_t_nostrata)
  em<- emmeans(m3_t_nostrata, ~ status, weights = "equal")
  summary(em)
  pairs(em)
  #DECISSION: GOOD residuals
  
  #SKIP, cannot include vegpresence as fixed effect
  if (F){
  #Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
  m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                             data = ri_ch4,
                             family = t_family(),
                             dispformula = ~1)
  #Evaluate basic model:
  check_convergence(m4_t_vegpresence)
  check_singularity(m4_t_vegpresence,tolerance = 1e-8)
  r2(m4_t_vegpresence, tolerance = 1e-10)
  Anova(m4_t_vegpresence)
  summary(m4_t_vegpresence)
  res<- simulateResiduals(m4_t_vegpresence)
  plotQQunif(res)
  plotResiduals(res)
  check_residuals(m4_t_vegpresence)
  em<- emmeans(m4_t_vegpresence, ~ status, weights = "proportional")
  summary(em)
  pairs(em)
  #DECISSION: 
}


#Compare models 
# anova(m1_gaus_nostrata,m3_t_nostrata)

#CANNOT HAVE VEGPRESENCE AS EFFECT FOR RIA DE AVEIRO, most basic model is best (and only allowed)

  #Save best models (simple and complex best)
ri_simplemodel_ch4<- m3_t_nostrata
# ri_complexmodel_ch4<- NA #NOt possible


#PLot observed vs predicted ONLY option: 
plot_obs_vs_pred_models(model1 = ri_simplemodel_ch4, 
                        # model2 = ri_complexmodel_ch4,
                        data=ri_ch4,
                        color_var = "vegpresence")

ggsave(filename = "RI_ch4_Observed_VS_predicted.png", 
       width = 160,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)

##2.1. VA_ch4 (ok) -------

#Is vegpresence random?
adonis_summary %>% filter(casepilot=="VA")
#YES, vegpresence is not linked to status, we can allow vegpresnce as a fixed effect

#Subset data and check transformation: 
va_ch4<- data4models %>% filter(casepilot=="VA"&ghgspecies=="ch4")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="VA"&
           ghgspecies=="ch4") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = va_ch4,
                           family = gaussian(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m1_gaus_nostrata) #OK
check_singularity(m1_gaus_nostrata,tolerance = 1e-8)
r2(m1_gaus_nostrata, tolerance = 1e-10) 
Anova(m1_gaus_nostrata)
res<- simulateResiduals(m1_gaus_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_nostrata)
em<- emmeans(m1_gaus_nostrata, ~ status,weights = "equal") #using equal weights as we are only averaging across season and all should have the same impact.
summary(em)
pairs(em)
#DECISSION: GOOD residuals 


#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                              data = va_ch4,
                              family = gaussian(),
                              dispformula = ~1)
#Evaluate model:
check_convergence(m2_gaus_vegpresence) #OK
check_singularity(m2_gaus_vegpresence,tolerance = 1e-8) #Singular
r2(m2_gaus_vegpresence, tolerance = 1e-10)
Anova(m2_gaus_vegpresence)
summary(m2_gaus_vegpresence)
res<- simulateResiduals(m2_gaus_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m2_gaus_vegpresence)
em<- emmeans(m2_gaus_vegpresence, ~ status, weights = "proportional")
summary(em)
pairs(em)
#DECISSION: GOOD residuals


#SKIP, Gaussian already good
if(F){
#T_family:
#Most basic: t_family, only status and season 
m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = va_ch4,
                        family = t_family(),
                        dispformula = ~1)
#Evaluate basic model:
check_convergence(m3_t_nostrata)
check_singularity(m3_t_nostrata,tolerance = 1e-8) 
r2(m3_t_nostrata, tolerance=1e-10) 
Anova(m3_t_nostrata)
res<- simulateResiduals(m3_t_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m3_t_nostrata)
em<- emmeans(m3_t_nostrata, ~ status, weights = "equal")
summary(em)
pairs(em)
#DECISSION:


#Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                           data = va_ch4,
                           family = t_family(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m4_t_vegpresence) #OK
check_singularity(m4_t_vegpresence,tolerance = 1e-8)
r2(m4_t_vegpresence, tolerance = 1e-10)
Anova(m4_t_vegpresence)
res<- simulateResiduals(m4_t_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m4_t_vegpresence)
em<- emmeans(m4_t_vegpresence, ~ status, weights = "proportional")
summary(em)
pairs(em)
#DECISSION:
}

#Compare models 
anova(m1_gaus_nostrata,m2_gaus_vegpresence)
#Adding vegpresence is significantly better (better fit worth complexity increase)

#Save best models (simple and complex best)
va_simplemodel_ch4<- m1_gaus_nostrata
va_complexmodel_ch4<- m2_gaus_vegpresence


#PLot observed vs predicted of two options: 
plot_obs_vs_pred_models(model1 = va_simplemodel_ch4, 
                        model2 = va_complexmodel_ch4,
                        data=va_ch4,
                        color_var = "vegpresence")

ggsave(filename = "VA_ch4_Observed_VS_predicted.png", 
       width = 160,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)




#GWP______-------

##2.1. CA_GWPco2andch4 (ok) -------

#Is vegpresence random?
adonis_summary %>% filter(casepilot=="CA")
#YES, vegpresence is not significantly linked to status, we can allow vegpresnce as a fixed effect

#Subset data and check transformation: 
ca_GWPco2andch4<- data4models %>% filter(casepilot=="CA"&ghgspecies=="GWPco2andch4")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="CA"&
           ghgspecies=="GWPco2andch4") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = ca_GWPco2andch4,
                           family = gaussian(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m1_gaus_nostrata) #OK
check_singularity(m1_gaus_nostrata,tolerance = 1e-8) #OK
r2(m1_gaus_nostrata, tolerance = 1e-10) 
Anova(m1_gaus_nostrata)
res<- simulateResiduals(m1_gaus_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_nostrata)
em<- emmeans(m1_gaus_nostrata, ~ status,weights = "equal") #using equal weights as we are only averaging across season and all should have the same impact.
summary(em)
pairs(em)
#DECISSION: Good residuals, keep



#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                              data = ca_GWPco2andch4,
                              family = gaussian(),
                              dispformula = ~1)
#Evaluate model:
check_convergence(m2_gaus_vegpresence) #OK
check_singularity(m2_gaus_vegpresence,tolerance = 1e-8) #Singular
r2(m2_gaus_vegpresence, tolerance = 1e-10)
Anova(m2_gaus_vegpresence)
res<- simulateResiduals(m2_gaus_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m2_gaus_vegpresence)
em<- emmeans(m2_gaus_vegpresence, ~ status, weights = "proportional")
summary(em)
pairs(em)
#DECISSION: BAD residuals


  #T_family:
  #Most basic: t_family, only status and season 
  m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                          data = ca_GWPco2andch4,
                          family = t_family(),
                          dispformula = ~1)
  #Evaluate basic model:
  check_convergence(m3_t_nostrata)
  check_singularity(m3_t_nostrata,tolerance = 1e-8) 
  r2(m3_t_nostrata, tolerance=1e-10) 
  Anova(m3_t_nostrata)
  res<- simulateResiduals(m3_t_nostrata)
  plotQQunif(res)
  plotResiduals(res)
  check_residuals(m3_t_nostrata)
  em<- emmeans(m3_t_nostrata, ~ status)
  summary(em)
  pairs(em)
  #DECISSION: Good residuals
  
  
  #Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
  m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                             data = ca_GWPco2andch4,
                             family = t_family(),
                             dispformula = ~1)
  #Evaluate basic model:
  check_convergence(m4_t_vegpresence) #OK
  check_singularity(m4_t_vegpresence,tolerance = 1e-8)
  r2(m4_t_vegpresence, tolerance = 1e-9)
  Anova(m4_t_vegpresence)
  res<- simulateResiduals(m4_t_vegpresence)
  plotQQunif(res)
  plotResiduals(res)
  check_residuals(m4_t_vegpresence)
  check_overdispersion(m4_t_vegpresence)
  em<- emmeans(m4_t_vegpresence, ~ status, weights = "proportional")
  summary(em)
  pairs(em)
  #DECISSION: Good  residuals
  

#Compare models 
anova(m1_gaus_nostrata,m4_t_vegpresence)
#Adding vegpresence does significantly improve results

#Save best models (simple and complex best)
ca_simplemodel_GWPco2andch4<- m1_gaus_nostrata
ca_complexmodel_GWPco2andch4<- m4_t_vegpresence


#PLot observed vs predicted of two options: 
plot_obs_vs_pred_models(model1 = ca_simplemodel_GWPco2andch4, 
                        model2 = ca_complexmodel_GWPco2andch4,
                        data=ca_GWPco2andch4,
                        color_var = "vegpresence")

ggsave(filename = "CA_GWPco2andch4_Observed_VS_predicted.png", 
       width = 160,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)




##2.2. CU_GWPco2andch4 (ok) -------

#Is vegpresence random?
adonis_summary %>% filter(casepilot=="CU")
#YES, by definition, vegpresence is not linked to status, we can allow vegpresnce as a fixed effect

#Subset data and check transformation: 
cu_GWPco2andch4<- data4models %>% filter(casepilot=="CU"&ghgspecies=="GWPco2andch4")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="CU"&
           ghgspecies=="GWPco2andch4") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = cu_GWPco2andch4,
                           family = gaussian(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m1_gaus_nostrata) #OK
check_singularity(m1_gaus_nostrata,tolerance = 1e-8)
r2(m1_gaus_nostrata, tolerance = 1e-10) 
Anova(m1_gaus_nostrata)
summary(m1_gaus_nostrata)
res<- simulateResiduals(m1_gaus_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_nostrata)
em<- emmeans(m1_gaus_nostrata, ~ status,weights = "equal") #using equal weights as we are only averaging across season and all should have the same impact.
summary(em)
pairs(em)
#DECISSION: BAD residuals


#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                              data = cu_GWPco2andch4,
                              family = gaussian(),
                              dispformula = ~1)
#Evaluate model:
check_convergence(m2_gaus_vegpresence) #OK
check_singularity(m2_gaus_vegpresence,tolerance = 1e-8) #Singular
r2(m2_gaus_vegpresence, tolerance = 1e-10)
Anova(m2_gaus_vegpresence)
summary(m2_gaus_vegpresence)
res<- simulateResiduals(m2_gaus_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m2_gaus_vegpresence)
em<- emmeans(m2_gaus_vegpresence, ~ status, weights = "proportional")
summary(em)
pairs(em)
#DECISSION: BAD residuals


  #T_family:
  #Most basic: t_family, only status and season 
  m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                          data = cu_GWPco2andch4,
                          family = t_family(),
                          dispformula = ~1)
  #Evaluate basic model:
  check_convergence(m3_t_nostrata)
  check_singularity(m3_t_nostrata,tolerance = 1e-8) 
  r2(m3_t_nostrata, tolerance=1e-10) 
  Anova(m3_t_nostrata)
  summary(m3_t_nostrata)
  res<- simulateResiduals(m3_t_nostrata)
  plotQQunif(res)
  plotResiduals(res)
  check_residuals(m3_t_nostrata)
  em<- emmeans(m3_t_nostrata, ~ status)
  summary(em)
  pairs(em)
  #DECISSION: Good enough residuals
  
  
  #Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
  m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                             data = cu_GWPco2andch4,
                             family = t_family(),
                             dispformula = ~1)
  #Evaluate basic model:
  check_convergence(m4_t_vegpresence) #OK
  check_singularity(m4_t_vegpresence,tolerance = 1e-8)
  r2(m4_t_vegpresence, tolerance = 1e-9)
  Anova(m4_t_vegpresence)
  summary(m4_t_vegpresence)
  res<- simulateResiduals(m4_t_vegpresence)
  plotQQunif(res)
  plotResiduals(res)
  check_residuals(m4_t_vegpresence)
  em<- emmeans(m4_t_vegpresence, ~ status, weights = "proportional")
  summary(em)
  pairs(em)
  #DECISSION: GOOD residuals

#Compare models 
anova(m3_t_nostrata,m4_t_vegpresence)
#Adding vegpresence is significantly better (better fit worth complexity increase)

#Save best models (simple and complex best)
cu_simplemodel_GWPco2andch4<- m3_t_nostrata
cu_complexmodel_GWPco2andch4<- m4_t_vegpresence


#PLot observed vs predicted of two options: 
plot_obs_vs_pred_models(model1 = cu_simplemodel_GWPco2andch4, 
                        model2 = cu_complexmodel_GWPco2andch4,
                        data=cu_GWPco2andch4,
                        color_var = "vegpresence")

ggsave(filename = "CU_GWPco2andch4_Observed_VS_predicted.png", 
       width = 160,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)

##2.3. DA_GWPco2andch4 (ok-ish) -------

#Is vegpresence random?
adonis_summary %>% filter(casepilot=="DA")
#YES, vegpresence is not linked to status, we can allow vegpresnce as a fixed effect

#Subset data and check transformation: 
da_GWPco2andch4<- data4models %>% filter(casepilot=="DA"&ghgspecies=="GWPco2andch4")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="DA"&
           ghgspecies=="GWPco2andch4") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = da_GWPco2andch4,
                           family = gaussian(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m1_gaus_nostrata) #OK
check_singularity(m1_gaus_nostrata,tolerance = 1e-8)
r2(m1_gaus_nostrata, tolerance = 1e-10) 
Anova(m1_gaus_nostrata)
summary(m1_gaus_nostrata)
res<- simulateResiduals(m1_gaus_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_nostrata)
em<- emmeans(m1_gaus_nostrata, ~ status,weights = "equal") #using equal weights as we are only averaging across season and all should have the same impact.
summary(em)
pairs(em)
#DECISSION: Bad residuals


#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                              data = da_GWPco2andch4,
                              family = gaussian(),
                              dispformula = ~1)
#Evaluate model:
check_convergence(m2_gaus_vegpresence) #OK
check_singularity(m2_gaus_vegpresence,tolerance = 1e-8) #Singular
r2(m2_gaus_vegpresence, tolerance = 1e-10)
Anova(m2_gaus_vegpresence)
summary(m2_gaus_vegpresence)
res<- simulateResiduals(m2_gaus_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m2_gaus_vegpresence)
em<- emmeans(m2_gaus_vegpresence, ~ status, weights = "proportional")
summary(em)
pairs(em)
#DECISSION: BAD residuals



#T_family:
#Most basic: t_family, only status and season 
m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = da_GWPco2andch4,
                        family = t_family(),
                        dispformula = ~1)
#Evaluate basic model:
check_convergence(m3_t_nostrata)
check_singularity(m3_t_nostrata,tolerance = 1e-8) 
r2(m3_t_nostrata, tolerance=1e-11) 
Anova(m3_t_nostrata)
summary(m3_t_nostrata)
res<- simulateResiduals(m3_t_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m3_t_nostrata)
em<- emmeans(m3_t_nostrata, ~ status, weights = "equal")
summary(em)
pairs(em)
#DECISSION: Better residuals, although they still fail


#Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                           data = da_GWPco2andch4,
                           family = t_family(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m4_t_vegpresence) #OK
check_singularity(m4_t_vegpresence,tolerance = 1e-8)
r2(m4_t_vegpresence, tolerance = 1e-11)
Anova(m4_t_vegpresence)
summary(m4_t_vegpresence)
res<- simulateResiduals(m4_t_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m4_t_vegpresence)
em<- emmeans(m4_t_vegpresence, ~ status, weights = "proportional")
summary(em)
pairs(em)
#Better residuals, although they present issues

#Compare models 
anova(m3_t_nostrata,m4_t_vegpresence)
#Adding vegpresence is significantly better (better fit worth complexity increase)

#Save best models (simple and complex best)
da_simplemodel_GWPco2andch4<- m3_t_nostrata
da_complexmodel_GWPco2andch4<- m4_t_vegpresence


#PLot observed vs predicted of two options: 
plot_obs_vs_pred_models(model1 = da_simplemodel_GWPco2andch4, 
                        model2 = da_complexmodel_GWPco2andch4,
                        data=da_GWPco2andch4,
                        color_var = "vegpresence")

ggsave(filename = "DA_GWPco2andch4_Observed_VS_predicted.png", 
       width = 160,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)


##2.1. DU_GWPco2andch4 (ok) -------
#Is vegpresence random?
adonis_summary %>% filter(casepilot=="DU")
#YES, vegpresence is not linked to status, we can allow vegpresnce as a fixed effect

#Subset data and check transformation: 
du_GWPco2andch4<- data4models %>% filter(casepilot=="DU"&ghgspecies=="GWPco2andch4")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="DU"&
           ghgspecies=="GWPco2andch4") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = du_GWPco2andch4,
                           family = gaussian(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m1_gaus_nostrata) #OK
check_singularity(m1_gaus_nostrata,tolerance = 1e-8)
r2(m1_gaus_nostrata, tolerance = 1e-10) 
Anova(m1_gaus_nostrata)
res<- simulateResiduals(m1_gaus_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_nostrata)
em<- emmeans(m1_gaus_nostrata, ~ status,weights = "equal") #using equal weights as we are only averaging across season and all should have the same impact.
summary(em)
pairs(em)
#DECISSION: Bad residuals


#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                              data = du_GWPco2andch4,
                              family = gaussian(),
                              dispformula = ~1)
#Evaluate model:
check_convergence(m2_gaus_vegpresence) #OK
check_singularity(m2_gaus_vegpresence,tolerance = 1e-8) #Singular
r2(m2_gaus_vegpresence, tolerance = 1e-10)
Anova(m2_gaus_vegpresence)
summary(m2_gaus_vegpresence)
res<- simulateResiduals(m2_gaus_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m2_gaus_vegpresence)
em<- emmeans(m2_gaus_vegpresence, ~ status, weights = "proportional")
summary(em)
pairs(em)
#DECISSION: GOOD-enough residuals


  #T_family:
  #Most basic: t_family, only status and season 
  m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                          data = du_GWPco2andch4,
                          family = t_family(),
                          dispformula = ~1)
  #Evaluate basic model:
  check_convergence(m3_t_nostrata)
  check_singularity(m3_t_nostrata,tolerance = 1e-8) 
  r2(m3_t_nostrata, tolerance=1e-10) 
  Anova(m3_t_nostrata)
  res<- simulateResiduals(m3_t_nostrata)
  plotQQunif(res)
  plotResiduals(res)
  check_residuals(m3_t_nostrata)
  em<- emmeans(m3_t_nostrata, ~ status, weights = "equal")
  summary(em)
  pairs(em)
  #DECISSION: DOES NOT CONVERGE
  
  
  #Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
  m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                             data = du_GWPco2andch4,
                             family = t_family(),
                             dispformula = ~1)
  #Evaluate basic model:
  check_convergence(m4_t_vegpresence) #OK
  check_singularity(m4_t_vegpresence,tolerance = 1e-8)
  r2(m4_t_vegpresence, tolerance = 1e-10)
  Anova(m4_t_vegpresence)
  res<- simulateResiduals(m4_t_vegpresence)
  plotQQunif(res)
  plotResiduals(res)
  check_residuals(m4_t_vegpresence)
  em<- emmeans(m4_t_vegpresence, ~ status, weights = "proportional")
  summary(em)
  pairs(em)
  #DECISSION: 

  

#Compare models 
anova(m1_gaus_nostrata,m2_gaus_vegpresence)
#Adding vegpresence does significantly improve fit (worth the extra complexity)

#Save best models (simple and complex best)
du_simplemodel_GWPco2andch4<- m1_gaus_nostrata
du_complexmodel_GWPco2andch4<- m2_gaus_vegpresence


#PLot observed vs predicted of two options: 
plot_obs_vs_pred_models(model1 = du_simplemodel_GWPco2andch4, 
                        model2 = du_complexmodel_GWPco2andch4,
                        data=du_GWPco2andch4,
                        color_var = "vegpresence")

ggsave(filename = "DU_GWPco2andch4_Observed_VS_predicted.png", 
       width = 160,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)

##2.1. RI_GWPco2andch4 (ok) -------
#Is vegpresence random?
adonis_summary %>% filter(casepilot=="RI")
#NO! vegpresence is SIGNIFICANTLY linked to status, we cannot allow vegpresnce as a fixed effect

#Subset data and check transformation: 
ri_GWPco2andch4<- data4models %>% filter(casepilot=="RI"&ghgspecies=="GWPco2andch4")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="RI"&
           ghgspecies=="GWPco2andch4") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = ri_GWPco2andch4,
                           family = gaussian(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m1_gaus_nostrata) #OK
check_singularity(m1_gaus_nostrata,tolerance = 1e-8)
r2(m1_gaus_nostrata, tolerance = 1e-10) 
Anova(m1_gaus_nostrata)
summary(m1_gaus_nostrata)
res<- simulateResiduals(m1_gaus_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_nostrata)
em<- emmeans(m1_gaus_nostrata, ~ status,weights = "equal") #using equal weights as we are only averaging across season and all should have the same impact.
summary(em)
pairs(em)
#DECISSION: GOOD residuals


#SKIP! cannot include vegpresence as fixed with status
if(F){
  #Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
  m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                                data = ri_GWPco2andch4,
                                family = gaussian(),
                                dispformula = ~1)
  #Evaluate model:
  check_convergence(m2_gaus_vegpresence) #OK
  check_singularity(m2_gaus_vegpresence,tolerance = 1e-8) #
  r2(m2_gaus_vegpresence, tolerance = 1e-10)
  Anova(m2_gaus_vegpresence)
  summary(m2_gaus_vegpresence)
  res<- simulateResiduals(m2_gaus_vegpresence)
  plotQQunif(res)
  plotResiduals(res)
  check_residuals(m2_gaus_vegpresence)
  em<- emmeans(m2_gaus_vegpresence, ~ status, weights = "proportional")
  summary(em)
  pairs(em)
  #DECISSION: 



#SKIP t_options: 

#T_family:
#Most basic: t_family, only status and season 
m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = ri_GWPco2andch4,
                        family = t_family(),
                        dispformula = ~1)
#Evaluate basic model:
check_convergence(m3_t_nostrata)
check_singularity(m3_t_nostrata,tolerance = 1e-8) 
r2(m3_t_nostrata, tolerance=1e-10) 
Anova(m3_t_nostrata)
res<- simulateResiduals(m3_t_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m3_t_nostrata)
em<- emmeans(m3_t_nostrata, ~ status, weights = "equal")
summary(em)
pairs(em)
#DECISSION: GOOD residuals

#SKIP, cannot include vegpresence as fixed effect
  #Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
  m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                             data = ri_GWPco2andch4,
                             family = t_family(),
                             dispformula = ~1)
  #Evaluate basic model:
  check_convergence(m4_t_vegpresence)
  check_singularity(m4_t_vegpresence,tolerance = 1e-8)
  r2(m4_t_vegpresence, tolerance = 1e-10)
  Anova(m4_t_vegpresence)
  summary(m4_t_vegpresence)
  res<- simulateResiduals(m4_t_vegpresence)
  plotQQunif(res)
  plotResiduals(res)
  check_residuals(m4_t_vegpresence)
  em<- emmeans(m4_t_vegpresence, ~ status, weights = "proportional")
  summary(em)
  pairs(em)
  #DECISSION: 
}

#Compare models: gaussian simple is already the best 
# anova(m1_gaus_nostrata,m3_t_nostrata)

#CANNOT HAVE VEGPRESENCE AS EFFECT FOR RIA DE AVEIRO, most basic model is best (and only allowed)

#Save best models (simple and complex best)
ri_simplemodel_GWPco2andch4<- m1_gaus_nostrata
# ri_complexmodel_GWPco2andch4<- NA # Complex not possible


#PLot observed vs predicted ONLY option: 
plot_obs_vs_pred_models(model1 = ri_simplemodel_GWPco2andch4, 
                        # model2 = ri_complexmodel_GWPco2andch4,
                        data=ri_GWPco2andch4,
                        color_var = "vegpresence")

ggsave(filename = "RI_GWPco2andch4_Observed_VS_predicted.png", 
       width = 160,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)

##2.1. VA_GWPco2andch4 (ok) -------

#Is vegpresence random?
adonis_summary %>% filter(casepilot=="VA")
#YES, vegpresence is not linked to status, we can allow vegpresnce as a fixed effect

#Subset data and check transformation: 
va_GWPco2andch4<- data4models %>% filter(casepilot=="VA"&ghgspecies=="GWPco2andch4")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="VA"&
           ghgspecies=="GWPco2andch4") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = va_GWPco2andch4,
                           family = gaussian(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m1_gaus_nostrata) #OK
check_singularity(m1_gaus_nostrata,tolerance = 1e-8)
r2(m1_gaus_nostrata, tolerance = 1e-10) 
Anova(m1_gaus_nostrata)
res<- simulateResiduals(m1_gaus_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_nostrata)
em<- emmeans(m1_gaus_nostrata, ~ status,weights = "equal") #using equal weights as we are only averaging across season and all should have the same impact.
summary(em)
pairs(em)
#DECISSION: BAD residuals 


#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                              data = va_GWPco2andch4,
                              family = gaussian(),
                              dispformula = ~1)
#Evaluate model:
check_convergence(m2_gaus_vegpresence) #OK
check_singularity(m2_gaus_vegpresence,tolerance = 1e-8) #Singular
r2(m2_gaus_vegpresence, tolerance = 1e-10)
Anova(m2_gaus_vegpresence)
summary(m2_gaus_vegpresence)
res<- simulateResiduals(m2_gaus_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m2_gaus_vegpresence)
em<- emmeans(m2_gaus_vegpresence, ~ status, weights = "proportional")
summary(em)
pairs(em)
#DECISSION: BAD residuals 


  #T_family:
  #Most basic: t_family, only status and season 
  m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                          data = va_GWPco2andch4,
                          family = t_family(),
                          dispformula = ~1)
  #Evaluate basic model:
  check_convergence(m3_t_nostrata)
  check_singularity(m3_t_nostrata,tolerance = 1e-8) 
  r2(m3_t_nostrata, tolerance=1e-10) 
  Anova(m3_t_nostrata)
  res<- simulateResiduals(m3_t_nostrata)
  plotQQunif(res)
  plotResiduals(res)
  check_residuals(m3_t_nostrata)
  em<- emmeans(m3_t_nostrata, ~ status, weights = "equal")
  summary(em)
  pairs(em)
  #DECISSION: DOES NOT CONVERGE
  
  
  #Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
  m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                             data = va_GWPco2andch4,
                             family = t_family(),
                             dispformula = ~1)
  #Evaluate basic model:
  check_convergence(m4_t_vegpresence) #OK
  check_singularity(m4_t_vegpresence,tolerance = 1e-8)
  r2(m4_t_vegpresence, tolerance = 1e-10)
  Anova(m4_t_vegpresence)
  res<- simulateResiduals(m4_t_vegpresence)
  plotQQunif(res)
  plotResiduals(res)
  check_residuals(m4_t_vegpresence)
  em<- emmeans(m4_t_vegpresence, ~ status, weights = "proportional")
  summary(em)
  pairs(em)
  #DECISSION:GOOD residuals

#Compare models 
anova(m1_gaus_nostrata,m4_t_vegpresence)
#Adding vegpresence is significantly better (better fit worth complexity increase)

#Save best models (simple and complex best)
va_simplemodel_GWPco2andch4<- m1_gaus_nostrata
va_complexmodel_GWPco2andch4<- m2_gaus_vegpresence


#PLot observed vs predicted of two options: 
plot_obs_vs_pred_models(model1 = va_simplemodel_GWPco2andch4, 
                        model2 = va_complexmodel_GWPco2andch4,
                        data=va_GWPco2andch4,
                        color_var = "vegpresence")

ggsave(filename = "VA_GWPco2andch4_Observed_VS_predicted.png", 
       width = 160,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)


#3. Save bestmodel-list -------

#Save lists of models: 1. simple_model_list (only status*season models), all casepilots, 2. complex_model_list (status*season*vegpresence), only casepilots with complex_model (i.e. no RI-models). Both lists named as: casepilot_ghgspecies


#List of simple models for each casepilot*ghgspecies 
simplemodel_list_allghg<- list(
  "CA_co2" = ca_simplemodel_co2,
  "CU_co2" = cu_simplemodel_co2, 
  "DA_co2" = da_simplemodel_co2,
  "DU_co2" = du_simplemodel_co2,
  "RI_co2" = ri_simplemodel_co2,
  "VA_co2" = va_simplemodel_co2,
  "CA_ch4" = ca_simplemodel_ch4,
  "CU_ch4" = cu_simplemodel_ch4, 
  "DA_ch4" = da_simplemodel_ch4,
  "DU_ch4" = du_simplemodel_ch4,
  "RI_ch4" = ri_simplemodel_ch4,
  "VA_ch4" = va_simplemodel_ch4,
  "CA_GWPco2andch4" = ca_simplemodel_GWPco2andch4,
  "CU_GWPco2andch4" = cu_simplemodel_GWPco2andch4, 
  "DA_GWPco2andch4" = da_simplemodel_GWPco2andch4,
  "DU_GWPco2andch4" = du_simplemodel_GWPco2andch4,
  "RI_GWPco2andch4" = ri_simplemodel_GWPco2andch4,
  "VA_GWPco2andch4" = va_simplemodel_GWPco2andch4
)


#List of complex models for each casepilot*ghgspecies 
complexmodel_list_allghg<- list(
  "CA_co2" = ca_complexmodel_co2,
  "CU_co2" = cu_complexmodel_co2, 
  "DA_co2" = da_complexmodel_co2,
  "DU_co2" = du_complexmodel_co2,
  # "RI_co2" = ri_complexmodel_co2, #NOT possible for RI
  "VA_co2" = va_complexmodel_co2,
  "CA_ch4" = ca_complexmodel_ch4,
  "CU_ch4" = cu_complexmodel_ch4, 
  "DA_ch4" = da_complexmodel_ch4,
  "DU_ch4" = du_complexmodel_ch4,
  # "RI_ch4" = ri_complexmodel_ch4, #NOT possible for RI
  "VA_ch4" = va_complexmodel_ch4,
  "CA_GWPco2andch4" = ca_complexmodel_GWPco2andch4,
  "CU_GWPco2andch4" = cu_complexmodel_GWPco2andch4, 
  "DA_GWPco2andch4" = da_complexmodel_GWPco2andch4,
  "DU_GWPco2andch4" = du_complexmodel_GWPco2andch4,
  # "RI_GWPco2andch4" = ri_complexmodel_GWPco2andch4, #NOT possible for RI
  "VA_GWPco2andch4" = va_complexmodel_GWPco2andch4
)


#Save last run environment (avoids recalculation of each model when modifying model result exports)
#Last saved: 
#simplemodels: only consideing status*season
#complexmodels: taking into account vegpresence (status*season*vegpresence). 
save.image(file = paste0(paper_path,"chambermodels.RData"))

#_______----
#Load bestmodels----

load(file = paste0(paper_path, "chambermodels.RData"))




#SUMMARY across GHGs-------


## 3.1. Overall fit -----
#Structure, pseudoR2 of best model:
simplemodel_fit<- get_pseudoR2s(simplemodel_list_allghg)
complexmodel_fit<- get_pseudoR2s(complexmodel_list_allghg)

#R2c and R2m of homocedastic model:
simplemodel_homoc_r2<- get_homoc_R2s(simplemodel_list_allghg)
complexmodel_homoc_r2<- get_homoc_R2s(complexmodel_list_allghg)

#DHARMa Residual diagnostics of best model:
simplemodel_resid_diag_summary <- summarize_dharma_diagnostics(simplemodel_list_allghg)
complexmodel_resid_diag_summary <- summarize_dharma_diagnostics(complexmodel_list_allghg)


##3.2. Significance of effects-----
#Effect of status will first be assessed via Anova (best_model,type=3): is there an effect averaging across all seasons? yes/no, is the effect dependent on season? 

#Obtain significance for main effects of best_models
simplemodel_results_anova<-test_get_all_anova_results(simplemodel_list_allghg)
complexmodel_results_anova<-test_get_all_anova_results(complexmodel_list_allghg)


all_simplemodel_outputs<- simplemodel_results_anova %>% 
  full_join(simplemodel_fit) %>% 
  full_join(simplemodel_homoc_r2 %>% dplyr::select(dataset, homoced_R2m, homoced_R2c))

all_complexmodel_outputs<- complexmodel_results_anova %>% 
  full_join(complexmodel_fit) %>% 
  full_join(complexmodel_homoc_r2 %>% dplyr::select(dataset, homoced_R2m, homoced_R2c))


#Save in dropbox tables with the models summary outputs: 

#Simple models
write.csv(x = all_simplemodel_outputs, file = paste0(paper_path,"Summary_bestsimplemodel_chambermodels.csv"),row.names = F)

write.csv(x = simplemodel_resid_diag_summary, file = paste0(paper_path,"Residualtests_bestsimplemodel_chambermodels.csv"),row.names = F)


#Complex models
write.csv(x = all_complexmodel_outputs, file = paste0(paper_path,"Summary_bestcomplexmodel_chambermodels.csv"),row.names = F)

write.csv(x = complexmodel_resid_diag_summary, file = paste0(paper_path,"Residualtests_bestcomplexmodel_chambermodels.csv"),row.names = F)


##3.3. EMEANs----

#Calculate relevant emmeans and post-hocs tests: separately for simple model list (status*season) and complex model list (status*season*vegpresence)

#NOTE on bootstrapping: there is the possibility to bootstrap the models (i.e. recalculate the model many times with re-sampled data, and calculate emmeans for each model, to be able to produce more robust Confidence Intervals for them, but this is not usually done and computationally intensive. We would have to take into account the model structure when re-sampling the data to ensure balanced resamplings (all factors represented by data). WE will not do this (for the moment).

#NOTE on non-gaussian models: upper and lower confidence limits are given as asymp.LCL/asymp.UCL (in our context we can treat them equally as the ones from gaussian models). Additionally, non-gaussian models will have Inf df (degrees of freedom are calculated differently for T-family models), not an issue. 


##Simple models-----
#No need to use weights, we use "equal" to average across seasons with the same weight (assume sampling is representative of every season). 

simple_comparison_list<- list()

# Loop to extract all comparisons for every casepilot best model
for (dataset in names(simplemodel_list_allghg)) {
  #get model
  cp_model <- simplemodel_list_allghg[[dataset]]
  
  #extract casepilot name and ghgspecies from model-list names
  casepilot_name<- sub("_.*", "", dataset)
  ghgspecies<- sub(paste0(casepilot_name,"_"),"",dataset)
  
  # Status comparisons:
  status_emmeans <- emmeans(cp_model, ~status, weights = "equal")
  status_comparisons <- cld(status_emmeans, 
                            alpha = 0.05,
                            Letters = letters,
                            adjust = "sidak") %>%
    mutate(comparison = "status",
           ghgspecies = ghgspecies,
           casepilot = casepilot_name,
           season = NA) %>%
    rename(lower.CL = any_of("asymp.LCL"),
           upper.CL = any_of("asymp.UCL"))
  
  # Season comparisons: 
  season_emmeans <- emmeans(cp_model, ~season, weights = "equal")
  season_comparisons <- cld(season_emmeans, 
                            alpha = 0.05,
                            Letters = letters,
                            adjust = "sidak") %>%
    mutate(comparison = "season",
           ghgspecies = ghgspecies,
           casepilot = casepilot_name,
           status = NA) %>%
    rename(lower.CL = any_of("asymp.LCL"),
           upper.CL = any_of("asymp.UCL"))
  
  #Status_within_season comparisons: 
  statuswithinseason_emmeans <- emmeans(cp_model, ~status|season, weights="equal")
  statuswithinseason_comparisons <- cld(statuswithinseason_emmeans, 
                                 alpha = 0.05,
                                 Letters = letters,
                                 adjust = "sidak") %>%
    mutate(comparison = "status_within_season",
           ghgspecies = ghgspecies,
           casepilot = casepilot_name) %>%
    rename(lower.CL = any_of("asymp.LCL"),
           upper.CL = any_of("asymp.UCL"))
  
  # Join and format CP comparisons
  all_comparisons <- status_comparisons %>%
    full_join(season_comparisons) %>%
    full_join(statuswithinseason_comparisons) %>%
    rename(group_letter = .group) %>%
    mutate(group_letter = gsub(" ", "", group_letter),
           # Create key for accessing BestNormalize object
           key_trans = paste(casepilot_name, ghgspecies, sep = "_")
    ) %>%
    rowwise() %>% 
    # BACK-TRANSFORMATION of emmean, SE and CL using BestNormalize object
    mutate(
      emmean_bt = if (!is.null(bn_list[[key_trans]])) predict(bn_list[[key_trans]], emmean, inverse = TRUE) else NA_real_,
      SE_bt = if (!is.null(bn_list[[key_trans]])) predict(bn_list[[key_trans]], SE, inverse = TRUE) else NA_real_,
      lower.CL_bt = if (!is.null(bn_list[[key_trans]])) predict(bn_list[[key_trans]], lower.CL, inverse = TRUE) else NA_real_,
      upper.CL_bt = if (!is.null(bn_list[[key_trans]])) predict(bn_list[[key_trans]], upper.CL, inverse = TRUE) else NA_real_
    ) %>%
    dplyr::select(any_of(c(
      "ghgspecies", "casepilot","comparison", "status", "season","strata", "df",
      "emmean", "SE", "lower.CL", "upper.CL", "group_letter",
      "emmean_bt", "SE_bt", "lower.CL_bt", "upper.CL_bt"
    ))) %>%
    arrange(ghgspecies, casepilot, comparison, status, season)
  
  # Store casepilot results
  simple_comparison_list[[dataset]] <- all_comparisons
}

#Remove within-loop objects
rm(all_comparisons, status_comparisons,season_comparisons,statuswithinseason_comparisons,
   statuswithinseason_emmeans,status_emmeans, season_emmeans)


# Unir todos los resultados en un solo data.frame
simplemodel_emmeans_withgroups <- bind_rows(simple_comparison_list)

#Save all comparisons:  
write.csv(x = simplemodel_emmeans_withgroups, file = paste0(paper_path,"Emmeans-posthoc_simplemodels_chambermodels.csv"),row.names = F)



##Complex models------

#We need to take into account the proportion of vegpresence in the field to estimate the status, season and status_within_season emmeans. For status_within_vegpresence, weights should not be applied. 

#NOTES: 
#Overall status: weights should be used to account for different vegpresence in different status and potentially in different season, but all seasons should have the same combined weight (to give the same importance to each season regardless of how many observation we have in each.)

formula(cu_best_ch4)
Anova(cu_best_ch4)
emmeans(cu_best_ch4, ~status,weights = "show.levels")
#Calling emmeans with weights="show.levels" give us the combination of levels over which each status emmean is averaging across. 

#IF  we simply use the default weights="equal", each combination of levels above will be treated with equal weight for all status, which will biass our estimates if for example, not every status have the same proportion of vegpresence, or if this vegproportion is different across seasons, or even just if the actual vegproportion in the field is not 50:50.  
emmeans(cu_best_ch4, ~status, weights = "equal")

#WE NEED to use WEIGHTS. 

#A passable alternative would be to use weights="proportional", this option scales each level combination to its proportion in the data with is ok-ish given our semi-proportional sampling strategy. HOWEVER, with this method we get the same "sampling biass" as we have in our data: strata with little cover will be over-represented due to the minimum 3-chambers per strata rule. Also, if the number of measurements differs by season, this imbalance will also be transferred to the estimates. 
emmeans(cu_best_ch4, ~status, weights = "proportional")

#HOWEVER, if we want to "revert" the sampling biasses that we have (by design or by accident), we shoudl use weights that inform of the actual distribution of vegpresence in each season for each status. 
#The optimal solution is to re-calculate the weights for each combination of status*season*vegpresence according to the independently assessed coverage of each site. Additionally, we can ensure that all seasons are given the same importance to account for cases where, by chance, we performed more measurements in one season or other. 



##Opt1. (NO) Using proportional weigths------
#Not ideal
#This method takes into account the different N of the actual data the model was fitted to. So it would preserve any proportionality and biasses of the dataset. 
#BUT it uses the same weight for all levels of the effect being estimated. Cross-status average for vegpresence and season.
summary(cp_model)
cp_model

# Initialize list for storing comparisons
complex_comparisons_list <- list()

# Loop to extract all comparisons for every casepilot best model
for (dataset in names(complexmodel_list_allghg)) {
  cp_model <- complexmodel_list_allghg[[dataset]]
  
  #extract list of fixed effects of model:
  fixed_effects_used<- attr(terms(cp_model), "term.labels")
  
  #extract casepilot name and ghgspecies from model-list names
  casepilot_name<- sub("_.*", "", dataset)
  ghgspecies<- sub(paste0(casepilot_name,"_"),"",dataset)
  
  
  # Status comparisons: 
  status_emmeans <- emmeans(cp_model, ~status, weights = "proportional") #Calculate overall status emmean, scaling results to database. This assumes perfect representativity of our sampling. Eg. If more samplings in winter--> more weight of winter in status emmean (same weight for all levels of status). If predominantly vegetated (across all status)--> more weight to vegetated (same weight for all). 
  status_comparisons <- cld(status_emmeans, 
                            alpha = 0.05,
                            Letters = letters,
                            adjust = "sidak") %>%
    mutate(comparison = "status",
           ghgspecies = ghgspecies,
           casepilot = casepilot_name,
           season = NA) %>%
    rename(lower.CL = any_of("asymp.LCL"),
           upper.CL = any_of("asymp.UCL"))
  
  # Season comparisons: 
  season_emmeans <- emmeans(cp_model, ~season, weights = "proportional") #Calculate overall season emmeans, using database weights. IF through the year there are more measurements in altered, altered will have a greater impact in the estimates.
  season_comparisons <- cld(season_emmeans, 
                            alpha = 0.05,
                            Letters = letters,
                            adjust = "sidak") %>%
    mutate(comparison = "season",
           ghgspecies = ghgspecies,
           casepilot = casepilot_name,
           status = NA) %>%
    rename(lower.CL = any_of("asymp.LCL"),
           upper.CL = any_of("asymp.UCL"))
  
  # Status within season comparisons: 
  statuswithinseason_emmeans <- emmeans(cp_model, ~status|season, weights = "proportional")
  statuswithinseason_comparisons <- cld(statuswithinseason_emmeans, 
                                        alpha = 0.05,
                                        Letters = letters,
                                        adjust = "sidak") %>%
    mutate(comparison = "status_within_season",
           ghgspecies = ghgspecies,
           casepilot = casepilot_name) %>%
    rename(lower.CL = any_of("asymp.LCL"),
           upper.CL = any_of("asymp.UCL"))
  
  
  #Status within vegpresence comparisons: 
  statuswithinvegpresence_emmeans<-  emmeans(cp_model, ~status|vegpresence, weights = "equal")#Using equal weights ensures that different proportions of vegpresence between different seasons are not taken into account. The effect does not take into account if spring had more vegetated class than winter to calcualte the overall effect. Is the equivalent to calculating the yearly effect giving equal importance to all seasons. 
  statuswithinvegpresence_comparisons<- cld(statuswithinvegpresence_emmeans,
                                            alpha=0.05,
                                            Letters=letters,
                                            adjust = "sidak") %>% 
    mutate(comparison="status_within_vegpresence",
           ghgspecies=ghgspecies,
           casepilot=casepilot_name) %>% 
    rename(lower.CL = any_of("asymp.LCL"),
           upper.CL = any_of("asymp.UCL"))
  
  
  # Join and format CP comparisons
  all_comparisons <- status_comparisons %>%
    full_join(season_comparisons) %>%
    full_join(statuswithinseason_comparisons) %>%
    full_join(statuswithinvegpresence_comparisons) %>% 
    rename(group_letter = .group) %>%
    mutate(group_letter = gsub(" ", "", group_letter),
           
           # Create key for accessing BestNormalize object
           key_trans = paste(casepilot_name, ghgspecies, sep = "_")
    ) %>%
    rowwise() %>% 
    # BACK-TRANSFORMATION of emmean, SE and CL using BestNormalize object
    mutate(
      emmean_bt = if (!is.null(bn_list[[key_trans]])) predict(bn_list[[key_trans]], emmean, inverse = TRUE) else NA_real_,
      SE_bt = if (!is.null(bn_list[[key_trans]])) predict(bn_list[[key_trans]], SE, inverse = TRUE) else NA_real_,
      lower.CL_bt = if (!is.null(bn_list[[key_trans]])) predict(bn_list[[key_trans]], lower.CL, inverse = TRUE) else NA_real_,
      upper.CL_bt = if (!is.null(bn_list[[key_trans]])) predict(bn_list[[key_trans]], upper.CL, inverse = TRUE) else NA_real_
    ) %>%
    dplyr::select(any_of(c(
      "ghgspecies", "casepilot","comparison", "status", "vegpresence","season", "df",
      "emmean", "SE", "lower.CL", "upper.CL", "group_letter",
      "emmean_bt", "SE_bt", "lower.CL_bt", "upper.CL_bt"
    ))) %>%
    arrange(ghgspecies, casepilot, comparison, status, season)
  
  # Store casepilot results
  complex_comparisons_list[[dataset]] <- all_comparisons
}

#Remove within-loop objects
rm(all_comparisons, 
   status_comparisons,season_comparisons,statuswithinseason_comparisons,statuswithinvegpresence_comparisons,
   status_emmeans, season_emmeans, statuswithinseason_emmeans,  statuswithinvegpresence_emmeans)


# Unir todos los resultados en un solo data.frame
complexmodel_emmeans_withgroups <- bind_rows(complex_comparisons_list)

#Save all comparisons:  
write.csv(x = complexmodel_emmeans_withgroups, file = paste0(paper_path,"Emmeans-posthoc_proportionalweights_complexmodels_chambermodels.csv"),row.names = F)





##Opt2.(YES) Using custom weights------
#WE NEED TO APPLY weights to scale the vegpresence effects according to their proportion in the field. This applies for comparisions: status, season, and status_within_season. For status_within_vegpresence, weights shoudl not be applied, use "equal".

#First, calculate vegpresence proportions at every level of casepilot, status and season (but combine vegpresence proportion of each subsite). FOr the moment use field-data, update when we have RS-retrieved proportions

#For CURONIAN: Calculate overall veg/noveg of curonian, by design this proportion should be the same for all subsites, additionally, we asume constant proportion across seasons. we will use this to override the seasonal chamber distribution
cu_weights<- data4models %>% 
  filter(casepilot=="CU") %>% 
  group_by(plotcode, casepilot,vegpresence) %>% 
  distinct() %>% #to remove duplicates (same plotcode, different ghgspecies)
  group_by(casepilot, vegpresence) %>% 
  #calculate overall vegpresence counts (across all seasons and all subsites)
  summarise(n_vegpresence=sum(!is.na(vegpresence)), .groups = "drop") %>% 
  pivot_wider(names_from = vegpresence, values_from = n_vegpresence,values_fill = 0) %>% 
  #Calculate vegpresence proportion 
  mutate(sum_all=`Non-vegetated`+Vegetated) %>% 
  mutate(Vegetated=Vegetated/sum_all,
         `Non-vegetated`=`Non-vegetated`/sum_all) %>% 
  dplyr::select(-sum_all)


#For the moment, use field-data, update when we have RS-retrieved proportions
vegpresence4weights<- data4models %>% 
  dplyr::select(plotcode, casepilot, season, status, vegpresence, subsite) %>% 
  #remove duplicates (same plotcode, different ghgspecies)
  distinct() %>% 
  dplyr::group_by(casepilot, season, status, vegpresence, subsite) %>%
  #Calculate vegpresence deployment counts for each (subsite) sampling and pivot_wider
  summarise(n_vegpresence=sum(!is.na(vegpresence)), .groups = "drop") %>% 
  pivot_wider(names_from = vegpresence, values_from = n_vegpresence,values_fill = 0) %>% 
  #Calculate vegpresence proportion for every (subsite) sampling.
  mutate(sum_all=`Non-vegetated`+Vegetated) %>% 
  mutate(Vegetated=Vegetated/sum_all,
         `Non-vegetated`=`Non-vegetated`/sum_all) %>% 
  dplyr::select(-sum_all) %>% 
  pivot_longer(cols = c(`Non-vegetated`,Vegetated), names_to = "vegpresence", values_to = "proportions") %>% 
  #Calculate the average of each status (giving equal weight to the two subsites)
  dplyr::group_by(casepilot, status, season, vegpresence) %>% 
  summarise(prop_weight=mean(proportions, na.rm = T), .groups = "drop") %>% 
  #Override CU composition: constant actual composition, differences in chamber deployment are due to systematic sampling biass in Restored sites
  mutate(prop_weight=if_else(casepilot=="CU"&vegpresence=="Vegetated", cu_weights %>% pull(`Vegetated`),prop_weight),
         prop_weight=if_else(casepilot=="CU"&vegpresence=="Non-vegetated",cu_weights %>% pull(`Non-vegetated`),prop_weight))


#Calculate emmeans for each relevant comparison, using weights when appropriate:

#MANUAL loop:
  #step 1. Estimate emmeans for all combinations status season vegpresence (equal weights)
  #step 2. Join custom weights and calculated weighted means (and SE) per comparison group
  #step 3. Perform post-hoc Wald test in transformed scale
  #step 4. Obtain CLDs letters to join with emmeans
  #step 5. Save estimates and results of post-hoc tests 


#OMIT FOR THE MOMENT: This loop  estimates SE via conservative SE formula assuming independence of EMMs 
#We could improve SE accuracy using the vcov() from emmeans and the wieghed formula: 
# SE = sqrt(t(w) %*% V %*% w) 
#AND CONSIDERING: 
# w <- custom_weights_vector  # same length as number of rows in emm_df
# Make sure w sums to 1 within each group
# 5. Get covariance matrix of EMM estimates
# V <- vcov(emm)


#POST HOC function: 
# Estimated marginal means (EMMs) were calculated for each level of the fixed effects (status, season, vegpresence) using model predictions from linear mixed-effects models (lmer). To account for sampling bias due to minimum representation constraints, custom weights reflecting the observed distribution of vegpresence within each status and season level were applied manually to the EMMs. Standard errors were conservatively estimated by propagating uncertainty assuming independence between strata.
# Post-hoc pairwise comparisons between EMMs were performed manually using Wald tests. For each comparison, the difference in EMMs was divided by the standard error of the difference (computed as the square root of the sum of squared standard errors), and two-sided p-values were obtained from the standard normal distribution. 
# To account for multiple comparisons within each grouping (e.g., status within season), p-values were adjusted using the Sidak correction. Letters were assigned to significantly different groups using the multcompLetters funcion (using sidak-adjusted pvalues and a significance threshold of 0.05).


# Function to compute post-hoc Wald tests of Emmeans within loop (with sidak adjustment for multiple comparisons)
compute_pairwise_comparisons_grouped <- function(df, compare_var, group_vars=character(0)) {
  df %>%
    group_by(across(all_of(group_vars))) %>%
    group_split() %>%
    purrr::map_dfr(function(group_df) {
      combs <- combn(unique(group_df[[compare_var]]), 2, simplify = FALSE)
      n_comparisons <- length(combs)
      
      purrr::map_dfr(combs, function(pair) {
        group1 <- group_df %>% filter(!!sym(compare_var) == pair[1])
        group2 <- group_df %>% filter(!!sym(compare_var) == pair[2])
        
        diff <- group1$emmean - group2$emmean
        se_diff <- sqrt(group1$SE^2 + group2$SE^2)
        z <- diff / se_diff
        p_value <- 2 * pnorm(-abs(z))
        
        # Sidak correction
        p_value_sidak <- 1 - (1 - p_value)^n_comparisons
        
        tibble(
          group1 = pair[1],
          group2 = pair[2],
          diff = diff,
          se_diff = se_diff,
          z = z,
          p_value = p_value,
          p_value_sidak = p_value_sidak,
          comparison = unique(group_df$comparison),
          ghgspecies = unique(group_df$ghgspecies),
          casepilot = unique(group_df$casepilot),
          !!!group_df[1, group_vars]
        )
      })
    })
}

#Function to generate group letters Compact letter display based on the output of the previous function: 
generate_cld <- function(df, group_var, p_col = "p_value_sidak") {
  # Create a named vector of p-values
  comparisons <- paste(df$group1, df$group2, sep = "-")
  pvals <- setNames(df[[p_col]], comparisons)
  
  # Generate CLDs
  cld <- multcompLetters(pvals, threshold = 0.05)$Letters
  
  # Return as data frame
  tibble(
    !!group_var := names(cld),
    group_letter = cld
  )
}


# Initialize list for storing comparisons
complex_comparisons_list <- list()

# Loop to extract all comparisons for every casepilot best model
for (dataset in names(complexmodel_list_allghg)) {
  cp_model <- complexmodel_list_allghg[[dataset]]
  
  #extract list of fixed effects of model:
  fixed_effects_used<- attr(terms(cp_model), "term.labels")
  
  #extract casepilot name and ghgspecies from model-list names
  casepilot_name<- sub("_.*", "", dataset)
  ghgspecies<- sub(paste0(casepilot_name,"_"),"",dataset)
  
  #Filter external weights
  ext_weights<- vegpresence4weights %>% filter(casepilot==casepilot_name)
  
  # STEP 1. Estimate emmeans for all combinations
  
  emm_df <- emmeans(cp_model, ~ status * season * vegpresence) %>% 
    as.data.frame()
  
  # STEP 2. Merge custom weights 
  emm_df <- emm_df %>%
    left_join(ext_weights, by = c("status", "season", "vegpresence"))
  
  
  # STATUS — weighted by custom weights
  status_emmeans_weighted <- emm_df %>%
    group_by(status) %>%
    summarise(
      emmean = sum(emmean * prop_weight) / sum(prop_weight),
      SE = sqrt(sum((SE^2) * prop_weight^2)) / sum(prop_weight),
      comparison = "status",
      ghgspecies = ghgspecies,
      casepilot = casepilot_name,
      season = NA,
      vegpresence = NA,
      .groups = "drop"
    )
  
  # SEASON — weighted by custom weights
  season_emmeans_weighted <- emm_df %>%
    group_by(season) %>%
    summarise(
      emmean = sum(emmean * prop_weight) / sum(prop_weight),
      SE = sqrt(sum((SE^2) * prop_weight^2)) / sum(prop_weight),
      comparison = "season",
      ghgspecies = ghgspecies,
      casepilot = casepilot_name,
      status = NA,
      vegpresence = NA,
      .groups = "drop"
    )
  
  # STATUS within SEASON — custom weights per group
  status_within_season_emmeans <- emm_df %>%
    group_by(season, status) %>%
    summarise(
      emmean = sum(emmean * prop_weight) / sum(prop_weight),
      SE = sqrt(sum((SE^2) * prop_weight^2)) / sum(prop_weight),
      comparison = "status_within_season",
      ghgspecies = ghgspecies,
      casepilot = casepilot_name,
      vegpresence = NA,
      .groups = "drop"
    )
  
  # STATUS within VEGPRESENCE — equal weights only
  # For equal weighting, we just average over seasons (we can obtain the same same output via emmeans(~status|vegpresence, weights= "equal"), potentially with more accurate SE)
  status_within_vegpresence_emmeans <- emm_df %>%
    group_by(vegpresence, status) %>%
    summarise(
      emmean = mean(emmean),
      SE = sqrt(mean(SE^2)),
      comparison = "status_within_vegpresence",
      ghgspecies = ghgspecies,
      casepilot = casepilot_name,
      season = NA,
      .groups = "drop"
    )
  
  # STEP 3: Manual Wald post-hoc (test differences in emmeans for every emmeans data frame above) 
  status_comparisons <- compute_pairwise_comparisons_grouped(status_emmeans_weighted, "status")
  season_comparisons <- compute_pairwise_comparisons_grouped(season_emmeans_weighted, "season")
  status_within_season_comparisons <- compute_pairwise_comparisons_grouped(status_within_season_emmeans, "status","season")
  status_within_vegpresence_comparisons <- compute_pairwise_comparisons_grouped(status_within_vegpresence_emmeans, "status","vegpresence")
  
  
  #STEP4: Obtain Compact Letter Display for different groups, based on previous post-hoc
  cld_status<- status_comparisons %>% 
    group_by(ghgspecies, casepilot) %>% 
    group_split() %>% 
    purrr::map_dfr(function(df) {
      ghgspecies_val <- df$ghgspecies[1]
      casepilot_val <- df$casepilot[1]
      cld <- generate_cld(df, "group1") %>%
        rename(status = group1) %>%
        mutate(
          ghgspecies = ghgspecies_val,
          casepilot = casepilot_val
        )
      return(cld)
    })
  
  cld_season<- season_comparisons %>% 
    group_by(ghgspecies, casepilot) %>% 
    group_split() %>% 
    purrr::map_dfr(function(df) {
      # Extract grouping values from the first row
      ghgspecies_val <- df$ghgspecies[1]
      casepilot_val <- df$casepilot[1]
      
      # Generate CLD
      cld <- generate_cld(df, "group1") %>%
        rename(season = group1) %>%
        mutate(
          ghgspecies = ghgspecies_val,
          casepilot = casepilot_val
        )
      return(cld)
    })
  

  cld_status_within_season <- status_within_season_comparisons %>%
    group_by(season, ghgspecies, casepilot) %>%
    group_split() %>%
    purrr::map_dfr(function(df) {
      # Extract grouping values from the first row
      season_val <- df$season[1]
      ghgspecies_val <- df$ghgspecies[1]
      casepilot_val <- df$casepilot[1]
      
      # Generate CLD
      cld <- generate_cld(df, "group1") %>%
        rename(status = group1) %>%
        mutate(
          season = season_val,
          ghgspecies = ghgspecies_val,
          casepilot = casepilot_val
        )
      return(cld)
    })
  
  cld_status_within_vegpresence <- status_within_vegpresence_comparisons %>%
    group_by(vegpresence, ghgspecies, casepilot) %>%
    group_split() %>%
    purrr::map_dfr(function(df) {
      # Extract grouping values from the first row
      vegpresence_val <- df$vegpresence[1]
      ghgspecies_val <- df$ghgspecies[1]
      casepilot_val <- df$casepilot[1]
      
      # Generate CLD
      cld <- generate_cld(df, "group1") %>%
        rename(status = group1) %>%
        mutate(
          vegpresence = vegpresence_val,
          ghgspecies = ghgspecies_val,
          casepilot = casepilot_val
        )
      return(cld)
    })
  
  
  
  # STEP 4: Combine outputs into 2 data frames.
  
  #Pairwise comparisons: 
  all_posthoc_comparisons <- bind_rows(
    status_comparisons,
    season_comparisons,
    status_within_season_comparisons,
    status_within_vegpresence_comparisons
  )
  
  # Merge CLD with emmeans
  status_emmeans_weighted<- status_emmeans_weighted %>% 
    left_join(cld_status, by=c("status","ghgspecies","casepilot"))
  
  season_emmeans_weighted<- season_emmeans_weighted %>% 
    left_join(cld_season, by=c("season","ghgspecies","casepilot"))
  
  status_within_season_emmeans <- status_within_season_emmeans %>%
    left_join(cld_status_within_season, by = c("status", "season", "ghgspecies", "casepilot"))
  
  status_within_vegpresence_emmeans <- status_within_vegpresence_emmeans %>%
    left_join(cld_status_within_vegpresence, by = c("status", "vegpresence", "ghgspecies", "casepilot"))
  
  #Join all emmeans: 
  all_emmeans_cld<- bind_rows(
    status_emmeans_weighted,
    season_emmeans_weighted,
    status_within_season_emmeans,
    status_within_vegpresence_emmeans
  )

  
# Store casepilot results as named list of 2 lists
  complex_comparisons_list[[dataset]] <- list(
    emmeans_with_cld = all_emmeans_cld,
    posthoc_comparisons = all_posthoc_comparisons
  )
  
}


# Extract all emmeans with CLDs
complexmodel_emmeans_withgroups<- purrr::map_dfr(complex_comparisons_list, "emmeans_with_cld") %>% 
  #Back-transform to molar units
  mutate(key_trans=paste0(casepilot, "_",ghgspecies)) %>% 
  rowwise() %>% 
  mutate(
    emmean_bt = if (!is.null(bn_list[[key_trans]])) predict(bn_list[[key_trans]], emmean, inverse = TRUE) else NA_real_,
    SE_bt = if (!is.null(bn_list[[key_trans]])) predict(bn_list[[key_trans]], SE, inverse = TRUE) else NA_real_
  )

#Save emmeans and groupletters (back-transformed):  
write.csv(x = complexmodel_emmeans_withgroups, file = paste0(paper_path,"Emmeans-groupletters_complexmodels_chambermodels.csv"),row.names = F)


# Extract all post-hoc comparisons (keep in model scale)
complexmodel_posthoc_tests <- purrr::map_dfr(complex_comparisons_list, "posthoc_comparisons")

#Save post-hoc comparisons (to go to supplementary table)
write.csv(x = complexmodel_posthoc_tests, file = paste0(paper_path,"Emmeans-posthoctests_complexmodels_chambermodels.csv"),row.names = F)









#UNTIL HERE----------
#Beyond here, the script contains old approaches that have been deemed not appropriate. 



##Opt0. (NO) Use equal weights------


# Initialize list for storing comparisons
comparisons_list <- list()

# Loop to extract all comparisons for every casepilot best model
for (dataset in names(specific_best_model_list_allghg)) {
  cp_model <- specific_best_model_list_allghg[[dataset]]
  
  #extract list of fixed effects of model:
  fixed_effects_used<- attr(terms(cp_model), "term.labels")
  
  #extract casepilot name and ghgspecies from model-list names
  casepilot_name<- sub("_.*", "", dataset)
  ghgspecies<- sub(paste0(casepilot_name,"_"),"",dataset)
  
  # Status comparisons: always
  status_emmeans <- emmeans(cp_model, ~status) #ADD weights="proportional" to use the chamber distribution as a weight to calcualte the overall emmean for status. It will also work for season (if we have differing number of observations each season)
  status_comparisons <- cld(status_emmeans, 
                            alpha = 0.05,
                            Letters = letters,
                            adjust = "sidak") %>%
    mutate(comparison = "status",
           ghgspecies = ghgspecies,
           casepilot = casepilot_name,
           season = NA) %>%
    rename(lower.CL = any_of("asymp.LCL"),
           upper.CL = any_of("asymp.UCL"))
  
  # Season comparisons: 
  season_emmeans <- emmeans(cp_model, ~season)
  season_comparisons <- cld(season_emmeans, 
                            alpha = 0.05,
                            Letters = letters,
                            adjust = "sidak") %>%
    mutate(comparison = "season",
           ghgspecies = ghgspecies,
           casepilot = casepilot_name,
           status = NA) %>%
    rename(lower.CL = any_of("asymp.LCL"),
           upper.CL = any_of("asymp.UCL"))
  
  # Interaction comparisons: 
  interaction_emmeans <- emmeans(cp_model, ~status|season)
  interaction_comparisons <- cld(interaction_emmeans, 
                                 alpha = 0.05,
                                 Letters = letters,
                                 adjust = "sidak") %>%
    mutate(comparison = "status_within_season",
           ghgspecies = ghgspecies,
           casepilot = casepilot_name) %>%
    rename(lower.CL = any_of("asymp.LCL"),
           upper.CL = any_of("asymp.UCL"))
  
  if("strata"%in%fixed_effects_used){
    status_within_strata_emmeans<-  emmeans(cp_model, ~status|strata)
    status_within_strata_comparisons<- cld(status_within_strata_emmeans,
                                           alpha=0.05,
                                           Letters=letters,
                                           adjust = "sidak") %>% 
      mutate(comparison="status_within_strata",
             ghgspecies=ghgspecies,
             casepilot=casepilot_name) %>% 
      rename(lower.CL = any_of("asymp.LCL"),
             upper.CL = any_of("asymp.UCL"))
  } else{ status_within_strata_comparisons<-data.frame(emmean=numeric(0))}
  
  # Join and format CP comparisons
  all_comparisons <- status_comparisons %>%
    full_join(season_comparisons) %>%
    full_join(interaction_comparisons) %>%
    full_join(status_within_strata_comparisons) %>% 
    rename(group_letter = .group) %>%
    mutate(group_letter = gsub(" ", "", group_letter),
           
           # Create key for accessing BestNormalize object
           key_trans = paste(casepilot_name, ghgspecies, sep = "_")
    ) %>%
    rowwise() %>% 
    # BACK-TRANSFORMATION of emmean, SE and CL using BestNormalize object
    mutate(
      emmean_bt = if (!is.null(bn_list[[key_trans]])) predict(bn_list[[key_trans]], emmean, inverse = TRUE) else NA_real_,
      SE_bt = if (!is.null(bn_list[[key_trans]])) predict(bn_list[[key_trans]], SE, inverse = TRUE) else NA_real_,
      lower.CL_bt = if (!is.null(bn_list[[key_trans]])) predict(bn_list[[key_trans]], lower.CL, inverse = TRUE) else NA_real_,
      upper.CL_bt = if (!is.null(bn_list[[key_trans]])) predict(bn_list[[key_trans]], upper.CL, inverse = TRUE) else NA_real_
    ) %>%
    
    dplyr::select(any_of(c(
      "ghgspecies", "casepilot","comparison", "status", "season","strata", "df",
      "emmean", "SE", "lower.CL", "upper.CL", "group_letter",
      "emmean_bt", "SE_bt", "lower.CL_bt", "upper.CL_bt"
    ))) %>%
    arrange(ghgspecies, casepilot, comparison, status, season)
  
  # Store casepilot results
  comparisons_list[[dataset]] <- all_comparisons
}

#Remove within-loop objects
rm(all_comparisons, status_comparisons,season_comparisons,interaction_comparisons,
   interaction_emmeans,status_emmeans, season_emmeans)


# Unir todos los resultados en un solo data.frame
specific_emmeans_withgroups <- bind_rows(comparisons_list)

#Save all comparisons:  
write.csv(x = specific_emmeans_withgroups, file = paste0(paper_path,"Allghg_emmeans-posthoc_specificchambermodels.csv"),row.names = F)



# Initialize list for storing comparisons
comparisons_list <- list()

# Loop to extract all comparisons for every casepilot best model
for (dataset in names(best_model_list_allghg)) {
  cp_model <- best_model_list_allghg[[dataset]]
  
  #extract casepilot name and ghgspecies from model-list names
  casepilot_name<- sub("_.*", "", dataset)
  ghgspecies<- sub(paste0(casepilot_name,"_"),"",dataset)
  
  # Status comparisons:
  status_emmeans <- emmeans(cp_model, ~status, weights = "proportional")
  status_comparisons <- cld(status_emmeans, 
                            alpha = 0.05,
                            Letters = letters,
                            adjust = "sidak") %>%
    mutate(comparison = "status",
           ghgspecies = ghgspecies,
           casepilot = casepilot_name,
           season = NA) %>%
    rename(lower.CL = any_of("asymp.LCL"),
           upper.CL = any_of("asymp.UCL"))
  
  # Season comparisons: 
  season_emmeans <- emmeans(cp_model, ~season, weights = "proportional")
  season_comparisons <- cld(season_emmeans, 
                            alpha = 0.05,
                            Letters = letters,
                            adjust = "sidak") %>%
    mutate(comparison = "season",
           ghgspecies = ghgspecies,
           casepilot = casepilot_name,
           status = NA) %>%
    rename(lower.CL = any_of("asymp.LCL"),
           upper.CL = any_of("asymp.UCL"))
  
  # Interaction comparisons: 
  interaction_emmeans <- emmeans(cp_model, ~status|season)
  interaction_comparisons <- cld(interaction_emmeans, 
                                 alpha = 0.05,
                                 Letters = letters,
                                 adjust = "sidak") %>%
    mutate(comparison = "status_within_season",
           ghgspecies = ghgspecies,
           casepilot = casepilot_name) %>%
    rename(lower.CL = any_of("asymp.LCL"),
           upper.CL = any_of("asymp.UCL"))
  
  # Join and format CP comparisons
  all_comparisons <- status_comparisons %>%
    full_join(season_comparisons) %>%
    full_join(interaction_comparisons) %>%
    rename(group_letter = .group) %>%
    mutate(
      comparison = factor(comparison, levels = c("status", "season", "status_within_season")),
      season = factor(season, levels = c("S1", "S2", "S3", "S4")),
      status = factor(status, levels = c("Altered", "Preserved", "Restored")),
      group_letter = gsub(" ", "", group_letter),
      
      # Create key for accessing BestNormalize object
      key_trans = paste(casepilot_name, ghgspecies, sep = "_")
    ) %>%
    rowwise() %>% 
    # BACK-TRANSFORMATION of emmean, SE and CL using BestNormalize object
    mutate(
      emmean_bt = if (!is.null(bn_list[[key_trans]])) predict(bn_list[[key_trans]], emmean, inverse = TRUE) else NA_real_,
      SE_bt = if (!is.null(bn_list[[key_trans]])) predict(bn_list[[key_trans]], SE, inverse = TRUE) else NA_real_,
      lower.CL_bt = if (!is.null(bn_list[[key_trans]])) predict(bn_list[[key_trans]], lower.CL, inverse = TRUE) else NA_real_,
      upper.CL_bt = if (!is.null(bn_list[[key_trans]])) predict(bn_list[[key_trans]], upper.CL, inverse = TRUE) else NA_real_
    ) %>%
    
    dplyr::select(any_of(c(
      "ghgspecies", "casepilot","comparison", "status", "season", "df",
      "emmean", "SE", "lower.CL", "upper.CL", "group_letter",
      "emmean_bt", "SE_bt", "lower.CL_bt", "upper.CL_bt"
    ))) %>%
    arrange(ghgspecies, casepilot, comparison, status, season)
  
  # Store casepilot results
  comparisons_list[[dataset]] <- all_comparisons
}

#Remove within-loop objects
rm(all_comparisons, status_comparisons,season_comparisons,interaction_comparisons,
   interaction_emmeans,status_emmeans, season_emmeans)


# Unir todos los resultados en un solo data.frame
all_emmeans_withgroups <- bind_rows(comparisons_list)

#Save all comparisons:  
write.csv(x = all_emmeans_withgroups, file = paste0(paper_path,"Allghg_emmeans-posthoc_chambermodels.csv"),row.names = F)


##Pairwise contrasts-----

# Initialize list for storing comparisons
comparisons_pairwise_list <- list()

# Loop to extract all comparisons for every casepilot best model
for (dataset in names(new_best_model_list_allghg)) {
  cp_model <- new_best_model_list_allghg[[dataset]]
  
  # Extract casepilot name and ghgspecies from model-list names
  casepilot_name <- sub("_.*", "", dataset)
  ghgspecies <- sub(paste0(casepilot_name, "_"), "", dataset)
  
  # Status pairwise comparisons
  status_emmeans <- emmeans(cp_model, ~status)
  status_pairs <- pairs(status_emmeans, adjust = "sidak") %>%
    as.data.frame() %>%
    mutate(comparison = "status",
           ghgspecies = ghgspecies,
           casepilot = casepilot_name)
  
  # Season pairwise comparisons
  season_emmeans <- emmeans(cp_model, ~season)
  season_pairs <- pairs(season_emmeans, adjust = "sidak") %>%
    as.data.frame() %>%
    mutate(comparison = "season",
           ghgspecies = ghgspecies,
           casepilot = casepilot_name)
  
  # Interaction pairwise comparisons (status within season)
  interaction_emmeans <- emmeans(cp_model, ~status | season)
  interaction_pairs <- pairs(interaction_emmeans, adjust = "sidak") %>%
    as.data.frame() %>%
    mutate(comparison = "status_within_season",
           ghgspecies = ghgspecies,
           casepilot = casepilot_name)
  
  # Combine all comparisons
  all_comparisons <- bind_rows(status_pairs, season_pairs, interaction_pairs) %>%
    mutate(
      comparison = factor(comparison, levels = c("status", "season", "status_within_season")),
      season = factor(season, levels = c("S1", "S2", "S3", "S4")),
      database = paste(casepilot_name, ghgspecies, sep = "_")
    ) 
  
  #Skip back-transformation, we are calculating for each comparison the difference in emmeans the SE of that difference and its significance, I dont think it should be back-transformed  
  # Store results
  comparisons_pairwise_list[[dataset]] <- all_comparisons
}

# Combine all results
all_pairwise_comparisons <- bind_rows(comparisons_pairwise_list) %>% 
  mutate(model_distribution=if_else(df==Inf, "t_family", "gaussian"),
         test_used=if_else(df==Inf, "Z.tests","T.tests"),
         stat.ratio=if_else(df==Inf, z.ratio, t.ratio)) %>% 
  dplyr::select(casepilot, ghgspecies, database, model_distribution,comparison, test_used, contrast, season, estimate, SE, df, stat.ratio, p.value)


#Save all comparisons:  
write.csv(x = all_pairwise_comparisons, file = paste0(paper_path,"Allghg_pairwise_diffemmeans-posthoc_chambermodels.csv"),row.names = F)




#_______-------


#CP- tailored models---------

#Here we will answer why there is or not an effect of status. 
#Each casepilot has different data-capabilities to show this, which will require specific models with different structure to best answer the question in each case.

#Main approaches: 

#TO-ADAPT: RI models ------
#Build veg/noveg model for RI (without status): dailyflux_trans~vegpresence*season +(1|subsite). This will be used to obtain the average flux for vegetated and not vegetated, as we cannot test for different strata-specific rates related to status. We could use these estimates Jointly with the strata coverage to answer the same question as in the other casepilots, plus we can give a reduction of CO2, CH4 GWP fluxes directly linked to zostera (and compare with literature). Doubt: remove tidal pools beforehand, remove non-zostera plants?



#OLD: 
#RI: only question that we can answer is whether the vegetated strata from restored areas and well-preserved areas behave the same. Filter only vegetated strata (re-transform) and repeat the modelling approach for the general model (status will only have 2 levels).

#This question is not important and it requires to re-transform the data (too complex to export results with the rest)


##specific_RI_co2-------
ri_vegonly_co2<- data4models %>% filter(strata=="vegetated",status%in%c("Preserved", "Restored"),casepilot=="RI", ghgspecies=="co2")

#specific data-tranformation
ri_vegonly_co2_bn_result<- bestNormalize(x = ri_vegonly_co2$dailyflux, new_transforms = pseudolog_transform, standardize = TRUE, warn = TRUE, k = 5,allow_orderNorm = F)

ri_vegonly_co2<-ri_vegonly_co2 %>% 
  mutate(dailyflux_trans=predict(ri_vegonly_co2_bn_result, dailyflux))

#See effect of transformation:
  plot_bestNormalize(ri_vegonly_co2_bn_result)


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = ri_vegonly_co2,
                           family = gaussian(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m1_gaus_nostrata) #OK
check_singularity(m1_gaus_nostrata,tolerance = 1e-8) #OK
r2(m1_gaus_nostrata) 
Anova(m1_gaus_nostrata)
summary(m1_gaus_nostrata)
res<- simulateResiduals(m1_gaus_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_nostrata)
em<- emmeans(m1_gaus_nostrata, ~ status)
summary(em)
pairs(em)
#DECISSION: GOOD MODEL


#SKIP T-family models, already gaussian works OK

# t_family options
#Most basic: t_family, only status and season 
m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = ri_vegonly_co2,
                        family = t_family(),
                        dispformula = ~1)
#Evaluate basic model:
check_convergence(m3_t_nostrata)
check_singularity(m3_t_nostrata,tolerance = 1e-8)
r2(m3_t_nostrata, tolerance=1e-10) 
Anova(m3_t_nostrata)
summary(m3_t_nostrata)
res<- simulateResiduals(m3_t_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m3_t_nostrata)
em<- emmeans(m3_t_nostrata, ~ status)
summary(em)
pairs(em)




#Set as best
specific_ri_best_co2<- m1_gaus_nostrata 


#check the emmean comparison graphically:
em<- emmeans(specific_ri_best_co2, ~ status)
ggplot(as.data.frame(em), aes(x = status, y = emmean, col=status)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  labs(
    title = "Estimated Marginal Means by Status for CO2",
    x = "Status",
    y = "EMM (model scale)"
  ) +
  theme_minimal()

#Plot obs vs pred 
data.frame(obs=specific_ri_best_co2$frame$dailyflux_trans,
           pred=predict(specific_ri_best_co2, type = "response"),
           strata=factor(ri_vegonly_co2$strata, levels = c("bare","vegetated", "open water"))) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=strata),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  stat_cor(aes(col=strata),label.x.npc = 0.4, label.y.npc = 0.2)+
  geom_smooth(method = "lm", aes(col=strata), se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  scale_color_manual(values = c(
    "bare" = "#d95f02",      # brown/orange
    "vegetated" = "#2db82d", # green
    "open water" = "blue" # blue
  ), drop = FALSE) +
  labs(title = "specific_RI_co2 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)",
       caption= "This model only contains vegetated strata from preserved and restored sites")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())


ggsave(filename = "specific_RI_co2_Observed_VS_predicted.png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)


##specific_RI_ch4-------
ri_vegonly_ch4<- data4models %>% filter(strata=="vegetated",status%in%c("Preserved", "Restored"),casepilot=="RI", ghgspecies=="ch4")

#specific data-tranformation
ri_vegonly_ch4_bn_result<- bestNormalize(x = ri_vegonly_ch4$dailyflux, new_transforms = pseudolog_transform, standardize = TRUE, warn = TRUE, k = 5,allow_orderNorm = F)

ri_vegonly_ch4<-ri_vegonly_ch4 %>% 
  mutate(dailyflux_trans=predict(ri_vegonly_ch4_bn_result, dailyflux))

#See effect of transformation:
plot_bestNormalize(ri_vegonly_ch4_bn_result)


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = ri_vegonly_ch4,
                           family = gaussian(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m1_gaus_nostrata) #OK
check_singularity(m1_gaus_nostrata,tolerance = 1e-8) #OK
r2(m1_gaus_nostrata) 
Anova(m1_gaus_nostrata)
summary(m1_gaus_nostrata)
res<- simulateResiduals(m1_gaus_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_nostrata)
em<- emmeans(m1_gaus_nostrata, ~ status)
summary(em)
pairs(em)
#DECISSION: RESIUDALS FAIL




# t_family options
#Most basic: t_family, only status and season 
m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = ri_vegonly_ch4,
                        family = t_family(),
                        dispformula = ~1)
#Evaluate basic model:
check_convergence(m3_t_nostrata)
check_singularity(m3_t_nostrata,tolerance = 1e-8)
r2(m3_t_nostrata, tolerance=1e-10) 
Anova(m3_t_nostrata)
summary(m3_t_nostrata)
res<- simulateResiduals(m3_t_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m3_t_nostrata)
em<- emmeans(m3_t_nostrata, ~ status)
summary(em)
pairs(em)




#Set as best
specific_ri_best_ch4<- m3_t_nostrata 


#check the emmean comparison graphically:
em<- emmeans(specific_ri_best_ch4, ~ status)
ggplot(as.data.frame(em), aes(x = status, y = emmean, col=status)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
  labs(
    title = "Estimated Marginal Means by Status for CH4",
    x = "Status",
    y = "EMM (model scale)"
  ) +
  theme_minimal()


#Plot obs vs pred 
data.frame(obs=specific_ri_best_ch4$frame$dailyflux_trans,
           pred=predict(specific_ri_best_ch4, type = "response"),
           strata=factor(ri_vegonly_ch4$strata, levels = c("bare","vegetated", "open water"))) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=strata),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  stat_cor(aes(col=strata),label.x.npc = 0.4, label.y.npc = 0.2)+
  geom_smooth(method = "lm", aes(col=strata), se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  scale_color_manual(values = c(
    "bare" = "#d95f02",      # brown/orange
    "vegetated" = "#2db82d", # green
    "open water" = "blue" # blue
  ), drop = FALSE) +
  labs(title = "specific_RI_ch4 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)",
       caption= "This model only contains vegetated strata from preserved and restored sites")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())


ggsave(filename = "specific_RI_ch4_Observed_VS_predicted.png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)


##specific_RI_GWPco2andch4-------
ri_vegonly_GWPco2andch4<- data4models %>% filter(strata=="vegetated",status%in%c("Preserved", "Restored"),casepilot=="RI", ghgspecies=="GWPco2andch4")

#specific data-tranformation
ri_vegonly_GWPco2andch4_bn_result<- bestNormalize(x = ri_vegonly_GWPco2andch4$dailyflux, new_transforms = pseudolog_transform, standardize = TRUE, warn = TRUE, k = 5,allow_orderNorm = F)

ri_vegonly_GWPco2andch4<-ri_vegonly_GWPco2andch4 %>% 
  mutate(dailyflux_trans=predict(ri_vegonly_GWPco2andch4_bn_result, dailyflux))

#See effect of transformation:
plot_bestNormalize(ri_vegonly_GWPco2andch4_bn_result)


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = ri_vegonly_GWPco2andch4,
                           family = gaussian(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m1_gaus_nostrata) #OK
check_singularity(m1_gaus_nostrata,tolerance = 1e-8) #OK
r2(m1_gaus_nostrata) 
Anova(m1_gaus_nostrata)
summary(m1_gaus_nostrata)
res<- simulateResiduals(m1_gaus_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_nostrata)
em<- emmeans(m1_gaus_nostrata, ~ status)
summary(em)
pairs(em)
#DECISSION: GOOD MODEL



#SKYP, gaussian works ok
# t_family options
#Most basic: t_family, only status and season 
m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = ri_vegonly_GWPco2andch4,
                        family = t_family(),
                        dispformula = ~1)
#Evaluate basic model:
check_convergence(m3_t_nostrata)
check_singularity(m3_t_nostrata,tolerance = 1e-8)
r2(m3_t_nostrata, tolerance=1e-10) 
Anova(m3_t_nostrata)
summary(m3_t_nostrata)
res<- simulateResiduals(m3_t_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m3_t_nostrata)
em<- emmeans(m3_t_nostrata, ~ status)
summary(em)
pairs(em)




#Set as best
specific_ri_best_GWPco2andch4<- m1_gaus_nostrata 



#check the emmean comparison graphically:
em<- emmeans(specific_ri_best_GWPco2andch4, ~ status)
ggplot(as.data.frame(em), aes(x = status, y = emmean, col=status)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  labs(
    title = "Estimated Marginal Means by Status for GWPco2andch4",
    x = "Status",
    y = "EMM (model scale)"
  ) +
  theme_minimal()

#Plot obs vs pred 
data.frame(obs=specific_ri_best_GWPco2andch4$frame$dailyflux_trans,
           pred=predict(specific_ri_best_GWPco2andch4, type = "response"),
           strata=factor(ri_vegonly_GWPco2andch4$strata, levels = c("bare","vegetated", "open water"))) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=strata),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  stat_cor(aes(col=strata),label.x.npc = 0.4, label.y.npc = 0.2)+
  geom_smooth(method = "lm", aes(col=strata), se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  scale_color_manual(values = c(
    "bare" = "#d95f02",      # brown/orange
    "vegetated" = "#2db82d", # green
    "open water" = "blue" # blue
  ), drop = FALSE) +
  labs(title = "specific_RI_GWPco2andch4 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)",
       caption= "This model only contains vegetated strata from preserved and restored sites")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())


ggsave(filename = "specific_RI_GWPco2andch4_Observed_VS_predicted.png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)




 


#CA models 2strata (ok)---------
#To decide: (CA?), build 3-way models in which we can see the effects of restoration strata-by-strata.
#We saw that strata is not trongly linked to status, we expect very differnet fluxes depending on strata. Try 3-way models. 

#Issues arose when using 3-level strata (veg/ow/bare) due to imbalanced strata (missing water in some combinations).  CANNOT CALCULATE RESIDUALS OR EMMEANS


#Potentially, reclassify strata into two separate T/F variables: veg/noveg and water/nowater, using veg/noveg as fixed factor.
#(OMIT) and water/nowater as random nested (to account for different baseline flux between bare and water). 



##specific_CA_co2-----

ca_co2_vegnoveg<- data4models %>% filter(casepilot=="CA", ghgspecies=="co2") %>% 
  mutate(strata=as.factor(if_else(strata=="vegetated","vegetated", "non-vegetated")))

#3-way model using gaussian family
m1_gaus_3way<- glmmTMB(formula = dailyflux_trans~status*season*strata + (1|subsite), 
                       data = ca_co2_vegnoveg,
                       family = gaussian(),
                       dispformula = ~1)


#Evaluate basic model:
check_convergence(m1_gaus_3way) #OK
check_singularity(m1_gaus_3way,tolerance = 1e-8) #OK
r2(m1_gaus_3way) 
Anova(m1_gaus_3way)
summary(m1_gaus_3way)
res<- simulateResiduals(m1_gaus_3way)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_3way)
em<- emmeans(m1_gaus_3way, ~ status|strata)
summary(em)
pairs(em)
#DECISSION: Residuals fail

#TRY t_family option
m2_t_3way<- glmmTMB(formula = dailyflux_trans~status*season*strata + (1|subsite), 
                    data = ca_co2_vegnoveg,
                    family = t_family(),
                    dispformula = ~1)

#Evaluate basic model:
check_convergence(m2_t_3way) #OK
check_singularity(m2_t_3way,tolerance = 1e-8) #Singular, ok non-issue
r2(m2_t_3way, tolerance=1e-11) 
Anova(m2_t_3way)
summary(m2_t_3way)
res<- simulateResiduals(m2_t_3way)
plotQQunif(res)
plotResiduals(res)
check_residuals(m2_t_3way)
em<- emmeans(m2_t_3way, ~ status|strataveg)
summary(em)
pairs(em)
#DECISSION: Better, residuals have issues bust mostly are OK. 

#Compare models 
anova(m1_gaus_3way,m2_t_3way) #gaus vs t_family

#Set as best
specific_ca_best_co2<- m2_t_3way

#check the emmean comparison graphically:
em<- emmeans(specific_ca_best_co2, ~ status|strataveg)
ggplot(as.data.frame(em), aes(x = status, y = emmean, col=strataveg)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
  labs(
    title = "Estimated Marginal Means by Status for CO2",
    x = "Status",
    y = "EMM (model scale)"
  ) +
  theme_minimal()



#Plot obs vs pred 
data.frame(obs=specific_ca_best_co2$frame$dailyflux_trans,
           pred=predict(specific_ca_best_co2, type = "response"),
           strataveg=as.factor(ca_co2_vegnoveg$strata)) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=strataveg),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  stat_cor(aes(col=strataveg),label.x.npc = 0.4, label.y.npc = 0.2)+
  geom_smooth(method = "lm", aes(col=strataveg), se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "specific_ca_co2 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)",
       caption = "This model relevels strata into veg/no-veg and takes it as fixed effect (season*status*strataveg)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())


ggsave(filename = "specific_ca_co2_Observed_VS_predicted.png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)



##specific_CA_ch4-----

ca_ch4_vegnoveg<- data4models %>% filter(casepilot=="CA", ghgspecies=="ch4") %>% 
  mutate(strata=as.factor(if_else(strata=="vegetated","vegetated", "non-vegetated")))

#3-way model using gaussian family
m1_gaus_3way<- glmmTMB(formula = dailyflux_trans~status*season*strata + (1|subsite), 
                       data = ca_ch4_vegnoveg,
                       family = gaussian(),
                       dispformula = ~1)


#Evaluate basic model:
check_convergence(m1_gaus_3way) #OK
check_singularity(m1_gaus_3way,tolerance = 1e-8) #OK
r2(m1_gaus_3way) 
Anova(m1_gaus_3way)
summary(m1_gaus_3way)
res<- simulateResiduals(m1_gaus_3way)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_3way)
em<- emmeans(m1_gaus_3way, ~ status|strata)
summary(em)
pairs(em)
#DECISSION: GOOD RESIDUALS


#SKYP, gaussian is good.
#TRY t_family option
m2_t_3way<- glmmTMB(formula = dailyflux_trans~status*season*strata + (1|subsite), 
                    data = ca_ch4_vegnoveg,
                    family = t_family(),
                    dispformula = ~1)

#Evaluate basic model:
check_convergence(m2_t_3way) #OK
check_singularity(m2_t_3way,tolerance = 1e-8) #Singular, ok non-issue
r2(m2_t_3way, tolerance=1e-11) 
Anova(m2_t_3way)
summary(m2_t_3way)
res<- simulateResiduals(m2_t_3way)
plotQQunif(res)
plotResiduals(res)
check_residuals(m2_t_3way)
em<- emmeans(m2_t_3way, ~ status|strata)
summary(em)
pairs(em)
#DECISSION: Better, residuals have issues bust mostly are OK. 

#Compare models 
anova(m1_gaus_3way,m2_t_3way) #gaus vs t_family

#Set as best
specific_ca_best_ch4<- m1_gaus_3way

#check the emmean comparison graphically:
em<- emmeans(specific_ca_best_ch4, ~ status|strata)
ggplot(as.data.frame(em), aes(x = status, y = emmean, col=strata)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  labs(
    title = "Estimated Marginal Means by Status for CH4",
    x = "Status",
    y = "EMM (model scale)"
  ) +
  theme_minimal()



#Plot obs vs pred 
data.frame(obs=specific_ca_best_ch4$frame$dailyflux_trans,
           pred=predict(specific_ca_best_ch4, type = "response"),
           strataveg=as.factor(ca_ch4_vegnoveg$strata)) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=strataveg),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  stat_cor(aes(col=strataveg),label.x.npc = 0.4, label.y.npc = 0.2)+
  geom_smooth(method = "lm", aes(col=strataveg), se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "specific_ca_ch4 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)",
       caption = "This model relevels strata into veg/no-veg and takes it as fixed effect (season*status*strataveg)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())


ggsave(filename = "specific_ca_ch4_Observed_VS_predicted.png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)



##specific_CA_GWPco2andch4-----

ca_GWPco2andch4_vegnoveg<- data4models %>% filter(casepilot=="CA", ghgspecies=="GWPco2andch4") %>% 
  mutate(strata=as.factor(if_else(strata=="vegetated","vegetated", "non-vegetated")))

#3-way model using gaussian family
m1_gaus_3way<- glmmTMB(formula = dailyflux_trans~status*season*strata + (1|subsite), 
                       data = ca_GWPco2andch4_vegnoveg,
                       family = gaussian(),
                       dispformula = ~1)


#Evaluate basic model:
check_convergence(m1_gaus_3way) #OK
check_singularity(m1_gaus_3way,tolerance = 1e-8) #OK
r2(m1_gaus_3way) 
Anova(m1_gaus_3way)
summary(m1_gaus_3way)
res<- simulateResiduals(m1_gaus_3way)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_3way)
em<- emmeans(m1_gaus_3way, ~ status|strata)
summary(em)
pairs(em)
#DECISSION: BAD RESIDUALS (almost good)

#TRY t_family option
m2_t_3way<- glmmTMB(formula = dailyflux_trans~status*season*strata + (1|subsite), 
                    data = ca_GWPco2andch4_vegnoveg,
                    family = t_family(),
                    dispformula = ~1)

#Evaluate basic model:
check_convergence(m2_t_3way) #OK
check_singularity(m2_t_3way,tolerance = 1e-8) #Singular, ok non-issue
r2(m2_t_3way, tolerance=1e-11) 
Anova(m2_t_3way)
summary(m2_t_3way)
res<- simulateResiduals(m2_t_3way)
plotQQunif(res)
plotResiduals(res)
check_residuals(m2_t_3way)
em<- emmeans(m2_t_3way, ~ status|strata)
summary(em)
pairs(em)
#DECISSION: Better, residuals have issues bust mostly are OK. 

#Compare models 
anova(m1_gaus_3way,m2_t_3way) #gaus vs t_family

#Set as best
specific_ca_best_GWPco2andch4<- m2_t_3way

#check the emmean comparison graphically:
em<- emmeans(specific_ca_best_GWPco2andch4, ~ status|strata)
ggplot(as.data.frame(em), aes(x = status, y = emmean, col=strata)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
  labs(
    title = "Estimated Marginal Means by Status for GWP",
    x = "Status",
    y = "EMM (model scale)"
  ) +
  theme_minimal()



#Plot obs vs pred 
data.frame(obs=specific_ca_best_GWPco2andch4$frame$dailyflux_trans,
           pred=predict(specific_ca_best_GWPco2andch4, type = "response"),
           strataveg=as.factor(ca_GWPco2andch4_vegnoveg$strata)) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=strataveg),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  stat_cor(aes(col=strataveg),label.x.npc = 0.4, label.y.npc = 0.2)+
  geom_smooth(method = "lm", aes(col=strataveg), se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "specific_ca_GWPco2andch4 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)",
       caption = "This model relevels strata into veg/no-veg and takes it as fixed effect (season*status*strataveg)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())


ggsave(filename = "specific_ca_GWPco2andch4_Observed_VS_predicted.png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)








#VA models (DOUBT)---------
#We saw that strata is not trongly linked to status, we expect very differnet fluxes depending on strata. Try 3-way models. 

#Issues may arise from unbalanced strata (missing open water-preserved-S4). Model works but cannot calculate emmean for water-preserved


#We could re-level strata to join bare and water. This would allow for year-round calcualtion of effect of restoration depending on vegetation presence/absence

#Potentially, reclassify strata into two separate T/F variables: veg/noveg and water/nowater, using veg/noveg as fixed factor and water/nowater as random nested (to account for different baseline flux between bare and water). 



##specific_VA_co2-----

va_co2<- data4models %>% filter(casepilot=="VA", ghgspecies=="co2")

#3-way model using gaussian family
m1_gaus_3way<- glmmTMB(formula = dailyflux_trans~status*season*strata + (1|subsite), 
                       data = va_co2,
                       family = gaussian(),
                       dispformula = ~1)

#ONE COMBINATION droped

#Evaluate basic model:
check_convergence(m1_gaus_3way) #OK
check_singularity(m1_gaus_3way,tolerance = 1e-8) #OK
r2(m1_gaus_3way) 
Anova(m1_gaus_3way)
summary(m1_gaus_3way)
res<- simulateResiduals(m1_gaus_3way)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_3way)
em<- emmeans(m1_gaus_3way, ~ status|strata)
summary(em)
pairs(em)
#DECISSION: Residuals fail

#TRY t_family option
m2_t_3way<- glmmTMB(formula = dailyflux_trans~status*season*strata + (1|subsite), 
                    data = va_co2,
                    family = t_family(),
                    dispformula = ~1)

#Evaluate basic model:
check_convergence(m2_t_3way) #OK
check_singularity(m2_t_3way,tolerance = 1e-8) #Singular, ok non-issue
r2(m2_t_3way, tolerance=1e-11) 
Anova(m2_t_3way)
summary(m2_t_3way)
res<- simulateResiduals(m2_t_3way)
plotQQunif(res)
plotResiduals(res)
check_residuals(m2_t_3way)
em<- emmeans(m2_t_3way, ~ status|strata)
summary(em)
pairs(em)
#DECISSION: Better, residuals have issues bust mostly are OK. 

#Compare models 
anova(m1_gaus_3way,m2_t_3way) #gaus vs t_family

#Set as best
specific_va_best_co2<- m2_t_3way

#check the emmean comparison graphically:
em<- emmeans(specific_va_best_co2, ~ status|strata)
ggplot(as.data.frame(em), aes(x = status, y = emmean, col=strata)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
  labs(
    title = "Estimated Marginal Means by Status for CO2",
    x = "Status",
    y = "EMM (model scale)"
  ) +
  theme_minimal()



#Plot obs vs pred 
data.frame(obs=specific_va_best_co2$frame$dailyflux_trans,
           pred=predict(specific_va_best_co2, type = "response"),
           strata=factor(va_co2$strata, levels = c("bare","vegetated", "open water"))) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=strata),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  stat_cor(aes(col=strata),label.x.npc = 0.4, label.y.npc = 0.2)+
  geom_smooth(method = "lm", aes(col=strata), se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  scale_color_manual(values = c(
    "bare" = "#d95f02",      # brown/orange
    "vegetated" = "#2db82d", # green
    "open water" = "blue" # blue
  ), drop = FALSE) +
  labs(title = "specific_va_co2 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)",
       caption = "This model uses strata as fixed effect (season*status*strata)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())


ggsave(filename = "specific_va_co2_Observed_VS_predicted.png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)



##specific_VA_ch4-----

va_ch4<- data4models %>% filter(casepilot=="VA", ghgspecies=="ch4")

#3-way model using gaussian family
m1_gaus_3way<- glmmTMB(formula = dailyflux_trans~status*season*strata + (1|subsite), 
                       data = va_ch4,
                       family = gaussian(),
                       dispformula = ~1)

#ONE COMBINATION droped

#Evaluate basic model:
check_convergence(m1_gaus_3way) #OK
check_singularity(m1_gaus_3way,tolerance = 1e-8) #OK
r2(m1_gaus_3way) 
Anova(m1_gaus_3way)
summary(m1_gaus_3way)
res<- simulateResiduals(m1_gaus_3way)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_3way)
em<- emmeans(m1_gaus_3way, ~ status|strata)
summary(em)
pairs(em)
#DECISSION: Residuals are ok-ish

#TRY t_family option
m2_t_3way<- glmmTMB(formula = dailyflux_trans~status*season*strata + (1|subsite), 
                    data = va_ch4,
                    family = t_family(),
                    dispformula = ~1)

#Evaluate basic model:
check_convergence(m2_t_3way) #OK
check_singularity(m2_t_3way,tolerance = 1e-8) #Singular, ok non-issue
r2(m2_t_3way, tolerance=1e-11) 
Anova(m2_t_3way)
summary(m2_t_3way)
res<- simulateResiduals(m2_t_3way)
plotQQunif(res)
plotResiduals(res)
check_residuals(m2_t_3way)
em<- emmeans(m2_t_3way, ~ status|strata)
summary(em)
pairs(em)
#DECISSION: Better, residuals have issues bust mostly are OK. 

#Compare models 
anova(m1_gaus_3way,m2_t_3way) #gaus vs t_family

#Set as best
specific_va_best_ch4<- m2_t_3way

#check the emmean comparison graphically:
em<- emmeans(specific_va_best_ch4, ~ status|strata)
ggplot(as.data.frame(em), aes(x = strata, y = emmean, col=status)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
  labs(
    title = "Estimated Marginal Means by Status for ch4",
    x = "Status",
    y = "EMM (model scale)"
  ) +
  theme_minimal()


#Plot obs vs pred 
data.frame(obs=specific_va_best_ch4$frame$dailyflux_trans,
           pred=predict(specific_va_best_ch4, type = "response"),
           strata=factor(va_ch4$strata, levels = c("bare","vegetated", "open water"))) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=strata),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  stat_cor(aes(col=strata),label.x.npc = 0.4, label.y.npc = 0.2)+
  geom_smooth(method = "lm", aes(col=strata), se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  scale_color_manual(values = c(
    "bare" = "#d95f02",      # brown/orange
    "vegetated" = "#2db82d", # green
    "open water" = "blue" # blue
  ), drop = FALSE) +
  labs(title = "specific_va_ch4 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)",
       caption = "This model uses strata as fixed effect (season*status*strata)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())


ggsave(filename = "specific_va_ch4_Observed_VS_predicted.png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)




##specific_VA_GWPco2andch4-----

va_GWPco2andch4<- data4models %>% filter(casepilot=="VA", ghgspecies=="GWPco2andch4")

#3-way model using gaussian family
m1_gaus_3way<- glmmTMB(formula = dailyflux_trans~status*season*strata + (1|subsite), 
                       data = va_GWPco2andch4,
                       family = gaussian(),
                       dispformula = ~1)

#ONE COMBINATION droped

#Evaluate basic model:
check_convergence(m1_gaus_3way) #OK
check_singularity(m1_gaus_3way,tolerance = 1e-8) #OK
r2(m1_gaus_3way) 
Anova(m1_gaus_3way)
summary(m1_gaus_3way)
res<- simulateResiduals(m1_gaus_3way)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_3way)
em<- emmeans(m1_gaus_3way, ~ status|strata)
summary(em)
pairs(em)
#DECISSION: BAD RESIDUALS

#TRY t_family option
m2_t_3way<- glmmTMB(formula = dailyflux_trans~status*season*strata + (1|subsite), 
                    data = va_GWPco2andch4,
                    family = t_family(),
                    dispformula = ~1)

#Evaluate basic model:
check_convergence(m2_t_3way) #OK
check_singularity(m2_t_3way,tolerance = 1e-8) #Singular, ok non-issue
r2(m2_t_3way, tolerance=1e-11) 
Anova(m2_t_3way)
summary(m2_t_3way)
res<- simulateResiduals(m2_t_3way)
plotQQunif(res)
plotResiduals(res)
check_residuals(m2_t_3way)
em<- emmeans(m2_t_3way, ~ status|strata)
summary(em)
pairs(em)
#DECISSION: Better, residuals have issues bust mostly are OK. 

#Compare models 
anova(m1_gaus_3way,m2_t_3way) #gaus vs t_family

#Set as best
specific_va_best_GWPco2andch4<- m2_t_3way

#check the emmean comparison graphically:
em<- emmeans(specific_va_best_GWPco2andch4, ~ status|strata)
ggplot(as.data.frame(em), aes(x = strata, y = emmean, col=status)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
  labs(
    title = "Estimated Marginal Means by Status for GWPco2andch4",
    x = "Status",
    y = "EMM (model scale)"
  ) +
  theme_minimal()


#Plot obs vs pred 
data.frame(obs=specific_va_best_GWPco2andch4$frame$dailyflux_trans,
           pred=predict(specific_va_best_GWPco2andch4, type = "response"),
           strata=factor(va_GWPco2andch4$strata, levels = c("bare","vegetated", "open water"))) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=strata),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  stat_cor(aes(col=strata),label.x.npc = 0.4, label.y.npc = 0.2)+
  geom_smooth(method = "lm", aes(col=strata), se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  scale_color_manual(values = c(
    "bare" = "#d95f02",      # brown/orange
    "vegetated" = "#2db82d", # green
    "open water" = "blue" # blue
  ), drop = FALSE) +
  labs(title = "specific_va_GWPco2andch4 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)",
       caption = "This model uses strata as fixed effect (season*status*strata)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())


ggsave(filename = "specific_va_GWPco2andch4_Observed_VS_predicted.png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)


#VA models 2strata (ok)---------
#We saw that strata is not trongly linked to status, we expect very differnet fluxes depending on strata.
#We tried 3-way models with trata, but due to rank deficiencies we cannot calcualte emmeans.
#Here we try re-leveling strata to vegetated or not vegetated (difference between bare and open water should be accounted by seasonality)


##specific_VA_co2-----

va_co2_vegnoveg<- data4models %>% filter(casepilot=="VA", ghgspecies=="co2") %>% 
  mutate(strata=as.factor(if_else(strata=="vegetated","vegetated", "non-vegetated")))

#3-way model using gaussian family
m1_gaus_3way<- glmmTMB(formula = dailyflux_trans~status*season*strata + (1|subsite), 
                       data = va_co2_vegnoveg,
                       family = gaussian(),
                       dispformula = ~1)


#Evaluate basic model:
check_convergence(m1_gaus_3way) #OK
check_singularity(m1_gaus_3way,tolerance = 1e-8) #OK
r2(m1_gaus_3way) 
Anova(m1_gaus_3way)
summary(m1_gaus_3way)
res<- simulateResiduals(m1_gaus_3way)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_3way)
em<- emmeans(m1_gaus_3way, ~ status|strata)
summary(em)
pairs(em)
#DECISSION: Residuals fail

#TRY t_family option
m2_t_3way<- glmmTMB(formula = dailyflux_trans~status*season*strata + (1|subsite), 
                    data = va_co2_vegnoveg,
                    family = t_family(),
                    dispformula = ~1)

#Evaluate basic model:
check_convergence(m2_t_3way) #OK
check_singularity(m2_t_3way,tolerance = 1e-8) #Singular, ok non-issue
r2(m2_t_3way, tolerance=1e-11) 
Anova(m2_t_3way)
summary(m2_t_3way)
res<- simulateResiduals(m2_t_3way)
plotQQunif(res)
plotResiduals(res)
check_residuals(m2_t_3way)
em<- emmeans(m2_t_3way, ~ status|strata)
summary(em)
pairs(em)
#DECISSION: Better, residuals have issues bust mostly are OK. 

#Compare models 
anova(m1_gaus_3way,m2_t_3way) #gaus vs t_family

#Set as best
specific_va_best_co2<- m2_t_3way

#check the emmean comparison graphically:
em<- emmeans(specific_va_best_co2, ~ status|strata)
ggplot(as.data.frame(em), aes(x = status, y = emmean, col=strata)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
  labs(
    title = "Estimated Marginal Means by Status for CO2",
    x = "Status",
    y = "EMM (model scale)"
  ) +
  theme_minimal()

#Only status
em<- emmeans(specific_va_best_co2, ~ status)
ggplot(as.data.frame(em), aes(x = status, y = emmean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
  labs(
    title = "Estimated Marginal Means by Status for CO2",
    x = "Status",
    y = "EMM (model scale)"
  ) +
  theme_minimal()


#Plot obs vs pred 
data.frame(obs=specific_va_best_co2$frame$dailyflux_trans,
           pred=predict(specific_va_best_co2, type = "response"),
           strataveg=as.factor(va_co2_vegnoveg$strata)) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=strataveg),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  stat_cor(aes(col=strataveg),label.x.npc = 0.4, label.y.npc = 0.2)+
  geom_smooth(method = "lm", aes(col=strataveg), se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "specific_va_co2 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)",
       caption = "This model relevels strata into veg/no-veg and takes it as fixed effect (season*status*strataveg)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())


ggsave(filename = "specific_va_co2_Observed_VS_predicted.png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)



##specific_VA_ch4-----

va_ch4_vegnoveg<- data4models %>% filter(casepilot=="VA", ghgspecies=="ch4") %>% 
  mutate(strata=as.factor(if_else(strata=="vegetated","vegetated", "non-vegetated")))

#3-way model using gaussian family
m1_gaus_3way<- glmmTMB(formula = dailyflux_trans~status*season*strata + (1|subsite), 
                       data = va_ch4_vegnoveg,
                       family = gaussian(),
                       dispformula = ~1)

#Evaluate basic model:
check_convergence(m1_gaus_3way) #OK
check_singularity(m1_gaus_3way,tolerance = 1e-8) #OK
r2(m1_gaus_3way) 
Anova(m1_gaus_3way)
summary(m1_gaus_3way)
res<- simulateResiduals(m1_gaus_3way)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_3way)
em<- emmeans(m1_gaus_3way, ~ status|strata)
summary(em)
pairs(em)
#DECISSION: GOOD RESIDUALS


#SKYP, gaussian is good.
#TRY t_family option
m2_t_3way<- glmmTMB(formula = dailyflux_trans~status*season*strata + (1|subsite), 
                    data = va_ch4_vegnoveg,
                    family = t_family(),
                    dispformula = ~1)

#Evaluate basic model:
check_convergence(m2_t_3way) #OK
check_singularity(m2_t_3way,tolerance = 1e-8) #Singular, ok non-issue
r2(m2_t_3way, tolerance=1e-11) 
Anova(m2_t_3way)
summary(m2_t_3way)
res<- simulateResiduals(m2_t_3way)
plotQQunif(res)
plotResiduals(res)
check_residuals(m2_t_3way)
em<- emmeans(m2_t_3way, ~ status|strata)
summary(em)
pairs(em)
#DECISSION: Better, residuals have issues bust mostly are OK. 

#Compare models 
anova(m1_gaus_3way,m2_t_3way) #gaus vs t_family

#Set as best
specific_va_best_ch4<- m1_gaus_3way


#check the emmean comparison graphically:
em<- emmeans(specific_va_best_ch4, ~ status|strata)
ggplot(as.data.frame(em), aes(x = status, y = emmean, col=strata)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  labs(
    title = "Estimated Marginal Means by Status for CH4",
    x = "Status",
    y = "EMM (model scale)"
  ) +
  theme_minimal()


#Plot obs vs pred 
data.frame(obs=specific_va_best_ch4$frame$dailyflux_trans,
           pred=predict(specific_va_best_ch4, type = "response"),
           strataveg=as.factor(va_ch4_vegnoveg$strata)) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=strataveg),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  stat_cor(aes(col=strataveg),label.x.npc = 0.4, label.y.npc = 0.2)+
  geom_smooth(method = "lm", aes(col=strataveg), se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "specific_va_ch4 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)",
       caption = "This model relevels strata into veg/no-veg and takes it as fixed effect (season*status*strataveg)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())


ggsave(filename = "specific_va_ch4_Observed_VS_predicted.png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)



##specific_VA_GWPco2andch4-----

va_GWPco2andch4_vegnoveg<- data4models %>% filter(casepilot=="VA", ghgspecies=="GWPco2andch4") %>% 
  mutate(strata=as.factor(if_else(strata=="vegetated","vegetated", "non-vegetated")))

#3-way model using gaussian family
m1_gaus_3way<- glmmTMB(formula = dailyflux_trans~status*season*strata + (1|subsite), 
                       data = va_GWPco2andch4_vegnoveg,
                       family = gaussian(),
                       dispformula = ~1)


#Evaluate basic model:
check_convergence(m1_gaus_3way) #OK
check_singularity(m1_gaus_3way,tolerance = 1e-8) #OK
r2(m1_gaus_3way) 
Anova(m1_gaus_3way)
summary(m1_gaus_3way)
res<- simulateResiduals(m1_gaus_3way)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_3way)
em<- emmeans(m1_gaus_3way, ~ status|strata)
summary(em)
pairs(em)
#DECISSION: BAD RESIDUASL


#TRY t_family option
m2_t_3way<- glmmTMB(formula = dailyflux_trans~status*season*strata + (1|subsite), 
                    data = va_GWPco2andch4_vegnoveg,
                    family = t_family(),
                    dispformula = ~1)

#Evaluate basic model:
check_convergence(m2_t_3way) #OK
check_singularity(m2_t_3way,tolerance = 1e-8) #Singular, ok non-issue
r2(m2_t_3way, tolerance=1e-11) 
Anova(m2_t_3way)
summary(m2_t_3way)
res<- simulateResiduals(m2_t_3way)
plotQQunif(res)
plotResiduals(res)
check_residuals(m2_t_3way)
em<- emmeans(m2_t_3way, ~ status|strata)
summary(em)
pairs(em)
#DECISSION: Better, residuals have issues bust mostly are OK. 

#Compare models 
anova(m1_gaus_3way,m2_t_3way) #gaus vs t_family

#Set as best
specific_va_best_GWPco2andch4<- m2_t_3way

#check the emmean comparison graphically:
em<- emmeans(specific_va_best_GWPco2andch4, ~ status|strata)
ggplot(as.data.frame(em), aes(x = status, y = emmean, col=strata)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
  labs(
    title = "Estimated Marginal Means by Status for GWP",
    x = "Status",
    y = "EMM (model scale)"
  ) +
  theme_minimal()



#Plot obs vs pred 
data.frame(obs=specific_va_best_GWPco2andch4$frame$dailyflux_trans,
           pred=predict(specific_va_best_GWPco2andch4, type = "response"),
           strataveg=as.factor(va_GWPco2andch4_vegnoveg$strata)) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=strataveg),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  stat_cor(aes(col=strataveg),label.x.npc = 0.4, label.y.npc = 0.2)+
  geom_smooth(method = "lm", aes(col=strataveg), se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "specific_va_GWPco2andch4 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)",
       caption = "This model relevels strata into veg/no-veg and takes it as fixed effect (season*status*strataveg)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())


ggsave(filename = "specific_va_GWPco2andch4_Observed_VS_predicted.png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)







#CU models 2strata (ok)------
#We know that the fluxes are very different depending on the strata, in curonian we have a constant strata compositon, it is very-well suited to test the effect of restoration accounting for strata.
#CU: build 3-way model with the two-level strata to observe strata-specific patterns. 


##specific_CU_co2-----

cu_co2<- data4models %>% filter(casepilot=="CU", ghgspecies=="co2")

#3-way model using gaussian family
m1_gaus_3way<- glmmTMB(formula = dailyflux_trans~status*season*strata + (1|subsite), 
                           data = cu_co2,
                           family = gaussian(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m1_gaus_3way) #OK
check_singularity(m1_gaus_3way,tolerance = 1e-8) #OK
r2(m1_gaus_3way) 
Anova(m1_gaus_3way)
summary(m1_gaus_3way)
res<- simulateResiduals(m1_gaus_3way)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_3way)
em<- emmeans(m1_gaus_3way, ~ status*strata)
summary(em)
pairs(em)
#DECISSION: Residuals fail

#TRY t_family option
m2_t_3way<- glmmTMB(formula = dailyflux_trans~status*season*strata + (1|subsite), 
                       data = cu_co2,
                       family = t_family(),
                       dispformula = ~1)
#Evaluate basic model:
check_convergence(m2_t_3way) #OK
check_singularity(m2_t_3way,tolerance = 1e-8) #Singular, ok non-issue
r2(m2_t_3way, tolerance=1e-11) 
Anova(m2_t_3way)
summary(m2_t_3way)
res<- simulateResiduals(m2_t_3way)
plotQQunif(res)
plotResiduals(res)
check_residuals(m2_t_3way)
em<- emmeans(m2_t_3way, ~ status|strata)
summary(em)
pairs(em)
#DECISSION: Better, residuals have issues bust mostly are OK. 

#Compare models 
anova(m1_gaus_3way,m2_t_3way) #gaus vs t_family

#Set as best
specific_cu_best_co2<- m2_t_3way

#check the emmean comparison graphically:
em<- emmeans(specific_cu_best_co2, ~ status|strata)
ggplot(as.data.frame(em), aes(x = status, y = emmean, col=strata)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
  labs(
    title = "Estimated Marginal Means by Status for CO2",
    x = "Status",
    y = "EMM (model scale)"
  ) +
  theme_minimal()



#Plot obs vs pred 
data.frame(obs=specific_cu_best_co2$frame$dailyflux_trans,
           pred=predict(specific_cu_best_co2, type = "response"),
           strata=factor(cu_co2$strata, levels = c("bare","vegetated", "open water"))) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=strata),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  stat_cor(aes(col=strata),label.x.npc = 0.4, label.y.npc = 0.2)+
  geom_smooth(method = "lm", aes(col=strata), se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  scale_color_manual(values = c(
    "bare" = "#d95f02",      # brown/orange
    "vegetated" = "#2db82d", # green
    "open water" = "blue" # blue
  ), drop = FALSE) +
  labs(title = "specific_CU_co2 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)",
       caption = "This model uses strata as fixed effect (season*status*strata)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())


ggsave(filename = "specific_CU_co2_Observed_VS_predicted.png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)




##specific_CU_ch4-----

cu_ch4<- data4models %>% filter(casepilot=="CU", ghgspecies=="ch4")

#3-way model using gaussian family
m1_gaus_3way<- glmmTMB(formula = dailyflux_trans~status*season*strata + (1|subsite), 
                       data = cu_ch4,
                       family = gaussian(),
                       dispformula = ~1)
#Evaluate basic model:
check_convergence(m1_gaus_3way) #OK
check_singularity(m1_gaus_3way,tolerance = 1e-8) #OK
r2(m1_gaus_3way) 
Anova(m1_gaus_3way)
summary(m1_gaus_3way)
res<- simulateResiduals(m1_gaus_3way)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_3way)
em<- emmeans(m1_gaus_3way, ~ status|strata)
summary(em)
pairs(em)
#DECISSION: Residuals are ok

#TRY t_family option
m2_t_3way<- glmmTMB(formula = dailyflux_trans~status*season*strata + (1|subsite), 
                    data = cu_ch4,
                    family = t_family(),
                    dispformula = ~1)
#Evaluate basic model:
check_convergence(m2_t_3way) #OK
check_singularity(m2_t_3way,tolerance = 1e-8) #OK
r2(m2_t_3way, tolerance=1e-11) 
Anova(m2_t_3way)
summary(m2_t_3way)
res<- simulateResiduals(m2_t_3way)
plotQQunif(res)
plotResiduals(res)
check_residuals(m2_t_3way)
em<- emmeans(m2_t_3way, ~ status|strata)
summary(em)
pairs(em)
#DECISSION: , residuals have issues bust mostly are OK. 

#Compare models 
anova(m1_gaus_3way,m2_t_3way) #gaus vs t_family

#Set as best
specific_cu_best_ch4<- m1_gaus_3way

#check the emmean comparison graphically:
em<- emmeans(specific_cu_best_ch4, ~ status|strata)
ggplot(as.data.frame(em), aes(x = status, y = emmean, col=strata)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  labs(
    title = "Estimated Marginal Means by Status for ch4",
    x = "Status",
    y = "EMM (model scale)"
  ) +
  theme_minimal()



#Plot obs vs pred 
data.frame(obs=specific_cu_best_ch4$frame$dailyflux_trans,
           pred=predict(specific_cu_best_ch4, type = "response"),
           strata=factor(cu_ch4$strata, levels = c("bare","vegetated", "open water"))) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=strata),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  stat_cor(aes(col=strata),label.x.npc = 0.4, label.y.npc = 0.2)+
  geom_smooth(method = "lm", aes(col=strata), se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  scale_color_manual(values = c(
    "bare" = "#d95f02",      # brown/orange
    "vegetated" = "#2db82d", # green
    "open water" = "blue" # blue
  ), drop = FALSE) +
  labs(title = "specific_CU_ch4 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)",
       caption = "This model uses strata as fixed effect (season*status*strata)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())


ggsave(filename = "specific_CU_ch4_Observed_VS_predicted.png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)




##specific_CU_GWPco2andch4-----

cu_GWPco2andch4<- data4models %>% filter(casepilot=="CU", ghgspecies=="GWPco2andch4")

#3-way model using gaussian family
m1_gaus_3way<- glmmTMB(formula = dailyflux_trans~status*season*strata + (1|subsite), 
                       data = cu_GWPco2andch4,
                       family = gaussian(),
                       dispformula = ~1)
#Evaluate basic model:
check_convergence(m1_gaus_3way) #OK
check_singularity(m1_gaus_3way,tolerance = 1e-8) #OK
r2(m1_gaus_3way) 
Anova(m1_gaus_3way)
summary(m1_gaus_3way)
res<- simulateResiduals(m1_gaus_3way)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_3way)
em<- emmeans(m1_gaus_3way, ~ status|strata)
summary(em)
pairs(em)
#DECISSION: BAD residuals

#TRY t_family option
m2_t_3way<- glmmTMB(formula = dailyflux_trans~status*season*strata + (1|subsite), 
                    data = cu_GWPco2andch4,
                    family = t_family(),
                    dispformula = ~1)
#Evaluate basic model:
check_convergence(m2_t_3way) #OK
check_singularity(m2_t_3way,tolerance = 1e-8) #OK
r2(m2_t_3way, tolerance=1e-11) 
Anova(m2_t_3way)
summary(m2_t_3way)
res<- simulateResiduals(m2_t_3way)
plotQQunif(res)
plotResiduals(res)
check_residuals(m2_t_3way)
em<- emmeans(m2_t_3way, ~ status|strata)
summary(em)
pairs(em)
#DECISSION: Better, residuals have issues bust mostly are OK. 

#Compare models 
anova(m1_gaus_3way,m2_t_3way) #gaus vs t_family

#Set as best
specific_cu_best_GWPco2andch4<- m2_t_3way

#check the emmean comparison graphically:
em<- emmeans(specific_cu_best_GWPco2andch4, ~ status|strata)
ggplot(as.data.frame(em), aes(x = status, y = emmean, col=strata)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
  labs(
    title = "Estimated Marginal Means by Status for GWPco2andch4",
    x = "Status",
    y = "EMM (model scale)"
  ) +
  theme_minimal()



#Plot obs vs pred 
data.frame(obs=specific_cu_best_GWPco2andch4$frame$dailyflux_trans,
           pred=predict(specific_cu_best_GWPco2andch4, type = "response"),
           strata=factor(cu_GWPco2andch4$strata, levels = c("bare","vegetated", "open water"))) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=strata),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  stat_cor(aes(col=strata),label.x.npc = 0.4, label.y.npc = 0.2)+
  geom_smooth(method = "lm", aes(col=strata), se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  scale_color_manual(values = c(
    "bare" = "#d95f02",      # brown/orange
    "vegetated" = "#2db82d", # green
    "open water" = "blue" # blue
  ), drop = FALSE) +
  labs(title = "specific_CU_GWPco2andch4 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)",
       caption = "This model uses strata as fixed effect (season*status*strata)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())


ggsave(filename = "specific_CU_GWPco2andch4_Observed_VS_predicted.png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)




#DU models 3strata(ok)------
#We know that the fluxes are very different depending on the strata, dutch delta we have a very balanced strata composition (almost constant), it is very-well suited to test the effect of restoration accounting for strata.

#DU: build 3-way model with the two-level strata to observe strata-specific patterns. 


##specific_DU_co2-----

du_co2<- data4models %>% filter(casepilot=="DU", ghgspecies=="co2")

#3-way model using gaussian family
m1_gaus_3way<- glmmTMB(formula = dailyflux_trans~status*season*strata + (1|subsite), 
                       data = du_co2,
                       family = gaussian(),
                       dispformula = ~1)
#Evaluate basic model:
check_convergence(m1_gaus_3way) #OK
check_singularity(m1_gaus_3way,tolerance = 1e-8) #OK
r2(m1_gaus_3way) 
Anova(m1_gaus_3way)
summary(m1_gaus_3way)
res<- simulateResiduals(m1_gaus_3way)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_3way)
em<- emmeans(m1_gaus_3way, ~ status*strata)
summary(em)
pairs(em)
#DECISSION: Residuals fail

#TRY t_family option
m2_t_3way<- glmmTMB(formula = dailyflux_trans~status*season*strata + (1|subsite), 
                    data = du_co2,
                    family = t_family(),
                    dispformula = ~1)
#Evaluate basic model:
check_convergence(m2_t_3way) #OK
check_singularity(m2_t_3way,tolerance = 1e-8) #Singular, ok non-issue
r2(m2_t_3way, tolerance=1e-11) 
Anova(m2_t_3way)
summary(m2_t_3way)
res<- simulateResiduals(m2_t_3way)
plotQQunif(res)
plotResiduals(res)
check_residuals(m2_t_3way)
em<- emmeans(m2_t_3way, ~ status|strata)
summary(em)
pairs(em)
#DECISSION: Better, residuals have issues bust mostly are OK. 

#Compare models 
anova(m1_gaus_3way,m2_t_3way) #gaus vs t_family

#Set as best
specific_du_best_co2<- m2_t_3way

#check the emmean comparison graphically:
em<- emmeans(specific_du_best_co2, ~ status|strata)
ggplot(as.data.frame(em), aes(x = status, y = emmean, col=strata)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
  labs(
    title = "Estimated Marginal Means by Status for CO2",
    x = "Status",
    y = "EMM (model scale)"
  ) +
  theme_minimal()



#Plot obs vs pred 
data.frame(obs=specific_du_best_co2$frame$dailyflux_trans,
           pred=predict(specific_du_best_co2, type = "response"),
           strata=factor(du_co2$strata, levels = c("bare","vegetated", "open water"))) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=strata),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  stat_cor(aes(col=strata),label.x.npc = 0.4, label.y.npc = 0.2)+
  geom_smooth(method = "lm", aes(col=strata), se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  scale_color_manual(values = c(
    "bare" = "#d95f02",      # brown/orange
    "vegetated" = "#2db82d", # green
    "open water" = "blue" # blue
  ), drop = FALSE) +
  labs(title = "specific_DU_co2 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)",
       caption = "This model uses strata as fixed effect (season*status*strata)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())


ggsave(filename = "specific_DU_co2_Observed_VS_predicted.png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)




##specific_DU_ch4-----

du_ch4<- data4models %>% filter(casepilot=="DU", ghgspecies=="ch4")

#3-way model using gaussian family
m1_gaus_3way<- glmmTMB(formula = dailyflux_trans~status*season*strata + (1|subsite), 
                       data = du_ch4,
                       family = gaussian(),
                       dispformula = ~1)
#Evaluate basic model:
check_convergence(m1_gaus_3way) #OK
check_singularity(m1_gaus_3way,tolerance = 1e-8) #OK
r2(m1_gaus_3way) 
Anova(m1_gaus_3way)
summary(m1_gaus_3way)
res<- simulateResiduals(m1_gaus_3way)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_3way)
em<- emmeans(m1_gaus_3way, ~ status|strata)
summary(em)
pairs(em)
#DECISSION: Residuals OK

#TRY t_family option
m2_t_3way<- glmmTMB(formula = dailyflux_trans~status*season*strata + (1|subsite), 
                    data = du_ch4,
                    family = t_family(),
                    dispformula = ~1)
#Evaluate basic model:
check_convergence(m2_t_3way) #OK
check_singularity(m2_t_3way,tolerance = 1e-8) #OK
r2(m2_t_3way, tolerance=1e-11) 
Anova(m2_t_3way)
summary(m2_t_3way)
res<- simulateResiduals(m2_t_3way)
plotQQunif(res)
plotResiduals(res)
check_residuals(m2_t_3way)
em<- emmeans(m2_t_3way, ~ status|strata)
summary(em)
pairs(em)
#DECISSION: Residuals OK

#Compare models 
anova(m1_gaus_3way,m2_t_3way) #gaus vs t_family

#Set as best
specific_du_best_ch4<- m1_gaus_3way

#check the emmean comparison graphically:
em<- emmeans(specific_du_best_ch4, ~ status|strata)
ggplot(as.data.frame(em), aes(x = status, y = emmean, col=strata)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  labs(
    title = "Estimated Marginal Means by Status for ch4",
    x = "Status",
    y = "EMM (model scale)"
  ) +
  theme_minimal()



#Plot obs vs pred 
data.frame(obs=specific_du_best_ch4$frame$dailyflux_trans,
           pred=predict(specific_du_best_ch4, type = "response"),
           strata=factor(du_ch4$strata, levels = c("bare","vegetated", "open water"))) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=strata),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  stat_cor(aes(col=strata),label.x.npc = 0.4, label.y.npc = 0.2)+
  geom_smooth(method = "lm", aes(col=strata), se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  scale_color_manual(values = c(
    "bare" = "#d95f02",      # brown/orange
    "vegetated" = "#2db82d", # green
    "open water" = "blue" # blue
  ), drop = FALSE) +
  labs(title = "specific_DU_ch4 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)",
       caption = "This model uses strata as fixed effect (season*status*strata)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())


ggsave(filename = "specific_DU_ch4_Observed_VS_predicted.png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)




##specific_DU_GWPco2andch4-----

du_GWPco2andch4<- data4models %>% filter(casepilot=="DU", ghgspecies=="GWPco2andch4")

#3-way model using gaussian family
m1_gaus_3way<- glmmTMB(formula = dailyflux_trans~status*season*strata + (1|subsite), 
                       data = du_GWPco2andch4,
                       family = gaussian(),
                       dispformula = ~1)
#Evaluate basic model:
check_convergence(m1_gaus_3way) #OK
check_singularity(m1_gaus_3way,tolerance = 1e-8) #OK
r2(m1_gaus_3way) 
Anova(m1_gaus_3way)
summary(m1_gaus_3way)
res<- simulateResiduals(m1_gaus_3way)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_3way)
em<- emmeans(m1_gaus_3way, ~ status|strata)
summary(em)
pairs(em)
#DECISSION: BAD residuals

#TRY t_family option
m2_t_3way<- glmmTMB(formula = dailyflux_trans~status*season*strata + (1|subsite), 
                    data = du_GWPco2andch4,
                    family = t_family(),
                    dispformula = ~1)
#Evaluate basic model:
check_convergence(m2_t_3way) #OK
check_singularity(m2_t_3way,tolerance = 1e-8) #OK
r2(m2_t_3way, tolerance=1e-11) 
Anova(m2_t_3way)
summary(m2_t_3way)
res<- simulateResiduals(m2_t_3way)
plotQQunif(res)
plotResiduals(res)
check_residuals(m2_t_3way)
em<- emmeans(m2_t_3way, ~ status|strata)
summary(em)
pairs(em)
#DECISSION: Better, residuals have issues bust mostly are OK. 

#Compare models 
anova(m1_gaus_3way,m2_t_3way) #gaus vs t_family

#Set as best
specific_du_best_GWPco2andch4<- m2_t_3way

#check the emmean comparison graphically:
em<- emmeans(specific_du_best_GWPco2andch4, ~ status|strata)
ggplot(as.data.frame(em), aes(x = status, y = emmean, col=strata)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
  labs(
    title = "Estimated Marginal Means by Status for GWPco2andch4",
    x = "Status",
    y = "EMM (model scale)"
  ) +
  theme_minimal()



#Plot obs vs pred 
data.frame(obs=specific_du_best_GWPco2andch4$frame$dailyflux_trans,
           pred=predict(specific_du_best_GWPco2andch4, type = "response"),
           strata=factor(du_GWPco2andch4$strata, levels = c("bare","vegetated", "open water"))) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=strata),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  stat_cor(aes(col=strata),label.x.npc = 0.4, label.y.npc = 0.2)+
  geom_smooth(method = "lm", aes(col=strata), se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  scale_color_manual(values = c(
    "bare" = "#d95f02",      # brown/orange
    "vegetated" = "#2db82d", # green
    "open water" = "blue" # blue
  ), drop = FALSE) +
  labs(title = "specific_DU_GWPco2andch4 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)",
       caption = "This model uses strata as fixed effect (season*status*strata)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())


ggsave(filename = "specific_DU_GWPco2andch4_Observed_VS_predicted.png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)


#TO ADAPT: DA models--------

#Build model with vegpresence without separating the two altered sites, 


#DA models 4status (ok)-----
#DA: we know that the alterations are very different in this casepilot. Likely with contrasting effects over different ghgs. We will separate status into 4 levels, differenciating the two altered subsites. 

#Prepare data: 
da_4status<- data4models %>% filter(casepilot=="DA") %>% 
  mutate(status=case_when(subsite=="DA-A1"~"Altered 1",
                                 subsite=="DA-A2"~"Altered 2",
                                 status=="Preserved"~"Preserved",
                                 status=="Restored"~"Restored")) %>% 
  mutate(status=as.factor(status))

#first check if strata is still related to status when we separate two alterations
  #Filter and format strata composition data:
da_stratacomp_4status<- stratacomp_formated %>% 
  filter(casepilot=="DA") %>% 
  mutate(status=case_when(subsite=="A1"~"Altered 1",
                                 subsite=="A2"~"Altered 2",
                                 status=="Preserved"~"Preserved",
                                 status=="Restored"~"Restored")) %>% 
  mutate(status=as.factor(status))
  #separate strata comp alone: 
da_comp_df<- da_stratacomp_4status %>% dplyr::select(bare,vegetated,water)

#CHECK: is strata composition dependent of status (separating 2 alterations)?
adonis2(da_comp_df ~ status + season,
  data = da_stratacomp_4status,
  method = "bray",
  permutations = 9999
)
#ANSWER: YES, we cannot use strata in our model (it is highly dependent of status and will mask its effect)

##specific_DA_co2---- 


#Subset data and check transformation: 
specific_da_co2<- da_4status %>% filter(casepilot=="DA"&ghgspecies=="co2")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="DA"&
           ghgspecies=="co2") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = specific_da_co2,
                           family = gaussian(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m1_gaus_nostrata) #OK
check_singularity(m1_gaus_nostrata,tolerance = 1e-8) #OK
r2(m1_gaus_nostrata) 
Anova(m1_gaus_nostrata)
summary(m1_gaus_nostrata)
res<- simulateResiduals(m1_gaus_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_nostrata)
em<- emmeans(m1_gaus_nostrata, ~ status)
summary(em)
pairs(em)
#DECISSION: BAD RESIDUALS 

#TRY t_family option

#Most basic: t_family, only status and season 
m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = specific_da_co2,
                        family = t_family(),
                        dispformula = ~1)
#Evaluate basic model:
check_convergence(m3_t_nostrata) #OK
check_singularity(m3_t_nostrata,tolerance = 1e-8) #Singular, ok non-issue
r2(m3_t_nostrata, tolerance=1e-11) 
Anova(m3_t_nostrata)
summary(m3_t_nostrata)
res<- simulateResiduals(m3_t_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m3_t_nostrata)
em<- emmeans(m3_t_nostrata, ~ status)
summary(em)
pairs(em)
#DECISSION: Better, residuals have issues bust mostly are OK. 

#Compare models 
anova(m1_gaus_nostrata,m3_t_nostrata) #gaus vs t_family

#Set as best
specific_da_best_co2<- m3_t_nostrata

#check the emmean comparison graphically:
em<- emmeans(specific_da_best_co2, ~ status)
ggplot(as.data.frame(em), aes(x = status, y = emmean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
  labs(
    title = "Estimated Marginal Means by Status for CO2",
    x = "Status",
    y = "EMM (model scale)"
  ) +
  theme_minimal()



#Plot obs vs pred 
data.frame(obs=specific_da_best_co2$frame$dailyflux_trans,
           pred=predict(specific_da_best_co2, type = "response"),
           strata=factor(specific_da_co2$strata, levels = c("bare","vegetated", "open water"))) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=strata),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  stat_cor(aes(col=strata),label.x.npc = 0.4, label.y.npc = 0.2)+
  geom_smooth(method = "lm", aes(col=strata), se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  scale_color_manual(values = c(
    "bare" = "#d95f02",      # brown/orange
    "vegetated" = "#2db82d", # green
    "open water" = "blue" # blue
  ), drop = FALSE) +
  labs(title = "specific_DA_co2 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)",caption = "This model separates altered1 and altered2")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())


ggsave(filename = "specific_DA_co2_Observed_VS_predicted.png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)


##specific_DA_ch4---- 

#Subset data and check transformation: 
specific_da_ch4<- da_4status %>% filter(casepilot=="DA"&ghgspecies=="ch4")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="DA"&
           ghgspecies=="ch4") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = specific_da_ch4,
                           family = gaussian(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m1_gaus_nostrata) #OK
check_singularity(m1_gaus_nostrata,tolerance = 1e-8) #OK
r2(m1_gaus_nostrata) 
Anova(m1_gaus_nostrata)
summary(m1_gaus_nostrata)
res<- simulateResiduals(m1_gaus_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_nostrata)
em<- emmeans(m1_gaus_nostrata, ~ status)
summary(em)
pairs(em)
#DECISSION: BAD RESIDUALS 

#TRY t_family options: 

#Most basic: t_family, only status and season 
m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = specific_da_ch4,
                        family = t_family(),
                        dispformula = ~1)
#Evaluate basic model:
check_convergence(m3_t_nostrata) #OK
check_singularity(m3_t_nostrata,tolerance = 1e-8) #Singular, ok non-issue
r2(m3_t_nostrata, tolerance=1e-11) 
Anova(m3_t_nostrata)
summary(m3_t_nostrata)
res<- simulateResiduals(m3_t_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m3_t_nostrata)
em<- emmeans(m3_t_nostrata, ~ status)
summary(em)
pairs(em)
#DECISSION: Better, residuals have issues bust mostly are OK. 

#Compare models 
anova(m1_gaus_nostrata,m3_t_nostrata) #gaus vs t_family

#Set as best
specific_da_best_ch4<- m3_t_nostrata

#check the emmean comparison graphically:
em<- emmeans(specific_da_best_ch4, ~ status)
ggplot(as.data.frame(em), aes(x = status, y = emmean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
  labs(
    title = "Estimated Marginal Means by Status for CH4",
    x = "Status",
    y = "EMM (model scale)"
  ) +
  theme_minimal()



#Plot obs vs pred 
data.frame(obs=specific_da_best_ch4$frame$dailyflux_trans,
           pred=predict(specific_da_best_ch4, type = "response"),
           strata=factor(specific_da_ch4$strata, levels = c("bare","vegetated", "open water"))) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=strata),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  stat_cor(aes(col=strata),label.x.npc = 0.4, label.y.npc = 0.2)+
  geom_smooth(method = "lm", aes(col=strata), se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  scale_color_manual(values = c(
    "bare" = "#d95f02",      # brown/orange
    "vegetated" = "#2db82d", # green
    "open water" = "blue" # blue
  ), drop = FALSE) +
  labs(title = "specific_DA_ch4 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)",caption = "This model separates altered1 and altered2")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())


ggsave(filename = "specific_DA_ch4_Observed_VS_predicted.png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)



##specific_DA_GWPco2andch4---- 


#Subset data and check transformation: 
specific_da_GWPco2andch4<- da_4status %>% filter(casepilot=="DA"&ghgspecies=="GWPco2andch4")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="DA"&
           ghgspecies=="GWPco2andch4") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = specific_da_GWPco2andch4,
                           family = gaussian(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m1_gaus_nostrata) #OK
check_singularity(m1_gaus_nostrata,tolerance = 1e-8) #OK
r2(m1_gaus_nostrata,tolerance = 1e-10) 
Anova(m1_gaus_nostrata)
summary(m1_gaus_nostrata)
res<- simulateResiduals(m1_gaus_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_nostrata)
em<- emmeans(m1_gaus_nostrata, ~ status)
summary(em)
pairs(em)
#DECISSION: BAD RESIDUALS 

#TRY t_family options: 

#Most basic: t_family, only status and season 
m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = specific_da_GWPco2andch4,
                        family = t_family(),
                        dispformula = ~1)
#Evaluate basic model:
check_convergence(m3_t_nostrata) #OK
check_singularity(m3_t_nostrata,tolerance = 1e-8) #Singular, ok non-issue
r2(m3_t_nostrata, tolerance=1e-11) 
Anova(m3_t_nostrata)
summary(m3_t_nostrata)
res<- simulateResiduals(m3_t_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m3_t_nostrata)
em<- emmeans(m3_t_nostrata, ~ status4levels)
summary(em)
pairs(em)
#DECISSION: Better, residuals have issues bust mostly are OK. 

#Compare models 
anova(m1_gaus_nostrata,m3_t_nostrata) #gaus vs t_family

#Set as best
specific_da_best_GWPco2andch4<- m3_t_nostrata

#check the emmean comparison graphically:
em<- emmeans(specific_da_best_GWPco2andch4, ~ status)
ggplot(as.data.frame(em), aes(x = status, y = emmean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
  labs(
    title = "Estimated Marginal Means by Status for GWP",
    x = "Status",
    y = "EMM (model scale)"
  ) +
  theme_minimal()



#Plot obs vs pred 
data.frame(obs=specific_da_best_GWPco2andch4$frame$dailyflux_trans,
           pred=predict(specific_da_best_GWPco2andch4, type = "response"),
           strata=factor(specific_da_GWPco2andch4$strata, levels = c("bare","vegetated", "open water"))) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=strata),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  stat_cor(aes(col=strata),label.x.npc = 0.4, label.y.npc = 0.2)+
  geom_smooth(method = "lm", aes(col=strata), se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  scale_color_manual(values = c(
    "bare" = "#d95f02",      # brown/orange
    "vegetated" = "#2db82d", # green
    "open water" = "blue" # blue
  ), drop = FALSE) +
  labs(title = "specific_DA_GWPco2andch4 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)",caption = "This model separates altered1 and altered2")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())


ggsave(filename = "specific_DA_GWPco2andch4_Observed_VS_predicted.png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)




#Save specific bestmodel-list -------

#Save lists with best-models and homocedastic models, named as casepilot_ghgspecies (needed in order to work appropriately with the summarising functions)

#List of best models for each casepilot (to get significances, pseudoR2s and residual tests)
specific_best_model_list_allghg<- list(
  "CA_co2" = specific_ca_best_co2, #considers veg/noveg strata
  "CU_co2" = specific_cu_best_co2, #considers strata
  "DA_co2" = specific_da_best_co2, #4-level status
  "DU_co2" = specific_du_best_co2, #considers strata
  # "RI_co2" = ri_best_co2, #SKIP 
  "VA_co2" = specific_va_best_co2, #considers veg/noveg strata
  "CA_ch4" = specific_ca_best_ch4, #considers veg/noveg strata
  "CU_ch4" = specific_cu_best_ch4, #considers strata
  "DA_ch4" = specific_da_best_ch4, #4-level status
  "DU_ch4" = specific_du_best_ch4, #considers strata
  # "RI_ch4" = ri_best_ch4, #SKIP
  "VA_ch4" = specific_va_best_ch4,#considers veg/noveg strata
  "CA_GWPco2andch4" = specific_ca_best_GWPco2andch4, #considers veg/noveg strata
  "CU_GWPco2andch4" = specific_cu_best_GWPco2andch4, #considers strata
  "DA_GWPco2andch4" = specific_da_best_GWPco2andch4, #4-level status
  "DU_GWPco2andch4" = specific_du_best_GWPco2andch4,#considers strata
  # "RI_GWPco2andch4" = ri_best_GWPco2andch4, #SKIP
  "VA_GWPco2andch4" = specific_va_best_GWPco2andch4 # consider veg/noveg strata
)


## 3.1. Overall fit -----
#Structure, pseudoR2 of best model:
specific_best_model_fit<- get_pseudoR2s(specific_best_model_list_allghg)


#R2c and R2m of homocedastic model:
specific_homoc_r2<- get_homoc_R2s(specific_best_model_list_allghg)


#DHARMa Residual diagnostics of best model:
specific_resid_diag_summary <- summarize_dharma_diagnostics(specific_best_model_list_allghg)

print(specific_resid_diag_summary)


##3.2. Significance of effects-----
#Effect of status will first be assessed via Anova (best_model,type=3): is there an effect averaging across all seasons? yes/no, is the effect dependent on season? 

#Obtain significance for main effects of best_models
specific_results_anova<-get_anova_results_3way(specific_best_model_list_allghg) 
print(specific_results_anova)



#Intercept refers to whether the grand mean (on the transformed scale) is significantly different than cero (Is the overall net exchange of CH4 balanced between uptake and emission)
#Status refers to whether the status (averaged across seasons) has a significant effect on flux. 
#Season refers to whether the season (averaged across status) has a significant effect on flux.
#Interaction refers to whether the effect of status depends on season — i.e., the difference between status levels changes across seasons.


specific_model_outputs<- specific_results_anova %>% 
  full_join(specific_best_model_fit) %>% 
  full_join(specific_homoc_r2 %>% dplyr::select(dataset, homoced_R2m, homoced_R2c))


#Save in dropbox tables with the models summary outputs: 

write.csv(x = specific_model_outputs, file = paste0(paper_path,"Allghg_summary_Specificchambermodels.csv"),row.names = F)

write.csv(x = specific_resid_diag_summary, file = paste0(paper_path,"Allghg_Residualtests_Specificchambermodels.csv"),row.names = F)


##3.3. EMEANs----

#Notes on non-gaussian models: upper and lower confidence limits are given as asymp.LCL/asymp.UCL (in our context we can treat them equally as the ones from gaussian models). Additionally, non-gaussian models will have Inf df (degrees of freedom are calculated differently for T-family models), not an issue. 

##Emmeans + groupletters-----

# Initialize list for storing comparisons
comparisons_list <- list()


#LOOP has been modified to perform status_within_strata comparisons when strata is a fixed effect

# Loop to extract all comparisons for every casepilot best model
for (dataset in names(specific_best_model_list_allghg)) {
  cp_model <- specific_best_model_list_allghg[[dataset]]
  
  #extract list of fixed effects of model:
  fixed_effects_used<- attr(terms(cp_model), "term.labels")
  
  #extract casepilot name and ghgspecies from model-list names
  casepilot_name<- sub("_.*", "", dataset)
  ghgspecies<- sub(paste0(casepilot_name,"_"),"",dataset)
  
  # Status comparisons: always
  status_emmeans <- emmeans(cp_model, ~status) #ADD weights="proportional" to use the chamber distribution as a weight to calcualte the overall emmean for status. It will also work for season (if we have differing number of observations each season)
  status_comparisons <- cld(status_emmeans, 
                            alpha = 0.05,
                            Letters = letters,
                            adjust = "sidak") %>%
    mutate(comparison = "status",
           ghgspecies = ghgspecies,
           casepilot = casepilot_name,
           season = NA) %>%
    rename(lower.CL = any_of("asymp.LCL"),
           upper.CL = any_of("asymp.UCL"))
  
  # Season comparisons: 
  season_emmeans <- emmeans(cp_model, ~season)
  season_comparisons <- cld(season_emmeans, 
                            alpha = 0.05,
                            Letters = letters,
                            adjust = "sidak") %>%
    mutate(comparison = "season",
           ghgspecies = ghgspecies,
           casepilot = casepilot_name,
           status = NA) %>%
    rename(lower.CL = any_of("asymp.LCL"),
           upper.CL = any_of("asymp.UCL"))
  
  # Interaction comparisons: 
  interaction_emmeans <- emmeans(cp_model, ~status|season)
  interaction_comparisons <- cld(interaction_emmeans, 
                                 alpha = 0.05,
                                 Letters = letters,
                                 adjust = "sidak") %>%
    mutate(comparison = "status_within_season",
           ghgspecies = ghgspecies,
           casepilot = casepilot_name) %>%
    rename(lower.CL = any_of("asymp.LCL"),
           upper.CL = any_of("asymp.UCL"))
  
  if("strata"%in%fixed_effects_used){
    status_within_strata_emmeans<-  emmeans(cp_model, ~status|strata)
    status_within_strata_comparisons<- cld(status_within_strata_emmeans,
                                           alpha=0.05,
                                           Letters=letters,
                                           adjust = "sidak") %>% 
      mutate(comparison="status_within_strata",
             ghgspecies=ghgspecies,
             casepilot=casepilot_name) %>% 
      rename(lower.CL = any_of("asymp.LCL"),
             upper.CL = any_of("asymp.UCL"))
  } else{ status_within_strata_comparisons<-data.frame(emmean=numeric(0))}
  
  # Join and format CP comparisons
  all_comparisons <- status_comparisons %>%
    full_join(season_comparisons) %>%
    full_join(interaction_comparisons) %>%
    full_join(status_within_strata_comparisons) %>% 
    rename(group_letter = .group) %>%
    mutate(group_letter = gsub(" ", "", group_letter),
      
      # Create key for accessing BestNormalize object
      key_trans = paste(casepilot_name, ghgspecies, sep = "_")
    ) %>%
    rowwise() %>% 
    # BACK-TRANSFORMATION of emmean, SE and CL using BestNormalize object
    mutate(
      emmean_bt = if (!is.null(bn_list[[key_trans]])) predict(bn_list[[key_trans]], emmean, inverse = TRUE) else NA_real_,
      SE_bt = if (!is.null(bn_list[[key_trans]])) predict(bn_list[[key_trans]], SE, inverse = TRUE) else NA_real_,
      lower.CL_bt = if (!is.null(bn_list[[key_trans]])) predict(bn_list[[key_trans]], lower.CL, inverse = TRUE) else NA_real_,
      upper.CL_bt = if (!is.null(bn_list[[key_trans]])) predict(bn_list[[key_trans]], upper.CL, inverse = TRUE) else NA_real_
    ) %>%
    
    dplyr::select(any_of(c(
      "ghgspecies", "casepilot","comparison", "status", "season","strata", "df",
      "emmean", "SE", "lower.CL", "upper.CL", "group_letter",
      "emmean_bt", "SE_bt", "lower.CL_bt", "upper.CL_bt"
    ))) %>%
    arrange(ghgspecies, casepilot, comparison, status, season)
  
  # Store casepilot results
  comparisons_list[[dataset]] <- all_comparisons
}

#Remove within-loop objects
rm(all_comparisons, status_comparisons,season_comparisons,interaction_comparisons,
   interaction_emmeans,status_emmeans, season_emmeans)


# Unir todos los resultados en un solo data.frame
specific_emmeans_withgroups <- bind_rows(comparisons_list)

#Save all comparisons:  
write.csv(x = specific_emmeans_withgroups, file = paste0(paper_path,"Allghg_emmeans-posthoc_specificchambermodels.csv"),row.names = F)


##Pairwise contrasts (not adapted)-----
#NOT ADAPTED YET TO WORK WITH SPECIFIC MODELS. 


# Initialize list for storing comparisons
comparisons_pairwise_list <- list()

# Loop to extract all comparisons for every casepilot best model
for (dataset in names(new_best_model_list_allghg)) {
  cp_model <- new_best_model_list_allghg[[dataset]]
  
  # Extract casepilot name and ghgspecies from model-list names
  casepilot_name <- sub("_.*", "", dataset)
  ghgspecies <- sub(paste0(casepilot_name, "_"), "", dataset)
  
  # Status pairwise comparisons
  status_emmeans <- emmeans(cp_model, ~status)
  status_pairs <- pairs(status_emmeans, adjust = "sidak") %>%
    as.data.frame() %>%
    mutate(comparison = "status",
           ghgspecies = ghgspecies,
           casepilot = casepilot_name)
  
  # Season pairwise comparisons
  season_emmeans <- emmeans(cp_model, ~season)
  season_pairs <- pairs(season_emmeans, adjust = "sidak") %>%
    as.data.frame() %>%
    mutate(comparison = "season",
           ghgspecies = ghgspecies,
           casepilot = casepilot_name)
  
  # Interaction pairwise comparisons (status within season)
  interaction_emmeans <- emmeans(cp_model, ~status | season)
  interaction_pairs <- pairs(interaction_emmeans, adjust = "sidak") %>%
    as.data.frame() %>%
    mutate(comparison = "status_within_season",
           ghgspecies = ghgspecies,
           casepilot = casepilot_name)
  
  # Combine all comparisons
  all_comparisons <- bind_rows(status_pairs, season_pairs, interaction_pairs) %>%
    mutate(
      comparison = factor(comparison, levels = c("status", "season", "status_within_season")),
      season = factor(season, levels = c("S1", "S2", "S3", "S4")),
      database = paste(casepilot_name, ghgspecies, sep = "_")
    ) 
  
  #Skip back-transformation, we are calculating for each comparison the difference in emmeans the SE of that difference and its significance, I dont think it should be back-transformed  
  # Store results
  comparisons_pairwise_list[[dataset]] <- all_comparisons
}

# Combine all results
all_pairwise_comparisons <- bind_rows(comparisons_pairwise_list) %>% 
  mutate(model_distribution=if_else(df==Inf, "t_family", "gaussian"),
         test_used=if_else(df==Inf, "Z.tests","T.tests"),
         stat.ratio=if_else(df==Inf, z.ratio, t.ratio)) %>% 
  dplyr::select(casepilot, ghgspecies, database, model_distribution,comparison, test_used, contrast, season, estimate, SE, df, stat.ratio, p.value)


#Save all comparisons:  
write.csv(x = all_pairwise_comparisons, file = paste0(paper_path,"Allghg_pairwise_diffemmeans-posthoc_chambermodels.csv"),row.names = F)







#________-------
#CHECKS from POL--------

#To decide: which is the basic model build (likely the one we already have, check improvement with strata as random nested in subsite)

#To decide: for which casepilots can we build a non-rank deficient model? does it help us in explaining patterns? 


testmodel<- glmmTMB(formula = dailyflux_trans~status*season*strata + (1|subsite), 
                            data = da_co2,
                            family = gaussian(),
                            dispformula = ~1)

Anova(testmodel)
r2(testmodel)
summary(testmodel)
restest<- simulateResiduals(testmodel)
plotQQunif(restest)
check_outliers(restest)
testOutliers(restest)
plotResiduals(restest)
check_residuals(testmodel)

#Options for rank-deficient models: 

#RANK deficient model: Problem due to missing combinations--> missing coeficient--> function tries to calcualte residual that does not exist (resid= obs-esp but in this case it is resid = NA - NA). LOOK FOR ALTERNATIVES TO CHECK RESIDUALS AND see if we can use such a model and trust the significances of effects when some coefficients are missing. 

#We can work around this by manually calculating and plotting the residuals: 

# Extract residuals and fitted values
res <- resid(ca_glmm_gaus_3fixed, type = "pearson")
fitted_vals <- fitted(ca_glmm_gaus_3fixed)

# Create a data frame for ggplot
res_df <- data.frame(
  residuals = res,
  fitted = fitted_vals,
  observed = ca_co2$dailyflux_trans,
  status = ca_co2$status,
  season = ca_co2$season,
  strata = ca_co2$strata
)

ggplot(res_df, aes(x = fitted, y = residuals)) +
  geom_point(alpha = 0.6, color = "gray30") +
  geom_smooth(method = "loess", se = FALSE, color = "blue") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Residuals vs Fitted Values",
       x = "Fitted Values",
       y = "Pearson Residuals") +
  theme_minimal()

ggplot(res_df, aes(sample = residuals)) +
  stat_qq() +
  stat_qq_line(color = "red") +
  labs(title = "Normal Q-Q Plot of Residuals",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles") +
  theme_minimal()

# residual normality
ks.test(res, "pnorm", mean(res, na.rm = TRUE), sd(res, na.rm = TRUE))
shapiro.test(res)

#Obesered vs theoretical quantiles 
obs_q <- quantile(res, probs = seq(0.01, 0.99, 0.01), na.rm = TRUE)

# Theoretical normal quantiles
theo_q <- qnorm(seq(0.01, 0.99, 0.01), mean = mean(res, na.rm = TRUE), sd = sd(res, na.rm = TRUE))

# Plot quantile deviation
library(ggplot2)

qdf <- data.frame(
  obs = obs_q,
  theo = theo_q
)

ggplot(qdf, aes(x = theo, y = obs)) +
  geom_point(color = "darkblue") +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(title = "Quantile Deviation (Observed vs Theoretical)",
       x = "Theoretical Quantiles (Normal)",
       y = "Observed Quantiles (Residuals)") +
  theme_minimal()


em <- emmeans(ca_glmm_gaus_3fixed, ~ status | strata)
summary(em)
pairs(em)

em_status<- emmeans(ca_glmm_gaus_3fixed, ~ status)
pairs(em_status)




res_df %>% 
  ggplot(aes(x=fitted, y=observed))+
  geom_point(aes(col=strata),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=T)+
  stat_cor(aes(col=strata),label.x.npc = 0.4, label.y.npc = 0.2)+
  geom_smooth(method = "lm", aes(col=strata), se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  scale_color_manual(values = c(
    "bare" = "#d95f02",      # brown/orange
    "vegetated" = "#2db82d", # green
    "open water" = "blue" # blue
  ), drop = FALSE) +
  labs(title = "CA_co2 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())



#CO2_______ ------

#Optimize GlMM for transformed data for each casepilot. 

#1. Optimice GLMMs-------

#STEPS: 
#subset data and check transformation effect
#fit model options: 
#1st gaussian, random effects (only subsite or, when appropriate, subsite/strata)
#2nd t_family, if residuals from fist model set fail. 

#check residuals for assumptions
#decide best model (AIC + RESIDUALS)
#store model (best and homocedastic option)

#MODEL options: Need to be tested at the same time. 
#Gausian or T-family
#With or Without heterocedasticity (season or season+status)

#Notes on singularity: to have a model with singularity means that some of the random effect terms likely have variance estimates of (almost) zero, or there is redundancy in the random structure. IN OUR CASE, we should not worry! the random effect is solely the subsite (to account for subsite-specific differences and repeated samplings across seasons). Our random effect structure is the most simple possible (one random effect for intercept), and we expect this random effect to account for variance when the 2 subsites of each status are fundamentally different. With this structure, we should not be too concerned about singuarity. The inclusion of the subsite as a random effect might not be statistically necesary (variance~0) but it is still conceptually justified. 




## 1.1. CA model (ok) ----
##CLEAN AFTER TESTING------
#Re-clean after testing options for model variants. 


#approach: compare fit between gausian models, without heterocedasticity, with for season only and with for season and status. 

ca_co2<- data4models %>% filter(casepilot=="CA"&ghgspecies=="co2")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="CA"&
           ghgspecies=="co2") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()

#Is normality achieved?
#YES, use gaussian family distribution only. 

#FIT MODEL OPTIONS: 

#Most basic: gaussian, only status and season 
ca_glmm_gaus0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = ca_co2,
                        family = gaussian(),
                        dispformula = ~1)
#Evaluate basic model:
check_convergence(ca_glmm_gaus0)
check_singularity(ca_glmm_gaus0,tolerance = 1e-8)
r2(ca_glmm_gaus0)
Anova(ca_glmm_gaus0)
summary(ca_glmm_gaus0)
res<- simulateResiduals(ca_glmm_gaus0)
plotQQunif(res)
plotResiduals(res)
check_residuals(ca_glmm_gaus0)
check_overdispersion(ca_glmm_gaus0)
em<- emmeans(ca_glmm_gaus0, ~ status)
summary(em)
pairs(em)
#DECISSION: Good model, can calculate everything

#Most complex: gaussian, 3fixed effects (status, season, strata)
ca_glmm_gaus_3fixed<- glmmTMB(formula = dailyflux_trans~status*season*strata + (1|subsite), 
                        data = ca_co2,
                        family = gaussian(),
                        dispformula = ~1)

#Evaluate 3fixed model: 
check_convergence(ca_glmm_gaus_3fixed)
check_singularity(ca_glmm_gaus_3fixed,tolerance = 1e-8)
r2(ca_glmm_gaus_3fixed)
Anova(ca_glmm_gaus_3fixed)
summary(ca_glmm_gaus_3fixed)
res<- simulateResiduals(ca_glmm_gaus_3fixed)
plotQQunif(res)
plotResiduals(res)
check_residuals(ca_glmm_gaus_3fixed)
em<- emmeans(ca_glmm_gaus_3fixed, ~ status)
summary(em)
pairs(em)
#RANK deficient model: Problem due to missing combinations--> missing coeficient--> 
#Can obtain some coefficients,
#Cannot calculate residuals using DHARMA or performance package: could do it manually
#CANNOT calcualte contrasts for status.
#DECISSION: 3-way interaction cannot be used: rank-deficiency prevents from calculating contrasts for status (even when averaging across season and strata)


#Intermediate: strata as nested random within subsite
ca_glmm_gaus_nested<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite/strata), 
                        data = ca_co2,
                        family = gaussian(),
                        dispformula = ~1)


#Evaluate nested strata model: 
check_convergence(ca_glmm_gaus_nested)
check_singularity(ca_glmm_gaus_nested,tolerance = 1e-8)
check_singularity(ca_glmm_gaus_nested,tolerance = 1e-9)
r2(ca_glmm_gaus_nested, tolerance=1e-8)
r2(ca_glmm_gaus_nested,tolerance = 1e-9)
summary(ca_glmm_gaus_nested)
Anova(ca_glmm_gaus_nested)
res<- simulateResiduals(ca_glmm_gaus_nested)
plotQQunif(res)
plotResiduals(res)
check_residuals(ca_glmm_gaus_nested)
em<- emmeans(ca_glmm_gaus_nested, ~ status)
summary(em)
pairs(em)
#MAYBE: model is borderline singular (depending on tolerance used), with default tolerance, model is singular and I cannot calculate R2c.Residuals look mostly OK. 


anova(ca_glmm_gaus0,ca_glmm_gaus_nested)
#Nested model is best. (issue with singularity)





ca_glmm_gaus1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = ca_co2,
                        family = gaussian(),
                        dispformula = ~season)

ca_glmm_gaus2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = ca_co2,
                        family = gaussian(),
                        dispformula = ~season+status)

ca_glmm_stu0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = ca_co2,
                       family = t_family,
                       dispformula = ~1) 

ca_glmm_stu1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = ca_co2,
                       family = t_family,
                       dispformula = ~season)

ca_glmm_stu2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = ca_co2,
                       family = t_family,
                       dispformula = ~season+status)

#CALCULATE RESIDUALS
resgaus0<- simulateResiduals(ca_glmm_gaus0)
resgaus1<- simulateResiduals(ca_glmm_gaus1)
resgaus2<- simulateResiduals(ca_glmm_gaus2)

restu0<- simulateResiduals(ca_glmm_stu0)
restu1<- simulateResiduals(ca_glmm_stu1)
restu2<- simulateResiduals(ca_glmm_stu2)


#Check gaussian model
plotQQunif(resgaus0) #check Q-Q
plotResiduals(resgaus0) #check residuals vs fitted, 
plotResiduals(resgaus0, ca_co2$season) #check heterocedasticity
plotResiduals(resgaus0, ca_co2$status) #check heterocedasticity

#Heterocedasticity

plotQQunif(resgaus1) #check Q-Q
plotResiduals(resgaus1) #check residuals vs fitted, 
plotResiduals(resgaus1, ca_co2$season) #check heterocedasticity
plotResiduals(resgaus1, ca_co2$status) #check heterocedasticity
#Heterocedasticity

plotQQunif(resgaus2) #check Q-Q
plotResiduals(resgaus2) #check residuals vs fitted, 
plotResiduals(resgaus2, ca_co2$season) #check heterocedasticity
plotResiduals(resgaus2, ca_co2$status) #check heterocedasticity
#Herocedasticity

#Check t_family model
check_convergence(ca_glmm_stu0)
plot(restu0) #check Q-Q, residuals vs fitted,  
plotResiduals(restu0, ca_co2$season) #check heterocedasticity
plotResiduals(restu0, ca_co2$status) #check heterocedasticity
#DOES NOT CONVERGE

check_convergence(ca_glmm_stu1)
plot(restu1) #check Q-Q, residuals vs fitted, 
plotResiduals(restu1, ca_co2$season) #check heterocedasticity
plotResiduals(restu1, ca_co2$status) #check heterocedasticity
#DOES NOT CONVERGE

check_convergence(ca_glmm_stu2)
plot(restu2) #check Q-Q, residuals vs fitted, 
plotResiduals(restu2, ca_co2$season) #check heterocedasticity
plotResiduals(restu2, ca_co2$status) #check heterocedasticity
#Heterocedasticity

#GAUssian comparison: dispformula options
anova(ca_glmm_gaus0,ca_glmm_gaus1,ca_glmm_gaus2)

#t_family comparison: dispformula options
anova(ca_glmm_stu0,ca_glmm_stu1,ca_glmm_stu2) #ONLY 1 model converges: DO NOT CONSIDER

#best vs best
# anova(ca_glmm_gaus1,ca_glmm_stu1) #GAUS is best, clear


#CA best model?: gaussian, allowing heterocedasticity across seasons
ca_best_co2<- ca_glmm_gaus_nested
Anova(ca_best_co2,type = 3) 

#R2 of homocedastic model
ca_best_homoc_co2<- ca_glmm_gaus0
r2(ca_best_homoc_co2)

#Singularity?
performance::check_singularity(ca_best_co2)
#Outliers?
testOutliers(ca_best_co2)


#plot obs vs predicted
data.frame(obs=ca_best_co2$frame$dailyflux_trans,
           pred=predict(ca_best_co2, type = "response"),
           strata=factor(ca_co2$strata, levels = c("bare","vegetated", "open water"))) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=strata),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  stat_cor(aes(col=strata),label.x.npc = 0.4, label.y.npc = 0.2)+
  geom_smooth(method = "lm", aes(col=strata), se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  scale_color_manual(values = c(
    "bare" = "#d95f02",      # brown/orange
    "vegetated" = "#2db82d", # green
    "open water" = "blue" # blue
  ), drop = FALSE) +
  labs(title = "CA_co2 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

ggsave(filename = "CA_co2_Observed_VS_predicted.png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
       )

#Remove non-best models
rm(ca_glmm_gaus0,ca_glmm_gaus1,ca_glmm_gaus2,
   ca_glmm_stu0,ca_glmm_stu1,ca_glmm_stu2,
   resgaus0,resgaus1,resgaus2,
   restu0, restu1, restu2)


## 1.2. CU model (ok)----

#approach: compare fit between gausian model and tstudent models without heterocedasticity, with for season only and with for season and status (not for interaction). 


cu_co2<- data4models %>% filter(casepilot=="CU"&ghgspecies=="co2")

#See effect of transformation: 
table_trans %>%
  filter(casepilot=="CU"&
           ghgspecies=="co2") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()

#Is normality achieved?
#How important is this? what implications does it have if it fails?

#FIT MODEL OPTIONS: 
cu_glmm_gaus0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = cu_co2,
                        family = gaussian(),
                        dispformula = ~1)

cu_glmm_gaus1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = cu_co2,
                        family = gaussian(),
                        dispformula = ~season)

cu_glmm_gaus2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = cu_co2,
                        family = gaussian(),
                        dispformula = ~season+status)

cu_glmm_stu0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = cu_co2,
                       family = t_family,
                       dispformula = ~1) 

cu_glmm_stu1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = cu_co2,
                       family = t_family,
                       dispformula = ~season)

cu_glmm_stu2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = cu_co2,
                       family = t_family,
                       dispformula = ~season+status)

#CALCULATE RESIDUALS
resgaus0<- simulateResiduals(cu_glmm_gaus0)
resgaus1<- simulateResiduals(cu_glmm_gaus1)
resgaus2<- simulateResiduals(cu_glmm_gaus2)

restu0<- simulateResiduals(cu_glmm_stu0)
restu1<- simulateResiduals(cu_glmm_stu1)
restu2<- simulateResiduals(cu_glmm_stu2)


#Check gaussian model
plotQQunif(resgaus0) #check Q-Q
plotResiduals(resgaus0) #check residuals vs fitted, 
plotResiduals(resgaus0, cu_co2$season) #check heterocedasticity
plotResiduals(resgaus0, cu_co2$status) #check heterocedasticity
#Assumpitons FAIL

plotQQunif(resgaus1) #check Q-Q
plotResiduals(resgaus1) #check residuals vs fitted, 
plotResiduals(resgaus1, cu_co2$season) #check heterocedasticity
plotResiduals(resgaus1, cu_co2$status) #check heterocedasticity
#Assumpitons FAIL

plotQQunif(resgaus2) #check Q-Q
plotResiduals(resgaus2) #check residuals vs fitted, 
plotResiduals(resgaus2, cu_co2$season) #check heterocedasticity
plotResiduals(resgaus2, cu_co2$status) #check heterocedasticity
#Assumpitons FAIL

#Check t_family model
check_convergence(cu_glmm_stu0)
plot(restu0) #check Q-Q, residuals vs fitted,  
plotResiduals(restu0, cu_co2$season) #check heterocedasticity
plotResiduals(restu0, cu_co2$status) #check heterocedasticity
##Assumpitons FAIL

check_convergence(cu_glmm_stu1)
plot(restu1) #check Q-Q, residuals vs fitted, 
plotResiduals(restu1, cu_co2$season) #check heterocedasticity
plotResiduals(restu1, cu_co2$status) #check heterocedasticity
#Assumpitons FAIL

check_convergence(cu_glmm_stu2)
plot(restu2) #check Q-Q, residuals vs fitted, 
plotResiduals(restu2, cu_co2$season) #check heterocedasticity
plotResiduals(restu2, cu_co2$status) #check heterocedasticity
#Heterocedasticity

#GAUssian comparison: dispformula options
anova(cu_glmm_gaus0,cu_glmm_gaus1,cu_glmm_gaus2)

#t_family comparison: dispformula options
anova(cu_glmm_stu0,cu_glmm_stu1,cu_glmm_stu2)

#best vs best
anova(cu_glmm_gaus1,cu_glmm_stu1) #T_family best, clear


#CU best model?: T_family, allowing heterocedasticity across season
cu_best_co2<- cu_glmm_stu1
Anova(cu_best_co2,type = 3) 


#R2 of homocedastic model
cu_best_homoc_co2<- cu_glmm_stu0
r2(cu_best_homoc_co2)

#Singularity?
performance::check_singularity(cu_best_co2)
#Outliers?
testOutliers(cu_best_co2)

#plot obs vs predicted

data.frame(obs=cu_best_co2$frame$dailyflux_trans,
             pred=predict(cu_best_co2, type = "response"),
           strata=factor(cu_co2$strata, levels = c("bare","vegetated", "open water"))) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=strata),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  stat_cor(aes(col=strata),label.x.npc = 0.4, label.y.npc = 0.2)+
  geom_smooth(method = "lm", aes(col=strata), se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  scale_color_manual(values = c(
    "bare" = "#d95f02",      # brown/orange
    "vegetated" = "#2db82d", # green
    "open water" = "blue" # blue
  ), drop = FALSE) +
  labs(title = "CU_co2 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

ggsave(filename = "CU_co2_Observed_VS_predicted.png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)



#Remove non-best models
rm(cu_glmm_gaus0,cu_glmm_gaus1,cu_glmm_gaus2,
   cu_glmm_stu0,cu_glmm_stu1,cu_glmm_stu2,
   resgaus0,resgaus1,resgaus2,
   restu0, restu1, restu2)


## 1.3. DA model (ok) ----
#approach: compare fit between gausian and tstudent models without heterocedasticity, with for season only and with for season and status (not for interaction). 

da_co2<- data4models %>% filter(casepilot=="DA"&ghgspecies=="co2")

#See effect of transformation: 
table_trans %>%
  filter(casepilot=="DA"&
           ghgspecies=="co2") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()

#Is normality achieved?
#How important is this? what implications does it have if it fails?


#FIT MODEL OPTIONS: 
da_glmm_gaus0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = da_co2,
                        family = gaussian(),
                        dispformula = ~1)

da_glmm_gaus1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = da_co2,
                        family = gaussian(),
                        dispformula = ~season)

da_glmm_gaus2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = da_co2,
                        family = gaussian(),
                        dispformula = ~season+status)

da_glmm_stu0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = da_co2,
                       family = t_family,
                       dispformula = ~1) 

da_glmm_stu1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = da_co2,
                       family = t_family,
                       dispformula = ~season)

da_glmm_stu2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = da_co2,
                       family = t_family,
                       dispformula = ~season+status)

#CALCULATE RESIDUALS
resgaus0<- simulateResiduals(da_glmm_gaus0)
resgaus1<- simulateResiduals(da_glmm_gaus1)
resgaus2<- simulateResiduals(da_glmm_gaus2)

restu0<- simulateResiduals(da_glmm_stu0)
restu1<- simulateResiduals(da_glmm_stu1)
restu2<- simulateResiduals(da_glmm_stu2)


#Check gaussian model
plotQQunif(resgaus0) #check Q-Q
plotResiduals(resgaus0) #check residuals vs fitted, 
plotResiduals(resgaus0, da_co2$season) #check heterocedasticity
plotResiduals(resgaus0, da_co2$status) #check heterocedasticity
#Assumpitons FAIL

plotQQunif(resgaus1) #check Q-Q
plotResiduals(resgaus1) #check residuals vs fitted, 
plotResiduals(resgaus1, da_co2$season) #check heterocedasticity
plotResiduals(resgaus1, da_co2$status) #check heterocedasticity
#Assumpitons FAIL

plotQQunif(resgaus2) #check Q-Q
plotResiduals(resgaus2) #check residuals vs fitted, 
plotResiduals(resgaus2, da_co2$season) #check heterocedasticity
plotResiduals(resgaus2, da_co2$status) #check heterocedasticity
#Assumpitons FAIL

#Check t_family model
check_convergence(da_glmm_stu0)
plot(restu0) #check Q-Q, residuals vs fitted,  
plotResiduals(restu0, da_co2$season) #check heterocedasticity
plotResiduals(restu0, da_co2$status) #check heterocedasticity
##Assumpitons FAIL

check_convergence(da_glmm_stu1)
plot(restu1) #check Q-Q, residuals vs fitted, 
plotResiduals(restu1, da_co2$season) #check heterocedasticity
plotResiduals(restu1, da_co2$status) #check heterocedasticity
#Assumpitons FAIL

check_convergence(da_glmm_stu2)
plot(restu2) #check Q-Q, residuals vs fitted, 
plotResiduals(restu2, da_co2$season) #check heterocedasticity
plotResiduals(restu2, da_co2$status) #check heterocedasticity
#All pass

#GAUssian comparison: dispformula options
anova(da_glmm_gaus0,da_glmm_gaus1,da_glmm_gaus2)

#t_family comparison: dispformula options
anova(da_glmm_stu0,da_glmm_stu1,da_glmm_stu2)

#best vs best
anova(da_glmm_gaus1,da_glmm_stu1) #T_family best, clear


#DA best model?: T_family, allowing heterocedasticity across season
da_best_co2<- da_glmm_stu1
Anova(da_best_co2,type = 3) 


#R2 of homocedastic model
da_best_homoc_co2<- da_glmm_stu0
r2(da_best_homoc_co2)

#Singularity?
performance::check_singularity(da_best_co2)
#Outliers?
testOutliers(da_best_co2)

check_overdispersion(da_best_co2)

#plot obs vs predicted
data.frame(obs=da_best_co2$frame$dailyflux_trans,
           pred=predict(da_best_co2, type = "response"),
           strata=factor(da_co2$strata, levels = c("bare","vegetated", "open water"))) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=strata),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  stat_cor(aes(col=strata),label.x.npc = 0.4, label.y.npc = 0.2)+
  geom_smooth(method = "lm", aes(col=strata), se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  scale_color_manual(values = c(
    "bare" = "#d95f02",      # brown/orange
    "vegetated" = "#2db82d", # green
    "open water" = "blue" # blue
  ), drop = FALSE) +
  labs(title = "DA_co2 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

ggsave(filename = "DA_co2_Observed_VS_predicted.png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)

data.frame(obs=da_best_co2$frame$dailyflux_trans,
           pred=predict(da_best_co2, type = "response"),
           strata=factor(da_co2$strata, levels = c("bare","vegetated", "open water")),
           subsite=da_best_co2$frame$subsite) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=subsite),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "DA_co2 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

ggsave(filename = "DA_co2_Observed_VS_predicted(subsite-color).png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)




#Remove non-best models
rm(da_glmm_gaus0,da_glmm_gaus1,da_glmm_gaus2,
   da_glmm_stu0,da_glmm_stu1,da_glmm_stu2,
   resgaus0,resgaus1,resgaus2,
   restu0, restu1, restu2)




## 1.4. DU model (ok)----

#approach: compare fit between gausian models: without heterocedasticity, with for season only and with for season and status. 


du_co2<- data4models %>% filter(casepilot=="DU"&ghgspecies=="co2")

#See effect of transformation: 
table_trans %>%
  filter(casepilot=="DU"&
           ghgspecies=="co2") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Is normality achieved?
#How important is this? what implications does it have if it fails?


#FIT MODEL OPTIONS: 
du_glmm_gaus0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = du_co2,
                        family = gaussian(),
                        dispformula = ~1)

du_glmm_gaus1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = du_co2,
                        family = gaussian(),
                        dispformula = ~season)

du_glmm_gaus2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = du_co2,
                        family = gaussian(),
                        dispformula = ~season+status)

du_glmm_stu0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = du_co2,
                       family = t_family,
                       dispformula = ~1) 

du_glmm_stu1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = du_co2,
                       family = t_family,
                       dispformula = ~season)

du_glmm_stu2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = du_co2,
                       family = t_family,
                       dispformula = ~season+status)

#CALCULATE RESIDUALS
resgaus0<- simulateResiduals(du_glmm_gaus0)
resgaus1<- simulateResiduals(du_glmm_gaus1)
resgaus2<- simulateResiduals(du_glmm_gaus2)

restu0<- simulateResiduals(du_glmm_stu0)
restu1<- simulateResiduals(du_glmm_stu1)
restu2<- simulateResiduals(du_glmm_stu2)


#Check gaussian model
plotQQunif(resgaus0) #check Q-Q
plotResiduals(resgaus0) #check residuals vs fitted, 
plotResiduals(resgaus0, du_co2$season) #check heterocedasticity
plotResiduals(resgaus0, du_co2$status) #check heterocedasticity
#Heterocedasticity

plotQQunif(resgaus1) #check Q-Q
plotResiduals(resgaus1) #check residuals vs fitted, 
plotResiduals(resgaus1, du_co2$season) #check heterocedasticity
plotResiduals(resgaus1, du_co2$status) #check heterocedasticity
#Good

plotQQunif(resgaus2) #check Q-Q
plotResiduals(resgaus2) #check residuals vs fitted, 
plotResiduals(resgaus2, du_co2$season) #check heterocedasticity
plotResiduals(resgaus2, du_co2$status) #check heterocedasticity
#Good

#Check t_family model
check_convergence(du_glmm_stu0)
plot(restu0) #check Q-Q, residuals vs fitted,  
plotResiduals(restu0, du_co2$season) #check heterocedasticity
plotResiduals(restu0, du_co2$status) #check heterocedasticity
#DOES NOT CONVERGE

check_convergence(du_glmm_stu1)
plot(restu1) #check Q-Q, residuals vs fitted, 
plotResiduals(restu1, du_co2$season) #check heterocedasticity
plotResiduals(restu1, du_co2$status) #check heterocedasticity
#DOES NOT CONVERGE

check_convergence(du_glmm_stu2)
plot(restu2) #check Q-Q, residuals vs fitted, 
plotResiduals(restu2, du_co2$season) #check heterocedasticity
plotResiduals(restu2, du_co2$status) #check heterocedasticity
#DOES NOT CONVERGE

#GAUssian comparison: dispformula options
anova(du_glmm_gaus0,du_glmm_gaus1,du_glmm_gaus2)

#t_family comparison: dispformula options
# anova(du_glmm_stu0,du_glmm_stu1,du_glmm_stu2)# SKIP NONE CONVERGES

#best vs best
# anova(du_glmm_gaus1,du_glmm_stu1) #GAUS best, clear


#DU best model?: GAUS, allowing heterocedasticity across season
du_best_co2<- du_glmm_gaus1
Anova(du_best_co2,type = 3) 


#R2 of homocedastic model
du_best_homoc_co2<- du_glmm_gaus0
r2(du_best_homoc_co2)

#Singularity?
performance::check_singularity(du_best_co2)
#Outliers?
testOutliers(du_best_co2)



#plot obs vs predicted
data.frame(obs=du_best_co2$frame$dailyflux_trans,
           pred=predict(du_best_co2, type = "response"),
           strata=factor(du_co2$strata, levels = c("bare","vegetated", "open water"))) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=strata),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  stat_cor(aes(col=strata),label.x.npc = 0.4, label.y.npc = 0.2)+
  geom_smooth(method = "lm", aes(col=strata), se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  scale_color_manual(values = c(
    "bare" = "#d95f02",      # brown/orange
    "vegetated" = "#2db82d", # green
    "open water" = "blue" # blue
  ), drop = FALSE) +
  labs(title = "DU_co2 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

ggsave(filename = "DU_co2_Observed_VS_predicted.png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)


#Remove non-best models
rm(du_glmm_gaus0,du_glmm_gaus1,du_glmm_gaus2,
   du_glmm_stu0,du_glmm_stu1,du_glmm_stu2,
   resgaus0,resgaus1,resgaus2,
   restu0, restu1, restu2)



## 1.5. RI model (ok)----


#approach: compare fit between gausian models: without heterocedasticity, with for season only and with for season and status. 


ri_co2<- data4models %>% filter(casepilot=="RI"&ghgspecies=="co2")

#See effect of transformation: 
table_trans %>%
  filter(casepilot=="RI"&
           ghgspecies=="co2") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()

#Is normality achieved? YES

#Most basic: gaussian, only status and season 
ca_glmm_gaus0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = ri_co2,
                        family = gaussian(),
                        dispformula = ~1)
#Evaluate basic model:
check_convergence(ca_glmm_gaus0)
check_singularity(ca_glmm_gaus0,tolerance = 1e-8)
r2(ca_glmm_gaus0)
Anova(ca_glmm_gaus0)
summary(ca_glmm_gaus0)
res<- simulateResiduals(ca_glmm_gaus0)
plotQQunif(res)
plotResiduals(res)
check_residuals(ca_glmm_gaus0)
check_overdispersion(ca_glmm_gaus0)
em<- emmeans(ca_glmm_gaus0, ~ status)
summary(em)
pairs(em)
#DECISSION: Good model, can calculate everything

#Most complex: gaussian, 3fixed effects (status, season, strata)
ca_glmm_gaus_3fixed<- glmmTMB(formula = dailyflux_trans~status*season*strata + (1|subsite), 
                              data = ri_co2,
                              family = gaussian(),
                              dispformula = ~1)

#Evaluate 3fixed model: 
check_convergence(ca_glmm_gaus_3fixed)
check_singularity(ca_glmm_gaus_3fixed,tolerance = 1e-8)
r2(ca_glmm_gaus_3fixed)
Anova(ca_glmm_gaus_3fixed)
summary(ca_glmm_gaus_3fixed)
res<- simulateResiduals(ca_glmm_gaus_3fixed)
plotQQunif(res)
plotResiduals(res)
check_residuals(ca_glmm_gaus_3fixed)
em<- emmeans(ca_glmm_gaus_3fixed, ~ status)
summary(em)
pairs(em)
#RANK deficient model: missing combinations--> missing coeficient
#Can obtain some coefficients,
#Cannot calculate residuals using DHARMA or performance package: could do it manually
#CANNOT calcualte contrasts for status.
#DECISSION: 3-way interaction cannot be used: rank-deficiency prevents from calculating contrasts for status (even when averaging across season and strata)


#Intermediate: strata as nested random within subsite
ca_glmm_gaus_nested<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite/strata), 
                              data = ri_co2,
                              family = gaussian(),
                              dispformula = ~1)


#Evaluate nested strata model: 
check_convergence(ca_glmm_gaus_nested)
check_singularity(ca_glmm_gaus_nested,tolerance = 1e-8)
check_singularity(ca_glmm_gaus_nested,tolerance = 1e-9)
r2(ca_glmm_gaus_nested, tolerance=1e-8)
r2(ca_glmm_gaus_nested,tolerance = 1e-9)
summary(ca_glmm_gaus_nested)
Anova(ca_glmm_gaus_nested)
res<- simulateResiduals(ca_glmm_gaus_nested)
plotQQunif(res)
plotResiduals(res)
check_residuals(ca_glmm_gaus_nested)
em<- emmeans(ca_glmm_gaus_nested, ~ status)
summary(em)
pairs(em)
#NO: model is borderline (depending on tolerance used), with default tolerance, model is singular and I cannot calculate R2c.Residuals look mostly OK. 
#NOT VALID FOR OUR QUESTION: due to almost-complete correlation of status and strata, strata s random effect negates the explanatory power of status. 


anova(ca_glmm_gaus0,ca_glmm_gaus_nested)
#Nested model is best. (BUT issue with singularity)






#FIT MODEL OPTIONS: 
ri_glmm_gaus0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = ri_co2,
                        family = gaussian(),
                        dispformula = ~1)

ri_glmm_gaus1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = ri_co2,
                        family = gaussian(),
                        dispformula = ~season)

ri_glmm_gaus2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = ri_co2,
                        family = gaussian(),
                        dispformula = ~season+status)

ri_glmm_stu0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = ri_co2,
                       family = t_family,
                       dispformula = ~1) 

ri_glmm_stu1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = ri_co2,
                       family = t_family,
                       dispformula = ~season)

ri_glmm_stu2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = ri_co2,
                       family = t_family,
                       dispformula = ~season+status)

#CALCULATE RESIDUALS
resgaus0<- simulateResiduals(ri_glmm_gaus0)
resgaus1<- simulateResiduals(ri_glmm_gaus1)
resgaus2<- simulateResiduals(ri_glmm_gaus2)

restu0<- simulateResiduals(ri_glmm_stu0)
restu1<- simulateResiduals(ri_glmm_stu1)
restu2<- simulateResiduals(ri_glmm_stu2)


#Check gaussian model
plotQQunif(resgaus0) #check Q-Q
plotResiduals(resgaus0) #check residuals vs fitted, 
plotResiduals(resgaus0, ri_co2$season) #check heterocedasticity
plotResiduals(resgaus0, ri_co2$status) #check heterocedasticity
#Heterocedasticity

plotQQunif(resgaus1) #check Q-Q
plotResiduals(resgaus1) #check residuals vs fitted, 
plotResiduals(resgaus1, ri_co2$season) #check heterocedasticity
plotResiduals(resgaus1, ri_co2$status) #check heterocedasticity
#Heterocedasticity

plotQQunif(resgaus2) #check Q-Q
plotResiduals(resgaus2) #check residuals vs fitted, 
plotResiduals(resgaus2, ri_co2$season) #check heterocedasticity
plotResiduals(resgaus2, ri_co2$status) #check heterocedasticity
#Good

#Check t_family model
check_convergence(ri_glmm_stu0)
plot(restu0) #check Q-Q, residuals vs fitted,  
plotResiduals(restu0, ri_co2$season) #check heterocedasticity
plotResiduals(restu0, ri_co2$status) #check heterocedasticity
#Heterocedasticity

check_convergence(ri_glmm_stu1)
plot(restu1) #check Q-Q, residuals vs fitted, 
plotResiduals(restu1, ri_co2$season) #check heterocedasticity
plotResiduals(restu1, ri_co2$status) #check heterocedasticity
#Heterocedasticity

check_convergence(ri_glmm_stu2)
plot(restu2) #check Q-Q, residuals vs fitted, 
plotResiduals(restu2, ri_co2$season) #check heterocedasticity
plotResiduals(restu2, ri_co2$status) #check heterocedasticity
#Good 


#GAUssian comparison: dispformula options
anova(ri_glmm_gaus0,ri_glmm_gaus1,ri_glmm_gaus2)

#t_family comparison: dispformula options
anova(ri_glmm_stu0,ri_glmm_stu1,ri_glmm_stu2)

#best vs best
anova(ri_glmm_gaus0,ri_glmm_stu0) #GAUS best, clear


#RI best model?: GAUS, Homocedastic
ri_best_co2<- ri_glmm_gaus0
Anova(ri_best_co2,type = 3) 


#R2 of homocedastic model
ri_best_homoc_co2<- ri_glmm_gaus0
r2(ri_best_homoc_co2)

#Singularity?
performance::check_singularity(ri_best_co2)
#Outliers?
testOutliers(ri_best_co2)


#plot obs vs predicted
data.frame(obs=ri_best_co2$frame$dailyflux_trans,
           pred=predict(ri_best_co2, type = "response"),
           strata=factor(ri_co2$strata, levels = c("bare","vegetated", "open water"))) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=strata),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  stat_cor(aes(col=strata),label.x.npc = 0.4, label.y.npc = 0.2)+
  geom_smooth(method = "lm", aes(col=strata), se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  scale_color_manual(values = c(
    "bare" = "#d95f02",      # brown/orange
    "vegetated" = "#2db82d", # green
    "open water" = "blue" # blue
  ), drop = FALSE) +
  labs(title = "RI_co2 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

ggsave(filename = "RI_co2_Observed_VS_predicted.png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)

data.frame(obs=ri_best_co2$frame$dailyflux_trans,
           pred=predict(ri_best_co2, type = "response"),
           strata=factor(ri_co2$strata, levels = c("bare","vegetated", "open water")),
           subsite=ri_best_co2$frame$subsite) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=subsite),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "RI_co2 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

ggsave(filename = "RI_co2_Observed_VS_predicted(subsite-color).png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)



#Remove non-best models
rm(ri_glmm_gaus0,ri_glmm_gaus1,ri_glmm_gaus2,
   ri_glmm_stu0,ri_glmm_stu1,ri_glmm_stu2,
   resgaus0,resgaus1,resgaus2,
   restu0, restu1, restu2)


## 1.6. VA model (ok, bad) ----

#approach: compare fit between gausian and tstudent models without heterocedasticity, with for season only and with for season and status (not for interaction). 


va_co2<- data4models %>% filter(casepilot=="VA"&ghgspecies=="co2")

#See effect of transformation: 
table_trans %>%
  filter(casepilot=="VA"&
           ghgspecies=="co2") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()

#Is normality achieved?
#How important is this? what implications does it have if it fails?


#FIT MODEL OPTIONS: 
va_glmm_gaus0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = va_co2,
                        family = gaussian(),
                        dispformula = ~1)

va_glmm_gaus1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = va_co2,
                        family = gaussian(),
                        dispformula = ~season)

va_glmm_gaus2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = va_co2,
                        family = gaussian(),
                        dispformula = ~season+status)

va_glmm_stu0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = va_co2,
                       family = t_family,
                       dispformula = ~1) 

va_glmm_stu1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = va_co2,
                       family = t_family,
                       dispformula = ~season)

va_glmm_stu2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = va_co2,
                       family = t_family,
                       dispformula = ~season+status)

#CALCULATE RESIDUALS
resgaus0<- simulateResiduals(va_glmm_gaus0)
resgaus1<- simulateResiduals(va_glmm_gaus1)
resgaus2<- simulateResiduals(va_glmm_gaus2)

restu0<- simulateResiduals(va_glmm_stu0)
restu1<- simulateResiduals(va_glmm_stu1)
restu2<- simulateResiduals(va_glmm_stu2)


#Check gaussian model
plotQQunif(resgaus0) #check Q-Q
plotResiduals(resgaus0) #check residuals vs fitted, 
plotResiduals(resgaus0, va_co2$season) #check heterocedasticity
plotResiduals(resgaus0, va_co2$status) #check heterocedasticity
#Assumptions fail

plotQQunif(resgaus1) #check Q-Q
plotResiduals(resgaus1) #check residuals vs fitted, 
plotResiduals(resgaus1, va_co2$season) #check heterocedasticity
plotResiduals(resgaus1, va_co2$status) #check heterocedasticity
#Heterocedasticity

plotQQunif(resgaus2) #check Q-Q
plotResiduals(resgaus2) #check residuals vs fitted, 
plotResiduals(resgaus2, va_co2$season) #check heterocedasticity
plotResiduals(resgaus2, va_co2$status) #check heterocedasticity
#Heterocedasticity

#Check t_family model
check_convergence(va_glmm_stu0)
plot(restu0) #check Q-Q, residuals vs fitted,  
plotResiduals(restu0, va_co2$season) #check heterocedasticity
plotResiduals(restu0, va_co2$status) #check heterocedasticity
#DOES NOT CONVERGE

check_convergence(va_glmm_stu1)
plot(restu1) #check Q-Q, residuals vs fitted, 
plotResiduals(restu1, va_co2$season) #check heterocedasticity
plotResiduals(restu1, va_co2$status) #check heterocedasticity
#Heterocedasticity

check_convergence(va_glmm_stu2)
plot(restu2) #check Q-Q, residuals vs fitted, 
plotResiduals(restu2, va_co2$season) #check heterocedasticity
plotResiduals(restu2, va_co2$status) #check heterocedasticity
#Assumptions fail


#GAUssian comparison: dispformula options
anova(va_glmm_gaus0,va_glmm_gaus1,va_glmm_gaus2)

#t_family comparison: dispformula options
anova(va_glmm_stu0,va_glmm_stu1,va_glmm_stu2) #SOME NOT CONVERGE; others funnel residuals

#best vs best
# anova(va_glmm_gaus1,va_glmm_stu0) #GAUS best, clear


#VA best model?: GAUS, Heteroced across season
va_best_co2<- va_glmm_gaus1
Anova(va_best_co2,type = 3) 


#R2 of homocedastic model
va_best_homoc_co2<- va_glmm_gaus0
r2(va_best_homoc_co2)

#Singularity?
performance::check_singularity(va_best_co2)
#Outliers?
testOutliers(va_best_co2)

#plot obs vs predicted
data.frame(obs=va_best_co2$frame$dailyflux_trans,
           pred=predict(va_best_co2, type = "response"),
           strata=factor(va_co2$strata, levels = c("bare","vegetated", "open water"))) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=strata),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  stat_cor(aes(col=strata),label.x.npc = 0.4, label.y.npc = 0.2)+
  geom_smooth(method = "lm", aes(col=strata), se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  scale_color_manual(values = c(
    "bare" = "#d95f02",      # brown/orange
    "vegetated" = "#2db82d", # green
    "open water" = "blue" # blue
  ), drop = FALSE) +
  labs(title = "VA_co2 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

ggsave(filename = "VA_co2_Observed_VS_predicted.png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)


#Remove non-best models
rm(va_glmm_gaus0,va_glmm_gaus1,va_glmm_gaus2,
   va_glmm_stu0,va_glmm_stu1,va_glmm_stu2,
   resgaus0,resgaus1,resgaus2,
   restu0, restu1, restu2)




#2. Save bestmodel-list -------

#Save lists with best-models and homocedastic models, named as casepilot_ghgspecies (needed in order to work appropriately with the summarising functions)

#Provide list of homocedastic model variants (to estimate R2m and R2c) Cannot be computed for heterocedastic models or singular models
homoced_model_list_co2<- list(
  "CA_co2" = ca_best_homoc_co2,
  "CU_co2" = cu_best_homoc_co2, 
  "DA_co2" = da_best_homoc_co2,
  "DU_co2" = du_best_homoc_co2,
  "RI_co2" = ri_best_homoc_co2, 
  "VA_co2" = va_best_homoc_co2)

#List of best models for each casepilot (to get significances, pseudoR2s and residual tests)
best_model_list_co2<- list(
  "CA_co2" = ca_best_co2,
  "CU_co2" = cu_best_co2, 
  "DA_co2" = da_best_co2,
  "DU_co2" = du_best_co2,
  "RI_co2" = ri_best_co2,
  "VA_co2" = va_best_co2)


#______________________-----

#CH4_______ ------

#FOR CH4, if we wanted to use gamma distributions for modelling, we would need to shift the dataset to avoid negative values. In old script there is a section to calculate the optimal shift to maximize normality when log-transformed (equivalent to what using gamma (link="log") would do).

#HOWEVER: using bestNormalize we are already including the possibility of this option with the a+log(b) transformation option. Additionally, a+log(b) or gamma family modelling is sensitive to the units used for ch4, while pseudolog can adjust for it by the sigma optimization. 

#DECISSION:   WE WONT USE GAMMA FAMILY. 

#1. Optimice GLMMs-------

#We will compare tranformed-gaussian vs t_family models potentially enabling heterocedasticity via dispformula. We will compare residual behaviours (via plots) and model fit (via AIC). 


#Notes on singularity: to have a model with singularity means that some of the random effect terms likely have variance estimates of (almost) zero, or there is redundancy in the random structure. IN our case, the random effect is solely the subsite (to account for subsite-specific differences and repeated samplings across seasons). Our random effect structure is the most simple possible (one random effect for intercept, not for slope), and we expect this random effect to account for variance when the 2 subsites of each status are fundamentally different. With this structure, we should not be too concerned about singuarity. The inclusion of the subsite as a random effect might not be statistically necessary (variance~0) but it is still conceptually justified.  



## 1.1. CA model (ok)----

#approach: compare fit between gausian, t_family models (pseudolog data), without heterocedasticity, with for season only and with for season and status. 


ca_ch4<- data4models %>% filter(casepilot=="CA"&ghgspecies=="ch4")

#See effect of transformation: 
table_trans %>%
  filter(casepilot=="CA"&
           ghgspecies=="ch4") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()

#Is normality achieved?
#How important is this? what implications does it have if it fails?

#FIT MODEL OPTIONS: 
ca_glmm_gaus0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = ca_ch4,
                        family = gaussian(),
                        dispformula = ~1)

ca_glmm_gaus1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = ca_ch4,
                        family = gaussian(),
                        dispformula = ~season)

ca_glmm_gaus2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = ca_ch4,
                        family = gaussian(),
                        dispformula = ~season+status)

ca_glmm_stu0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = ca_ch4,
                       family = t_family,
                       dispformula = ~1) 

ca_glmm_stu1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = ca_ch4,
                       family = t_family,
                       dispformula = ~season)

ca_glmm_stu2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = ca_ch4,
                       family = t_family,
                       dispformula = ~season+status)

#CALCULATE RESIDUALS
resgaus0<- simulateResiduals(ca_glmm_gaus0)
resgaus1<- simulateResiduals(ca_glmm_gaus1)
resgaus2<- simulateResiduals(ca_glmm_gaus2)

restu0<- simulateResiduals(ca_glmm_stu0)
restu1<- simulateResiduals(ca_glmm_stu1)
restu2<- simulateResiduals(ca_glmm_stu2)


#Check gaussian model
plotQQunif(resgaus0) #check Q-Q
plotResiduals(resgaus0) #check residuals vs fitted, 
plotResiduals(resgaus0, ca_ch4$season) #check heterocedasticity
plotResiduals(resgaus0, ca_ch4$status) #check heterocedasticity
#Heterocedasticity

plotQQunif(resgaus1) #check Q-Q
plotResiduals(resgaus1) #check residuals vs fitted, 
plotResiduals(resgaus1, ca_ch4$season) #check heterocedasticity
plotResiduals(resgaus1, ca_ch4$status) #check heterocedasticity
#ALL PASS

plotQQunif(resgaus2) #check Q-Q
plotResiduals(resgaus2) #check residuals vs fitted, 
plotResiduals(resgaus2, ca_ch4$season) #check heterocedasticity
plotResiduals(resgaus2, ca_ch4$status) #check heterocedasticity
#ALL PASS

#Check t_family model
check_convergence(ca_glmm_stu0)
plot(restu0) #check Q-Q, residuals vs fitted,  
plotResiduals(restu0, ca_ch4$season) #check heterocedasticity
plotResiduals(restu0, ca_ch4$status) #check heterocedasticity
#ALL PASS

check_convergence(ca_glmm_stu1)
plot(restu1) #check Q-Q, residuals vs fitted, 
plotResiduals(restu1, ca_ch4$season) #check heterocedasticity
plotResiduals(restu1, ca_ch4$status) #check heterocedasticity
#ALL PASS

check_convergence(ca_glmm_stu2)
plot(restu2) #check Q-Q, residuals vs fitted, 
plotResiduals(restu2, ca_ch4$season) #check heterocedasticity
plotResiduals(restu2, ca_ch4$status) #check heterocedasticity
#ALL PASS

#GAUssian comparison: dispformula options
anova(ca_glmm_gaus0,ca_glmm_gaus1,ca_glmm_gaus2)

#t_family comparison: dispformula options
anova(ca_glmm_stu0,ca_glmm_stu1,ca_glmm_stu2) 

#best vs best
anova(ca_glmm_gaus2,ca_glmm_stu2) #T_family is best, clear


#CA best model?: gaussian, allowing heterocedasticity across seasons
ca_best_ch4<- ca_glmm_stu2
Anova(ca_best_ch4,type = 3) 


#R2 of homocedastic model
ca_best_homoc_ch4<- ca_glmm_stu0
r2(ca_best_homoc_ch4)

#Singularity?
performance::check_singularity(ca_best_ch4)
#Outliers?
testOutliers(ca_best_ch4)

#plot obs vs predicted
data.frame(obs=ca_best_ch4$frame$dailyflux_trans,
           pred=predict(ca_best_ch4, type = "response"),
           strata=factor(ca_ch4$strata, levels = c("bare","vegetated", "open water"))) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=strata),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  stat_cor(aes(col=strata),label.x.npc = 0.4, label.y.npc = 0.2)+
  geom_smooth(method = "lm", aes(col=strata), se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  scale_color_manual(values = c(
    "bare" = "#d95f02",      # brown/orange
    "vegetated" = "#2db82d", # green
    "open water" = "blue" # blue
  ), drop = FALSE) +
  labs(title = "CA_ch4 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

ggsave(filename = "CA_ch4_Observed_VS_predicted.png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)


#Remove non-best models
rm(ca_glmm_gaus0,ca_glmm_gaus1,ca_glmm_gaus2,
   ca_glmm_stu0,ca_glmm_stu1,ca_glmm_stu2,
   resgaus0,resgaus1,resgaus2,
   restu0, restu1, restu2)





## 1.2. CU model (ok)----
#approach: compare fit between gausian, t_family models (trans data), without heterocedasticity, with for season only and with for season and status. 



cu_ch4<- data4models %>% filter(casepilot=="CU"&ghgspecies=="ch4")

#See effect of transformation: 
table_trans %>%
  filter(casepilot=="CU"&
           ghgspecies=="ch4") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()

#Is normality achieved?
#How important is this? what implications does it have if it fails?

#FIT MODEL OPTIONS: 
cu_glmm_gaus0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = cu_ch4,
                        family = gaussian(),
                        dispformula = ~1)

cu_glmm_gaus1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = cu_ch4,
                        family = gaussian(),
                        dispformula = ~season)

cu_glmm_gaus2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = cu_ch4,
                        family = gaussian(),
                        dispformula = ~season+status)

cu_glmm_stu0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = cu_ch4,
                       family = t_family,
                       dispformula = ~1) 

cu_glmm_stu1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = cu_ch4,
                       family = t_family,
                       dispformula = ~season)

cu_glmm_stu2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = cu_ch4,
                       family = t_family,
                       dispformula = ~season+status)

#CALCULATE RESIDUALS
resgaus0<- simulateResiduals(cu_glmm_gaus0)
resgaus1<- simulateResiduals(cu_glmm_gaus1)
resgaus2<- simulateResiduals(cu_glmm_gaus2)

restu0<- simulateResiduals(cu_glmm_stu0)
restu1<- simulateResiduals(cu_glmm_stu1)
restu2<- simulateResiduals(cu_glmm_stu2)


#Check gaussian model
plotQQunif(resgaus0) #check Q-Q
plotResiduals(resgaus0) #check residuals vs fitted, 
plotResiduals(resgaus0, cu_ch4$season) #check heterocedasticity
plotResiduals(resgaus0, cu_ch4$status) #check heterocedasticity
#FAIL assumptions

plotQQunif(resgaus1) #check Q-Q
plotResiduals(resgaus1) #check residuals vs fitted, 
plotResiduals(resgaus1, cu_ch4$season) #check heterocedasticity
plotResiduals(resgaus1, cu_ch4$status) #check heterocedasticity
#FAIL assumptions

plotQQunif(resgaus2) #check Q-Q
plotResiduals(resgaus2) #check residuals vs fitted, 
plotResiduals(resgaus2, cu_ch4$season) #check heterocedasticity
plotResiduals(resgaus2, cu_ch4$status) #check heterocedasticity
#FAIL assumptions

#Check t_family model
check_convergence(cu_glmm_stu0)
plot(restu0) #check Q-Q, residuals vs fitted,  
plotResiduals(restu0, cu_ch4$season) #check heterocedasticity
plotResiduals(restu0, cu_ch4$status) #check heterocedasticity
#FAIL assumptions

check_convergence(cu_glmm_stu1)
plot(restu1) #check Q-Q, residuals vs fitted, 
plotResiduals(restu1, cu_ch4$season) #check heterocedasticity
plotResiduals(restu1, cu_ch4$status) #check heterocedasticity
#FAIL assumptions

check_convergence(cu_glmm_stu2)
plot(restu2) #check Q-Q, residuals vs fitted, 
plotResiduals(restu2, cu_ch4$season) #check heterocedasticity
plotResiduals(restu2, cu_ch4$status) #check heterocedasticity
#FAIL assumptions

#GAUssian comparison: dispformula options
anova(cu_glmm_gaus0,cu_glmm_gaus1,cu_glmm_gaus2)

#t_family comparison: dispformula options
anova(cu_glmm_stu0,cu_glmm_stu1,cu_glmm_stu2) 

#best vs best
anova(cu_glmm_gaus2,cu_glmm_stu2) #T_family is best, clear


#CU best model?: gaussian, allowing heterocedasticity across seasons and status
cu_best_ch4<- cu_glmm_stu2
Anova(cu_best_ch4,type = 3) 


#R2 of homocedastic model
cu_best_homoc_ch4<- cu_glmm_stu0
r2(cu_best_homoc_ch4)

#Singularity?
performance::check_singularity(cu_best_ch4)
#Outliers?
testOutliers(cu_best_ch4)


#plot obs vs predicted
data.frame(obs=cu_best_ch4$frame$dailyflux_trans,
           pred=predict(cu_best_ch4, type = "response"),
           strata=factor(cu_ch4$strata, levels = c("bare","vegetated", "open water"))) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=strata),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  stat_cor(aes(col=strata),label.x.npc = 0.4, label.y.npc = 0.2)+
  geom_smooth(method = "lm", aes(col=strata), se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  scale_color_manual(values = c(
    "bare" = "#d95f02",      # brown/orange
    "vegetated" = "#2db82d", # green
    "open water" = "blue" # blue
  ), drop = FALSE) +
  labs(title = "CU_ch4 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

ggsave(filename = "CU_ch4_Observed_VS_predicted.png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)


#Remove non-best models
rm(cu_glmm_gaus0,cu_glmm_gaus1,cu_glmm_gaus2,
   cu_glmm_stu0,cu_glmm_stu1,cu_glmm_stu2,
   resgaus0,resgaus1,resgaus2,
   restu0, restu1, restu2)




## 1.3. DA model (ok)----
#approach: compare fit between gausian, t_family models (trans data), without heterocedasticity, with for season only and with for season and status. 


da_ch4<- data4models %>% filter(casepilot=="DA"&ghgspecies=="ch4")

#See effect of transformation: 
table_trans %>%
  filter(casepilot=="DA"&
           ghgspecies=="ch4") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()

#Is normality achieved?
#How important is this? what implications does it have if it fails?

#FIT MODEL OPTIONS: 
da_glmm_gaus0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = da_ch4,
                        family = gaussian(),
                        dispformula = ~1)

da_glmm_gaus1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = da_ch4,
                        family = gaussian(),
                        dispformula = ~season)

da_glmm_gaus2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = da_ch4,
                        family = gaussian(),
                        dispformula = ~season+status)

da_glmm_stu0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = da_ch4,
                       family = t_family,
                       dispformula = ~1) 

da_glmm_stu1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = da_ch4,
                       family = t_family,
                       dispformula = ~season)

da_glmm_stu2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = da_ch4,
                       family = t_family,
                       dispformula = ~season+status)

#CALCULATE RESIDUALS
resgaus0<- simulateResiduals(da_glmm_gaus0)
resgaus1<- simulateResiduals(da_glmm_gaus1)
resgaus2<- simulateResiduals(da_glmm_gaus2)

restu0<- simulateResiduals(da_glmm_stu0)
restu1<- simulateResiduals(da_glmm_stu1)
restu2<- simulateResiduals(da_glmm_stu2)

#Check gaussian model
plot(resgaus0) #check Q-Q, residuals vs fitted, 
plotResiduals(resgaus0, da_ch4$season) #check heterocedasticity
plotResiduals(resgaus0, da_ch4$status) #check heterocedasticity
#Assumptions fail, funnel shape residuals, heterocedasticity for both season and status

plot(resgaus1) #check Q-Q, residuals vs fitted, 
plotResiduals(resgaus1, da_ch4$season) #check heterocedasticity
plotResiduals(resgaus1, da_ch4$status) #check heterocedasticity
#Assumptions fail

plot(resgaus2) #check Q-Q, residuals vs fitted, 
plotResiduals(resgaus2, da_ch4$season) #check heterocedasticity
plotResiduals(resgaus2, da_ch4$status) #check heterocedasticity
#Assumptions fail


#Check t_family model
plot(restu0) #check Q-Q, residuals vs fitted,  
plotResiduals(restu0, da_ch4$season) #check heterocedasticity
plotResiduals(restu0, da_ch4$status) #check heterocedasticity
#Assumptions fail

plot(restu1) #check Q-Q, residuals vs fitted, 
plotResiduals(restu1, da_ch4$season) #check heterocedasticity
plotResiduals(restu1, da_ch4$status) #check heterocedasticity
#Assumptions fail

plot(restu2) #check Q-Q, residuals vs fitted, 
plotResiduals(restu2, da_ch4$season) #check heterocedasticity
plotResiduals(restu2, da_ch4$status) #check heterocedasticity
#Assumptions fail


#GAUssian comparison: dispformula options
anova(da_glmm_gaus0,da_glmm_gaus1,da_glmm_gaus2)

#t_family comparison: dispformula options
anova(da_glmm_stu0,da_glmm_stu1,da_glmm_stu2)

#best vs best
anova(da_glmm_gaus2,da_glmm_stu2) #DOUBT t-family model is best, should I use it?

check_singularity(da_glmm_stu2)

#DA best model?: t_family, allowing heterocedasticity across seasons and status (expected, altered will have more dispersion due to higher impact of ebullition). 
da_best_ch4<- da_glmm_stu2
Anova(da_best_ch4,type = 3) 

#R2 of homocedastic model
da_best_homoc_ch4<- da_glmm_stu0
r2(da_best_homoc_ch4)

#Singularity?
performance::check_singularity(da_best_ch4)


#plot obs vs predicted
data.frame(obs=da_best_ch4$frame$dailyflux_trans,
           pred=predict(da_best_ch4, type = "response"),
           strata=factor(da_ch4$strata, levels = c("bare","vegetated", "open water"))) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=strata),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  stat_cor(aes(col=strata),label.x.npc = 0.4, label.y.npc = 0.2)+
  geom_smooth(method = "lm", aes(col=strata), se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  scale_color_manual(values = c(
    "bare" = "#d95f02",      # brown/orange
    "vegetated" = "#2db82d", # green
    "open water" = "blue" # blue
  ), drop = FALSE) +
  labs(title = "DA_ch4 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

ggsave(filename = "DA_ch4_Observed_VS_predicted.png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)


data.frame(obs=da_best_ch4$frame$dailyflux_trans,
           pred=predict(da_best_ch4, type = "response"),
           strata=factor(da_ch4$strata, levels = c("bare","vegetated", "open water")),
           subsite=da_best_ch4$frame$subsite) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=subsite),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "DA_ch4 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

ggsave(filename = "DA_ch4_Observed_VS_predicted(subsite-color).png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)

#Remove non-best models
rm(da_glmm_gaus0,da_glmm_gaus1,da_glmm_gaus2,
   da_glmm_stu0,da_glmm_stu1,da_glmm_stu2,
   resgaus0,resgaus1,resgaus2,
   restu0, restu1, restu2)



## 1.4. DU model (ok)----
#approach: compare fit between gausian, t_family models (trans data), without heterocedasticity, with for season only and with for season and status. 


du_ch4<- data4models %>% filter(casepilot=="DU"&ghgspecies=="ch4")

#See effect of transformation: 
table_trans %>%
  filter(casepilot=="DU"&
           ghgspecies=="ch4") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()

#FIT MODEL OPTIONS: 
du_glmm_gaus0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = du_ch4,
                        family = gaussian(),
                        dispformula = ~1)

du_glmm_gaus1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = du_ch4,
                        family = gaussian(),
                        dispformula = ~season)

du_glmm_gaus2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = du_ch4,
                        family = gaussian(),
                        dispformula = ~season+status)

du_glmm_stu0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = du_ch4,
                       family = t_family,
                       dispformula = ~1) 

du_glmm_stu1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = du_ch4,
                       family = t_family,
                       dispformula = ~season)

du_glmm_stu2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = du_ch4,
                       family = t_family,
                       dispformula = ~season+status)

#CALCULATE RESIDUALS
resgaus0<- simulateResiduals(du_glmm_gaus0)
resgaus1<- simulateResiduals(du_glmm_gaus1)
resgaus2<- simulateResiduals(du_glmm_gaus2)

restu0<- simulateResiduals(du_glmm_stu0)
restu1<- simulateResiduals(du_glmm_stu1)
restu2<- simulateResiduals(du_glmm_stu2)

#Check gaussian model
plot(resgaus0) #check Q-Q, residuals vs fitted, 
plotResiduals(resgaus0, du_ch4$season) #check heterocedasticity
plotResiduals(resgaus0, du_ch4$status) #check heterocedasticity
#Heterocedasticity only 

plot(resgaus1) #check Q-Q, residuals vs fitted, 
plotResiduals(resgaus1, du_ch4$season) #check heterocedasticity
plotResiduals(resgaus1, du_ch4$status) #check heterocedasticity
#Assumptions fail

plot(resgaus2) #check Q-Q, residuals vs fitted, 
plotResiduals(resgaus2, du_ch4$season) #check heterocedasticity
plotResiduals(resgaus2, du_ch4$status) #check heterocedasticity
#Heterocedasticity only


#Check t_family model
plot(restu0) #check Q-Q, residuals vs fitted,  
plotResiduals(restu0, du_ch4$season) #check heterocedasticity
plotResiduals(restu0, du_ch4$status) #check heterocedasticity
#Heterocedasticity only

plot(restu1) #check Q-Q, residuals vs fitted, 
plotResiduals(restu1, du_ch4$season) #check heterocedasticity
plotResiduals(restu1, du_ch4$status) #check heterocedasticity
#Assumptions fail

plot(restu2) #check Q-Q, residuals vs fitted, 
plotResiduals(restu2, du_ch4$season) #check heterocedasticity
plotResiduals(restu2, du_ch4$status) #check heterocedasticity
#OK assumptions


#GAUssian comparison: dispformula options
anova(du_glmm_gaus0,du_glmm_gaus1,du_glmm_gaus2)

#t_family comparison: dispformula options
anova(du_glmm_stu0,du_glmm_stu1,du_glmm_stu2)

#best vs best
anova(du_glmm_gaus2,du_glmm_stu2) #DOUBT t-family model is best, should I use it?

check_singularity(du_glmm_gaus2)

#DU best model: gaussian, allowing heterocedasticity across seasons and status (expected, altered will have more dispersion due to higher impact of ebullition). 
du_best_ch4<- du_glmm_gaus2
Anova(du_best_ch4,type = 3) 

#R2 of homocedastic model
du_best_homoc_ch4<- du_glmm_gaus0
r2(du_best_homoc_ch4)

#Singularity?
performance::check_singularity(du_best_ch4)


#plot obs vs predicted
data.frame(obs=du_best_ch4$frame$dailyflux_trans,
           pred=predict(du_best_ch4, type = "response"),
           strata=factor(du_ch4$strata, levels = c("bare","vegetated", "open water"))) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=strata),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  stat_cor(aes(col=strata),label.x.npc = 0.4, label.y.npc = 0.2)+
  geom_smooth(method = "lm", aes(col=strata), se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  scale_color_manual(values = c(
    "bare" = "#d95f02",      # brown/orange
    "vegetated" = "#2db82d", # green
    "open water" = "blue" # blue
  ), drop = FALSE) +
  labs(title = "DU_ch4 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

ggsave(filename = "DU_ch4_Observed_VS_predicted.png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)


data.frame(obs=du_best_ch4$frame$dailyflux_trans,
           pred=predict(du_best_ch4, type = "response"),
           strata=factor(du_ch4$strata, levels = c("bare","vegetated", "open water")),
           subsite=du_best_ch4$frame$subsite) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=subsite),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "DU_ch4 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

# ggsave(filename = "DU_ch4_Observed_VS_predicted(subsite-color).png", 
#        width = 115,height = 120,units = "mm",device = "png",dpi = 400,
#        path = extraplots_path
# )


#Remove non-best models
rm(du_glmm_gaus0,du_glmm_gaus1,du_glmm_gaus2,
   du_glmm_stu0,du_glmm_stu1,du_glmm_stu2,
   resgaus0,resgaus1,resgaus2,
   restu0, restu1, restu2)



## 1.5. RI model (ok)----
#Without removing outliers, t_family is best. Gaussian models detect 7 outliers. But after, refitting without them, gaussian model becomes singular. 
#DECISSION: keep all data and use t_family

#approach: compare fit between gausian, t_family models (trans data), without heterocedasticity, with for season only and with for season and status. 

ri_ch4<- data4models %>% filter(casepilot=="RI"&ghgspecies=="ch4")

#See effect of transformation: 
table_trans %>%
  filter(casepilot=="RI"&
           ghgspecies=="ch4") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()

#FIT MODEL OPTIONS: 
ri_glmm_gaus0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = ri_ch4,
                        family = gaussian(),
                        dispformula = ~1)

ri_glmm_gaus1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = ri_ch4,
                        family = gaussian(),
                        dispformula = ~season)

ri_glmm_gaus2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = ri_ch4,
                        family = gaussian(),
                        dispformula = ~season+status)

ri_glmm_stu0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = ri_ch4,
                       family = t_family,
                       dispformula = ~1) 

ri_glmm_stu1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = ri_ch4,
                       family = t_family,
                       dispformula = ~season)

ri_glmm_stu2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = ri_ch4,
                       family = t_family,
                       dispformula = ~season+status)

#CALCULATE RESIDUALS
resgaus0<- simulateResiduals(ri_glmm_gaus0)
resgaus1<- simulateResiduals(ri_glmm_gaus1)
resgaus2<- simulateResiduals(ri_glmm_gaus2)

restu0<- simulateResiduals(ri_glmm_stu0)
restu1<- simulateResiduals(ri_glmm_stu1)
restu2<- simulateResiduals(ri_glmm_stu2)

#Check gaussian model
plotQQunif(resgaus0)
plot(resgaus0) #check Q-Q, residuals vs fitted, 
plotResiduals(resgaus0, ri_ch4$season) #check heterocedasticity
plotResiduals(resgaus0, ri_ch4$status) #check heterocedasticity
testOutliers(resgaus0) #7 outliers detected 
#Assumptions fail

plot(resgaus1) #check Q-Q, residuals vs fitted, 
plotResiduals(resgaus1, ri_ch4$season) #check heterocedasticity
plotResiduals(resgaus1, ri_ch4$status) #check heterocedasticity
testOutliers(resgaus1) #7 outliers detected 
#Assumptions fail

plot(resgaus2) #check Q-Q, residuals vs fitted, 
plotResiduals(resgaus2, ri_ch4$season) #check heterocedasticity
plotResiduals(resgaus2, ri_ch4$status) #check heterocedasticity
testOutliers(resgaus2) #6 outliers detected 
#Heterocedasticity only


#Check t_family model
plot(restu0) #check Q-Q, residuals vs fitted,  
plotResiduals(restu0, ri_ch4$season) #check heterocedasticity
plotResiduals(restu0, ri_ch4$status) #check heterocedasticity
testOutliers(restu0) #NO outliers detected (more permisive distribution)
#ALL pass

plot(restu1) #check Q-Q, residuals vs fitted, 
plotResiduals(restu1, ri_ch4$season) #check heterocedasticity
plotResiduals(restu1, ri_ch4$status) #check heterocedasticity
#ALL pass

plot(restu2) #check Q-Q, residuals vs fitted, 
plotResiduals(restu2, ri_ch4$season) #check heterocedasticity
plotResiduals(restu2, ri_ch4$status) #check heterocedasticity
#ALL pass

#GAUssian comparison: dispformula options
anova(ri_glmm_gaus0,ri_glmm_gaus1,ri_glmm_gaus2)

#t_family comparison: dispformula options
anova(ri_glmm_stu0,ri_glmm_stu1,ri_glmm_stu2)

#best vs best
anova(ri_glmm_gaus2,ri_glmm_stu0) #DOUBT t-family model is best, should I use it?

check_singularity(ri_glmm_stu0)

#RI best model: t-family, homocedastic
ri_best_ch4<- ri_glmm_stu0
Anova(ri_best_ch4,type = 3) 

#R2 of homocedastic model
ri_best_homoc_ch4<- ri_glmm_stu0
r2(ri_best_homoc_ch4)

#Singularity?
performance::check_singularity(ri_best_ch4)

#plot obs vs predicted
data.frame(obs=ri_best_ch4$frame$dailyflux_trans,
           pred=predict(ri_best_ch4, type = "response"),
           strata=factor(ri_ch4$strata, levels = c("bare","vegetated", "open water"))) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=strata),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  stat_cor(aes(col=strata),label.x.npc = 0.4, label.y.npc = 0.2)+
  geom_smooth(method = "lm", aes(col=strata), se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  scale_color_manual(values = c(
    "bare" = "#d95f02",      # brown/orange
    "vegetated" = "#2db82d", # green
    "open water" = "blue" # blue
  ), drop = FALSE) +
  labs(title = "RI_ch4 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

ggsave(filename = "RI_ch4_Observed_VS_predicted.png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)


data.frame(obs=ri_best_ch4$frame$dailyflux_trans,
           pred=predict(ri_best_ch4, type = "response"),
           strata=factor(ri_ch4$strata, levels = c("bare","vegetated", "open water")),
           subsite=ri_best_ch4$frame$subsite) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=subsite),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "RI_ch4 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

ggsave(filename = "RI_ch4_Observed_VS_predicted(subsite-color).png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)


#CHECK RE-fit gaussian removing outliers:

out_ch4_ri<-outliers(resgaus0)

ri_ch4_no_outliers <- ri_ch4[-out_ch4_ri, ]

# Refit the model without outliers
ri_glmm_gaus0_no_outliers <- glmmTMB(
  formula = dailyflux_trans ~ status * season + (1 | subsite),
  data = ri_ch4_no_outliers,
  family = gaussian(),
  dispformula = ~1
)

#GREAT IMPROVEMENT OF AIC from original model:
AIC(ri_glmm_gaus0)
AIC(ri_glmm_gaus0_no_outliers)

#GREAT IMPROVEMENT OF AIC from best_model
AIC(ri_best_ch4)
AIC(ri_glmm_gaus0_no_outliers)

#Inspect Residuals model_without outliers: 
resgaus0_noout<- simulateResiduals(ri_glmm_gaus0_no_outliers)
plotQQunif(resgaus0_noout)
plot(resgaus0_noout) #check Q-Q, residuals vs fitted, 
plotResiduals(resgaus0_noout, ri_ch4_no_outliers$season) #check heterocedasticity
plotResiduals(resgaus0_noout, ri_ch4_no_outliers$status) #check heterocedasticity
testOutliers(resgaus0_noout) #No further detected outliers

#BUT model without outliers is singular....
check_singularity(ri_glmm_gaus0_no_outliers)

data.frame(obs=ri_glmm_gaus0_no_outliers$frame$dailyflux_trans,
           pred=predict(ri_glmm_gaus0_no_outliers, type = "response"),
           strata=ri_ch4_no_outliers$strata) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=strata))+
  stat_cor()+
  geom_smooth(method="lm")+
  geom_abline(intercept = 0, slope = 1, col="red")

#DECISSION: we keep outliers and use t_family for modelling. 
#Remove objects from outlier inspection
rm(ri_glmm_gaus0_no_outliers,out_ch4_ri,ri_ch4_no_outliers,resgaus0_noout)

#Remove non-best models
rm(ri_glmm_gaus0,ri_glmm_gaus1,ri_glmm_gaus2,
   ri_glmm_stu0,ri_glmm_stu1,ri_glmm_stu2,
   resgaus0,resgaus1,resgaus2,
   restu0, restu1, restu2)



## 1.6. VA model (ok) ----
#approach: compare fit between gausian, t_family models (trans data), without heterocedasticity, with for season only and with for season and status. 


va_ch4<- data4models %>% filter(casepilot=="VA"&ghgspecies=="ch4")

#See effect of transformation: 
table_trans %>%
  filter(casepilot=="VA"&
           ghgspecies=="ch4") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#FIT MODEL OPTIONS: 
va_glmm_gaus0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = va_ch4,
                        family = gaussian(),
                        dispformula = ~1)

va_glmm_gaus1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = va_ch4,
                        family = gaussian(),
                        dispformula = ~season)

va_glmm_gaus2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = va_ch4,
                        family = gaussian(),
                        dispformula = ~season+status)

va_glmm_stu0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = va_ch4,
                       family = t_family,
                       dispformula = ~1) 

va_glmm_stu1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = va_ch4,
                       family = t_family,
                       dispformula = ~season)

va_glmm_stu2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = va_ch4,
                       family = t_family,
                       dispformula = ~season+status)

#CALCULATE RESIDUALS
resgaus0<- simulateResiduals(va_glmm_gaus0)
resgaus1<- simulateResiduals(va_glmm_gaus1)
resgaus2<- simulateResiduals(va_glmm_gaus2)

restu0<- simulateResiduals(va_glmm_stu0)
restu1<- simulateResiduals(va_glmm_stu1)
restu2<- simulateResiduals(va_glmm_stu2)

#Check gaussian model
plot(resgaus0) #check Q-Q, residuals vs fitted, 
plotResiduals(resgaus0, va_ch4$season) #check heterocedasticity
plotResiduals(resgaus0, va_ch4$status) #check heterocedasticity
#Heterocedasticity

plot(resgaus1) #check Q-Q, residuals vs fitted, 
plotResiduals(resgaus1, va_ch4$season) #check heterocedasticity
plotResiduals(resgaus1, va_ch4$status) #check heterocedasticity
check_convergence(va_glmm_gaus1)
#Heterocedasticity, and step failure in fitting (unclear)

plot(resgaus2) #check Q-Q, residuals vs fitted, 
plotResiduals(resgaus2, va_ch4$season) #check heterocedasticity
plotResiduals(resgaus2, va_ch4$status) #check heterocedasticity
#ALL pass

#Check t_family model
plot(restu0) #check Q-Q, residuals vs fitted,  
plotResiduals(restu0, va_ch4$season) #check heterocedasticity
plotResiduals(restu0, va_ch4$status) #check heterocedasticity
#Heterocedasticity

plot(restu1) #check Q-Q, residuals vs fitted, 
plotResiduals(restu1, va_ch4$season) #check heterocedasticity
plotResiduals(restu1, va_ch4$status) #check heterocedasticity
#Heterocedasticity

plot(restu2) #check Q-Q, residuals vs fitted, 
plotResiduals(restu2, va_ch4$season) #check heterocedasticity
plotResiduals(restu2, va_ch4$status) #check heterocedasticity
#ALL pass


#GAUssian comparison: dispformula options
anova(va_glmm_gaus0,va_glmm_gaus1,va_glmm_gaus2)

#t_family comparison: dispformula options
anova(va_glmm_stu0,va_glmm_stu1,va_glmm_stu2)

#best vs best
anova(va_glmm_gaus2,va_glmm_stu1) #DOUBT t-family model is best, should I use it?

check_singularity(va_glmm_stu1)

#VA best model?: t-family, allowing heterocedasticity across seasons only. 
va_best_ch4<- va_glmm_stu1
Anova(va_best_ch4,type = 3) 

#R2 of homocedastic model
va_best_homoc_ch4<- va_glmm_stu0
r2(va_best_homoc_ch4)

#Singularity?
performance::check_singularity(va_best_ch4)



#plot obs vs predicted
data.frame(obs=va_best_ch4$frame$dailyflux_trans,
           pred=predict(va_best_ch4, type = "response"),
           strata=factor(va_ch4$strata, levels = c("bare","vegetated", "open water"))) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=strata),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  stat_cor(aes(col=strata),label.x.npc = 0.4, label.y.npc = 0.2)+
  geom_smooth(method = "lm", aes(col=strata), se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  scale_color_manual(values = c(
    "bare" = "#d95f02",      # brown/orange
    "vegetated" = "#2db82d", # green
    "open water" = "blue" # blue
  ), drop = FALSE) +
  labs(title = "VA_ch4 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

ggsave(filename = "VA_ch4_Observed_VS_predicted.png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)



#Remove non-best models
rm(va_glmm_gaus0,va_glmm_gaus1,va_glmm_gaus2,
   va_glmm_stu0,va_glmm_stu1,va_glmm_stu2,
   resgaus0,resgaus1,resgaus2,
   restu0, restu1, restu2)



#2. Save bestmodel-list -------

#Save lists with best-models and homocedastic models, named as casepilot_ghgspecies (needed in order to work appropriately with the summarising functions)

#Provide list of homocedastic model variants (to estimate R2m and R2c) Cannot be computed for heterocedastic models or singular models
homoced_model_list_ch4<- list(
  "CA_ch4" = ca_best_homoc_ch4,
  "CU_ch4" = cu_best_homoc_ch4, 
  "DA_ch4" = da_best_homoc_ch4,
  "DU_ch4" = du_best_homoc_ch4,
  "RI_ch4" = ri_best_homoc_ch4, 
  "VA_ch4" = va_best_homoc_ch4)

#List of best models for each casepilot (to get significances, pseudoR2s and residual tests)
best_model_list_ch4<- list(
  "CA_ch4" = ca_best_ch4,
  "CU_ch4" = cu_best_ch4, 
  "DA_ch4" = da_best_ch4,
  "DU_ch4" = du_best_ch4,
  "RI_ch4" = ri_best_ch4,
  "VA_ch4" = va_best_ch4)


#______________________-----

#GWPco2andch4_______ ------
#IMPORTANT, crossing of NAs for Co2 and Ch4 might cause strata distributions to become biased or non-balanced between subsites/status in different casepilots. Especially critical for systems with a lot of ebullition (which often causes artefacts in CO2 leading to NA fluxes)

#1. Optimice GLMMs-------


## 1.1. CA model (ok)----

#approach: compare fit between gausian models, without heterocedasticity, with for season only and with for season and status. 

ca_GWPco2andch4<- data4models %>% filter(casepilot=="CA"&ghgspecies=="GWPco2andch4")

#See effect of transformation: 
table_trans %>%
  filter(casepilot=="CA"&
           ghgspecies=="GWPco2andch4") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#FIT MODEL OPTIONS: 
ca_glmm_gaus0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = ca_GWPco2andch4,
                        family = gaussian(),
                        dispformula = ~1)

ca_glmm_gaus1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = ca_GWPco2andch4,
                        family = gaussian(),
                        dispformula = ~season)

ca_glmm_gaus2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = ca_GWPco2andch4,
                        family = gaussian(),
                        dispformula = ~season+status)

ca_glmm_stu0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = ca_GWPco2andch4,
                       family = t_family,
                       dispformula = ~1) 

ca_glmm_stu1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = ca_GWPco2andch4,
                       family = t_family,
                       dispformula = ~season)

ca_glmm_stu2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = ca_GWPco2andch4,
                       family = t_family,
                       dispformula = ~season+status)

#CALCULATE RESIDUALS
resgaus0<- simulateResiduals(ca_glmm_gaus0)
resgaus1<- simulateResiduals(ca_glmm_gaus1)
resgaus2<- simulateResiduals(ca_glmm_gaus2)

restu0<- simulateResiduals(ca_glmm_stu0)
restu1<- simulateResiduals(ca_glmm_stu1)
restu2<- simulateResiduals(ca_glmm_stu2)


#Check gaussian model
plot(resgaus0) #check Q-Q, residuals vs fitted, 
plotResiduals(resgaus0, ca_GWPco2andch4$season) #check heterocedasticity
plotResiduals(resgaus0, ca_GWPco2andch4$status) #check heterocedasticity
#Heterocedasticity

plot(resgaus1) #check Q-Q, residuals vs fitted, 
plotResiduals(resgaus1, ca_GWPco2andch4$season) #check heterocedasticity
plotResiduals(resgaus1, ca_GWPco2andch4$status) #check heterocedasticity
#Heterocedasticity

plot(resgaus2) #check Q-Q, residuals vs fitted, 
plotResiduals(resgaus2, ca_GWPco2andch4$season) #check heterocedasticity
plotResiduals(resgaus2, ca_GWPco2andch4$status) #check heterocedasticity
#ALL pass

#Check t_family model
plot(restu0) #check Q-Q, residuals vs fitted,  
plotResiduals(restu0, ca_GWPco2andch4$season) #check heterocedasticity
plotResiduals(restu0, ca_GWPco2andch4$status) #check heterocedasticity
#Heterocedasticity

plot(restu1) #check Q-Q, residuals vs fitted, 
plotResiduals(restu1, ca_GWPco2andch4$season) #check heterocedasticity
plotResiduals(restu1, ca_GWPco2andch4$status) #check heterocedasticity
#ALL pass

plot(restu2) #check Q-Q, residuals vs fitted, 
plotResiduals(restu2, ca_GWPco2andch4$season) #check heterocedasticity
plotResiduals(restu2, ca_GWPco2andch4$status) #check heterocedasticity
check_convergence(ca_glmm_stu2)#DOES NOT CONVERGE


#GAUssian comparison: dispformula options
anova(ca_glmm_gaus0,ca_glmm_gaus1,ca_glmm_gaus2)

#t_family comparison: dispformula options
anova(ca_glmm_stu0,ca_glmm_stu1,ca_glmm_stu2)

#best vs best
anova(ca_glmm_gaus1,ca_glmm_stu1) #GAUS is best, clear


#CA best model?: gaussian, allowing heterocedasticity across seasons and status (expected, altered will have more dispersion due to higher impact of ebullition). 
ca_best_GWPco2andch4<- ca_glmm_gaus1
Anova(ca_best_GWPco2andch4,type = 3) 

#R2 of homocedastic model
ca_best_homoc_GWPco2andch4<- ca_glmm_gaus0
r2(ca_best_homoc_GWPco2andch4)

#Singularity?
performance::check_singularity(ca_best_GWPco2andch4)



#plot obs vs predicted
data.frame(obs=ca_best_GWPco2andch4$frame$dailyflux_trans,
           pred=predict(ca_best_GWPco2andch4, type = "response"),
           strata=factor(ca_GWPco2andch4$strata, levels = c("bare","vegetated", "open water"))) %>%
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=strata),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  stat_cor(aes(col=strata),label.x.npc = 0.4, label.y.npc = 0.2)+
  geom_smooth(method = "lm", aes(col=strata), se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  scale_color_manual(values = c(
    "bare" = "#d95f02",      # brown/orange
    "vegetated" = "#2db82d", # green
    "open water" = "blue" # blue
  ), drop = FALSE) +
  labs(title = "CA_gwp Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

ggsave(filename = "CA_gwp_Observed_VS_predicted.png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)


data.frame(obs=ca_best_GWPco2andch4$frame$dailyflux_trans,
           pred=predict(ca_best_GWPco2andch4, type = "response"),
           strata=factor(ca_GWPco2andch4$strata, levels = c("bare","vegetated", "open water")),
           subsite=ca_best_GWPco2andch4$frame$subsite) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=subsite),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "CA_gwp Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

# ggsave(filename = "CA_gwp_Observed_VS_predicted(subsite-color).png", 
#        width = 115,height = 120,units = "mm",device = "png",dpi = 400,
#        path = extraplots_path
# )

#Remove non-best models
rm(ca_glmm_gaus0,ca_glmm_gaus1,ca_glmm_gaus2,
   ca_glmm_stu0,ca_glmm_stu1,ca_glmm_stu2,
   resgaus0,resgaus1,resgaus2,
   restu0, restu1, restu2)


## 1.2. CU model (ok)----
#approach: compare fit between gausian models, without heterocedasticity, with for season only and with for season and status. 

cu_GWPco2andch4<- data4models %>% filter(casepilot=="CU"&ghgspecies=="GWPco2andch4")

#See effect of transformation: 
table_trans %>%
  filter(casepilot=="CU"&
           ghgspecies=="GWPco2andch4") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#FIT MODEL OPTIONS: 
cu_glmm_gaus0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = cu_GWPco2andch4,
                        family = gaussian(),
                        dispformula = ~1)

cu_glmm_gaus1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = cu_GWPco2andch4,
                        family = gaussian(),
                        dispformula = ~season)

cu_glmm_gaus2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = cu_GWPco2andch4,
                        family = gaussian(),
                        dispformula = ~season+status)

cu_glmm_stu0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = cu_GWPco2andch4,
                       family = t_family,
                       dispformula = ~1) 

cu_glmm_stu1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = cu_GWPco2andch4,
                       family = t_family,
                       dispformula = ~season)

cu_glmm_stu2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = cu_GWPco2andch4,
                       family = t_family,
                       dispformula = ~season+status)

#CALCULATE RESIDUALS
resgaus0<- simulateResiduals(cu_glmm_gaus0)
resgaus1<- simulateResiduals(cu_glmm_gaus1)
resgaus2<- simulateResiduals(cu_glmm_gaus2)

restu0<- simulateResiduals(cu_glmm_stu0)
restu1<- simulateResiduals(cu_glmm_stu1)
restu2<- simulateResiduals(cu_glmm_stu2)

#Check gaussian model
plot(resgaus0) #check Q-Q, residuals vs fitted, 
plotResiduals(resgaus0, cu_GWPco2andch4$season) #check heterocedasticity
plotResiduals(resgaus0, cu_GWPco2andch4$status) #check heterocedasticity
#Assumptions fail

plot(resgaus1) #check Q-Q, residuals vs fitted, 
plotResiduals(resgaus1, cu_GWPco2andch4$season) #check heterocedasticity
plotResiduals(resgaus1, cu_GWPco2andch4$status) #check heterocedasticity
#Assumptions fail

plot(resgaus2) #check Q-Q, residuals vs fitted, 
plotResiduals(resgaus2, cu_GWPco2andch4$season) #check heterocedasticity
plotResiduals(resgaus2, cu_GWPco2andch4$status) #check heterocedasticity
#Assumptions fail

#Check t_family model
plot(restu0) #check Q-Q, residuals vs fitted,  
plotResiduals(restu0, cu_GWPco2andch4$season) #check heterocedasticity
plotResiduals(restu0, cu_GWPco2andch4$status) #check heterocedasticity
#Heterocedasticity + funnel resid

plot(restu1) #check Q-Q, residuals vs fitted, 
plotResiduals(restu1, cu_GWPco2andch4$season) #check heterocedasticity
plotResiduals(restu1, cu_GWPco2andch4$status) #check heterocedasticity
#Heterocedasticity + funnel resid

plot(restu2) #check Q-Q, residuals vs fitted, 
plotResiduals(restu2, cu_GWPco2andch4$season) #check heterocedasticity
plotResiduals(restu2, cu_GWPco2andch4$status) #check heterocedasticity
#ALL PASS

#GAUssian comparison: dispformula options
anova(cu_glmm_gaus0,cu_glmm_gaus1,cu_glmm_gaus2)

#t_family comparison: dispformula options
anova(cu_glmm_stu0,cu_glmm_stu1,cu_glmm_stu2)

#best vs best
anova(cu_glmm_gaus2,cu_glmm_stu2) #GAUS is best clear


#CU best model: T-family, allowing heterocedasticity across seasons and status (expected, altered will have more dispersion due to higher impact of ebullition). 
cu_best_GWPco2andch4<- cu_glmm_stu2
Anova(cu_best_GWPco2andch4,type = 3) 

#R2 of homocedastic model
cu_best_homoc_GWPco2andch4<- cu_glmm_stu0
r2(cu_best_homoc_GWPco2andch4)

#Singularity?
performance::check_singularity(cu_best_GWPco2andch4)



#plot obs vs predicted
data.frame(obs=cu_best_GWPco2andch4$frame$dailyflux_trans,
           pred=predict(cu_best_GWPco2andch4, type = "response"),
           strata=factor(cu_GWPco2andch4$strata, levels = c("bare","vegetated", "open water"))) %>%
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=strata),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  stat_cor(aes(col=strata),label.x.npc = 0.4, label.y.npc = 0.2)+
  geom_smooth(method = "lm", aes(col=strata), se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  scale_color_manual(values = c(
    "bare" = "#d95f02",      # brown/orange
    "vegetated" = "#2db82d", # green
    "open water" = "blue" # blue
  ), drop = FALSE) +
  labs(title = "CU_gwp Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

ggsave(filename = "CU_gwp_Observed_VS_predicted.png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)


data.frame(obs=cu_best_GWPco2andch4$frame$dailyflux_trans,
           pred=predict(cu_best_GWPco2andch4, type = "response"),
           strata=factor(cu_GWPco2andch4$strata, levels = c("bare","vegetated", "open water")),
           subsite=cu_best_GWPco2andch4$frame$subsite) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=subsite),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "CU_gwp Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

# ggsave(filename = "CU_gwp_Observed_VS_predicted(subsite-color).png", 
#        width = 115,height = 120,units = "mm",device = "png",dpi = 400,
#        path = extraplots_path
# )


#Remove non-best models
rm(cu_glmm_gaus0,cu_glmm_gaus1,cu_glmm_gaus2,
   cu_glmm_stu0,cu_glmm_stu1,cu_glmm_stu2,
   resgaus0,resgaus1,resgaus2,
   restu0, restu1, restu2)



## 1.3. DA model (ok)----
#approach: compare fit between gausian models, without heterocedasticity, with for season only and with for season and status. 

da_GWPco2andch4<- data4models %>% filter(casepilot=="DA"&ghgspecies=="GWPco2andch4")

#See effect of transformation: 
table_trans %>%
  filter(casepilot=="DA"&
           ghgspecies=="GWPco2andch4") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#FIT MODEL OPTIONS: 
da_glmm_gaus0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = da_GWPco2andch4,
                        family = gaussian(),
                        dispformula = ~1)

da_glmm_gaus1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = da_GWPco2andch4,
                        family = gaussian(),
                        dispformula = ~season)

da_glmm_gaus2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = da_GWPco2andch4,
                        family = gaussian(),
                        dispformula = ~season+status)

da_glmm_stu0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = da_GWPco2andch4,
                       family = t_family,
                       dispformula = ~1) 

da_glmm_stu1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = da_GWPco2andch4,
                       family = t_family,
                       dispformula = ~season)

da_glmm_stu2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = da_GWPco2andch4,
                       family = t_family,
                       dispformula = ~season+status)

#CALCULATE RESIDUALS
resgaus0<- simulateResiduals(da_glmm_gaus0)
resgaus1<- simulateResiduals(da_glmm_gaus1)
resgaus2<- simulateResiduals(da_glmm_gaus2)

restu0<- simulateResiduals(da_glmm_stu0)
restu1<- simulateResiduals(da_glmm_stu1)
restu2<- simulateResiduals(da_glmm_stu2)

#Check gaussian model
plot(resgaus0) #check Q-Q, residuals vs fitted, 
plotResiduals(resgaus0, da_GWPco2andch4$season) #check heterocedasticity
plotResiduals(resgaus0, da_GWPco2andch4$status) #check heterocedasticity
#Assumptions fail

plot(resgaus1) #check Q-Q, residuals vs fitted, 
plotResiduals(resgaus1, da_GWPco2andch4$season) #check heterocedasticity
plotResiduals(resgaus1, da_GWPco2andch4$status) #check heterocedasticity
#Assumptions fail

plot(resgaus2) #check Q-Q, residuals vs fitted, 
plotResiduals(resgaus2, da_GWPco2andch4$season) #check heterocedasticity
plotResiduals(resgaus2, da_GWPco2andch4$status) #check heterocedasticity
#Assumptions fail

#Check t_family model
plot(restu0) #check Q-Q, residuals vs fitted,  
plotResiduals(restu0, da_GWPco2andch4$season) #check heterocedasticity
plotResiduals(restu0, da_GWPco2andch4$status) #check heterocedasticity
#Heterocedasticity + funnel resid

plot(restu1) #check Q-Q, residuals vs fitted, 
plotResiduals(restu1, da_GWPco2andch4$season) #check heterocedasticity
plotResiduals(restu1, da_GWPco2andch4$status) #check heterocedasticity
#ALL pass

plot(restu2) #check Q-Q, residuals vs fitted, 
plotResiduals(restu2, da_GWPco2andch4$season) #check heterocedasticity
plotResiduals(restu2, da_GWPco2andch4$status) #check heterocedasticity
#ALL PASS

#GAUssian comparison: dispformula options
anova(da_glmm_gaus0,da_glmm_gaus1,da_glmm_gaus2) #none pass

#t_family comparison: dispformula options
anova(da_glmm_stu0,da_glmm_stu1,da_glmm_stu2)

#best vs best
anova(da_glmm_gaus1,da_glmm_stu1) #t_family is the best.


#DA best model: T-family, allowing heterocedasticity across seasons
da_best_GWPco2andch4<- da_glmm_stu1
Anova(da_best_GWPco2andch4,type = 3) 

#R2 of homocedastic model
da_best_homoc_GWPco2andch4<- da_glmm_stu0
r2(da_best_homoc_GWPco2andch4) #homocedastic t-family is singular....

#Singularity?
performance::check_singularity(da_best_GWPco2andch4) #Best-model is also singular

#plot obs vs predicted
data.frame(obs=da_best_GWPco2andch4$frame$dailyflux_trans,
           pred=predict(da_best_GWPco2andch4, type = "response"),
           strata=factor(da_GWPco2andch4$strata, levels = c("bare","vegetated", "open water"))) %>%
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=strata),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  stat_cor(aes(col=strata),label.x.npc = 0.4, label.y.npc = 0.2)+
  geom_smooth(method = "lm", aes(col=strata), se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  scale_color_manual(values = c(
    "bare" = "#d95f02",      # brown/orange
    "vegetated" = "#2db82d", # green
    "open water" = "blue" # blue
  ), drop = FALSE) +
  labs(title = "DA_gwp Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

ggsave(filename = "DA_gwp_Observed_VS_predicted.png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)


data.frame(obs=da_best_GWPco2andch4$frame$dailyflux_trans,
           pred=predict(da_best_GWPco2andch4, type = "response"),
           strata=factor(da_GWPco2andch4$strata, levels = c("bare","vegetated", "open water")),
           subsite=da_best_GWPco2andch4$frame$subsite) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=subsite),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "DA_gwp Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

# ggsave(filename = "DA_gwp_Observed_VS_predicted(subsite-color).png", 
#        width = 115,height = 120,units = "mm",device = "png",dpi = 400,
#        path = extraplots_path
# )


#Remove non-best models
rm(da_glmm_gaus0,da_glmm_gaus1,da_glmm_gaus2,
   da_glmm_stu0,da_glmm_stu1,da_glmm_stu2,
   resgaus0,resgaus1,resgaus2,
   restu0, restu1, restu2)



## 1.4. DU model (ok) ----
#approach: compare fit between gausian models: without heterocedasticity, with for season only and with for season and status. 

du_GWPco2andch4<- data4models %>% filter(casepilot=="DU"&ghgspecies=="GWPco2andch4")

#See effect of transformation: 
table_trans %>%
  filter(casepilot=="DU"&
           ghgspecies=="GWPco2andch4") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#FIT MODEL OPTIONS: 
du_glmm_gaus0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = du_GWPco2andch4,
                        family = gaussian(),
                        dispformula = ~1)

du_glmm_gaus1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = du_GWPco2andch4,
                        family = gaussian(),
                        dispformula = ~season)

du_glmm_gaus2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = du_GWPco2andch4,
                        family = gaussian(),
                        dispformula = ~season+status)

du_glmm_stu0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = du_GWPco2andch4,
                       family = t_family,
                       dispformula = ~1) 

du_glmm_stu1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = du_GWPco2andch4,
                       family = t_family,
                       dispformula = ~season)

du_glmm_stu2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = du_GWPco2andch4,
                       family = t_family,
                       dispformula = ~season+status)

#CALCULATE RESIDUALS
resgaus0<- simulateResiduals(du_glmm_gaus0)
resgaus1<- simulateResiduals(du_glmm_gaus1)
resgaus2<- simulateResiduals(du_glmm_gaus2)

restu0<- simulateResiduals(du_glmm_stu0)
restu1<- simulateResiduals(du_glmm_stu1)
restu2<- simulateResiduals(du_glmm_stu2)

#Check gaussian model
plot(resgaus0) #check Q-Q, residuals vs fitted, 
plotResiduals(resgaus0, du_GWPco2andch4$season) #check heterocedasticity
plotResiduals(resgaus0, du_GWPco2andch4$status) #check heterocedasticity
#Assumptions fail

plot(resgaus1) #check Q-Q, residuals vs fitted, 
plotResiduals(resgaus1, du_GWPco2andch4$season) #check heterocedasticity
plotResiduals(resgaus1, du_GWPco2andch4$status) #check heterocedasticity
#Assumptions fail

plot(resgaus2) #check Q-Q, residuals vs fitted, 
plotResiduals(resgaus2, du_GWPco2andch4$season) #check heterocedasticity
plotResiduals(resgaus2, du_GWPco2andch4$status) #check heterocedasticity
#Assumptions fail

#Check t_family model
check_convergence(du_glmm_stu0)
plot(restu0) #check Q-Q, residuals vs fitted,  
plotResiduals(restu0, du_GWPco2andch4$season) #check heterocedasticity
plotResiduals(restu0, du_GWPco2andch4$status) #check heterocedasticity
#DOES NOT CONVERGE

check_convergence(du_glmm_stu1)
plot(restu1) #check Q-Q, residuals vs fitted, 
plotResiduals(restu1, du_GWPco2andch4$season) #check heterocedasticity
plotResiduals(restu1, du_GWPco2andch4$status) #check heterocedasticity
#DOES NOT CONVERGE

check_convergence(du_glmm_stu2)
plot(restu2) #check Q-Q, residuals vs fitted, 
plotResiduals(restu2, du_GWPco2andch4$season) #check heterocedasticity
plotResiduals(restu2, du_GWPco2andch4$status) #check heterocedasticity
#DOES NOT CONVERGE

#GAUssian comparison: dispformula options
anova(du_glmm_gaus0,du_glmm_gaus1,du_glmm_gaus2)

#t_family comparison: dispformula options
anova(du_glmm_stu0,du_glmm_stu1,du_glmm_stu2) #NONE CONVERGE

#best vs best
#anova(du_glmm_gaus1,du_glmm_stu1) # only gaussian converge


#DU best model: GAUSSIAN, allowing heterocedasticity across seasons
du_best_GWPco2andch4<- du_glmm_gaus1
Anova(du_best_GWPco2andch4,type = 3) 

#R2 of homocedastic model
du_best_homoc_GWPco2andch4<- du_glmm_gaus0
r2(du_best_homoc_GWPco2andch4)

#Singularity?
performance::check_singularity(du_best_GWPco2andch4)


#plot obs vs predicted
data.frame(obs=du_best_GWPco2andch4$frame$dailyflux_trans,
           pred=predict(du_best_GWPco2andch4, type = "response"),
           strata=factor(du_GWPco2andch4$strata, levels = c("bare","vegetated", "open water"))) %>%
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=strata),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  stat_cor(aes(col=strata),label.x.npc = 0.4, label.y.npc = 0.2)+
  geom_smooth(method = "lm", aes(col=strata), se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  scale_color_manual(values = c(
    "bare" = "#d95f02",      # brown/orange
    "vegetated" = "#2db82d", # green
    "open water" = "blue" # blue
  ), drop = FALSE) +
  labs(title = "DA_gwp Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

ggsave(filename = "DU_gwp_Observed_VS_predicted.png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)


data.frame(obs=du_best_GWPco2andch4$frame$dailyflux_trans,
           pred=predict(du_best_GWPco2andch4, type = "response"),
           strata=factor(du_GWPco2andch4$strata, levels = c("bare","vegetated", "open water")),
           subsite=du_best_GWPco2andch4$frame$subsite) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=subsite),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "DU_gwp Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

# ggsave(filename = "DU_gwp_Observed_VS_predicted(subsite-color).png", 
#        width = 115,height = 120,units = "mm",device = "png",dpi = 400,
#        path = extraplots_path
# )

#Remove non-best models
rm(du_glmm_gaus0,du_glmm_gaus1,du_glmm_gaus2,
   du_glmm_stu0,du_glmm_stu1,du_glmm_stu2,
   resgaus0,resgaus1,resgaus2,
   restu0, restu1, restu2)


## 1.5. RI model (ok)----
#approach: compare fit between gausian models: without heterocedasticity, with for season only and with for season and status. 

ri_GWPco2andch4<- data4models %>% filter(casepilot=="RI"&ghgspecies=="GWPco2andch4")

#See effect of transformation: 
table_trans %>%
  filter(casepilot=="RI"&
           ghgspecies=="GWPco2andch4") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#FIT MODEL OPTIONS: 
ri_glmm_gaus0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = ri_GWPco2andch4,
                        family = gaussian(),
                        dispformula = ~1)

ri_glmm_gaus1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = ri_GWPco2andch4,
                        family = gaussian(),
                        dispformula = ~season)

ri_glmm_gaus2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = ri_GWPco2andch4,
                        family = gaussian(),
                        dispformula = ~season+status)

ri_glmm_stu0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = ri_GWPco2andch4,
                       family = t_family,
                       dispformula = ~1) 

ri_glmm_stu1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = ri_GWPco2andch4,
                       family = t_family,
                       dispformula = ~season)

ri_glmm_stu2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = ri_GWPco2andch4,
                       family = t_family,
                       dispformula = ~season+status)

#CALCULATE RESIDUALS
resgaus0<- simulateResiduals(ri_glmm_gaus0)
resgaus1<- simulateResiduals(ri_glmm_gaus1)
resgaus2<- simulateResiduals(ri_glmm_gaus2)

restu0<- simulateResiduals(ri_glmm_stu0)
restu1<- simulateResiduals(ri_glmm_stu1)
restu2<- simulateResiduals(ri_glmm_stu2)

#Check gaussian model
plot(resgaus0) #check Q-Q, residuals vs fitted, 
plotResiduals(resgaus0, ri_GWPco2andch4$season) #check heterocedasticity
plotResiduals(resgaus0, ri_GWPco2andch4$status) #check heterocedasticity
#Heterocedasticity status

plot(resgaus1) #check Q-Q, residuals vs fitted, 
plotResiduals(resgaus1, ri_GWPco2andch4$season) #check heterocedasticity
plotResiduals(resgaus1, ri_GWPco2andch4$status) #check heterocedasticity
#Heterocedasticity

plot(resgaus2) #check Q-Q, residuals vs fitted, 
plotResiduals(resgaus2, ri_GWPco2andch4$season) #check heterocedasticity
plotResiduals(resgaus2, ri_GWPco2andch4$status) #check heterocedasticity
#ALL PASS

#Check t_family model
check_convergence(ri_glmm_stu0)
plot(restu0) #check Q-Q, residuals vs fitted,  
plotResiduals(restu0, ri_GWPco2andch4$season) #check heterocedasticity
plotResiduals(restu0, ri_GWPco2andch4$status) #check heterocedasticity
#Heterocedasticity

check_convergence(ri_glmm_stu1)
plot(restu1) #check Q-Q, residuals vs fitted, 
plotResiduals(restu1, ri_GWPco2andch4$season) #check heterocedasticity
plotResiduals(restu1, ri_GWPco2andch4$status) #check heterocedasticity
#Heterocedasticity

check_convergence(ri_glmm_stu2)
plot(restu2) #check Q-Q, residuals vs fitted, 
plotResiduals(restu2, ri_GWPco2andch4$season) #check heterocedasticity
plotResiduals(restu2, ri_GWPco2andch4$status) #check heterocedasticity
#ALL PASS

#GAUssian comparison: dispformula options
anova(ri_glmm_gaus0,ri_glmm_gaus1,ri_glmm_gaus2)

#t_family comparison: dispformula options
anova(ri_glmm_stu0,ri_glmm_stu1,ri_glmm_stu2) 

#best vs best
anova(ri_glmm_gaus0,ri_glmm_stu2) 


#RI best model: GAUSSIAN, Homocedastic
ri_best_GWPco2andch4<- ri_glmm_gaus0
Anova(ri_best_GWPco2andch4,type = 3) 

#R2 of homocedastic model
ri_best_homoc_GWPco2andch4<- ri_glmm_gaus0
r2(ri_best_homoc_GWPco2andch4)

#Singularity?
performance::check_singularity(ri_best_GWPco2andch4)

#plot obs vs predicted
data.frame(obs=ri_best_GWPco2andch4$frame$dailyflux_trans,
           pred=predict(ri_best_GWPco2andch4, type = "response"),
           strata=factor(ri_GWPco2andch4$strata, levels = c("bare","vegetated", "open water"))) %>%
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=strata),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  stat_cor(aes(col=strata),label.x.npc = 0.4, label.y.npc = 0.2)+
  geom_smooth(method = "lm", aes(col=strata), se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  scale_color_manual(values = c(
    "bare" = "#d95f02",      # brown/orange
    "vegetated" = "#2db82d", # green
    "open water" = "blue" # blue
  ), drop = FALSE) +
  labs(title = "RI_gwp Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

ggsave(filename = "RI_gwp_Observed_VS_predicted.png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)


data.frame(obs=ri_best_GWPco2andch4$frame$dailyflux_trans,
           pred=predict(ri_best_GWPco2andch4, type = "response"),
           strata=factor(ri_GWPco2andch4$strata, levels = c("bare","vegetated", "open water")),
           subsite=ri_best_GWPco2andch4$frame$subsite) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=subsite),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "RI_gwp Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

# ggsave(filename = "RI_gwp_Observed_VS_predicted(subsite-color).png", 
#        width = 115,height = 120,units = "mm",device = "png",dpi = 400,
#        path = extraplots_path
# )


#Remove non-best models
rm(ri_glmm_gaus0,ri_glmm_gaus1,ri_glmm_gaus2,
   ri_glmm_stu0,ri_glmm_stu1,ri_glmm_stu2,
   resgaus0,resgaus1,resgaus2,
   restu0, restu1, restu2)



## 1.6. VA model (ok, bad) ----
#approach: compare fit between gausian and tstudent models without heterocedasticity, with for season only and with for season and status (not for interaction). 

va_GWPco2andch4<- data4models %>% filter(casepilot=="VA"&ghgspecies=="GWPco2andch4")

#See effect of transformation: 
table_trans %>%
  filter(casepilot=="VA"&
           ghgspecies=="GWPco2andch4") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#FIT MODEL OPTIONS: 
va_glmm_gaus0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = va_GWPco2andch4,
                        family = gaussian(),
                        dispformula = ~1)

va_glmm_gaus1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = va_GWPco2andch4,
                        family = gaussian(),
                        dispformula = ~season)

va_glmm_gaus2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = va_GWPco2andch4,
                        family = gaussian(),
                        dispformula = ~season+status)

va_glmm_stu0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = va_GWPco2andch4,
                       family = t_family,
                       dispformula = ~1) 

va_glmm_stu1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = va_GWPco2andch4,
                       family = t_family,
                       dispformula = ~season)

va_glmm_stu2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = va_GWPco2andch4,
                       family = t_family,
                       dispformula = ~season+status)

#CALCULATE RESIDUALS
resgaus0<- simulateResiduals(va_glmm_gaus0)
resgaus1<- simulateResiduals(va_glmm_gaus1)
resgaus2<- simulateResiduals(va_glmm_gaus2)

restu0<- simulateResiduals(va_glmm_stu0)
restu1<- simulateResiduals(va_glmm_stu1)
restu2<- simulateResiduals(va_glmm_stu2)

#Check gaussian model
plot(resgaus0) #check Q-Q, residuals vs fitted, 
plotResiduals(resgaus0, va_GWPco2andch4$season) #check heterocedasticity
plotResiduals(resgaus0, va_GWPco2andch4$status) #check heterocedasticity
#Assumptions fail

plot(resgaus1) #check Q-Q, residuals vs fitted, 
plotResiduals(resgaus1, va_GWPco2andch4$season) #check heterocedasticity
plotResiduals(resgaus1, va_GWPco2andch4$status) #check heterocedasticity
#Heterocedasticity

plot(resgaus2) #check Q-Q, residuals vs fitted, 
plotResiduals(resgaus2, va_GWPco2andch4$season) #check heterocedasticity
plotResiduals(resgaus2, va_GWPco2andch4$status) #check heterocedasticity
#Heterocedasticity

#Check t_family model
check_convergence(va_glmm_stu0) #DOES NOT CONVEGE
plot(restu0) #check Q-Q, residuals vs fitted,  
plotResiduals(restu0, va_GWPco2andch4$season) #check heterocedasticity
plotResiduals(restu0, va_GWPco2andch4$status) #check heterocedasticity
#DOES NOT CONVEGE

check_convergence(va_glmm_stu1)
plot(restu1) #check Q-Q, residuals vs fitted, 
plotResiduals(restu1, va_GWPco2andch4$season) #check heterocedasticity
plotResiduals(restu1, va_GWPco2andch4$status) #check heterocedasticity
#Heterocedasticity

check_convergence(va_glmm_stu2)
plot(restu2) #check Q-Q, residuals vs fitted, 
plotResiduals(restu2, va_GWPco2andch4$season) #check heterocedasticity
plotResiduals(restu2, va_GWPco2andch4$status) #check heterocedasticity
#Heterocedasticity

#GAUssian comparison: dispformula options
anova(va_glmm_gaus0,va_glmm_gaus1,va_glmm_gaus2)

#t_family comparison: dispformula options
anova(va_glmm_stu0,va_glmm_stu1,va_glmm_stu2) 

#best vs best
anova(va_glmm_gaus1,va_glmm_stu1) 


#VA best model: GAUSSIAN, Heterocedasticity for season 
va_best_GWPco2andch4<- va_glmm_gaus1
Anova(va_best_GWPco2andch4,type = 3) 

#R2 of homocedastic model
va_best_homoc_GWPco2andch4<- va_glmm_gaus0
r2(va_best_homoc_GWPco2andch4)

#Singularity?
performance::check_singularity(va_best_GWPco2andch4)


#plot obs vs predicted
data.frame(obs=va_best_GWPco2andch4$frame$dailyflux_trans,
           pred=predict(va_best_GWPco2andch4, type = "response"),
           strata=factor(va_GWPco2andch4$strata, levels = c("bare","vegetated", "open water"))) %>%
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=strata),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  stat_cor(aes(col=strata),label.x.npc = 0.4, label.y.npc = 0.2)+
  geom_smooth(method = "lm", aes(col=strata), se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  scale_color_manual(values = c(
    "bare" = "#d95f02",      # brown/orange
    "vegetated" = "#2db82d", # green
    "open water" = "blue" # blue
  ), drop = FALSE) +
  labs(title = "VA_gwp Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

ggsave(filename = "VA_gwp_Observed_VS_predicted.png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)


data.frame(obs=va_best_GWPco2andch4$frame$dailyflux_trans,
           pred=predict(va_best_GWPco2andch4, type = "response"),
           strata=factor(va_GWPco2andch4$strata, levels = c("bare","vegetated", "open water")),
           subsite=va_best_GWPco2andch4$frame$subsite) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=subsite),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "VA_gwp Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

# ggsave(filename = "VA_gwp_Observed_VS_predicted(subsite-color).png", 
#        width = 115,height = 120,units = "mm",device = "png",dpi = 400,
#        path = extraplots_path
# )


#Remove non-best models
rm(va_glmm_gaus0,va_glmm_gaus1,va_glmm_gaus2,
   va_glmm_stu0,va_glmm_stu1,va_glmm_stu2,
   resgaus0,resgaus1,resgaus2,
   restu0, restu1, restu2)


#2. Save bestmodel-list -------

#Save lists with best-models and homocedastic models, named as casepilot_ghgspecies (needed in order to work appropriately with the summarising functions)

#Provide list of homocedastic model variants (to estimate R2m and R2c) Cannot be computed for heterocedastic models or singular models
homoced_model_list_GWPco2andch4<- list(
  "CA_GWPco2andch4" = ca_best_homoc_GWPco2andch4,
  "CU_GWPco2andch4" = cu_best_homoc_GWPco2andch4, 
  "DA_GWPco2andch4" = da_best_homoc_GWPco2andch4,
  "DU_GWPco2andch4" = du_best_homoc_GWPco2andch4,
  "RI_GWPco2andch4" = ri_best_homoc_GWPco2andch4, 
  "VA_GWPco2andch4" = va_best_homoc_GWPco2andch4)

#List of best models for each casepilot (to get significances, pseudoR2s and residual tests)
best_model_list_GWPco2andch4<- list(
  "CA_GWPco2andch4" = ca_best_GWPco2andch4,
  "CU_GWPco2andch4" = cu_best_GWPco2andch4, 
  "DA_GWPco2andch4" = da_best_GWPco2andch4,
  "DU_GWPco2andch4" = du_best_GWPco2andch4,
  "RI_GWPco2andch4" = ri_best_GWPco2andch4,
  "VA_GWPco2andch4" = va_best_GWPco2andch4)


#______________________-----
#SUMMARY across GHGs-------

#Combine model-lists
all_best_models<- c(best_model_list_co2,best_model_list_ch4, best_model_list_GWPco2andch4)

all_homoced_models<- c(homoced_model_list_co2,homoced_model_list_ch4, homoced_model_list_GWPco2andch4)
  
  ## 2.1. Overall fit -----
#Structure, pseudoR2 of best model:
all_best_model_fit<- get_pseudoR2s(all_best_models)
#TWO warnings of "false convergence"



#R2c and R2m of homocedastic model:
all_homoc_r2<- get_homoc_R2s(all_homoced_models)
#2 Warnings:  homocedastic models with singularity: CU_co2 and DA_gwp


#DHARMa Residual diagnostics of best model:
all_resid_diag_summary <- summarize_dharma_diagnostics(all_best_models)

print(all_resid_diag_summary)


##2.2. Significance of effects-----
#Effect of status will first be assessed via Anova (best_model,type=3): is there an effect averaging across all seasons? yes/no, is the effect dependent on season? 

#Obtain significance for main effects of best_models
all_results_anova<-get_anova_results(all_best_models) 
print(all_results_anova)

#Add symbol on top of p_value

#Function to get symbols
pval_to_symbol <- function(p) {
  ifelse(p < 0.001, "***",
         ifelse(p < 0.01, "**",
                ifelse(p < 0.05, "*",
                       ifelse(p < 0.1, ".", "ns")
                )
         )
  )
}

#Add symbols from significance pvalues
all_results_anova<- all_results_anova %>% 
  mutate(across(c(intercept_pval,status_pval,season_pval,interaction_pval),pval_to_symbol,.names = "symbol_{.col}"))


#Intercept refers to whether the grand mean (on the transformed scale) is significantly different than cero (Is the overall net exchange of CH4 balanced between uptake and emission)
#Status refers to whether the status (averaged across seasons) has a significant effect on flux. 
#Season refers to whether the season (averaged across status) has a significant effect on flux.
#Interaction refers to whether the effect of status depends on season — i.e., the difference between status levels changes across seasons.


all_model_outputs<- all_results_anova %>% 
  full_join(all_best_model_fit) %>% 
  full_join(all_homoc_r2 %>% dplyr::select(dataset, homoced_R2m, homoced_R2c))


#Look for significance at either status or status:season---> there is an effect (that may or may not be evident when averaging across seasons)
#IF status:season significant---> effect depends on season. 

#Does status have an effect (significance of status or interaction)
all_results_anova %>% filter(status_pval<0.05|interaction_pval<0.05)


#Save in dropbox tables with the models summary outputs: 

write.csv(x = all_model_outputs, file = paste0(paper_path,"Allghg_summary_chambermodels.csv"),row.names = F)

write.csv(x = all_resid_diag_summary, file = paste0(paper_path,"Allghg_Residualtests_chambermodels.csv"),row.names = F)


##2.3. Effect direction and magnitude----

#Use model scale for statistics, back-transform for interpretation.  

#NOTE on bootstrapping: there is the possibility to bootstrap the models (i.e. recalculate the model many times with re-sampled data, and calculate emmeans for each model, to be able to produce more robust Confidence Intervals for them, but this is not usually done and computationally intensive. We would have to take into account the model structure when re-sampling the data to ensure balanced resamplings (all factors represented by data). WE will not do this (for the moment).

#Some guides to plot and use emmeans: https://rcompanion.org/handbook/G_06.html


#Then emmeans comparisons to get direction of effect and magnitude (averaging over season interaction if any)
plot(emmeans(ri_best_GWPco2andch4, ~status), comparisons=TRUE) #Unclear interpretation of arrows, not the best way of presenting it. 


#Notes on non-gaussian models: upper and lower confidence limits are given as asymp.LCL/asymp.UCL (in our context we can treat them equally as the ones from gaussian models). Additionally, non-gaussian models will have Inf df (degrees of freedom are calculated differently for T-family models), not an issue. 

# Initialize list for storing comparisons
comparisons_list <- list()

# Loop to extract all comparisons for every casepilot best model
for (dataset in names(all_best_models)) {
  cp_model <- all_best_models[[dataset]]
  
  #extract casepilot name and ghgspecies from model-list names
  casepilot_name<- sub("_.*", "", dataset)
  ghgspecies<- sub(paste0(casepilot_name,"_"),"",dataset)
  
  # Status comparisons:
  status_emmeans <- emmeans(cp_model, ~status)
  status_comparisons <- cld(status_emmeans, 
                            alpha = 0.05,
                            Letters = letters,
                            adjust = "sidak") %>%
    mutate(comparison = "status",
           ghgspecies = ghgspecies,
           casepilot = casepilot_name,
           season = NA) %>%
    rename(lower.CL = any_of("asymp.LCL"),
           upper.CL = any_of("asymp.UCL"))
  
  # Season comparisons: 
  season_emmeans <- emmeans(cp_model, ~season)
  season_comparisons <- cld(season_emmeans, 
                            alpha = 0.05,
                            Letters = letters,
                            adjust = "sidak") %>%
    mutate(comparison = "season",
           ghgspecies = ghgspecies,
           casepilot = casepilot_name,
           status = NA) %>%
    rename(lower.CL = any_of("asymp.LCL"),
           upper.CL = any_of("asymp.UCL"))
  
  # Interaction comparisons: 
  interaction_emmeans <- emmeans(cp_model, ~status:season)
  interaction_comparisons <- cld(interaction_emmeans, 
                                 alpha = 0.05,
                                 Letters = letters,
                                 adjust = "sidak") %>%
    mutate(comparison = "interaction",
           ghgspecies = ghgspecies,
           casepilot = casepilot_name) %>%
    rename(lower.CL = any_of("asymp.LCL"),
           upper.CL = any_of("asymp.UCL"))
  
  # Join and format CP comparisons
  all_comparisons <- status_comparisons %>%
    full_join(season_comparisons) %>%
    full_join(interaction_comparisons) %>%
    rename(group_letter = .group) %>%
    mutate(
      comparison = factor(comparison, levels = c("status", "season", "interaction")),
      season = factor(season, levels = c("S1", "S2", "S3", "S4")),
      status = factor(status, levels = c("Altered", "Preserved", "Restored")),
      group_letter = gsub(" ", "", group_letter),
      
      # Create key for accessing BestNormalize object
      key_trans = paste(casepilot_name, ghgspecies, sep = "_")
    ) %>%
    rowwise() %>% 
    # BACK-TRANSFORMATION of emmean, SE and CL using BestNormalize object
    mutate(
      emmean_bt = if (!is.null(bn_list[[key_trans]])) predict(bn_list[[key_trans]], emmean, inverse = TRUE) else NA_real_,
      SE_bt = if (!is.null(bn_list[[key_trans]])) predict(bn_list[[key_trans]], SE, inverse = TRUE) else NA_real_,
      lower.CL_bt = if (!is.null(bn_list[[key_trans]])) predict(bn_list[[key_trans]], lower.CL, inverse = TRUE) else NA_real_,
      upper.CL_bt = if (!is.null(bn_list[[key_trans]])) predict(bn_list[[key_trans]], upper.CL, inverse = TRUE) else NA_real_
    ) %>%
    
    dplyr::select(any_of(c(
      "ghgspecies", "casepilot","comparison", "status", "season", "df",
      "emmean", "SE", "lower.CL", "upper.CL", "group_letter",
      "emmean_bt", "SE_bt", "lower.CL_bt", "upper.CL_bt"
    ))) %>%
    arrange(ghgspecies, casepilot, comparison, status, season)
  
  # Store casepilot results
  comparisons_list[[dataset]] <- all_comparisons
}

#Remove within-loop objects
rm(all_comparisons, status_comparisons,season_comparisons,interaction_comparisons,
   interaction_emmeans,status_emmeans, season_emmeans)


# Unir todos los resultados en un solo data.frame
all_emmeans_comparisons <- bind_rows(comparisons_list)


#Save all comparisons:  
write.csv(x = all_emmeans_comparisons, file = paste0(paper_path,"Allghg_emmeans-posthoc_chambermodels.csv"),row.names = F)


# save.image(paste0(plots_path,"/model_results_last.Rdata"))

#____________________--------

#INSPECT Bad models-------
#Keeping all strata together (not differentiating) and subsites with specific behaviors under the same status can lead to model issues and lack of detection of restoration effects. 
#Here we inspect the performance of the models to identify "bad" models.

fit_best<- read.csv(paste0(paper_path,"Allghg_summary_chambermodels.csv"))
resid_best<- read.csv(paste0(paper_path,"Allghg_Residualtests_chambermodels.csv"))

#Convergence issue: 
convergence_flag<- resid_best %>% 
  filter(!convergence) %>% 
  mutate(flag="Model does not converge") %>% 
  dplyr::select(dataset,flag)

#Singular models: 
singular_flag<- resid_best %>% 
  filter(singular) %>%
  mutate(flag="Singular model") %>% 
  dplyr::select(dataset,flag)

#Uniformity fail: 
uniformity_flag<-resid_best %>% 
  filter(uniformity_pval<0.05) %>%
  mutate(flag="Uniformity fail (p<0.05)") %>% 
  dplyr::select(dataset,flag)

#Dispersion fail: 
dispersion_flag<-resid_best %>% 
  filter(dispersion_pval<0.05) %>%
  mutate(flag="Dispersion fail (p<0.05)") %>% 
  dplyr::select(dataset,flag)
         
#VERY LOW pseudoR2: obsVSpred
low_predict_flag<-fit_best %>%
  filter(pseudoR2_Efron<0.15) %>% 
  mutate(flag="Low Obs vs pred power (pseudoR2_Efron <0.15)") %>% 
  dplyr::select(dataset,flag)

#Great improvement R2c vs R2m
r2c_r2m_flag<- fit_best %>% 
  mutate(rel_improve=homoced_R2c/homoced_R2m) %>% 
  filter(rel_improve>1.3) %>% 
  mutate(flag="Big rel. improvement (R2c/R2m)>1.2") %>% 
  dplyr::select(dataset,flag)

all_flags<- convergence_flag%>% 
  bind_rows(singular_flag, dispersion_flag,uniformity_flag,low_predict_flag, r2c_r2m_flag)

all_flags

unique(all_flags$dataset)


#TESTS including strata as factor--------

#We need to have all combinations of factor levels in order to create a model that accounts for these 3 factors and that can say if any of these factors is signifcant in determining the fluxes. 






#____________________--------

#STRATA-Specific MODELS-------



#Here we test the approach of modelling different strata independently (when appropriate).
#Rationale: different strata can have different emission profiles and respond differently to alterations, so that modelling each independently can: (I) improve overall explanatory power AND (II) provide strata-specific effects (and magnitudes) for different status that will help with management and biogeochemical interpretation.


#CU: ow model, vegetation model, for co2 and ch4 OK!
#RI: overall model is enough, we cannot make strata-specific models

#Other casepilots


#Strata-specific transformations-----
#RE-check bestNormalization: separating data for each ghg*casepilot*strata combination. We will only use some of this, as others are not appropriate (e.g. RI cannot be separated into strata, as alteration and strata composition are strictly linked and are not fully shared, i.e. no bare in restored for example)

#USE bestNormalize (allowing for pseudo-log trans) to obtain for each dataset the most-Normal transformation: 


rm(table_trans, table_trans_i)
for (cp in unique(data4models$casepilot)) {
  for (ghg in  unique(data4models$ghgspecies)) {
    for (stratum in unique(data4models$strata)){
    x<- data4models %>% filter(casepilot==cp&ghgspecies==ghg&stratum==strata) %>%
      pull(dailyflux)
    
    #Check if there is enough data:
    if (length(x)>10) {
      #If more than 10 values for a particular dataset(ghg*cp*strata), perform transformation, 
    bn_result<- bestNormalize(x = x, new_transforms = pseudolog_transform, standardize = TRUE, warn = TRUE, k = 5,allow_orderNorm = F)
    #calculate normality test (shapiro)
    x_trans <- bn_result$chosen_transform$x.t
    pval <- shapiro.test(x_trans)$p.value
    
    table_trans_i <- tibble(
      dataset = paste(cp, ghg, stratum, sep = "_"),
      casepilot = cp,
      ghgspecies = ghg,
      strata = stratum,
      nobs= length(x),
      trans_name = attr(x = bn_result$chosen_transform,which = "class")[1],
      trans_Pearson_p=bn_result$chosen_transform$norm_stat, #Pearson's P (the lower, the more-normal data, only used for ranking our transformation options)
      shapiro_pval=pval, #P value of Normality
      bn_result = list(bn_result)  # wrap BestNormalize object in list to store as list-column
    )
    } else{
      #If less than 10 values for a particular dataset, skip transformation and create dummy row. 
      table_trans_i <- tibble(
        dataset = paste(cp, ghg, stratum, sep = "_"),
        casepilot = cp,
        ghgspecies = ghg,
        strata = stratum,
        nobs= length(x),
        trans_name = "Not enough data",
        trans_Pearson_p=as.numeric(NA_real_), 
        shapiro_pval=as.numeric(NA_real_), 
        bn_result = list(NA_real_))
      
    }
    #Create table_trans for first round of loop, then apend to table
    if(cp==unique(data4models$casepilot)[1]&
       ghg==unique(data4models$ghgspecies)[1]&
       stratum==unique(data4models$strata)[1]){
      table_trans<- table_trans_i
    }else{ table_trans<- bind_rows(table_trans, table_trans_i)}
    }
  }
}

#View all bestNormalizing transformations:
table_trans

rm(cp, bn_result, table_trans_i, ghg, x, x_trans, pval,stratum)

# Create a named list of BestNormalize objects for easy lookup
bn_list_strata <- table_trans %>%
  dplyr::select(dataset, bn_result) %>%
  deframe()

# Add a key column to data4models for matching
data4models_strata <- data4models %>%
  mutate(dataset = paste(casepilot, ghgspecies,strata, sep = "_"))

# Apply the corresponding BestNormalize transformation
data4models_strata <- data4models_strata %>%
  rowwise() %>%
  mutate(
    trans_name = table_trans$trans_name[match(dataset, paste(table_trans$casepilot, table_trans$ghgspecies,table_trans$strata, sep = "_"))],
    dailyflux_trans = if (!is.null(bn_list_strata[[dataset]])) predict(bn_list_strata[[dataset]], dailyflux) else NA_real_
  ) %>%
  ungroup()



#CURONIAN  ------

  #CU_co2_ow---- 

cu_co2_ow<- data4models_strata %>% filter(casepilot=="CU"&ghgspecies=="co2"&strata=="open water")

#See effect of transformation: 
table_trans %>%
  filter(casepilot=="CU",
         ghgspecies=="co2", strata=="open water") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#FIT MODEL OPTIONS: 
cu_glmm_gaus0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = cu_co2_ow,
                        family = gaussian(),
                        dispformula = ~1)

cu_glmm_gaus1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = cu_co2_ow,
                        family = gaussian(),
                        dispformula = ~season)

cu_glmm_gaus2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = cu_co2_ow,
                        family = gaussian(),
                        dispformula = ~season+status)

cu_glmm_stu0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = cu_co2_ow,
                       family = t_family,
                       dispformula = ~1) 

cu_glmm_stu1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = cu_co2_ow,
                       family = t_family,
                       dispformula = ~season)

cu_glmm_stu2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = cu_co2_ow,
                       family = t_family,
                       dispformula = ~season+status)

#CALCULATE RESIDUALS
resgaus0<- simulateResiduals(cu_glmm_gaus0)
resgaus1<- simulateResiduals(cu_glmm_gaus1)
resgaus2<- simulateResiduals(cu_glmm_gaus2)

restu0<- simulateResiduals(cu_glmm_stu0)
restu1<- simulateResiduals(cu_glmm_stu1)
restu2<- simulateResiduals(cu_glmm_stu2)


#Check gaussian model
plotQQunif(resgaus0) #check Q-Q
plotResiduals(resgaus0) #check residuals vs fitted, 
plotResiduals(resgaus0, cu_co2_ow$season) #check heterocedasticity
plotResiduals(resgaus0, cu_co2_ow$status) #check heterocedasticity
#season heterocedasticity

plotQQunif(resgaus1) #check Q-Q
plotResiduals(resgaus1) #check residuals vs fitted, 
plotResiduals(resgaus1, cu_co2_ow$season) #check heterocedasticity
plotResiduals(resgaus1, cu_co2_ow$status) #check heterocedasticity
#ALL PASS

plotQQunif(resgaus2) #check Q-Q
plotResiduals(resgaus2) #check residuals vs fitted, 
plotResiduals(resgaus2, cu_co2_ow$season) #check heterocedasticity
plotResiduals(resgaus2, cu_co2_ow$status) #check heterocedasticity
#ALL PASS

#Check t_family model
check_convergence(cu_glmm_stu0)
plot(restu0) #check Q-Q, residuals vs fitted,  
plotResiduals(restu0, cu_co2_ow$season) #check heterocedasticity
plotResiduals(restu0, cu_co2_ow$status) #check heterocedasticity
#season heterocedasticity

check_convergence(cu_glmm_stu1)
plot(restu1) #check Q-Q, residuals vs fitted, 
plotResiduals(restu1, cu_co2_ow$season) #check heterocedasticity
plotResiduals(restu1, cu_co2_ow$status) #check heterocedasticity
#ALL PASS

check_convergence(cu_glmm_stu2)
plot(restu2) #check Q-Q, residuals vs fitted, 
plotResiduals(restu2, cu_co2_ow$season) #check heterocedasticity
plotResiduals(restu2, cu_co2_ow$status) #check heterocedasticity
#ALL PASS

#GAUssian comparison: dispformula options
anova(cu_glmm_gaus0,cu_glmm_gaus1,cu_glmm_gaus2)

#t_family comparison: dispformula options
anova(cu_glmm_stu0,cu_glmm_stu1,cu_glmm_stu2)

#best vs best
anova(cu_glmm_gaus2,cu_glmm_stu1) #T_family, with season heterocedasticity


#CU_co2_ow best model?: #T_family, with season heterocedasticity
cu_co2_ow_best<- cu_glmm_stu1
Anova(cu_co2_ow_best,type = 3) 


#R2 of homocedastic model
cu_co2_ow_best_homoc<- cu_glmm_stu0
r2(cu_co2_ow_best_homoc)

#Singularity?
performance::check_singularity(cu_co2_ow_best)

#Outliers?
testOutliers(cu_co2_ow_best)

#plot obs vs predicted

data.frame(obs=cu_co2_ow_best$frame$dailyflux_trans,
           pred=predict(cu_co2_ow_best, type = "response"),
           strata=factor(cu_co2_ow$strata, levels = c("bare","vegetated", "open water")),
           status=factor(cu_co2_ow$status)) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=status),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "CU_co2_ow Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

ggsave(filename = "CU_co2_ow_Observed_VS_predicted.png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)



#Remove non-best models
rm(cu_glmm_gaus0,cu_glmm_gaus1,cu_glmm_gaus2,
   cu_glmm_stu0,cu_glmm_stu1,cu_glmm_stu2,
   resgaus0,resgaus1,resgaus2,
   restu0, restu1, restu2)




#CU_co2_v---- 

cu_co2_v<- data4models_strata %>% filter(casepilot=="CU"&ghgspecies=="co2"&strata=="vegetated")

#See effect of transformation: 
table_trans %>%
  filter(casepilot=="CU",
         ghgspecies=="co2", strata=="vegetated") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#FIT MODEL OPTIONS: 
cu_glmm_gaus0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = cu_co2_v,
                        family = gaussian(),
                        dispformula = ~1)

cu_glmm_gaus1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = cu_co2_v,
                        family = gaussian(),
                        dispformula = ~season)

cu_glmm_gaus2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = cu_co2_v,
                        family = gaussian(),
                        dispformula = ~season+status)

cu_glmm_stu0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = cu_co2_v,
                       family = t_family,
                       dispformula = ~1) 

cu_glmm_stu1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = cu_co2_v,
                       family = t_family,
                       dispformula = ~season)

cu_glmm_stu2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = cu_co2_v,
                       family = t_family,
                       dispformula = ~season+status)

#CALCULATE RESIDUALS
resgaus0<- simulateResiduals(cu_glmm_gaus0)
resgaus1<- simulateResiduals(cu_glmm_gaus1)
resgaus2<- simulateResiduals(cu_glmm_gaus2)

restu0<- simulateResiduals(cu_glmm_stu0)
restu1<- simulateResiduals(cu_glmm_stu1)
restu2<- simulateResiduals(cu_glmm_stu2)


#Check gaussian model
plotQQunif(resgaus0) #check Q-Q
plotResiduals(resgaus0) #check residuals vs fitted, 
plotResiduals(resgaus0, cu_co2_v$season) #check heterocedasticity
plotResiduals(resgaus0, cu_co2_v$status) #check heterocedasticity
#season heterocedasticity, residuals funnel

plotQQunif(resgaus1) #check Q-Q
plotResiduals(resgaus1) #check residuals vs fitted, 
plotResiduals(resgaus1, cu_co2_v$season) #check heterocedasticity
plotResiduals(resgaus1, cu_co2_v$status) #check heterocedasticity
#ALL PASS

plotQQunif(resgaus2) #check Q-Q
plotResiduals(resgaus2) #check residuals vs fitted, 
plotResiduals(resgaus2, cu_co2_v$season) #check heterocedasticity
plotResiduals(resgaus2, cu_co2_v$status) #check heterocedasticity
#ALL PASS

#Check t_family model
check_convergence(cu_glmm_stu0)
plot(restu0) #check Q-Q, residuals vs fitted,  
plotResiduals(restu0, cu_co2_v$season) #check heterocedasticity
plotResiduals(restu0, cu_co2_v$status) #check heterocedasticity
#Assumptions fail

check_convergence(cu_glmm_stu1)
plot(restu1) #check Q-Q, residuals vs fitted, 
plotResiduals(restu1, cu_co2_v$season) #check heterocedasticity
plotResiduals(restu1, cu_co2_v$status) #check heterocedasticity
#DOES NOT CONVERGE

check_convergence(cu_glmm_stu2)
plot(restu2) #check Q-Q, residuals vs fitted, 
plotResiduals(restu2, cu_co2_v$season) #check heterocedasticity
plotResiduals(restu2, cu_co2_v$status) #check heterocedasticity
#ALL PASS

#GAUssian comparison: dispformula options
anova(cu_glmm_gaus0,cu_glmm_gaus1,cu_glmm_gaus2)

#t_family comparison: dispformula options
# anova(cu_glmm_stu0,cu_glmm_stu1,cu_glmm_stu2)#OMIT, FAIL CONVERGENCE

#best vs best
# anova(cu_glmm_gaus1,cu_glmm_stu1) #OMIT, T_family not appropriate


#CU_co2_v best model?: #Gaussian, with season heterocedasticity
cu_co2_v_best<- cu_glmm_gaus1
Anova(cu_co2_v_best,type = 3) 


#R2 of homocedastic model
cu_co2_v_best_homoc<- cu_glmm_gaus0
r2(cu_co2_v_best_homoc)

#Singularity?
performance::check_singularity(cu_co2_v_best)
#Outliers?
testOutliers(cu_co2_ow_best)

#plot obs vs predicted

data.frame(obs=cu_co2_v_best$frame$dailyflux_trans,
           pred=predict(cu_co2_v_best, type = "response"),
           strata=factor(cu_co2_v$strata, levels = c("bare","vegetated", "open water")),
           status=factor(cu_co2_v$status)) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=status),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "CU_co2_v Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

ggsave(filename = "CU_co2_v_Observed_VS_predicted.png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)

#Remove non-best models
rm(cu_glmm_gaus0,cu_glmm_gaus1,cu_glmm_gaus2,
   cu_glmm_stu0,cu_glmm_stu1,cu_glmm_stu2,
   resgaus0,resgaus1,resgaus2,
   restu0, restu1, restu2)




#CU_ch4_ow---- 

cu_ch4_ow<- data4models_strata %>% filter(casepilot=="CU"&ghgspecies=="ch4"&strata=="open water")

#See effect of transformation: 
table_trans %>%
  filter(casepilot=="CU",
         ghgspecies=="ch4", strata=="open water") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#FIT MODEL OPTIONS: 
cu_glmm_gaus0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = cu_ch4_ow,
                        family = gaussian(),
                        dispformula = ~1)

cu_glmm_gaus1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = cu_ch4_ow,
                        family = gaussian(),
                        dispformula = ~season)

cu_glmm_gaus2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = cu_ch4_ow,
                        family = gaussian(),
                        dispformula = ~season+status)

cu_glmm_stu0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = cu_ch4_ow,
                       family = t_family,
                       dispformula = ~1) 

cu_glmm_stu1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = cu_ch4_ow,
                       family = t_family,
                       dispformula = ~season)

cu_glmm_stu2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = cu_ch4_ow,
                       family = t_family,
                       dispformula = ~season+status)

#CALCULATE RESIDUALS
resgaus0<- simulateResiduals(cu_glmm_gaus0)
resgaus1<- simulateResiduals(cu_glmm_gaus1)
resgaus2<- simulateResiduals(cu_glmm_gaus2)

restu0<- simulateResiduals(cu_glmm_stu0)
restu1<- simulateResiduals(cu_glmm_stu1)
restu2<- simulateResiduals(cu_glmm_stu2)


#Check gaussian model
plotQQunif(resgaus0) #check Q-Q
plotResiduals(resgaus0) #check residuals vs fitted, 
plotResiduals(resgaus0, cu_ch4_ow$season) #check heterocedasticity
plotResiduals(resgaus0, cu_ch4_ow$status) #check heterocedasticity
#Status heterocedasticity

plotQQunif(resgaus1) #check Q-Q
plotResiduals(resgaus1) #check residuals vs fitted, 
plotResiduals(resgaus1, cu_ch4_ow$season) #check heterocedasticity
plotResiduals(resgaus1, cu_ch4_ow$status) #check heterocedasticity
#Status heterocedasticity

plotQQunif(resgaus2) #check Q-Q
plotResiduals(resgaus2) #check residuals vs fitted, 
plotResiduals(resgaus2, cu_ch4_ow$season) #check heterocedasticity
plotResiduals(resgaus2, cu_ch4_ow$status) #check heterocedasticity
#Status heterocedasticity

#Check t_family model
check_convergence(cu_glmm_stu0)
plot(restu0) #check Q-Q, residuals vs fitted,  
plotResiduals(restu0, cu_ch4_ow$season) #check heterocedasticity
plotResiduals(restu0, cu_ch4_ow$status) #check heterocedasticity
#status heterocedasticity

check_convergence(cu_glmm_stu1)
plot(restu1) #check Q-Q, residuals vs fitted, 
plotResiduals(restu1, cu_ch4_ow$season) #check heterocedasticity
plotResiduals(restu1, cu_ch4_ow$status) #check heterocedasticity
#status heterocedasticity

check_convergence(cu_glmm_stu2)
plot(restu2) #check Q-Q, residuals vs fitted, 
plotResiduals(restu2, cu_ch4_ow$season) #check heterocedasticity
plotResiduals(restu2, cu_ch4_ow$status) #check heterocedasticity
#Status heterocedasticity

#GAUssian comparison: dispformula options
anova(cu_glmm_gaus0,cu_glmm_gaus1,cu_glmm_gaus2)

#t_family comparison: dispformula options
anova(cu_glmm_stu0,cu_glmm_stu1,cu_glmm_stu2)

#best vs best
anova(cu_glmm_gaus2,cu_glmm_stu2) 


#cu_ch4_ow best model?: #T_family, with season and status heterocedasticity
cu_ch4_ow_best<- cu_glmm_stu2
Anova(cu_ch4_ow_best,type = 3) 


#R2 of homocedastic model
cu_ch4_ow_best_homoc<- cu_glmm_stu0
r2(cu_ch4_ow_best_homoc)

#Singularity?
performance::check_singularity(cu_ch4_ow_best)

#Outliers?
testOutliers(cu_ch4_ow_best)

#plot obs vs predicted

data.frame(obs=cu_ch4_ow_best$frame$dailyflux_trans,
           pred=predict(cu_ch4_ow_best, type = "response"),
           strata=factor(cu_ch4_ow$strata, levels = c("bare","vegetated", "open water")),
           status=factor(cu_ch4_ow$status)) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=status),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "cu_ch4_ow Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

ggsave(filename = "CU_ch4_ow_Observed_VS_predicted.png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)



#Remove non-best models
rm(cu_glmm_gaus0,cu_glmm_gaus1,cu_glmm_gaus2,
   cu_glmm_stu0,cu_glmm_stu1,cu_glmm_stu2,
   resgaus0,resgaus1,resgaus2,
   restu0, restu1, restu2)





#CU_ch4_v---- 

cu_ch4_v<- data4models_strata %>% filter(casepilot=="CU"&ghgspecies=="ch4"&strata=="vegetated")

#See effect of transformation: 
table_trans %>%
  filter(casepilot=="CU",
         ghgspecies=="ch4", strata=="vegetated") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#FIT MODEL OPTIONS: 
cu_glmm_gaus0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = cu_ch4_v,
                        family = gaussian(),
                        dispformula = ~1)

cu_glmm_gaus1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = cu_ch4_v,
                        family = gaussian(),
                        dispformula = ~season)

cu_glmm_gaus2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = cu_ch4_v,
                        family = gaussian(),
                        dispformula = ~season+status)

cu_glmm_stu0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = cu_ch4_v,
                       family = t_family,
                       dispformula = ~1) 

cu_glmm_stu1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = cu_ch4_v,
                       family = t_family,
                       dispformula = ~season)

cu_glmm_stu2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = cu_ch4_v,
                       family = t_family,
                       dispformula = ~season+status)

#CALCULATE RESIDUALS
resgaus0<- simulateResiduals(cu_glmm_gaus0)
resgaus1<- simulateResiduals(cu_glmm_gaus1)
resgaus2<- simulateResiduals(cu_glmm_gaus2)

restu0<- simulateResiduals(cu_glmm_stu0)
restu1<- simulateResiduals(cu_glmm_stu1)
restu2<- simulateResiduals(cu_glmm_stu2)


#Check gaussian model
plotQQunif(resgaus0) #check Q-Q
plotResiduals(resgaus0) #check residuals vs fitted, 
plotResiduals(resgaus0, cu_ch4_v$season) #check heterocedasticity
plotResiduals(resgaus0, cu_ch4_v$status) #check heterocedasticity
#ALL pass

plotQQunif(resgaus1) #check Q-Q
plotResiduals(resgaus1) #check residuals vs fitted, 
plotResiduals(resgaus1, cu_ch4_v$season) #check heterocedasticity
plotResiduals(resgaus1, cu_ch4_v$status) #check heterocedasticity
#ALL PASS

plotQQunif(resgaus2) #check Q-Q
plotResiduals(resgaus2) #check residuals vs fitted, 
plotResiduals(resgaus2, cu_ch4_v$season) #check heterocedasticity
plotResiduals(resgaus2, cu_ch4_v$status) #check heterocedasticity
#ALL PASS

#Check t_family model
check_convergence(cu_glmm_stu0)
plot(restu0) #check Q-Q, residuals vs fitted,  
plotResiduals(restu0, cu_ch4_v$season) #check heterocedasticity
plotResiduals(restu0, cu_ch4_v$status) #check heterocedasticity
#ALL pass

check_convergence(cu_glmm_stu1)
plot(restu1) #check Q-Q, residuals vs fitted, 
plotResiduals(restu1, cu_ch4_v$season) #check heterocedasticity
plotResiduals(restu1, cu_ch4_v$status) #check heterocedasticity
#All pass

check_convergence(cu_glmm_stu2)
plot(restu2) #check Q-Q, residuals vs fitted, 
plotResiduals(restu2, cu_ch4_v$season) #check heterocedasticity
plotResiduals(restu2, cu_ch4_v$status) #check heterocedasticity
#DOES NOT CONVERGE

#GAUssian comparison: dispformula options
anova(cu_glmm_gaus0,cu_glmm_gaus1,cu_glmm_gaus2)

#t_family comparison: dispformula options
anova(cu_glmm_stu0,cu_glmm_stu1,cu_glmm_stu2)

#best vs best
anova(cu_glmm_gaus1,cu_glmm_stu1) 


#CU_ch4_v best model?: #Gaussian, with season heterocedasticity
cu_ch4_v_best<- cu_glmm_gaus1
Anova(cu_ch4_v_best,type = 3) 


#R2 of homocedastic model
cu_ch4_v_best_homoc<- cu_glmm_gaus0
r2(cu_ch4_v_best_homoc)

#Singularity?
performance::check_singularity(cu_ch4_v_best) #ALL MODELS ARE SINGULAR.
#Outliers?
testOutliers(cu_ch4_v_best)

#plot obs vs predicted

data.frame(obs=cu_ch4_v_best$frame$dailyflux_trans,
           pred=predict(cu_ch4_v_best, type = "response"),
           strata=factor(cu_ch4_v$strata, levels = c("bare","vegetated", "open water")),
           status=factor(cu_ch4_v$status)) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=status),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "CU_ch4_v Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

ggsave(filename = "CU_ch4_v_Observed_VS_predicted.png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)

#Remove non-best models
rm(cu_glmm_gaus0,cu_glmm_gaus1,cu_glmm_gaus2,
   cu_glmm_stu0,cu_glmm_stu1,cu_glmm_stu2,
   resgaus0,resgaus1,resgaus2,
   restu0, restu1, restu2)






#CU_GWPco2andch4_ow---- 

cu_GWPco2andch4_ow<- data4models_strata %>% filter(casepilot=="CU"&ghgspecies=="GWPco2andch4"&strata=="open water")

#See effect of transformation: 
table_trans %>%
  filter(casepilot=="CU",
         ghgspecies=="GWPco2andch4", strata=="open water") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#FIT MODEL OPTIONS: 
cu_glmm_gaus0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = cu_GWPco2andch4_ow,
                        family = gaussian(),
                        dispformula = ~1)

cu_glmm_gaus1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = cu_GWPco2andch4_ow,
                        family = gaussian(),
                        dispformula = ~season)

cu_glmm_gaus2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = cu_GWPco2andch4_ow,
                        family = gaussian(),
                        dispformula = ~season+status)

cu_glmm_stu0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = cu_GWPco2andch4_ow,
                       family = t_family,
                       dispformula = ~1) 

cu_glmm_stu1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = cu_GWPco2andch4_ow,
                       family = t_family,
                       dispformula = ~season)

cu_glmm_stu2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = cu_GWPco2andch4_ow,
                       family = t_family,
                       dispformula = ~season+status)

#CALCULATE RESIDUALS
resgaus0<- simulateResiduals(cu_glmm_gaus0)
resgaus1<- simulateResiduals(cu_glmm_gaus1)
resgaus2<- simulateResiduals(cu_glmm_gaus2)

restu0<- simulateResiduals(cu_glmm_stu0)
restu1<- simulateResiduals(cu_glmm_stu1)
restu2<- simulateResiduals(cu_glmm_stu2)


#Check gaussian model
plotQQunif(resgaus0) #check Q-Q
plotResiduals(resgaus0) #check residuals vs fitted, 
plotResiduals(resgaus0, cu_GWPco2andch4_ow$season) #check heterocedasticity
plotResiduals(resgaus0, cu_GWPco2andch4_ow$status) #check heterocedasticity
#Status and season heterocedasticity

plotQQunif(resgaus1) #check Q-Q
plotResiduals(resgaus1) #check residuals vs fitted, 
plotResiduals(resgaus1, cu_GWPco2andch4_ow$season) #check heterocedasticity
plotResiduals(resgaus1, cu_GWPco2andch4_ow$status) #check heterocedasticity
#Status and season heterocedasticity

plotQQunif(resgaus2) #check Q-Q
plotResiduals(resgaus2) #check residuals vs fitted, 
plotResiduals(resgaus2, cu_GWPco2andch4_ow$season) #check heterocedasticity
plotResiduals(resgaus2, cu_GWPco2andch4_ow$status) #check heterocedasticity
#Status heterocedasticity

#Check t_family model
check_convergence(cu_glmm_stu0)
plot(restu0) #check Q-Q, residuals vs fitted,  
plotResiduals(restu0, cu_GWPco2andch4_ow$season) #check heterocedasticity
plotResiduals(restu0, cu_GWPco2andch4_ow$status) #check heterocedasticity
#Status and season heterocedasticity

check_convergence(cu_glmm_stu1)
plot(restu1) #check Q-Q, residuals vs fitted, 
plotResiduals(restu1, cu_GWPco2andch4_ow$season) #check heterocedasticity
plotResiduals(restu1, cu_GWPco2andch4_ow$status) #check heterocedasticity
#Status and season heterocedasticity

check_convergence(cu_glmm_stu2)
plot(restu2) #check Q-Q, residuals vs fitted, 
plotResiduals(restu2, cu_GWPco2andch4_ow$season) #check heterocedasticity
plotResiduals(restu2, cu_GWPco2andch4_ow$status) #check heterocedasticity
#Status heterocedasticity

#GAUssian comparison: dispformula options
anova(cu_glmm_gaus0,cu_glmm_gaus1,cu_glmm_gaus2)

#t_family comparison: dispformula options
anova(cu_glmm_stu0,cu_glmm_stu1,cu_glmm_stu2)

#best vs best
anova(cu_glmm_gaus2,cu_glmm_stu2) 


#cu_GWPco2andch4_ow best model?: #T_family, with season and status heterocedasticity
cu_GWPco2andch4_ow_best<- cu_glmm_stu2
Anova(cu_GWPco2andch4_ow_best,type = 3) 


#R2 of homocedastic model
cu_GWPco2andch4_ow_best_homoc<- cu_glmm_stu0
r2(cu_GWPco2andch4_ow_best_homoc)

#Singularity?
performance::check_singularity(cu_GWPco2andch4_ow_best)

#Outliers?
testOutliers(cu_GWPco2andch4_ow_best)

#plot obs vs predicted

data.frame(obs=cu_GWPco2andch4_ow_best$frame$dailyflux_trans,
           pred=predict(cu_GWPco2andch4_ow_best, type = "response"),
           strata=factor(cu_GWPco2andch4_ow$strata, levels = c("bare","vegetated", "open water")),
           status=factor(cu_GWPco2andch4_ow$status)) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=status),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "cu_GWPco2andch4_ow Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

ggsave(filename = "cu_GWPco2andch4_ow_Observed_VS_predicted.png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)



#Remove non-best models
rm(cu_glmm_gaus0,cu_glmm_gaus1,cu_glmm_gaus2,
   cu_glmm_stu0,cu_glmm_stu1,cu_glmm_stu2,
   resgaus0,resgaus1,resgaus2,
   restu0, restu1, restu2)





#CU_GWPco2andch4_v---- 

cu_GWPco2andch4_v<- data4models_strata %>% filter(casepilot=="CU"&ghgspecies=="GWPco2andch4"&strata=="vegetated")

#See effect of transformation: 
table_trans %>%
  filter(casepilot=="CU",
         ghgspecies=="GWPco2andch4", strata=="vegetated") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#FIT MODEL OPTIONS: 
cu_glmm_gaus0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = cu_GWPco2andch4_v,
                        family = gaussian(),
                        dispformula = ~1)

cu_glmm_gaus1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = cu_GWPco2andch4_v,
                        family = gaussian(),
                        dispformula = ~season)

cu_glmm_gaus2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = cu_GWPco2andch4_v,
                        family = gaussian(),
                        dispformula = ~season+status)

cu_glmm_stu0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = cu_GWPco2andch4_v,
                       family = t_family,
                       dispformula = ~1) 

cu_glmm_stu1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = cu_GWPco2andch4_v,
                       family = t_family,
                       dispformula = ~season)

cu_glmm_stu2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = cu_GWPco2andch4_v,
                       family = t_family,
                       dispformula = ~season+status)

#CALCULATE RESIDUALS
resgaus0<- simulateResiduals(cu_glmm_gaus0)
resgaus1<- simulateResiduals(cu_glmm_gaus1)
resgaus2<- simulateResiduals(cu_glmm_gaus2)

restu0<- simulateResiduals(cu_glmm_stu0)
restu1<- simulateResiduals(cu_glmm_stu1)
restu2<- simulateResiduals(cu_glmm_stu2)


#Check gaussian model
plotQQunif(resgaus0) #check Q-Q
plotResiduals(resgaus0) #check residuals vs fitted, 
plotResiduals(resgaus0, cu_GWPco2andch4_v$season) #check heterocedasticity
plotResiduals(resgaus0, cu_GWPco2andch4_v$status) #check heterocedasticity
#Assumptions fail

plotQQunif(resgaus1) #check Q-Q
plotResiduals(resgaus1) #check residuals vs fitted, 
plotResiduals(resgaus1, cu_GWPco2andch4_v$season) #check heterocedasticity
plotResiduals(resgaus1, cu_GWPco2andch4_v$status) #check heterocedasticity
#ALL PASS

plotQQunif(resgaus2) #check Q-Q
plotResiduals(resgaus2) #check residuals vs fitted, 
plotResiduals(resgaus2, cu_GWPco2andch4_v$season) #check heterocedasticity
plotResiduals(resgaus2, cu_GWPco2andch4_v$status) #check heterocedasticity
#ALL PASS

#Check t_family model
check_convergence(cu_glmm_stu0)
plot(restu0) #check Q-Q, residuals vs fitted,  
plotResiduals(restu0, cu_GWPco2andch4_v$season) #check heterocedasticity
plotResiduals(restu0, cu_GWPco2andch4_v$status) #check heterocedasticity
#Assumptions fail

check_convergence(cu_glmm_stu1)
plot(restu1) #check Q-Q, residuals vs fitted, 
plotResiduals(restu1, cu_ch4_v$season) #check heterocedasticity
plotResiduals(restu1, cu_ch4_v$status) #check heterocedasticity
#Does not converge

check_convergence(cu_glmm_stu2)
plot(restu2) #check Q-Q, residuals vs fitted, 
plotResiduals(restu2, cu_GWPco2andch4_v$season) #check heterocedasticity
plotResiduals(restu2, cu_GWPco2andch4_v$status) #check heterocedasticity
#All pass 

#GAUssian comparison: dispformula options
anova(cu_glmm_gaus0,cu_glmm_gaus1,cu_glmm_gaus2)

#t_family comparison: dispformula options
anova(cu_glmm_stu0,cu_glmm_stu1,cu_glmm_stu2)

#best vs best
anova(cu_glmm_gaus1,cu_glmm_stu0) #GAUssian chosen


#CU_ch4_v best model?: #Gaussian, with season heterocedasticity
cu_GWPco2andch4_v_best<- cu_glmm_gaus1
Anova(cu_GWPco2andch4_v_best,type = 3) 


#R2 of homocedastic model
cu_GWPco2andch4_v_best_homoc<- cu_glmm_gaus0
r2(cu_GWPco2andch4_v_best_homoc)

#Singularity?
performance::check_singularity(cu_GWPco2andch4_v_best) #ALL MODELS ARE SINGULAR.
#Outliers?
testOutliers(cu_GWPco2andch4_v_best)

#plot obs vs predicted

data.frame(obs=cu_GWPco2andch4_v_best$frame$dailyflux_trans,
           pred=predict(cu_GWPco2andch4_v_best, type = "response"),
           strata=factor(cu_GWPco2andch4_v$strata, levels = c("bare","vegetated", "open water")),
           status=factor(cu_GWPco2andch4_v$status)) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=status),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "CU_GWPco2andch4_v Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

ggsave(filename = "CU_GWPco2andch4_v_Observed_VS_predicted.png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)

#Remove non-best models
rm(cu_glmm_gaus0,cu_glmm_gaus1,cu_glmm_gaus2,
   cu_glmm_stu0,cu_glmm_stu1,cu_glmm_stu2,
   resgaus0,resgaus1,resgaus2,
   restu0, restu1, restu2)



#VALENCIA  ------

#See interaction-level nobs 
data4models_strata %>%
  filter(casepilot == "VA", ghgspecies == "co2") %>%
  filter(!is.na(dailyflux_trans)) %>%
  group_by(status, season, strata) %>%
  summarise(nobs = n(), .groups = "drop") %>%
  complete(status, season, strata) %>%
  arrange(nobs)


#MISSING Preserved-S4-open water. WE cannot build the full model specifically for Open water. Additionally, only 1 observation for Altered-S3-bare: possible issues.  

#Potentially, we could combine open water and bare into a strata2 category: non-vegetated vs vegetated. 


#VA_co2_v---- 

#Models are always singular, they only capture seasonal variability as effects, not status. 

va_co2_v<- data4models_strata %>% filter(casepilot=="VA"&ghgspecies=="co2"&strata=="vegetated")

#See effect of transformation: 
table_trans %>%
  filter(casepilot=="VA",
         ghgspecies=="co2", strata=="vegetated") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#FIT MODEL OPTIONS: 
va_glmm_gaus0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = va_co2_v,
                        family = gaussian(),
                        dispformula = ~1)

va_glmm_gaus1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = va_co2_v,
                        family = gaussian(),
                        dispformula = ~season)

va_glmm_gaus2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = va_co2_v,
                        family = gaussian(),
                        dispformula = ~season+status)

va_glmm_stu0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = va_co2_v,
                       family = t_family,
                       dispformula = ~1) 

va_glmm_stu1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = va_co2_v,
                       family = t_family,
                       dispformula = ~season)

va_glmm_stu2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = va_co2_v,
                       family = t_family,
                       dispformula = ~season+status)

#CALCULATE RESIDUALS
resgaus0<- simulateResiduals(va_glmm_gaus0)
resgaus1<- simulateResiduals(va_glmm_gaus1)
resgaus2<- simulateResiduals(va_glmm_gaus2)

restu0<- simulateResiduals(va_glmm_stu0)
restu1<- simulateResiduals(va_glmm_stu1)
restu2<- simulateResiduals(va_glmm_stu2)


#Check gaussian model
plotQQunif(resgaus0) #check Q-Q
plotResiduals(resgaus0) #check residuals vs fitted, 
plotResiduals(resgaus0, va_co2_v$season) #check heterocedasticity
plotResiduals(resgaus0, va_co2_v$status) #check heterocedasticity
#Heterocedasticity

plotQQunif(resgaus1) #check Q-Q
plotResiduals(resgaus1) #check residuals vs fitted, 
plotResiduals(resgaus1, va_co2_v$season) #check heterocedasticity
plotResiduals(resgaus1, va_co2_v$status) #check heterocedasticity
#ALL pass

plotQQunif(resgaus2) #check Q-Q
plotResiduals(resgaus2) #check residuals vs fitted, 
plotResiduals(resgaus2, va_co2_v$season) #check heterocedasticity
plotResiduals(resgaus2, va_co2_v$status) #check heterocedasticity
#ALL pass

#Check t_family model
check_convergence(cu_glmm_stu0)
plot(restu0) #check Q-Q, residuals vs fitted,  
plotResiduals(restu0, va_co2_v$season) #check heterocedasticity
plotResiduals(restu0, va_co2_v$status) #check heterocedasticity
#heterocedasticity

check_convergence(cu_glmm_stu1)
plot(restu1) #check Q-Q, residuals vs fitted, 
plotResiduals(restu1, va_co2_v$season) #check heterocedasticity
plotResiduals(restu1, va_co2_v$status) #check heterocedasticity
#Heterocedasticity

check_convergence(cu_glmm_stu2)
plot(restu2) #check Q-Q, residuals vs fitted, 
plotResiduals(restu2, va_co2_v$season) #check heterocedasticity
plotResiduals(restu2, va_co2_v$status) #check heterocedasticity
#ALL pass

#GAUssian comparison: dispformula options
anova(va_glmm_gaus0,va_glmm_gaus1,va_glmm_gaus2)

#t_family comparison: dispformula options
anova(va_glmm_stu0,va_glmm_stu1,va_glmm_stu2)

#best vs best
anova(va_glmm_gaus1,va_glmm_stu1) 

#VA_co2_v best model?: #GAUSSIAN, with season heterocedasticity
va_co2_v_best<- va_glmm_gaus1
Anova(va_co2_v_best,type = 3) 


#R2 of homocedastic model
va_co2_v_best_homoc<- va_glmm_gaus0
r2(va_co2_v_best_homoc)
r2(va_co2_v_best_homoc,tolerance = 1e-10)

#Singularity?
performance::check_singularity(va_co2_v_best)

#Outliers?
testOutliers(va_co2_v_best)

#plot obs vs predicted

data.frame(obs=va_co2_v_best$frame$dailyflux_trans,
           pred=predict(va_co2_v_best, type = "response"),
           strata=factor(va_co2_v$strata, levels = c("bare","vegetated", "open water")),
           status=factor(va_co2_v$status)) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=status),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "VA_co2_v Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

ggsave(filename = "VA_co2_v_Observed_VS_predicted.png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)



#Remove non-best models
rm(cu_glmm_gaus0,cu_glmm_gaus1,cu_glmm_gaus2,
   cu_glmm_stu0,cu_glmm_stu1,cu_glmm_stu2,
   resgaus0,resgaus1,resgaus2,
   restu0, restu1, restu2)



#VA_co2_nonveg-----

#Check distribution across open-water vs bare (are those two comparable, in order to model them together?)

va_co2_nonveg<- data4models_strata %>% filter(casepilot=="VA"&ghgspecies=="co2"&strata!="vegetated")


bn_co2_nonveg<- bestNormalize(va_co2_nonveg$dailyflux,new_transforms = pseudolog_transform, standardize = TRUE, warn = TRUE, k = 5,allow_orderNorm = F)

plot_bestNormalize(bn_co2_nonveg)

va_co2_nonveg<- va_co2_nonveg %>% 
  mutate(dailyflux_trans=predict(bn_co2_nonveg, dailyflux))


data.frame(x=va_co2_nonveg$dailyflux,
           x_t=bn_co2_nonveg$x.t,
           strata=va_co2_nonveg$strata) %>% 
  ggplot(aes(x = strata, y = x_t)) +
  geom_violin(trim = FALSE, fill = "lightblue") +
  geom_boxplot(width = 0.1, fill = "white") +
  labs(title = "Distribution of dailyflux by strata",
       x = "Group",
       y = "Value") +
  theme_minimal()


#FIT MODEL OPTIONS: 
va_glmm_gaus0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = va_co2_nonveg,
                        family = gaussian(),
                        dispformula = ~1)

va_glmm_gaus1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = va_co2_nonveg,
                        family = gaussian(),
                        dispformula = ~season)

va_glmm_gaus2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = va_co2_nonveg,
                        family = gaussian(),
                        dispformula = ~season+status)

va_glmm_stu0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = va_co2_nonveg,
                       family = t_family,
                       dispformula = ~1) 

va_glmm_stu1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = va_co2_nonveg,
                       family = t_family,
                       dispformula = ~season)

va_glmm_stu2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = va_co2_nonveg,
                       family = t_family,
                       dispformula = ~season+status)

#CALCULATE RESIDUALS
resgaus0<- simulateResiduals(va_glmm_gaus0)
resgaus1<- simulateResiduals(va_glmm_gaus1)
resgaus2<- simulateResiduals(va_glmm_gaus2)

restu0<- simulateResiduals(va_glmm_stu0)
restu1<- simulateResiduals(va_glmm_stu1)
restu2<- simulateResiduals(va_glmm_stu2)


#Check gaussian model
plotQQunif(resgaus0) #check Q-Q
plotResiduals(resgaus0) #check residuals vs fitted, 
plotResiduals(resgaus0, va_co2_nonveg$season) #check heterocedasticity
plotResiduals(resgaus0, va_co2_nonveg$status) #check heterocedasticity
#Heterocedasticity

plotQQunif(resgaus1) #check Q-Q
plotResiduals(resgaus1) #check residuals vs fitted, 
plotResiduals(resgaus1, va_co2_nonveg$season) #check heterocedasticity
plotResiduals(resgaus1, va_co2_nonveg$status) #check heterocedasticity
#Heterocedasticity

plotQQunif(resgaus2) #check Q-Q
plotResiduals(resgaus2) #check residuals vs fitted, 
plotResiduals(resgaus2, va_co2_nonveg$season) #check heterocedasticity
plotResiduals(resgaus2, va_co2_nonveg$status) #check heterocedasticity
#Heterocedasticity

#Check t_family model
check_convergence(cu_glmm_stu0)
plot(restu0) #check Q-Q, residuals vs fitted,  
plotResiduals(restu0, va_co2_nonveg$season) #check heterocedasticity
plotResiduals(restu0, va_co2_nonveg$status) #check heterocedasticity
#heterocedasticity

check_convergence(cu_glmm_stu1)
plot(restu1) #check Q-Q, residuals vs fitted, 
plotResiduals(restu1, va_co2_nonveg$season) #check heterocedasticity
plotResiduals(restu1, va_co2_nonveg$status) #check heterocedasticity
#Heterocedasticity

check_convergence(cu_glmm_stu2)
plot(restu2) #check Q-Q, residuals vs fitted, 
plotResiduals(restu2, va_co2_nonveg$season) #check heterocedasticity
plotResiduals(restu2, va_co2_nonveg$status) #check heterocedasticity
#ALL pass

#GAUssian comparison: dispformula options
anova(va_glmm_gaus0,va_glmm_gaus1,va_glmm_gaus2)

#t_family comparison: dispformula options
anova(va_glmm_stu0,va_glmm_stu1,va_glmm_stu2)

#best vs best
anova(va_glmm_gaus2,va_glmm_stu2) 

#VA_co2_v best model?: #GAUSSIAN, with season and status heterocedasticity
va_co2_nonveg_best<- va_glmm_gaus2
Anova(va_co2_nonveg_best,type = 3) 


#R2 of homocedastic model
va_co2_nonveg_best_homoc<- va_glmm_gaus0
r2(va_co2_nonveg_best_homoc)

#Singularity?
performance::check_singularity(va_co2_nonveg_best)

#Outliers?
testOutliers(va_co2_nonveg_best)

#plot obs vs predicted

data.frame(obs=va_co2_nonveg_best$frame$dailyflux_trans,
           pred=predict(va_co2_nonveg_best, type = "response"),
           strata=factor(va_co2_nonveg$strata, levels = c("bare","vegetated", "open water")),
           status=factor(va_co2_nonveg$status)) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=status),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "va_co2_nonveg Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

ggsave(filename = "va_co2_nonveg_Observed_VS_predicted.png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)



#Remove non-best models
rm(cu_glmm_gaus0,cu_glmm_gaus1,cu_glmm_gaus2,
   cu_glmm_stu0,cu_glmm_stu1,cu_glmm_stu2,
   resgaus0,resgaus1,resgaus2,
   restu0, restu1, restu2)



#VA_co2_ow---- 

va_co2_ow<- data4models_strata %>% filter(casepilot=="VA"&ghgspecies=="co2"&strata=="open water")

#See effect of transformation: 
table_trans %>%
  filter(casepilot=="VA",
         ghgspecies=="co2", strata=="open water") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#FIT MODEL OPTIONS: 
va_glmm_gaus0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = va_co2_ow,
                        family = gaussian(),
                        dispformula = ~1)

va_glmm_gaus1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = va_co2_ow,
                        family = gaussian(),
                        dispformula = ~season)

va_glmm_gaus2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = va_co2_ow,
                        family = gaussian(),
                        dispformula = ~season+status)

va_glmm_stu0<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = va_co2_ow,
                       family = t_family,
                       dispformula = ~1) 

va_glmm_stu1<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = va_co2_ow,
                       family = t_family,
                       dispformula = ~season)

va_glmm_stu2<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                       data = va_co2_ow,
                       family = t_family,
                       dispformula = ~season+status)

#CALCULATE RESIDUALS
resgaus0<- simulateResiduals(va_glmm_gaus0)
resgaus1<- simulateResiduals(va_glmm_gaus1)
resgaus2<- simulateResiduals(va_glmm_gaus2)

restu0<- simulateResiduals(va_glmm_stu0)
restu1<- simulateResiduals(va_glmm_stu1)
restu2<- simulateResiduals(va_glmm_stu2)


#Check gaussian model
plotQQunif(resgaus0) #check Q-Q
plotResiduals(resgaus0) #check residuals vs fitted, 
plotResiduals(resgaus0, va_co2_ow$season) #check heterocedasticity
plotResiduals(resgaus0, va_co2_ow$status) #check heterocedasticity
#Assumptions fail

plotQQunif(resgaus1) #check Q-Q
plotResiduals(resgaus1) #check residuals vs fitted, 
plotResiduals(resgaus1, va_co2_ow$season) #check heterocedasticity
plotResiduals(resgaus1, va_co2_ow$status) #check heterocedasticity
#Heterocedasticity

plotQQunif(resgaus2) #check Q-Q
plotResiduals(resgaus2) #check residuals vs fitted, 
plotResiduals(resgaus2, va_co2_ow$season) #check heterocedasticity
plotResiduals(resgaus2, va_co2_ow$status) #check heterocedasticity
#Heterocedasticity

#Check t_family model
check_convergence(cu_glmm_stu0)
plot(restu0) #check Q-Q, residuals vs fitted,  
plotResiduals(restu0, va_co2_ow$season) #check heterocedasticity
plotResiduals(restu0, va_co2_ow$status) #check heterocedasticity
#heterocedasticity

check_convergence(cu_glmm_stu1)
plot(restu1) #check Q-Q, residuals vs fitted, 
plotResiduals(restu1, va_co2_ow$season) #check heterocedasticity
plotResiduals(restu1, va_co2_ow$status) #check heterocedasticity
#Heterocedasticity

check_convergence(cu_glmm_stu2)
plot(restu2) #check Q-Q, residuals vs fitted, 
plotResiduals(restu2, va_co2_ow$season) #check heterocedasticity
plotResiduals(restu2, va_co2_ow$status) #check heterocedasticity
#Heterocedasticity

#GAUssian comparison: dispformula options
anova(va_glmm_gaus0,va_glmm_gaus1,va_glmm_gaus2)

#t_family comparison: dispformula options
anova(va_glmm_stu0,va_glmm_stu1,va_glmm_stu2)

#best vs best
anova(va_glmm_gaus2,va_glmm_stu0) #T_family, with season heterocedasticity


#CU_co2_ow best model?: #T_family, with season heterocedasticity
va_co2_ow_best<- va_glmm_gaus0
Anova(va_co2_ow_best,type = 3) 


#R2 of homocedastic model
va_co2_ow_best_homoc<- va_glmm_stu0
r2(va_co2_ow_best_homoc)

#Singularity?
performance::check_singularity(va_co2_ow_best)

#Outliers?
testOutliers(va_co2_ow_best)

#plot obs vs predicted

data.frame(obs=va_co2_ow_best$frame$dailyflux_trans,
           pred=predict(va_co2_ow_best, type = "response"),
           strata=factor(va_co2_ow$strata, levels = c("bare","vegetated", "open water")),
           status=factor(va_co2_ow$status)) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=status),size = 1)+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "VA_co2_ow Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

ggsave(filename = "VA_co2_ow_Observed_VS_predicted.png", 
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)



#Remove non-best models
rm(cu_glmm_gaus0,cu_glmm_gaus1,cu_glmm_gaus2,
   cu_glmm_stu0,cu_glmm_stu1,cu_glmm_stu2,
   resgaus0,resgaus1,resgaus2,
   restu0, restu1, restu2)





#Save strata-best models--------

#Homocedastic list: casepilot_ghgspecies_strata
homoced_model_list_strata<- list(
  "CU_co2_open water" = cu_co2_ow_best_homoc,
  "CU_co2_vegetated" = cu_co2_v_best_homoc,
  "CU_ch4_open water" = cu_ch4_ow_best_homoc,
  "CU_ch4_vegetated" = cu_ch4_v_best_homoc,
  "CU_GWPco2andch4_open water" = cu_GWPco2andch4_ow_best_homoc,
  "CU_GWPco2andch4_vegetated" = cu_GWPco2andch4_v_best_homoc)

#List of best models for each casepilot (to get significances, pseudoR2s and residual tests)
best_model_list_strata<- list(
  "CU_co2_open water" = cu_co2_ow_best,
  "CU_co2_vegetated" = cu_co2_v_best,
  "CU_ch4_open water" = cu_ch4_ow_best,
  "CU_ch4_vegetated" = cu_ch4_v_best,
  "CU_GWPco2andch4_open water" = cu_GWPco2andch4_ow_best,
  "CU_GWPco2andch4_vegetated" = cu_GWPco2andch4_v_best)



#Overall fit -----
#Structure, pseudoR2 of best model:
strata_best_model_fit<- get_pseudoR2s(best_model_list_strata)
#MANY warnings of "false convergence"



#R2c and R2m of homocedastic model:
strata_homoc_r2<- get_homoc_R2s(homoced_model_list_strata)


#DHARMa Residual diagnostics of best model:
strata_resid_diag_summary <- summarize_dharma_diagnostics(best_model_list_strata)

print(strata_resid_diag_summary)



#Obtain significance for main effects of best_models
strata_results_anova<-get_anova_results(best_model_list_strata) 
print(strata_results_anova)

#Add symbol on top of p_value

#Function to get symbols
pval_to_symbol <- function(p) {
  ifelse(p < 0.001, "***",
         ifelse(p < 0.01, "**",
                ifelse(p < 0.05, "*",
                       ifelse(p < 0.1, ".", "ns")
                )
         )
  )
}

#Add symbols from significance pvalues
strata_results_anova<- strata_results_anova %>% 
  mutate(across(c(intercept_pval,status_pval,season_pval,interaction_pval),pval_to_symbol,.names = "symbol_{.col}"))


#Intercept refers to whether the grand mean (on the transformed scale) is significantly different than cero (Is the overall net exchange of CH4 balanced between uptake and emission)
#Status refers to whether the status (averaged across seasons) has a significant effect on flux. 
#Season refers to whether the season (averaged across status) has a significant effect on flux.
#Interaction refers to whether the effect of status depends on season — i.e., the difference between status levels changes across seasons.


strata_model_outputs<- strata_results_anova %>% 
  full_join(strata_best_model_fit) %>% 
  full_join(strata_homoc_r2 %>% dplyr::select(dataset, homoced_R2m, homoced_R2c))




#CALCULATE strata-EMMEANS-------

# Initialize list for storing comparisons
comparisons_list <- list()

# Loop to extract all comparisons for every dataset best model
for (dataset in names(best_model_list_strata)) {
  cp_model <- best_model_list_strata[[dataset]]
  
  #extract casepilot name and ghgspecies from model-list names
  casepilot_name<- sub("_.*", "", dataset)
  ghgspecies<- sub("_.*", "", sub(paste0(casepilot_name,"_"),"",dataset))
  stratum<- sub(paste0(casepilot_name,"_",ghgspecies,"_"),"", dataset)
  
  # Status comparisons:
  status_emmeans <- emmeans(cp_model, ~status)
  status_comparisons <- cld(status_emmeans, 
                            alpha = 0.05,
                            Letters = letters,
                            adjust = "sidak") %>%
    mutate(comparison = "status",
           ghgspecies = ghgspecies,
           casepilot = casepilot_name,
           strata= stratum, 
           season = NA) %>%
    rename(lower.CL = any_of("asymp.LCL"),
           upper.CL = any_of("asymp.UCL"))
  
  # Season comparisons: 
  season_emmeans <- emmeans(cp_model, ~season)
  season_comparisons <- cld(season_emmeans, 
                            alpha = 0.05,
                            Letters = letters,
                            adjust = "sidak") %>%
    mutate(comparison = "season",
           ghgspecies = ghgspecies,
           casepilot = casepilot_name,
           strata= stratum, 
           status = NA) %>%
    rename(lower.CL = any_of("asymp.LCL"),
           upper.CL = any_of("asymp.UCL"))
  
  # Interaction comparisons: 
  interaction_emmeans <- emmeans(cp_model, ~status:season)
  interaction_comparisons <- cld(interaction_emmeans, 
                                 alpha = 0.05,
                                 Letters = letters,
                                 adjust = "sidak") %>%
    mutate(comparison = "interaction",
           ghgspecies = ghgspecies,
           casepilot = casepilot_name,
           strata= stratum) %>%
    rename(lower.CL = any_of("asymp.LCL"),
           upper.CL = any_of("asymp.UCL"))
  
  # Join and format CP comparisons
  all_comparisons <- status_comparisons %>%
    full_join(season_comparisons) %>%
    full_join(interaction_comparisons) %>%
    rename(group_letter = .group) %>%
    mutate(
      comparison = factor(comparison, levels = c("status", "season", "interaction")),
      season = factor(season, levels = c("S1", "S2", "S3", "S4")),
      status = factor(status, levels = c("Altered", "Preserved", "Restored")),
      group_letter = gsub(" ", "", group_letter),
      
      # Create key for accessing BestNormalize object
      key_trans = paste(casepilot_name, ghgspecies,stratum, sep = "_")
    ) %>%
    rowwise() %>% 
    # BACK-TRANSFORMATION of emmean, SE and CL using BestNormalize object
    mutate(
      emmean_bt = if (!is.null(bn_list_strata[[key_trans]])) predict(bn_list_strata[[key_trans]], emmean, inverse = TRUE) else NA_real_,
      SE_bt = if (!is.null(bn_list_strata[[key_trans]])) predict(bn_list_strata[[key_trans]], SE, inverse = TRUE) else NA_real_,
      lower.CL_bt = if (!is.null(bn_list_strata[[key_trans]])) predict(bn_list_strata[[key_trans]], lower.CL, inverse = TRUE) else NA_real_,
      upper.CL_bt = if (!is.null(bn_list_strata[[key_trans]])) predict(bn_list_strata[[key_trans]], upper.CL, inverse = TRUE) else NA_real_
    ) %>%
    
    dplyr::select(any_of(c(
      "ghgspecies", "casepilot","strata","comparison", "status", "season", "df",
      "emmean", "SE", "lower.CL", "upper.CL", "group_letter",
      "emmean_bt", "SE_bt", "lower.CL_bt", "upper.CL_bt"
    ))) %>%
    arrange(ghgspecies, casepilot,strata, comparison, status, season)
  
  # Store casepilot results
  comparisons_list[[dataset]] <- all_comparisons
}

#Remove within-loop objects
rm(all_comparisons, status_comparisons,season_comparisons,interaction_comparisons,
   interaction_emmeans,status_emmeans, season_emmeans)


# Unir todos los resultados en un solo data.frame
all_emmeans_comparisons_strata <- bind_rows(comparisons_list)



#PLOT strata-EMMmeans------

#Plots emmeans in molar scale:  Emm +- SE
all_emmeans_comparisons_strata %>% 
  filter(casepilot == "CU", comparison == "status") %>%
  ggplot(., aes(x = status, y = emmean_bt, label = group_letter)) +
  geom_errorbar(aes(ymin = emmean_bt - SE_bt, ymax = emmean_bt + SE_bt),
                linewidth = 0.2, width = 0.2) +
  geom_point(shape = 15, size = 4, aes(col=status)) +
  geom_text(nudge_x = 0.2, color = "black") +
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw() +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  labs(title = expression(Status~effect~on~GHG~NEE~(per~strata)),
       subtitle = paste0("EMM with SE"), 
       y=expression(EMM~(mmol~GHG~d^-1~m^-2)),
       x="Conservation status",
       col=paste0("Status"))+
  facet_wrap(facets=vars(ghgspecies, strata),scales="free")



#Plots emmeans in molar scale:  Emm +- CL
all_emmeans_comparisons_strata %>% 
  filter(casepilot == "CU", comparison == "status") %>%
  ggplot(., aes(x = status, y = emmean_bt, label = group_letter)) +
  geom_errorbar(aes(ymin = lower.CL_bt, ymax =upper.CL_bt),
                linewidth = 0.2, width = 0.2) +
  geom_point(shape = 15, size = 4, aes(col=status)) +
  geom_text(nudge_x = 0.2, color = "black") +
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw() +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  labs(title = expression(CURONIAN:Status~effect~on~GHG~NEE~(per~strata)),
       subtitle = paste0("EMM with Confidence Limits (alpha = 0.05)"), 
       y=expression(EMM~(mmol~GHG~d^-1~m^-2)),
       x="Conservation status",
       col=paste0("Status"))+
  facet_wrap(facets=vars(ghgspecies, strata),scales="free")




#FUTHER CHECKS: 

# check for autorcorrelation structure: 
data4models_strata %>%
  filter(casepilot=="RI") %>% 
  # filter(strata == "open water") %>%
  mutate(season_num=parse_number(as.character(season))) %>% 
  arrange(subsite, season_num) %>%
  group_by(subsite) %>%
  summarise(acf = acf(dailyflux_trans, plot = FALSE)$acf[2]) %>%
  ggplot(aes(x = subsite, y = acf)) +
  geom_col() +
  labs(title = "Lag-1 Autocorrelation by Site")

acf()

#Specify temporal autocorrelation structure in model
library(glmmTMB)

glmmTMB(dailyflux_trans ~ status * season + (1|site),
        data = data4models_strata,
        family = gaussian(),
        control = glmmTMBControl(profile = TRUE),
        dispformula = ~1,
        ziformula = ~0,
        map = NULL,
        start = NULL,
        REML = FALSE,
        covstruct = list(ar1 = ~date | site)) #date should be numeric and be oredered (eg. days since start of sampling)
        
acf(resgaus0$fittedResiduals)

acf(residuals(ri_best_co2), main = "ACF of residuals")

