#GHGpaper_modelCoresData.R

#Author: Miguel Cabrera-Brufau
#Date: October 2025
#Project: Restore4cs



#Description----
#This scrip is used to model the effect of restoration in each casepilot using data from the incubation of cores. The analysis is repeated for every GHGspecies available (co2, ch4, n2o) and for the combined GWPtotal. 

#For N2O and GWPtotal we will only have data for seasons 2,3 and 4 (data for S1 was lost).

#DECISSIONS: 
#Modelling approach: GLMMtmb with general formula dailyflux~status*season + (1|subsite), subsite as random to account for repeated samplings and subsite-specific patterns

#Contrast: Set contrast options to "contr. sum", so that we are not using any one level of a factor as the reference level, but rather the tests will asses whether there is an overall average effect of one factor across all levels of the other factors. For example, is there an effect of status across all levels of season?




#STEPS: 
#1. Import and format data
#2. Transform data: Using BestNormalize package, find and apply for each casepilot, the transformation that maximizes normality. Allow pseudo-log transformation (with optimized sigma). 
#3. Optimize and select best model for each casepilot*GHGspecies combo (co2, ch4, n2o, GWPtotal).
# 3.1. Decide final structure: gaussian vs t-family distribution
# 3.2. Detect and treat outliers (if present).

#4. Export model results:
# 4.1. Allghg_summary_coresmodels.csv(normalization, structure, distributionfamily, dispformula, significances, Homoced_R2c, Homoced_R2m, PredObs_cor, nobs, extra(AIC?, BIC?, check common practices for reporting).   
# 4.2. Allghg_Residualtests_Coresmodels.csv: pvalues of residual tests for models (model assumptions)
# 4.3. Allghg_emmeans-posthoc_coresmodels.csv: estimated marginal means for each casepilot*GHGspecies and post-hoc groups. comparison for interaction is status|season (do not make all-pairwise) 



rm(list = ls()) # clear workspace
cat("/014") # clear console


# ---- Packages ----
#For data-handling & plotting
library(tidyverse)
library(readxl)
library(lubridate)
library(zoo)
library(ggpubr)
library(ggtext)#for plot captions
require(purrr)
require(data.table)
require(tools)
library(hms)
library(suncalc)



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


#Repo-functions
repo_root <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
files.sources = list.files(path = paste0(repo_root,"/functions"), full.names = T)
for (f in files.sources){source(f)}



# ---- Directories ----
dropbox_root <- "C:/Users/Miguel/Dropbox/RESTORE4Cs - Fieldwork/Data" # You have to make sure this is pointing to the write folder on your local machine


#Path to datasets and to save model outputs: 
paper_path<- paste0(dropbox_root,"/GHG/Main_paper_results/")
#Path to exploratory figures:
extraplots_path<- paste0(paper_path,"Exploratory_figures/")



#0.Contrast options-----
#NOTES on contrasts: in R by default, contrast is "contr.treatment" wich uses the first level of each factor as the "reference" level and all subsequent as "treatment" levels. This implies that, with default contrast type, when looking at the main effects of my model, what is shown for each main effect is whether there is a significant effect of 1 factor (eg. status) at the reference level of all other factors (i.e. season). What we want is to asses whether there is an overall average effect of status across all levels of season. For this purpose we need to set contrasts to "contr. sum", and always call Anova (model, type="III").

#Set contrasts to contr.sum (for the whole R session)
options(contrasts = c("contr.sum", "contr.poly"))



#0. Import data--------

#Import formated data (see formatting in GHGpaper_prepData.R), remove NAs and factor relevant variables: 

data4models<- read.csv(paste0(paper_path, "CoresData4paper.csv")) %>% 
  #Remove NAs
  filter(!is.na(dailyflux)) %>%
  #Factors 
  mutate(season=factor(season, ordered = F),
         status=factor(status, ordered = F),
         status_difalter=factor(status_difalter, ordered=F))

str(data4models)


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



#1.Transformations------

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
#See effect of transformation:
table_trans %>%
  filter(casepilot=="VA"&
           ghgspecies=="n2o") %>%
  pull(bn_result) %>%
  pluck(1) %>%
  plot_bestNormalize()


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


##best_CASEPILOT_GHGSPECIES -------
#REPEAT for each CASEPILOT*GHGSPECIES combo. 

#Subset data and check transformation: 
CASEPILOT_GHGSPECIES<- data4models %>% filter(casepilot=="CASEPILOT"&ghgspecies=="GHGSPECIES")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="CASEPILOT"&
           ghgspecies=="GHGSPECIES") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#GAUSSIAN OPTION, status*season
m1_gaus<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = CASEPILOT_GHGSPECIES,
                           family = gaussian(),
                           dispformula = ~1)
#Evaluate model:
check_convergence(m1_gaus) #OK
check_singularity(m1_gaus,tolerance = 1e-8) #OK
r2(m1_gaus, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m1_gaus)
summary(m1_gaus)
res<- simulateResiduals(m1_gaus)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus)
em<- emmeans(m1_gaus, ~ status)
summary(em)
pairs(em)
DHARMa::testOutliers(res)
#DECISSION: 


#T_family OPTION, status*season
m2_t<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = CASEPILOT_GHGSPECIES,
                        family = t_family(),
                        dispformula = ~1)
#Evaluate basic model:
check_convergence(m2_t)
check_singularity(m2_t,tolerance = 1e-8)
r2(m2_t, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m2_t)
summary(m2_t)
res<- simulateResiduals(m2_t)
plotQQunif(res)
plotResiduals(res)
DHARMa::testOutliers(res)
check_residuals(m2_t)
em<- emmeans(m2_t, ~ status)
summary(em)
pairs(em)
#DECISSION: 


#Compare models (if needed):AIC & pval
anova(m1_gaus,m2_t)


#Set as best
best_CASEPILOT_GHGSPECIES<- m2_t


#Plot obs vs pred 
data.frame(obs=best_CASEPILOT_GHGSPECIES$frame$dailyflux_trans,
           pred=predict(best_CASEPILOT_GHGSPECIES, type = "response"),
           status=best_CASEPILOT_GHGSPECIES$frame$status) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=status))+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "CASEPILOT_GHGSPECIES Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

#Save obs vs pred plot
# ggsave(filename = "coremodels_CASEPILOT_GHGSPECIES_Observed_VS_predicted.png", 
#        width = 115,height = 120,units = "mm",device = "png",dpi = 400,
#        path = extraplots_path
# )






#CO2 models------

## best_CA_co2 -------

#Subset data and check transformation: 
CA_co2<- data4models %>% filter(casepilot=="CA"&ghgspecies=="co2")

#TRANSFORMATION via bestNormalize creates issues


CA_co2 %>% 
  ggplot(aes(x=log10(dailyflux+4)))+
  geom_histogram()

#See effect of transformation:
table_trans %>%
  filter(casepilot=="CA"&
           ghgspecies=="co2") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()

#See effect of transformation:
ca_bn<- table_trans %>%
  filter(casepilot=="CA"&
           ghgspecies=="co2") %>%
  pull(bn_result) %>%
  pluck(1)


#Pseudolog tranformation (not chosen by bestNormalize) achieves more normal-looking distribution
ca_bn$other_transforms$best_pseudo_log$x.t %>% 
  gghistogram()
#even when log_b (x+a) has a lower Pearson's P (thats why bestNormalize chooses the distribution)
ca_bn$other_transforms$best_pseudo_log$x.t %>% shapiro.test()
ks.test(ca_bn$other_transforms$best_pseudo_log$x.t, "pnorm")



#GAUSSIAN OPTION, status*season
m1_gaus<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                  data = CA_co2,
                  family = gaussian(),
                  dispformula = ~1)
#Evaluate model:
check_convergence(m1_gaus) #OK
check_singularity(m1_gaus,tolerance = 1e-8) #OK
r2(m1_gaus, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m1_gaus)
summary(m1_gaus)
res<- simulateResiduals(m1_gaus)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus)
em<- emmeans(m1_gaus, ~ status)
summary(em)
pairs(em)
DHARMa::testOutliers(res)
#DECISSION: RESIDUALS FAIL


#T_family OPTION, status*season
m2_t<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
               data = CA_co2,
               family = t_family(),
               dispformula = ~1)
#Evaluate basic model:
check_convergence(m2_t)
check_singularity(m2_t,tolerance = 1e-8)
r2(m2_t, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m2_t)
summary(m2_t)
res<- simulateResiduals(m2_t)
plotQQunif(res)
plotResiduals(res)
DHARMa::testOutliers(res)
check_residuals(m2_t)
em<- emmeans(m2_t, ~ status)
summary(em)
pairs(em)
#DECISSION: OK model


#Compare models (if needed):AIC & pval
anova(m1_gaus,m2_t)


#Set as best
best_CA_co2<- m2_t


#Plot obs vs pred 
data.frame(obs=best_CA_co2$frame$dailyflux_trans,
           pred=predict(best_CA_co2, type = "response"),
           status=best_CA_co2$frame$status) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=status))+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "CA_co2 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

#Save obs vs pred plot
ggsave(filename = "coremodels_CA_co2_Observed_VS_predicted.png",
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)



## best_CU_co2 -------
#REPEAT for each CU*co2 combo. 

#Subset data and check transformation: 
CU_co2<- data4models %>% filter(casepilot=="CU"&ghgspecies=="co2")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="CU"&
           ghgspecies=="co2") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#GAUSSIAN OPTION, status*season
m1_gaus<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                  data = CU_co2,
                  family = gaussian(),
                  dispformula = ~1)
#Evaluate model:
check_convergence(m1_gaus) #OK
check_singularity(m1_gaus,tolerance = 1e-8) #OK
r2(m1_gaus, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m1_gaus)
summary(m1_gaus)
res<- simulateResiduals(m1_gaus)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus)
em<- emmeans(m1_gaus, ~ status)
summary(em)
pairs(em)
DHARMa::testOutliers(res)
#DECISSION: Model is OK


#T_family OPTION, status*season
m2_t<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
               data = CU_co2,
               family = t_family(),
               dispformula = ~1)
#Evaluate basic model:
check_convergence(m2_t)
check_singularity(m2_t,tolerance = 1e-8)
r2(m2_t, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m2_t)
summary(m2_t)
res<- simulateResiduals(m2_t)
plotQQunif(res)
plotResiduals(res)
DHARMa::testOutliers(res)
check_residuals(m2_t)
em<- emmeans(m2_t, ~ status)
summary(em)
pairs(em)
#DECISSION: 


#Compare models (if needed):AIC & pval
anova(m1_gaus,m2_t)


#Set as best
best_CU_co2<- m1_gaus


#Plot obs vs pred 
data.frame(obs=best_CU_co2$frame$dailyflux_trans,
           pred=predict(best_CU_co2, type = "response"),
           status=best_CU_co2$frame$status) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=status))+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "CU_co2 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

#Save obs vs pred plot
ggsave(filename = "coremodels_CU_co2_Observed_VS_predicted.png",
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)



## best_DA_co2 -------

#Subset data and check transformation: 
DA_co2<- data4models %>% filter(casepilot=="DA"&ghgspecies=="co2")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="DA"&
           ghgspecies=="co2") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#GAUSSIAN OPTION, status*season
m1_gaus<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                  data = DA_co2,
                  family = gaussian(),
                  dispformula = ~1)
#Evaluate model:
check_convergence(m1_gaus) #OK
check_singularity(m1_gaus,tolerance = 1e-8) #OK
r2(m1_gaus, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m1_gaus)
summary(m1_gaus)
res<- simulateResiduals(m1_gaus)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus)
em<- emmeans(m1_gaus, ~ status)
summary(em)
pairs(em)
DHARMa::testOutliers(res)
#DECISSION: OK-ish residuals


#T_family OPTION, status*season
m2_t<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
               data = DA_co2,
               family = t_family(),
               dispformula = ~1)
#Evaluate basic model:
check_convergence(m2_t)
check_singularity(m2_t,tolerance = 1e-8)
r2(m2_t, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m2_t)
summary(m2_t)
res<- simulateResiduals(m2_t)
plotQQunif(res)
plotResiduals(res)
DHARMa::testOutliers(res)
check_residuals(m2_t)
em<- emmeans(m2_t, ~ status)
summary(em)
pairs(em)
#DECISSION: OK-ish residuals


#Compare models (if needed):AIC & pval
anova(m1_gaus,m2_t)


#Set as best
best_DA_co2<- m1_gaus


#Plot obs vs pred 
data.frame(obs=best_DA_co2$frame$dailyflux_trans,
           pred=predict(best_DA_co2, type = "response"),
           status=best_DA_co2$frame$status) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=status))+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "DA_co2 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

#Save obs vs pred plot
ggsave(filename = "coremodels_DA_co2_Observed_VS_predicted.png",
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)


## best_DU_co2 -------

#Subset data and check transformation: 
DU_co2<- data4models %>% filter(casepilot=="DU"&ghgspecies=="co2")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="DU"&
           ghgspecies=="co2") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#GAUSSIAN OPTION, status*season
m1_gaus<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                  data = DU_co2,
                  family = gaussian(),
                  dispformula = ~1)
#Evaluate model:
check_convergence(m1_gaus) #OK
check_singularity(m1_gaus,tolerance = 1e-8) #OK
r2(m1_gaus, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m1_gaus)
summary(m1_gaus)
res<- simulateResiduals(m1_gaus)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus)
em<- emmeans(m1_gaus, ~ status)
summary(em)
pairs(em)
DHARMa::testOutliers(res)
#DECISSION: GOOD Model


#T_family OPTION, status*season
m2_t<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
               data = DU_co2,
               family = t_family(),
               dispformula = ~1)
#Evaluate basic model:
check_convergence(m2_t)
check_singularity(m2_t,tolerance = 1e-8)
r2(m2_t, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m2_t)
summary(m2_t)
res<- simulateResiduals(m2_t)
plotQQunif(res)
plotResiduals(res)
DHARMa::testOutliers(res)
check_residuals(m2_t)
em<- emmeans(m2_t, ~ status)
summary(em)
pairs(em)
#DECISSION: 


#Compare models (if needed):AIC & pval
anova(m1_gaus,m2_t)


#Set as best
best_DU_co2<- m1_gaus


#Plot obs vs pred 
data.frame(obs=best_DU_co2$frame$dailyflux_trans,
           pred=predict(best_DU_co2, type = "response"),
           status=best_DU_co2$frame$status) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=status))+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "DU_co2 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

#Save obs vs pred plot
ggsave(filename = "coremodels_DU_co2_Observed_VS_predicted.png",
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)



## best_RI_co2-------

#Subset data and check transformation: 
RI_co2<- data4models %>% filter(casepilot=="RI"&ghgspecies=="co2")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="RI"&
           ghgspecies=="co2") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#GAUSSIAN OPTION, status*season
m1_gaus<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                  data = RI_co2,
                  family = gaussian(),
                  dispformula = ~1)
#Evaluate model:
check_convergence(m1_gaus) #OK
check_singularity(m1_gaus,tolerance = 1e-8) #OK
r2(m1_gaus, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m1_gaus)
summary(m1_gaus)
res<- simulateResiduals(m1_gaus)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus)
em<- emmeans(m1_gaus, ~ status)
summary(em)
pairs(em)
DHARMa::testOutliers(res)
#DECISSION: Bad residuals


#T_family OPTION, status*season
m2_t<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
               data = RI_co2,
               family = t_family(),
               dispformula = ~1)
#Evaluate basic model:
check_convergence(m2_t)
check_singularity(m2_t,tolerance = 1e-8)
r2(m2_t, tolerance=1e-11) #Decrease tolerance to force R2 calculation when singular
Anova(m2_t)
summary(m2_t)
res<- simulateResiduals(m2_t)
plotQQunif(res)
plotResiduals(res)
DHARMa::testOutliers(res)
check_residuals(m2_t)
em<- emmeans(m2_t, ~ status)
summary(em)
pairs(em)
#DECISSION: GOOD MODEL


#Compare models (if needed):AIC & pval
anova(m1_gaus,m2_t)


#Set as best
best_RI_co2<- m2_t


#Plot obs vs pred 
data.frame(obs=best_RI_co2$frame$dailyflux_trans,
           pred=predict(best_RI_co2, type = "response"),
           status=best_RI_co2$frame$status) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=status))+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "RI_co2 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

#Save obs vs pred plot
ggsave(filename = "coremodels_RI_co2_Observed_VS_predicted.png",
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)


## best_VA_co2 -------


#Subset data and check transformation: 
VA_co2<- data4models %>% filter(casepilot=="VA"&ghgspecies=="co2")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="VA"&
           ghgspecies=="co2") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#GAUSSIAN OPTION, status*season
m1_gaus<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                  data = VA_co2,
                  family = gaussian(),
                  dispformula = ~1)
#Evaluate model:
check_convergence(m1_gaus) #OK
check_singularity(m1_gaus,tolerance = 1e-8) #OK
r2(m1_gaus, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m1_gaus)
summary(m1_gaus)
res<- simulateResiduals(m1_gaus)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus)
em<- emmeans(m1_gaus, ~ status)
summary(em)
pairs(em)
DHARMa::testOutliers(res)
#DECISSION: OK-ish


#T_family OPTION, status*season
m2_t<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
               data = VA_co2,
               family = t_family(),
               dispformula = ~1)
#Evaluate basic model:
check_convergence(m2_t)
check_singularity(m2_t,tolerance = 1e-8)
r2(m2_t, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m2_t)
summary(m2_t)
res<- simulateResiduals(m2_t)
plotQQunif(res)
plotResiduals(res)
DHARMa::testOutliers(res)
check_residuals(m2_t)
em<- emmeans(m2_t, ~ status)
summary(em)
pairs(em)
#DECISSION: 


#Compare models (if needed):AIC & pval
anova(m1_gaus,m2_t)


#Set as best
best_VA_co2<- m1_gaus


#Plot obs vs pred 
data.frame(obs=best_VA_co2$frame$dailyflux_trans,
           pred=predict(best_VA_co2, type = "response"),
           status=best_VA_co2$frame$status) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=status))+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "VA_co2 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

# Save obs vs pred plot
ggsave(filename = "coremodels_VA_co2_Observed_VS_predicted.png",
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)


#CH4 models------


##best_CA_ch4-------


#Subset data and check transformation: 
CA_ch4<- data4models %>% filter(casepilot=="CA"&ghgspecies=="ch4")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="CA"&
           ghgspecies=="ch4") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#GAUSSIAN OPTION, status*season
m1_gaus<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                  data = CA_ch4,
                  family = gaussian(),
                  dispformula = ~1)
#Evaluate model:
check_convergence(m1_gaus) #OK
check_singularity(m1_gaus,tolerance = 1e-8) #OK
r2(m1_gaus, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m1_gaus)
summary(m1_gaus)
res<- simulateResiduals(m1_gaus)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus)
em<- emmeans(m1_gaus, ~ status)
summary(em)
pairs(em)
DHARMa::testOutliers(res)
#DECISSION: OK model


#T_family OPTION, status*season
m2_t<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
               data = CA_ch4,
               family = t_family(),
               dispformula = ~1)
#Evaluate basic model:
check_convergence(m2_t)
check_singularity(m2_t,tolerance = 1e-8)
r2(m2_t, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m2_t)
summary(m2_t)
res<- simulateResiduals(m2_t)
plotQQunif(res)
plotResiduals(res)
DHARMa::testOutliers(res)
check_residuals(m2_t)
em<- emmeans(m2_t, ~ status)
summary(em)
pairs(em)
#DECISSION: BAD residuals


#Compare models (if needed):AIC & pval
anova(m1_gaus,m2_t)


#Set as best
best_CA_ch4<- m1_gaus


#Plot obs vs pred 
data.frame(obs=best_CA_ch4$frame$dailyflux_trans,
           pred=predict(best_CA_ch4, type = "response"),
           status=best_CA_ch4$frame$status) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=status))+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "CA_ch4 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

#Save obs vs pred plot
ggsave(filename = "coremodels_CA_ch4_Observed_VS_predicted.png",
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)




##best_CU_ch4 -------

#Subset data and check transformation: 
CU_ch4<- data4models %>% filter(casepilot=="CU"&ghgspecies=="ch4")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="CU"&
           ghgspecies=="ch4") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#GAUSSIAN OPTION, status*season
m1_gaus<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                  data = CU_ch4,
                  family = gaussian(),
                  dispformula = ~1)
#Evaluate model:
check_convergence(m1_gaus) #OK
check_singularity(m1_gaus,tolerance = 1e-8) #OK
r2(m1_gaus, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m1_gaus)
summary(m1_gaus)
res<- simulateResiduals(m1_gaus)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus)
em<- emmeans(m1_gaus, ~ status)
summary(em)
pairs(em)
DHARMa::testOutliers(res)
#DECISSION: GOOD MODEL


#T_family OPTION, status*season
m2_t<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
               data = CU_ch4,
               family = t_family(),
               dispformula = ~1)
#Evaluate basic model:
check_convergence(m2_t)
check_singularity(m2_t,tolerance = 1e-8)
r2(m2_t, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m2_t)
summary(m2_t)
res<- simulateResiduals(m2_t)
plotQQunif(res)
plotResiduals(res)
DHARMa::testOutliers(res)
check_residuals(m2_t)
em<- emmeans(m2_t, ~ status)
summary(em)
pairs(em)
#DECISSION: Not relevant (gaussian is good already)


#Compare models (if needed):AIC & pval
anova(m1_gaus,m2_t)


#Set as best
best_CU_ch4<- m1_gaus


#Plot obs vs pred 
data.frame(obs=best_CU_ch4$frame$dailyflux_trans,
           pred=predict(best_CU_ch4, type = "response"),
           status=best_CU_ch4$frame$status) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=status))+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "CU_ch4 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

#Save obs vs pred plot
ggsave(filename = "coremodels_CU_ch4_Observed_VS_predicted.png",
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)



##best_DA_ch4 -------


#Subset data and check transformation: 
DA_ch4<- data4models %>% filter(casepilot=="DA"&ghgspecies=="ch4")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="DA"&
           ghgspecies=="ch4") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#GAUSSIAN OPTION, status*season
m1_gaus<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                  data = DA_ch4,
                  family = gaussian(),
                  dispformula = ~1)
#Evaluate model:
check_convergence(m1_gaus) #OK
check_singularity(m1_gaus,tolerance = 1e-8) #OK
r2(m1_gaus, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m1_gaus)
summary(m1_gaus)
res<- simulateResiduals(m1_gaus)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus)
em<- emmeans(m1_gaus, ~ status)
summary(em)
pairs(em)
DHARMa::testOutliers(res)
#DECISSION: OK residuals


#T_family OPTION, status*season
m2_t<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
               data = DA_ch4,
               family = t_family(),
               dispformula = ~1)
#Evaluate basic model:
check_convergence(m2_t)
check_singularity(m2_t,tolerance = 1e-8)
r2(m2_t, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m2_t)
summary(m2_t)
res<- simulateResiduals(m2_t)
plotQQunif(res)
plotResiduals(res)
DHARMa::testOutliers(res)
check_residuals(m2_t)
em<- emmeans(m2_t, ~ status)
summary(em)
pairs(em)
#DECISSION: 


#Compare models (if needed):AIC & pval
anova(m1_gaus,m2_t)


#Set as best
best_DA_ch4<- m1_gaus


#Plot obs vs pred 
data.frame(obs=best_DA_ch4$frame$dailyflux_trans,
           pred=predict(best_DA_ch4, type = "response"),
           status=best_DA_ch4$frame$status) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=status))+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "DA_ch4 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

#Save obs vs pred plot
ggsave(filename = "coremodels_DA_ch4_Observed_VS_predicted.png",
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)


##best_DU_ch4 -------

#Subset data and check transformation: 
DU_ch4<- data4models %>% filter(casepilot=="DU"&ghgspecies=="ch4")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="DU"&
           ghgspecies=="ch4") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#GAUSSIAN OPTION, status*season
m1_gaus<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                  data = DU_ch4,
                  family = gaussian(),
                  dispformula = ~1)
#Evaluate model:
check_convergence(m1_gaus) #OK
check_singularity(m1_gaus,tolerance = 1e-8) #OK
r2(m1_gaus, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m1_gaus)
summary(m1_gaus)
res<- simulateResiduals(m1_gaus)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus)
em<- emmeans(m1_gaus, ~ status)
summary(em)
pairs(em)
DHARMa::testOutliers(res)
#DECISSION: OK-ish residuals


#T_family OPTION, status*season
m2_t<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
               data = DU_ch4,
               family = t_family(),
               dispformula = ~1)
#Evaluate basic model:
check_convergence(m2_t)
check_singularity(m2_t,tolerance = 1e-8)
r2(m2_t, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m2_t)
summary(m2_t)
res<- simulateResiduals(m2_t)
plotQQunif(res)
plotResiduals(res)
DHARMa::testOutliers(res)
check_residuals(m2_t)
em<- emmeans(m2_t, ~ status)
summary(em)
pairs(em)
#DECISSION: 


#Compare models (if needed):AIC & pval
anova(m1_gaus,m2_t)


#Set as best
best_DU_ch4<- m2_t


#Plot obs vs pred 
data.frame(obs=best_DU_ch4$frame$dailyflux_trans,
           pred=predict(best_DU_ch4, type = "response"),
           status=best_DU_ch4$frame$status) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=status))+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "DU_ch4 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

#Save obs vs pred plot
ggsave(filename = "coremodels_DU_ch4_Observed_VS_predicted.png",
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)



##best_RI_ch4 -------

#Subset data and check transformation: 
RI_ch4<- data4models %>% filter(casepilot=="RI"&ghgspecies=="ch4")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="RI"&
           ghgspecies=="ch4") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#GAUSSIAN OPTION, status*season
m1_gaus<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                  data = RI_ch4,
                  family = gaussian(),
                  dispformula = ~1)
#Evaluate model:
check_convergence(m1_gaus) #OK
check_singularity(m1_gaus,tolerance = 1e-8) #OK
r2(m1_gaus, tolerance=1e-10) #Decrease tolerance to force R2 calculation when singular
Anova(m1_gaus)
summary(m1_gaus)
res<- simulateResiduals(m1_gaus)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus)
em<- emmeans(m1_gaus, ~ status)
summary(em)
pairs(em)
DHARMa::testOutliers(res)
#DECISSION: good residuals


#T_family OPTION, status*season
m2_t<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
               data = RI_ch4,
               family = t_family(),
               dispformula = ~1)
#Evaluate basic model:
check_convergence(m2_t)
check_singularity(m2_t,tolerance = 1e-8)
r2(m2_t, tolerance=1e-10) #Decrease tolerance to force R2 calculation when singular
Anova(m2_t)
summary(m2_t)
res<- simulateResiduals(m2_t)
plotQQunif(res)
plotResiduals(res)
DHARMa::testOutliers(res)
check_residuals(m2_t)
em<- emmeans(m2_t, ~ status)
summary(em)
pairs(em)
#DECISSION: 


#Compare models (if needed):AIC & pval
anova(m1_gaus,m2_t)


#Set as best
best_RI_ch4<- m1_gaus


#Plot obs vs pred 
data.frame(obs=best_RI_ch4$frame$dailyflux_trans,
           pred=predict(best_RI_ch4, type = "response"),
           status=best_RI_ch4$frame$status) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=status))+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "RI_ch4 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

#Save obs vs pred plot
ggsave(filename = "coremodels_RI_ch4_Observed_VS_predicted.png",
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)


##best_VA_ch4 -------

#Subset data and check transformation: 
VA_ch4<- data4models %>% filter(casepilot=="VA"&ghgspecies=="ch4")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="VA"&
           ghgspecies=="ch4") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#GAUSSIAN OPTION, status*season
m1_gaus<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                  data = VA_ch4,
                  family = gaussian(),
                  dispformula = ~1)
#Evaluate model:
check_convergence(m1_gaus) #OK
check_singularity(m1_gaus,tolerance = 1e-8) #OK
r2(m1_gaus, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m1_gaus)
summary(m1_gaus)
res<- simulateResiduals(m1_gaus)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus)
em<- emmeans(m1_gaus, ~ status)
summary(em)
pairs(em)
DHARMa::testOutliers(res)
#DECISSION: OK-ish resiudals


#T_family OPTION, status*season
m2_t<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
               data = VA_ch4,
               family = t_family(),
               dispformula = ~1)
#Evaluate basic model:
check_convergence(m2_t)
check_singularity(m2_t,tolerance = 1e-8)
r2(m2_t, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m2_t)
summary(m2_t)
res<- simulateResiduals(m2_t)
plotQQunif(res)
plotResiduals(res)
DHARMa::testOutliers(res)
check_residuals(m2_t)
em<- emmeans(m2_t, ~ status)
summary(em)
pairs(em)
#DECISSION: 


#Compare models (if needed):AIC & pval
anova(m1_gaus,m2_t)


#Set as best
best_VA_ch4<- m1_gaus


#Plot obs vs pred 
data.frame(obs=best_VA_ch4$frame$dailyflux_trans,
           pred=predict(best_VA_ch4, type = "response"),
           status=best_VA_ch4$frame$status) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=status))+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "VA_ch4 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

#Save obs vs pred plot
ggsave(filename = "coremodels_VA_ch4_Observed_VS_predicted.png",
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)





  

#N2O models------


##best_CA_n2o-------


#Subset data and check transformation: 
CA_n2o<- data4models %>% filter(casepilot=="CA"&ghgspecies=="n2o")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="CA"&
           ghgspecies=="n2o") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#GAUSSIAN OPTION, status*season
m1_gaus<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                  data = CA_n2o,
                  family = gaussian(),
                  dispformula = ~1)
#Evaluate model:
check_convergence(m1_gaus) #OK
check_singularity(m1_gaus,tolerance = 1e-8) #OK
r2(m1_gaus, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m1_gaus)
summary(m1_gaus)
res<- simulateResiduals(m1_gaus)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus)
em<- emmeans(m1_gaus, ~ status)
summary(em)
pairs(em)
DHARMa::testOutliers(res)
#DECISSION: Bad residuals


#T_family OPTION, status*season
m2_t<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
               data = CA_n2o,
               family = t_family(),
               dispformula = ~1)
#Evaluate basic model:
check_convergence(m2_t)
check_singularity(m2_t,tolerance = 1e-8)
r2(m2_t, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m2_t)
summary(m2_t)
res<- simulateResiduals(m2_t)
plotQQunif(res)
plotResiduals(res)
DHARMa::testOutliers(res)
check_residuals(m2_t)
em<- emmeans(m2_t, ~ status)
summary(em)
pairs(em)
#DECISSION: OK residuals


#Compare models (if needed):AIC & pval
anova(m1_gaus,m2_t)


#Set as best
best_CA_n2o<- m2_t


#Plot obs vs pred 
data.frame(obs=best_CA_n2o$frame$dailyflux_trans,
           pred=predict(best_CA_n2o, type = "response"),
           status=best_CA_n2o$frame$status) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=status))+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "CA_n2o Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

#Save obs vs pred plot
ggsave(filename = "coremodels_CA_n2o_Observed_VS_predicted.png",
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)




##best_CU_n2o -------

#Subset data and check transformation: 
CU_n2o<- data4models %>% filter(casepilot=="CU"&ghgspecies=="n2o")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="CU"&
           ghgspecies=="n2o") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#GAUSSIAN OPTION, status*season
m1_gaus<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                  data = CU_n2o,
                  family = gaussian(),
                  dispformula = ~1)
#Evaluate model:
check_convergence(m1_gaus) #OK
check_singularity(m1_gaus,tolerance = 1e-8) #OK
r2(m1_gaus, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m1_gaus)
summary(m1_gaus)
res<- simulateResiduals(m1_gaus)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus)
em<- emmeans(m1_gaus, ~ status)
summary(em)
pairs(em)
DHARMa::testOutliers(res)
#DECISSION: GOOD MODEL


#T_family OPTION, status*season
m2_t<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
               data = CU_n2o,
               family = t_family(),
               dispformula = ~1)
#Evaluate basic model:
check_convergence(m2_t)
check_singularity(m2_t,tolerance = 1e-8)
r2(m2_t, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m2_t)
summary(m2_t)
res<- simulateResiduals(m2_t)
plotQQunif(res)
plotResiduals(res)
DHARMa::testOutliers(res)
check_residuals(m2_t)
em<- emmeans(m2_t, ~ status)
summary(em)
pairs(em)
#DECISSION: Good residuals


#Compare models (if needed):AIC & pval
anova(m1_gaus,m2_t)


#Set as best
best_CU_n2o<- m1_gaus


#Plot obs vs pred 
data.frame(obs=best_CU_n2o$frame$dailyflux_trans,
           pred=predict(best_CU_n2o, type = "response"),
           status=best_CU_n2o$frame$status) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=status))+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "CU_n2o Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

#Save obs vs pred plot
ggsave(filename = "coremodels_CU_n2o_Observed_VS_predicted.png",
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)



##best_DA_n2o -------


#Subset data and check transformation: 
DA_n2o<- data4models %>% filter(casepilot=="DA"&ghgspecies=="n2o")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="DA"&
           ghgspecies=="n2o") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#GAUSSIAN OPTION, status*season
m1_gaus<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                  data = DA_n2o,
                  family = gaussian(),
                  dispformula = ~1)
#Evaluate model:
check_convergence(m1_gaus) #OK
check_singularity(m1_gaus,tolerance = 1e-8) #OK
r2(m1_gaus, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m1_gaus)
summary(m1_gaus)
res<- simulateResiduals(m1_gaus)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus)
em<- emmeans(m1_gaus, ~ status)
summary(em)
pairs(em)
DHARMa::testOutliers(res)
#DECISSION: BADresiduals


#T_family OPTION, status*season
m2_t<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
               data = DA_n2o,
               family = t_family(),
               dispformula = ~1)
#Evaluate basic model:
check_convergence(m2_t)
check_singularity(m2_t,tolerance = 1e-8)
r2(m2_t, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m2_t)
summary(m2_t)
res<- simulateResiduals(m2_t)
plotQQunif(res)
plotResiduals(res)
DHARMa::testOutliers(res)
check_residuals(m2_t)
em<- emmeans(m2_t, ~ status)
summary(em)
pairs(em)
#DECISSION: OKish residuals


#Compare models (if needed):AIC & pval
anova(m1_gaus,m2_t)


#Set as best
best_DA_n2o<- m2_t


#Plot obs vs pred 
data.frame(obs=best_DA_n2o$frame$dailyflux_trans,
           pred=predict(best_DA_n2o, type = "response"),
           status=best_DA_n2o$frame$status) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=status))+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "DA_n2o Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

#Save obs vs pred plot
ggsave(filename = "coremodels_DA_n2o_Observed_VS_predicted.png",
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)


##best_DU_n2o -------

#Subset data and check transformation: 
DU_n2o<- data4models %>% filter(casepilot=="DU"&ghgspecies=="n2o")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="DU"&
           ghgspecies=="n2o") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#GAUSSIAN OPTION, status*season
m1_gaus<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                  data = DU_n2o,
                  family = gaussian(),
                  dispformula = ~1)
#Evaluate model:
check_convergence(m1_gaus) #OK
check_singularity(m1_gaus,tolerance = 1e-8) #OK
r2(m1_gaus, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m1_gaus)
summary(m1_gaus)
res<- simulateResiduals(m1_gaus)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus)
em<- emmeans(m1_gaus, ~ status)
summary(em)
pairs(em)
DHARMa::testOutliers(res)
#DECISSION: OK-ish residuals


#T_family OPTION, status*season
m2_t<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
               data = DU_n2o,
               family = t_family(),
               dispformula = ~1)
#Evaluate basic model:
check_convergence(m2_t)
check_singularity(m2_t,tolerance = 1e-8)
r2(m2_t, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m2_t)
summary(m2_t)
res<- simulateResiduals(m2_t)
plotQQunif(res)
plotResiduals(res)
DHARMa::testOutliers(res)
check_residuals(m2_t)
em<- emmeans(m2_t, ~ status)
summary(em)
pairs(em)
#DECISSION: 


#Compare models (if needed):AIC & pval
anova(m1_gaus,m2_t)


#Set as best
best_DU_n2o<- m2_t


#Plot obs vs pred 
data.frame(obs=best_DU_n2o$frame$dailyflux_trans,
           pred=predict(best_DU_n2o, type = "response"),
           status=best_DU_n2o$frame$status) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=status))+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "DU_n2o Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

#Save obs vs pred plot
ggsave(filename = "coremodels_DU_n2o_Observed_VS_predicted.png",
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)



##best_RI_n2o -------

#Subset data and check transformation: 
RI_n2o<- data4models %>% filter(casepilot=="RI"&ghgspecies=="n2o")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="RI"&
           ghgspecies=="n2o") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#GAUSSIAN OPTION, status*season
m1_gaus<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                  data = RI_n2o,
                  family = gaussian(),
                  dispformula = ~1)
#Evaluate model:
check_convergence(m1_gaus) #OK
check_singularity(m1_gaus,tolerance = 1e-8) #OK
r2(m1_gaus, tolerance=1e-10) #Decrease tolerance to force R2 calculation when singular
Anova(m1_gaus)
summary(m1_gaus)
res<- simulateResiduals(m1_gaus)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus)
em<- emmeans(m1_gaus, ~ status)
summary(em)
pairs(em)
DHARMa::testOutliers(res)
#DECISSION: OK-ish residuals


#T_family OPTION, status*season
m2_t<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
               data = RI_n2o,
               family = t_family(),
               dispformula = ~1)
#Evaluate basic model:
check_convergence(m2_t)
check_singularity(m2_t,tolerance = 1e-8)
r2(m2_t, tolerance=1e-10) #Decrease tolerance to force R2 calculation when singular
Anova(m2_t)
summary(m2_t)
res<- simulateResiduals(m2_t)
plotQQunif(res)
plotResiduals(res)
DHARMa::testOutliers(res)
check_residuals(m2_t)
em<- emmeans(m2_t, ~ status)
summary(em)
pairs(em)
#DECISSION: 


#Compare models (if needed):AIC & pval
anova(m1_gaus,m2_t)


#Set as best
best_RI_n2o<- m1_gaus


#Plot obs vs pred 
data.frame(obs=best_RI_n2o$frame$dailyflux_trans,
           pred=predict(best_RI_n2o, type = "response"),
           status=best_RI_n2o$frame$status) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=status))+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "RI_n2o Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

#Save obs vs pred plot
ggsave(filename = "coremodels_RI_n2o_Observed_VS_predicted.png",
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)


##best_VA_n2o -------

#Subset data and check transformation: 
VA_n2o<- data4models %>% filter(casepilot=="VA"&ghgspecies=="n2o")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="VA"&
           ghgspecies=="n2o") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#GAUSSIAN OPTION, status*season
m1_gaus<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                  data = VA_n2o,
                  family = gaussian(),
                  dispformula = ~1)
#Evaluate model:
check_convergence(m1_gaus) #OK
check_singularity(m1_gaus,tolerance = 1e-8) #OK
r2(m1_gaus, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m1_gaus)
summary(m1_gaus)
res<- simulateResiduals(m1_gaus)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus)
em<- emmeans(m1_gaus, ~ status)
summary(em)
pairs(em)
DHARMa::testOutliers(res)
#DECISSION: OK-ish resiudals


#T_family OPTION, status*season
m2_t<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
               data = VA_n2o,
               family = t_family(),
               dispformula = ~1)
#Evaluate basic model:
check_convergence(m2_t)
check_singularity(m2_t,tolerance = 1e-8)
r2(m2_t, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m2_t)
summary(m2_t)
res<- simulateResiduals(m2_t)
plotQQunif(res)
plotResiduals(res)
DHARMa::testOutliers(res)
check_residuals(m2_t)
em<- emmeans(m2_t, ~ status)
summary(em)
pairs(em)
#DECISSION: OK-ish residuals


#Compare models (if needed):AIC & pval
anova(m1_gaus,m2_t)


#Set as best
best_VA_n2o<- m2_t


#Plot obs vs pred 
data.frame(obs=best_VA_n2o$frame$dailyflux_trans,
           pred=predict(best_VA_n2o, type = "response"),
           status=best_VA_n2o$frame$status) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=status))+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "VA_n2o Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

#Save obs vs pred plot
ggsave(filename = "coremodels_VA_n2o_Observed_VS_predicted.png",
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)



#GWPco2andch4 models------


##best_CA_GWPco2andch4-------


#Subset data and check transformation: 
CA_GWPco2andch4<- data4models %>% filter(casepilot=="CA"&ghgspecies=="GWPco2andch4")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="CA"&
           ghgspecies=="GWPco2andch4") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#GAUSSIAN OPTION, status*season
m1_gaus<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                  data = CA_GWPco2andch4,
                  family = gaussian(),
                  dispformula = ~1)
#Evaluate model:
check_convergence(m1_gaus) #OK
check_singularity(m1_gaus,tolerance = 1e-8) #OK
r2(m1_gaus, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m1_gaus)
summary(m1_gaus)
res<- simulateResiduals(m1_gaus)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus)
em<- emmeans(m1_gaus, ~ status)
summary(em)
pairs(em)
DHARMa::testOutliers(res)
#DECISSION: OK residuals


#T_family OPTION, status*season
m2_t<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
               data = CA_GWPco2andch4,
               family = t_family(),
               dispformula = ~1)
#Evaluate basic model:
check_convergence(m2_t)
check_singularity(m2_t,tolerance = 1e-8)
r2(m2_t, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m2_t)
summary(m2_t)
res<- simulateResiduals(m2_t)
plotQQunif(res)
plotResiduals(res)
DHARMa::testOutliers(res)
check_residuals(m2_t)
em<- emmeans(m2_t, ~ status)
summary(em)
pairs(em)
#DECISSION: OK residuals


#Compare models (if needed):AIC & pval
anova(m1_gaus,m2_t)


#Set as best
best_CA_GWPco2andch4<- m2_t


#Plot obs vs pred 
data.frame(obs=best_CA_GWPco2andch4$frame$dailyflux_trans,
           pred=predict(best_CA_GWPco2andch4, type = "response"),
           status=best_CA_GWPco2andch4$frame$status) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=status))+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "CA_GWPco2andch4 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

#Save obs vs pred plot
ggsave(filename = "coremodels_CA_GWPco2andch4_Observed_VS_predicted.png",
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)




##best_CU_GWPco2andch4 -------

#Subset data and check transformation: 
CU_GWPco2andch4<- data4models %>% filter(casepilot=="CU"&ghgspecies=="GWPco2andch4")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="CU"&
           ghgspecies=="GWPco2andch4") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#GAUSSIAN OPTION, status*season
m1_gaus<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                  data = CU_GWPco2andch4,
                  family = gaussian(),
                  dispformula = ~1)
#Evaluate model:
check_convergence(m1_gaus) #OK
check_singularity(m1_gaus,tolerance = 1e-8) #OK
r2(m1_gaus, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m1_gaus)
summary(m1_gaus)
res<- simulateResiduals(m1_gaus)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus)
em<- emmeans(m1_gaus, ~ status)
summary(em)
pairs(em)
DHARMa::testOutliers(res)
#DECISSION: GOOD MODEL


#T_family OPTION, status*season
m2_t<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
               data = CU_GWPco2andch4,
               family = t_family(),
               dispformula = ~1)
#Evaluate basic model:
check_convergence(m2_t)
check_singularity(m2_t,tolerance = 1e-8)
r2(m2_t, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m2_t)
summary(m2_t)
res<- simulateResiduals(m2_t)
plotQQunif(res)
plotResiduals(res)
DHARMa::testOutliers(res)
check_residuals(m2_t)
em<- emmeans(m2_t, ~ status)
summary(em)
pairs(em)
#DECISSION: Good residuals


#Compare models (if needed):AIC & pval
anova(m1_gaus,m2_t)


#Set as best
best_CU_GWPco2andch4<- m1_gaus


#Plot obs vs pred 
data.frame(obs=best_CU_GWPco2andch4$frame$dailyflux_trans,
           pred=predict(best_CU_GWPco2andch4, type = "response"),
           status=best_CU_GWPco2andch4$frame$status) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=status))+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "CU_GWPco2andch4 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

#Save obs vs pred plot
ggsave(filename = "coremodels_CU_GWPco2andch4_Observed_VS_predicted.png",
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)



##best_DA_GWPco2andch4 -------


#Subset data and check transformation: 
DA_GWPco2andch4<- data4models %>% filter(casepilot=="DA"&ghgspecies=="GWPco2andch4")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="DA"&
           ghgspecies=="GWPco2andch4") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#GAUSSIAN OPTION, status*season
m1_gaus<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                  data = DA_GWPco2andch4,
                  family = gaussian(),
                  dispformula = ~1)
#Evaluate model:
check_convergence(m1_gaus) #OK
check_singularity(m1_gaus,tolerance = 1e-8) #OK
r2(m1_gaus, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m1_gaus)
summary(m1_gaus)
res<- simulateResiduals(m1_gaus)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus)
em<- emmeans(m1_gaus, ~ status)
summary(em)
pairs(em)
DHARMa::testOutliers(res)
#DECISSION: OK-ish 


#T_family OPTION, status*season
m2_t<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
               data = DA_GWPco2andch4,
               family = t_family(),
               dispformula = ~1)
#Evaluate basic model:
check_convergence(m2_t)
check_singularity(m2_t,tolerance = 1e-8)
r2(m2_t, tolerance=1e-10) #Decrease tolerance to force R2 calculation when singular
Anova(m2_t)
summary(m2_t)
res<- simulateResiduals(m2_t)
plotQQunif(res)
plotResiduals(res)
DHARMa::testOutliers(res)
check_residuals(m2_t)
em<- emmeans(m2_t, ~ status)
summary(em)
pairs(em)
#DECISSION: BAD residuals


#Compare models (if needed):AIC & pval
anova(m1_gaus,m2_t)


#Set as best
best_DA_GWPco2andch4<- m1_gaus


#Plot obs vs pred 
data.frame(obs=best_DA_GWPco2andch4$frame$dailyflux_trans,
           pred=predict(best_DA_GWPco2andch4, type = "response"),
           status=best_DA_GWPco2andch4$frame$status) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=status))+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "DA_GWPco2andch4 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

#Save obs vs pred plot
ggsave(filename = "coremodels_DA_GWPco2andch4_Observed_VS_predicted.png",
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)


##best_DU_GWPco2andch4 -------

#Subset data and check transformation: 
DU_GWPco2andch4<- data4models %>% filter(casepilot=="DU"&ghgspecies=="GWPco2andch4")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="DU"&
           ghgspecies=="GWPco2andch4") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#GAUSSIAN OPTION, status*season
m1_gaus<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                  data = DU_GWPco2andch4,
                  family = gaussian(),
                  dispformula = ~1)
#Evaluate model:
check_convergence(m1_gaus) #OK
check_singularity(m1_gaus,tolerance = 1e-8) #OK
r2(m1_gaus, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m1_gaus)
summary(m1_gaus)
res<- simulateResiduals(m1_gaus)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus)
em<- emmeans(m1_gaus, ~ status)
summary(em)
pairs(em)
DHARMa::testOutliers(res)
#DECISSION: OK-ish residuals


#T_family OPTION, status*season
m2_t<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
               data = DU_GWPco2andch4,
               family = t_family(),
               dispformula = ~1)
#Evaluate basic model:
check_convergence(m2_t)
check_singularity(m2_t,tolerance = 1e-8)
r2(m2_t, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m2_t)
summary(m2_t)
res<- simulateResiduals(m2_t)
plotQQunif(res)
plotResiduals(res)
DHARMa::testOutliers(res)
check_residuals(m2_t)
em<- emmeans(m2_t, ~ status)
summary(em)
pairs(em)
#DECISSION: non-appropriate, gaussian is ok


#Compare models (if needed):AIC & pval
anova(m1_gaus,m2_t)


#Set as best
best_DU_GWPco2andch4<- m1_gaus


#Plot obs vs pred 
data.frame(obs=best_DU_GWPco2andch4$frame$dailyflux_trans,
           pred=predict(best_DU_GWPco2andch4, type = "response"),
           status=best_DU_GWPco2andch4$frame$status) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=status))+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "DU_GWPco2andch4 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

#Save obs vs pred plot
ggsave(filename = "coremodels_DU_GWPco2andch4_Observed_VS_predicted.png",
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)



##best_RI_GWPco2andch4 -------

#Subset data and check transformation: 
RI_GWPco2andch4<- data4models %>% filter(casepilot=="RI"&ghgspecies=="GWPco2andch4")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="RI"&
           ghgspecies=="GWPco2andch4") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#GAUSSIAN OPTION, status*season
m1_gaus<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                  data = RI_GWPco2andch4,
                  family = gaussian(),
                  dispformula = ~1)
#Evaluate model:
check_convergence(m1_gaus) #OK
check_singularity(m1_gaus,tolerance = 1e-8) #OK
r2(m1_gaus, tolerance=1e-10) #Decrease tolerance to force R2 calculation when singular
Anova(m1_gaus)
summary(m1_gaus)
res<- simulateResiduals(m1_gaus)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus)
em<- emmeans(m1_gaus, ~ status)
summary(em)
pairs(em)
DHARMa::testOutliers(res)
#DECISSION: OK-ish residuals


#T_family OPTION, status*season
m2_t<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
               data = RI_GWPco2andch4,
               family = t_family(),
               dispformula = ~1)
#Evaluate basic model:
check_convergence(m2_t)
check_singularity(m2_t,tolerance = 1e-8)
r2(m2_t, tolerance=1e-11) #Decrease tolerance to force R2 calculation when singular
Anova(m2_t)
summary(m2_t)
res<- simulateResiduals(m2_t)
plotQQunif(res)
plotResiduals(res)
DHARMa::testOutliers(res)
check_residuals(m2_t)
em<- emmeans(m2_t, ~ status)
summary(em)
pairs(em)
#DECISSION: 


#Compare models (if needed):AIC & pval
anova(m1_gaus,m2_t)


#Set as best
best_RI_GWPco2andch4<- m2_t


#Plot obs vs pred 
data.frame(obs=best_RI_GWPco2andch4$frame$dailyflux_trans,
           pred=predict(best_RI_GWPco2andch4, type = "response"),
           status=best_RI_GWPco2andch4$frame$status) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=status))+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "RI_GWPco2andch4 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

#Save obs vs pred plot
ggsave(filename = "coremodels_RI_GWPco2andch4_Observed_VS_predicted.png",
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)


##best_VA_GWPco2andch4 -------

#Subset data and check transformation: 
VA_GWPco2andch4<- data4models %>% filter(casepilot=="VA"&ghgspecies=="GWPco2andch4")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="VA"&
           ghgspecies=="GWPco2andch4") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#GAUSSIAN OPTION, status*season
m1_gaus<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                  data = VA_GWPco2andch4,
                  family = gaussian(),
                  dispformula = ~1)
#Evaluate model:
check_convergence(m1_gaus) #OK
check_singularity(m1_gaus,tolerance = 1e-8) #OK
r2(m1_gaus, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m1_gaus)
summary(m1_gaus)
res<- simulateResiduals(m1_gaus)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus)
em<- emmeans(m1_gaus, ~ status)
summary(em)
pairs(em)
DHARMa::testOutliers(res)
#DECISSION: OK-ish resiudals


#T_family OPTION, status*season
m2_t<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
               data = VA_GWPco2andch4,
               family = t_family(),
               dispformula = ~1)
#Evaluate basic model:
check_convergence(m2_t)
check_singularity(m2_t,tolerance = 1e-8)
r2(m2_t, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m2_t)
summary(m2_t)
res<- simulateResiduals(m2_t)
plotQQunif(res)
plotResiduals(res)
DHARMa::testOutliers(res)
check_residuals(m2_t)
em<- emmeans(m2_t, ~ status)
summary(em)
pairs(em)
#DECISSION: OK-ish residuals


#Compare models (if needed):AIC & pval
anova(m1_gaus,m2_t)


#Set as best
best_VA_GWPco2andch4<- m1_gaus


#Plot obs vs pred 
data.frame(obs=best_VA_GWPco2andch4$frame$dailyflux_trans,
           pred=predict(best_VA_GWPco2andch4, type = "response"),
           status=best_VA_GWPco2andch4$frame$status) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=status))+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "VA_GWPco2andch4 Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

#Save obs vs pred plot
ggsave(filename = "coremodels_VA_GWPco2andch4_Observed_VS_predicted.png",
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)




#GWPtotal models------


##best_CA_GWPtotal-------


#Subset data and check transformation: 
CA_GWPtotal<- data4models %>% filter(casepilot=="CA"&ghgspecies=="GWPtotal")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="CA"&
           ghgspecies=="GWPtotal") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#GAUSSIAN OPTION, status*season
m1_gaus<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                  data = CA_GWPtotal,
                  family = gaussian(),
                  dispformula = ~1)
#Evaluate model:
check_convergence(m1_gaus) #OK
check_singularity(m1_gaus,tolerance = 1e-8) #OK
r2(m1_gaus, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m1_gaus)
summary(m1_gaus)
res<- simulateResiduals(m1_gaus)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus)
em<- emmeans(m1_gaus, ~ status)
summary(em)
pairs(em)
DHARMa::testOutliers(res)
#DECISSION: OK residuals


#T_family OPTION, status*season
m2_t<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
               data = CA_GWPtotal,
               family = t_family(),
               dispformula = ~1)
#Evaluate basic model:
check_convergence(m2_t)
check_singularity(m2_t,tolerance = 1e-8)
r2(m2_t, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m2_t)
summary(m2_t)
res<- simulateResiduals(m2_t)
plotQQunif(res)
plotResiduals(res)
DHARMa::testOutliers(res)
check_residuals(m2_t)
em<- emmeans(m2_t, ~ status)
summary(em)
pairs(em)
#DECISSION: OK residuals


#Compare models (if needed):AIC & pval
anova(m1_gaus,m2_t)


#Set as best
best_CA_GWPtotal<- m2_t


#Plot obs vs pred 
data.frame(obs=best_CA_GWPtotal$frame$dailyflux_trans,
           pred=predict(best_CA_GWPtotal, type = "response"),
           status=best_CA_GWPtotal$frame$status) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=status))+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "CA_GWPtotal Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

#Save obs vs pred plot
ggsave(filename = "coremodels_CA_GWPtotal_Observed_VS_predicted.png",
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)




##best_CU_GWPtotal -------

#Subset data and check transformation: 
CU_GWPtotal<- data4models %>% filter(casepilot=="CU"&ghgspecies=="GWPtotal")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="CU"&
           ghgspecies=="GWPtotal") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#GAUSSIAN OPTION, status*season
m1_gaus<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                  data = CU_GWPtotal,
                  family = gaussian(),
                  dispformula = ~1)
#Evaluate model:
check_convergence(m1_gaus) #OK
check_singularity(m1_gaus,tolerance = 1e-8) #OK
r2(m1_gaus, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m1_gaus)
summary(m1_gaus)
res<- simulateResiduals(m1_gaus)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus)
em<- emmeans(m1_gaus, ~ status)
summary(em)
pairs(em)
DHARMa::testOutliers(res)
#DECISSION: GOOD MODEL


#T_family OPTION, status*season
m2_t<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
               data = CU_GWPtotal,
               family = t_family(),
               dispformula = ~1)
#Evaluate basic model:
check_convergence(m2_t)
check_singularity(m2_t,tolerance = 1e-8)
r2(m2_t, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m2_t)
summary(m2_t)
res<- simulateResiduals(m2_t)
plotQQunif(res)
plotResiduals(res)
DHARMa::testOutliers(res)
check_residuals(m2_t)
em<- emmeans(m2_t, ~ status)
summary(em)
pairs(em)
#DECISSION: Good residuals


#Compare models (if needed):AIC & pval
anova(m1_gaus,m2_t)


#Set as best
best_CU_GWPtotal<- m1_gaus


#Plot obs vs pred 
data.frame(obs=best_CU_GWPtotal$frame$dailyflux_trans,
           pred=predict(best_CU_GWPtotal, type = "response"),
           status=best_CU_GWPtotal$frame$status) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=status))+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "CU_GWPtotal Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

#Save obs vs pred plot
ggsave(filename = "coremodels_CU_GWPtotal_Observed_VS_predicted.png",
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)



##best_DA_GWPtotal -------


#Subset data and check transformation: 
DA_GWPtotal<- data4models %>% filter(casepilot=="DA"&ghgspecies=="GWPtotal")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="DA"&
           ghgspecies=="GWPtotal") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#GAUSSIAN OPTION, status*season
m1_gaus<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                  data = DA_GWPtotal,
                  family = gaussian(),
                  dispformula = ~1)
#Evaluate model:
check_convergence(m1_gaus) #OK
check_singularity(m1_gaus,tolerance = 1e-8) #OK
r2(m1_gaus, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m1_gaus)
summary(m1_gaus)
res<- simulateResiduals(m1_gaus)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus)
em<- emmeans(m1_gaus, ~ status)
summary(em)
pairs(em)
DHARMa::testOutliers(res)
#DECISSION: OK residuals


#T_family OPTION, status*season
m2_t<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
               data = DA_GWPtotal,
               family = t_family(),
               dispformula = ~1)
#Evaluate basic model:
check_convergence(m2_t)
check_singularity(m2_t,tolerance = 1e-8)
r2(m2_t, tolerance=1e-10) #Decrease tolerance to force R2 calculation when singular
Anova(m2_t)
summary(m2_t)
res<- simulateResiduals(m2_t)
plotQQunif(res)
plotResiduals(res)
DHARMa::testOutliers(res)
check_residuals(m2_t)
em<- emmeans(m2_t, ~ status)
summary(em)
pairs(em)
#DECISSION: BAD residuals


#Compare models (if needed):AIC & pval
anova(m1_gaus,m2_t)


#Set as best
best_DA_GWPtotal<- m1_gaus


#Plot obs vs pred 
data.frame(obs=best_DA_GWPtotal$frame$dailyflux_trans,
           pred=predict(best_DA_GWPtotal, type = "response"),
           status=best_DA_GWPtotal$frame$status) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=status))+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "DA_GWPtotal Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

#Save obs vs pred plot
ggsave(filename = "coremodels_DA_GWPtotal_Observed_VS_predicted.png",
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)


##best_DU_GWPtotal -------

#Subset data and check transformation: 
DU_GWPtotal<- data4models %>% filter(casepilot=="DU"&ghgspecies=="GWPtotal")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="DU"&
           ghgspecies=="GWPtotal") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#GAUSSIAN OPTION, status*season
m1_gaus<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                  data = DU_GWPtotal,
                  family = gaussian(),
                  dispformula = ~1)
#Evaluate model:
check_convergence(m1_gaus) #OK
check_singularity(m1_gaus,tolerance = 1e-8) #OK
r2(m1_gaus, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m1_gaus)
summary(m1_gaus)
res<- simulateResiduals(m1_gaus)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus)
em<- emmeans(m1_gaus, ~ status)
summary(em)
pairs(em)
DHARMa::testOutliers(res)
#DECISSION: OK-ish residuals


#T_family OPTION, status*season
m2_t<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
               data = DU_GWPtotal,
               family = t_family(),
               dispformula = ~1)
#Evaluate basic model:
check_convergence(m2_t)
check_singularity(m2_t,tolerance = 1e-8)
r2(m2_t, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m2_t)
summary(m2_t)
res<- simulateResiduals(m2_t)
plotQQunif(res)
plotResiduals(res)
DHARMa::testOutliers(res)
check_residuals(m2_t)
em<- emmeans(m2_t, ~ status)
summary(em)
pairs(em)
#DECISSION: OK residuals


#Compare models (if needed):AIC & pval
anova(m1_gaus,m2_t)


#Set as best
best_DU_GWPtotal<- m1_gaus


#Plot obs vs pred 
data.frame(obs=best_DU_GWPtotal$frame$dailyflux_trans,
           pred=predict(best_DU_GWPtotal, type = "response"),
           status=best_DU_GWPtotal$frame$status) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=status))+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "DU_GWPtotal Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

#Save obs vs pred plot
ggsave(filename = "coremodels_DU_GWPtotal_Observed_VS_predicted.png",
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)



##best_RI_GWPtotal -------

#Subset data and check transformation: 
RI_GWPtotal<- data4models %>% filter(casepilot=="RI"&ghgspecies=="GWPtotal")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="RI"&
           ghgspecies=="GWPtotal") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#GAUSSIAN OPTION, status*season
m1_gaus<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                  data = RI_GWPtotal,
                  family = gaussian(),
                  dispformula = ~1)
#Evaluate model:
check_convergence(m1_gaus) #OK
check_singularity(m1_gaus,tolerance = 1e-8) #OK
r2(m1_gaus, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m1_gaus)
summary(m1_gaus)
res<- simulateResiduals(m1_gaus)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus)
em<- emmeans(m1_gaus, ~ status)
summary(em)
pairs(em)
DHARMa::testOutliers(res)
#DECISSION: bad residuals


#T_family OPTION, status*season
m2_t<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
               data = RI_GWPtotal,
               family = t_family(),
               dispformula = ~1)
#Evaluate basic model:
check_convergence(m2_t)
check_singularity(m2_t,tolerance = 1e-8)
r2(m2_t, tolerance=1e-11) #Decrease tolerance to force R2 calculation when singular
Anova(m2_t)
summary(m2_t)
res<- simulateResiduals(m2_t)
plotQQunif(res)
plotResiduals(res)
DHARMa::testOutliers(res)
check_residuals(m2_t)
em<- emmeans(m2_t, ~ status)
summary(em)
pairs(em)
#DECISSION: Ok residuals


#Compare models (if needed):AIC & pval
anova(m1_gaus,m2_t)


#Set as best
best_RI_GWPtotal<- m2_t


#Plot obs vs pred 
data.frame(obs=best_RI_GWPtotal$frame$dailyflux_trans,
           pred=predict(best_RI_GWPtotal, type = "response"),
           status=best_RI_GWPtotal$frame$status) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=status))+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "RI_GWPtotal Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

#Save obs vs pred plot
ggsave(filename = "coremodels_RI_GWPtotal_Observed_VS_predicted.png",
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)


##best_VA_GWPtotal -------

#Subset data and check transformation: 
VA_GWPtotal<- data4models %>% filter(casepilot=="VA"&ghgspecies=="GWPtotal")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="VA"&
           ghgspecies=="GWPtotal") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#GAUSSIAN OPTION, status*season
m1_gaus<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                  data = VA_GWPtotal,
                  family = gaussian(),
                  dispformula = ~1)
#Evaluate model:
check_convergence(m1_gaus) #OK
check_singularity(m1_gaus,tolerance = 1e-8) #OK
r2(m1_gaus, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m1_gaus)
summary(m1_gaus)
res<- simulateResiduals(m1_gaus)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus)
em<- emmeans(m1_gaus, ~ status)
summary(em)
pairs(em)
DHARMa::testOutliers(res)
#DECISSION: OK-ish resiudals


#T_family OPTION, status*season
m2_t<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
               data = VA_GWPtotal,
               family = t_family(),
               dispformula = ~1)
#Evaluate basic model:
check_convergence(m2_t)
check_singularity(m2_t,tolerance = 1e-8)
r2(m2_t, tolerance=1e-8) #Decrease tolerance to force R2 calculation when singular
Anova(m2_t)
summary(m2_t)
res<- simulateResiduals(m2_t)
plotQQunif(res)
plotResiduals(res)
DHARMa::testOutliers(res)
check_residuals(m2_t)
em<- emmeans(m2_t, ~ status)
summary(em)
pairs(em)
#DECISSION: OK-ish residuals


#Compare models (if needed):AIC & pval
anova(m1_gaus,m2_t)


#Set as best
best_VA_GWPtotal<- m2_t


#Plot obs vs pred 
data.frame(obs=best_VA_GWPtotal$frame$dailyflux_trans,
           pred=predict(best_VA_GWPtotal, type = "response"),
           status=best_VA_GWPtotal$frame$status) %>% 
  ggplot(aes(x=pred, y=obs))+
  geom_point(aes(col=status))+
  stat_cor()+
  geom_smooth(method="lm", col="black", se=F)+
  geom_abline(intercept = 0, slope = 1, col="red",linetype = 2,linewidth = 1)+
  labs(title = "VA_GWPtotal Observed vs Predicted",
       subtitle = "Best-fit lines and correations in model scale (transformed)")+
  theme_bw()+
  theme(plot.caption = element_textbox_simple())

#Save obs vs pred plot
ggsave(filename = "coremodels_VA_GWPtotal_Observed_VS_predicted.png",
       width = 115,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)



#3. Save bestmodel-list -------

#Save list with best-models, named as casepilot_ghgspecies (needed in order to work appropriately with the summarising custom functions)

#List of best models for each casepilot
cores_best_model_list_allghg<- list(
  "CA_co2" = best_CA_co2,
  "CU_co2" = best_CU_co2,
  "DA_co2" = best_DA_co2,
  "DU_co2" = best_DU_co2,
  "RI_co2" = best_RI_co2,
  "VA_co2" = best_VA_co2,
  "CA_ch4" = best_CA_ch4,
  "CU_ch4" = best_CU_ch4,
  "DA_ch4" = best_DA_ch4,
  "DU_ch4" = best_DU_ch4,
  "RI_ch4" = best_RI_ch4,
  "VA_ch4" = best_VA_ch4,
  "CA_n2o" = best_CA_n2o,
  "CU_n2o" = best_CU_n2o,
  "DA_n2o" = best_DA_n2o,
  "DU_n2o" = best_DU_n2o,
  "RI_n2o" = best_RI_n2o,
  "VA_n2o" = best_VA_n2o,
  "CA_GWPco2andch4" = best_CA_GWPco2andch4,
  "CU_GWPco2andch4" = best_CU_GWPco2andch4,
  "DA_GWPco2andch4" = best_DA_GWPco2andch4,
  "DU_GWPco2andch4" = best_DU_GWPco2andch4,
  "RI_GWPco2andch4" = best_RI_GWPco2andch4,
  "VA_GWPco2andch4" = best_VA_GWPco2andch4,
  "CA_GWPtotal" = best_CA_GWPtotal,
  "CU_GWPtotal" = best_CU_GWPtotal,
  "DA_GWPtotal" = best_DA_GWPtotal,
  "DU_GWPtotal" = best_DU_GWPtotal,
  "RI_GWPtotal" = best_RI_GWPtotal,
  "VA_GWPtotal" = best_VA_GWPtotal
)



## 3.1. Overall fit -----
#Structure, and pseudoR2s of best model:
all_best_model_fit<- get_pseudoR2s(cores_best_model_list_allghg)


#R2c and R2m of best model:
all_homoc_r2<- get_homoc_R2s(cores_best_model_list_allghg)
print(all_homoc_r2)

#DHARMa Residual diagnostics of best model:
all_resid_diag_summary <- summarize_dharma_diagnostics(cores_best_model_list_allghg)
print(all_resid_diag_summary)


##3.2. Significance of effects-----
#via Anova (best_model,type=3): 
#intercept pvalue: across all levels of effects, is the flux_trans different from cero? yes/no
#status pvalue: does status have an effect across all seasons? yes/no, 
#season pvalue: does season have an effect across all status? yes/no
#interaction pvalue: does the effect of status depend on the season? yes/no

#Obtain significance for main effects of best_models
all_results_anova<-get_anova_results(cores_best_model_list_allghg) 
print(all_results_anova)


#Combine model summary: 
all_model_outputs<- all_results_anova %>% 
  full_join(all_best_model_fit) %>% 
  full_join(all_homoc_r2 %>% dplyr::select(dataset, homoced_R2m, homoced_R2c))


#Save in dropbox tables with the models summary outputs: 
write.csv(x = all_model_outputs, file = paste0(paper_path,"Allghg_summary_coresmodels.csv"),row.names = F)

write.csv(x = all_resid_diag_summary, file = paste0(paper_path,"Allghg_Residualtests_coresmodels.csv"),row.names = F)


##3.3. Effect direction and magnitude----
#TO-ADAPT:------
#Here we are calculating group-letters for each emmean. For the interaction plots and for tables, we will want to exctract the p-value for each comparions, create similar loop using "pairs(emmeans)"


#Emmeans comparisons give us the direction of effect and its magnitude for each level (averaging over the other effects). Also obtain emmeans for status within each season. 


#Notes on t-family models: upper and lower confidence limits are given as asymp.LCL/asymp.UCL (in our context we can treat them equally as the ones from gaussian models). Additionally, non-gaussian models will have Inf df (degrees of freedom are calculated differently for T-family models), not an issue. 


# Initialize list for storing comparisons
comparisons_list <- list()

# Loop to extract all comparisons for every casepilot best model
for (dataset in names(cores_best_model_list_allghg)) {
  cp_model <- cores_best_model_list_allghg[[dataset]]
  
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
all_emmeans_comparisons <- bind_rows(comparisons_list)


#Save all emmeans post-hocs:  
write.csv(x = all_emmeans_comparisons, 
          file = paste0(paper_path,"Allghg_emmeans-posthoc_coresmodels.csv"),
          row.names = F)

