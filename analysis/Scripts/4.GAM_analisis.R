# 3. Analyse GAM
# Explanation of this script ------------------------------------------------

# combination fo the habitat monitoring and the reef monitoring
# 
#
#


# Clean and load packages ------------------------------------------------------
cat("\014"); rm(list = ls())#; dev.off()
#sapply(.packages(), unloadNamespace)

#Set working derectory
setwd("~/PhD/Data/3. Natural habitat monitoring")

# Define the list of packages
packages <- c("gratia","mgcv","tidyverse", "dplyr", "ggplot2","readxl", "writexl","readr")

# Load packages if not already installed
lapply(packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, dependencies = TRUE)
  library(pkg, character.only = TRUE)})


# Import the data sets ---------------------------------------------------------
Monitoring_CPUE_data <- read_csv("Data_mod/Monitoring_CPUE_data.csv")

M_C_data <- Monitoring_CPUE_data 

names(M_C_data)
str(M_C_data)
summary(M_C_data)

# select variabeals most suitable ----------------------------------------------
GAM_data <- M_C_data %>%
  Weighted_CPUE_Kōura, Date_Time_Numeric , Mean_depth_m, Elevation_m 

# GAM (Generalized Additive Model) ---------------------------------------------
# Full model
model_full <- gam(Weighted_CPUE_Kōura ~ 
                    s(Date_Time_Numeric, k = 3) + 
                    s(Mean_depth_m, k = 3) + 
                    s(Elevation_m, k = 3) + 
                    s(Slope, k = 3) + 
                    s(Distance_deep, k = 3) + 
                    s(Temperature, k = 3) + 
                    s(TLI, k = 3) + 
                    s(Lake_water_clarity_m, k = 3) + 
                    s(Boulders, k = 3) + 
                    s(Mud, k = 3) + 
                    s(Overhanging_trees, k = 3) + 
                    Presence_Morihana +
                    Presence_Eel, 
                  data = M_C_data, family = gaussian)




# List of predictor variables
variables2 <- c(
  "s(Elevation_m, k = 3)", 
  "s(Date_Time_Numeric, k = 3)", 
  "s(Lake_water_clarity_m, k = 3)", 
  "s(Slope, k = 3)", 
  "s(Boulders, k = 3)", 
  "s(Overhanging_trees, k = 3)", 
  "s(Mean_depth_m, k = 3)", 
  "s(Mud, k = 3)", 
  "Riparian_vegetation", 
  "s(Temperature, k = 3)", 
  "s(TLI, k = 3)", 
  "s(Distance_deep, k = 3)", 
  "Wood_cover", 
  "Presence_Kōaro", 
  "Presence_rocks", 
  "Presence_Morihana", 
  "Presence_Eel")

# Response variable
response <-"Presence_Kōura"# "Weighted_CPUE_Kōura"

# Initialize a data frame to store results
results <- data.frame(
  Model = character(),
  Variables = character(),
  AIC = numeric(),
  stringsAsFactors = FALSE)

# Loop through all combinations of variables
for (i in 1:length(variables2)) {
  combos <- combn(variables2, i, simplify = FALSE) # Generate combinations of size 'i'
  for (combo in combos) {
    # Construct the formula
    formula <- as.formula(paste(response, "~", paste(combo, collapse = " + ")))
    
    # Fit the model
    model <- tryCatch(
      gam(formula, data = M_C_data, family = binomial), #binomial, gaussian poisson(link = "log")
      error = function(e) NULL)
    
    # If the model is successfully fitted, record the AIC
    if (!is.null(model)) {
      results <- rbind(
        results, 
        data.frame(
          Model = paste("Model", nrow(results) + 1),
          Variables = paste(combo, collapse = " + "),
          AIC = AIC(model)))}}}

# Sort results by AIC
results <- results[order(results$AIC), ]

# Display the top 10 models with the lowest AIC
head(results, 10)


best_model <- gam(
  Weighted_CPUE_Kōura ~ 
    s(Elevation_m, k = 3)+
    s(Date_Time_Numeric, k = 3) +
    s(Boulders, k = 3) + 
    s(Mud, k = 3) +
    Presence_Morihana,
  data = M_C_data,
  family = gaussian,
  method = "REML"
)
summary(best_model)
plot(best_model, pages = 1)
gam.check(best_model)

draw(best_model)





# Fit the GAM model with adjusted degrees of freedom (k)
model <- gam(Presence_Kōura ~ 
             #  s(Date_Time_Numeric, k = 3) + 
              # s(Mean_depth_m, k = 3) + 
            #   s(Elevation_m, k = 3) + 
               s(Slope, k = 3) + 
               s(Distance_deep, k = 3) + 
              # s(Temperature, k = 3) + 
           #    s(TLI, k = 3) + 
               s(Lake_water_clarity_m) + 
               s(Boulders, k = 3) + 
               s(Mud, k = 3) + 
               s(Overhanging_trees, k = 3) + 
               Presence_Morihana, 
              # Presence_Eel, 
             data = M_C_data, family = binomial)


summary(model)
plot(model, pages = 1)

# Scale numeric predictors
GAM_data <- M_C_data %>%
  mutate(across(c(Date_Time_Numeric, Mean_depth_m, Elevation_m, Slope, 
                  Distance_deep, Temperature, TLI, Lake_water_clarity_m, 
                  Boulders, Mud, Overhanging_trees), scale))

# Fit GAM model with adjusted settings
model2 <- gam(Weighted_CPUE_Kōura ~ 
                s(Date_Time_Numeric, k = 5) + 
                s(Mean_depth_m, k = 5) + 
                s(Elevation_m, k = 5) + 
                s(Slope, k = 5) + 
                s(Distance_deep, k = 5) + 
                s(Temperature, k = 5) + 
                s(TLI, k = 5) + 
                s(Lake_water_clarity_m, k = 5) + 
                s(Boulders, k = 5) + 
                s(Mud, k = 5) + 
                s(Overhanging_trees, k = 5) + 
                Presence_Morihana +
                Presence_Eel, 
              data = GAM_data, family = gaussian)

# Model diagnostics
summary(model2)    # View model summary
gam.check(model2)  # Check smoothness adequacy

# Plot results
plot(model2, pages = 1)








gam_model <- gam(Weighted_CPUE_Kōura ~ s(Distance_deep, k = 3) + s(Boulders, k = 3) + 
                   s(Overhanging_trees, k = 4),
                 family = poisson(link = "log"),
                 data = M_C_data)

gam_model <- gam(Weighted_CPUE_Kōura ~ Selected_variables, 
                 family = poisson(link = "log"), 
                 data = M_C_data)

plot(gam_model, pages = 1)

# Summarize the model
summary(gam_model)

plot(gam_model, pages = 1)












# GAM (Generalized Additive Model) -----------------------------------------------
library(mgcv)

# Fit the GAM model with adjusted degrees of freedom (k)
gam_model <- gam(Weighted_CPUE_Kōura ~ s(Weighted_CPUE_Kōura, k = 3) + s(Overhanging_trees, k = 3) + 
                   s(Weighted_CPUE_Morihana, k = 4),
                 family = poisson(link = "log"),
                 data = M_C_data)

gam_model <- gam(Weighted_CPUE_Kōura ~ Selected_variables, 
                 family = poisson(link = "log"), 
                 data = M_C_data)

plot(gam_model, pages = 1)

# Summarize the model
summary(gam_model)

plot(gam_model, pages = 1)

# GLM --------------------------------------------------------------------------
# List of predictor variables
variables_glm <- c(
  "Elevation_m", 
  "Date_Time_Numeric", 
  "Lake_water_clarity_m", 
  "Slope", 
  "Boulders", 
  "Overhanging_trees", 
  "Mean_depth_m", 
  "Mud", 
  "Riparian_vegetation", 
  "Temperature", 
  "TLI", 
  "Distance_deep", 
  "Wood_cover", 
  "Presence_Kōaro", 
  "Presence_rocks", 
  "Presence_Morihana", 
  "Presence_Eel"
)

# Response variable
response_glm <- "Presence_Kōura" # "Weighted_CPUE_Kōura"

# Initialize a data frame to store results
glm_results <- data.frame(
  Model = character(),
  Variables = character(),
  AIC = numeric(),
  BIC = numeric(),
  stringsAsFactors = FALSE
)

# Loop through all combinations of variables
for (i in 1:length(variables_glm)) {
  combos <- combn(variables_glm, i, simplify = FALSE) # Generate combinations of size 'i'
  for (combo in combos) {
    # Construct the formula
    formula <- as.formula(paste(response_glm, "~", paste(combo, collapse = " + ")))
    
    # Fit the model
    model <- tryCatch(
      glm(formula, data = M_C_data, family = binomial), # binomial, gaussian, poisson(link = "log")
      error = function(e) NULL
    )
    
    # If the model is successfully fitted, record the AIC and BIC
    if (!is.null(model)) {
      glm_results <- rbind(
        glm_results, 
        data.frame(
          Model = paste("Model", nrow(glm_results) + 1),
          Variables = paste(combo, collapse = " + "),
          AIC = AIC(model),
          BIC = BIC(model)
        )
      )
    }
  }
}

# Sort results by AIC (or BIC, depending on preference)
glm_results <- glm_results[order(glm_results$AIC), ]

# Display the top 10 models with the lowest AIC
head(glm_results, 10)


# PCA (Principal Component Analysis) -------------------------------------------
predictor_data <- M_C_data[, variables]
predictor_data_scaled <- scale(predictor_data)
pca_result <- prcomp(predictor_data_scaled, scale. = TRUE)
summary(pca_result)
biplot(pca_result)




# NMDS (Non-Metric Multidimensional Scaling) ---------------------------------------------------------------
library(vegan)

# Prepare data for NMDS: Exclude the response variable CPUE_Kōura if it’s not needed
predictor_data_clean <- predictor_data[complete.cases(predictor_data), ]

# Perform NMDS using Bray-Curtis dissimilarity (use other measures if suitable)
nmds_result <- metaMDS(predictor_data_clean, distance = "bray", k = 2, trymax = 100)

# Summarize and plot NMDS results
print(nmds_result)
plot(nmds_result, type = "t", main = "NMDS Plot")




# CCA (Canonical Correspondence Analysis) ---------------------------------------------------------------
# Separate response and predictor data (assuming CPUE data as response)
response_data <- Data1[, c("CPUE_Kōura", "CPUE_Bullies", "CPUE_Morihana", 
                           "CPUE_Eel", "CPUE_Common_smelt", "CPUE_Kōaro", "CPUE_Trout")]
environmental_data <- Data1[, c("Distance_deep", "Riparian_vegetation","Overhanging_trees",
                                "Boulders","Gravel","Sand","Mud","Organic_matter",
                                "Temperature","DO_mgl","DO_percent","Connectivity","pH",
                                "Emergent_weed","Submerged_weed", "Wood_cover")]

# Ensure no missing values in response or predictor data
response_data <- response_data[complete.cases(response_data), ]
environmental_data <- environmental_data[complete.cases(environmental_data), ]

# Perform CCA
cca_result <- cca(response_data ~ ., data = environmental_data)

# Summarize and plot CCA results
summary(cca_result)
plot(cca_result, main = "Canonical Correspondence Analysis (CCA)")




# RDA (Redundancy Analysis) ---------------------------------------------------------------
rda_result <- rda(response_data ~ ., data = environmental_data)

# Summarize and plot RDA results
summary(rda_result)
plot(rda_result, main = "Redundancy Analysis (RDA)")


# db-RDA (Constrained Ordination) ---------------------------------------------------------------
species_data <- Data1[, c("CPUE_Kōura", "CPUE_Bullies", "CPUE_Morihana", "CPUE_Eel", "CPUE_Common_smelt", "CPUE_Kōaro", "CPUE_Trout")]
environmental_data <- Data1[, c("Distance_deep", "Riparian_vegetation", "Overhanging_trees", 
                                "Bedrock", "Boulders", "Cobble", "Gravel", "Sand", "Mud", 
                                "Organic_matter", "Temperature", "DO_mgl", "DO_percent", 
                                "Connectivity", "pH", "Emergent_weed", "Submerged_weed", 
                                "Wood_cover")]

# Compute the dissimilarity matrix on the response variables (species data)
species_dist <- vegdist(species_data, method = "bray")

# Perform db-RDA using capscale (distance-based constrained ordination)
db_rda_model <- capscale(species_dist ~ ., data = environmental_data)

# Summarize and inspect the model
summary(db_rda_model)

# Visualize the db-RDA results
plot(db_rda_model, main = "db-RDA of Species Data Constrained by Environmental Predictors")





