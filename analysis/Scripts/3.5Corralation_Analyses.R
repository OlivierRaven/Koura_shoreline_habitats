# 3.5 Corralation Analyses 
# Explanation of this script ------------------------------------------------

# 
# 
#
#


# Clean and load packages ------------------------------------------------------
cat("\014"); rm(list = ls())#; dev.off()
#sapply(.packages(), unloadNamespace)

#Set working derectory
setwd("~/PhD/Data/3. Natural habitat monitoring")

# Define the list of packages
packages <- c("lmerTest","lme4","gratia","ggmosaic","randomForest","logistf","detectseparation","reshape2","ggridges","naniar","car","patchwork", "sf", "mgcv","tidyverse", "dplyr", "ggplot2","readxl", "writexl","readr")

# Load packages if not already installed
lapply(packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, dependencies = TRUE)
  library(pkg, character.only = TRUE)})


# Import the data sets ---------------------------------------------------------
Monitoring_CPUE_data <- read_csv("Data_mod/Monitoring_CPUE_data.csv")

M_C_data <- Monitoring_CPUE_data 

names(M_C_data)
summary(M_C_data)
str(M_C_data)
colSums(is.na(M_C_data))

# List of variables to test
Numeric_variables <- M_C_data %>% select_if(is.numeric) %>% names()

variables <- c("Presence_Kōura",#"Weighted_CPUE_Kōura",#"Total_Individuals_Kōura","", "Weighted_BCUE_Kōura",  #"Total_Weight_Kōura",
               "LID","Size_km2", "Catchment_size_km2","Mean_depth_m","Max_depth_m", "Elevation_m","Slope_5m","Distance_5m",
               "Date_Time_Numeric","Temperature", "DO_mgl", "DO_percent", "pH","Connectivity", "TLI", "Lake_water_clarity_m",
               "Bedrock","Boulders", "Cobble", "Gravel", "Sand", "Mud", "Organic_matter","Turf" ,"Presence_rocks", #"Substrate_index", 
               "Riparian_vegetation", "Overhanging_trees", "Emergent_Native","Emergent_Non_Native","Submerged_Native","Submerged_Non_Native", "Wood_cover",
               "Presence_Kōaro","Presence_Common_smelt", "Presence_Trout", "Presence_Morihana", "Presence_Eel","Presence_Catfish",
               "Dragonfly larvae Small","Diving beetle","Damselfly larvae","Dragonfly larvae Large","Snails")

# select variables based on common knowledge of variables that will have high corralition
Selected_variables <- c("Presence_Kōura","Mean_depth_m", "Elevation_m","Slope_5m","Temperature","DO_mgl", "pH","Bedrock", "Cobble", "Gravel", "Mud", "Organic_matter","Turf" ,"Presence_rocks", 
                        "Riparian_vegetation", "Overhanging_trees","Emergent_Non_Native","Submerged_Native","Submerged_Non_Native", "Wood_cover","Presence_Kōaro","Presence_Common_smelt", "Presence_Trout", "Presence_Morihana", "Presence_Eel","Presence_Catfish","Dragonfly larvae Small","Damselfly larvae","Snails")

variables1 <- c("Presence_Kōura", "Total_Individuals_Kōura","Weighted_CPUE_Kōura","Weighted_BCUE_Kōura","Date_Time_Numeric" ,"Mean_depth_m","Elevation_m","Slope_5m","Slope","Distance_5m", "Distance_20m","Calculated_slope", "Distance_deep","Temperature", "TLI", "Lake_water_clarity_m","Boulders", "Mud","Overhanging_trees" ,"Presence_Morihana" , "Presence_Eel")
variables2 <- c("Presence_Kōura", "Date_Time_Numeric","Mean_depth_m","Elevation_m","Slope_5m","Slope","Distance_5m", "Distance_20m","Calculated_slope", "Distance_deep","Temperature", "TLI", "Lake_water_clarity_m","Boulders", "Mud","Overhanging_trees" ,"Presence_Morihana" , "Presence_Eel")
variables3 <- c("Presence_Kōura", "Date_Time_Numeric","Mean_depth_m","Elevation_m","Slope_5m","Temperature", "TLI", "Lake_water_clarity_m","Boulders", "Mud","Overhanging_trees" ,"Presence_Morihana" , "Presence_Eel")
variables4 <- c("Presence_Kōura", "Date_Time_Numeric","Elevation_m","Slope_5m","Temperature","Boulders", "Mud","Overhanging_trees" ,"Presence_Morihana" , "Presence_Eel")
variables5 <- c("Presence_Kōura", "Elevation_m","Slope_5m","Temperature","Boulders","Overhanging_trees" ,"Presence_Morihana" , "Presence_Eel")

variables9 <- c("Presence_Kōura", "Distance_5m","Temperature","Boulders","Overhanging_trees" ,"Presence_Morihana" , "Presence_Eel")

M_C_data_numeric <- M_C_data[, colnames(M_C_data) %in% Numeric_variables]
M_C_data0 <- M_C_data[, colnames(M_C_data) %in% variables]
M_C_data1 <- M_C_data[, colnames(M_C_data) %in% variables1]
M_C_data2 <- M_C_data[, colnames(M_C_data) %in% variables2]
M_C_data3 <- M_C_data[, colnames(M_C_data) %in% variables3]
M_C_data4 <- M_C_data[, colnames(M_C_data) %in% variables4]
M_C_data5 <- M_C_data[, colnames(M_C_data) %in% variables5]
M_C_data6 <- M_C_data[, colnames(M_C_data) %in% variables6]

#
# Step 1: Check normality ------------------------------------------------------
normally_distributed_vars <- c()

# Loop through numeric variables in the dataset
for (var in variables) {
  if (var %in% colnames(M_C_data) && is.numeric(M_C_data[[var]])) {
    hist(M_C_data[[var]], main = paste("Histogram of", var), xlab = var, breaks = 20)
    qqnorm(M_C_data[[var]], main = paste("Q-Q plot of", var)); qqline(M_C_data[[var]], col = "red")
    shapiro_p <- shapiro.test(M_C_data[[var]])$p.value
    print(paste("Shapiro-Wilk test for", var, ": p-value =", shapiro_p))
    if (shapiro_p >= 0.05) normally_distributed_vars <- c(normally_distributed_vars, var)} 
  else {print(paste("Variable", var, "is not numeric or does not exist in the dataset"))}}


# Step 2: Perform Mann–Whitney U test for each variable ------------------------
results <- lapply(variables[-1], function(var) {
  wilcox.test(M_C_data[[var]] ~ M_C_data$Presence_Kōura   , data = M_C_data)})

# Format results into a summary
results_summary <- data.frame(
  Variable = variables[-1],
  p_value = sapply(results, function(x) x$p.value),
  W_statistic = sapply(results, function(x) x$statistic))

# Display significant results (e.g., p < 0.05)
significant_results <- results_summary %>%
  filter(p_value < 0.05)

print(significant_results)

# Step 3: Scatter plot ---------------------------------------------------------
# Create the pair plot
pairs(M_C_data5[, variables5])

pairs(M_C_data[, Selected_variables])

plot(M_C_data$Weighted_CPUE_Kōura, M_C_data$TLI)

# Step 4: Correlation matrix ---------------------------------------------------
cor_matrix <- cor(M_C_data[, Selected_variables], method = "spearman")  # "pearson", "spearman", "kendall" # spearman chosen as
cor_matrix_melted <- melt(cor_matrix)

# Create a heatmap
ggplot(cor_matrix_melted, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 2)), color = "black", size = 3) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Set threshold for strong correlations
threshold <- 0.5

# Filter for correlations above the threshold (absolute value)
strong_correlations <- cor_matrix_melted %>%
  filter(abs(value) > threshold, Var1 != Var2) %>%
  arrange(desc(abs(value)))


# Calculate VIF
model <- glm(Presence_Kōura ~ ., family = binomial, data = M_C_data5)
vif_values <- vif(model)
print(vif_values)

# explore imporatanc with randomForest)
rf_model <- randomForest(Presence_Kōura ~ ., data = M_C_data0, importance = TRUE)
varImpPlot(rf_model)


M_C_data5_scaled <- M_C_data5 %>%
  mutate(across(c("Calculated_slope", "Overhanging_trees", "Boulders", "Temperature"), scale))



# Create the design matrix and response variable
X <- model.matrix(Presence_Kōura ~ ., data = M_C_data5)
y <- M_C_data5$Presence_Kōura

# Perform the separation test
separation_test <- detect_separation(X, y, family = binomial())

# Display the results
print(separation_test)

library(logistf)
firth_model <- logistf(Presence_Kōura ~ ., data = M_C_data5)
summary(firth_model)

glm_model <- glm(Presence_Kōura ~ ., family = binomial, data = M_C_data5)
summary(glm_model)

step(glm_model, direction = "both")


# Presence_Kōura analysis -----------------------------------------------------
# Define relevant variables
variables <- c("Presence_Kōura", "Elevation_m","Temperature",
                "LID", "Mean_depth_m","Max_depth_m","Slope_5m",
                 "DO_mgl", "DO_percent", "pH", "TLI","Connectivity","Date_Time_Numeric",
                "Bedrock","Boulders", "Cobble", "Gravel", "Sand", "Mud", "Organic_matter","Turf" ,"Presence_rocks",  
                "Riparian_vegetation", "Overhanging_trees", "Emergent_Native","Emergent_Non_Native","Submerged_Native","Submerged_Non_Native", "Wood_cover",
                "Presence_Kōaro","Presence_Common_smelt", "Presence_Trout", "Presence_Morihana", "Presence_Eel","Presence_Catfish",
                "Dragonfly larvae Small","Diving beetle","Damselfly larvae","Dragonfly larvae Large","Snails")

variables <- c("Presence_Kōura",
                "LID", "Mean_depth_m","Max_depth_m", "Elevation_m","Slope_5m", "TLI","Temperature",
                "Date_Time_Numeric", "DO_mgl", "DO_percent", "pH","Connectivity", 
                "Bedrock","Boulders", "Cobble", "Gravel", "Sand", "Mud", "Organic_matter","Turf" ,"Presence_rocks",  
                "Riparian_vegetation", "Overhanging_trees", "Emergent_Native","Emergent_Non_Native","Submerged_Native","Submerged_Non_Native", "Wood_cover",
                "Presence_Kōaro","Presence_Common_smelt", "Presence_Trout", "Presence_Morihana", "Presence_Eel","Presence_Catfish",
                "Dragonfly larvae Small","Diving beetle","Damselfly larvae","Dragonfly larvae Large","Snails")


#variables removed round 1:,"Size_km2", "TLI","Elevation_m","Temperature","Catchment_size_km2", "Presence_Morihana", "Lake_water_clarity_m"
#variables removed round 2:,"Size_km2", "TLI", "Elevation_m","Temperature","Catchment_size_km2"
#variables removed round 3 based on VIF and heatmap:,"TLI","Connectivity","Lake_water_clarity_m","Date_Time_Numeric",

variables <- c("Presence_Kōura", "LID", "Mean_depth_m", "Slope_5m", "Temperature", 
               "DO_mgl", "DO_percent", "pH", "Bedrock", "Boulders", "Cobble", 
               "Gravel", "Sand", "Mud", "Organic_matter", "Turf", "Presence_rocks", 
               "Riparian_vegetation", "Overhanging_trees", "Emergent_Native", 
               "Emergent_Non_Native", "Submerged_Native", "Submerged_Non_Native", 
               "Wood_cover", "Presence_Kōaro", "Presence_Common_smelt", 
               "Presence_Trout", "Presence_Morihana", "Presence_Eel", 
               "Presence_Catfish", "Dragonfly larvae Small", "Diving beetle", 
               "Damselfly larvae", "Dragonfly larvae Large", "Snails")

# Subset data to relevant variables
M_C_data0 <- M_C_data[, colnames(M_C_data) %in% variables]

# Check for missing values
colSums(is.na(M_C_data0))

# 3. Initial Multicollinearity Check
cor_matrix_initial <- cor(M_C_data0[, sapply(M_C_data0, is.numeric)], use = "complete.obs", method = "pearson")
cor_melted_initial <- reshape2::melt(cor_matrix_initial)
ggplot(cor_melted_initial, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 2)), color = "black", size = 3) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Spearman Correlation Heatmap", x = "Variables", y = "Variables")

# Significant Variable Selection 
# Calculate correlations and filter significant variables
numeric_vars <- M_C_data0[, sapply(M_C_data0, is.numeric)]
#numeric_vars <- numeric_vars[, colnames(numeric_vars) != "Presence_Kōura"]
cor_results <- sapply(names(numeric_vars), function(var) {
  cor.test(numeric_vars[[var]], M_C_data0$Presence_Kōura, method = "pearson")$p.value})

adjusted_p <- p.adjust(cor_results, method = "fdr") # or method = "bonferroni"
significant_vars <- names(adjusted_p[adjusted_p < 0.05])
print(significant_vars)

# Visualize Significant Variables 
# Correlation heatmap
cor_matrix <- cor(M_C_data0[, significant_vars], use = "complete.obs", method = "pearson")
cor_melted <- melt(cor_matrix)
ggplot(cor_melted, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 2)), color = "black", size = 3) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Correlation Heatmap of Significant Variables")

# Preliminary GLM Model 
# Fit a GLM with significant variables
glm_formula <- as.formula(paste("Presence_Kōura ~", 
                                paste(significant_vars[!significant_vars %in% "Presence_Kōura"], collapse = " + ")))
glm_model <- glm(glm_formula, data = M_C_data0, family = binomial)

alias(glm_model)
# Calculate Variance Inflation Factor (VIF)
vif_values <- vif(glm_model)
print(vif_values)

# GAM Model Development 
# Fit initial GAM model
gam_model <- gam(Presence_Kōura ~ 
                   Elevation_m + TLI + Boulders + Temperature + s(Slope_5m, k = 10),
                    #Mean_depth_m + Elevation_m +Overhanging_trees + Boulders + Mud + s(Temperature, k = 10) + Presence_Morihana + Presence_Eel + s(Slope_5m, k = 10), 
                 data = M_C_data0, family = binomial)
summary(gam_model)

# Refine GAM model
gam_model_refined <- gam(Presence_Kōura ~ 
                         #  Presence_Morihana + Presence_Eel+ s(Slope_5m, k = 10),
                           Boulders  + s(Slope_5m, k = 10), 
                         data = M_C_data0, family = binomial)
summary(gam_model_refined)


# Validate and Visualize GAM Model 
# GAM diagnostics
gam.check(gam_model_refined)

# Draw smooth terms and residuals
draw(gam_model_refined, residuals = TRUE)

# Bar Plot with Annotations 
# Summarize data for plotting
data_with_counts <- M_C_data0 %>%
  group_by(Presence_Morihana, Presence_Eel, Presence_Kōura) %>%
  summarise(n = n(), .groups = "drop")

# Plot annotated bar plot
ggplot(data_with_counts, aes(x = factor(Presence_Morihana), y = n, fill = factor(Presence_Kōura))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = n), position = position_dodge(width = 0.9), vjust = -0.5, size = 3.5) +
  facet_wrap(~ Presence_Eel, labeller = label_both) +
  labs(x = "Presence of Morihana", y = "Count of Kōura Presence", fill = "Kōura Presence") +
  theme_minimal() +
  ggtitle("Effect of Morihana and Eel on Kōura Presence")

# Presence_Kōura analysis automated -----------------------------------------------------
# Define relevant variables
variables <- c("Presence_Kōura",
               "LID","Mean_depth_m","Slope_5m",
               "Date_Time_Numeric","Temperature", "DO_mgl", "pH","Connectivity", "Lake_water_clarity_m",
               "Bedrock","Boulders", "Cobble", "Gravel", "Sand", "Mud", "Organic_matter","Turf" ,"Presence_rocks", 
               "Riparian_vegetation", "Overhanging_trees", "Emergent_Native","Emergent_Non_Native","Submerged_Native","Submerged_Non_Native", "Wood_cover",
               "Presence_Kōaro","Presence_Common_smelt", "Presence_Trout", "Presence_Morihana", "Presence_Eel","Presence_Catfish",
               "Dragonfly larvae Small","Diving beetle","Damselfly larvae","Dragonfly larvae Large","Snails")

# Subset data to relevant variables
M_C_data0 <- M_C_data[, colnames(M_C_data) %in% variables]

# Check for missing values
print(colSums(is.na(M_C_data0)))

# Significant Variable Selection 
# Calculate correlations and filter significant variables
numeric_vars <- M_C_data0[, sapply(M_C_data0, is.numeric)]
cor_results <- sapply(names(numeric_vars), function(var) {
  cor.test(numeric_vars[[var]], M_C_data0$Presence_Kōura, method = "spearman")$p.value})

significant_vars <- names(cor_results[cor_results < 0.05])
print(significant_vars)

# Iterative VIF Reduction Process
reduce_vif <- function(data, response, threshold = 5) {
  # Separate response and predictors
  remaining_vars <- setdiff(names(data), response)
  
  repeat {
    # Construct the formula for the GLM
    glm_formula <- as.formula(paste(response, "~", paste(remaining_vars, collapse = " + ")))
    
    # Fit the GLM and handle aliasing issues
    glm_model <- tryCatch({
      glm(glm_formula, data = data, family = binomial)
    }, error = function(e) {
      stop("GLM fitting failed: ", e$message)
    })
    
    # Check for aliased coefficients
    aliased <- names(alias(glm_model)$Complete)
    if (length(aliased) > 0) {
      message("Removed aliased variable(s): ", paste(aliased, collapse = ", "))
      remaining_vars <- setdiff(remaining_vars, aliased)
      next  # Restart loop with updated variables
    }
    
    # Calculate VIF values and handle exceptions
    vif_values <- tryCatch({
      vif(glm_model)
    }, error = function(e) {
      stop("VIF calculation failed: ", e$message)
    })
    
    print(vif_values)  # Show the current VIF values
    
    # Check if all VIF values are below the threshold
    if (all(vif_values < threshold)) {
      message("All VIF values are below the threshold.")
      break
    }
    
    # Identify the variable with the highest VIF
    max_vif_var <- names(vif_values)[which.max(vif_values)]
    message("Removed variable due to high VIF: ", max_vif_var)
    remaining_vars <- setdiff(remaining_vars, max_vif_var)
    
    # Stop if no variables are left
    if (length(remaining_vars) == 0) {
      warning("No variables left after VIF reduction.")
      break
    }
  }
  
  # Return final results
  return(list(
    final_model = glm_model,
    remaining_vars = remaining_vars,
    vif_values = vif_values
  ))
}


# Run the iterative VIF reduction
vif_results <- reduce_vif(data = M_C_data0, response = "Presence_Kōura", threshold = 5)

# Final selected variables
print(vif_results$remaining_vars)

# Final GLM Model with reduced variables
final_glm_model <- vif_results$final_model
summary(final_glm_model)

# GAM Model Development 
# Example GAM with reduced variables
reduced_vars <- vif_results$remaining_vars
gam_formula <- as.formula(paste("Presence_Kōura ~", paste(reduced_vars, collapse = " + ")))
gam_model <- gam(gam_formula, data = M_C_data0, family = binomial)

# Validate and Visualize GAM Model
summary(gam_model)
gam.check(gam_model)
draw(gam_model, residuals = TRUE)




# Weighted_CPUE_Kōura analysis -------------------------------------------------
# 1. Exploratory Data Analysis (EDA)
#hist(M_C_data$Weighted_CPUE_Kōura), main = "Histogram of Weighted CPUE Kōura", xlab = "Weighted CPUE Kōura", col = "skyblue", border = "white")

# Additional EDA: Density plot and boxplot for Weighted_CPUE_Kōura
plot(density(M_C_data$Weighted_CPUE_Kōura, na.rm = TRUE), main = "Density Plot of Weighted CPUE Kōura", xlab = "Weighted CPUE Kōura")
boxplot(M_C_data$Weighted_CPUE_Kōura, main = "Boxplot of Weighted CPUE Kōura", ylab = "Weighted CPUE Kōura", col = "lightgreen")

# Summary statistics for variables
summary(M_C_data)

# 2. Define relevant variables
variables1 <- c("Weighted_CPUE_Kōura",
                "LID", "Mean_depth_m","Max_depth_m", "Elevation_m","Slope_5m", "TLI",
                "Date_Time_Numeric", "DO_mgl", "DO_percent", "pH","Connectivity", "Temperature",
                "Bedrock","Boulders", "Cobble", "Gravel", "Sand", "Mud", "Organic_matter","Turf" ,"Presence_rocks",  
                "Riparian_vegetation", "Overhanging_trees", "Emergent_Native","Emergent_Non_Native","Submerged_Native","Submerged_Non_Native", "Wood_cover",
                "Presence_Kōaro","Presence_Common_smelt", "Presence_Trout", "Presence_Morihana", "Presence_Eel","Presence_Catfish",
                "Dragonfly larvae Small","Diving beetle","Damselfly larvae","Dragonfly larvae Large","Snails")

#variables removed: "Size_km2","Lake_water_clarity_m","Catchment_size_km2","Temperature"

variables1 <- c("Weighted_CPUE_Kōura", "Boulders", "Overhanging_trees","Slope_5m", "Presence_Catfish","LID", "TLI")

# Subset data to relevant variables
M_C_data1 <- M_C_data[, colnames(M_C_data) %in% variables1]
summary(M_C_data1)
# Check for missing values
colSums(is.na(M_C_data1))

# 2.5. 
# Standardize numeric predictors to aid convergence and interpretation
num_vars <- c("Weighted_CPUE_Kōura", "Boulders", "Overhanging_trees","Slope_5m", "Presence_Catfish","LID", "TLI")

M_C_data1_scaled <- M_C_data1 %>%
  mutate(across(all_of(num_vars), scale))

# Fit the linear mixed effects model
# Random effect: LID (Lake ID) to account for grouping by lake
mixed_model <- lmer(log1p(Weighted_CPUE_Kōura) ~ 
                      Boulders+Overhanging_trees+Slope_5m+Presence_Catfish+
                      (1|LID),
                    data = M_C_data1_scaled)

# Summary of the model with p-values
summary(mixed_model)

vif_values <- vif(mixed_model)
print(vif_values)

# Optional: Check model diagnostics
par(mfrow = c(2, 2))
plot(mixed_model)

# 3. Initial Multicollinearity Check
cor_matrix_initial <- cor(M_C_data1[, sapply(M_C_data1, is.numeric)], use = "complete.obs", method = "spearman")
cor_melted_initial <- reshape2::melt(cor_matrix_initial)
ggplot(cor_melted_initial, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 2)), color = "black", size = 3) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Spearman Correlation Heatmap", x = "Variables", y = "Variables")

# 4. Significant Variable Selection 
numeric_vars <- M_C_data1[, sapply(M_C_data1, is.numeric)]
cor_results <- sapply(names(numeric_vars), function(var) {
  cor.test(numeric_vars[[var]], M_C_data1$Weighted_CPUE_Kōura, method = "spearman")$p.value})

adjusted_p <- p.adjust(cor_results, method = "fdr") # or method = "bonferroni"
significant_vars <- names(adjusted_p[adjusted_p < 0.05])
print(significant_vars)


# 4.5. PCA 
library(FactoMineR) 
library(factoextra)

### PCA for all variables
names(M_C_data1)

# reduce vriables
data_for_pca <- M_C_data1#[, !(names(M_C_data1) %in% c("LID", "Date_Time_Numeric"))]
PCA_data_sig <- M_C_data1[, significant_vars]

# Scale the data
data_scaled <- scale(data_for_pca)

# Perform PCA
pca_result <- prcomp(data_scaled, center = TRUE, scale. = TRUE)

# Summary of PCA
summary(pca_result)


# Make plots
fviz_eig(pca_result)                 # Variance explained
fviz_pca_biplot(pca_result)         # Biplot
fviz_pca_ind(pca_result)             # Samples on PCA plane
fviz_pca_var(pca_result)             # Variables contributions


# 5. GLM Model as Baseline
M_C_data1$Adjusted_CPUE <- M_C_data1$Weighted_CPUE_Kōura + 0.01
glm_formula <- as.formula(paste("Adjusted_CPUE ~", 
                                paste(significant_vars[!significant_vars %in% "Adjusted_CPUE"], collapse = " + ")))
glm_model <- glm(glm_formula, data = M_C_data1, family = Gamma(link = "log"))

summary(glm_model)

# 6. Multicollinearity Check (VIF)
vif_values <- vif(glm_model)
print(vif_values)



# 7. GAM Model Development
gam_model <- gam(Adjusted_CPUE ~ 
                   Elevation_m + Overhanging_trees + Boulders + TLI +
                   s(Connectivity, k = 10) +
                   s(Slope_5m, k = 3), 
                 data = M_C_data1, 
                 family = Gamma(link = "log"))
summary(gam_model)
gam.check(gam_model)

# 8. Refined GAM Model
gam_model_refined <- gam(Adjusted_CPUE ~ 
                           Overhanging_trees + s(Slope_5m, k = 10, select = TRUE), 
                         data = M_C_data1, family = Gamma(link = "log"))
summary(gam_model_refined)

# AIC Comparison
AIC(gam_model, gam_model_refined)

# GAM Validation
gam.check(gam_model_refined)

# Visualizations
plot(M_C_data1$Adjusted_CPUE ~ M_C_data1$Overhanging_trees)
plot(gam_model_refined, pages = 1, residuals = TRUE)
draw(gam_model_refined, residuals = TRUE)
