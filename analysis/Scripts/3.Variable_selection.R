# 3. Variable selection
# Explanation of this script ------------------------------------------------

# 
# 
#
#


# Clean and load packages ------------------------------------------------------
cat("\014"); rm(list = ls())
#dev.off();sapply(.packages(), unloadNamespace)

#Set working derectory
setwd("~/PhD/Data/3. Natural habitat monitoring")

# Define the list of packages
packages <- c("gratia", "mgcv", "corrplot","glmmTMB", "performance", "dplyr", "tidyverse", "writexl")

# Load packages if not already installed
lapply(packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, dependencies = TRUE)
  library(pkg, character.only = TRUE)})


# Import the data set ----------------------------------------------------------
Monitoring_CPUE_data <- read_csv("Data_mod/Monitoring_CPUE_data.csv")

# Filter to ignor all sites that have been sampled double.
M_C_data <- Monitoring_CPUE_data %>%
  filter(Monitoring==0)

names(M_C_data)
#summary(M_C_data)
#str(M_C_data)
#colSums(is.na(M_C_data))

# Built Functions --------------------------------------------------------------
# remove high vif in glmmTMB function
remove_high_vif_glmmTMB <- function(data, response, predictors, threshold = 5) {
  removed <- character()
  repeat {
    formula <- as.formula(paste(response, "~", paste(predictors, collapse = " + "))) # , "+ (1 | LID)"
    model <- try(glmmTMB(formula, data = data, family = tweedie(link = "log")), silent = TRUE)
    if (inherits(model, "try-error")) break
    vif_data <- performance::check_collinearity(model)
    vif_data <- vif_data[!grepl("\\|", vif_data$Term), ]
    if (all(vif_data$VIF < threshold)) break
    to_remove <- vif_data$Term[which.max(vif_data$VIF)]
    cat("Removing:", to_remove, "\n")
    predictors <- setdiff(predictors, to_remove)
    removed <- c(removed, to_remove)
    if (length(predictors) == 0) stop("All predictors removed.")}
  cat("Removed variables:\n", paste(removed, collapse = ", "), "\n")
  return(model)
}

# remove high vif in glm function for binominal data
remove_high_vif_glmBI <- function(data, response, predictors, threshold = 5) {
  removed <- character()
  repeat {
    formula <- as.formula(paste(response, "~", paste(predictors, collapse = " + "))) # , "+ (1 | LID)"
    model <- try(glm(formula, data = data, family = binomial(link = "logit")), silent = TRUE)
    if (inherits(model, "try-error")) break
    vif_data <- performance::check_collinearity(model)
    vif_data <- vif_data[!grepl("\\|", vif_data$Term), ]
    if (all(vif_data$VIF < threshold)) break
    to_remove <- vif_data$Term[which.max(vif_data$VIF)]
    cat("Removing:", to_remove, "\n")
    predictors <- setdiff(predictors, to_remove)
    removed <- c(removed, to_remove)
    if (length(predictors) == 0) stop("All predictors removed.")}
  cat("Removed variables:\n", paste(removed, collapse = ", "), "\n")
  return(model)
}

# Build smooitng funcion
build_smooth <- function(var) {
  n <- length(unique(M_C_data_scaled3[[var]]))
  if (n >= 5) {
    k <- if (!is.null(custom_k[[var]])) paste0(", k=", custom_k[[var]]) else ""
    paste0("s(", var, ", bs='ts'", k, ")")
  } else var}

# automatic remove the variable with highest p value untill only p < 0.1
stepwise_gam_vars <- function(response, vars, data, family, method = "ML", p_cutoff = 0.1) {
  remaining_vars <- vars
  m <- NULL  # Initialize model object outside loop
  
  repeat {
    # Build formula with smooth terms
    rhs_terms <- sapply(remaining_vars, build_smooth)
    formula_str <- paste(response, "~", paste(rhs_terms, collapse = " + "))
    current_formula <- as.formula(formula_str)
    
    # Fit GAM model
    m <- gam(current_formula, data = data, family = family, method = method)
    summ <- summary(m)
    
    # Extract parametric terms p-values, excluding intercept
    if (!is.null(summ$p.table) && nrow(summ$p.table) > 0) {
      pcol <- intersect(c("Pr(>|t|)", "Pr(>|z|)"), colnames(summ$p.table))
      if (length(pcol) == 0) stop("No recognised p-value column in parametric table")
      param_p <- summ$p.table[, pcol[1]]
      names(param_p) <- rownames(summ$p.table)
      param_p <- param_p[names(param_p) != "(Intercept)"]
    } else {
      param_p <- numeric(0)
    }
    
    # Extract smooth terms p-values
    smooth_p <- summ$s.table[, "p-value"]
    names(smooth_p) <- rownames(summ$s.table)
    
    # Combine all p-values
    all_p <- c(param_p, smooth_p)
    
    # Filter p-values above cutoff
    all_p_filtered <- all_p[all_p > p_cutoff]
    
    # If no terms have p-value above cutoff, stop
    if (length(all_p_filtered) == 0) {
      message("All terms have p <= ", p_cutoff, ". Stopping.")
      break
    }
    
    # Find term with max p-value above cutoff
    term_to_remove <- names(which.max(all_p_filtered))
    max_pval <- max(all_p_filtered)
    
    message("Removing term: ", term_to_remove, " with p-value = ", max_pval)
    
    # Identify the variable to remove:
    if (grepl("^s\\(", term_to_remove)) {
      var_to_remove <- sub("^s\\(([^,]+).*\\)$", "\\1", term_to_remove)
    } else {
      var_to_remove <- term_to_remove
    }
    
    # Remove variable from remaining_vars
    if (var_to_remove %in% remaining_vars) {
      remaining_vars <- setdiff(remaining_vars, var_to_remove)
    } else {
      message("Warning: variable '", var_to_remove, "' not found in remaining_vars!")
      break  # avoid infinite loop
    }
    
    if (length(remaining_vars) == 0) {
      message("No variables left, stopping.")
      break
    }
  }
  
  return(m)
}

#
# Test weighted cpue kōura (abundances)-----------------------------------------
# 1. Preprocessing 
vars <- c("Weighted_CPUE_Kōura", #"Lake",
          "Slope_5m", "Riparian_vegetation", "Overhanging_trees",
          "Wood_cover",# "Substrate_index", 
          "Bedrock", "Boulders", "Cobble", "Gravel", "Sand", "Organic_matter",
          "Temperature", "DO_percent", "pH", "Specific_conductivity", "DO_mgl",
          "Emergent_Native", "Submerged_Native", "Submerged_Non_Native", "Emergent_Non_Native", #"Turf_Native",
          "Presence_Kōaro", "Presence_Common_smelt", "Presence_Morihana", 
          "Presence_Eel", "Presence_Catfish"
         )

M_C_data_scaled <- M_C_data %>%
  dplyr::select(all_of(vars)) %>%
  na.omit() %>%
  mutate(across(where(is.numeric) & !c(Weighted_CPUE_Kōura))) #, scale

# 2. Correlation check 
cor_matrix <- cor(dplyr::select(M_C_data_scaled, where(is.numeric)), use = "complete.obs")
corrplot(cor_matrix)

# 3. VIF reduction 
response <- "Weighted_CPUE_Kōura"
predictors <- setdiff(names(M_C_data_scaled), response)

model <- remove_high_vif_glmmTMB(M_C_data_scaled, response, predictors)
summary(model)
performance::check_collinearity(model)

# 4. Post-VIF correlation filtering 
clean_vars <- c(setdiff(names(fixef(model)$cond), "(Intercept)"), response)
cor_filtered <- M_C_data_scaled %>% dplyr::select(all_of(clean_vars))
corrplot(cor(cor_filtered, use = "complete.obs"))
clean_vars <- setdiff(clean_vars, response)



pca <- prcomp(M_C_data_scaled[, clean_vars], scale. = TRUE)

library(factoextra)
fviz_eig(pca)
fviz_pca_biplot(pca)
fviz_pca_ind(pca)
fviz_pca_var(pca)


# spearman_test
spearman_test <- sapply(clean_vars, function(var) {cor.test(M_C_data_scaled[[var]], M_C_data_scaled$Weighted_CPUE_Kōura, method = "spearman", exact = FALSE)$p.value})
significant <- names(spearman_test)[spearman_test < 0.05]
significant

# 5. Build GAM formulas 
# Custom smooth parameters
sapply(M_C_data_scaled[clean_vars], function(x) length(unique(x)))
custom_k <- list(Overhanging_trees=5, Riparian_vegetation=7,
                 Bedrock=3, Boulders=5, Cobble=4, Gravel=6,
                 Emergent_Native=5, Submerged_Native=4, Submerged_Non_Native=9, Turf_Native=6)


m_final <- stepwise_gam_vars(
  response = response,
  vars = clean_vars,
  data = M_C_data_scaled,
  family = tw(link = "log"),
  method = "ML",
  p_cutoff = 0.1)

summary(m_final)
formula_final <- formula(m_final)

formula_final <- Weighted_CPUE_Kōura ~ Lake+ s(Boulders, bs = "ts", k = 5) + s(Cobble,bs = "ts", k = 4) + 
  #s(Temperature, bs = "ts") + s(pH, bs = "ts") + 
  s(Emergent_Native, bs = "ts", k = 5) + 
  s(Submerged_Non_Native, bs = "ts", k = 9) + Presence_Common_smelt


m1 <- gam(formula_final, family = tw(link = "log"), method = "REML", data = M_C_data_scaled)
summary(m1)
gam.check(m1)
concurvity(m1)
influence.gam(m1)

# Visual diagnostics
plot(m1, pages = 1, shade = TRUE)
draw(m1)

# Build new data frame 
new_data <- expand.grid(
  Presence_Common_smelt = c(0, 1),
  Temperature = mean(M_C_data_scaled$Temperature, na.rm = TRUE),
  pH = mean(M_C_data_scaled$pH, na.rm = TRUE),
  Emergent_Native = mean(M_C_data_scaled$Emergent_Native, na.rm = TRUE),
  Substrate_index = mean(M_C_data_scaled$Substrate_index, na.rm = TRUE),
  Boulders = mean(M_C_data_scaled$Boulders, na.rm = TRUE),
  Submerged_Non_Native = mean(M_C_data_scaled$Submerged_Non_Native, na.rm = TRUE),
  Cobble = mean(M_C_data_scaled$Cobble, na.rm = TRUE))

# Predict using your model
new_data$predicted_CPUE <- predict(m1, newdata = new_data, type = "response")

p1 <- draw(m1)
p2<-ggplot(new_data, aes(factor(Presence_Common_smelt), predicted_CPUE)) +
  geom_col(width = 0.5) +
  labs(x = "Presence Common Smelt", y = "Predicted CPUE") 

library(patchwork)
p1 + p2 + plot_layout(ncol = 2)


M_C_data %>%
  #filter(pH >= 6) %>%
  ggplot(aes(pH, Weighted_CPUE_Kōura)) +
  geom_point(aes(col = Lake)) +
  geom_smooth()


boxplot(M_C_data_scaled$Weighted_CPUE_Kōura~M_C_data_scaled$LID)

ggplot(M_C_data, aes(y=Weighted_CPUE_Kōura,x= Slope_5m))+
  geom_point(aes(col=Lake))

# Testing lake effect
library(lme4)
model_mixed <- lmer(
  Weighted_CPUE_Kōura ~ Boulders + Cobble + Temperature + Emergent_Native + 
    Presence_Common_smelt + (1 | LID),
  data = M_C_data_scaled,
  REML = FALSE)

model_mixed_null <- lm(
  Weighted_CPUE_Kōura ~ Boulders + Cobble + Temperature + Emergent_Native + 
    Presence_Common_smelt,
  data = M_C_data_scaled,
  REML = FALSE)

lapply(list(model_mixed_null, model_mixed), summary)
AIC(model_mixed_null, model_mixed)




#------
formula_0 <- as.formula(paste(response, "~", 
                              paste(sapply(setdiff(remaining_vars, response), build_smooth), collapse = " + ")))

formula_0.1 <- Weighted_CPUE_Kōura ~
  s(Slope_5m, bs = "ts") + 
  s(Riparian_vegetation, bs = "ts", k = 7) + 
  s(Overhanging_trees, bs = "ts", k = 5) + 
  s(Wood_cover, bs = "ts") + 
  s(Bedrock, bs = "ts", k = 3) +
  s(Boulders, bs = "ts", k = 5) + 
  s(Cobble, bs = "ts", k = 4) + 
  s(Gravel, bs = "ts", k = 6) + 
  s(Sand, bs = "ts") + 
  s(Temperature, bs= "ts") +
  s(pH, bs = "ts") + 
  s(Specific_conductivity, bs = "ts") +
  s(DO_mgl, bs = "ts") + 
  s(Emergent_Native, bs = "ts", k = 5) +
  s(Submerged_Native, bs = "ts", k = 4) + 
  s(Submerged_Non_Native, bs = "ts", k = 9) + 
  Emergent_Non_Native + Presence_Kōaro + Presence_Common_smelt + Presence_Catfish +  Presence_Morihana + Presence_Eel
m0.1 <- gam(formula_0.1, family = tw(link = "log"), method = "ML", data = M_C_data_scaled)
summary(m0.1)


# Manually curated formulas (based on p-values)
formula_1 <- Weighted_CPUE_Kōura ~
  s(Temperature, bs = "ts") +
  s(pH, bs = "ts") +
  s(Emergent_Native, bs = "ts", k = 5) +
  Presence_Common_smelt +
  s(Boulders, bs = "ts", k = 5) +
  s(Cobble, bs = "ts", k = 4) 


# backwards selection remove variable with the highest p value run the model and do again until all are below 0.1  
formula_8 <- Weighted_CPUE_Kōura ~
  #s(Slope_5m, bs = "ts") + 
  #s(Riparian_vegetation, bs = "ts", k = 7) + 
  #s(Overhanging_trees, bs = "ts", k = 5) + 
  #s(Wood_cover, bs = "ts") + 
  #s(Bedrock, bs = "ts", k = 3) +
  s(Boulders, bs = "ts", k = 5) + 
  s(Cobble, bs = "ts", k = 4) + 
  #s(Gravel, bs = "ts", k = 6) + 
  #s(Sand, bs = "ts") + 
  s(Temperature, bs= "ts") +
  s(pH, bs = "ts") + 
  #s(Specific_conductivity, bs = "ts") +
  #s(DO_mgl, bs = "ts") + 
  s(Emergent_Native, bs = "ts", k = 5) +
  #s(Submerged_Native, bs = "ts", k = 4) + 
  s(Submerged_Non_Native, bs = "ts", k = 9) + 
  #Emergent_Non_Native + #Presence_Kōaro + 
  Presence_Common_smelt #+ #Presence_Catfish +  Presence_Morihana #+ Presence_Eel
m8 <- gam(formula_8, family = tw(link = "log"), method = "ML", data = M_C_data_scaled)
summary(m8)


formula_9 <- Weighted_CPUE_Kōura ~
  s(LID, bs = "ts", k = 4)+
  #s(Slope_5m, bs = "ts") + 
  #s(Riparian_vegetation, bs = "ts", k = 7) + 
  #s(Overhanging_trees, bs = "ts", k = 5) + 
  #s(Wood_cover, bs = "ts") + 
  #s(Bedrock, bs = "ts", k = 3) +
  s(Boulders, bs = "ts", k = 5) + 
  s(Cobble, bs = "ts", k = 4) + 
  #s(Gravel, bs = "ts", k = 6) + 
  #s(Sand, bs = "ts") + 
  s(Temperature, bs= "ts") +
  s(pH, bs = "ts") + 
  #s(Specific_conductivity, bs = "ts") +
  #s(DO_mgl, bs = "ts") + 
  s(Emergent_Native, bs = "ts", k = 5) +
  #s(Submerged_Native, bs = "ts", k = 4) + 
  s(Submerged_Non_Native, bs = "ts", k = 9) + 
  #Emergent_Non_Native + #Presence_Kōaro + 
  Presence_Common_smelt #+ #Presence_Catfish +  Presence_Morihana #+ Presence_Eel
m9 <- gam(formula_9, family = tw(link = "log"), method = "ML", data = M_C_data_scaled)
summary(m9)


formula_2 <- update(formula_1, . ~ . - s(Cobble, bs = "ts", k = 4))
formula_3 <- update(formula_1, . ~ . - s(Boulders, bs = "ts", k = 5))
formula_4 <- update(formula_1, . ~ . - Presence_Common_smelt)
formula_5 <- update(formula_1, . ~ . - s(Emergent_Native, bs = "ts", k = 5))
formula_6 <- update(formula_1, . ~ . - s(pH, bs = "ts"))
formula_7 <- update(formula_1, . ~ . - s(Temperature, bs = "ts"))


# 6. Fit GAMs 
m0 <- gam(formula_0, family = tw(link = "log"), method = "ML", data = M_C_data_scaled)
m1 <- gam(formula_1, family = tw(link = "log"), method = "ML", data = M_C_data_scaled)
m2 <- gam(formula_2, family = tw(link = "log"), method = "ML", data = M_C_data_scaled)
m3 <- gam(formula_3, family = tw(link = "log"), method = "ML", data = M_C_data_scaled)
m4 <- gam(formula_4, family = tw(link = "log"), method = "ML", data = M_C_data_scaled)
m5 <- gam(formula_5, family = tw(link = "log"), method = "ML", data = M_C_data_scaled)
m6 <- gam(formula_6, family = tw(link = "log"), method = "ML", data = M_C_data_scaled)
m7 <- gam(formula_7, family = tw(link = "log"), method = "ML", data = M_C_data_scaled)
m8 <- gam(formula_8, family = tw(link = "log"), method = "ML", data = M_C_data_scaled)


# 7. Model summaries & comparison 
lapply(list(m0, m1, m2, m3, m4, m5, m6, m7, m8), summary)
anova(m1, m0, test = "Chisq")
anova(m2, m1, test = "Chisq")
anova(m3, m1, test = "Chisq")
anova(m4, m1, test = "Chisq")
anova(m5, m1, test = "Chisq")
anova(m6, m1, test = "Chisq")
anova(m7, m1, test = "Chisq")
anova(m1, m8, test = "Chisq")
anova(m8, m9, test = "Chisq")
AIC(m0, m1, m2, m3, m4, m5, m6, m7, m8)

# 8. Final model diagnostics & plotting 
m8.1 <- gam(formula_8, family = tw(link = "log"), method = "REML", data = M_C_data_scaled)
summary(m8.1)
gam.check(m8.1)
concurvity(m8.1)
influence.gam(m8.1)

# Visual diagnostics
plot(m8.1, pages = 1, shade = TRUE)
draw(m8.1)

# Predicted vs observed
plot(predict(m8.1, type = "response"), m8.1$y,
     xlab = "Predicted CPUE", ylab = "Observed CPUE")
abline(0, 1, col = "red")

# Cook’s Distance-style
plot(influence.gam(m8.1), type = "h")
abline(h = 0.2, col = "red", lty = 2)

# plot Presence_Common_smelt
plot(M_C_data$Weighted_CPUE_Kōura, M_C_data$Presence_Common_smelt)
plot(M_C_data$Site_ID, M_C_data$Weighted_CPUE_Common_smelt, col=M_C_data$Weighted_CPUE_Kōura)


# Build new data frame 
new_data <- expand.grid(
  Presence_Common_smelt = c(0, 1),
  Temperature = mean(M_C_data_scaled$Temperature, na.rm = TRUE),
  pH = mean(M_C_data_scaled$pH, na.rm = TRUE),
  Emergent_Native = mean(M_C_data_scaled$Emergent_Native, na.rm = TRUE),
  Boulders = mean(M_C_data_scaled$Boulders, na.rm = TRUE),
  Submerged_Non_Native = mean(M_C_data_scaled$Submerged_Non_Native, na.rm = TRUE),
  Cobble = mean(M_C_data_scaled$Cobble, na.rm = TRUE))

# Predict using your model
new_data$predicted_CPUE <- predict(m8.1, newdata = new_data, type = "response")

p1 <- draw(m8.1)
p2<-ggplot(new_data, aes(factor(Presence_Common_smelt), predicted_CPUE)) +
  geom_col(width = 0.5) +
  labs(x = "Presence Common Smelt", y = "Predicted CPUE") 

library(patchwork)
p1 + p2 + plot_layout(ncol = 4)























library(rpart)
library(rpart.plot)
# Formula with all predictors from m0
formula <- as.formula(paste(response, "~", paste(setdiff(remaining_vars, response), collapse = " + ")))
tree_model <- rpart(formula,data = M_C_data_scaled,method = "anova")

# Plot the tree
rpart.plot(tree_model)
summary(tree_model)

#
# OLD stuff ####
# Select the best model function
get_final_model <- function(model, data) {
  options(na.action = "na.fail")
  dredged <- dredge(model, trace = FALSE)
  avg <- model.avg(dredged, subset = delta < 4)
  vars <- names(sw(avg)[sw(avg) > 0.5])  
  vars_clean <- gsub("^cond\\(|\\)$", "", vars)
  form <- as.formula(paste("Weighted_CPUE_Kōura ~", paste(vars_clean, collapse = " + ")))
  refit <- glmmTMB(form, data = data, family = tweedie(link = "log"), ziformula = ~0, dispformula = ~1)
  return(list(model = refit, averaged = avg, used_models = subset(dredged, delta < 4)))
}


# 1. Preprocess
M_C_data_scaled <- M_C_data %>%
  dplyr::select(all_of(c(
    "Weighted_CPUE_Kōura", #"LID", #"Mean_depth_m", #"TLI", 
    "Slope_5m", "Riparian_vegetation", "Overhanging_trees",
    "Wood_cover", "Bedrock", "Boulders", "Cobble", "Gravel", "Sand", "Organic_matter", #"Mud",#Mud: Estimate = -35.67, SE = 22,180#"Presence_rocks", # binominal from boulders and cobble
    "Temperature", "DO_percent", "pH","Specific_conductivity", "DO_mgl", 
    "Emergent_Native", "Submerged_Native", "Submerged_Non_Native", "Emergent_Non_Native", 
    "Presence_Kōaro","Presence_Common_smelt","Presence_Morihana","Presence_Eel", "Presence_Catfish"
  ))) %>%
  na.omit() %>%
  mutate(across(where(is.numeric) & !c(Weighted_CPUE_Kōura), scale))#,LID = factor(LID))

# Explore correlations
cor_data <- M_C_data_scaled %>% dplyr::select(where(is.numeric))#, -LID
cor_matrix <- cor(cor_data, use = "complete.obs")
corrplot(cor_matrix)

# 2. VIF-based reduction
response <- "Weighted_CPUE_Kōura"
predictors <- setdiff(names(M_C_data_scaled), c("Weighted_CPUE_Kōura")) #, "LID"
model <- remove_high_vif_glmmTMB(M_C_data_scaled, response, predictors)
performance::check_collinearity(model)

# Correlation filtering after VIF
remaining_vars <- c(setdiff(names(fixef(model)$cond), "(Intercept)"),"Weighted_CPUE_Kōura")
cor_vars <- M_C_data_scaled %>% 
  dplyr::select(all_of(remaining_vars)) 
cor_matrix <- cor(cor_vars, use = "complete.obs")
corrplot(cor_matrix)
cor_df <- as.data.frame(as.table(cor_matrix)) %>% 
  filter(Var1 != Var2) %>% arrange(desc(abs(Freq)))

# Gam modeling
sapply(remaining_vars, function(var) length(unique(M_C_data_scaled[[var]])))

library(mgcv)

# Custom k values for specific smooths
custom_k <- list(Mean_depth_m=4, Boulders=5, Cobble=4, Emergent_Native=5, Overhanging_trees=5,
                 Riparian_vegetation=7, Bedrock=3, Gravel=6, Submerged_Native=4,
                 Submerged_Non_Native=9)

# Auto formula
formula_1 <- as.formula(paste(
  "Weighted_CPUE_Kōura ~",
  paste(sapply(setdiff(remaining_vars, "Weighted_CPUE_Kōura"), function(v) {
    n <- length(unique(M_C_data_scaled[[v]]))
    if (n >= 5) {
      k <- if (!is.null(custom_k[[v]])) paste0(", k=", custom_k[[v]]) else ""
      paste0("s(", v, ", bs='ts'", k, ")")
    } else v
  }), collapse = " + ")))

# only keep variables with p < 0.1
formula_2 <- (Weighted_CPUE_Kōura ~
    s(Temperature, bs = "ts")+
    s(Boulders, bs = "ts", k = 5) +
    s(Emergent_Native, bs = "ts", k = 5) +
    s(pH, bs = "ts") +
    s(Cobble, bs = "ts", k = 4)+
    Presence_Common_smelt)

# remove coubles 
formula_3 <- (Weighted_CPUE_Kōura ~ 
  s(Temperature, bs = "ts") + 
  s(Boulders, bs = "ts", k = 5) + 
  s(Emergent_Native, bs = "ts", k = 5) + 
  s(pH, bs = "ts") + 
  #s(Cobble, bs = "ts", k = 4) +  
  Presence_Common_smelt)

# remove 
formula_4 <- (Weighted_CPUE_Kōura ~ 
                #s(Temperature, bs = "ts") + 
                s(Boulders,bs = "ts", k = 5) + 
                s(Emergent_Native, bs = "ts", k = 5) + 
                s(pH, bs = "ts") + 
                s(Cobble, bs = "ts", k = 4) +  
                Presence_Common_smelt)


m1 <- gam(formula_1, family = tw(link = "log"), method = "ML", data = M_C_data_scaled)
m2 <- gam(formula_2, family = tw(link = "log"), method = "ML", data = M_C_data_scaled)
m3 <- gam(formula_3, family = tw(link = "log"), method = "ML", data = M_C_data_scaled)
m4 <- gam(formula_4, family = tw(link = "log"), method = "ML", data = M_C_data_scaled)
#m5 <- gam(formula_5, family = tw(link = "log"), method = "ML", data = M_C_data_scaled)
#m6 <- gam(formula_6, family = tw(link = "log"), method = "ML", data = M_C_data_scaled)

summary(m1)
summary(m2)
summary(m3)
summary(m4)

# Compare models using ANOVA
anova(m1, m2, m3, m4, test = "Chisq")
AIC(m1, m2, m3, m4)

m2.1 <- gam(formula_2, family = tw(link = "log"), data = M_C_data_scaled)
summary(m2.1)
gam.check(m2.1)
concurvity(m2.1)
influence.gam(m2.1)
plot(m2.1, pages = 1, shade = TRUE)

library(gratia)
draw(m2.1)
plot(predict(m2.1, type = "response"), m2.1$y,
     xlab = "Predicted CPUE", ylab = "Observed CPUE")
abline(0, 1, col = "red")




# 3. select variables based on RF and adjusted P values
# Adjusted p-values
coefs_df <- as.data.frame(summary(model)$coefficients$cond)
coefs_df$adj_p <- p.adjust(coefs_df$`Pr(>|z|)`, method = "fdr")
summary(model)
sig_vars <- rownames(coefs_df)[coefs_df$`Pr(>|z|)` < 0.05 & rownames(coefs_df) != "(Intercept)"]
sig_vars_adj <- rownames(coefs_df)[coefs_df$adj_p < 0.5 & rownames(coefs_df) != "(Intercept)"]

# RF
library(Boruta)
set.seed(123)
boruta_result <- Boruta(Weighted_CPUE_Kōura ~ ., data = M_C_data_scaled, doTrace = 1)
confirmed_vars <- getSelectedAttributes(boruta_result, withTentative = T)


# 4. Reduced model
relavant_var <- sort(unique(c(confirmed_vars, sig_vars)))
final_formula <- as.formula(paste(response, "~", paste(relavant_var, collapse = " + "))) #, "+ (1|LID)"
reduced_model <- glmmTMB(final_formula, data = M_C_data_scaled, family = tweedie(link = "log"))
summary(reduced_model)
performance::check_collinearity(reduced_model)

# Perform PCA
pca <- prcomp(M_C_data_scaled[, relavant_var], scale. = TRUE)

library(factoextra)
fviz_eig(pca)
fviz_pca_biplot(pca)
fviz_pca_ind(pca)
fviz_pca_var(pca)

# 5. Model selection
result <- get_final_model(reduced_model, M_C_data_scaled)
summary(result$model)           # Final refit model
summary(result$averaged)        # Model-averaged output
result$used_models              # Models included in averaging
performance::check_collinearity(result$model)




# 6. Model visualization
# Simulate residuals with DHARMa and plot diagnostics
sim_res <- simulateResiduals(fittedModel = result$model, plot = TRUE)
check_model(result$model)


# check best model vs basic model
library(cplm)
null_model <- cpglm(Weighted_CPUE_Kōura ~ 1, 
                    data = M_C_data_scaled)
AIC(null_model, result$model)


# Plot predicted vs observed
predicted <- predict(result$model, type = "response")  
observed <- M_C_data_scaled$Weighted_CPUE_Kōura     

plot(predicted, observed)
abline(a = 0, b = 1, col = "red", lwd = 2)

plot(M_C_data_scaled$Weighted_CPUE_Kōura, M_C_data_scaled$pH)


# p- adjust AIC 209.1
# RF model AIC 208.1

# Weighted_CPUE_Kōura ~ Boulders + Temperature, AIC: 208.1
# Weighted_CPUE_Kōura ~ Boulders + pH + Temperature, AIC: 206.3
# Weighted_CPUE_Kōura ~ Boulders + Emergent_Native + Temperature, AIC: 206.2
# Weighted_CPUE_Kōura ~ Emergent_Native + Overhanging_trees + Temperature, AIC: 205.5
# Weighted_CPUE_Kōura ~ pH + Temperature + Overhanging_trees + Presence_Morihana + Boulders, AIC: 204.4  based on only RF
# Weighted_CPUE_Kōura ~ Boulders + Emergent_Native + pH + Presence_Morihana + Temperature, AIC: 202.84755
# Weighted_CPUE_Kōura ~ Temperature + Emergent_Native + Presence_Common_smelt + pH + Cobble + Overhanging_trees, AIC 199.7
# Weighted_CPUE_Kōura ~ Temperature + Emergent_Native + Presence_Common_smelt + pH + Boulders + Cobble, AIC: 198.7


#

# 1. Keep only relevant variables 
M_C_data_subset <- M_C_data %>%
  dplyr::select(all_of(c(
    "Weighted_CPUE_Kōura", # response "Weighted_BCUE_Kōura", "Presence_Kōura"
    "LID",                 
    "Mean_depth_m","TLI",# "Elevation_m",
    "Slope_5m","Riparian_vegetation","Overhanging_trees","Erosion",
    "Wood_cover","Bedrock","Boulders","Cobble","Gravel","Sand","Organic_matter","Mud", "Presence_rocks",#"Turf",
    "Temperature","DO_mgl","DO_percent","Specific_conductivity","pH",
    "Emergent_Native","Emergent_Non_Native","Submerged_Native","Submerged_Non_Native",
    "Presence_Kōaro","Presence_Common_smelt","Presence_Trout","Presence_Morihana","Presence_Eel","Presence_Catfish"
  ))) %>%  na.omit()

# 2. Scale the data 
M_C_data_scaled <- M_C_data_subset
M_C_data_scaled <- M_C_data_subset %>%
  mutate(across(where(is.numeric) & !c(Weighted_CPUE_Kōura), scale)) %>%
  mutate(LID = factor(LID))

# Make sure LID is a factor
M_C_data_scaled$LID <- as.factor(M_C_data_scaled$LID)






response <- "Weighted_CPUE_Kōura"
glmmTMB_reduced <- remove_high_vif_glmmTMB(M_C_data_scaled, response, predictors)

model_summary <- summary(glmmTMB_reduced)


significant_vars <- rownames(coef(summary(glmmTMB_reduced))$cond)[coef(summary(glmmTMB_reduced))$cond[, "Pr(>|z|)"] < 0.10 & rownames(coef(summary(glmmTMB_reduced))$cond) != "(Intercept)"]
print(significant_vars)

coefs <- model_summary$coefficients$cond
raw_pvalues <- coefs[, "Pr(>|z|)"]
adjusted_p <- p.adjust(raw_pvalues, method = "fdr")
coefs_df <- as.data.frame(coefs)
coefs_df$adjusted_p <- adjusted_p
coefs_df_sorted <- coefs_df[order(coefs_df$adjusted_p), ]
print(coefs_df_sorted)

significant_vars_p <- rownames(coefs_df)[coefs_df$adjusted_p < 0.1]
print(significant_vars_p)

significant_formula <- as.formula(paste("Weighted_CPUE_Kōura ~", paste(significant_vars_p, collapse = " + "), "+ (1|LID)"))

reduced_model <- glmmTMB(
  significant_formula,
  data = M_C_data_scaled,
  family = tweedie(link = "log"))

summary(reduced_model)
check_collinearity(reduced_model)



best_model <- get_best_model(reduced_model)
summary(best_model)




# Fit GAM model
numeric_vars <- c("Mean_depth_m", "Slope_5m", "DO_mgl", "pH")
substrate_vars <- c("Wood_cover", "Bedrock", "Boulders", "Cobble", "Gravel", "Sand")
vegetation_vars <- c("Riparian_vegetation", "Overhanging_trees", "Erosion", "Emergent_Native", "Emergent_Non_Native", "Submerged_Native", "Submerged_Non_Native")
presence_vars <- c("Presence_rocks", "Presence_Kōaro", "Presence_Common_smelt", "Presence_Trout", "Presence_Morihana", "Presence_Eel", "Presence_Catfish")

smooth_terms <- paste0("s(", numeric_vars, ", bs = 'ts', k = 5)", collapse = " + ")
linear_terms <- paste(c(substrate_vars, vegetation_vars, presence_vars), collapse = " + ")

# Combine into full formula
gam_formula <- as.formula(paste("Weighted_CPUE_Kōura ~", smooth_terms, "+", linear_terms))

# Fit GAM
gam_model <- gam(gam_formula, data = M_C_data_scaled, family = tw(), select = TRUE)
summary(gam_model)

gam_model<- gam(Weighted_CPUE_Kōura ~ Cobble + s(Slope_5m, k = 3)+
                  Boulders #  Overhanging_trees +   +
                #s(Elevation_m,k = 3) + s(TLI, k = 3)
                , data = M_C_data_scaled, family = tw(), select = TRUE)

gam_formula2 <- as.formula(
  paste("Weighted_CPUE_Kōura ~ s(Slope_5m, k=3) + Cobble + Wood_cover + Riparian_vegetation + Presence_rocks")
)
gam_model2 <- gam(gam_formula2, data = M_C_data_scaled, family = tw(), select = TRUE)
summary(gam_model2)



# Fit LASSO mixed model
lasso_model <- glmmLasso(
  fix = fixed_formula,
  rnd = list(LID = ~1),
  data = M_C_data_scaled,
  lambda = 10,
  family = poisson(link = "log"),
  switch.NR = TRUE)

# Extract non-zero fixed effect coefficients
coef_fixed <- lasso_model$coefficients[-1]  # remove intercept
selected_vars <- names(coef_fixed)[coef_fixed != 0]

# Print selected variables
print(selected_vars)

# Build final formula
final_formula <- reformulate(
  selected_vars[selected_vars != "(Intercept)"], response = "log1p(Weighted_CPUE_Kōura)")

# Build the fixed effects part of the formula
fixed_part <- paste(selected_vars, collapse = " + ")

# Build the full formula with random effect
full_formula <- as.formula(
  paste("log1p(Weighted_CPUE_Kōura) ~", fixed_part, "+ (1 | LID)"))

# Fit final LMM
final_model <- lmer(full_formula, data = M_C_data_scaled)

summary(final_model)

check_collinearity(final_model)

# Test random effect
ranova(final_model)

# test the best models
final_model <- lmer(log1p(Weighted_CPUE_Kōura) ~ 
                         Boulders +  Elevation_m + #Temperature + 
                         TLI +
                         Slope_5m + Overhanging_trees + Riparian_vegetation +
                         Presence_Morihana+
                         (1 | LID), data = M_C_data_scaled)
summary(final_model)

best_final_model <- dredge(final_model)

# make model average of the top best models
averaged_model <- model.avg(best_final_model, subset = delta < 4) # average models within delta < 4
summary(averaged_model)

final_model_reduced <- lmer(log1p(Weighted_CPUE_Kōura) ~ 
                         Boulders +  #Elevation_m + Temperature  + TLI +
                         #Slope_5m + Overhanging_trees +   
                        # Presence_Morihana+
                         (1 | LID), data = M_C_data_scaled)

# 1. Refit models with ML (not REML) for valid comparison
final_model_ml <- update(final_model, REML = FALSE)
final_model_reduced_ml <- update(final_model_reduced, REML = FALSE)

# Compare
anova(final_model_ml, final_model_reduced_ml)

r.squaredGLMM(final_model_ml)
r.squaredGLMM(final_model_reduced_ml)

# Test random effect
ranova(final_model)

# Test the LM models
final_model_fixed <- lm(log1p(Weighted_CPUE_Kōura) ~ 
                          Boulders + Elevation_m + Temperature + TLI + Slope_5m +
                          Overhanging_trees + Riparian_vegetation + Presence_Morihana,
                        data = M_C_data_scaled)
summary(final_model_fixed)

best_fixed_model <- dredge(final_model_fixed)
averaged_model <- model.avg(best_fixed_model, subset = delta < 4) # average models within delta < 4
summary(averaged_model)

plot(final_model_fixed)

final_model_fixed_reduced <- lm(log1p(Weighted_CPUE_Kōura) ~ 
                           Boulders, #Elevation_m + Temperature + TLI + Slope_5m +
                          #Overhanging_trees + Riparian_vegetation +Presence_Morihana,
                        data = M_C_data_scaled)
summary(final_model_fixed_reduced)
plot(final_model_fixed_reduced)

anova(final_model_fixed, final_model_fixed_reduced)
r.squaredGLMM(final_model_fixed)
r.squaredGLMM(final_model_fixed_reduced)


model_summary <- summary(final_model_fixed_reduced)

# Extract R^2 and p-value
r2 <- round(model_summary$r.squared, 3)
pval <- signif(coef(model_summary)[2, 4], 3)

# Plot with annotation
ggplot(M_C_data_scaled, aes(Boulders, log1p(Weighted_CPUE_Kōura))) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  labs(x = "Boulders (scaled)",
    y = "log(Kōura CPUE + 1)") +
  annotate("text", x = 1, y = 1.5,
    label = paste0("R² = ", r2, "\nP = ", pval),
    size = 5, hjust = 0)+
  theme_bw()

ggplot(M_C_data_subset, aes(Weighted_CPUE_Kōura, Slope_5m )) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE)


# Fit a GLM 
library(tweedie)
library(statmod)  
library(cplm)
library(MASS)
glm_model <- glm(Weighted_CPUE_Kōura ~ Boulders + Elevation_m + TLI + Slope_5m + Overhanging_trees + Riparian_vegetation + Presence_Morihana,
                 data = M_C_data_scaled,
                 family = tweedie(var.power = 1.5, link.power = 0),
                 control = glm.control(maxit = 50))

summary(glm_model)
vif(glm_model)

tweedie.profile(Weighted_CPUE_Kōura ~ Boulders, 
                data = M_C_data_scaled, 
                p.vec = seq(1.1, 1.9, 0.05), do.plot = TRUE)

best_glm_model <- dredge(glm_model)
averaged_model <- model.avg(best_glm_model, subset = delta < 4) # average models within delta < 4
summary(averaged_model)

library(mgcv)
# Fit a GAM with a smooth term for Boulders
gam_model <- gam(log1p(Weighted_CPUE_Kōura) ~ s(Boulders, k = 3), 
                 data = M_C_data_scaled, 
                 method = "REML")

summary(gam_model)
plot(gam_model, residuals = TRUE, pch = 16)


# Step 2: Perform Spearman’s Rank Correlation test for each variable 
results_df <- purrr::map_dfr(variables, function(var) {
  test <- cor.test(M_C_data_subset[[var]], M_C_data_subset$Weighted_CPUE_Kōura, method = "spearman")
  tibble(Variable = var,Spearman_rho = test$estimate,p_value = test$p.value)})

# View sorted results
results_df %>% arrange(p_value)

# 2. Correlation matrix 
#numeric_vars <- M_C_data_subset[, sapply(M_C_data_subset, is.numeric)]
cor_results <- sapply(names(M_C_data_subset), function(var) {
  cor.test(M_C_data_subset[[var]], M_C_data_subset$Weighted_CPUE_Kōura, method = "spearman")$p.value})

adjusted_p <- p.adjust(cor_results, method = "fdr") # or method = "bonferroni"
significant_vars <- names(adjusted_p[adjusted_p < 0.1])
print(significant_vars)

cor_matrix_initial <- cor(M_C_data_subset[, sapply(M_C_data_subset, is.numeric)], use = "complete.obs", method = "spearman")
cor_melted_initial <- reshape2::melt(cor_matrix_initial)
ggplot(cor_melted_initial, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 2)), color = "black", size = 3) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Spearman Correlation Heatmap", x = "Variables", y = "Variables")


hist(M_C_data_scaled$Weighted_CPUE_Kōura)
hist(log1p(M_C_data_scaled$Weighted_CPUE_Kōura))


# make Test mix effect models
# test all variables
glmm_Presence0 <- lmer(log1p(Weighted_CPUE_Kōura) ~ 
                          TLI + Elevation_m + Mean_depth_m + 
                          Slope_5m + Riparian_vegetation + Overhanging_trees + Wood_cover + Erosion +
                          Bedrock + Boulders + Cobble + Gravel + Sand + Organic_matter +  Mud + Presence_rocks + 
                          Temperature + DO_mgl + pH + Specific_conductivity + DO_percent +
                          Emergent_Native + Emergent_Non_Native + Submerged_Native + Submerged_Non_Native + 
                          Presence_Morihana + Presence_Eel + Presence_Catfish +
                          Presence_Kōaro + Presence_Common_smelt + Presence_Trout + Presence_Mosquitofish +
                          (1 | LID), data = M_C_data_scaled)

# remove variables with high Correlation
glmm_Presence1 <- lmer(log1p(Weighted_CPUE_Kōura) ~ 
                         #TLI + 
                         Elevation_m + Mean_depth_m + 
                         Slope_5m + Riparian_vegetation + Overhanging_trees + Wood_cover + Erosion +
                         Bedrock + Boulders + Cobble + Gravel + Sand + Organic_matter +  Mud + Presence_rocks + 
                         Temperature + DO_mgl + pH + Specific_conductivity + #DO_percent +
                         Emergent_Native + Emergent_Non_Native + Submerged_Native + Submerged_Non_Native + 
                         Presence_Morihana + Presence_Eel + Presence_Catfish +
                         Presence_Kōaro + Presence_Common_smelt + Presence_Trout + Presence_Mosquitofish +
                         (1 | LID), data = M_C_data_scaled)


# remove variables with high Correlation
glmm_Presence2 <- lmer(log1p(Weighted_CPUE_Kōura) ~ 
                          TLI + 
                          Elevation_m + 
                          Mean_depth_m + 
                          Slope_5m +  Overhanging_trees + 
                          Boulders + Mud + 
                          Temperature  + 
                          DO_percent + 
                          Presence_Morihana + Presence_Eel + 
                          (1 | LID), data = M_C_data_scaled)

# test the significant results < 0.1
glmm_Presence3 <- lmer(log1p(Weighted_CPUE_Kōura) ~ 
                         Elevation_m + TLI +
                         Slope_5m + Overhanging_trees + Boulders + Temperature  + Presence_Morihana +
                         (1 | LID), data = M_C_data_scaled)

# test the Spearman’s Rank variables 
glmm_Presence4 <- lmer(log1p(Weighted_CPUE_Kōura) ~ 
                         Boulders +  #Elevation_m + 
                         Temperature  + #TLI +
                         #Slope_5m + Overhanging_trees +   
                         Presence_Morihana+
                          (1 | LID), data = M_C_data_scaled)

# test the significant results < 0.1 redefined
glmm_Presence5 <- lmer(log1p(Weighted_CPUE_Kōura) ~ 
                         Boulders +  Elevation_m + Temperature  + TLI +
                         Slope_5m + Overhanging_trees + Riparian_vegetation +
                         Presence_Morihana+
                         (1 | LID), data = M_C_data_scaled)

summary(glmm_Presence0)
summary(glmm_Presence1)
summary(glmm_Presence2)
summary(glmm_Presence3)
summary(glmm_Presence4)
summary(glmm_Presence5)

# check collinearity
library(performance)
check_collinearity(glmm_Presence0)
check_collinearity(glmm_Presence1)
check_collinearity(glmm_Presence2)
check_collinearity(glmm_Presence3)
check_collinearity(glmm_Presence4)
check_collinearity(glmm_Presence5)

# find the most parsimonious model
library(MuMIn)
options(na.action = "na.fail") # required for dredge
best_model0 <- dredge(glmm_Presence0)
best_model1 <- dredge(glmm_Presence1)
best_model2 <- dredge(glmm_Presence2)
best_model3 <- dredge(glmm_Presence3)
best_model4 <- dredge(glmm_Presence4)
best_model5 <- dredge(glmm_Presence5)

# make model average of the top best models
averaged_model <- model.avg(best_model4, subset = delta < 4) # average models within delta < 4
summary(averaged_model)


# Compare
anova(glmm_Presence4, glmm_Presence5)
AIC(glmm_Presence4, glmm_Presence5)

# Fit model without random effect
glmm_no_RE <- lm(log1p(Weighted_CPUE_Kōura) ~ Boulders, data = M_C_data_scaled)
# Fit model with random effect
glmm_with_RE <- lmer(log1p(Weighted_CPUE_Kōura) ~ Boulders + (1 | LID), data = M_C_data_scaled)
# Compare models
anova(glmm_with_RE, glmm_no_RE)

# Bayesian version
library(brms)
brm_Presence3 <- brm(
  formula = log1p(Weighted_CPUE_Kōura) ~ 
    Boulders + #Presence_Morihana+ #Elevation_m + TLI +
    #Slope_5m + Overhanging_trees +  
    Temperature  + #Riparian_vegetation +
    (1 | LID),
  data = M_C_data_scaled,
  family = gaussian(),
  chains = 4,
  iter = 4000,
  warmup = 2000,
  cores = 4,
  seed = 123,
  control = list(adapt_delta = 0.9999, max_treedepth = 15)
)


# Compute LOO for each model
loo1 <- brms::loo(brm_Presence1)
loo2 <- brms::loo(brm_Presence2)
loo3 <- brms::loo(brm_Presence3)

# Compare models
loo_compare(loo1, loo2, loo3)



# Load libraries
library(optimx)
library(projpred)

# Fit the full model
full_model <- brm(
  formula = log1p(Weighted_CPUE_Kōura) ~ 
    Boulders + Presence_Morihana + Elevation_m + TLI +
    Slope_5m + Overhanging_trees + Riparian_vegetation + Temperature + 
    (1 | LID),
  data = M_C_data_scaled,
  family = gaussian(),
  chains = 4,
  iter = 4000,
  warmup = 2000,
  cores = 4,
  seed = 123,
  control = list(adapt_delta = 0.99999, max_treedepth = 15))

# Run variable selection
vs <- varsel(full_model, method = "forward") # or "backward"

# Inspect selection path
plot(vs)  # Shows how predictive accuracy changes as predictors are added

# Find the suggested model size
best_size <- suggest_size(vs)
print(best_size)

# Project onto the submodel with the suggested size
submodel <- project(vs, nv = best_size)

# Look at summary of reduced model
summary(submodel)

# Get the selected predictor terms
selected_terms <- predictor_terms(submodel)
print(selected_terms)

# Build the reduced formula manually
reduced_formula <- reformulate(
  termlabels = selected_terms,
  response = "log1p(Weighted_CPUE_Kōura)"
)

cat("Reduced model formula:\n")
print(reduced_formula)



brm_reduced <- brm(
  formula = reduced_formula,
  data = M_C_data_scaled,
  family = gaussian(),
  chains = 4,
  iter = 4000,
  warmup = 2000,
  cores = 4,
  seed = 123,
  control = list(adapt_delta = 0.9999, max_treedepth = 15)
)

summary(brm_reduced)


# Compare full and reduced models
loo_full <- loo(full_model)
loo_reduced <- loo(brm_reduced)

# Model comparison
loo_compare(loo_full, loo_reduced)








# explanation of dropping mixed effect model:
library(lmerTest)
library(lmtest)
# Fixed-effects model
fixed_model <- lm(log1p(Weighted_CPUE_Kōura) ~ Elevation_m + Boulders + Cobble + Mud, data = M_C_data_scaled)

# Mixed-effects model
mixed_model <- lmer(log1p(Weighted_CPUE_Kōura) ~ Elevation_m + Boulders + Cobble + Mud + (1|LID), data = M_C_data_scaled)

# Check if random effect variance is zero
summary(mixed_model)
isSingular(mixed_model)

# Likelihood ratio test
lrtest(mixed_model, fixed_model)

ranova(mixed_model)


# 3.5 Filter for Multicollinearity Check (VIF) using mixed effect model
fixed_formula <- as.formula(paste("Weighted_CPUE_Kōura ~", paste(num_cols, collapse = " + "), "+ (1 | LID)"))

full_model <- lmer(fixed_formula, data = M_C_data_scaled)
summary(full_model)

vif_values <- vif(full_model)
print(vif_values)

# test several mixed models
model1 <- log1p(Weighted_CPUE_Kōura) ~ 
  LID +
  Elevation_m + 
  Slope_5m + 
  Boulders + 
  Cobble + 
  Gravel +
  Mud + 
  Emergent_Native + 
  Emergent_Non_Native +
  Presence_Morihana + 
  Temperature +
  (1 | LID)

model2 <- log1p(Weighted_CPUE_Kōura) ~ 
  Elevation_m + 
  Slope_5m + 
  Boulders + 
  Cobble + 
  Gravel +
  Mud + 
  Emergent_Native + 
  #Emergent_Non_Native +
  #Presence_Morihana +
  (1 | LID)

model3 <- log1p(Weighted_CPUE_Kōura) ~ 
  Elevation_m + 
  Slope_5m + 
  Boulders + 
  Cobble + 
  Gravel +
  Mud + 
  Emergent_Native + 
  #Emergent_Non_Native +
  #Presence_Morihana +
  Temperature +
  (1|LID)


lmm_CPUE1 <- lmer(model1, data = M_C_data_scaled)
lmm_CPUE2 <- lmer(model2, data = M_C_data_scaled)
lmm_CPUE3 <- lmer(model3, data = M_C_data_scaled)

summary(lmm_CPUE1)
summary(lmm_CPUE2)
summary(lmm_CPUE3)

# Compare models using likelihood ratio test
anova(lmm_CPUE1,lmm_CPUE2, lmm_CPUE3, test = "Chisq")


lm_CPUE <- lm(model1, data = M_C_data_scaled)
summary(lm_CPUE)

par(mfrow=c(2,2))
plot(lm_CPUE)


# 4. PCA
selected_vars
predictors <- M_C_data_scaled[, c("Elevation_m","TLI","Slope_5m", "Riparian_vegetation", 
                                  "Overhanging_trees","Boulders", "Temperature", "Presence_Morihana")]
# Perform PCA
pca_result <- prcomp(predictors, center = TRUE, scale. = TRUE)

# Summary of PCA
summary(pca_result)

# Make plots
library(factoextra)
fviz_eig(pca_result)                 # Variance explained
fviz_pca_biplot(pca_result)         # Biplot
fviz_pca_ind(pca_result)             # Samples on PCA plane
fviz_pca_var(pca_result)



#
# Test for weighted BCUE kōura (biomass) --------------------------------------
# 1. Preprocessing 
vars <- c("Weighted_BCUE_Kōura", "Slope_5m", "Riparian_vegetation", "Overhanging_trees",
          "Wood_cover", "Bedrock", "Boulders", "Cobble", "Gravel", "Sand", "Organic_matter",
          "Temperature", "DO_percent", "pH", "Specific_conductivity", "DO_mgl",
          "Emergent_Native", "Submerged_Native", "Submerged_Non_Native", "Emergent_Non_Native",
          "Presence_Kōaro", "Presence_Common_smelt", "Presence_Morihana", 
          "Presence_Eel", "Presence_Catfish")

M_C_data_scaled2 <- M_C_data %>%
  dplyr::select(all_of(vars)) %>%
  na.omit() %>%
  mutate(across(where(is.numeric) & !c(Weighted_BCUE_Kōura), scale))

# 2. Correlation check 
cor_matrix <- cor(dplyr::select(M_C_data_scaled2, where(is.numeric)), use = "complete.obs")
corrplot(cor_matrix)

# 3. VIF reduction 
response <- "Weighted_BCUE_Kōura"
predictors <- setdiff(names(M_C_data_scaled2), response)

model <- remove_high_vif_glmmTMB(M_C_data_scaled2, response, predictors)
performance::check_collinearity(model)

# 4. Post-VIF correlation filtering 
remaining_vars <- c(setdiff(names(fixef(model)$cond), "(Intercept)"), response)
cor_filtered <- M_C_data_scaled2 %>% dplyr::select(all_of(remaining_vars))
corrplot(cor(cor_filtered, use = "complete.obs"))

# 5. Build GAM formulas 
# Custom smooth parameters
custom_k <- list(Boulders=5, Cobble=4, Emergent_Native=5, Overhanging_trees=5,
                 Riparian_vegetation=7, Bedrock=3, Gravel=6, Submerged_Native=4,
                 Submerged_Non_Native=9)

# Auto formula generator
build_smooth <- function(var) {
  n <- length(unique(M_C_data_scaled2[[var]]))
  if (n >= 5) {
    k <- if (!is.null(custom_k[[var]])) paste0(", k=", custom_k[[var]]) else ""
    paste0("s(", var, ", bs='ts'", k, ")")
  } else var
}
formula_1 <- as.formula(paste(response, "~", 
                              paste(sapply(setdiff(remaining_vars, response), build_smooth), 
                                    collapse = " + ")))

# Manually curated formulas (based on p-values)
formula_2 <- Weighted_BCUE_Kōura ~
  s(Overhanging_trees, bs = "ts", k = 5) +
  s(Cobble, bs = "ts", k = 4) +
  s(pH, bs = "ts") +
  Presence_Kōaro +
  Presence_Eel

formula_3 <- update(formula_2, . ~ . - Presence_Kōaro)
formula_4 <- update(formula_3, . ~ . - s(pH, bs = "ts"))
#formula_5 <- update(formula_3, . ~ . - s(Cobble, bs = "ts", k = 4))

# 6. Fit GAMs 
m1 <- gam(formula_1, family = tw(link = "log"), method = "ML", data = M_C_data_scaled2)
m2 <- gam(formula_2, family = tw(link = "log"), method = "ML", data = M_C_data_scaled2)
m3 <- gam(formula_3, family = tw(link = "log"), method = "ML", data = M_C_data_scaled2)
m4 <- gam(formula_4, family = tw(link = "log"), method = "ML", data = M_C_data_scaled2)
#m5 <- gam(formula_5, family = tw(link = "log"), method = "ML", data = M_C_data_scaled2)

summary(m2)

# 7. Model summaries & comparison 
lapply(list(m1, m2, m3, m4), summary)
anova(m1, m2, m3, m4, test = "Chisq")
AIC(m1, m2, m3, m4)

# 8. Final model diagnostics & plotting 
m2.1 <- gam(formula_2, family = tw(link = "log"), method = "REML", data = M_C_data_scaled2)
summary(m2.1)
gam.check(m2.1)
concurvity(m2.1)
influence.gam(m2.1)

# Visual diagnostics
plot(m2.1, pages = 1, shade = TRUE)
draw(m2.1)

# Predicted vs observed
plot(predict(m2.1, type = "response"), m2.1$y,
     xlab = "Predicted CPUE", ylab = "Observed CPUE")
abline(0, 1, col = "red")





#
# OLD ####
# 1. Preprocess
M_C_data_scaled1 <- M_C_data %>%
  dplyr::select(all_of(c(
    "Weighted_BCUE_Kōura", "LID",
    "Mean_depth_m", #"TLI", 
    "Slope_5m", "Riparian_vegetation", "Overhanging_trees", "Erosion",
    "Wood_cover", "Bedrock", "Boulders", "Cobble", "Gravel", "Sand", "Organic_matter", "Mud",
    "Presence_rocks", "Temperature", "DO_mgl", "DO_percent", "Specific_conductivity", "pH",
    "Emergent_Native", "Emergent_Non_Native", "Submerged_Native", "Submerged_Non_Native",
    "Presence_Kōaro", "Presence_Common_smelt", "Presence_Trout", "Presence_Morihana",
    "Presence_Eel", "Presence_Catfish"
  ))) %>%
  na.omit() %>%
  mutate(across(where(is.numeric) & !c(Weighted_BCUE_Kōura), scale),
         LID = factor(LID))

# Explore correlations
cor_data <- M_C_data_scaled1 %>% select(where(is.numeric), -LID)
cor_matrix <- cor(cor_data, use = "complete.obs")
corrplot(cor_matrix)

# 2. VIF-based reduction
response <- "Weighted_BCUE_Kōura"
predictors <- setdiff(names(M_C_data_scaled1), c("Weighted_BCUE_Kōura", "LID"))
model <- remove_high_vif_glmmTMB(M_C_data_scaled1, response, predictors)

# 3. Adjusted p-values
coefs_df <- as.data.frame(summary(model)$coefficients$cond)
coefs_df$adj_p <- p.adjust(coefs_df$`Pr(>|z|)`, method = "fdr")
sig_vars <- rownames(coefs_df)[coefs_df$`Pr(>|z|)` < 0.05 & rownames(coefs_df) != "(Intercept)"]
sig_vars_adj <- rownames(coefs_df)[coefs_df$adj_p < 0.2 & rownames(coefs_df) != "(Intercept)"]

# 4. Reduced model
final_formula <- as.formula(paste(response, "~", paste(sig_vars, collapse = " + "), "+ (1|LID)"))
reduced_model <- glmmTMB(final_formula, data = M_C_data_scaled1, family = tweedie(link = "log"))

# 5. Model selection
best_model <- get_best_model(reduced_model)
summary(best_model)


# 6. Model visualization
# Simulate residuals with DHARMa and plot diagnostics
sim_res <- simulateResiduals(fittedModel = best_model, plot = TRUE)

# Extract and print R2 values (marginal and conditional)
r2_vals <- r2(best_model, by_group = TRUE)
print(r2_vals)

# Plot predicted vs observed
predicted <- predict(best_model, type = "response")  
observed <- M_C_data_scaled1$Weighted_BCUE_Kōura     

plot(predicted, observed)
abline(a = 0, b = 1, col = "red", lwd = 2)





# 1. Keep only relevant variables 
M_C_data_subset2 <- M_C_data %>%
  select(all_of(c(
    "Weighted_BCUE_Kōura", # response "Weighted_CPUE_Kōura", "Presence_Kōura"
    "LID",                 
    "Mean_depth_m","Elevation_m","TLI",
    "Slope_5m","Riparian_vegetation","Overhanging_trees","Erosion",
    "Wood_cover","Bedrock","Boulders","Cobble","Gravel","Sand","Organic_matter","Mud","Presence_rocks",
    "Temperature","DO_mgl","DO_percent","Conductivity","Specific_conductivity","pH",
    "Emergent_Native","Emergent_Non_Native","Submerged_Native","Submerged_Non_Native",
    "Presence_Kōaro","Presence_Common_smelt","Presence_Trout","Presence_Morihana","Presence_Eel","Presence_Catfish","Presence_Mosquitofish"
  ))) %>%  na.omit()

variables2<-c("Weighted_BCUE_Kōura","LID",
             "TLI","Elevation_m","Mean_depth_m",
             "Slope_5m","Riparian_vegetation","Overhanging_trees","Wood_cover","Erosion",
             "Bedrock","Boulders","Cobble","Gravel","Sand","Organic_matter","Mud","Presence_rocks",
             "Temperature","DO_mgl", "Specific_conductivity","pH","DO_percent",
             "Emergent_Native","Emergent_Non_Native","Submerged_Native","Submerged_Non_Native",
             "Presence_Kōaro","Presence_Common_smelt","Presence_Trout","Presence_Morihana","Presence_Eel","Presence_Catfish","Presence_Mosquitofish")


# 2. Scale the data 
M_C_data_scaled2 <- M_C_data_subset2
num_cols <- setdiff(names(M_C_data_scaled2), c("Weighted_BCUE_Kōura", "LID"))
num_cols <- num_cols[sapply(M_C_data_scaled2[num_cols], is.numeric)]
M_C_data_scaled2[num_cols] <- lapply(M_C_data_scaled2[num_cols], scale)

# Make sure LID is a factor
M_C_data_scaled2$LID <- as.factor(M_C_data_scaled2$LID)

# Get predictor names
predictors <- setdiff(names(M_C_data_scaled2), c("Weighted_BCUE_Kōura", "LID"))

# Build formula
fixed_formula <- as.formula(
  paste("log1p(Weighted_BCUE_Kōura) ~", paste(predictors, collapse = " + ")))

# Fit LASSO mixed model
lasso_model <- glmmLasso(
  fix = fixed_formula,
  rnd = list(LID = ~1),
  data = M_C_data_scaled2,
  lambda = 10,
  family = gaussian(link = "identity"),
  switch.NR = TRUE)

# Extract non-zero fixed effect coefficients
coef_fixed <- lasso_model$coefficients[-1]  # remove intercept
selected_vars <- names(coef_fixed)[coef_fixed != 0]

# Print selected variables
print(selected_vars)

# Build final formula
final_formula <- reformulate(
  selected_vars[selected_vars != "(Intercept)"], response = "log1p(Weighted_BCUE_Kōura)")

# Build the fixed effects part of the formula
fixed_part <- paste(selected_vars, collapse = " + ")

# Build the full formula with random effect
full_formula <- as.formula(
  paste("log1p(Weighted_BCUE_Kōura) ~", fixed_part, "+ (1 | LID)"))

# Fit final LMM
final_model <- lmer(full_formula, data = M_C_data_scaled2)

summary(final_model)

check_collinearity(final_model)




# test the best models
final_model <- lmer(log1p(Weighted_BCUE_Kōura) ~ 
                      Elevation_m +  Overhanging_trees + Boulders + Cobble +
                      Mud + Presence_Kōaro + Presence_Morihana +
                      Presence_Eel + Presence_Mosquitofish +
                      (1 | LID), data = M_C_data_scaled2)

best_final_model <- dredge(final_model)

# make model average of the top best models
averaged_model <- model.avg(best_final_model, subset = delta < 4) # average models within delta < 4
summary(averaged_model)

final_model_reduced <- lmer(log1p(Weighted_BCUE_Kōura) ~ 
                              Overhanging_trees + Presence_Morihana + Presence_Eel + 
                              (1 | LID), data = M_C_data_scaled2)
summary(final_model_reduced)

# 1. Refit models with ML (not REML) for valid comparison
final_model_ml <- update(final_model, REML = FALSE)
final_model_reduced_ml <- update(final_model_reduced, REML = FALSE)

# Compare
anova(final_model_ml, final_model_reduced_ml)


#
# Test for kōura occupancy (presence/absence) ---------------------------------
# 1. Preprocessing 
vars <- c("Presence_Kōura", "Lake",
          "Slope_5m", "Riparian_vegetation", "Overhanging_trees","Wood_cover", #"Substrate_index",
          "Bedrock","Boulders", "Cobble", "Gravel", "Sand", "Organic_matter",
          "Temperature", "DO_percent", "pH", "Specific_conductivity", "DO_mgl",
          "Emergent_Native","Submerged_Non_Native", "Turf_Native", #"Submerged_Native","Emergent_Non_Native",
          "Presence_Kōaro", "Presence_Common_smelt","Presence_Morihana","Presence_Eel") #"Presence_Catfish")

M_C_data_scaled3 <- M_C_data %>%
  dplyr::select(all_of(vars)) %>%
  na.omit() %>%
  mutate(across(where(is.numeric) & !c(Presence_Kōura))) #, scale

library(caret)
nearZeroVar(M_C_data_scaled3, saveMetrics = TRUE)

# 2. Correlation check 
cor_matrix <- cor(dplyr::select(M_C_data_scaled3, where(is.numeric)), use = "complete.obs")
corrplot(cor_matrix)

# 3. VIF reduction 
response <- "Presence_Kōura"
predictors <- setdiff(names(M_C_data_scaled3), response)

model <- remove_high_vif_glmBI(M_C_data_scaled3, response, predictors)
#summary(model)
performance::check_collinearity(model)

# 4. Post-VIF correlation filtering 
coef_names <- names(coef(model))
clean_vars <- setdiff(coef_names, "(Intercept)")
clean_vars <- c(clean_vars, response)
cor_filtered <- M_C_data_scaled3 %>% dplyr::select(all_of(clean_vars))
corrplot(cor(cor_filtered, use = "complete.obs"))
clean_vars <- setdiff(clean_vars, response)

# Mann–Whitney test
MW_test <- sapply(clean_vars, function(var) {
  wilcox.test(M_C_data_scaled3[[var]] ~ M_C_data_scaled3$Presence_Kōura, exact = FALSE)$p.value})
significant <- names(MW_test)[MW_test < 0.05]
#significant

# 5. Build GAM formulas 
# Custom smooth parameters
custom_k <- list(#LID=4, Mean_depth_m=4,
                 Overhanging_trees=5, Riparian_vegetation=7,Overhanging_trees=5,
                 Bedrock=3, Boulders=5, Cobble=4, Gravel=6,
                 Emergent_Native=5, Submerged_Native=4, Submerged_Non_Native=9, Turf_Native=6)

m_final2 <- stepwise_gam_vars(
  response = response,
  vars = clean_vars,
  data = M_C_data_scaled3,
  family = binomial,
  method = "ML",
  p_cutoff = 0.1)

summary(m_final2)
formula_final2 <- formula(m_final2)

formula_final2 <- response ~ Lake+ s(Overhanging_trees, bs = "ts", k = 5) + 
  s(Temperature, bs = "ts") + Presence_Morihana + Presence_Eel


m1 <- gam(formula_final2, family = binomial, method = "REML", data = M_C_data_scaled3)
summary(m1)
gam.check(m1)
concurvity(m1)
influence.gam(m1)

# Predicted vs observed
plot(predict(m1, type = "response"), m1$y,xlab = "Predicted presence", ylab = "Observed presence")
abline(0, 1, col = "red")

# Build new data frame 
newdata_param <- expand.grid(
  Presence_Morihana = c(0, 1),
  Presence_Eel = c(0, 1),
  #Riparian_vegetation = mean(M_C_data_scaled3$Riparian_vegetation, na.rm = TRUE),
  Overhanging_trees = mean(M_C_data_scaled3$Overhanging_trees, na.rm = TRUE),
  Temperature = mean(M_C_data_scaled3$Temperature, na.rm = TRUE))

newdata_param$predicted_prob <- predict(m1, newdata = newdata_param, type = "response")

p1 <- draw(m1)
p2 <- ggplot(newdata_param, aes(x = interaction(Presence_Morihana, Presence_Eel), y = predicted_prob)) +
  geom_col(width = 0.6) +
  labs(x = "Presence Morihana and Eel", y = "Predicted Probability of Presence_Kōura") 

library(patchwork)
p1 + p2 + plot_layout(ncol = 3)




#----
formula_0 <- as.formula(paste(response, "~", 
                              paste(sapply(setdiff(remaining_vars, response), build_smooth), 
                                    collapse = " + ")))

# Manually curated formulas (based on p-values)
formula_1 <- Presence_Kōura ~
  s(Riparian_vegetation, bs = "ts", k = 7) +
  s(Boulders, bs = "ts", k = 5) +
  s(DO_percent, bs = "ts") +
  Presence_Morihana +
  Presence_Eel

formula_7 <- Presence_Kōura ~
  s(Slope_5m, bs = "ts") +
  s(Boulders, bs = "ts", k = 5) +
  Presence_Morihana 


formula_8 <- Presence_Kōura ~
  #s(Slope_5m, bs = "ts") + 
  #s(Riparian_vegetation, bs = "ts", k = 7) + 
  #s(Overhanging_trees, bs = "ts", k = 5) + 
  #s(Wood_cover, bs = "ts") + 
  s(Boulders, bs = "ts", k = 5) + 
  #s(Cobble, bs = "ts", k = 4) + 
  #s(Gravel, bs = "ts", k = 6) + 
  #s(Sand, bs = "ts") + 
  s(DO_percent, bs = "ts") + 
  #s(pH, bs = "ts") + 
  #s(Specific_conductivity, bs = "ts") + 
  #s(Emergent_Native, bs = "ts", k = 5) +
  #s(Submerged_Native, bs = "ts", k = 4) + 
  #s(Submerged_Non_Native, bs = "ts", k = 9) + 
  #Emergent_Non_Native + #Presence_Kōaro + #Presence_Common_smelt + 
  Presence_Morihana + Presence_Eel
m8 <- gam(formula_8, family = binomial, method = "ML", data = M_C_data_scaled3)
summary(m8)

formula_2 <- update(formula_1, . ~ . - s(Riparian_vegetation, bs = "ts", k = 7))
formula_3 <- update(formula_2, . ~ . - s(Boulders, bs = "ts", k = 5))
formula_4 <- update(formula_2, . ~ . - s(DO_percent, bs = "ts"))
formula_5 <- update(formula_2, . ~ . - Presence_Morihana)
formula_6 <- update(formula_2, . ~ . - Presence_Eel)


# 6. Fit GAMs 
m0 <- gam(formula_0, family = binomial, method = "ML", data = M_C_data_scaled3)
m1 <- gam(formula_1, family = binomial, method = "ML", data = M_C_data_scaled3)
m2 <- gam(formula_2, family = binomial, method = "ML", data = M_C_data_scaled3)
m3 <- gam(formula_3, family = binomial, method = "ML", data = M_C_data_scaled3)
m4 <- gam(formula_4, family = binomial, method = "ML", data = M_C_data_scaled3)
m5 <- gam(formula_5, family = binomial, method = "ML", data = M_C_data_scaled3)
m6 <- gam(formula_6, family = binomial, method = "ML", data = M_C_data_scaled3)
m7 <- gam(formula_7, family = binomial, method = "ML", data = M_C_data_scaled3)

summary(m0)

# 7. Model summaries & comparison 
lapply(list(m0, m1, m2, m3, m4, m5, m6), summary)
anova(m0, m1, test = "Chisq")
anova(m1, m2, test = "Chisq")
anova(m2, m3, test = "Chisq")
anova(m2, m4, test = "Chisq")
anova(m2, m5, test = "Chisq")
anova(m2, m6, test = "Chisq")
AIC(m1, m2, m3, m4, m5, m6)

# 8. Final model diagnostics & plotting 
m2.1 <- gam(formula_2, family = binomial, method = "REML", data = M_C_data_scaled3)
summary(m2.1)
gam.check(m2.1)
concurvity(m2.1)
influence.gam(m2.1)

# Visual diagnostics
plot(m2.1, pages = 1, shade = TRUE)
draw(m2.1)

# Predicted vs observed
plot(predict(m2.1, type = "response"), m2.1$y,
     xlab = "Predicted presence", ylab = "Observed presence")
abline(0, 1, col = "red")


ggplot(M_C_data_scaled3, aes(as.factor(Presence_Morihana), fill = as.factor(Presence_Kōura))) +
  geom_bar(position = "dodge") +
  facet_wrap(~ as.factor(Presence_Eel), labeller = label_both) 

# Build new data frame 
newdata_param <- expand.grid(
  Presence_Morihana = c(0, 1),
  Presence_Eel = c(0, 1),
  #Riparian_vegetation = mean(M_C_data_scaled3$Riparian_vegetation, na.rm = TRUE),
  Boulders = mean(M_C_data_scaled3$Boulders, na.rm = TRUE),
  DO_percent = mean(M_C_data_scaled3$DO_percent, na.rm = TRUE))

newdata_param$predicted_prob <- predict(m2.1, newdata = newdata_param, type = "response")

p1 <- draw(m2.1)
p2 <- ggplot(newdata_param, aes(x = interaction(Presence_Morihana, Presence_Eel), y = predicted_prob)) +
  geom_col(width = 0.6) +
  labs(x = "Presence Morihana and Eel", y = "Predicted Probability of Presence_Kōura") 

library(patchwork)
p1 + p2 + plot_layout(ncol = 3)


ggplot(M_C_data, aes(DO_percent, Weighted_CPUE_Kōura))+
  geom_point(aes(col=Lake))+
  geom_smooth()












#
# OLD ----

# 1. Keep only relevant variables 
M_C_data_subset3 <- M_C_data %>%
  select(all_of(c("Presence_Kōura", "LID",
    "TLI","Elevation_m","Mean_depth_m",
    "Slope_5m","Riparian_vegetation","Overhanging_trees","Wood_cover","Erosion",
    "Bedrock","Boulders","Cobble","Gravel","Sand","Organic_matter","Mud","Presence_rocks",
    "Temperature","DO_mgl", "Specific_conductivity","pH","DO_percent",
    "Emergent_Native","Emergent_Non_Native","Submerged_Native","Submerged_Non_Native",
    "Presence_Kōaro","Presence_Common_smelt","Presence_Trout","Presence_Morihana","Presence_Eel","Presence_Catfish","Presence_Mosquitofish"))) %>%  na.omit()


# 2. Scale the data 
M_C_data_scaled3 <- M_C_data_subset3
num_cols <- setdiff(names(M_C_data_scaled3), c("Presence_Kōura", "LID"))
num_cols <- num_cols[sapply(M_C_data_scaled3[num_cols], is.numeric)]
M_C_data_scaled3[num_cols] <- lapply(M_C_data_scaled3[num_cols], scale)

# Make sure LID is a factor
M_C_data_scaled3$LID <- as.factor(M_C_data_scaled3$LID)

# Get predictor names
predictors <- setdiff(names(M_C_data_scaled3), c("Presence_Kōura", "LID"))

# Build formula
fixed_formula <- as.formula(
  paste("Presence_Kōura ~", paste(predictors, collapse = " + ")))

# Fit LASSO mixed model
lasso_model <- glmmLasso(
  fix = fixed_formula,
  rnd = list(LID = ~1),
  data = M_C_data_scaled3,
  lambda = 10,
  family = binomial(link = "logit"),
  switch.NR = TRUE)

# Extract non-zero fixed effect coefficients
coef_fixed <- lasso_model$coefficients[-1]  # remove intercept
selected_vars <- names(coef_fixed)[coef_fixed != 0]

# Print selected variables
print(selected_vars)


# Build final formula
final_formula <- reformulate(
  selected_vars[selected_vars != "(Intercept)"], response = "Presence_Kōura")

# Build the fixed effects part of the formula
fixed_part <- paste(selected_vars, collapse = " + ")

# Build the full formula with random effect
full_formula <- as.formula(
  paste("Presence_Kōura ~", fixed_part, "+ (1 | LID)"))

# Fit final LMM
final_model <- glmer(full_formula, data = M_C_data_scaled3, family = binomial(link = "logit"))

summary(final_model)

check_collinearity(final_model)










variables3 <-c("Presence_Kōura","LID",
               "TLI","Elevation_m","Mean_depth_m",
               "Slope_5m","Riparian_vegetation","Overhanging_trees","Wood_cover","Erosion",
               "Bedrock","Boulders","Cobble","Gravel","Sand","Organic_matter","Mud","Presence_rocks",
               "Temperature","DO_mgl", "Specific_conductivity","pH","DO_percent",
               "Emergent_Native","Emergent_Non_Native","Submerged_Native","Submerged_Non_Native",
               "Presence_Kōaro","Presence_Common_smelt","Presence_Trout","Presence_Morihana","Presence_Eel","Presence_Catfish","Presence_Mosquitofish")

# Step 2: Perform Mann–Whitney U test for each variable 
results <- lapply(variables3[-1], function(var) {
  wilcox.test(M_C_data_subset3[[var]] ~ M_C_data_subset3$Presence_Kōura   , data = M_C_data_subset3)})

# Format results into a summary
results_summary <- data.frame(
  Variable = variables3[-1],
  p_value = sapply(results, function(x) x$p.value),
  W_statistic = sapply(results, function(x) x$statistic))

# Display significant results (e.g., p < 0.05)
Wilcoxon_results <- results_summary %>%
  filter(p_value < 0.05)

print(Wilcoxon_results)


# scale variables
M_C_data_scaled3 <- M_C_data_subset3
num_cols <- setdiff(names(M_C_data_scaled3), c("Presence_Kōura", "LID"))
num_cols <- num_cols[sapply(M_C_data_scaled3[num_cols], is.numeric)]
M_C_data_scaled3[num_cols] <- lapply(M_C_data_scaled3[num_cols], scale)


# Check for collinearity
num_cols <- names(M_C_data_scaled3)[sapply(M_C_data_scaled3, is.numeric)]
cor_matrix <- cor(M_C_data_scaled3[num_cols], use = "pairwise.complete.obs")
cor_melted <- reshape2::melt(cor_matrix)
high_corr <- findCorrelation(cor_matrix, cutoff = 0.7, names = TRUE)  # From caret package

ggplot(cor_melted, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 2)), color = "black", size = 3) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Spearman Correlation Heatmap", x = "Variables", y = "Variables")


# make Test mix effect models
# test all variables
glmm_Presence0 <- glmer(Presence_Kōura ~ 
  TLI + Elevation_m + Mean_depth_m + 
  Slope_5m + Riparian_vegetation + Overhanging_trees + Wood_cover + Erosion +
  Bedrock + Boulders + Cobble + Gravel + Sand + Organic_matter +  Mud + Presence_rocks + 
  Temperature + DO_mgl + pH + Specific_conductivity + DO_percent + Conductivity +
  Emergent_Native + Emergent_Non_Native + Submerged_Native + Submerged_Non_Native + 
  Presence_Morihana + Presence_Eel + Presence_Catfish +
  Presence_Kōaro + Presence_Common_smelt + Presence_Trout + Presence_Mosquitofish +
  (1 | LID), data = M_C_data_scaled3, family = binomial(link = "logit"))

# test the Wilcoxon_results < .1
glmm_Presence1 <- glmer(Presence_Kōura ~ 
  TLI + Elevation_m + Mean_depth_m + 
  Slope_5m +  Overhanging_trees + #Riparian_vegetation +Wood_cover + Erosion +
  Boulders + Mud + #Bedrock +Cobble + Gravel + Sand + Organic_matter + Presence_rocks + 
  Temperature + Conductivity + DO_percent + #DO_mgl + pH + Specific_conductivity +  
  #Emergent_Native + Emergent_Non_Native + Submerged_Native + Submerged_Non_Native + 
  Presence_Morihana + Presence_Eel + #Presence_Catfish + Presence_Kōaro + Presence_Common_smelt + Presence_Trout + Presence_Mosquitofish +
  (1 | LID), data = M_C_data_scaled3, family = binomial(link = "logit"))


# remove variables with high Correlation
glmm_Presence2 <- glmer(Presence_Kōura ~ 
  #TLI + Elevation_m + 
  Mean_depth_m + 
  Slope_5m +  Overhanging_trees + 
  Boulders + Mud + 
  Temperature + #Conductivity + 
  DO_percent + 
  Presence_Morihana + Presence_Eel + 
  (1 | LID), data = M_C_data_scaled3, family = binomial(link = "logit"))

# remove diffrent variables than in model2
glmm_Presence3 <- glmer(Presence_Kōura ~ 
  #TLI + 
  #Elevation_m + 
  #Mean_depth_m + 
  Slope_5m +  Overhanging_trees + 
  Boulders + #Mud + 
  Temperature + 
  #Conductivity + DO_percent + 
  Presence_Morihana + Presence_Eel + 
  (1 | LID), data = M_C_data_scaled3, family = binomial(link = "logit"))

# test the Wilcoxon results < .5
glmm_Presence4 <- glmer(Presence_Kōura ~ 
  #TLI + Elevation_m + 
  Mean_depth_m + 
  Slope_5m +  Boulders +  
  #Temperature +  
  Presence_Morihana + 
  (1 | LID), data = M_C_data_scaled3, family = binomial(link = "logit"))

summary(glmm_Presence0)
summary(glmm_Presence1)
summary(glmm_Presence2)
summary(glmm_Presence3)
summary(glmm_Presence4)
summary(glmm_Presence5)

# check collinearity
library(performance)
check_collinearity(glmm_Presence0)
check_collinearity(glmm_Presence1)
check_collinearity(glmm_Presence2)
check_collinearity(glmm_Presence3)
check_collinearity(glmm_Presence4)
check_collinearity(glmm_Presence5)

# find the most parsimonious model
library(MuMIn)
options(na.action = "na.fail") # required for dredge
best_model0 <- dredge(glmm_Presence0)
best_model1 <- dredge(glmm_Presence1)
best_model2 <- dredge(glmm_Presence2)
best_model3 <- dredge(glmm_Presence3)
best_model4 <- dredge(glmm_Presence4)
best_model5 <- dredge(glmm_Presence5)

# make model average of the top best models
averaged_model <- model.avg(best_model4, subset = delta < 4) # average models within delta < 4
summary(averaged_model)


# Compare models using likelihood ratio test
anova(glmm_Presence0, glmm_Presence1, glmm_Presence2, glmm_Presence3, test = "Chisq")


# Bayesian version
library(brms)
brm_Presence1 <- brm(
  formula = Presence_Kōura| trials(1) ~ 
    TLI + Elevation_m + 
    Mean_depth_m + 
    Slope_5m +  Boulders +  
    Temperature +  
    Presence_Morihana + 
    (1 | LID),
  data = M_C_data_scaled3,
  family = binomial(link = "logit"),
  chains = 4,
  iter = 4000,
  warmup = 2000,
  cores = 4,
  seed = 123,
  control = list(adapt_delta = 0.9999, max_treedepth = 15)
)

plot(brm_Presence1)




for (var in predictors) {
  formula <- as.formula(paste("Presence_Kōura ~", var, "+ (1 | LID)"))
  model <- glmer(formula, data = M_C_data_scaled3, family = binomial)
  cat("\n", var, "\n")
  print(summary(model))
}

(Presence_Kōura ~ TLI + 
  Elevation_m + 
  Riparian_vegetation + 
  Overhanging_trees +  
  Boulders + 
  Temperature + 
  Submerged_Non_Native + 
  Presence_Morihana + 
  Presence_Eel + 
  Slope_5m + (1 | LID),
data = M_C_data_scaled3, family = binomial)






glm_reduced <- glmer(Presence_Kōura ~ #TLI + 
                     # Elevation_m + 
                      Mean_depth_m +
                      #Riparian_vegetation + 
                      #Overhanging_trees +  
                      Boulders + 
                      #Temperature + 
                      #Submerged_Non_Native + 
                      Presence_Morihana + 
                      Presence_Eel + 
                      Slope_5m + (1 | LID),
                    data = M_C_data_scaled3, family = binomial)
summary(glm_reduced)
vif(glm_reduced)


glm_fixed <- glm(Presence_Kōura ~ TLI + Elevation_m + Boulders + 
                   Presence_Morihana + Presence_Eel + Slope_5m,
                 data = M_C_data_scaled3, family = binomial)
anova(glm_fixed, glm_reduced, test = "Chisq")





glm_step <- step(glm_reduced, direction = "both")
summary(glm_step)
vif(glm_step)

# Make diffrent models
modelall <- as.formula(paste("Presence_Kōura ~", paste(num_cols, collapse = " + "), "+ (1 | LID)"))

full_model3 <- glmer(Presence_Kōura ~ 
                     Slope_5m + TLI + #LID + 
                     Elevation_m + Mean_depth_m + 
                     Riparian_vegetation + Overhanging_trees + Wood_cover + Erosion +
                     Bedrock + Boulders + Cobble + Gravel + Sand + Organic_matter +  Mud + #Presence_rocks + 
                     Temperature + DO_mgl + pH + Specific_conductivity + # DO_percent + Conductivity +
                     Emergent_Native + Emergent_Non_Native + Submerged_Native + Submerged_Non_Native + 
                     Presence_Morihana + Presence_Eel + Presence_Catfish +
                     #Presence_Kōaro + Presence_Common_smelt + Presence_Trout + Presence_Mosquitofish,
                     (1 | LID),
                   data = M_C_data_scaled3, family = binomial(link = "logit"))

options(na.action = "na.fail") # required by dredge
model_set <- dredge(full_model3)

# View the top models
head(model_set)

# Get the best model
best_model <- get.models(model_set, subset = 1)[[1]]
summary(best_model)

# 
model0 <- Presence_Kōura ~ 
  TLI + Elevation_m + Mean_depth_m + 
  Slope_5m + Riparian_vegetation + Overhanging_trees + Wood_cover + Erosion +
  Bedrock + Boulders + Cobble + Gravel + Sand + Organic_matter +  Mud + #Presence_rocks + 
  Temperature + DO_mgl + pH + Specific_conductivity + # DO_percent + Conductivity +
  Emergent_Native + Emergent_Non_Native + Submerged_Native + Submerged_Non_Native + 
  Presence_Morihana + Presence_Eel + Presence_Catfish +
  #Presence_Kōaro + Presence_Common_smelt + Presence_Trout + Presence_Mosquitofish,
  (1 | LID)

model1 <- Presence_Kōura ~
  TLI + Elevation_m + Mean_depth_m + 
  Slope_5m + Riparian_vegetation + Overhanging_trees + Wood_cover + Erosion +
  #Bedrock + Boulders + Cobble + Gravel + Sand + Organic_matter +  Mud + 
  #Temperature + DO_mgl + pH + Specific_conductivity + 
  #Emergent_Native + Emergent_Non_Native + Submerged_Native + Submerged_Non_Native + 
  Presence_Morihana + Presence_Eel + Presence_Catfish +
  (1 | LID)

model2 <- Presence_Kōura ~  
  Slope_5m + Overhanging_trees + 
  Temperature + 
  Presence_Morihana + Presence_Eel + 
  (1 | LID)

model3 <- Presence_Kōura ~  
  Slope_5m + Overhanging_trees + 
  Temperature + Boulders +
  Presence_Morihana + Presence_Eel + 
  (1 | LID)

glmm_Presence0 <- glmer(model0, data = M_C_data_scaled3, family = binomial(link = "logit"))
glmm_Presence1 <- glmer(model1, data = M_C_data_scaled3, family = binomial(link = "logit"))
glmm_Presence2 <- glmer(model2, data = M_C_data_scaled3, family = binomial(link = "logit"))
glmm_Presence3 <- glmer(model3, data = M_C_data_scaled3, family = binomial(link = "logit"))

summary(glmm_Presence0)
summary(glmm_Presence1)
summary(glmm_Presence2)
summary(glmm_Presence3)

# Compare models using likelihood ratio test
anova(glmm_Presence0, glmm_Presence1, glmm_Presence2, glmm_Presence3, test = "Chisq")

fixed_formula1 <- Presence_Kōura ~  
  TLI + #Elevation_m + Mean_depth_m + 
  Slope_5m + Overhanging_trees + #Wood_cover + 
  #Bedrock + Boulders + Cobble + Gravel + Sand + Organic_matter + Mud + 
  #Temperature + #DO_mgl + Specific_conductivity + pH + 
  #Emergent_Native + Emergent_Non_Native + Submerged_Native + Submerged_Non_Native  + 
  Presence_Morihana + Presence_Eel + #Presence_Catfish  + 
  (1 | LID)
