# 3. Data Modeling using GAM
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
packages <- c("patchwork", "gratia", "mgcv", "corrplot","glmmTMB", "performance", "tidyverse", "writexl")

# Load packages if not already installed
lapply(packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, dependencies = TRUE)
  library(pkg, character.only = TRUE)})


# Import the data set ----------------------------------------------------------
Monitoring_CPUE_data <- read_csv("Data_mod/Monitoring_CPUE_data.csv")

# Filter to ignor all sites that have been sampled double.
Modeling_data <- Monitoring_CPUE_data %>%
  filter(Monitoring==0)

names(Modeling_data)
summary(Modeling_data)
#str(Modeling_data)
#colSums(is.na(Modeling_data))

plot(Modeling_data$CPUE_Kōura, Modeling_data$CPUE_Catfish)
plot(Modeling_data$Substrate_index, Modeling_data$CPUE_Catfish)


# Built Functions --------------------------------------------------------------
# remove high VIF for Bionomial data
remove_high_vif_glmBI <- function(data, response, predictors, threshold = 5) {
  # Load brglm2 for bias-reduced logistic regression
  if (!requireNamespace("brglm2", quietly = TRUE)) {
    install.packages("brglm2")}
  library(brglm2)
  removed <- character()
  # Step 1: Remove near-zero variance predictors
  nzv <- caret::nearZeroVar(data[, predictors, drop = FALSE])
  if (length(nzv) > 0) {
    cat("Removing near-zero variance predictors:", predictors[nzv], "\n")
    removed <- c(removed, predictors[nzv])
    predictors <- predictors[-nzv]}
  # Step 2: Iteratively remove high-VIF predictors
  repeat {formula <- as.formula(paste(response, "~", paste(predictors, collapse = " + ")))
    # Use bias-reduced logistic regression to avoid separation issues
    model <- try(glm(formula, data = data, family = binomial(link = "logit"),method = "brglmFit"), silent = TRUE)
    if (inherits(model, "try-error")) break
    vif_data <- performance::check_collinearity(model)
    vif_data <- vif_data[!grepl("\\|", vif_data$Term), ] # Remove random effects if any
    print(vif_data) # Debugging: show VIFs each iteration
    if (all(vif_data$VIF < threshold)) break
    to_remove <- vif_data$Term[which.max(vif_data$VIF)]
    cat("Removing:", to_remove, "\n")
    predictors <- setdiff(predictors, to_remove)
    removed <- c(removed, to_remove)
    if (length(predictors) == 0) stop("All predictors removed.")}
  cat("Removed variables:\n", paste(removed, collapse = ", "), "\n")
  return(model)}

# remove high VIF in glmmTMB function with Tweedie
remove_high_vif_glmmTMB <- function(data, response, predictors, threshold = 5) {
  removed <- character()
  repeat {
    formula <- as.formula(paste(response, "~", paste(predictors, collapse = " + "))) # , "+ (1 | LID)"
    model <- try(glmmTMB(formula, data = data, family = tweedie()), silent = TRUE)
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

# automatic remove the variable with highest p value include LID via build_smooth()
stepwise_gam_vars_with_LID <- function(response, vars, data, family, method = "ML", p_cutoff = 0.1) {
  remaining_vars <- vars
  m <- NULL
  repeat {
    rhs_terms <- vapply(remaining_vars, \(v) build_smooth(v, include_re_for_LID = TRUE), character(1))
    current_formula <- as.formula(paste(response, "~", paste(rhs_terms, collapse = " + ")))
    m <- gam(current_formula, data = data, family = family, method = method, select = TRUE)
    summ <- summary(m)
    
    # parametric p
    if (!is.null(summ$p.table) && nrow(summ$p.table) > 0) {
      pcol <- intersect(c("Pr(>|t|)", "Pr(>|z|)"), colnames(summ$p.table))
      param_p <- summ$p.table[, pcol[1]]
      names(param_p) <- rownames(summ$p.table)
      param_p <- param_p[names(param_p) != "(Intercept)"]
    } else param_p <- numeric(0)
    
    # smooth p
    smooth_p <- if (!is.null(summ$s.table)) {
      p <- summ$s.table[, "p-value"]; names(p) <- rownames(summ$s.table); p
    } else numeric(0)
    
    all_p <- c(param_p, smooth_p)
    all_p_filtered <- all_p[all_p > p_cutoff]
    if (length(all_p_filtered) == 0) {
      message("All terms have p <= ", p_cutoff, ". Stopping.")
      break
    }
    term_to_remove <- names(which.max(all_p_filtered))
    max_pval <- max(all_p_filtered)
    message("Removing term: ", term_to_remove, " with p-value = ", max_pval)
    
    # map smooth label back to variable name
    var_to_remove <- if (grepl("^s\\(", term_to_remove)) sub("^s\\(([^,]+).*\\)$", "\\1", term_to_remove) else term_to_remove
    
    if (var_to_remove %in% remaining_vars) {
      remaining_vars <- setdiff(remaining_vars, var_to_remove)
    } else {
      message("Warning: variable '", var_to_remove, "' not found in remaining_vars!")
      break
    }
    if (length(remaining_vars) == 0) {
      message("No variables left, stopping.")
      break
    }
  }
  m
}

###############################################################################.
# Kōura occupancy (presence/absence) -------------------------------------------
###############################################################################.
## 1) Prep: select vars and coerce types
vars <- c(
  "Presence_Kōura","LID",
  "Slope_5m","Riparian_vegetation","Overhanging_trees","Wood_cover","Substrate_index",
  "Temperature","pH","DO_mgl","DO_percent","Specific_conductivity",
  "Emergent_Native","Submerged_Non_Native","Turf_Native","Submerged_Native","Emergent_Non_Native",
  "Presence_Common_smelt","Presence_Morihana","Presence_Eel","Presence_Catfish")

PData <- Modeling_data |>
  dplyr::select(all_of(vars)) |>
  mutate(LID = factor(LID),
         across(c(Presence_Kōura, Presence_Common_smelt, Presence_Morihana,
                  Presence_Eel, Presence_Catfish), ~ as.numeric(.)))

## 2) VIF pruning with your function — EXCLUDE LID here
response   <- "Presence_Kōura"
pred_fixed <- setdiff(names(PData), c(response, "LID"))  # only fixed effects into VIF
vif_model  <- remove_high_vif_glmBI(PData, response, pred_fixed, threshold = 5)
kept_fixed <- setdiff(names(coef(vif_model)), "(Intercept)")

## 3) Smooth builder that respects variable type (your stepwise uses this)
custom_k <- list(Slope_5m=10, Riparian_vegetation=7, Overhanging_trees=5, Wood_cover=10,
                 Substrate_index=10, Temperature=10, pH=10, DO_mgl=10, Emergent_Native=5, Submerged_Non_Native=9, Turf_Native=6)

is_cont <- function(x) is.numeric(x) && dplyr::n_distinct(x, na.rm = TRUE) >= 5
build_smooth <- function(var, include_re_for_LID = FALSE) {
  if (identical(var, "LID") && include_re_for_LID) return("s(LID, bs='re')")
  x <- PData[[var]]
  if (is_cont(x)) {
    k <- if (!is.null(custom_k[[var]])) paste0(", k=", custom_k[[var]]) else ""
    return(paste0("s(", var, ", bs='ts'", k, ")"))
  }
  var  # binary/low-level → parametric
}

## 4) Full models (± lake RE) using the same fixed set
rhs_fixed <- paste(vapply(kept_fixed, build_smooth, character(1)), collapse = " + ")
form_noLID   <- as.formula(paste0(response, " ~ ", rhs_fixed))
form_withLID <- as.formula(paste0(response, " ~ ", rhs_fixed, " + s(LID, bs='re')"))

mp_full    <- gam(form_noLID,   data = PData, family = binomial(link = "logit"), method = "ML", select = TRUE)
mp_full_re <- gam(form_withLID, data = PData, family = binomial(link = "logit"), method = "ML", select = TRUE)

summary(mp_full)
summary(mp_full_re)
AIC(mp_full, mp_full_re)
anova(mp_full, mp_full_re, test = "Chisq")

## 5) Reduced models via stepwise at two cut-offs
vars_for_stepwise <- c(kept_fixed, "LID")

mp.05 <- stepwise_gam_vars_with_LID(
  response = response, vars = vars_for_stepwise,
  data = PData, family = binomial(link = "logit"),
  method = "ML", p_cutoff = 0.05)

mp.10 <- stepwise_gam_vars_with_LID(
  response = response, vars = vars_for_stepwise,
  data = PData, family = binomial(link = "logit"),
  method = "ML", p_cutoff = 0.10)


## 6) Compare models
summary(mp_full)
summary(mp_full_re)
summary(mp.05)
summary(mp.10)
anova(mp_full, mp_full_re, test = "Chisq")
anova(mp.05, mp_full, test = "Chisq")
anova(mp.05, mp_full_re, test = "Chisq")
anova(mp.10, mp_full, test = "Chisq")
anova(mp.10, mp_full_re, test = "Chisq")
anova(mp.05, mp.10, test = "Chisq")
AIC(mp_full, mp_full_re, mp.05, mp.10)


## 7) Final refit in REML for reporting (best model only)
bestp_model <- formula(mp.10)  
mp_final  <- gam(bestp_model, data = PData, family = binomial(link = "logit"), method = "REML", select = TRUE)
summary(mp_final)
gam.check(mp_final)
concurvity(mp_final)
influence.gam(mp_final)

# Predicted vs observed
plot(predict(mp_final, type = "response"), mp_final$y,xlab = "Predicted presence", ylab = "Observed presence")
abline(0, 1, col = "red")

# helper
logit2prob <- function(x) exp(x) / (1 + exp(x))

make_fish_plot <- function(fish_var, label, model, PData){
  # base grid with other predictors fixed, vary LID
  newdat <- expand.grid(
    LID                 = levels(PData$LID),
    Riparian_vegetation = mean(PData$Riparian_vegetation, na.rm = TRUE),
    Substrate_index     = mean(PData$Substrate_index,     na.rm = TRUE),
    Temperature         = mean(PData$Temperature,         na.rm = TRUE),
    Presence_Morihana   = 0,
    Presence_Eel        = 0,
    Presence_Catfish    = 0
  )
  
  # duplicate rows: absent (0) and present (1) for this fish
  newdat0 <- newdat
  newdat1 <- newdat
  newdat0[[fish_var]] <- 0
  newdat1[[fish_var]] <- 1
  newdat <- rbind(newdat0, newdat1)
  
  # predict excluding RE
  pr <- predict(model, newdata = newdat, type = "link",
                se.fit = TRUE, exclude = "s(LID)")
  
  logit2prob <- function(x) exp(x) / (1 + exp(x))
  newdat$prob <- logit2prob(pr$fit)
  newdat$lwr  <- logit2prob(pr$fit - 1.96*pr$se.fit)
  newdat$upr  <- logit2prob(pr$fit + 1.96*pr$se.fit)
  
  # summarise across lakes (mean over LID)
  plot_df <- newdat |>
    dplyr::group_by(.data[[fish_var]]) |>
    dplyr::summarise(prob = mean(prob),
                     lwr  = mean(lwr),
                     upr  = mean(upr),
                     .groups = "drop")
  
  ggplot(plot_df,
         aes(x = factor(.data[[fish_var]], labels = c("0","1")),
             y = prob)) +
    geom_col(width = 0.5, fill = "grey70") +
    geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2) +
    ylim(0,1) +
    labs(x = label, y = "Predicted kōura presence") 
}

# Now call it:
p_mori <- make_fish_plot("Presence_Morihana", "Presence Morihana", mp_final, PData) +labs(y = NULL)
p_eel  <- make_fish_plot("Presence_Eel", "Presence Eel", mp_final, PData)
p_cat  <- make_fish_plot("Presence_Catfish", "Presence Catfish", mp_final, PData) +labs(y = NULL)

p1 <- draw(mp_final)
p2 <- p_mori/ p_eel/ p_cat

finalp_plot <- p1 | p2
finalp_plot <- finalp_plot + plot_layout(widths = c(2, 1))
finalp_plot

ggsave("Figures/precence_GAM.png", finalp_plot, width = 8, height = 5, dpi = 300)




# Build new data frame 
newdata_param <- expand.grid(
  Presence_Morihana = c(0, 1),
  Presence_Eel = c(0, 1),
  Presence_Catfish = c(0, 1),
  LID = levels(PData$LID),
  Riparian_vegetation = mean(PData$Riparian_vegetation, na.rm = TRUE),
  Substrate_index = mean(PData$Substrate_index, na.rm = TRUE),
  Temperature = mean(PData$Temperature, na.rm = TRUE))

# Predict on link (logit) scale
pred <- predict(mp_final, newdata = newdata_param, type = "link", se.fit = TRUE)

# Convert log-odds to probabilities
logit2prob <- function(x) exp(x) / (1 + exp(x))
newdata_param$predicted_prob <- logit2prob(pred$fit)
newdata_param$lwr <- logit2prob(pred$fit - 1.96 * pred$se.fit)
newdata_param$upr <- logit2prob(pred$fit + 1.96 * pred$se.fit)

# Plots
p1 <- draw(mp_final)
p2 <- ggplot(newdata_param, aes(x = interaction(Presence_Morihana, Presence_Eel, Presence_Catfish), y = predicted_prob)) +
  geom_col(width = 0.5, fill = "grey70") +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2) +
  labs(x = "Presence Morihana and Eel", y = "Predicted Kōura Presence") 



# Get predicted probabilities
pred_probs <- predict(mp_final, type = "response")

# Equal-frequency binning (10 bins)
n_bins <- 10
calib_plot_data <- tibble(
  pred = pred_probs,
  obs = mp_final$y) %>%
  mutate(bin = ntile(pred, n_bins)) %>%  # create equal-frequency bins
  group_by(bin) %>%
  summarise(
    mean_pred = mean(pred),
    obs_rate = mean(obs),
    n = n(),
    .groups = "drop")

# Calibration plot
Calibration_plotp <- ggplot(calib_plot_data, aes(mean_pred, obs_rate)) +
  geom_point(size = 3) +
  geom_line() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "Mean predicted probability",
       y = "Observed proportion") 
Calibration_plotp

ggsave("Figures/Precence_Predictions.png", Calibration_plotp, width = 5, height = 4, dpi = 300)



###############################################################################.
# CPUE kōura (abundances)-------------------------------------------------------
###############################################################################.
## 1) Prep: select vars and coerce types
vars <- c("CPUE_Kōura","LID",
  "Slope_5m","Riparian_vegetation","Overhanging_trees","Wood_cover","Substrate_index",
  "Temperature","pH","DO_mgl","DO_percent","Specific_conductivity",
  "Emergent_Native","Submerged_Non_Native","Turf_Native","Submerged_Native","Emergent_Non_Native",
  "Presence_Common_smelt","Presence_Morihana","Presence_Eel","Presence_Catfish")

CData <- Modeling_data |>
  dplyr::select(all_of(vars)) |>
  mutate(LID = factor(LID),
         across(c(Presence_Common_smelt, Presence_Morihana,
                  Presence_Eel, Presence_Catfish), ~ as.numeric(.)))

## 2) VIF pruning with your function — EXCLUDE LID here
response   <- "CPUE_Kōura"
pred_fixed <- setdiff(names(CData), c(response, "LID"))  # only fixed effects into VIF
vif_model  <- remove_high_vif_glmmTMB(CData, response, pred_fixed, threshold = 5)
kept_fixed <- names(fixef(vif_model)$cond)
kept_fixed <- setdiff(kept_fixed, "(Intercept)")

## 3) Smooth builder that respects variable type (your stepwise uses this)
custom_k <- list(Slope_5m=10, Riparian_vegetation=7, Overhanging_trees=5, Wood_cover=10,
                 Substrate_index=10, Temperature=10, pH=10, DO_mgl=10, Emergent_Native=9, 
                 Submerged_Native = 4, Submerged_Non_Native=9, Turf_Native=6)

is_cont <- function(x) is.numeric(x) && dplyr::n_distinct(x, na.rm = TRUE) >= 5
build_smooth <- function(var, include_re_for_LID = FALSE) {
  if (identical(var, "LID") && include_re_for_LID) return("s(LID, bs='re')")
  x <- CData[[var]]
  if (is_cont(x)) {
    k <- if (!is.null(custom_k[[var]])) paste0(", k=", custom_k[[var]]) else ""
    return(paste0("s(", var, ", bs='ts'", k, ")"))
  }
  var  # binary/low-level → parametric
}

## 4) Full models (± lake RE) using the same fixed set
rhs_fixed <- paste(vapply(kept_fixed, build_smooth, character(1)), collapse = " + ")
form_noLID   <- as.formula(paste0(response, " ~ ", rhs_fixed))
form_withLID <- as.formula(paste0(response, " ~ ", rhs_fixed, " + s(LID, bs='re')"))

mc_full    <- gam(form_noLID,   data = CData, family = tw(link = "log"), method = "ML", select = TRUE)
mc_full_re <- gam(form_withLID, data = CData, family = tw(link = "log"), method = "ML", select = TRUE)

## 5) Reduced models via stepwise at two cut-offs
vars_for_stepwise <- c(kept_fixed, "LID")

mc.05 <- stepwise_gam_vars_with_LID(
  response = response, vars = vars_for_stepwise,
  data = CData, family = tw(link = "log"),
  method = "ML", p_cutoff = 0.05)


## 6) Compare models
summary(mc_full)
summary(mc_full_re)
summary(mc.05)
anova(mc_full, mc_full_re, test = "Chisq")
anova(mc.05, mc_full, test = "Chisq")
anova(mc.05, mc_full_re, test = "Chisq")
AIC(mc_full, mc_full_re, mc.05)


## 7) Final refit in REML for reporting (best model only)
bestc_model <- formula(mc.05)  
mc_final  <- gam(bestc_model, data = CData, family = tw(link = "log"), method = "REML", select = TRUE)
summary(mc_final)
gam.check(mc_final)
concurvity(mc_final)
influence.gam(mc_final)


## 8) Make it visual into plots
new_datac <- expand.grid(
  Presence_Common_smelt = c(0, 1),
  Temperature = mean(CData$Temperature, na.rm = TRUE),
  pH = mean(CData$pH, na.rm = TRUE),
  Emergent_Native = mean(CData$Emergent_Native, na.rm = TRUE),
  Substrate_index = mean(CData$Substrate_index, na.rm = TRUE))

# Predict with SE for CI
pred <- predict(mc_final, newdata = new_datac, type = "link", se.fit = TRUE)
new_datac$predicted_CPUE <- exp(pred$fit)
new_datac$lwr <- exp(pred$fit - 1.96 * pred$se.fit)
new_datac$upr <- exp(pred$fit + 1.96 * pred$se.fit)

p1 <- draw(mc_final)

p2 <- ggplot(new_datac, aes(factor(Presence_Common_smelt), predicted_CPUE)) +
  geom_col(width = 0.5, fill = "grey70") +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2) +
  labs(x = "Presence Common Smelt",y = "Predicted kōura CPUE") 

final_plotc <- p1 | p2
final_plotc <- final_plotc + plot_layout(widths = c(2, 1))
final_plotc

ggsave("Figures/CPUE_GAM.png", final_plotc, width = 8, height = 5, dpi = 300)

# Predicted and observed
pred <- predict(mc_final, type = "response")
obs <- mc_final$y
R2 <- cor(pred, obs)^2
RMSE <- sqrt(mean((pred - obs)^2))

df_pred_obs <- data.frame(ID = 1:length(obs), Predicted = pred, Observed = obs)

# Plot
library(ggrepel)
Predictions_plotc = ggplot(df_pred_obs, aes(x = Predicted, y = Observed)) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  #geom_text_repel(aes(label = ID), size = 3, max.overlaps = 5) +  
  annotate("text", x = max(df_pred_obs$Predicted) * 0.7, y = max(df_pred_obs$Observed) * 0.8,
           label = paste0("Adj. R² = ", round(R2, 2), "\nRMSE = ", round(RMSE, 2)),
           hjust = 0, size = 4, colour = "black") +
  labs(x = "Predicted CPUE", y = "Observed CPUE")

Predictions_plotc

ggsave("Figures/CPUE_Predictions.png", Predictions_plotc, width = 5, height = 4, dpi = 300)





###############################################################################.
# BCUE kōura (abundances)-------------------------------------------------------
###############################################################################.
## 1) Prep: select vars and coerce types
vars <- c("BCUE_Kōura","LID",
          "Slope_5m","Riparian_vegetation","Overhanging_trees","Wood_cover","Substrate_index",
          "Temperature","pH","DO_mgl","DO_percent","Specific_conductivity",
          "Emergent_Native","Submerged_Non_Native","Turf_Native","Submerged_Native","Emergent_Non_Native",
          "Presence_Common_smelt","Presence_Morihana","Presence_Eel","Presence_Catfish")

BData <- Modeling_data |>
  dplyr::select(all_of(vars)) |>
  mutate(LID = factor(LID),
         across(c(Presence_Common_smelt, Presence_Morihana,
                  Presence_Eel, Presence_Catfish), ~ as.numeric(.)))

## 2) VIF pruning with your function — EXCLUDE LID here
response   <- "BCUE_Kōura"
pred_fixed <- setdiff(names(BData), c(response, "LID"))  # only fixed effects into VIF
vif_model  <- remove_high_vif_glmmTMB(BData, response, pred_fixed, threshold = 5)
kept_fixed <- names(fixef(vif_model)$cond)
kept_fixed <- setdiff(kept_fixed, "(Intercept)")

## 3) Smooth builder that respects variable type (your stepwise uses this)
custom_k <- list(Slope_5m=10, Riparian_vegetation=7, Overhanging_trees=5, Wood_cover=10,
                 Substrate_index=10, Temperature=10, pH=10, DO_mgl=10, Emergent_Native=9, 
                 Submerged_Native = 4, Submerged_Non_Native=9, Turf_Native=6)

is_cont <- function(x) is.numeric(x) && dplyr::n_distinct(x, na.rm = TRUE) >= 5
build_smooth <- function(var, include_re_for_LID = FALSE) {
  if (identical(var, "LID") && include_re_for_LID) return("s(LID, bs='re')")
  x <- BData[[var]]
  if (is_cont(x)) {
    k <- if (!is.null(custom_k[[var]])) paste0(", k=", custom_k[[var]]) else ""
    return(paste0("s(", var, ", bs='ts'", k, ")"))
  }
  var  # binary/low-level → parametric
}

## 4) Full models (± lake RE) using the same fixed set
rhs_fixed <- paste(vapply(kept_fixed, build_smooth, character(1)), collapse = " + ")
form_noLID   <- as.formula(paste0(response, " ~ ", rhs_fixed))
form_withLID <- as.formula(paste0(response, " ~ ", rhs_fixed, " + s(LID, bs='re')"))

mb_full    <- gam(form_noLID,   data = BData, family = tw(link = "log"), method = "ML", select = TRUE)
mb_full_re <- gam(form_withLID, data = BData, family = tw(link = "log"), method = "ML", select = TRUE)

## 5) Reduced models via stepwise at two cut-offs
vars_for_stepwise <- c(kept_fixed, "LID")

mb.05 <- stepwise_gam_vars_with_LID(
  response = response, vars = vars_for_stepwise,
  data = BData, family = tw(link = "log"),
  method = "ML", p_cutoff = 0.05)


## 6) Compare models
summary(mb_full)
summary(mb_full_re)
summary(mb.05)
anova(mb_full, mb_full_re, test = "Chisq")
anova(mb.05, mb_full, test = "Chisq")
anova(mb.05, mb_full_re, test = "Chisq")
AIC(mb_full, mb_full_re, mb.05)


## 7) Final refit in REML for reporting (best model only)
bestb_model <- formula(mb.05)  
mb_final  <- gam(bestb_model, data = BData, family = tw(link = "log"), method = "REML", select = TRUE)
summary(mb_final)
gam.check(mb_final)
concurvity(mb_final)
influence.gam(mb_final)


## 8) Make it visual into plots
newdata <- expand.grid(
  Presence_Morihana = c(0, 1),
  Presence_Eel = c(0, 1),
  Overhanging_trees = mean(BData$Overhanging_trees, na.rm = TRUE),
  pH = mean(BData$pH, na.rm = TRUE))

# Predict on link (logit) scale
pred <- predict(mb_final, newdata = newdata, type = "link", se.fit = TRUE)

# Convert log-odds to probabilities
logit2prob <- function(x) exp(x) / (1 + exp(x))
newdata$predicted_prob <- logit2prob(pred$fit)
newdata$lwr <- logit2prob(pred$fit - 1.96 * pred$se.fit)
newdata$upr <- logit2prob(pred$fit + 1.96 * pred$se.fit)

# Plots
p1 <- draw(mb_final, ncol = 1)
p2 <- ggplot(newdata, aes(x = interaction(Presence_Morihana, Presence_Eel), y = predicted_prob)) +
  geom_col(width = 0.5, fill = "grey70") +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2) +
  labs(x = "Presence Morihana and Eel", y = "Predicted kōura BCUE") 

final_plot <- p1 | p2
final_plot <- final_plot + plot_layout(widths = c(1.5, 2))
final_plot

ggsave("Figures/BCUE_GAM.png", final_plot, width = 8, height = 5, dpi = 300)




# Predicted and observed
pred <- predict(mb_final, type = "response")
obs <- mb_final$y
R2 <- cor(pred, obs)^2
RMSE <- sqrt(mean((pred - obs)^2))

df_pred_obs <- data.frame(ID = 1:length(obs), Predicted = pred, Observed = obs)

# Plot
library(ggrepel)
Predictions_plot=ggplot(df_pred_obs, aes(x = Predicted, y = Observed)) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  #geom_text_repel(aes(label = ID), size = 3, max.overlaps = 5) +  
  annotate("text", x = max(df_pred_obs$Predicted) * 0.7, y = max(df_pred_obs$Observed) * 0.8,
           label = paste0("Adj. R² = ", round(R2, 2), "\nRMSE = ", round(RMSE, 2)),
           hjust = 0, size = 4, colour = "black") +
  labs(x = "Predicted BCUE", y = "Observed BCUE")

Predictions_plot

ggsave("Figures/BCUE_Predictions.png", Predictions_plot, width = 5, height = 4, dpi = 300)







#
# OLD####
# 1. Preprocessing 
vars <- c("CPUE_Kōura", "Lake", "LID",
          "Slope_5m", "Riparian_vegetation", "Overhanging_trees",
          "Wood_cover", "Substrate_index", 
         # "Bedrock", "Boulders", "Cobble", "Gravel", "Sand", "Organic_matter",
          "Temperature", "DO_percent", "pH", "Specific_conductivity", "DO_mgl",
          "Emergent_Native", "Submerged_Native", "Submerged_Non_Native", "Emergent_Non_Native", "Turf_Native",
          "Presence_Kōaro", "Presence_Common_smelt", "Presence_Morihana",#"Presence_Eel","Presence_Catfish", 
         "Predator_Fish_Presence")

Modeling_data1 <- Modeling_data %>%
  dplyr::select(all_of(vars)) %>%
  na.omit() 

# 2. Correlation check 
cor_matrix <- cor(dplyr::select(Modeling_data1, where(is.numeric)), use = "complete.obs")
corrplot(cor_matrix)

# 3. VIF reduction 
response <- "CPUE_Kōura"
predictors <- setdiff(names(Modeling_data1), response)

model <- remove_high_vif_glmmTMB(Modeling_data1, response, predictors)
#summary(model)
performance::check_collinearity(model)

# 4. Post-VIF correlation filtering 
clean_vars <- c(setdiff(names(fixef(model)$cond), "(Intercept)"), response)
cor_filtered <- Modeling_data1 %>% dplyr::select(all_of(clean_vars))
corrplot(cor(cor_filtered, use = "complete.obs"))
clean_vars <- setdiff(clean_vars, response)


# preform PCA to explore vomon variables.
pca <- prcomp(Modeling_data1[, clean_vars], scale. = TRUE)
library(factoextra)
fviz_eig(pca)
fviz_pca_biplot(pca)
fviz_pca_ind(pca)
fviz_pca_var(pca)

# spearman_test
#spearman_test <- sapply(clean_vars, function(var) {cor.test(Modeling_data1[[var]], Modeling_data1$CPUE_Kōura, method = "spearman", exact = FALSE)$p.value})
#significant <- names(spearman_test)[spearman_test < 0.05]
#significant

# 5. Build GAM formulas 
# Custom smooth parameters
sapply(Modeling_data1[clean_vars], function(x) length(unique(x)))
custom_k <- list(Overhanging_trees=5, Riparian_vegetation=7,
                 Bedrock=3, Boulders=5, Cobble=4, Gravel=6,
                 Emergent_Native=5, Submerged_Native=4, Submerged_Non_Native=9, Turf_Native=6)

build_smooth <- function(var) {
  n <- length(unique(Modeling_data1[[var]]))
  if (n >= 5) {
    k <- if (!is.null(custom_k[[var]])) paste0(", k=", custom_k[[var]]) else ""
    paste0("s(", var, ", bs='ts'", k, ")")
  } else var}

# 1) Which predictors are numeric with enough unique values to smooth?
is_num <- sapply(Modeling_data1[clean_vars], function(x) is.numeric(x) && length(unique(x)) >= 5)
num_vars <- names(is_num[is_num])
non_num  <- setdiff(clean_vars, num_vars)  # binary 0/1 or small-level vars → parametric

# 2) Your build_smooth() already defined; use it across num_vars
smooth_terms <- sapply(num_vars, build_smooth, USE.NAMES = FALSE)

# 3) Combine smooths + parametric terms
rhs_terms <- c(smooth_terms, non_num)

# 4) Build a proper formula (response has a macron → backticks)
form_str <- paste0("`", response, "` ~ ", paste(rhs_terms, collapse = " + "))
full_formula <- as.formula(form_str)

# 5) Fit shrinkage GAM
m_full <- gam(full_formula,
              data   = Modeling_data1,
              family = tw(link = "log"),
              method = "ML",
              select = TRUE)

summary(m_full)
gam.check(m_full)
smooth_df <- as.data.frame(summary(m_full)$s.table) %>%
  tibble::rownames_to_column("term")
smooth_df %>%
  filter(edf > 0.1)

# try with lake effect
Modeling_data1 <- Modeling_data1 |>
  dplyr::mutate(LID = factor(LID))

# if you build terms programmatically, just append the RE smooth:
rhs_terms <- c(smooth_terms, non_num, "s(LID, bs='re')")

form_str <- paste0("`", response, "` ~ ", paste(rhs_terms, collapse = " + "))
full_formula <- as.formula(form_str)

m_full_re <- gam(full_formula,
                 data   = Modeling_data1,
                 family = tw(link="log"),
                 method = "ML",
                 select = TRUE)
summary(m_full_re)
gam.check(m_full_re)

# use stepwise to find best model
m.05 <- stepwise_gam_vars(
  response = response,
  vars = clean_vars,
  data = Modeling_data1,
  family = tw(link = "log"),
  method = "ML",
  p_cutoff = 0.05)


summary(m_full)
summary(m_full_re)
summary(m.05)

anova(m.05, m_full, test = "Chisq")
anova(m.05, m_full_re, test = "Chisq")
anova(m_full, m_full_re, test = "Chisq")
AIC(m_full, m_full_re, m.05)

formula_m.05 <- formula(m.05)

#formula_final <- CPUE_Kōura ~ Lake + s(Substrate_index, bs = "ts") +
#  s(Temperature, bs = "ts") + s(pH, bs = "ts")

formula_m.11 <- CPUE_Kōura ~ s(Substrate_index, bs = "ts") +
  s(Temperature, bs = "ts") + s(pH, bs = "ts") + Presence_Common_smelt

m.11 <- gam(formula_m.11, family = tw(link = "log"), method = "ML", data = Modeling_data1)


formula_m.12 <- CPUE_Kōura ~ s(Substrate_index, bs = "ts") + s(Temperature,bs = "ts") + 
  s(pH, bs = "ts") + s(Emergent_Native, bs = "ts",k = 5)

m.12 <- gam(formula_m.12, family = tw(link = "log"), method = "ML", data = Modeling_data1)

summary(m.1)
summary(m.11)
summary(m.12)
anova(m.11, m.1, test = "Chisq")
anova(m.12, m.1, test = "Chisq")
anova(m.11, m.12, test = "Chisq")
AIC(m.1, m.11, m.12)

m05<- gam(formula_m.05, family = tw(link = "log"), method = "REML", data = Modeling_data1)
summary(m05)
gam.check(m05)
concurvity(m05)
influence.gam(m05)

# Visual diagnostics
plot(m05, pages = 1, shade = TRUE)
draw(m05)



new_data <- expand.grid(
  Presence_Common_smelt = c(0, 1),
  Temperature = mean(Modeling_data1$Temperature, na.rm = TRUE),
  pH = mean(Modeling_data1$pH, na.rm = TRUE),
  Emergent_Native = mean(Modeling_data1$Emergent_Native, na.rm = TRUE),
  Substrate_index = mean(Modeling_data1$Substrate_index, na.rm = TRUE))

# Predict with SE for CI
pred <- predict(m05, newdata = new_data, type = "link", se.fit = TRUE)
new_data$predicted_CPUE <- exp(pred$fit)
new_data$lwr <- exp(pred$fit - 1.96 * pred$se.fit)
new_data$upr <- exp(pred$fit + 1.96 * pred$se.fit)
 
p1 <- draw(m05)

p2 <- ggplot(new_data, aes(factor(Presence_Common_smelt), predicted_CPUE)) +
  geom_col(width = 0.5, fill = "grey70") +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2) +
  labs(x = "Presence of Common Smelt",y = "Predicted CPUE (kōura)") 

final_plot <- p1 | p2
final_plot <- final_plot + plot_layout(widths = c(2, 1))
final_plot

ggsave("Figures/CPUE_GAM.png", final_plot, width = 8, height = 5, dpi = 300)

# Predicted and observed
pred <- predict(m05, type = "response")
obs <- m05$y
R2 <- cor(pred, obs)^2
RMSE <- sqrt(mean((pred - obs)^2))

df_pred_obs <- data.frame(ID = 1:length(obs), Predicted = pred, Observed = obs)

# Plot
library(ggrepel)
Predictions_plot=ggplot(df_pred_obs, aes(x = Predicted, y = Observed)) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  #geom_text_repel(aes(label = ID), size = 3, max.overlaps = 5) +  
  annotate("text", x = max(df_pred_obs$Predicted) * 0.7, y = max(df_pred_obs$Observed) * 0.8,
           label = paste0("Adj. R² = ", round(R2, 2), "\nRMSE = ", round(RMSE, 2)),
           hjust = 0, size = 4, colour = "black") +
  labs(x = "Predicted CPUE", y = "Observed CPUE")

ggsave("Figures/CPUE_Predictions.png", Predictions_plot, width = 5, height = 4, dpi = 300)


Modeling_data %>%
  #filter(pH >= 6) %>%
  ggplot(aes(pH, CPUE_Kōura)) +
  geom_point(aes(col = Lake)) +
  geom_smooth()
Modeling_data %>%
  #filter(pH >= 6) %>%
  ggplot(aes(Temperature, CPUE_Kōura)) +
  geom_point(aes(col = Lake)) +
  geom_smooth()
Modeling_data %>%
  #filter(pH >= 6) %>%
  ggplot(aes(Substrate_index, CPUE_Kōura)) +
  geom_point(aes(col = Lake)) +
  geom_smooth()

boxplot(Modeling_data1$CPUE_Kōura~Modeling_data1$Lake)

ggplot(Modeling_data1, aes(y=CPUE_Kōura,x= Slope_5m))+
  geom_point(aes(col=Lake))+
  geom_smooth()

# Testing lake effect
library(lme4)
model_mixed <- lmer(
  CPUE_Kōura ~ Boulders + Cobble + Temperature + Emergent_Native + 
    Presence_Common_smelt + (1 | LID),
  data = Modeling_data1,
  REML = FALSE)

model_mixed_null <- lm(
  CPUE_Kōura ~ Boulders + Cobble + Temperature + Emergent_Native + 
    Presence_Common_smelt,
  data = Modeling_data1,
  REML = FALSE)

lapply(list(model_mixed_null, model_mixed), summary)
AIC(model_mixed_null, model_mixed)




normality_results <- lapply(clean_vars, function(var) {
  # Fit a simple model for each variable with Lake as factor
  model <- lm(as.formula(paste(var, "~ Lake")), data = Modeling_data1)
  
  # Extract residuals
  res <- residuals(model)
  
  # Shapiro-Wilk normality test
  shapiro <- shapiro.test(res)
  
  data.frame(
    Variable = var,
    W = shapiro$statistic,
    p_value = shapiro$p.value,
    Normal = ifelse(shapiro$p.value > 0.05, TRUE, FALSE)
  )
})

normality_df <- bind_rows(normality_results)
#normality_df

results <- lapply(clean_vars, function(var) {
  model <- lm(as.formula(paste(var, "~ Lake")), data = Modeling_data1)
  anova_res <- anova(model)
  data.frame(
    Variable = var,
    F_value = anova_res$`F value`[1],
    p_value = anova_res$`Pr(>F)`[1]
  )
})
results_df <- do.call(rbind, results)
#results_df

kruskal_results <- lapply(clean_vars, function(var) {
  kruskal_res <- kruskal.test(as.formula(paste(var, "~ Lake")), data = Modeling_data1)
  data.frame(
    Variable = var,
    chi_sq = kruskal_res$statistic,
    p_value = kruskal_res$p.value
  )
})
kruskal_df <- do.call(rbind, kruskal_results)
#kruskal_df

combined_df <- normality_df %>%
  left_join(results_df, by = "Variable") %>%
  left_join(kruskal_df, by = "Variable") %>%
  mutate(Test_used = ifelse(Normal, "ANOVA", "Kruskal-Wallis"), p_value_used = ifelse(Normal, p_value.x, p_value.y)) %>%
  select(Variable, Normal, Test_used, p_value_used) %>%
  arrange(p_value_used)
combined_df

pairwise.wilcox.test(Modeling_data1$Temperature, Modeling_data1$Lake, p.adjust.method = "BH", exact = FALSE)

lake_levels <- levels(Chemical_data$Lake)
pairwise_comparisons <- combn(lake_levels, 2, simplify = FALSE)

ggboxplot(Chemical_data, x = "Lake", y = "Values", color = "Lake", facet.by = "Variable") +
  stat_compare_means(comparisons = pairwise_comparisons, method = "wilcox.test", 
                     label = "p.signif", tip.length = 0.01)
#

# 1. Preprocessing 
vars <- c("Presence_Kōura", "LID",
          "Slope_5m", "Riparian_vegetation", "Overhanging_trees","Wood_cover", "Substrate_index",
          # "Bedrock","Boulders", "Cobble", "Gravel", "Sand", "Organic_matter",
          "Temperature", "DO_percent", "pH", "Specific_conductivity", "DO_mgl",
          "Emergent_Native","Submerged_Non_Native", "Turf_Native", "Submerged_Native","Emergent_Non_Native","Predator_Fish_Presence",
          "Presence_Kōaro", "Presence_Common_smelt","Presence_Morihana","Presence_Eel","Presence_Catfish")

Modeling_data2 <- Modeling_data %>%
  dplyr::select(all_of(vars))

library(caret)
nearZeroVar(Modeling_data2, saveMetrics = TRUE)

# 2. Correlation check 
cor_matrix <- cor(dplyr::select(Modeling_data2, where(is.numeric)), use = "complete.obs")
corrplot(cor_matrix)

# 3. VIF reduction 
response <- "Presence_Kōura"
predictors <- setdiff(names(Modeling_data2), response)

model <- remove_high_vif_glmBI(Modeling_data2, response, predictors)
#summary(model)
performance::check_collinearity(model)

# 4. Post-VIF correlation filtering 
coef_names <- names(coef(model))
clean_vars <- setdiff(coef_names, "(Intercept)")
clean_vars <- c(clean_vars, response)
cor_filtered <- Modeling_data2 %>% dplyr::select(all_of(clean_vars))
corrplot(cor(cor_filtered, use = "complete.obs"))
clean_vars <- setdiff(clean_vars, response)

# Mann–Whitney test
MW_test <- sapply(clean_vars, function(var) {
  wilcox.test(Modeling_data2[[var]] ~ Modeling_data2$Presence_Kōura, exact = FALSE)$p.value})
significant <- names(MW_test)[MW_test < 0.05]
significant



# 5. Build GAM formulas 
# test full models with and without lake as ramdom effect.
# 0) make sure types are right
bin_vars <- c("Presence_Common_smelt","Presence_Morihana","Presence_Eel","Presence_Catfish")
Modeling_data2 <- Modeling_data2 |>
  mutate(
    LID = factor(LID),
    across(all_of(bin_vars), ~ as.numeric(.)))

# 1) custom k (optional)
custom_k <- list(
  Slope_5m=10, Riparian_vegetation=7, Overhanging_trees=5, Wood_cover=10,
  Substrate_index=10, Temperature=10, pH=10, DO_mgl=10,
  Emergent_Native=5, Submerged_Non_Native=9, Turf_Native=6
)

# 2) helper to decide term type
is_continuous <- function(x) is.numeric(x) && dplyr::n_distinct(x, na.rm=TRUE) >= 5

build_smooth <- function(var, include_re = TRUE) {
  x <- Modeling_data2[[var]]
  # random effect for LID only if requested
  if (var == "LID" && include_re) return("s(LID, bs='re')")
  # continuous → smooth
  if (is.numeric(x) && dplyr::n_distinct(x, na.rm=TRUE) >= 5) {
    k <- if (!is.null(custom_k[[var]])) paste0(", k=", custom_k[[var]]) else ""
    return(paste0("s(", var, ", bs='ts'", k, ")"))
  }
  # otherwise parametric
  return(var)
}

# 3) split variables by type 
cont_vars <- clean_vars[vapply(clean_vars, \(v) is_continuous(Modeling_data2[[v]]), logical(1))] 
param_vars <- setdiff(clean_vars, c(cont_vars, "LID"))

# without lake effect
rhs_terms_noLID <- c(
  vapply(setdiff(cont_vars, "LID"), build_smooth, character(1), include_re = FALSE),
  param_vars)
# with lake effect
rhs_terms_withLID <- c(
  vapply(c(cont_vars, "LID"), build_smooth, character(1), include_re = TRUE),
  param_vars)

# 4) Assemble formulas
response <- "Presence_Kōura"

form_noLID   <- as.formula(paste0("`", response, "` ~ ", paste(rhs_terms_noLID, collapse = " + ")))
form_withLID <- as.formula(paste0("`", response, "` ~ ", paste(rhs_terms_withLID, collapse = " + ")))

# 5) Fit shrinkage GAM
m_full <- gam(form_noLID,
              data   = Modeling_data2,
              family = binomial,
              method = "ML",
              select = TRUE)
summary(m_full)
gam.check(m_full)

m_full_re <- gam(form_withLID,
                 data   = Modeling_data2,
                 family = binomial,
                 method = "ML",
                 select = TRUE)
summary(m_full_re)
gam.check(m_full_re)

# 7) Stepwise model selection
mp.05 <- stepwise_gam_vars(
  response = response,
  vars     = clean_vars,        
  data     = Modeling_data2,
  family   = binomial(link = "logit"),
  method   = "ML",              
  p_cutoff = 0.05)

mp.1 <- stepwise_gam_vars(
  response = response,
  vars     = clean_vars,        
  data     = Modeling_data2,
  family   = binomial(link = "logit"),
  method   = "ML",              
  p_cutoff = 0.1)

summary(m_full)
summary(m_full_re)
summary(mp.05)
summary(mp.1)
anova(m_full, m_full_re, test = "Chisq")
anova(mp.05, m_full, test = "Chisq")
anova(mp.1, m_full, test = "Chisq")
anova(mp.05, m_full_re, test = "Chisq")
anova(mp.1, m_full_re, test = "Chisq")
anova(mp.05, mp.1, test = "Chisq")
AIC(m_full,m_full_re ,mp.05, mp.1)



# Make best model into plots 
formula_mp.1 <- formula(mp.1)

mp1 <- gam(formula_mp.1, family = binomial, method = "REML", data = Modeling_data2)
summary(mp1)
gam.check(mp1)
concurvity(mp1)
influence.gam(mp1)

# Predicted vs observed
plot(predict(mp1, type = "response"), mp1$y,xlab = "Predicted presence", ylab = "Observed presence")
abline(0, 1, col = "red")

# Build new data frame 
newdata_param <- expand.grid(
  Presence_Morihana = c(0, 1),
  Presence_Eel = c(0, 1),
  Presence_Catfish = c(0, 1),   
  LID = levels(Modeling_data2$LID),
  Riparian_vegetation = mean(Modeling_data2$Riparian_vegetation, na.rm = TRUE),
  Substrate_index = mean(Modeling_data2$Substrate_index, na.rm = TRUE),
  Temperature = mean(Modeling_data2$Temperature, na.rm = TRUE))

# Predict on link (logit) scale
pred <- predict(mp1, newdata = newdata_param, type = "link", se.fit = TRUE)

# Convert log-odds to probabilities
logit2prob <- function(x) exp(x) / (1 + exp(x))
newdata_param$predicted_prob <- logit2prob(pred$fit)
newdata_param$lwr <- logit2prob(pred$fit - 1.96 * pred$se.fit)
newdata_param$upr <- logit2prob(pred$fit + 1.96 * pred$se.fit)

# Plots
p1 <- draw(mp1)
p2 <- ggplot(newdata_param, aes(x = interaction(Presence_Morihana, Presence_Eel), y = predicted_prob)) +
  geom_col(width = 0.5, fill = "grey70") +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2) +
  labs(x = "Presence Morihana and Eel", y = "Predicted Kōura Presence") 

final_plot <- p1 | p2
final_plot <- final_plot + plot_layout(widths = c(2, 1))
final_plot




























# Custom smooth parameters
custom_k <- list(LID=4, Mean_depth_m=4,
  Overhanging_trees=5, Riparian_vegetation=7,Overhanging_trees=5,
  Bedrock=3, Boulders=5, Cobble=4, Gravel=6,
  Emergent_Native=5, Submerged_Native=4, Submerged_Non_Native=9, Turf_Native=6)

# Build smooitng funcion
build_smooth <- function(var) {
  n <- length(unique(Modeling_data2[[var]]))
  if (n >= 5) {
    k <- if (!is.null(custom_k[[var]])) paste0(", k=", custom_k[[var]]) else ""
    paste0("s(", var, ", bs='ts'", k, ")")
  } else var}

# 1) Which predictors are numeric with enough unique values to smooth?
is_num <- sapply(Modeling_data2[clean_vars], function(x) is.numeric(x) && length(unique(x)) >= 5)
num_vars <- names(is_num[is_num])
non_num  <- setdiff(clean_vars, num_vars)  # binary 0/1 or small-level vars → parametric

# 2) Your build_smooth() already defined; use it across num_vars
smooth_terms <- sapply(num_vars, build_smooth, USE.NAMES = FALSE)

# 3) Combine smooths + parametric terms
rhs_terms <- c(smooth_terms, non_num)

# 4) Build a proper formula (response has a macron → backticks)
form_str <- paste0("`", response, "` ~ ", paste(rhs_terms, collapse = " + "))
full_formula <- as.formula(form_str)

# 5) Fit shrinkage GAM
m_full <- gam(full_formula,
              data   = Modeling_data2,
              family = binomial,
              method = "ML",
              select = TRUE)

summary(m_full)
gam.check(m_full)

# try with lake effect
Modeling_data2 <- Modeling_data2 |>
  dplyr::mutate(LID = factor(LID))

# if you build terms programmatically, just append the RE smooth:
rhs_terms <- c(smooth_terms, non_num, "s(LID, bs='re')")

form_str <- paste0("`", response, "` ~ ", paste(rhs_terms, collapse = " + "))
full_formula <- as.formula(form_str)

m_full_re <- gam(full_formula,
                 data   = Modeling_data2,
                 family = binomial,
                 method = "ML",
                 select = TRUE)
summary(m_full_re)
gam.check(m_full_re)

mp.05 <- stepwise_gam_vars(
  response = response,
  vars = clean_vars,
  data = Modeling_data2,
  family = binomial,
  method = "ML",
  p_cutoff = 0.05)

mp.1 <- stepwise_gam_vars(
  response = response,
  vars = clean_vars,
  data = Modeling_data2,
  family = binomial,
  method = "ML",
  p_cutoff = 0.1)

summary(m_full)
summary(m_full_re)
summary(mp.05)
summary(mp.1)
anova(m_full, m_full_re, test = "Chisq")
anova(mp.05, m_full, test = "Chisq")
anova(mp.1, m_full, test = "Chisq")
anova(mp.05, m_full_re, test = "Chisq")
anova(mp.1, m_full_re, test = "Chisq")
anova(mp.05, mp.1, test = "Chisq")
AIC(m_full,m_full_re ,mp.05, mp.1)

formula_mp.1 <- formula(mp.1)


#formula_final2 <- response ~ Lake+ s(Overhanging_trees, bs = "ts", k = 5) + 
#  s(Temperature, bs = "ts") + Presence_Morihana + Presence_Eel

formula_mp.11 <- Presence_Kōura ~ s(Substrate_index, bs = "ts") +  s(Temperature, bs = "ts") +
  s(Riparian_vegetation, bs = "ts", k = 7) 
mp.11 <- gam(formula_mp.11, family = binomial, method = "ML", data = Modeling_data2)

formula_mp.12 <- Presence_Kōura ~ s(Substrate_index, bs = "ts") + s(Temperature,bs = "ts") + 
  Presence_Morihana 
mp.12 <- gam(formula_mp.12, family = binomial, method = "ML", data = Modeling_data2)

formula_mp.13 <- Presence_Kōura ~ s(Substrate_index, bs = "ts") + s(Temperature,bs = "ts") + 
  Presence_Eel 
mp.13 <- gam(formula_mp.13, family = binomial, method = "ML", data = Modeling_data2)


formula_mp.14 <- Presence_Kōura ~ s(Substrate_index, bs = "ts") +  s(Temperature, bs = "ts") +
  s(Riparian_vegetation, bs = "ts", k = 7) + Presence_Morihana
mp.14 <- gam(formula_mp.14, family = binomial, method = "ML", data = Modeling_data2)

formula_mp.15 <- Presence_Kōura ~ s(Substrate_index, bs = "ts") + s(Temperature,bs = "ts") + 
  s(Riparian_vegetation, bs = "ts", k = 7) + Presence_Eel 
mp.15 <- gam(formula_mp.15, family = binomial, method = "ML", data = Modeling_data2)

formula_mp.16 <- Presence_Kōura ~ s(Substrate_index, bs = "ts") + s(Temperature,bs = "ts") + 
  Presence_Morihana + Presence_Eel 
mp.16 <- gam(formula_mp.16, family = binomial, method = "ML", data = Modeling_data2)



summary(mp.1)
summary(mp.11)
summary(mp.12)
summary(mp.13)
summary(mp.14)
summary(mp.15)
summary(mp.16)
anova(mp.11, mp.1, test = "Chisq")
anova(mp.12, mp.1, test = "Chisq")
anova(mp.13, mp.1, test = "Chisq")
anova(mp.14, mp.1, test = "Chisq")
anova(mp.15, mp.1, test = "Chisq")
anova(mp.16, mp.1, test = "Chisq")
#anova(mp.11, mp.12, test = "Chisq")
#anova(mp.11, mp.13, test = "Chisq")
#anova(mp.12, mp.13, test = "Chisq")
AIC(mp.1, mp.11, mp.12, mp.13, mp.14, mp.15, mp.16)


mp1 <- gam(formula_mp.1, family = binomial, method = "REML", data = Modeling_data2)
summary(mp1)
gam.check(mp1)
concurvity(mp1)
influence.gam(mp1)

# Predicted vs observed
plot(predict(mp1, type = "response"), mp1$y,xlab = "Predicted presence", ylab = "Observed presence")
abline(0, 1, col = "red")

# Build new data frame 
newdata_param <- expand.grid(
  Presence_Morihana = c(0, 1),
  Presence_Eel = c(0, 1),
  LID = mean(Modeling_data2$LID, na.rm = TRUE),
  Riparian_vegetation = mean(Modeling_data2$Riparian_vegetation, na.rm = TRUE),
  Substrate_index = mean(Modeling_data2$Substrate_index, na.rm = TRUE),
  Temperature = mean(Modeling_data2$Temperature, na.rm = TRUE))

# Predict on link (logit) scale
pred <- predict(mp1, newdata = newdata_param, type = "link", se.fit = TRUE)

# Convert log-odds to probabilities
logit2prob <- function(x) exp(x) / (1 + exp(x))
newdata_param$predicted_prob <- logit2prob(pred$fit)
newdata_param$lwr <- logit2prob(pred$fit - 1.96 * pred$se.fit)
newdata_param$upr <- logit2prob(pred$fit + 1.96 * pred$se.fit)

# Plots
p1 <- draw(mp1)
p2 <- ggplot(newdata_param, aes(x = interaction(Presence_Morihana, Presence_Eel), y = predicted_prob)) +
  geom_col(width = 0.5, fill = "grey70") +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2) +
  labs(x = "Presence Morihana and Eel", y = "Predicted Kōura Presence") 

final_plot <- p1 | p2
final_plot <- final_plot + plot_layout(widths = c(2, 1))
final_plot

ggsave("Figures/precence_GAM.png", final_plot, width = 8, height = 5, dpi = 300)


# Get predicted probabilities
pred_probs <- predict(mp1, type = "response")

# Equal-frequency binning (10 bins)
n_bins <- 10
calib_plot_data <- tibble(
  pred = pred_probs,
  obs = mp1$y) %>%
  mutate(bin = ntile(pred, n_bins)) %>%  # create equal-frequency bins
  group_by(bin) %>%
  summarise(
    mean_pred = mean(pred),
    obs_rate = mean(obs),
    n = n(),
    .groups = "drop")

# Calibration plot
Calibration_plot <- ggplot(calib_plot_data, aes(mean_pred, obs_rate)) +
  geom_point(size = 3) +
  geom_line() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "Mean predicted probability",
    y = "Observed proportion") 
Calibration_plot

ggsave("Figures/Precence_Predictions.png", Calibration_plot, width = 5, height = 4, dpi = 300)















