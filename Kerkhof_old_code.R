# Koura plot without 0 sites
```{r fig.width=6, fig.height=6, dpi=100}
plot_koura_stats <- function(data, y_var, y_label, 
                             test_type = c("wilcox", "fisher"),
                             show_x_title = FALSE) {
  test_type <- match.arg(test_type)
  lakes <- levels(data$Lake)
  xlab_text <- if (show_x_title) "Lake" else NULL
  
  if (test_type == "wilcox") {
    data_nz <- data %>% dplyr::filter(.data[[y_var]] > 0)
    keep_lakes <- data_nz %>%
      dplyr::group_by(Lake) %>%
      dplyr::summarise(n = sum(!is.na(.data[[y_var]])), .groups = "drop") %>%
      dplyr::filter(n > 0) %>% dplyr::pull(Lake) %>% as.character()
    level_order <- as.character(lakes)[as.character(lakes) %in% keep_lakes]
    data_nz <- data_nz %>% dplyr::mutate(Lake = factor(Lake, levels = level_order))
    
    stats <- data_nz %>%
      rstatix::wilcox_test(as.formula(paste0(y_var, " ~ Lake")), p.adjust.method = "BH") %>%
      { multcompView::multcompLetters(setNames(.$p.adj, paste(.$group1, .$group2, sep = "-")))$Letters } %>%
      tibble::enframe("Lake", "Letter")
    
    ggplot2::ggplot(data_nz, ggplot2::aes(Lake, .data[[y_var]])) +
      ggplot2::geom_boxplot(fill = "grey70") +
      ggplot2::geom_text(data = stats, ggplot2::aes(x = Lake, y = Inf, label = Letter), vjust = 1.1, inherit.aes = FALSE) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.15))) +
      ggplot2::labs(y = y_label, x = xlab_text)
    
  } else {
    presence_summary <- data %>%
      dplyr::group_by(Lake) %>%
      dplyr::summarise(Presence_Rate = mean(.data[[y_var]], na.rm = TRUE), .groups = "drop") %>%
      dplyr::mutate(Lake = factor(Lake, levels = lakes))
    
    lake_pairs <- combn(lakes, 2, simplify = FALSE)
    fisher_pair_test <- function(pair, df) {
      sub <- df %>% dplyr::filter(Lake %in% pair)
      tab <- table(sub[[y_var]], sub$Lake)
      pval <- fisher.test(tab)$p.value
      tibble::tibble(group1 = pair[1], group2 = pair[2], p.value = pval)
    }
    pairwise_results <- purrr::map_dfr(lake_pairs, fisher_pair_test, df = data) %>%
      dplyr::mutate(p.adj = p.adjust(p.value, method = "BH"))
    
    letters <- multcompView::multcompLetters(
      setNames(pairwise_results$p.adj, paste(pairwise_results$group1, pairwise_results$group2, sep = "-"))
    )$Letters
    
    stats_col <- tibble::tibble(
      Lake   = factor(names(letters), levels = lakes),
      Letter = unname(letters)
    )
    
    ggplot2::ggplot(presence_summary, ggplot2::aes(x = Lake, y = Presence_Rate)) +
      ggplot2::geom_col(fill = "grey50", col = "black") +
      ggplot2::geom_text(data = stats_col, ggplot2::aes(x = Lake, y = Inf, label = Letter), vjust = 1.1, inherit.aes = FALSE) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.15))) +
      ggplot2::labs(x = xlab_text, y = y_label) +
      ggplot2::theme(legend.position = "none")
  }
}

KPRES_plot <- plot_koura_stats(M_C_data, "Presence_Kōura", "Kōura presence", test_type = "fisher", show_x_title = FALSE)
KCPUE_plot <- plot_koura_stats(M_C_data, "Weighted_CPUE_Kōura",     "Kōura CPUE",     test_type = "wilcox",  show_x_title = FALSE)
KBPUE_plot <- plot_koura_stats(M_C_data, "Weighted_BPUE_Kōura",     "Kōura BPUE",     test_type = "wilcox",  show_x_title = TRUE)

Koura_plots_0 <- KPRES_plot / KCPUE_plot / KBPUE_plot
Koura_plots_0

ggsave(file.path(fig_dir, "Koura_plots_0.png"), Koura_plots_0, width = 6, height = 6, dpi = 300)


```


# all collected variables plot
```{r fig.width=10, fig.height=10, dpi=100}
Env_data <- M_C_data %>%
  mutate(
    Emergent_vegetation  = Emergent_Native + Emergent_Non_Native,
    Submerged_vegetation = Submerged_Native + Submerged_Non_Native) %>%
  select(
    Lake, Slope_5m, Riparian_vegetation, Overhanging_trees, Wood_cover,
    Substrate_index, Temperature, DO_mgl, pH, Specific_conductivity,
    Emergent_vegetation, Submerged_vegetation) %>%
  pivot_longer(
    cols = -Lake,
    names_to = "Variable",
    values_to = "Value") %>%
  mutate(
    Lake     = factor(Lake, levels = lake_order),
    Var_type = "Env")

# 2. Fish presence (bars)
presence_data2 <- M_C_data %>%
  select(
    Monitoring_ID, Lake,
    Presence_Common_smelt, Presence_Goldfish,
    Presence_Eel, Presence_Catfish, Presence_Bullies,
    Presence_Trout, Presence_Mosquitofish,
    `Dragonfly larvae Small`, `Diving beetle`,
    Caddishfly, `Damselfly larvae`,
    `Dragonfly larvae Large`, Snails) %>%
  pivot_longer(
    cols = -c(Monitoring_ID, Lake),
    names_to  = "Fish_Type",
    values_to = "Presence") %>%
  mutate(
    Lake      = factor(Lake, levels = lake_order),
    Fish_Type = factor(Fish_Type, levels = fish_order))

Fish_presence_summary2 <- presence_data2 %>%
  group_by(Lake, Fish_Type) %>%
  summarise(
    Presence_prop = mean(Presence, na.rm = TRUE),
    .groups = "drop") %>%
  mutate(
    Variable = as.character(Fish_Type),
    Value    = Presence_prop,
    Var_type = "Fish") %>%
  select(Lake, Variable, Value, Var_type)

# 3. Combine into one dataset
EnvBio_data <- bind_rows(
  Env_data %>% select(Lake, Variable, Value, Var_type),
  Fish_presence_summary2) %>%
  mutate(
    Variable = factor(
      Variable,
      levels = c(
        "Substrate_index","Slope_5m","Riparian_vegetation",
        "Overhanging_trees","Wood_cover","Temperature","DO_mgl","pH",
        "Specific_conductivity","Emergent_vegetation","Submerged_vegetation",
        as.character(fish_order))))

# 4. Single plot with conditional geoms
EnvBio_plot <- ggplot(EnvBio_data, aes(Lake, Value)) +
  # boxplots only for environmental vars
  geom_boxplot( data = subset(EnvBio_data, Var_type == "Env")) +
  geom_col(data = subset(EnvBio_data, Var_type == "Fish"),fill = "grey40") +
  facet_wrap(~ Variable, scales = "free_y") 

EnvBio_plot

ggsave(
  file.path(fig_dir, "EnvBio_plot.png"),
  EnvBio_plot,
  width = 10,
  height = 10,
  dpi = 300)
```












# Physical parameters
```{r fig.width=14, fig.height=15, dpi=100}
# Filter to ignor all sites that have been sampled double.


# Lake order and colors
lake_order <- c("Rotorua", "Rotoiti", "Rotoehu", "Rotomā", "Ōkāreka")
M_C_data$Lake <- factor(M_C_data$Lake, levels = lake_order)

sediment_colors <- c(
  "Bedrock" = "black","Boulders" = "gray25","Cobble" = "gray55",
  "Gravel" = "gray80","Sand" = "gold","Mud" = "saddlebrown","Organic_matter" = "darkgreen"
)
weed_colors <- c(
  "Emergent_Native" = "darkgreen","Emergent_Non_Native" = "green3",
  "Submerged_Native" = "skyblue","Submerged_Non_Native" = "royalblue4",
  "Turf_Native" = "#00bfc4","Wood_cover" = "saddlebrown"
)
habitat_colors <- c("Rocky" = "gray25", "Sandy" = "gold", "Muddy" = "saddlebrown","Emergent Macrophyte" = "darkgreen")
habitat_order  <- c("Rocky", "Sandy", "Muddy", "Emergent Macrophyte")
M_C_data$Habitat_Type <- factor(M_C_data$Habitat_Type, levels = habitat_order)

# Tidy physical layers
sediment_data <- M_C_data %>%
  select(Site_ID_, DHT, Habitat_Type, Lake, Bedrock, Boulders, Cobble, Gravel, Sand, Mud, Organic_matter) %>%
  pivot_longer(
    cols = c(Bedrock, Boulders, Cobble, Gravel, Sand, Mud, Organic_matter),
    names_to = "Sediment_Type", values_to = "Percentage"
  ) %>%
  mutate(
    Lake = factor(Lake, levels = lake_order),
    Sediment_Type = factor(Sediment_Type, levels = c("Bedrock","Boulders","Cobble","Gravel","Sand","Mud","Organic_matter"))
  )

Substrate_index_data <- M_C_data

weed_data <- M_C_data %>%
  select(Site_ID_, DHT, Habitat_Type, Lake, Emergent_Native, Emergent_Non_Native,
         Submerged_Native, Submerged_Non_Native, Turf_Native, Wood_cover) %>%
  pivot_longer(
    cols = c(Emergent_Native,Emergent_Non_Native,Submerged_Native,Submerged_Non_Native,Turf_Native,Wood_cover),
    names_to = "Weeds", values_to = "Percentage"
  ) %>%
  mutate(Lake = factor(Lake, levels = lake_order))

# Statistics for substrate index (Wilcoxon pairwise + compact letters)
substrate_stats <- Substrate_index_data %>%
  rstatix::wilcox_test(Substrate_index ~ Lake, p.adjust.method = "BH") %>%
  select(group1, group2, p.adj) %>%
  { multcompView::multcompLetters(setNames(.$p.adj, paste(.$group1, .$group2, sep = "-")))$Letters } %>%
  enframe(name = "Lake", value = "Letter") %>%
  mutate(Label = paste0(Lake, " (", Letter, ")"))

# Plots
Sediment_plot <- ggplot(sediment_data, aes(factor(Site_ID_), Percentage, fill = Sediment_Type)) +
  geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.2) +
  scale_fill_manual(values = sediment_colors) +
  facet_grid(~ Lake, scales = "free_x") +
  labs(x = "Site", y = "Sediment Cover (%)", fill = "Sediment Type") +
  theme_bw() +
  theme(axis.title.x = element_blank())

Substrate_index_plot <- Substrate_index_data %>%
  mutate(Lake = factor(Lake, levels = substrate_stats$Lake)) %>%
  ggplot(aes(factor(Site_ID_), Substrate_index, fill = Habitat_Type)) +
  geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.2) +
  scale_fill_manual(values = habitat_colors) +
  facet_grid(~ Lake, scales = "free_x") +
  labs(x = "Site", y = "Substrate_index", fill = "Habitat Type") +
  theme_bw() +
  theme(axis.title.x = element_blank())

Weed_plot <- ggplot(weed_data, aes(factor(Site_ID_), Percentage, fill = Weeds)) +
  geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.2) +
  scale_fill_manual(values = weed_colors) +
  facet_grid(~ Lake, scales = "free_x") +
  labs(x = "Site", y = "Cover (%)", fill = "Cover Type") +
  theme_bw()+
  theme(axis.title.x = element_blank())

Koura_plot <- M_C_data %>%
  ggplot(aes(factor(Site_ID_), Weighted_CPUE_Kōura, fill = Habitat_Type)) +
  geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.2) +
  scale_fill_manual(values = habitat_colors) +
  facet_grid(~ Lake, scales = "free_x") +
  labs(x = "Site", y = "Kōura CPUE", fill = "Habitat Type") +
  theme_bw() +
  theme(axis.title.x = element_blank())

Catfish_plot <- M_C_data %>%
  ggplot(aes(factor(Site_ID_), Weighted_CPUE_Catfish, fill = Habitat_Type)) +
  geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.2) +
  scale_fill_manual(values = habitat_colors) +
  facet_grid(~ Lake, scales = "free_x") +
  labs(x = "Site", y = "Catfish CPUE", fill = "Habitat Type") +
  theme_bw() +
  theme(axis.title.x = element_blank())

Smelt_plot <- M_C_data %>%
  ggplot(aes(factor(Site_ID_), Presence_Common_smelt, fill = Habitat_Type)) +
  geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.2) +
  scale_fill_manual(values = habitat_colors) +
  facet_grid(~ Lake, scales = "free_x") +
  labs(x = "Site", y = "Smelt", fill = "Habitat Type") +
  theme_bw() +
  theme(axis.title.x = element_blank())

Temperature_plot <- M_C_data %>%
  ggplot(aes(factor(Site_ID_),Temperature, fill = Habitat_Type)) +
  geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.2) +
  scale_fill_manual(values = habitat_colors) +
  facet_grid(~ Lake, scales = "free_x") +
  theme_bw() +
  theme(axis.title.x = element_blank())

Specific_conductivity_plot <- M_C_data %>%
  ggplot(aes(factor(Site_ID_), Specific_conductivity, fill = Habitat_Type)) +
  geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.2) +
  scale_fill_manual(values = habitat_colors) +
  facet_grid(~ Lake, scales = "free_x") +
  theme_bw() +
  theme(axis.title.x = element_blank())

pH_plot <- M_C_data %>%
  ggplot(aes(factor(Site_ID_), pH, fill = Habitat_Type)) +
  geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.2) +
  scale_fill_manual(values = habitat_colors) +
  facet_grid(~ Lake, scales = "free_x") +
  labs(x = "Site")+
  theme_bw() 


Physical_plot <- Sediment_plot / Substrate_index_plot / Weed_plot / Koura_plot / Catfish_plot / Smelt_plot / Temperature_plot / Specific_conductivity_plot / pH_plot
Physical_plot

ggsave(file.path(fig_dir, "Physical_plot.png"), Physical_plot, width = 14, height = 10, dpi = 300)

```
# Chemical Parameters
```{r fig.width=12, fig.height=5, dpi=100}
Chemical_data <- M_C_data %>%
  select(Site_ID_, DHT, Habitat_Type, Lake, DO_mgl, DO_percent, Conductivity, Specific_conductivity, pH, Temperature) %>%
  pivot_longer(
    cols = c(DO_mgl, DO_percent, Conductivity, Specific_conductivity, pH, Temperature),
    names_to = "Variable", values_to = "Values"
  ) %>%
  mutate(Lake = factor(Lake, levels = lake_order))

# Compact letter display per variable (Wilcoxon + BH)
cld_df <- Chemical_data %>%
  group_by(Variable) %>%
  rstatix::wilcox_test(Values ~ Lake, p.adjust.method = "BH") %>%
  select(Variable, group1, group2, p.adj) %>%
  group_by(Variable) %>%
  summarise(Letters = list(multcompView::multcompLetters(setNames(p.adj, paste(group1, group2, sep = "-")))$Letters), .groups = "drop") %>%
  unnest_wider(Letters) %>%
  pivot_longer(-Variable, names_to = "Lake", values_to = "Letter") %>%
  left_join(
    Chemical_data %>% group_by(Variable, Lake) %>% summarise(y = max(Values, na.rm = TRUE), .groups = "drop"),
    by = c("Variable","Lake")
  ) %>%
  mutate(y = y * 1.05)

Chemical_plot <- ggplot(Chemical_data, aes(Lake, Values, fill = Lake)) +
  geom_boxplot() +
  facet_wrap(~ Variable, scales = "free", nrow = 1) +
  geom_text(data = cld_df, aes(Lake, y, label = Letter), vjust = 0) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

Chemical_plot

ggsave(file.path(fig_dir, "Chemical_plot.png"), Chemical_plot, width = 12, height = 5, dpi = 300)

```

# Biological parameters
```{r fig.width=12, fig.height=5, dpi=100}




CPUE_data <- M_C_data %>%
  select(Monitoring_ID, DHT, Habitat_Type, Lake,
         Weighted_CPUE_Kōura, Weighted_CPUE_Eel, Weighted_CPUE_Kōaro, Weighted_CPUE_Common_smelt, Weighted_CPUE_Bullies, Weighted_CPUE_Catfish, Weighted_CPUE_Goldfish, Weighted_CPUE_Mosquitofish, Weighted_CPUE_Trout) %>%
  pivot_longer(
    cols = c(Weighted_CPUE_Kōura, Weighted_CPUE_Eel, Weighted_CPUE_Kōaro, Weighted_CPUE_Common_smelt, Weighted_CPUE_Bullies, Weighted_CPUE_Catfish, Weighted_CPUE_Goldfish, Weighted_CPUE_Mosquitofish, Weighted_CPUE_Trout),
    names_to = "Fish_Type", values_to = "CPUE"
  ) %>%
  mutate(Lake = factor(Lake, levels = lake_order),
         Fish_Type = factor(Fish_Type, levels = fish_order_CPUE))

BCUE_data <- M_C_data %>%
  select(Monitoring_ID, DHT, Habitat_Type, Lake,
         Weighted_BCUE_Kōura, Weighted_BCUE_Eel, Weighted_BCUE_Kōaro, Weighted_BCUE_Common_smelt, Weighted_BCUE_Bullies, Weighted_BCUE_Catfish, Weighted_BCUE_Goldfish, Weighted_BCUE_Mosquitofish, Weighted_BCUE_Trout) %>%
  pivot_longer(
    cols = c(Weighted_BCUE_Kōura, Weighted_BCUE_Eel, Weighted_BCUE_Kōaro, Weighted_BCUE_Common_smelt, Weighted_BCUE_Bullies, Weighted_BCUE_Catfish, Weighted_BCUE_Goldfish, Weighted_BCUE_Mosquitofish, Weighted_BCUE_Trout),
    names_to = "Fish_Type", values_to = "BCUE"
  ) %>%
  mutate(Lake = factor(Lake, levels = lake_order),
         Fish_Type = factor(Fish_Type, levels = fish_order_BCUE))

BCUE_summary <- M_C_data %>%
  select(Monitoring_ID, DHT, Habitat_Type, Lake,
         Weighted_BCUE_Kōura, Weighted_BCUE_Eel, Weighted_BCUE_Kōaro, Weighted_BCUE_Common_smelt, Weighted_BCUE_Bullies, Weighted_BCUE_Catfish, Weighted_BCUE_Goldfish, Weighted_BCUE_Mosquitofish, Weighted_BCUE_Trout) %>%
  pivot_longer(
    cols = c(Weighted_BCUE_Kōura, Weighted_BCUE_Eel, Weighted_BCUE_Kōaro, Weighted_BCUE_Common_smelt, Weighted_BCUE_Bullies, Weighted_BCUE_Catfish, Weighted_BCUE_Goldfish, Weighted_BCUE_Mosquitofish, Weighted_BCUE_Trout),
    names_to = "Fish_Type", values_to = "BCUE"
  ) %>%
  mutate(Lake = factor(Lake, levels = lake_order),
         Fish_Type = factor(Fish_Type, levels = fish_order_BCUE)) %>%
  group_by(Lake, Habitat_Type, Fish_Type) %>%
  summarise(mean_BCUE = mean(BCUE, na.rm = TRUE),
            se_BCUE   = sd(BCUE,   na.rm = TRUE) / sqrt(n()),
            .groups = "drop")

Total_Weight_summary <- M_C_data %>%
  select(Monitoring_ID, DHT, Habitat_Type, Lake,
         Total_Weight_Bullies, Total_Weight_Goldfish, Total_Weight_Kōura,
         Total_Weight_Common_smelt, Total_Weight_Eel, Total_Weight_Kōaro,
         Total_Weight_Trout, Total_Weight_Mosquitofish, Total_Weight_Catfish) %>%
  pivot_longer(
    cols = c(Total_Weight_Bullies, Total_Weight_Goldfish, Total_Weight_Kōura,
             Total_Weight_Common_smelt, Total_Weight_Eel, Total_Weight_Kōaro,
             Total_Weight_Trout, Total_Weight_Mosquitofish, Total_Weight_Catfish),
    names_to = "Fish_Type", values_to = "Total_Weight"
  ) %>%
  mutate(Lake = factor(Lake, levels = lake_order)) %>%
  group_by(Fish_Type) %>%
  summarise(Total_Weight = sum(Total_Weight, na.rm = TRUE), .groups = "drop")

# Richness & Abundance
plot_abundance <- ggplot(presence_data, aes(Lake, Abundance, fill = Lake)) +
  geom_boxplot() +
  facet_grid(~ Habitat_Type, scales = "free")

plot_richness <- ggplot(presence_data, aes(Lake, Richness, fill = Lake)) +
  geom_boxplot() +
  facet_grid(~ Habitat_Type, scales = "free")

plot_abundance / plot_richness

```


































# CPUE kōura (abundances)
```{r fig.width=14, fig.height=6, dpi=100, echo=FALSE, message=FALSE, warning=FALSE}
## 1) Prep: select vars and coerce types
vars <- c("CPUE_Kōura","LID",
          "Slope_5m","Riparian_vegetation","Overhanging_trees","Wood_cover","Substrate_index",
          "Temperature","pH","DO_mgl","DO_percent","Specific_conductivity",
          "Emergent_Native","Submerged_Non_Native","Turf_Native","Submerged_Native","Emergent_Non_Native",
          "Presence_Common_smelt","Presence_Goldfish","Presence_Eel","Presence_Catfish")

CData <- Modeling_data |>
  dplyr::select(all_of(vars)) |>
  mutate(LID = factor(LID),
         across(c(Presence_Common_smelt, Presence_Goldfish,
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

mc.1 <- stepwise_gam_vars_with_LID(
  response = response, vars = vars_for_stepwise,
  data = CData, family = tw(link = "log"),
  method = "ML", p_cutoff = 0.1)

## 6) Compare models
summary(mc_full)
summary(mc_full_re)
summary(mc.05)
summary(mc.1)
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

final_plotc <- p1 | p2 | Predictions_plotc 
final_plotc

ggsave(file.path(fig_dir, "CPUE_Predictions.png"), final_plotc, width = 14, height = 6, dpi = 300)


```
# BCUE kōura (abundances)
```{r fig.width=14, fig.height=6, dpi=100, echo=FALSE, message=FALSE, warning=FALSE}
## 1) Prep: select vars and coerce types
vars <- c("BCUE_Kōura","LID",
          "Slope_5m","Riparian_vegetation","Overhanging_trees","Wood_cover","Substrate_index",
          "Temperature","pH","DO_mgl","DO_percent","Specific_conductivity",
          "Emergent_Native","Submerged_Non_Native","Turf_Native","Submerged_Native","Emergent_Non_Native",
          "Presence_Common_smelt","Presence_Goldfish","Presence_Eel","Presence_Catfish")

BData <- Modeling_data |>
  dplyr::select(all_of(vars)) |>
  mutate(LID = factor(LID),
         across(c(Presence_Common_smelt, Presence_Goldfish,
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

mb.1 <- stepwise_gam_vars_with_LID(
  response = response, vars = vars_for_stepwise,
  data = BData, family = tw(link = "log"),
  method = "ML", p_cutoff = 0.1)

## 6) Compare models
summary(mb_full)
summary(mb_full_re)
summary(mb.05)
summary(mb.1)
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
  Presence_Goldfish = c(0, 1),
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
p2 <- ggplot(newdata, aes(x = interaction(Presence_Goldfish, Presence_Eel), y = predicted_prob)) +
  geom_col(width = 0.5, fill = "grey70") +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2) +
  labs(x = "Presence Goldfish and Eel", y = "Predicted kōura BCUE") 

# Predicted and observed
pred <- predict(mb_final, type = "response")
obs <- mb_final$y
R2 <- cor(pred, obs)^2
RMSE <- sqrt(mean((pred - obs)^2))

R2
RMSE

df_pred_obs <- data.frame(ID = 1:length(obs), Predicted = pred, Observed = obs)

# Plot
library(ggrepel)
Predictions_plotb=ggplot(df_pred_obs, aes(x = Predicted, y = Observed)) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  #geom_text_repel(aes(label = ID), size = 3, max.overlaps = 5) +  
  annotate("text", x = max(df_pred_obs$Predicted) * 0.7, y = max(df_pred_obs$Observed) * 0.8,
           label = paste0("Adj. R² = ", round(R2, 2), "\nRMSE = ", round(RMSE, 2)),
           hjust = 0, size = 4, colour = "black") +
  labs(x = "Predicted BCUE", y = "Observed BCUE")

Predictions_plotb

final_plotb <- p1 | p2 | Predictions_plotb 
final_plotb

ggsave(file.path(fig_dir, "BCUE_Predictions.png"), final_plotb, width = 14, height = 6, dpi = 300)

```



# remake the CPUE and BCUE
```{r fig.width=14, fig.height=6, dpi=100, echo=FALSE, message=FALSE, warning=FALSE}
# ---------- Shared helpers ----------
is_cont <- function(x) is.numeric(x) && dplyr::n_distinct(x, na.rm = TRUE) >= 5

#build_smooth <- function(var, data, custom_k = list(), include_re_for_LID = FALSE) {
if (identical(var, "LID") && include_re_for_LID) return("s(LID, bs='re')")

x <- data[[var]]

# Numeric? Try a smooth with safe k; else parametric.
if (is.numeric(x)) {
  nuniq <- dplyr::n_distinct(x, na.rm = TRUE)
  
  # If too few uniques for a smooth, use parametric (linear) term
  if (nuniq < 5) return(var)
  
  # Requested k (default 10), but cap to <= nuniq - 1 and >= 3
  k_req <- if (!is.null(custom_k[[var]])) custom_k[[var]] else 10
  k_cap <- max(3, min(k_req, nuniq - 1))
  
  return(paste0("s(", var, ", bs='ts', k=", k_cap, ")"))
}

# Binary/low-level factors -> parametric
var
}

mk_rhs <- function(vars, data, custom_k = list(), allow_re = TRUE) {
  paste(vapply(vars, function(v) build_smooth(v, data, custom_k, include_re_for_LID = allow_re), character(1)),
        collapse = " + ")
}

fit_binom <- function(formula, data) {
  mgcv::gam(formula, data = data, family = binomial(link = "logit"),
            method = "REML", select = TRUE)
}

fit_gamma <- function(formula, data) {
  mgcv::gam(formula, data = data, family = Gamma(link = "log"),
            method = "REML", select = TRUE)
}

logit2prob <- function(x) exp(x)/(1+exp(x))

stepwise_gam_vars_with_LID <- function(response, vars, data, family,
                                       method = "ML", p_cutoff = 0.1,
                                       custom_k = list()) {
  remaining_vars <- vars
  m <- NULL
  repeat {
    rhs_terms <- vapply(
      remaining_vars,
      function(v) build_smooth(v, data = data, custom_k = custom_k, include_re_for_LID = TRUE),
      character(1)
    )
    current_formula <- as.formula(paste(response, "~", paste(rhs_terms, collapse = " + ")))
    m <- mgcv::gam(current_formula, data = data, family = family, method = method, select = TRUE)
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
    var_to_remove <- if (grepl("^s\\(", term_to_remove))
      sub("^s\\(([^,]+).*\\)$", "\\1", term_to_remove) else term_to_remove
    
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


mk_rhs <- function(vars, data, custom_k = list(), allow_re = TRUE) {
  paste(vapply(vars, function(v) build_smooth(v, data = data, custom_k = custom_k,
                                              include_re_for_LID = allow_re),
               character(1)),
        collapse = " + ")
}

# Extract significant smooth and parametric terms from an mgcv::gam
get_sig_terms <- function(m, alpha = 0.05){
  sm <- summary(m)
  sig_smooths <- if (!is.null(sm$s.table)) {
    rn <- rownames(sm$s.table)
    rn[sm$s.table[,"p-value"] < alpha]
  } else character(0)
  
  sig_params <- if (!is.null(sm$p.table) && nrow(sm$p.table) > 0) {
    pcol <- intersect(colnames(sm$p.table), c("Pr(>|t|)","Pr(>|z|)"))[1]
    rn <- rownames(sm$p.table)
    rn[rn != "(Intercept)" & sm$p.table[, pcol] < alpha]
  } else character(0)
  
  list(smooths = sig_smooths, params = sig_params)
}

# Draw only the significant smooths (if any)
draw_sig_smooths <- function(m, alpha = 0.05){
  sig <- get_sig_terms(m, alpha)$smooths
  if (!length(sig)) return(NULL)
  
  # map smooth labels in summary to draw()'s "select"
  labs <- unique(smooth_estimates(m)$smooth)   # e.g., "s(Temperature)"
  idx  <- which(labs %in% sig)
  if (!length(idx)) return(NULL)
  
  draw(m, select = idx)
}

# Build a simple base grid with everything at means and all LIDs to average over
# 'req' are term labels from both models you need to supply (union of terms)
make_base_grid <- function(df, req){
  out <- list()
  for (v in req) {
    if (!v %in% names(df)) next
    x <- df[[v]]
    if (is.numeric(x)) {
      out[[v]] <- mean(x, na.rm = TRUE)
    } else if (is.factor(x)) {
      out[[v]] <- factor(levels(x)[1], levels = levels(x))
    } else {
      out[[v]] <- 0
    }
  }
  # Keep all lakes to average over RE when excluded
  if ("LID" %in% names(df)) out$LID <- levels(df$LID)
  as.data.frame(out)
}

# Exclude s(LID) only if present in the model (avoids warnings)
exclude_if_RE <- function(mod){
  if (!length(mod$smooth)) return(NULL)
  incl <- vapply(mod$smooth, function(s) "LID" %in% s$term, logical(1))
  if (any(incl)) "s(LID)" else NULL
}

# Plot overall effect (p * mu) for a single binary var using the two-part models
plot_overall_binary <- function(var, pres_model, pos_model, df, ylab = "Overall expected value"){
  # figure out which terms were used
  req_pres <- attr(terms(pres_model), "term.labels")
  req_pos  <- attr(terms(pos_model),  "term.labels")
  req_all  <- union(req_pres, req_pos)
  
  base <- make_base_grid(df, req_all)
  
  # make two grids differing only in var (0 vs 1)
  new0 <- new1 <- base
  if (!var %in% names(new0)) new0[[var]] <- 0
  if (!var %in% names(new1)) new1[[var]] <- 1
  
  # predictions (averaged over LID by excluding s(LID))
  p0 <- predict(pres_model, newdata = new0, type = "response", exclude = exclude_if_RE(pres_model))
  p1 <- predict(pres_model, newdata = new1, type = "response", exclude = exclude_if_RE(pres_model))
  mu0 <- predict(pos_model,  newdata = new0, type = "response", exclude = exclude_if_RE(pos_model))
  mu1 <- predict(pos_model,  newdata = new1, type = "response", exclude = exclude_if_RE(pos_model))
  
  dfp <- tibble(level = factor(c(0,1)),
                overall = c(mean(p0 * mu0, na.rm = TRUE), mean(p1 * mu1, na.rm = TRUE)))
  
  ggplot(dfp, aes(level, overall)) +
    geom_col(width = 0.6, fill = "grey70") +
    labs(x = var, y = ylab)
}


# ---------- CPUE data ----------
vars_c <- c("CPUE_Kōura","LID",
            "Slope_5m","Riparian_vegetation","Overhanging_trees","Wood_cover","Substrate_index",
            "Temperature","pH","DO_mgl","DO_percent","Specific_conductivity",
            "Emergent_Native","Submerged_Non_Native","Turf_Native","Submerged_Native","Emergent_Non_Native",
            "Presence_Common_smelt","Presence_Goldfish","Presence_Eel","Presence_Catfish")

CData <- Modeling_data |>
  dplyr::select(all_of(vars_c)) |>
  mutate(
    LID = factor(LID),
    across(c(Presence_Common_smelt, Presence_Goldfish, Presence_Eel, Presence_Catfish), ~ as.numeric(.)),
    CPUE_pos = as.integer(CPUE_Kōura > 0)
  )

# Custom k (same spirit as yours)
custom_k <- list(Slope_5m=10, Riparian_vegetation=7, Overhanging_trees=5, Wood_cover=10,
                 Substrate_index=10, Temperature=10, pH=10, DO_mgl=10,
                 Emergent_Native=9, Submerged_Native=4, Submerged_Non_Native=9, Turf_Native=6)

# ---------- VIF pruning on predictors (exclude response & LID) ----------
resp_posflag <- "CPUE_pos"
pred_fixed_c <- setdiff(names(CData), c("CPUE_Kōura", resp_posflag, "LID"))
vif_c <- remove_high_vif_glmmTMB(CData, "CPUE_Kōura", pred_fixed_c, threshold = 5)
kept_fixed_c <- setdiff(names(fixef(vif_c)$cond), "(Intercept)")

# Use same candidate set for both hurdle parts:
vars_for_stepwise_c <- c(kept_fixed_c, "LID")

# ---------- Part 1: presence (CPUE > 0) ----------
rhs1 <- mk_rhs(kept_fixed_c, CData, custom_k, allow_re = FALSE)
form1_noRE <- as.formula(paste0(resp_posflag, " ~ ", rhs1))
form1_RE   <- as.formula(paste0(resp_posflag, " ~ ", rhs1, " + s(LID, bs='re')"))

# Stepwise on binomial with potential RE
mc.pres <- stepwise_gam_vars_with_LID(
  response = resp_posflag, vars = vars_for_stepwise_c,
  data = CData, family = binomial(link = "logit"),
  method = "ML", p_cutoff = 0.05, custom_k = custom_k)

mc.pres_final <- fit_binom(formula(mc.pres), CData)
summary(mc.pres_final); gam.check(mc.pres_final); concurvity(mc.pres_final)

# ---------- Part 2: positive abundance (CPUE | CPUE>0), Gamma(log) ----------
Cpos <- CData |> filter(CPUE_pos == 1)

rhs2 <- mk_rhs(kept_fixed_c, Cpos, custom_k, allow_re = FALSE)
form2_noRE <- as.formula(paste0("CPUE_Kōura ~ ", rhs2))
form2_RE   <- as.formula(paste0("CPUE_Kōura ~ ", rhs2, " + s(LID, bs='re')"))

mc.pos <- stepwise_gam_vars_with_LID(
  response = "CPUE_Kōura", vars = vars_for_stepwise_c,
  data = Cpos, family = Gamma(link = "log"),
  method = "ML", p_cutoff = 0.05)

mc.pos_final <- fit_gamma(formula(mc.pos), Cpos)
summary(mc.pos_final); gam.check(mc.pos_final); concurvity(mc.pos_final)

# ---------- Combine to overall expected CPUE ----------
# p_pos = Pr(CPUE>0), mu_pos = E[CPUE | CPUE>0]
# Overall mean = p_pos * mu_pos
# Example prediction grid for fish presence effect (hold other terms at means; exclude RE for clarity)
new_cpue <- expand.grid(
  Presence_Goldfish = 0:1,
  Presence_Eel      = 0:1,
  Presence_Catfish  = 0:1,
  Presence_Common_smelt = 0:1,
  LID = levels(CData$LID),
  Riparian_vegetation = mean(CData$Riparian_vegetation, na.rm = TRUE),
  Substrate_index     = mean(CData$Substrate_index,     na.rm = TRUE),
  Temperature         = mean(CData$Temperature,         na.rm = TRUE),
  pH                  = mean(CData$pH,                  na.rm = TRUE),
  DO_mgl              = mean(CData$DO_mgl,              na.rm = TRUE),
  Specific_conductivity = mean(CData$Specific_conductivity, na.rm = TRUE),
  Slope_5m            = mean(CData$Slope_5m,            na.rm = TRUE),
  Overhanging_trees   = mean(CData$Overhanging_trees,   na.rm = TRUE),
  Wood_cover          = mean(CData$Wood_cover,          na.rm = TRUE),
  Emergent_Native     = mean(CData$Emergent_Native,     na.rm = TRUE),
  Submerged_Native    = mean(CData$Submerged_Native,    na.rm = TRUE),
  Submerged_Non_Native= mean(CData$Submerged_Non_Native,na.rm = TRUE),
  Turf_Native         = mean(CData$Turf_Native,         na.rm = TRUE)
)

p_link <- predict(mc.pres_final, newdata = new_cpue, type = "link", se.fit = TRUE, exclude = "s(LID)")
new_cpue$p_pos  <- logit2prob(p_link$fit)

mu_pos <- predict(mc.pos_final, newdata = new_cpue, type = "response", se.fit = TRUE, exclude = "s(LID)")
new_cpue$mu_pos <- mu_pos$fit

new_cpue$overall_mean <- new_cpue$p_pos * new_cpue$mu_pos

# Summarise across lakes to plot clean contrasts
eff_cpue <- new_cpue |>
  dplyr::group_by(Presence_Goldfish, Presence_Eel, Presence_Catfish, Presence_Common_smelt) |>
  dplyr::summarise(overall = mean(overall_mean), .groups = "drop")

pred_pos  <- predict(mc.pos_final, type = "response")
obs_pos   <- mc.pos_final$y
R2_pos    <- cor(pred_pos, obs_pos)^2
RMSE_pos  <- sqrt(mean((pred_pos - obs_pos)^2))




# ---------- BCUE data ----------
vars_b <- c("BCUE_Kōura","LID",
            "Slope_5m","Riparian_vegetation","Overhanging_trees","Wood_cover","Substrate_index",
            "Temperature","pH","DO_mgl","DO_percent","Specific_conductivity",
            "Emergent_Native","Submerged_Non_Native","Turf_Native","Submerged_Native","Emergent_Non_Native",
            "Presence_Common_smelt","Presence_Goldfish","Presence_Eel","Presence_Catfish")

BData <- Modeling_data |>
  dplyr::select(all_of(vars_b)) |>
  mutate(
    LID = factor(LID),
    across(c(Presence_Common_smelt, Presence_Goldfish, Presence_Eel, Presence_Catfish), ~ as.numeric(.)),
    BCUE_pos = as.integer(BCUE_Kōura > 0)
  )

# VIF on predictors
pred_fixed_b <- setdiff(names(BData), c("BCUE_Kōura", "BCUE_pos", "LID"))
vif_b <- remove_high_vif_glmmTMB(BData, "BCUE_Kōura", pred_fixed_b, threshold = 5)
kept_fixed_b <- setdiff(names(fixef(vif_b)$cond), "(Intercept)")
vars_for_stepwise_b <- c(kept_fixed_b, "LID")

# Part 1: presence (BCUE > 0)
mb.pres <- stepwise_gam_vars_with_LID(
  response = "BCUE_pos", vars = vars_for_stepwise_b,
  data = BData, family = binomial(link = "logit"),
  method = "ML", p_cutoff = 0.05)

mb.pres_final <- fit_binom(formula(mb.pres), BData)
summary(mb.pres_final); gam.check(mb.pres_final); concurvity(mb.pres_final)

# Part 2: positive biomass (Gamma log) OR lognormal sensitivity
Bpos <- BData |> filter(BCUE_pos == 1)

mb.pos <- stepwise_gam_vars_with_LID(
  response = "BCUE_Kōura", vars = vars_for_stepwise_b,
  data = Bpos, family = Gamma(link = "log"),
  method = "ML", p_cutoff = 0.05)

mb.pos_final <- fit_gamma(formula(mb.pos), Bpos)
#summary(mb.pos_final); gam.check(mb.pos_final); concurvity(mb.pos_final)

# Combine to overall expected BCUE
new_bcue <- expand.grid(
  Presence_Goldfish = 0:1,
  Presence_Eel      = 0:1,
  Presence_Catfish  = 0:1,
  Presence_Common_smelt = 0:1,
  LID = levels(BData$LID),
  Riparian_vegetation = mean(BData$Riparian_vegetation, na.rm = TRUE),
  Substrate_index     = mean(BData$Substrate_index,     na.rm = TRUE),
  Temperature         = mean(BData$Temperature,         na.rm = TRUE),
  pH                  = mean(BData$pH,                  na.rm = TRUE),
  DO_mgl              = mean(BData$DO_mgl,              na.rm = TRUE),
  Overhanging_trees   = mean(BData$Overhanging_trees,   na.rm = TRUE),
  Wood_cover          = mean(BData$Wood_cover,          na.rm = TRUE),
  Emergent_Native     = mean(BData$Emergent_Native,     na.rm = TRUE),
  Emergent_Non_Native   = mean(BData$Emergent_Non_Native,   na.rm = TRUE),
  Submerged_Native    = mean(BData$Submerged_Native,    na.rm = TRUE),
  Submerged_Non_Native= mean(BData$Submerged_Non_Native,na.rm = TRUE),
  Turf_Native         = mean(BData$Turf_Native,         na.rm = TRUE)
)

p_link_b <- predict(mb.pres_final, newdata = new_bcue, type = "link", se.fit = TRUE, exclude = "s(LID)")
new_bcue$p_pos  <- logit2prob(p_link_b$fit)

mu_pos_b <- predict(mb.pos_final, newdata = new_bcue, type = "response", se.fit = TRUE, exclude = "s(LID)")
new_bcue$mu_pos <- mu_pos_b$fit

new_bcue$overall_mean <- new_bcue$p_pos * new_bcue$mu_pos

eff_bcue <- new_bcue |>
  dplyr::group_by(Presence_Goldfish, Presence_Eel, Presence_Catfish, Presence_Common_smelt) |>
  dplyr::summarise(overall = mean(overall_mean), .groups = "drop")






sig_smooth_labels <- function(m, alpha = 0.05){
  sm <- summary(m)
  if (is.null(sm$s.table)) return(character(0))
  # rownames look like "s(var)" already — exactly what draw(select=) accepts
  labs <- rownames(sm$s.table)
  labs[ sm$s.table[, "p-value"] < alpha ]
}

draw_sig_smooths <- function(m, alpha = 0.05){
  labs <- sig_smooth_labels(m, alpha)
  if (!length(labs)) return(NULL)     # no significant smooths
  draw(m, select = labs)
}

# Optional: see what will be drawn
list_sig_terms <- function(..., alpha = 0.05){
  mods <- list(...)
  lapply(mods, function(m){
    list(
      model = deparse(formula(m)),
      sig_smooths = sig_smooth_labels(m, alpha),
      sig_params  = {
        sm <- summary(m)
        if (is.null(sm$p.table) || nrow(sm$p.table) == 0) character(0) else {
          pcol <- intersect(colnames(sm$p.table), c("Pr(>|t|)", "Pr(>|z|)"))[1]
          rn <- rownames(sm$p.table)
          rn[rn != "(Intercept)" & sm$p.table[, pcol] < alpha]
        }
      }
    )
  })
}


alpha <- 0.05

# See what's significant first (helps debug)
list_sig_terms(mc.pos_final, mc.pres_final, mb.pos_final, mb.pres_final, alpha = alpha)

# Draw only significant smooths, if any
p_cpue_smooths <- {
  up <- draw_sig_smooths(mc.pos_final, alpha)
  lo <- draw_sig_smooths(mc.pres_final, alpha)
  if (is.null(up) && is.null(lo)) NULL else up / lo
}

p_bcue_smooths <- {
  up <- draw_sig_smooths(mb.pos_final, alpha)
  lo <- draw_sig_smooths(mb.pres_final, alpha)
  if (is.null(up) && is.null(lo)) NULL else up / lo
}
# ---- Significant fish binaries only (presence vars) ----
fish_vars <- c("Presence_Goldfish","Presence_Eel","Presence_Catfish","Presence_Common_smelt")

sig_cpue <- get_sig_terms(mc.pos_final, alpha);  # positive part
sig_cpue_p <- get_sig_terms(mc.pres_final, alpha) # presence part
sig_fish_cpue <- intersect(fish_vars, unique(c(sig_cpue$params, sig_cpue_p$params)))

sig_bcue <- get_sig_terms(mb.pos_final, alpha);
sig_bcue_p <- get_sig_terms(mb.pres_final, alpha)
sig_fish_bcue <- intersect(fish_vars, unique(c(sig_bcue$params, sig_bcue_p$params)))

# Build plots only for significant fish vars
plots_cpue <- lapply(sig_fish_cpue, function(v) plot_overall_binary(v, mc.pres_final, mc.pos_final, CData,
                                                                    ylab = "Overall expected CPUE"))
plots_bcue <- lapply(sig_fish_bcue, function(v) plot_overall_binary(v, mb.pres_final, mb.pos_final, BData,
                                                                    ylab = "Overall expected BCUE"))

# If you use patchwork:
# library(patchwork)
p_cpue_final <- wrap_plots(plots_cpue) | p_cpue_smooths
p_bcue_final <- wrap_plots(plots_bcue) | p_bcue_smooths










# GAM Helpers
```{r}
Modeling_data <- Monitoring_CPUE_data %>%
  filter(Monitoring==0)

alpha_sig   <- 0.05
p_cutoff_ml <- 0.05        
vif_thresh  <- 5
INCLUDE_RE   <- TRUE          
PROTECT_FISH <- FALSE 

# Optional: per-variable k (capped safely by unique-1)
custom_k <- list(
  Slope_5m=10, Riparian_vegetation=7, Overhanging_trees=5, Wood_cover=10,
  Substrate_index=10, Temperature=10, pH=10, DO_mgl=10,
  Emergent_Native=9, Submerged_Non_Native=9, Turf_Native=6, Submerged_Native=6)

fish_vars <- c("Presence_Goldfish","Presence_Eel","Presence_Catfish","Presence_Common_smelt")


# Variables per analysis
vars_common <- c(
  "LID",
  "Slope_5m","Riparian_vegetation","Overhanging_trees","Wood_cover","Substrate_index",
  "Temperature","pH","DO_mgl","DO_percent","Specific_conductivity",
  "Emergent_Native","Submerged_Non_Native","Turf_Native","Submerged_Native","Emergent_Non_Native",
  fish_vars)


is_cont <- function(x) is.numeric(x) && dplyr::n_distinct(x, na.rm = TRUE) >= 5
logit2prob <- function(x) exp(x) / (1 + exp(x))

# Smooth builder
build_smooth_gam <- function(var, data, custom_k = list(), include_re_for_LID = FALSE) {
  if (identical(var, "LID") && include_re_for_LID) return("s(LID, bs='re')")
  x <- data[[var]]
  if (is.numeric(x)) {
    nuniq <- dplyr::n_distinct(x, na.rm = TRUE)
    if (nuniq < 5) return(var)
    k_req <- if (!is.null(custom_k[[var]])) custom_k[[var]] else 10
    k_cap <- max(3, min(k_req, nuniq - 1))
    return(paste0("s(", var, ", bs='ts', k=", k_cap, ")"))
  }
  var
}

# Build FULL RHS from a set of variables
make_full_rhs <- function(vars, data, include_re = TRUE, custom_k = list()){
  paste(
    vapply(vars, function(v) build_smooth_gam(v, data, custom_k, include_re_for_LID = include_re),
           character(1)),
    collapse = " + "
  )
}

# Exclude s(LID) only if present
exclude_if_RE <- function(mod){
  if (!length(mod$smooth)) return(NULL)
  has_re <- vapply(mod$smooth, function(s) "LID" %in% s$term, logical(1))
  if (any(has_re)) "s(LID)" else NULL
}

# VIF pruning — with protect option
remove_high_vif_glmBI <- function(data, response, predictors, threshold = 5, protect_vars = character(0)) {
  if (!requireNamespace("brglm2", quietly = TRUE))
    stop("Package 'brglm2' is required for bias-reduced logistic VIF pruning.")
  nzv <- caret::nearZeroVar(data[, predictors, drop = FALSE])
  if (length(nzv)) predictors <- predictors[-nzv]
  
  rpt <- list(removed = character(), start = predictors)
  repeat {
    if (!length(predictors)) stop("All predictors removed during VIF pruning.")
    fml <- as.formula(paste(response, "~", paste(predictors, collapse = " + ")))
    model <- try(glm(fml, data = data, family = binomial(link = "logit"), method = "brglmFit"), silent = TRUE)
    if (inherits(model, "try-error")) break
    vif_data <- performance::check_collinearity(model)
    vif_data <- vif_data[!grepl("\\|", vif_data$Term), ]
    if (!nrow(vif_data) || all(vif_data$VIF < threshold)) break
    
    ord <- order(vif_data$VIF, decreasing = TRUE)
    to_remove <- NA_character_
    for (i in ord) {
      cand <- vif_data$Term[i]
      if (!(cand %in% protect_vars)) { to_remove <- cand; break }
    }
    if (is.na(to_remove)) break
    predictors <- setdiff(predictors, to_remove)
    rpt$removed <- c(rpt$removed, to_remove)
  }
  attr(model, "vif_report") <- rpt
  model
}

remove_high_vif_glmmTMB <- function(data, response, predictors, threshold = 5, protect_vars = character(0)) {
  nzv <- caret::nearZeroVar(data[, predictors, drop = FALSE])
  if (length(nzv)) predictors <- predictors[-nzv]
  
  rpt <- list(removed = character(), start = predictors)
  repeat {
    if (!length(predictors)) stop("All predictors removed during VIF pruning.")
    fml <- as.formula(paste(response, "~", paste(predictors, collapse = " + ")))
    model <- try(glmmTMB::glmmTMB(fml, data = data, family = glmmTMB::tweedie()), silent = TRUE)
    if (inherits(model, "try-error")) break
    vif_data <- performance::check_collinearity(model)
    vif_data <- vif_data[!grepl("\\|", vif_data$Term), ]
    if (!nrow(vif_data) || all(vif_data$VIF < threshold)) break
    
    ord <- order(vif_data$VIF, decreasing = TRUE)
    to_remove <- NA_character_
    for (i in ord) {
      cand <- vif_data$Term[i]
      if (!(cand %in% protect_vars)) { to_remove <- cand; break }
    }
    if (is.na(to_remove)) break
    predictors <- setdiff(predictors, to_remove)
    rpt$removed <- c(rpt$removed, to_remove)
  }
  attr(model, "vif_report") <- rpt
  model
}

# Stepwise (ML) with optional mandatory vars protected from dropping
stepwise_gam_vars_with_LID <- function(response, vars, data, family,
                                       method = "ML", p_cutoff = 0.1,
                                       custom_k = list(),
                                       mandatory_vars = character(0)) {
  remaining_vars <- vars
  drops <- character()
  m <- NULL
  repeat {
    if (!length(remaining_vars)) { message("No variables left; stopping."); break }
    
    rhs_terms <- vapply(
      remaining_vars,
      function(v) build_smooth_gam(v, data = data, custom_k = custom_k, include_re_for_LID = TRUE),
      character(1)
    )
    current_formula <- as.formula(paste(response, "~", paste(rhs_terms, collapse = " + ")))
    m <- mgcv::gam(current_formula, data = data, family = family, method = method, select = TRUE)
    sm <- summary(m)
    
    # Collect p-values
    param_p <- if (!is.null(sm$p.table) && nrow(sm$p.table) > 0) {
      pcol <- intersect(colnames(sm$p.table), c("Pr(>|t|)","Pr(>|z|)"))[1]
      pp <- sm$p.table[, pcol]; names(pp) <- rownames(sm$p.table); pp[names(pp) != "(Intercept)"]
    } else numeric(0)
    
    smooth_p <- if (!is.null(sm$s.table)) {
      p <- sm$s.table[, "p-value"]; names(p) <- rownames(sm$s.table); p
    } else numeric(0)
    
    all_p <- c(param_p, smooth_p)
    drop <- all_p[all_p > p_cutoff]
    if (!length(drop)) { message("All terms have p <= ", p_cutoff, "; stopping."); break }
    
    # choose worst non-mandatory
    ordered <- names(sort(drop, decreasing = TRUE))
    cand <- setdiff(ordered, mandatory_vars)
    if (!length(cand)) { message("Only protected variables remain; stopping."); break }
    
    term_to_remove <- cand[1]
    var_to_remove  <- if (grepl("^s\\(", term_to_remove)) sub("^s\\(([^,]+).*\\)$", "\\1", term_to_remove) else term_to_remove
    
    if (var_to_remove %in% remaining_vars) {
      remaining_vars <- setdiff(remaining_vars, var_to_remove)
      drops <- c(drops, var_to_remove)
    } else {
      message("Variable not found in current set. Stopping."); break
    }
  }
  attr(m, "dropped_terms") <- drops
  m
}

# Extraction helpers
term_table <- function(m){
  sm <- summary(m)
  param <- if (!is.null(sm$p.table) && nrow(sm$p.table) > 0) {
    pcol <- intersect(colnames(sm$p.table), c("Pr(>|t|)","Pr(>|z|)"))[1]
    tibble(term = rownames(sm$p.table),
           type = "parametric",
           p    = sm$p.table[, pcol]) %>%
      filter(term != "(Intercept)")
  } else tibble(term=character(), type=character(), p=numeric())
  
  smooth <- if (!is.null(sm$s.table)) {
    tibble(term = rownames(sm$s.table),
           type = "smooth",
           edf  = sm$s.table[, "edf"],
           p    = sm$s.table[, "p-value"])
  } else tibble(term=character(), type=character(), edf=numeric(), p=numeric())
  
  suppressWarnings(full_join(param, smooth, by = c("term","type","p")))
}

has_term <- function(mod, var){
  tl <- attr(terms(mod), "term.labels")
  any(tl == var) || any(grepl(paste0("^s\\(", var, "(,|\\))"), tl))
}

# Plot helpers
as_plot <- function(p) {
  if (is.null(p)) return(patchwork::plot_spacer())
  if (inherits(p, c("gg","ggplot","patchwork"))) return(p)
  if (is.list(p)) {
    if (length(p) == 0) return(patchwork::plot_spacer())
    return(patchwork::wrap_plots(p))
  }
  patchwork::plot_spacer()
}

placeholder <- function(txt) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = txt, size = 4) +
    xlim(0,1) + ylim(0,1) + theme_void()
}

plot_binary_single_component <- function(model, var, data,
                                         family = c("binomial","gamma"),
                                         S = 2000, alpha = 0.05,
                                         exclude_RE = TRUE,
                                         hold_binaries = c("mean","zero","one"),
                                         id = "LID",
                                         fish_vars = fish_vars,
                                         seed = 1) {
  family <- match.arg(family)
  hold_binaries <- match.arg(hold_binaries)
  
  if (!has_term(model, var))
    return(placeholder(paste(var, "not retained in model")))
  
  # Base at numeric means
  num_means <- data %>% summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))
  base <- as.list(num_means)
  if (id %in% names(data) && is.factor(data[[id]]))
    base[[id]] <- factor(levels(data[[id]])[1], levels = levels(data[[id]]))
  
  # Hold other binaries
  bin_fish <- intersect(fish_vars, names(data))
  others <- setdiff(bin_fish, var)
  set_bin <- function(x){
    if (hold_binaries == "mean") round(mean(data[[x]], na.rm = TRUE))
    else if (hold_binaries == "zero") 0L else 1L
  }
  for (x in others) base[[x]] <- set_bin(x)
  
  # newdata for 0/1 replicated across LID
  mk_nd <- function(level){
    nd <- as.data.frame(base, stringsAsFactors = FALSE)
    nd[[var]] <- as.integer(level)
    if (id %in% names(data) && is.factor(data[[id]])) {
      nd <- do.call(rbind, lapply(levels(data[[id]]), function(lv){
        r <- nd; r[[id]] <- factor(lv, levels = levels(data[[id]])); r
      }))
    }
    # Align factor levels to model frame if available
    if (!is.null(model$model)) {
      common <- intersect(names(nd), names(model$model))
      for (nm in common) {
        if (is.factor(model$model[[nm]])) {
          nd[[nm]] <- factor(nd[[nm]], levels = levels(model$model[[nm]]))
        }
      }
    }
    nd
  }
  nd0 <- mk_nd(0L); nd1 <- mk_nd(1L)
  nd <- bind_rows(mutate(nd0, .level = 0L), mutate(nd1, .level = 1L))
  
  excl <- if (exclude_RE) exclude_if_RE(model) else NULL
  lp <- predict(model, newdata = nd, type = "link", se.fit = TRUE, exclude = excl)
  
  set.seed(seed)
  n <- nrow(nd)
  Z <- matrix(rnorm(n * S, lp$fit, pmax(lp$se.fit, .Machine$double.eps)), nrow = n, ncol = S)
  
  if (family == "binomial") {
    Y <- plogis(Z); ylab <- "Presence probability"
  } else {
    Y <- exp(Z);    ylab <- "Mean given presence (μ)"
  }
  
  sim_by_level <- function(level_flag){
    sims <- Y[nd$.level == level_flag, , drop = FALSE]
    colMeans(sims)
  }
  sim0 <- sim_by_level(0L); sim1 <- sim_by_level(1L)
  
  qlo <- alpha/2; qhi <- 1 - alpha/2
  df <- data.frame(
    level = factor(c(0,1), levels = c(0,1), labels = c("Absent","Present")),
    mean  = c(mean(sim0), mean(sim1)),
    lwr   = c(quantile(sim0, qlo), quantile(sim1, qlo)),
    upr   = c(quantile(sim0, qhi), quantile(sim1, qhi))
  )
  
  ggplot(df, aes(x = level, y = mean)) +
    geom_col(width = 0.6) +
    geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2) +
    labs(x = var, y = ylab)
}

diag_roc <- function(obs, pred){
  roc_curve <- pROC::roc(obs, pred)
  auc_val <- pROC::auc(roc_curve)
  roc_df <- data.frame(tpr = roc_curve$sensitivities, fpr = 1 - roc_curve$specificities)
  p <- ggplot(roc_df, aes(fpr, tpr)) +
    geom_path() + geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    annotate("text", x = 0.6, y = 0.1, label = paste0("AUC = ", round(as.numeric(auc_val), 3))) +
    labs(x = "False positive rate", y = "True positive rate")
  list(plot = p, auc = as.numeric(auc_val))
}

diag_scatter <- function(obs, pred, xlab, ylab){
  R2   <- cor(pred, obs)^2
  RMSE <- sqrt(mean((pred - obs)^2))
  df   <- data.frame(Predicted = pred, Observed = obs)
  p <- ggplot(df, aes(Predicted, Observed)) +
    geom_point(alpha = 0.6) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    annotate("text",
             x = max(df$Predicted, na.rm = TRUE) * 0.7,
             y = max(df$Observed,  na.rm = TRUE) * 0.8,
             label = paste0("R² = ", round(R2, 2), "\nRMSE = ", round(RMSE, 2))) +
    labs(x = xlab, y = ylab)
  list(plot = p, R2 = R2, RMSE = RMSE)
}

# Selection summary printer
print_selection_report <- function(title, vif_report, full_formula, full_ml, red_ml, final_reml){
  cat("\n==================== ", title, " ====================\n", sep = "")
  cat("\nVIF pruning\n-----------\n")
  cat("Started with: ", paste(vif_report$start, collapse=", "), "\n", sep="")
  if (length(vif_report$removed)) {
    cat("Removed by VIF (>", vif_thresh, "): ", paste(vif_report$removed, collapse=", "), "\n", sep="")
  } else cat("No removals by VIF.\n")
  
  cat("\nFULL (ML) formula\n-----------------\n", deparse(full_formula), "\n", sep="")
  cat("\nFULL vs REDUCED (ML)\n---------------------\n")
  print(AIC(full_ml, red_ml))
  suppressWarnings(print(anova(full_ml, red_ml, test="Chisq")))
  
  cat("\nREDUCED (final REML) summary\n----------------------------\n")
  print(summary(final_reml))
  
  cat("\nTerm table (final REML)\n-----------------------\n")
  print(term_table(final_reml))
}

# ======================================================================
# Data prep from Modeling_data
# ======================================================================
prepare_block <- function(Modeling_data, response, vars_common, id = "LID"){
  vars <- c(response, vars_common)
  out <- Modeling_data %>%
    select(all_of(vars)) %>%
    mutate(
      "{id}" := factor(.data[[id]]),
      across(any_of(fish_vars), ~ as.numeric(.x))
    )
  out
}

# ======================================================================
# Core selection routine (presence-only or component of hurdle)
# ======================================================================
run_component_selection <- function(data, response, family,
                                    include_re = TRUE,
                                    p_cutoff = 0.05,
                                    vif_threshold = 5,
                                    protect_vars = character(0)) {
  
  id <- "LID"
  # 1) VIF pruning on fixed effects only
  pred_fixed <- setdiff(names(data), c(response, id))
  if (identical(family$family, "binomial")) {
    vif_mod <- remove_high_vif_glmBI(data, response, pred_fixed, threshold = vif_threshold,
                                     protect_vars = protect_vars)
  } else {
    vif_mod <- remove_high_vif_glmmTMB(data, response, pred_fixed, threshold = vif_threshold,
                                       protect_vars = protect_vars)
  }
  vif_report <- attr(vif_mod, "vif_report")
  kept_fixed <- setdiff(names(if (identical(family$family, "binomial")) coef(vif_mod) else glmmTMB::fixef(vif_mod)$cond),
                        "(Intercept)")
  # 2) Build FULL (ML)
  vars_step <- c(kept_fixed, if (include_re) id)
  rhs_full  <- make_full_rhs(vars_step, data, include_re = include_re, custom_k = custom_k)
  form_full <- as.formula(paste(response, "~", rhs_full))
  full_ml   <- mgcv::gam(form_full, data = data, family = family, method = "ML", select = TRUE)
  
  # 3) REDUCED via stepwise (ML)
  red_ml <- stepwise_gam_vars_with_LID(
    response = response, vars = vars_step, data = data,
    family = family, method = "ML", p_cutoff = p_cutoff,
    custom_k = custom_k, mandatory_vars = protect_vars
  )
  # 4) Final REML refit
  final_reml <- mgcv::gam(formula(red_ml), data = data, family = family, method = "REML", select = TRUE)
  
  list(vif_report = vif_report,
       full_formula = form_full,
       full_ml = full_ml,
       reduced_ml = red_ml,
       final_reml = final_reml)
}

# Build a stacked fish panel (same behaviour as CPUE section)
build_fish_panel <- function(model, data, family = c("binomial","gamma"),
                             fish_vars = c("Presence_Goldfish","Presence_Eel","Presence_Catfish","Presence_Common_smelt"),
                             only_retained = TRUE,  # if TRUE, show placeholder for dropped vars
                             hold_binaries = "mean") {
  family <- match.arg(family)
  
  # local wrapper using your robust plot function
  plot_one <- function(v){
    if (only_retained) {
      # show placeholder when var not in model
      plot_binary_single_component(
        model, v, data, family = family,
        exclude_RE = TRUE, hold_binaries = hold_binaries,
        fish_vars = fish_vars
      )
    } else {
      # force plotting even if not retained (will still show placeholder text if not in model)
      plot_binary_single_component(
        model, v, data, family = family,
        exclude_RE = TRUE, hold_binaries = hold_binaries,
        fish_vars = fish_vars
      )
    }
  }
  
  plots <- lapply(fish_vars, plot_one)
  if (!length(plots)) return(ggplot2::ggplot() + ggplot2::theme_void())
  patchwork::wrap_plots(plots, ncol = 1)
}

run_component_selection <- function(data, response, family,
                                    include_re = TRUE,
                                    p_cutoff = 0.05,
                                    vif_threshold = 5,
                                    protect_vars = character(0),
                                    exclude_predictors = character(0)) {
  
  id <- "LID"
  
  # Auto-exclude base variable if response is *_pos (e.g., CPUE_pos -> drop CPUE_Kōura/CPUE)
  base_from_pos <- character(0)
  if (grepl("_pos$", response)) {
    base_from_pos <- sub("_pos$", "", response)
    # also consider typical naming with macron or without
    # if a column exactly matches base_from_pos or common variants, exclude them
    variants <- c(base_from_pos,
                  gsub("ō", "o", base_from_pos),                # macron → plain o
                  paste0(base_from_pos, "_Kōura"),              # CPUE_pos -> CPUE_Kōura
                  paste0(gsub("ō", "o", base_from_pos), "_Koura"))
    variants <- intersect(variants, names(data))
    exclude_predictors <- unique(c(exclude_predictors, variants))
  }
  
  # 1) VIF pruning on fixed effects only
  pred_fixed <- setdiff(names(data), c(response, id, exclude_predictors))
  if (identical(family$family, "binomial")) {
    vif_mod <- remove_high_vif_glmBI(data, response, pred_fixed, threshold = vif_threshold,
                                     protect_vars = protect_vars)
  } else {
    vif_mod <- remove_high_vif_glmmTMB(data, response, pred_fixed, threshold = vif_threshold,
                                       protect_vars = protect_vars)
  }
  vif_report <- attr(vif_mod, "vif_report")
  kept_fixed <- setdiff(names(if (identical(family$family, "binomial")) coef(vif_mod) else glmmTMB::fixef(vif_mod)$cond),
                        "(Intercept)")
  
  # 2) Build FULL (ML) from VIF-kept variables (+ optional RE)
  vars_step <- c(kept_fixed, if (include_re) id)
  rhs_full  <- make_full_rhs(vars_step, data, include_re = include_re, custom_k = custom_k)
  form_full <- as.formula(paste(response, "~", rhs_full))
  full_ml   <- mgcv::gam(form_full, data = data, family = family, method = "ML", select = TRUE)
  
  # 3) REDUCED via stepwise (ML)
  red_ml <- stepwise_gam_vars_with_LID(
    response = response, vars = vars_step, data = data,
    family = family, method = "ML", p_cutoff = p_cutoff,
    custom_k = custom_k, mandatory_vars = protect_vars
  )
  
  # 4) Final REML refit
  final_reml <- mgcv::gam(formula(red_ml), data = data, family = family, method = "REML", select = TRUE)
  
  list(vif_report = vif_report,
       full_formula = form_full,
       full_ml = full_ml,
       reduced_ml = red_ml,
       final_reml = final_reml)
}

# Safer: no self-referential defaults, explicit local name
plot_binary_single_component <- function(model, var, data,
                                         family = c("binomial","gamma"),
                                         S = 2000, alpha = 0.05,
                                         exclude_RE = TRUE,
                                         hold_binaries = c("mean","zero","one"),
                                         id = "LID",
                                         fish_covars = c("Presence_Goldfish","Presence_Eel","Presence_Catfish","Presence_Common_smelt"),
                                         seed = 1) {
  family <- match.arg(family)
  hold_binaries <- match.arg(hold_binaries)
  
  # Helper to check presence of a term in a model
  has_term <- function(mod, v){
    tl <- attr(terms(mod), "term.labels")
    any(tl == v) || any(grepl(paste0("^s\\(", v, "(,|\\))"), tl))
  }
  
  # If variable isn't in the final model, return a labelled placeholder
  if (!has_term(model, var)) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5,
                          label = paste(var, "not retained in model"), size = 4) +
        ggplot2::xlim(0,1) + ggplot2::ylim(0,1) + ggplot2::theme_void()
    )
  }
  
  # Base at numeric means (+ default id level)
  num_means <- data %>% dplyr::summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))
  base <- as.list(num_means)
  if (id %in% names(data) && is.factor(data[[id]]))
    base[[id]] <- factor(levels(data[[id]])[1], levels = levels(data[[id]]))
  
  # Hold other binaries at mean/zero/one
  bin_fish <- intersect(fish_covars, names(data))
  others <- setdiff(bin_fish, var)
  set_bin <- function(x){
    if (hold_binaries == "mean") round(mean(data[[x]], na.rm = TRUE))
    else if (hold_binaries == "zero") 0L else 1L
  }
  for (x in others) base[[x]] <- set_bin(x)
  
  # newdata for 0/1, replicated across id
  mk_nd <- function(level){
    nd <- as.data.frame(base, stringsAsFactors = FALSE)
    nd[[var]] <- as.integer(level)
    if (id %in% names(data) && is.factor(data[[id]])) {
      nd <- do.call(rbind, lapply(levels(data[[id]]), function(lv){
        r <- nd; r[[id]] <- factor(lv, levels = levels(data[[id]])); r
      }))
    }
    # Align factor levels to model frame if available
    if (!is.null(model$model)) {
      common <- intersect(names(nd), names(model$model))
      for (nm in common) if (is.factor(model$model[[nm]]))
        nd[[nm]] <- factor(nd[[nm]], levels = levels(model$model[[nm]]))
    }
    nd
  }
  nd0 <- mk_nd(0L); nd1 <- mk_nd(1L)
  nd <- dplyr::bind_rows(dplyr::mutate(nd0, .level = 0L),
                         dplyr::mutate(nd1, .level = 1L))
  
  excl <- if (exclude_RE) {
    if (!length(model$smooth)) NULL else {
      has_re <- vapply(model$smooth, function(s) "LID" %in% s$term, logical(1))
      if (any(has_re)) "s(LID)" else NULL
    }
  } else NULL
  
  lp <- predict(model, newdata = nd, type = "link", se.fit = TRUE, exclude = excl)
  
  set.seed(seed)
  n <- nrow(nd)
  Z <- matrix(rnorm(n * S, lp$fit, pmax(lp$se.fit, .Machine$double.eps)), nrow = n, ncol = S)
  
  if (family == "binomial") {
    Y <- plogis(Z); ylab <- "Presence probability"
  } else {
    Y <- exp(Z);    ylab <- "Mean given presence (μ)"
  }
  
  sim_by_level <- function(level_flag){
    sims <- Y[nd$.level == level_flag, , drop = FALSE]
    colMeans(sims)
  }
  sim0 <- sim_by_level(0L); sim1 <- sim_by_level(1L)
  
  qlo <- alpha/2; qhi <- 1 - alpha/2
  df <- data.frame(
    level = factor(c(0,1), levels = c(0,1), labels = c("Absent","Present")),
    mean  = c(mean(sim0), mean(sim1)),
    lwr   = c(quantile(sim0, qlo), quantile(sim1, qlo)),
    upr   = c(quantile(sim0, qhi), quantile(sim1, qhi))
  )
  
  ggplot2::ggplot(df, ggplot2::aes(x = level, y = mean)) +
    ggplot2::geom_col(width = 0.6) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lwr, ymax = upr), width = 0.2) +
    ggplot2::labs(x = var, y = ylab)
}

# Consistent builder for stacked fish panels
build_fish_panel <- function(model, data,
                             family = c("binomial","gamma"),
                             fish_covars = c("Presence_Goldfish","Presence_Eel","Presence_Catfish","Presence_Common_smelt"),
                             hold_binaries = "mean") {
  family <- match.arg(family)
  plots <- lapply(fish_covars, function(v)
    plot_binary_single_component(model, v, data,
                                 family = family,
                                 hold_binaries = hold_binaries,
                                 fish_covars = fish_covars))
  if (!length(plots)) return(ggplot2::ggplot() + ggplot2::theme_void())
  patchwork::wrap_plots(plots, ncol = 1)
}


```

# Occupancy (presence/absence)
```{r fig.width=14, fig.height=5, dpi=100, echo=FALSE, message=FALSE, warning=FALSE}
PData <- prepare_block(Modeling_data, "Presence_Kōura", vars_common)

occ <- run_component_selection(
  data = PData, response = "Presence_Kōura",
  family = binomial(link = "logit"),
  include_re = INCLUDE_RE,
  p_cutoff = p_cutoff_ml,
  vif_threshold = vif_thresh,
  protect_vars = if (PROTECT_FISH) fish_vars else character(0))

print_selection_report(
  "OCCUPANCY (Presence/Absence)", occ$vif_report, occ$full_formula, occ$full_ml, occ$reduced_ml, occ$final_reml)

summary(occ$full_ml)
mp_final <- occ$final_reml
mgcv::gam.check(mp_final)
if (length(mp_final$smooth)) try(print(mgcv::concurvity(mp_final)), silent = TRUE)

# Occupancy plots
p_occ_smooths <- gratia::draw(mp_final)
p_occ_fish    <- build_fish_panel(
  model = mp_final, data = PData, family = "binomial",
  fish_vars = fish_vars, only_retained = TRUE, hold_binaries = "mean"
)

# Simple grouped CV by lake (5-fold) for calibration
set.seed(1)
k <- min(5, length(unique(PData$LID)))
fold_ids <- sample(rep(1:k, length.out = nrow(PData)))
pred_cv <- numeric(nrow(PData))
for (fold in 1:k) {
  idx <- which(fold_ids == fold)
  m <- mgcv::gam(formula(mp_final), data = PData[-idx, , drop = FALSE],
                 family = binomial(link="logit"), method = "REML", select = TRUE)
  pred_cv[idx] <- predict(m, newdata = PData[idx, , drop = FALSE], type = "response")
}
calib_data <- tibble(pred = pred_cv, obs = mp_final$y) %>%
  mutate(bin = ntile(pred, 5)) %>%
  group_by(bin) %>%
  summarise(mean_pred = mean(pred), obs_rate = mean(obs), n = n(), .groups = "drop")
Calibration_plotp <- ggplot(calib_data, aes(mean_pred, obs_rate)) +
  geom_point(size = 3) + geom_line() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "Mean predicted probability", y = "Observed proportion")

final_plot_occupancy <- patchwork::wrap_plots(
  list(as_plot(p_occ_smooths), as_plot(p_occ_fish), as_plot(Calibration_plotp)), ncol = 3
)
final_plot_occupancy
ggsave(file.path(fig_dir, "Occupancy_panels.png"), final_plot_occupancy, width = 14, height = 6, dpi = 300)


```

# CPUE — hurdle
```{r fig.width=14, fig.height=5, dpi=100, echo=FALSE, message=FALSE, warning=FALSE}
CData <- prepare_block(Modeling_data, "CPUE_Kōura", vars_common) %>%
  dplyr::mutate(CPUE_pos = as.integer(CPUE_Kōura > 0))


# Presence component

cpue_pres <- run_component_selection(
  data = CData,
  response = "CPUE_pos",
  family   = binomial(link = "logit"),
  include_re = INCLUDE_RE,
  p_cutoff = p_cutoff_ml,
  vif_threshold = vif_thresh,
  protect_vars = if (PROTECT_FISH) fish_vars else character(0),
  exclude_predictors = c("CPUE_Kōura","CPUE_Koura","CPUE")  # safety
)


summary(cpue_pres$full_ml)
print_selection_report( "CPUE — Presence (CPUE>0)",  cpue_pres$vif_report, cpue_pres$full_formula, cpue_pres$full_ml, cpue_pres$reduced_ml, cpue_pres$final_reml)

# Positives component
Cpos <- filter(CData, CPUE_pos == 1L)
cpue_pos <- run_component_selection(
  data = Cpos,
  response = "CPUE_Kōura",
  family   = Gamma(link = "log"),
  include_re = INCLUDE_RE,
  p_cutoff = p_cutoff_ml,
  vif_threshold = vif_thresh,
  protect_vars = if (PROTECT_FISH) fish_vars else character(0)
)

print_selection_report(
  "CPUE — Positives (Gamma log)",
  cpue_pos$vif_report, cpue_pos$full_formula, cpue_pos$full_ml, cpue_pos$reduced_ml, cpue_pos$final_reml
)

# Panels
p_smooth_pres <- gratia::draw(cpue_pres$final_reml)
p_smooth_pos  <- gratia::draw(cpue_pos$final_reml)

fish_pres <- lapply(fish_vars, function(v) plot_binary_single_component(cpue_pres$final_reml, v, CData, family="binomial"))
fish_pos  <- lapply(fish_vars, function(v) plot_binary_single_component(cpue_pos$final_reml,  v, Cpos,  family="gamma"))
p_fish_pres <- if (length(fish_pres)) patchwork::wrap_plots(fish_pres, ncol = 1) else ggplot() + theme_void()
p_fish_pos  <- if (length(fish_pos))  patchwork::wrap_plots(fish_pos,  ncol = 1) else ggplot() + theme_void()

pred_pres <- predict(cpue_pres$final_reml, type = "response", exclude = exclude_if_RE(cpue_pres$final_reml))
obs_pres  <- CData$CPUE_pos
roc_out   <- diag_roc(obs_pres, pred_pres)

pred_pos <- predict(cpue_pos$final_reml, type = "response", exclude = exclude_if_RE(cpue_pos$final_reml))
obs_pos  <- cpue_pos$final_reml$y
scat_out <- diag_scatter(obs_pos, pred_pos, "Predicted (Gamma; positives only — CPUE)", "Observed CPUE (positives)")

CPUE_panel_presence  <- patchwork::wrap_plots(list(as_plot(p_smooth_pres), as_plot(p_fish_pres), as_plot(roc_out$plot)), ncol = 3)
CPUE_panel_positives <- patchwork::wrap_plots(list(as_plot(p_smooth_pos),  as_plot(p_fish_pos),  as_plot(scat_out$plot)), ncol = 3)

CPUE_panel_presence
CPUE_panel_positives

ggsave(file.path(fig_dir, "CPUE_panel_presence.png"),  CPUE_panel_presence,  width = 14, height = 5, dpi = 300)
ggsave(file.path(fig_dir, "CPUE_panel_positives.png"), CPUE_panel_positives, width = 14, height = 5, dpi = 300)

```

# BCUE — hurdle
```{r fig.width=14, fig.height=5, dpi=100, echo=FALSE, message=FALSE, warning=FALSE}
BData <- prepare_block(Modeling_data, "BCUE_Kōura", vars_common) %>%
  dplyr::mutate(BCUE_pos = as.integer(BCUE_Kōura > 0))

# Presence component
bcue_pres <- run_component_selection(
  data = BData,
  response = "BCUE_pos",
  family   = binomial(link = "logit"),
  include_re = INCLUDE_RE,
  p_cutoff = p_cutoff_ml,
  vif_threshold = vif_thresh,
  protect_vars = if (PROTECT_FISH) fish_vars else character(0),
  exclude_predictors = c("BCUE_Kōura","BCUE_Koura","BCUE"))

print_selection_report("BCUE — Presence (BCUE>0)",bcue_pres$vif_report, bcue_pres$full_formula, bcue_pres$full_ml, bcue_pres$reduced_ml, bcue_pres$final_reml)

# Positives component
Bpos <- filter(BData, BCUE_pos == 1L)
bcue_pos <- run_component_selection(
  data = Bpos,
  response = "BCUE_Kōura",
  family   = Gamma(link = "log"),
  include_re = INCLUDE_RE,
  p_cutoff = p_cutoff_ml,
  vif_threshold = vif_thresh,
  protect_vars = if (PROTECT_FISH) fish_vars else character(0))

print_selection_report("BCUE — Positives (Gamma log)",bcue_pos$vif_report, bcue_pos$full_formula, bcue_pos$full_ml, bcue_pos$reduced_ml, bcue_pos$final_reml)

# Panels
p_smooth_pres_b <- gratia::draw(bcue_pres$final_reml)
p_smooth_pos_b  <- gratia::draw(bcue_pos$final_reml)

fish_covars <- c("Presence_Goldfish","Presence_Eel","Presence_Catfish","Presence_Common_smelt")

fish_pres_b <- lapply(fish_covars,function(v) plot_binary_single_component(bcue_pres$final_reml, v, BData, family = "binomial",fish_covars = fish_covars))
p_fish_pres_b <- patchwork::wrap_plots(fish_pres_b, ncol = 1)
fish_pos_b <- lapply(fish_covars,function(v) plot_binary_single_component(bcue_pos$final_reml, v, Bpos,family = "gamma",fish_covars = fish_covars))
p_fish_pos_b <- patchwork::wrap_plots(fish_pos_b, ncol = 1)

pred_pres_b <- predict(bcue_pres$final_reml, type = "response", exclude = exclude_if_RE(bcue_pres$final_reml))
obs_pres_b  <- BData$BCUE_pos
roc_out_b   <- diag_roc(obs_pres_b, pred_pres_b)

pred_pos_b <- predict(bcue_pos$final_reml, type = "response", exclude = exclude_if_RE(bcue_pos$final_reml))
obs_pos_b  <- bcue_pos$final_reml$y
scat_out_b <- diag_scatter(obs_pos_b, pred_pos_b, "Predicted (Gamma; positives only — BCUE)", "Observed BCUE (positives)")

BCUE_panel_presence  <- patchwork::wrap_plots(list(as_plot(p_smooth_pres_b), as_plot(p_fish_pres_b), as_plot(roc_out_b$plot)), ncol = 3)
BCUE_panel_positives <- patchwork::wrap_plots(list(as_plot(p_smooth_pos_b),  as_plot(p_fish_pos_b),  as_plot(scat_out_b$plot)), ncol = 3)

BCUE_panel_presence
BCUE_panel_positives

ggsave(file.path(fig_dir, "BCUE_panel_presence.png"),  BCUE_panel_presence,  width = 14, height = 5, dpi = 300)
ggsave(file.path(fig_dir, "BCUE_panel_positives.png"), BCUE_panel_positives, width = 14, height = 5, dpi = 300)

```





























# GAM model helpers
```{r}
# ---------- Helpers (drop-in) ----------
# Basic
is_cont <- function(x) is.numeric(x) && dplyr::n_distinct(x, na.rm = TRUE) >= 5
logit2prob <- function(x) exp(x) / (1 + exp(x))

# NEW builder: dataset-agnostic, caps k safely
build_smooth_gam <- function(var, data, custom_k = list(), include_re_for_LID = FALSE) {
  if (identical(var, "LID") && include_re_for_LID) return("s(LID, bs='re')")
  x <- data[[var]]
  if (is.numeric(x)) {
    nuniq <- dplyr::n_distinct(x, na.rm = TRUE)
    if (nuniq < 5) return(var)               # too few uniques → parametric
    k_req <- if (!is.null(custom_k[[var]])) custom_k[[var]] else 10
    k_cap <- max(3, min(k_req, nuniq - 1))   # cap k safely
    return(paste0("s(", var, ", bs='ts', k=", k_cap, ")"))
  }
  var
}

# RHS builder that uses a passed-in builder function
mk_rhs <- function(vars, data, custom_k = list(), allow_re = TRUE, builder = build_smooth_gam) {
  paste(
    vapply(vars, function(v) builder(v, data = data, custom_k = custom_k, include_re_for_LID = allow_re),
           character(1)),
    collapse = " + "
  )
}

# Model fit wrappers
fit_binom <- function(formula, data) {
  mgcv::gam(formula, data = data, family = binomial(link = "logit"), method = "REML", select = TRUE)
}
fit_gamma <- function(formula, data) {
  mgcv::gam(formula, data = data, family = Gamma(link = "log"), method = "REML", select = TRUE)
}

# ----- VIF pruning -----

# Bias-reduced logistic (binomial) VIF pruning (handles separation better)
remove_high_vif_glmBI <- function(data, response, predictors, threshold = 5) {
  if (!requireNamespace("brglm2", quietly = TRUE)) install.packages("brglm2")
  # Step 1: remove near-zero variance predictors up front (robustness)
  nzv <- caret::nearZeroVar(data[, predictors, drop = FALSE])
  if (length(nzv)) predictors <- predictors[-nzv]
  
  removed <- character()
  repeat {
    if (!length(predictors)) stop("All predictors removed.")
    fml <- as.formula(paste(response, "~", paste(predictors, collapse = " + ")))
    model <- try(glm(fml, data = data, family = binomial(link = "logit"), method = "brglmFit"), silent = TRUE)
    if (inherits(model, "try-error")) break
    vif_data <- performance::check_collinearity(model)
    vif_data <- vif_data[!grepl("\\|", vif_data$Term), ]
    if (!nrow(vif_data) || all(vif_data$VIF < threshold)) break
    to_remove <- vif_data$Term[which.max(vif_data$VIF)]
    predictors <- setdiff(predictors, to_remove)
    removed <- c(removed, to_remove)
  }
  model
}

# Tweedie GLMM (via glmmTMB) VIF pruning for continuous responses
remove_high_vif_glmmTMB <- function(data, response, predictors, threshold = 5) {
  # Step 1: remove near-zero variance predictors up front (robustness)
  nzv <- caret::nearZeroVar(data[, predictors, drop = FALSE])
  if (length(nzv)) predictors <- predictors[-nzv]
  
  removed <- character()
  repeat {
    if (!length(predictors)) stop("All predictors removed.")
    fml <- as.formula(paste(response, "~", paste(predictors, collapse = " + ")))
    model <- try(glmmTMB::glmmTMB(fml, data = data, family = glmmTMB::tweedie()), silent = TRUE)
    if (inherits(model, "try-error")) break
    vif_data <- performance::check_collinearity(model)
    vif_data <- vif_data[!grepl("\\|", vif_data$Term), ]
    if (!nrow(vif_data) || all(vif_data$VIF < threshold)) break
    to_remove <- vif_data$Term[which.max(vif_data$VIF)]
    predictors <- setdiff(predictors, to_remove)
    removed <- c(removed, to_remove)
  }
  model
}

# ----- Stepwise GAM with optional LID random effect via build_smooth('LID') -----
# Stepwise GAM that uses a passed-in builder function
stepwise_gam_vars_with_LID <- function(response, vars, data, family,
                                       method = "ML", p_cutoff = 0.1,
                                       custom_k = list(),
                                       builder = build_smooth_gam) {
  remaining_vars <- vars
  m <- NULL
  repeat {
    if (!length(remaining_vars)) { message("No variables left → stop."); break }
    
    rhs_terms <- vapply(
      remaining_vars,
      function(v) builder(v, data = data, custom_k = custom_k, include_re_for_LID = TRUE),
      character(1)
    )
    current_formula <- as.formula(paste(response, "~", paste(rhs_terms, collapse = " + ")))
    m <- mgcv::gam(current_formula, data = data, family = family, method = method, select = TRUE)
    sm <- summary(m)
    
    # collect p-values
    param_p <- if (!is.null(sm$p.table) && nrow(sm$p.table) > 0) {
      pcol <- intersect(colnames(sm$p.table), c("Pr(>|t|)","Pr(>|z|)"))[1]
      pp <- sm$p.table[, pcol]; names(pp) <- rownames(sm$p.table); pp[names(pp) != "(Intercept)"]
    } else numeric(0)
    
    smooth_p <- if (!is.null(sm$s.table)) {
      p <- sm$s.table[, "p-value"]; names(p) <- rownames(sm$s.table); p
    } else numeric(0)
    
    all_p <- c(param_p, smooth_p)
    drop <- all_p[all_p > p_cutoff]
    if (!length(drop)) { message("All terms have p <= ", p_cutoff, " → stop."); break }
    
    term_to_remove <- names(which.max(drop))
    var_to_remove <- if (grepl("^s\\(", term_to_remove)) sub("^s\\(([^,]+).*\\)$", "\\1", term_to_remove) else term_to_remove
    if (var_to_remove %in% remaining_vars) {
      remaining_vars <- setdiff(remaining_vars, var_to_remove)
    } else {
      message("Warning: variable '", var_to_remove, "' not in current set. Stopping."); break
    }
  }
  m
}

# ----- Significance helpers for plotting -----

get_sig_terms <- function(m, alpha = 0.05){
  sm <- summary(m)
  sig_smooths <- if (!is.null(sm$s.table)) {
    rn <- rownames(sm$s.table); rn[ sm$s.table[, "p-value"] < alpha ]
  } else character(0)
  sig_params  <- if (!is.null(sm$p.table) && nrow(sm$p.table) > 0) {
    pcol <- intersect(colnames(sm$p.table), c("Pr(>|t|)","Pr(>|z|)"))[1]
    rn <- rownames(sm$p.table); rn[rn != "(Intercept)" & sm$p.table[, pcol] < alpha]
  } else character(0)
  list(smooths = sig_smooths, params = sig_params)
}

sig_smooth_labels <- function(m, alpha = 0.05){
  sm <- summary(m)
  if (is.null(sm$s.table)) return(character(0))
  labs <- rownames(sm$s.table)
  labs[ sm$s.table[, "p-value"] < alpha ]
}

draw_sig_smooths <- function(m, alpha = 0.05){
  labs <- sig_smooth_labels(m, alpha)
  if (!length(labs)) return(NULL)
  gratia::draw(m, select = labs)
}

# ----- Prediction utilities -----

# Build a base grid at means (and all LIDs to average over RE via exclude)
make_base_grid <- function(df, req){
  out <- list()
  for (v in req) {
    if (!v %in% names(df)) next
    x <- df[[v]]
    if (is.numeric(x)) out[[v]] <- mean(x, na.rm = TRUE)
    else if (is.factor(x)) out[[v]] <- factor(levels(x)[1], levels = levels(x))
    else out[[v]] <- 0
  }
  if ("LID" %in% names(df)) out$LID <- levels(df$LID)
  as.data.frame(out)
}

# Exclude s(LID) only if present
exclude_if_RE <- function(mod){
  if (!length(mod$smooth)) return(NULL)
  has_re <- vapply(mod$smooth, function(s) "LID" %in% s$term, logical(1))
  if (any(has_re)) "s(LID)" else NULL
}

# Overall (p * mu) bar plot for a binary predictor
plot_overall_binary <- function(var, pres_model, pos_model, df, ylab = "Overall expected value"){
  req_pres <- attr(terms(pres_model), "term.labels")
  req_pos  <- attr(terms(pos_model),  "term.labels")
  req_all  <- union(req_pres, req_pos)
  base <- make_base_grid(df, req_all)
  
  new0 <- new1 <- base
  if (!var %in% names(new0)) new0[[var]] <- 0
  if (!var %in% names(new1)) new1[[var]] <- 1
  
  p0  <- predict(pres_model, newdata = new0, type = "response", exclude = exclude_if_RE(pres_model))
  p1  <- predict(pres_model, newdata = new1, type = "response", exclude = exclude_if_RE(pres_model))
  mu0 <- predict(pos_model,  newdata = new0, type = "response", exclude = exclude_if_RE(pos_model))
  mu1 <- predict(pos_model,  newdata = new1, type = "response", exclude = exclude_if_RE(pos_model))
  
  dfp <- tibble::tibble(level = factor(c(0,1)), overall = c(mean(p0*mu0, na.rm=TRUE), mean(p1*mu1, na.rm=TRUE)))
  ggplot2::ggplot(dfp, ggplot2::aes(level, overall)) +
    ggplot2::geom_col(width = 0.6, fill = "grey70") +
    ggplot2::labs(x = var, y = ylab)
}

# ----- Reporting -----

# Print selected formula and summary of the (REML) refit; returns the refit model
report_selected_model <- function(step_obj, data, family, fitter, title){
  cat("\n--- ", title, " ---\n", sep = "")
  cat("Selected formula:\n  ", deparse(formula(step_obj)), "\n\n", sep = "")
  final <- fitter(formula(step_obj), data)
  print(summary(final))
  invisible(final)
}

# Concurvity only if there are smooths
safe_concurvity <- function(m){
  if (length(m$smooth)) {
    try(print(mgcv::concurvity(m)), silent = TRUE)
  } else {
    message("No smooths → concurvity not applicable.")
  }
}


# ---- Hurdle-aware binary plot with error bars ----
plot_overall_binary_hurdle <- function(var, pres_model, pos_model, data,
                                       S = 2000, alpha = 0.05,
                                       exclude_RE = TRUE,
                                       hold_binaries = c("mean","zero","one")) {
  hold_binaries <- match.arg(hold_binaries)
  
  # 1) Identify predictors in models
  all_terms <- unique(names(model.frame(pres_model)))
  # Keep only columns present in data
  all_terms <- intersect(all_terms, names(data))
  
  # 2) Build a base row with numeric means; binaries handled below
  num_means <- data |>
    dplyr::summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))
  
  base <- as.list(num_means)
  # ensure non-numeric columns that models may expect
  if ("LID" %in% names(data)) base$LID <- factor(levels(data$LID)[1], levels = levels(data$LID))
  
  # 3) Set other binary fish covariates
  bin_fish <- c("Presence_Goldfish","Presence_Eel","Presence_Catfish","Presence_Common_smelt")
  bin_fish <- intersect(bin_fish, names(data))
  others <- setdiff(bin_fish, var)
  set_bin <- function(x){
    if (hold_binaries == "mean") round(mean(data[[x]], na.rm = TRUE))
    else if (hold_binaries == "zero") 0L else 1L
  }
  for (x in others) base[[x]] <- set_bin(x)
  
  # 4) Build newdata for var = 0 and 1 (duplicate across all LID levels if present)
  mk_nd <- function(level){
    nd <- as.data.frame(base, stringsAsFactors = FALSE)
    nd[[var]] <- level
    # replicate across LID levels (doesn't change predictions if we exclude RE)
    if ("LID" %in% names(data)) {
      nd <- do.call(rbind, lapply(levels(data$LID), function(lv){
        row <- nd
        row$LID <- factor(lv, levels = levels(data$LID))
        row
      }))
    }
    # Keep only columns the models can use
    keep <- unique(c(names(pres_model$model), names(pos_model$model)))
    keep <- intersect(keep, names(nd))
    nd[, keep, drop = FALSE]
  }
  
  nd0 <- mk_nd(0L)
  nd1 <- mk_nd(1L)
  nd  <- rbind(nd0, nd1)
  nd$.level <- rep(c(0L,1L), each = nrow(nd) / 2)
  
  # 5) Predictions on LINK scales with SEs
  excl_pres <- if (exclude_RE) exclude_if_RE(pres_model) else NULL
  excl_pos  <- if (exclude_RE) exclude_if_RE(pos_model)  else NULL
  
  lp_pres <- predict(pres_model, newdata = nd, type = "link", se.fit = TRUE, exclude = excl_pres)
  lp_pos  <- predict(pos_model,  newdata = nd, type = "link", se.fit = TRUE, exclude = excl_pos)
  
  # 6) Parametric sims on link scales → response
  set.seed(1)
  n  <- nrow(nd)
  Zp <- matrix(rnorm(n * S, mean = lp_pres$fit, sd = lp_pres$se.fit), nrow = n, ncol = S)
  Zm <- matrix(rnorm(n * S, mean = lp_pos$fit,  sd = lp_pos$se.fit),  nrow = n, ncol = S)
  
  P  <- plogis(Zp)     # presence prob
  MU <- exp(Zm)        # positives mean
  OVERALL <- P * MU    # unconditional mean = p * mu
  
  # 7) Aggregate by binary level (average over LID rows)
  idx0 <- nd$.level == 0L
  idx1 <- nd$.level == 1L
  sim0 <- colMeans(OVERALL[idx0, , drop = FALSE])
  sim1 <- colMeans(OVERALL[idx1, , drop = FALSE])
  
  qlo <- alpha/2; qhi <- 1 - alpha/2
  df <- data.frame(
    level = factor(c(0,1), levels = c(0,1), labels = c("Absent","Present")),
    mean  = c(mean(sim0), mean(sim1)),
    lwr   = c(quantile(sim0, qlo), quantile(sim1, qlo)),
    upr   = c(quantile(sim0, qhi), quantile(sim1, qhi))
  )
  
  # 8) Plot
  ggplot(df, aes(x = level, y = mean)) +
    geom_col(width = 0.6, col="grey70") +
    geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2) 
}

plot_binary_single_component <- function(model, var, data,
                                         family = c("binomial","gamma"),
                                         S = 2000, alpha = 0.05,
                                         exclude_RE = TRUE,
                                         hold_binaries = c("mean","zero","one")) {
  family <- match.arg(family)
  hold_binaries <- match.arg(hold_binaries)
  
  # Base numeric means, plus a default LID
  num_means <- data %>% dplyr::summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))
  base <- as.list(num_means)
  if ("LID" %in% names(data)) base$LID <- factor(levels(data$LID)[1], levels = levels(data$LID))
  
  # Hold other fish binaries
  bin_fish <- intersect(fish_vars, names(data))
  others <- setdiff(bin_fish, var)
  set_bin <- function(x){
    if (hold_binaries == "mean") round(mean(data[[x]], na.rm = TRUE))
    else if (hold_binaries == "zero") 0L else 1L
  }
  for (x in others) base[[x]] <- set_bin(x)
  
  # Make newdata for var = 0/1, optionally replicate across LID levels
  mk_nd <- function(level){
    nd <- as.data.frame(base, stringsAsFactors = FALSE)
    nd[[var]] <- as.integer(level)
    if ("LID" %in% names(data)) {
      nd <- do.call(rbind, lapply(levels(data$LID), function(lv){
        r <- nd; r$LID <- factor(lv, levels = levels(data$LID)); r
      }))
    }
    keep <- intersect(names(model$model), names(nd))
    nd[, keep, drop = FALSE]
  }
  nd0 <- mk_nd(0L); nd1 <- mk_nd(1L)
  nd <- dplyr::bind_rows(
    dplyr::mutate(nd0, .level = 0L),
    dplyr::mutate(nd1, .level = 1L)
  )
  
  excl <- if (exclude_RE) exclude_if_RE(model) else NULL
  
  # Link-scale sims → response
  lp <- predict(model, newdata = nd, type = "link", se.fit = TRUE, exclude = excl)
  set.seed(1)
  n <- nrow(nd)
  Z <- matrix(rnorm(n * S, lp$fit, lp$se.fit), nrow = n, ncol = S)
  
  if (family == "binomial") {
    Y <- plogis(Z); ylab <- "Presence probability"
  } else {
    Y <- exp(Z);    ylab <- "Mean given presence (μ)"
  }
  
  sim_by_level <- function(level_flag){
    sims <- Y[nd$.level == level_flag, , drop = FALSE]
    colMeans(sims)  # average over LID rows → population-level
  }
  sim0 <- sim_by_level(0L)
  sim1 <- sim_by_level(1L)
  
  qlo <- alpha/2; qhi <- 1 - alpha/2
  df <- data.frame(
    level = factor(c(0,1), levels = c(0,1), labels = c("Absent","Present")),
    mean  = c(mean(sim0), mean(sim1)),
    lwr   = c(quantile(sim0, qlo), quantile(sim1, qlo)),
    upr   = c(quantile(sim0, qhi), quantile(sim1, qhi))
  )
  
  ggplot(df, aes(x = level, y = mean)) +
    geom_col(width = 0.6) +
    geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2) +
    labs(x = var, y = ylab) 
}

```

# Kōura occupancy (presence/absence)
```{r fig.width=14, fig.height=5, dpi=100, echo=FALSE, message=FALSE, warning=FALSE}
custom_k <- list(Slope_5m=10, Riparian_vegetation=7, Overhanging_trees=5, Wood_cover=10,
                 Substrate_index=10, Temperature=10, pH=10, DO_mgl=10,
                 Emergent_Native=9, Submerged_Native=4, Submerged_Non_Native=9, Turf_Native=6)

## 1) Prep: select vars and coerce types
vars <- c(
  "Presence_Kōura","LID",
  "Slope_5m","Riparian_vegetation","Overhanging_trees","Wood_cover","Substrate_index",
  "Temperature","pH","DO_mgl","DO_percent","Specific_conductivity",
  "Emergent_Native","Submerged_Non_Native","Turf_Native","Submerged_Native","Emergent_Non_Native",
  "Presence_Common_smelt","Presence_Goldfish","Presence_Eel","Presence_Catfish")

PData <- Modeling_data |>
  dplyr::select(all_of(vars)) |>
  mutate(LID = factor(LID),
         across(c(Presence_Kōura, Presence_Common_smelt, Presence_Goldfish,
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
  method = "ML", p_cutoff = 0.05,custom_k = custom_k)

mp.10 <- stepwise_gam_vars_with_LID(
  response = response, vars = vars_for_stepwise,
  data = PData, family = binomial(link = "logit"),
  method = "ML", p_cutoff = 0.10,custom_k = custom_k)


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
    Presence_Goldfish   = 0,
    Presence_Eel        = 0,
    Presence_Catfish    = 0
  )
  
  # duplicate rows: absent (0) and present (1) for this fish
  newdat0 <- newdat
  newdat1 <- newdat
  newdat0[[fish_var]] <- "Absent"
  newdat1[[fish_var]] <- "Present"
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
         aes(x = factor(.data[[fish_var]], labels = c("Absent","Present")),
             y = prob)) +
    geom_col(width = 0.5, fill = "grey50") +
    geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2) +
    ylim(0,1) +
    labs(x = label, y = "Predicted kōura presence") 
}

# Now call it:
p_mori <- make_fish_plot("Presence_Goldfish", "Presence Goldfish", mp_final, PData) +labs(y = NULL)
p_eel  <- make_fish_plot("Presence_Eel", "Presence Eel", mp_final, PData)
p_cat  <- make_fish_plot("Presence_Catfish", "Presence Catfish", mp_final, PData) +labs(y = NULL)

p1 <- draw(mp_final)
p2 <- p_mori/ p_eel/ p_cat


# Build new data frame 
newdata_param <- expand.grid(
  Presence_Goldfish = c(0, 1),
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
p2 <- ggplot(newdata_param, aes(x = interaction(Presence_Goldfish, Presence_Eel, Presence_Catfish), y = predicted_prob)) +
  geom_col(width = 0.5, fill = "grey70") +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2) +
  labs(x = "Presence Goldfish and Eel", y = "Predicted Kōura Presence") 

# Stack the fish-effect plots first
p2 <- p_mori / p_eel / p_cat + plot_layout(tag_level = "new")


# Get predicted probabilities
#  Helpers (append) 
logit2prob <- function(x) exp(x)/(1+exp(x))

# Fit function wrapper: same formula and settings you used for mp_final
fit_gam_bin <- function(formula, data) {
  mgcv::gam(formula, data = data, family = binomial(link="logit"),
            method = "REML", select = TRUE)
}

#  Leave-one-site-out (within-lake)
cv_preds_site <- function(formula, data, id_col = "RowID") {
  df <- data
  if (!id_col %in% names(df)) df[[id_col]] <- seq_len(nrow(df))
  out <- vector("numeric", nrow(df))
  for (i in seq_len(nrow(df))) {
    train <- df[-i, , drop = FALSE]
    test  <- df[i,  , drop = FALSE]
    m <- fit_gam_bin(formula, train)
    out[i] <- predict(m, newdata = test, type = "response") # includes s(LID)
  }
  out
}

# Using your final occupancy model formula
form_final <- formula(mp_final)

#pred_in_sample <- predict(mp_final, type="response")  
pred_cv_site   <- cv_preds_site(form_final, PData)     # within-lake CV

# Equal-frequency binning (5 bins)
n_bins <- 5
calib_plot_data <- tibble(
  pred = pred_cv_site,
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



# Combine all three horizontally and tag A–C
final_plot_1 <- p1 | p2 | Calibration_plotp 

final_plot_1

ggsave(file.path(fig_dir, "Precence_Predictions.png"), final_plot_1, width = 14, height = 5, dpi = 300)
```

# CPUE — Hurdle GAM Script
```{r fig.width=14, fig.height=5, dpi=100, echo=FALSE, message=FALSE, warning=FALSE}
vars_c <- c(
  "CPUE_Kōura","LID",
  "Slope_5m","Riparian_vegetation","Overhanging_trees","Wood_cover","Substrate_index",
  "Temperature","pH","DO_mgl","DO_percent","Specific_conductivity",
  "Emergent_Native","Submerged_Non_Native","Turf_Native","Submerged_Native","Emergent_Non_Native",
  "Presence_Common_smelt","Presence_Goldfish","Presence_Eel","Presence_Catfish")



CData <- Modeling_data |>
  dplyr::select(all_of(vars_c)) |>
  dplyr::mutate(
    LID = factor(LID),
    across(c(Presence_Common_smelt, Presence_Goldfish, Presence_Eel, Presence_Catfish), ~ as.numeric(.)),
    CPUE_pos = as.integer(CPUE_Kōura > 0))

resp_posflag <- "CPUE_pos"


# 1) VIF pruning (exclude response & LID)
pred_fixed_c <- setdiff(names(CData), c("CPUE_Kōura", resp_posflag, "LID"))
vif_c <- remove_high_vif_glmmTMB(CData, "CPUE_Kōura", pred_fixed_c, threshold = 5)
kept_fixed_c <- setdiff(names(fixef(vif_c)$cond), "(Intercept)")
vars_for_stepwise_c <- c(kept_fixed_c, "LID")

# Utilities for formula building
is_cont <- function(x) is.numeric(x) && length(unique(na.omit(x))) > 5
k_for <- function(v) if (exists("custom_k") && !is.null(custom_k[[v]])) custom_k[[v]] else 10


# 2) FULL presence (binomial) formula
preds_pres <- setdiff(vars_for_stepwise_c, "LID")
sm_pres <- preds_pres[sapply(CData[preds_pres], is_cont)]
fx_pres <- setdiff(preds_pres, sm_pres)
re_term <- if ("LID" %in% vars_for_stepwise_c) "s(LID, bs='re')" else NULL
sm_terms_pres <- if (length(sm_pres)) paste0("s(", sm_pres, ", bs='ts', k=", sapply(sm_pres, k_for), ")") else NULL
full_terms_pres <- c(sm_terms_pres, fx_pres, re_term)
full_formula_pres <- as.formula(paste(resp_posflag, "~", paste(full_terms_pres, collapse = " + ")))


# 3) Stepwise presence (ML) → compare FULL vs REDUCED (ML) → Refit final (REML)
mc.pres <- stepwise_gam_vars_with_LID(
  response = resp_posflag, vars = vars_for_stepwise_c, data = CData,
  family = binomial(link = "logit"), method = "ML", p_cutoff = 0.05, custom_k = custom_k)

mc.pres_final <- report_selected_model(mc.pres, CData, binomial(link="logit"), fit_binom, "CPUE presence (CPUE>0) — final")

gam.check(mc.pres_final); safe_concurvity(mc.pres_final)

full_pres_ML <- mgcv::gam(full_formula_pres, data = CData, family = binomial(link="logit"), method = "ML")
red_pres_ML  <- mgcv::gam(formula(mc.pres_final), data = CData, family = binomial(link="logit"), method = "ML")
anova(full_pres_ML, red_pres_ML, test = "Chisq")
AIC(full_pres_ML, red_pres_ML)

mc.pres_final_REML <- mgcv::gam(formula(mc.pres_final), data = CData, family = binomial(link = "logit"), method = "REML")
summary(mc.pres_final_REML)


# 4) Positives subset + FULL positives (Gamma/log) formula
Cpos <- dplyr::filter(CData, CPUE_pos == 1)

preds_pos <- setdiff(vars_for_stepwise_c, "LID")
sm_pos <- preds_pos[sapply(Cpos[preds_pos], is_cont)]
fx_pos <- setdiff(preds_pos, sm_pos)
re_term_pos  <- if ("LID" %in% vars_for_stepwise_c) "s(LID, bs='re')" else NULL
sm_terms_pos <- if (length(sm_pos)) paste0("s(", sm_pos, ", bs='ts', k=", sapply(sm_pos, k_for), ")") else NULL
full_terms_pos <- c(sm_terms_pos, fx_pos, re_term_pos)
full_formula_pos <- as.formula(paste("CPUE_Kōura ~", paste(full_terms_pos, collapse = " + ")))


# 5) Stepwise positives (ML) → compare FULL vs REDUCED (ML) → Refit final (REML)
mc.pos <- stepwise_gam_vars_with_LID(
  response = "CPUE_Kōura", vars = vars_for_stepwise_c, data = Cpos,
  family = Gamma(link = "log"), method = "ML", p_cutoff = 0.05, custom_k = custom_k)

mc.pos_final <- report_selected_model(mc.pos, Cpos, Gamma(link="log"), fit_gamma, "CPUE positives (Gamma log) — final")

gam.check(mc.pos_final); safe_concurvity(mc.pos_final)

full_pos_ML <- mgcv::gam(full_formula_pos, data = Cpos, family = Gamma(link="log"), method = "ML")
red_pos_ML  <- mgcv::gam(formula(mc.pos_final), data = Cpos, family = Gamma(link="log"), method = "ML")
anova(full_pos_ML, red_pos_ML, test = "Chisq")
AIC(full_pos_ML, red_pos_ML)

mc.pos_final_REML <- mgcv::gam(formula(mc.pos_final), data = Cpos, family = Gamma(link = "log"), method = "REML")
summary(mc.pos_final_REML)


# 6) Overall expected CPUE = p * μ
new_cpue <- expand.grid(
  Presence_Goldfish = 0:1,
  Presence_Eel      = 0:1,
  Presence_Catfish  = 0:1,
  Presence_Common_smelt = 0:1,
  LID = levels(CData$LID),
  Riparian_vegetation   = mean(CData$Riparian_vegetation,   na.rm = TRUE),
  Substrate_index       = mean(CData$Substrate_index,       na.rm = TRUE),
  Temperature           = mean(CData$Temperature,           na.rm = TRUE),
  pH                    = mean(CData$pH,                    na.rm = TRUE),
  DO_mgl                = mean(CData$DO_mgl,                na.rm = TRUE),
  Specific_conductivity = mean(CData$Specific_conductivity, na.rm = TRUE),
  Slope_5m              = mean(CData$Slope_5m,              na.rm = TRUE),
  Overhanging_trees     = mean(CData$Overhanging_trees,     na.rm = TRUE),
  Wood_cover            = mean(CData$Wood_cover,            na.rm = TRUE),
  Emergent_Native       = mean(CData$Emergent_Native,       na.rm = TRUE),
  Submerged_Native      = mean(CData$Submerged_Native,      na.rm = TRUE),
  Submerged_Non_Native  = mean(CData$Submerged_Non_Native,  na.rm = TRUE),
  Turf_Native           = mean(CData$Turf_Native,           na.rm = TRUE))

p_link <- predict(
  mc.pres_final_REML, newdata = new_cpue, type = "response", se.fit = FALSE,
  exclude = exclude_if_RE(mc.pres_final_REML))

mu_pos <- predict(
  mc.pos_final_REML, newdata = new_cpue, type = "response", se.fit = FALSE,
  exclude = exclude_if_RE(mc.pos_final_REML))

new_cpue$p_pos <- p_link
new_cpue$mu_pos <- mu_pos
new_cpue$overall_mean <- new_cpue$p_pos * new_cpue$mu_pos

eff_cpue <- new_cpue |>
  dplyr::group_by(Presence_Goldfish, Presence_Eel, Presence_Catfish, Presence_Common_smelt) |>
  dplyr::summarise(overall = mean(overall_mean), .groups = "drop")


# 7) Positives-only fit metrics (REML objects)
pred_pos <- predict(mc.pos_final_REML, type = "response")
obs_pos  <- mc.pos_final_REML$y
R2_pos   <- cor(pred_pos, obs_pos)^2
RMSE_pos <- sqrt(mean((pred_pos - obs_pos)^2))


# 8) Plots 
alpha <- 0.05

# --- Smooths (significant only) ---
p_smooth_pres <- draw_sig_smooths(mc.pres_final_REML, alpha)   # presence (binomial)
p_smooth_pos  <- draw_sig_smooths(mc.pos_final_REML,  alpha)   # positives (Gamma)

# --- Fish binaries per component (95% CI bars) ---
fish_vars <- c("Presence_Goldfish","Presence_Eel","Presence_Catfish","Presence_Common_smelt")
sig_terms_pres <- get_sig_terms(mc.pres_final_REML, alpha)$params
sig_terms_pos  <- get_sig_terms(mc.pos_final_REML,  alpha)$params

fish_pres <- intersect(fish_vars, sig_terms_pres)  # or set fish_pres <- fish_vars to force all
fish_pos  <- intersect(fish_vars, sig_terms_pos)   # or set fish_pos  <- fish_vars to force all

# Build fish plots
plots_fish_pres <- lapply(fish_pres, function(v)
  plot_binary_single_component(mc.pres_final_REML, v, CData,
                               family = "binomial", S = 2000, alpha = 0.05,
                               exclude_RE = TRUE, hold_binaries = "mean"))
plots_fish_pos <- lapply(fish_pos, function(v)
  plot_binary_single_component(mc.pos_final_REML, v, Cpos,
                               family = "gamma", S = 2000, alpha = 0.05,
                               exclude_RE = TRUE, hold_binaries = "mean"))

# If none were significant, show a blank panel rather than NULL
p_fish_pres <- if (length(plots_fish_pres)) {patchwork::wrap_plots(plots_fish_pres, ncol = 1) } else { ggplot() + theme_void()}
p_fish_pos <- if (length(plots_fish_pos)) {patchwork::wrap_plots(plots_fish_pos, ncol = 1)} else { ggplot() + theme_void()}

# --- Diagnostics per component (kept inside the respective pane) ---

## Presence ROC
pred_pres <- predict(mc.pres_final_REML, type = "response", exclude = exclude_if_RE(mc.pres_final_REML))
obs_pres  <- CData$CPUE_pos
library(pROC)
roc_curve <- roc(obs_pres, pred_pres); auc_val <- auc(roc_curve)
roc_df <- data.frame(tpr = roc_curve$sensitivities, fpr = 1 - roc_curve$specificities)
p_roc <- ggplot(roc_df, aes(fpr, tpr)) +
  geom_path() + geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.6, y = 0.1, label = paste0("AUC = ", round(as.numeric(auc_val), 3))) +
  labs(x = "False positive rate", y = "True positive rate") 

## Positives observed vs predicted
pred_pos <- predict(mc.pos_final_REML, type = "response", exclude = exclude_if_RE(mc.pos_final_REML))
obs_pos  <- mc.pos_final_REML$y
R2_pos   <- cor(pred_pos, obs_pos)^2
RMSE_pos <- sqrt(mean((pred_pos - obs_pos)^2))
df_pred_obs <- data.frame(Predicted = pred_pos, Observed = obs_pos)
p_pos_scatter <- ggplot(df_pred_obs, aes(Predicted, Observed)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text",
           x = max(df_pred_obs$Predicted, na.rm = TRUE) * 0.7,
           y = max(df_pred_obs$Observed,  na.rm = TRUE) * 0.8,
           label = paste0("R² = ", round(R2_pos, 2), "\nRMSE = ", round(RMSE_pos, 2))) +
  labs(x = "Predicted (Gamma; positives only)", y = "Observed CPUE (positives)") 

# --- Two panes exactly: Presence and Positives ---
as_plot <- function(p) {
  if (is.null(p)) return(plot_spacer())
  if (inherits(p, c("gg", "ggplot", "patchwork"))) return(p)
  if (is.list(p)) {
    if (length(p) == 0) return(plot_spacer())
    return(wrap_plots(p))
  }
  plot_spacer()
}

# placeholder with a message
placeholder <- function(txt) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.5, y = 0.5, label = txt, size = 5) +
    ggplot2::xlim(0,1) + ggplot2::ylim(0,1) + ggplot2::theme_void()
}

#BCUE_panel_presence
p_smooth_pres_n  <- as_plot(p_smooth_pres)
p_fish_pres_n  <- as_plot(p_fish_pres)
p_roc_n <- as_plot(p_roc)
CPUE_panel_presence <- wrap_plots(list(p_smooth_pres_n, p_fish_pres_n, p_roc_n), ncol = 3)

# BCUE_panel_positives
p_smooth_pos_n  <- as_plot(p_smooth_pos)
p_fish_pos_n  <- as_plot(p_fish_pos)
p_pos_scatter_n <- as_plot(p_pos_scatter)
CPUE_panel_positives <- wrap_plots(list(p_smooth_pos, p_pos_scatter_n), ncol = 2)

# Show them
CPUE_panel_presence
CPUE_panel_positives

# Save 
ggsave(file.path(fig_dir, "CPUE_panel_presence.png"),  CPUE_panel_presence,  width = 14, height = 5, dpi = 300)
ggsave(file.path(fig_dir, "CPUE_panel_positives.png"), CPUE_panel_positives, width = 10, height = 5, dpi = 300)
```

# BCUE — Hurdle GAM Script
```{r fig.width=14, fig.height=5, dpi=100, echo=FALSE, message=FALSE, warning=FALSE}
vars_b <- c(
  "BCUE_Kōura","LID",
  "Slope_5m","Riparian_vegetation","Overhanging_trees","Wood_cover","Substrate_index",
  "Temperature","pH","DO_mgl","DO_percent","Specific_conductivity",
  "Emergent_Native","Submerged_Non_Native","Turf_Native","Submerged_Native","Emergent_Non_Native",
  "Presence_Common_smelt","Presence_Goldfish","Presence_Eel","Presence_Catfish")

BData <- Modeling_data |>
  dplyr::select(all_of(vars_b)) |>
  dplyr::mutate(
    LID = factor(LID),
    across(c(Presence_Common_smelt, Presence_Goldfish, Presence_Eel, Presence_Catfish), ~ as.numeric(.)),
    BCUE_pos = as.integer(BCUE_Kōura > 0))

resp_posflag <- "BCUE_pos"


# 1) VIF pruning (exclude response & LID)
pred_fixed_b <- setdiff(names(BData), c("BCUE_Kōura", resp_posflag, "LID"))
vif_b <- remove_high_vif_glmmTMB(BData, "BCUE_Kōura", pred_fixed_b, threshold = 5)
kept_fixed_b <- setdiff(names(fixef(vif_b)$cond), "(Intercept)")
vars_for_stepwise_b <- c(kept_fixed_b, "LID")

# Utilities for formula building
is_cont <- function(x) is.numeric(x) && length(unique(na.omit(x))) > 5
k_for <- function(v) if (exists("custom_k") && !is.null(custom_k[[v]])) custom_k[[v]] else 10

# 2) FULL presence (binomial) formula
preds_pres <- setdiff(vars_for_stepwise_b, "LID")
sm_pres <- preds_pres[sapply(BData[preds_pres], is_cont)]
fx_pres <- setdiff(preds_pres, sm_pres)
re_term <- if ("LID" %in% vars_for_stepwise_b) "s(LID, bs='re')" else NULL
sm_terms_pres <- if (length(sm_pres)) paste0("s(", sm_pres, ", bs='ts', k=", sapply(sm_pres, k_for), ")") else NULL
full_terms_pres <- c(sm_terms_pres, fx_pres, re_term)
full_formula_pres <- as.formula(paste(resp_posflag, "~", paste(full_terms_pres, collapse = " + ")))

# 3) Stepwise presence (ML) → compare FULL vs REDUCED (ML) → Refit final (REML)
mb.pres <- stepwise_gam_vars_with_LID(
  response = resp_posflag, vars = vars_for_stepwise_b, data = BData,
  family = binomial(link = "logit"), method = "ML", p_cutoff = 0.05, custom_k = custom_k)

mb.pres_final <- report_selected_model(mb.pres, BData, binomial(link="logit"), fit_binom, "BCUE presence (BCUE>0) — final")

gam.check(mb.pres_final); safe_concurvity(mb.pres_final)


full_pres_ML <- mgcv::gam(full_formula_pres, data = BData, family = binomial(link="logit"), method = "ML")
red_pres_ML  <- mgcv::gam(formula(mb.pres_final), data = BData, family = binomial(link="logit"), method = "ML")
summary(full_pres_ML)
summary(red_pres_ML)
anova(full_pres_ML, red_pres_ML, test = "Chisq")
AIC(full_pres_ML, red_pres_ML)


mb.pres_final_REML <- mgcv::gam(formula(mb.pres_final), data = BData, family = binomial(link = "logit"), method = "REML")
summary(mb.pres_final_REML)


# 4) Positives subset + FULL positives (Gamma/log) formula
Bpos <- dplyr::filter(BData, BCUE_pos == 1)

preds_pos <- setdiff(vars_for_stepwise_b, "LID")
sm_pos <- preds_pos[sapply(Bpos[preds_pos], is_cont)]
fx_pos <- setdiff(preds_pos, sm_pos)
re_term_pos  <- if ("LID" %in% vars_for_stepwise_b) "s(LID, bs='re')" else NULL
sm_terms_pos <- if (length(sm_pos)) paste0("s(", sm_pos, ", bs='ts', k=", sapply(sm_pos, k_for), ")") else NULL
full_terms_pos <- c(sm_terms_pos, fx_pos, re_term_pos)
full_formula_pos <- as.formula(paste("BCUE_Kōura ~", paste(full_terms_pos, collapse = " + ")))


# 5) Stepwise positives (ML) → compare FULL vs REDUCED (ML) → Refit final (REML)
mb.pos <- stepwise_gam_vars_with_LID(
  response = "BCUE_Kōura", vars = vars_for_stepwise_c, data = Bpos,
  family = Gamma(link = "log"), method = "ML", p_cutoff = 0.05, custom_k = custom_k)

mb.pos_final <- report_selected_model(mb.pos, Bpos, Gamma(link="log"), fit_gamma, "BCUE positives (Gamma log) — final")

gam.check(mb.pos_final); safe_concurvity(mb.pos_final)

full_pos_ML <- mgcv::gam(full_formula_pos, data = Bpos, family = Gamma(link="log"), method = "ML")
red_pos_ML  <- mgcv::gam(formula(mb.pos_final), data = Bpos, family = Gamma(link="log"), method = "ML")
summary(full_pos_ML)
summary(red_pos_ML)
anova(full_pos_ML, red_pos_ML, test = "Chisq")
AIC(full_pos_ML, red_pos_ML)

mb.pos_final_REML <- mgcv::gam(formula(mb.pos_final), data = Bpos, family = Gamma(link = "log"), method = "REML")
summary(mb.pos_final_REML)


# 6) Overall expected CPUE = p * μ
new_bcue <- expand.grid(
  Presence_Goldfish = 0:1,
  Presence_Eel      = 0:1,
  Presence_Catfish  = 0:1,
  Presence_Common_smelt = 0:1,
  LID = levels(BData$LID),
  Riparian_vegetation   = mean(BData$Riparian_vegetation,   na.rm = TRUE),
  Substrate_index       = mean(BData$Substrate_index,       na.rm = TRUE),
  Temperature           = mean(BData$Temperature,           na.rm = TRUE),
  pH                    = mean(BData$pH,                    na.rm = TRUE),
  DO_mgl                = mean(BData$DO_mgl,                na.rm = TRUE),
  Specific_conductivity = mean(BData$Specific_conductivity, na.rm = TRUE),
  Slope_5m              = mean(BData$Slope_5m,              na.rm = TRUE),
  Overhanging_trees     = mean(BData$Overhanging_trees,     na.rm = TRUE),
  Wood_cover            = mean(BData$Wood_cover,            na.rm = TRUE),
  Emergent_Native       = mean(BData$Emergent_Native,       na.rm = TRUE),
  Submerged_Native      = mean(BData$Submerged_Native,      na.rm = TRUE),
  Submerged_Non_Native  = mean(BData$Submerged_Non_Native,  na.rm = TRUE),
  Turf_Native           = mean(BData$Turf_Native,           na.rm = TRUE))

p_link <- predict(
  mb.pres_final_REML, newdata = new_bcue, type = "response", se.fit = FALSE,
  exclude = exclude_if_RE(mb.pres_final_REML))

mu_pos <- predict(
  mb.pos_final_REML, newdata = new_bcue, type = "response", se.fit = FALSE,
  exclude = exclude_if_RE(mb.pos_final_REML))

new_bcue$p_pos <- p_link
new_bcue$mu_pos <- mu_pos
new_bcue$overall_mean <- new_bcue$p_pos * new_bcue$mu_pos

eff_cpue <- new_bcue |>
  dplyr::group_by(Presence_Goldfish, Presence_Eel, Presence_Catfish, Presence_Common_smelt) |>
  dplyr::summarise(overall = mean(overall_mean), .groups = "drop")


# 7) Positives-only fit metrics (REML objects)
pred_pos <- predict(mb.pos_final_REML, type = "response")
obs_pos  <- mb.pos_final_REML$y
R2_pos   <- cor(pred_pos, obs_pos)^2
RMSE_pos <- sqrt(mean((pred_pos - obs_pos)^2))


# 8) Plots 
alpha <- 0.05

# --- Smooths (significant only) ---
p_smooth_pres <- draw_sig_smooths(mb.pres_final_REML, alpha)   # presence (binomial)
p_smooth_pos  <- draw_sig_smooths(mb.pos_final_REML,  alpha)   # positives (Gamma)

# --- Fish binaries per component (95% CI bars) ---
fish_vars <- c("Presence_Goldfish","Presence_Eel","Presence_Catfish","Presence_Common_smelt")
sig_terms_pres <- get_sig_terms(mb.pres_final_REML, alpha)$params
sig_terms_pos  <- get_sig_terms(mb.pos_final_REML,  alpha)$params

fish_pres <- intersect(fish_vars, sig_terms_pres)  # or set fish_pres <- fish_vars to force all
fish_pos  <- intersect(fish_vars, sig_terms_pos)   # or set fish_pos  <- fish_vars to force all



# Build fish plots
plots_fish_pres <- lapply(fish_pres, function(v)
  plot_binary_single_component(mb.pres_final_REML, v, BData,
                               family = "binomial", S = 2000, alpha = 0.05,
                               exclude_RE = TRUE, hold_binaries = "mean")
)
plots_fish_pos <- lapply(fish_pos, function(v)
  plot_binary_single_component(mb.pos_final_REML, v, Bpos,
                               family = "gamma", S = 2000, alpha = 0.05,
                               exclude_RE = TRUE, hold_binaries = "mean")
)

# If none were significant, show a blank panel rather than NULL
p_fish_pres <- if (length(plots_fish_pres)) {patchwork::wrap_plots(plots_fish_pres, ncol = 1) } else { ggplot() + theme_void()}
p_fish_pos <- if (length(plots_fish_pos)) {patchwork::wrap_plots(plots_fish_pos, ncol = 1)} else { ggplot() + theme_void()}

# --- Diagnostics per component (kept inside the respective pane) ---

## Presence ROC
pred_pres <- predict(mb.pres_final_REML, type = "response", exclude = exclude_if_RE(mb.pres_final_REML))
obs_pres  <- BData$BCUE_pos
library(pROC)
roc_curve <- roc(obs_pres, pred_pres); auc_val <- auc(roc_curve)
roc_df <- data.frame(tpr = roc_curve$sensitivities, fpr = 1 - roc_curve$specificities)
p_roc <- ggplot(roc_df, aes(fpr, tpr)) +
  geom_path() + geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.6, y = 0.1, label = paste0("AUC = ", round(as.numeric(auc_val), 3))) +
  labs(x = "False positive rate", y = "True positive rate") 

## Positives observed vs predicted
pred_pos <- predict(mb.pos_final_REML, type = "response", exclude = exclude_if_RE(mb.pos_final_REML))
obs_pos  <- mb.pos_final_REML$y
R2_pos   <- cor(pred_pos, obs_pos)^2
RMSE_pos <- sqrt(mean((pred_pos - obs_pos)^2))
df_pred_obs <- data.frame(Predicted = pred_pos, Observed = obs_pos)
p_pos_scatter <- ggplot(df_pred_obs, aes(Predicted, Observed)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text",
           x = max(df_pred_obs$Predicted, na.rm = TRUE) * 0.7,
           y = max(df_pred_obs$Observed,  na.rm = TRUE) * 0.8,
           label = paste0("R² = ", round(R2_pos, 2), "\nRMSE = ", round(RMSE_pos, 2))) +
  labs(x = "Predicted (Gamma; positives only)", y = "Observed BCUE (positives)") 

# --- Two panes exactly: Presence and Positives ---
as_plot <- function(p) {
  if (is.null(p)) return(plot_spacer())
  if (inherits(p, c("gg", "ggplot", "patchwork"))) return(p)
  if (is.list(p)) {
    if (length(p) == 0) return(plot_spacer())
    return(wrap_plots(p))
  }
  plot_spacer()
}

# (Optional) placeholder with a message
placeholder <- function(txt) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.5, y = 0.5, label = txt, size = 5) +
    ggplot2::xlim(0,1) + ggplot2::ylim(0,1) + ggplot2::theme_void()
}

#BCUE_panel_presence
p_smooth_pres_n  <- as_plot(p_smooth_pres)
p_fish_pres_n  <- as_plot(p_fish_pres)
p_roc_n <- as_plot(p_roc)
BCUE_panel_presence <- wrap_plots(list(p_smooth_pres_n, p_fish_pres_n, p_roc_n), ncol = 3)

# BCUE_panel_positives
p_smooth_pos_n  <- as_plot(p_smooth_pos)
p_fish_pos_n  <- as_plot(p_fish_pos)
p_pos_scatter_n <- as_plot(p_pos_scatter)
BCUE_panel_positives <- wrap_plots(list(p_fish_pos_n, p_pos_scatter_n), ncol = 2)

# Show them
BCUE_panel_presence
BCUE_panel_positives


# Save (optional)
ggsave(file.path(fig_dir, "BCUE_panel_presence.png"),  BCUE_panel_presence,  width = 14, height = 5, dpi = 300)
ggsave(file.path(fig_dir, "BCUE_panel_positives.png"), BCUE_panel_positives, width = 8, height = 5, dpi = 300)

```
# transform the Fish_data
```{r}
# Safe helpers
safe_min  <- function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
safe_max  <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
safe_mean <- function(x) { m <- mean(x, na.rm = TRUE); if (is.nan(m)) NA_real_ else m }

# Base fish data
fish <- Fish_data %>%
  filter(!is.na(Species)) %>%
  filter(!is.na(Monitoring_ID))

# 1) Deployment (effort) per site × net type
deploy <- fish %>%
  filter(!is.na(Net_type), !is.na(Amount_nets)) %>%
  distinct(Monitoring_ID, Net_type, Amount_nets) %>%
  group_by(Monitoring_ID, Net_type) %>%
  summarise(Effort = max(Amount_nets, na.rm = TRUE), .groups = "drop") %>%
  mutate(Effort = ifelse(is.finite(Effort), Effort, 0))

# 2) Catches per site × net type × species (totals)
catches <- fish %>%
  group_by(Monitoring_ID, Net_type, Species) %>%
  summarise(
    Total_Individuals = sum(Amount,   na.rm = TRUE),
    Total_Weight      = sum(Weight_g, na.rm = TRUE),
    .groups = "drop"
  )

# 3) Ensure zero-catch combinations exist
catches_full <- deploy %>%
  left_join(catches, by = c("Monitoring_ID", "Net_type")) %>%
  tidyr::complete(Monitoring_ID, Net_type, Species,
                  fill = list(Total_Individuals = 0, Total_Weight = 0)) %>%
  group_by(Monitoring_ID, Net_type) %>%
  mutate(Effort = dplyr::first(Effort)) %>%
  ungroup() %>%
  mutate(Effort = coalesce(Effort, 0))

# 4) Correct CPUE/BCUE per site × species (effort-weighted across nets)
CPUE_BCUE <- catches_full %>%
  group_by(Monitoring_ID, Species) %>%
  summarise(
    Total_Individuals = sum(Total_Individuals, na.rm = TRUE),
    Total_Weight      = sum(Total_Weight,      na.rm = TRUE),
    Total_Effort      = sum(Effort,            na.rm = TRUE),
    CPUE              = ifelse(Total_Effort > 0, Total_Individuals / Total_Effort, NA_real_),
    BCUE              = ifelse(Total_Effort > 0, Total_Weight      / Total_Effort, NA_real_),
    .groups = "drop"
  )

# 5) Size stats per site × species
size_stats <- fish %>%
  group_by(Monitoring_ID, Species) %>%
  summarise(
    Mean_Length = safe_mean(Length_mm),
    Min_Length  = safe_min(Length_mm),
    Max_Length  = safe_max(Length_mm),
    Mean_Weight = safe_mean(Weight_g),
    Min_Weight  = safe_min(Weight_g),
    Max_Weight  = safe_max(Weight_g),
    .groups = "drop"
  )

# 6) Presence/absence (+ Predator_Fish_Presence)
species_presence <- fish %>%
  distinct(Monitoring_ID, Species) %>%
  mutate(Presence = 1L) %>%
  pivot_wider(
    names_from  = Species,
    values_from = Presence,
    values_fill = list(Presence = 0L),
    names_prefix = "Presence_"
  )

for (pred in c("Trout","Eel","Catfish")) {
  nm <- paste0("Presence_", pred)
  if (!nm %in% names(species_presence)) species_presence[[nm]] <- 0L
}

species_presence <- species_presence %>%
  mutate(
    Predator_Fish_Presence = pmax(
      coalesce(Presence_Trout,   0L),
      coalesce(Presence_Eel,     0L),
      coalesce(Presence_Catfish, 0L)
    )
  )

# 7) LONG table with CPUE/BCUE + size stats
site_species <- CPUE_BCUE %>%
  left_join(size_stats, by = c("Monitoring_ID","Species"))

# 8) WIDE per site summary
species_wide <- site_species %>%
  select(Monitoring_ID, Species,
         Total_Individuals, Total_Weight, Total_Effort,
         CPUE, BCUE,
         Mean_Length, Min_Length, Max_Length,
         Mean_Weight, Min_Weight, Max_Weight) %>%
  pivot_wider(
    names_from  = Species,
    values_from = c(Total_Individuals, Total_Weight, Total_Effort, CPUE, BCUE,
                    Mean_Length, Min_Length, Max_Length, Mean_Weight, Min_Weight, Max_Weight),
    names_sep   = "_",
    values_fill = 0
  ) %>%
  left_join(species_presence, by = "Monitoring_ID") %>%
  mutate(
    Richness = rowSums(across(starts_with("Total_Individuals_"), ~ . > 0), na.rm = TRUE),
    Abundance = rowSums(across(starts_with("Total_Individuals_"), ~ .), na.rm = TRUE) -
      rowSums(across(matches("^Total_Individuals_(Bullies|Common_smelt)$"), ~ .), na.rm = TRUE)
  )

























# GAM Modeling setup & helpers
```{r}
Modeling_data <- Monitoring_CPUE_data %>% dplyr::filter(Monitoring == 0)

alpha_sig    <- 0.05      # significance threshold for showing fish plots
p_cutoff_ml  <- 0.05     # stepwise ML drop threshold
vif_thresh   <- 5
INCLUDE_RE   <- TRUE
PROTECT_FISH <- FALSE

custom_k <- list(
  Slope_5m=10, Riparian_vegetation=7, Overhanging_trees=5, Wood_cover=10,
  Substrate_index=10, Temperature=10, pH=10, DO_mgl=10,
  Emergent_Native=9, Submerged_Non_Native=9, Turf_Native=6, Submerged_Native=6)

fish_vars <- c("Presence_Goldfish","Presence_Eel","Presence_Catfish","Presence_Common_smelt")

vars_common <- c("LID",
                 "Slope_5m","Riparian_vegetation","Overhanging_trees","Wood_cover","Substrate_index",
                 "Temperature","pH","DO_mgl","DO_percent","Specific_conductivity",
                 "Emergent_Native","Submerged_Non_Native","Turf_Native","Submerged_Native","Emergent_Non_Native",
                 fish_vars)

# Helpers ---
is_cont <- function(x) is.numeric(x) && dplyr::n_distinct(x, na.rm = TRUE) >= 5

# Construct a smooth or parametric term string for a single variable
build_smooth_gam <- function(var, data, custom_k = list(), include_re_for_LID = FALSE) {
  if (identical(var, "LID") && include_re_for_LID) return("s(LID, bs='re')")
  x <- data[[var]]
  if (is.numeric(x)) {
    nuniq <- dplyr::n_distinct(x, na.rm = TRUE)
    if (nuniq < 5) return(var)  # treat binary/few-level numerics as parametric
    k_req <- if (!is.null(custom_k[[var]])) custom_k[[var]] else 10
    k_cap <- max(3, min(k_req, nuniq - 1))
    return(paste0("s(", var, ", bs='ts', k=", k_cap, ")"))
  }
  var
}

exclude_if_RE <- function(mod){
  if (!length(mod$smooth)) return(NULL)
  has_re <- vapply(mod$smooth, function(s) "LID" %in% s$term, logical(1))
  if (any(has_re)) "s(LID)" else NULL
}

# VIF pruning (binomial) using brglm2 to reduce separation issues
remove_high_vif_glmBI <- function(data, response, predictors, threshold = 5, protect_vars = character(0)) {
  nzv <- caret::nearZeroVar(data[, predictors, drop = FALSE])
  if (length(nzv)) predictors <- predictors[-nzv]
  
  rpt <- list(removed = character(), start = predictors)
  repeat {
    if (!length(predictors)) stop("All predictors removed during VIF pruning.")
    fml <- as.formula(paste(response, "~", paste(predictors, collapse = " + ")))
    model <- try(glm(fml, data = data, family = binomial(link = "logit"), method = "brglmFit"), silent = TRUE)
    if (inherits(model, "try-error")) break
    vif_data <- performance::check_collinearity(model)
    vif_data <- vif_data[!grepl("\\|", vif_data$Term), ]
    if (!nrow(vif_data) || all(vif_data$VIF < threshold)) break
    
    ord <- order(vif_data$VIF, decreasing = TRUE)
    to_remove <- NA_character_
    for (i in ord) {
      cand <- vif_data$Term[i]
      if (!(cand %in% protect_vars)) { to_remove <- cand; break }
    }
    if (is.na(to_remove)) break
    predictors <- setdiff(predictors, to_remove)
    rpt$removed <- c(rpt$removed, to_remove)
  }
  list(predictors = predictors, report = rpt)
}

# VIF pruning (Gamma/Tweedie proxy) via glmmTMB
remove_high_vif_glmmTMB <- function(data, response, predictors, threshold = 5, protect_vars = character(0)) {
  nzv <- caret::nearZeroVar(data[, predictors, drop = FALSE])
  if (length(nzv)) predictors <- predictors[-nzv]
  
  rpt <- list(removed = character(), start = predictors)
  repeat {
    if (!length(predictors)) stop("All predictors removed during VIF pruning.")
    fml <- as.formula(paste(response, "~", paste(predictors, collapse = " + ")))
    model <- try(glmmTMB::glmmTMB(fml, data = data, family = glmmTMB::tweedie()), silent = TRUE)
    if (inherits(model, "try-error")) break
    vif_data <- performance::check_collinearity(model)
    vif_data <- vif_data[!grepl("\\|", vif_data$Term), ]
    if (!nrow(vif_data) || all(vif_data$VIF < threshold)) break
    
    ord <- order(vif_data$VIF, decreasing = TRUE)
    to_remove <- NA_character_
    for (i in ord) {
      cand <- vif_data$Term[i]
      if (!(cand %in% protect_vars)) { to_remove <- cand; break }
    }
    if (is.na(to_remove)) break
    predictors <- setdiff(predictors, to_remove)
    rpt$removed <- c(rpt$removed, to_remove)
  }
  list(predictors = predictors, report = rpt)
}

# Extract p-value for either parametric var or smooth s(var)
get_term_p <- function(m, var, kind = c("any","param","smooth")){
  kind <- match.arg(kind)
  sm <- summary(m)
  
  if (kind %in% c("any","param")) {
    if (!is.null(sm$p.table) && nrow(sm$p.table) > 0) {
      pcol <- intersect(colnames(sm$p.table), c("Pr(>|t|)","Pr(>|z|)"))[1]
      if (!is.na(pcol) && var %in% rownames(sm$p.table)) {
        return(as.numeric(sm$p.table[var, pcol]))
      }
    }
  }
  
  if (kind %in% c("any","smooth")) {
    if (!is.null(sm$s.table) && nrow(sm$s.table) > 0) {
      sname <- paste0("s(", var, ")")
      if (sname %in% rownames(sm$s.table)) {
        return(as.numeric(sm$s.table[sname, "p-value"]))
      }
    }
  }
  
  NA_real_
}

# Plot a binary covariate effect (0 vs 1) with uncertainty, holding other fish at mean/zero/one
plot_binary_single_component <- function(model, var, data,
                                         family = c("binomial","gamma"),
                                         S = 2000, alpha = 0.05,
                                         exclude_RE = TRUE,
                                         hold_binaries = c("mean","zero","one"),
                                         id = "LID",
                                         fish_covars = c("Presence_Goldfish","Presence_Eel","Presence_Catfish","Presence_Common_smelt"),
                                         seed = 1) {
  family <- match.arg(family)
  hold_binaries <- match.arg(hold_binaries)
  
  # check term presence
  tl <- attr(terms(model), "term.labels")
  has_term <- any(tl == var) || any(grepl(paste0("^s\\(", var, "(,|\\))"), tl))
  if (!has_term) return(ggplot() + theme_void())
  
  # base at numeric means (+ default id level)
  num_means <- data %>% dplyr::summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))
  base <- as.list(num_means)
  if (id %in% names(data) && is.factor(data[[id]]))
    base[[id]] <- factor(levels(data[[id]])[1], levels = levels(data[[id]]))
  
  # hold other fish binaries
  others <- setdiff(intersect(fish_covars, names(data)), var)
  set_bin <- function(x){
    if (hold_binaries == "mean") round(mean(data[[x]], na.rm = TRUE))
    else if (hold_binaries == "zero") 0L else 1L
  }
  for (x in others) base[[x]] <- set_bin(x)
  
  mk_nd <- function(level){
    nd <- as.data.frame(base, stringsAsFactors = FALSE)
    nd[[var]] <- as.integer(level)
    if (id %in% names(data) && is.factor(data[[id]])) {
      nd <- do.call(rbind, lapply(levels(data[[id]]), function(lv){
        r <- nd; r[[id]] <- factor(lv, levels = levels(data[[id]])); r
      }))
    }
    if (!is.null(model$model)) {
      common <- intersect(names(nd), names(model$model))
      for (nm in common) if (is.factor(model$model[[nm]]))
        nd[[nm]] <- factor(nd[[nm]], levels = levels(model$model[[nm]]))
    }
    nd
  }
  nd0 <- mk_nd(0L); nd1 <- mk_nd(1L)
  nd  <- dplyr::bind_rows(dplyr::mutate(nd0, .level = 0L),
                          dplyr::mutate(nd1, .level = 1L))
  
  excl <- if (exclude_RE) exclude_if_RE(model) else NULL
  lp <- predict(model, newdata = nd, type = "link", se.fit = TRUE, exclude = excl)
  
  set.seed(seed)
  n <- nrow(nd)
  Z <- matrix(rnorm(n * S, lp$fit, pmax(lp$se.fit, .Machine$double.eps)), nrow = n, ncol = S)
  
  if (family == "binomial") {
    Y <- plogis(Z); ylab <- "Presence probability"
  } else {
    Y <- exp(Z);    ylab <- "Mean given presence (μ)"
  }
  
  sim_by_level <- function(level_flag){
    sims <- Y[nd$.level == level_flag, , drop = FALSE]
    colMeans(sims)
  }
  sim0 <- sim_by_level(0L); sim1 <- sim_by_level(1L)
  
  qlo <- alpha/2; qhi <- 1 - alpha/2
  df <- data.frame(
    level = factor(c(0,1), levels = c(0,1), labels = c("Absent","Present")),
    mean  = c(mean(sim0), mean(sim1)),
    lwr   = c(quantile(sim0, qlo), quantile(sim1, qlo)),
    upr   = c(quantile(sim0, qhi), quantile(sim1, qhi))
  )
  
  ggplot(df, aes(x = level, y = mean)) +
    geom_col(width = 0.6) +
    geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2) +
    labs(x = var, y = ylab)
}

as_plot <- function(p) if (inherits(p, c("gg","ggplot","patchwork"))) p else patchwork::plot_spacer()

# DATA PREP HELPERS

prepare_block <- function(Modeling_data, response, vars_common, id = "LID"){
  vars <- c(response, vars_common)
  out <- Modeling_data %>%
    dplyr::select(dplyr::all_of(vars)) %>%
    dplyr::mutate(
      "{id}" := factor(.data[[id]]),
      dplyr::across(dplyr::any_of(fish_vars), ~ as.numeric(.x))
    )
  out
}


# Row of neutral reference values
ref_row <- function(data){
  dplyr::summarise(
    data,
    dplyr::across(
      dplyr::everything(),
      \(x){
        if (is.numeric(x)) stats::median(x, na.rm = TRUE)
        else if (is.factor(x)) levels(x)[1L]
        else if (is.logical(x)) FALSE
        else if (is.character(x)) unique(stats::na.omit(x))[1L]
        else x[1L]
      }
    )
  )
}

is_binary_numeric <- function(x) is.numeric(x) && dplyr::n_distinct(x, na.rm = TRUE) == 2

# Get non-smooth model terms that exist in data (parametric)
get_parametric_terms <- function(model, data){
  tl <- attr(stats::terms(model), "term.labels") # e.g. "Presence_Eel", "s(pH)"
  # keep non-smooth, drop s(), te(), ti(), t2()
  param_labels <- tl[!grepl("^(s|te|ti|t2)\\(", tl)]
  # only those present in data
  param_labels[param_labels %in% names(data)]
}

# Plot parametric numeric (continuous) partial on response scale
plot_param_numeric <- function(model, var, data, n = 100){
  r <- range(data[[var]], na.rm = TRUE)
  grid <- data.frame(x = seq(r[1], r[2], length.out = n)); names(grid) <- var
  
  base <- ref_row(data)
  newd <- base[rep(1, nrow(grid)), , drop = FALSE]
  newd[[var]] <- grid[[var]]
  
  pr  <- predict(model, newdata = newd, se.fit = TRUE, type = "link")
  il  <- model$family$linkinv
  fit <- il(pr$fit); lo <- il(pr$fit - 1.96*pr$se.fit); hi <- il(pr$fit + 1.96*pr$se.fit)
  
  ggplot2::ggplot(data.frame(x = newd[[var]], fit, lo, hi),
                  ggplot2::aes(x = .data[[var]], y = fit)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), alpha = 0.2) +
    ggplot2::geom_line() +
    ggplot2::labs(x = var, y = "Partial effect (response scale)")
}

# Plot parametric binary/factor partial on response scale
plot_param_factor <- function(model, var, data){
  # Detect how the model was fit
  mod_is_factor <- !is.null(model$model) && is.factor(model$model[[var]])
  mod_is_numeric <- !is.null(model$model) && is.numeric(model$model[[var]])
  
  # Decide the levels we want to show
  if (mod_is_numeric || is_binary_numeric(data[[var]])) {
    # binary numeric 0/1
    lev_vals <- c(0, 1)              # numeric for newdata
    lev_labs <- c("0","1")           # labels for plotting
  } else if (is.factor(data[[var]])) {
    lev_vals <- levels(data[[var]])  # factor levels
    lev_labs <- lev_vals
  } else {
    # fallback: unique sorted values (rare for true parametric factors)
    lev_vals <- sort(unique(stats::na.omit(data[[var]])))
    lev_labs <- as.character(lev_vals)
  }
  
  # Build reference newdata row(s)
  base <- ref_row(data)
  newd <- base[rep(1, length(lev_vals)), , drop = FALSE]
  
  # Critically: keep the same TYPE as used in the model
  if (mod_is_numeric) {
    newd[[var]] <- as.numeric(lev_vals)
  } else {
    newd[[var]] <- factor(lev_vals, levels = lev_vals)
  }
  
  # Keep factor columns compatible with model frame
  if (!is.null(model$model)) {
    common <- intersect(names(newd), names(model$model))
    for (nm in common) if (is.factor(model$model[[nm]]))
      newd[[nm]] <- factor(newd[[nm]], levels = levels(model$model[[nm]]))
  }
  
  pr  <- predict(model, newdata = newd, se.fit = TRUE, type = "link")
  il  <- model$family$linkinv
  fit <- il(pr$fit); lo <- il(pr$fit - 1.96*pr$se.fit); hi <- il(pr$fit + 1.96*pr$se.fit)
  
  # For plotting we can show labels as factors, independent of newdata type
  x_for_plot <- factor(lev_labs, levels = lev_labs)
  
  ggplot2::ggplot(data.frame(x = x_for_plot, fit, lo, hi),
                  ggplot2::aes(x = x, y = fit)) +
    ggplot2::geom_point() +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lo, ymax = hi), width = 0.15) +
    ggplot2::labs(x = var, y = "Partial effect (response scale)")
}

# Dispatch: prefer gratia for smooths, otherwise parametric plot (factor vs numeric)
plot_single_component <- function(model, var, data){
  sm_lab <- paste0("s(", var, ")")
  if (sm_lab %in% gratia::smooths(model)) {
    return(gratia::draw(model, select = sm_lab))
  }
  if (is.factor(data[[var]]) || is_binary_numeric(data[[var]])) {
    return(plot_param_factor(model, var, data))
  }
  plot_param_numeric(model, var, data)
}

# Significant parametric terms, optionally excluding LID (random effect)
significant_parametric_terms <- function(model, data, alpha = 0.1, exclude = c("LID")){
  params <- setdiff(get_parametric_terms(model, data), exclude)
  keep <- vapply(params, \(v){
    p <- get_term_p(model, v)
    !is.na(p) && p <= alpha
  }, logical(1))
  params[keep]
}

plot_linear_single_component <- function(model, var, data,
                                         family = c("binomial","gamma"),
                                         n = 100, alpha = 0.05,
                                         exclude_RE = TRUE) {
  family <- match.arg(family)
  
  # base at numeric means + factor baselines
  num_means <- data %>% dplyr::summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))
  base <- as.list(num_means)
  
  # align factor levels
  if (!is.null(model$model)) {
    for (nm in names(model$model)) {
      if (is.factor(model$model[[nm]])) {
        base[[nm]] <- factor(levels(model$model[[nm]])[1], levels = levels(model$model[[nm]]))
      }
    }
  }
  
  # grid over the param var
  x <- data[[var]]
  xr <- range(x, na.rm = TRUE)
  grid <- seq(xr[1], xr[2], length.out = n)
  
  nd <- as.data.frame(base)
  nd <- nd[rep(1, n), , drop = FALSE]
  nd[[var]] <- grid
  
  # predict (population level if RE present)
  excl <- if (exclude_RE) exclude_if_RE(model) else NULL
  pr <- predict(model, newdata = nd, type = "link", se.fit = TRUE, exclude = excl)
  
  # inverse link
  if (family == "binomial") {
    fit <- plogis(pr$fit)
    lwr <- plogis(pr$fit - qnorm(1 - alpha/2) * pr$se.fit)
    upr <- plogis(pr$fit + qnorm(1 - alpha/2) * pr$se.fit)
    ylab <- "Response (probability)"
  } else {
    fit <- exp(pr$fit)
    lwr <- exp(pr$fit - qnorm(1 - alpha/2) * pr$se.fit)
    upr <- exp(pr$fit + qnorm(1 - alpha/2) * pr$se.fit)
    ylab <- "Response (μ)"
  }
  
  df <- data.frame(x = grid, fit = fit, lwr = lwr, upr = upr)
  ggplot(df, aes(x, fit)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
    geom_line() +
    labs(x = var, y = ylab)
}


# Build a panel of parametric plots (significant-only by default)
parametric_panel <- function(model, data, alpha = 0.1, exclude = character(0),
                             family = c("binomial","gamma"),
                             fish_covars = c("Presence_Goldfish","Presence_Eel","Presence_Catfish","Presence_Common_smelt"),
                             hold_binaries = "mean") {
  family <- match.arg(family)
  
  sm <- summary(model)
  if (is.null(sm$p.table) || nrow(sm$p.table) == 0) return(ggplot() + theme_void())
  
  pcol <- intersect(colnames(sm$p.table), c("Pr(>|t|)","Pr(>|z|)"))[1]
  if (is.na(pcol)) return(ggplot() + theme_void())
  
  # all parametric term names (no intercept)
  par_names <- setdiff(rownames(sm$p.table), "(Intercept)")
  
  # remove excluded and any that also have a smooth s(var)
  smooth_rows <- if (!is.null(sm$s.table)) rownames(sm$s.table) else character(0)
  smooth_vars <- sub("^s\\(([^,\\)]+).*$", "\\1", smooth_rows)
  par_names <- setdiff(par_names, union(exclude, smooth_vars))
  
  # keep only significant parametric terms
  par_sig <- par_names[sm$p.table[par_names, pcol, drop = TRUE] <= alpha]
  if (!length(par_sig)) return(ggplot() + theme_void())
  
  # plot each parametric term appropriately
  plots <- lapply(par_sig, function(v) {
    x <- data[[v]]
    if (is.null(x)) return(ggplot() + theme_void())
    
    # binary (0/1 numeric) or 2-level factor -> use binary plotter
    is_bin_num <- is.numeric(x) && dplyr::n_distinct(x, na.rm = TRUE) <= 2
    is_bin_fac <- is.factor(x)  && nlevels(x) == 2
    
    if (is_bin_num || is_bin_fac) {
      plot_binary_single_component(model, v, data,
                                   family = family,
                                   hold_binaries = hold_binaries,
                                   fish_covars = fish_covars)
    } else {
      # continuous linear parametric -> linear effect plot
      plot_linear_single_component(model, v, data, family = family)
    }
  })
  
  if (!length(plots)) return(ggplot() + theme_void())
  patchwork::wrap_plots(plots, ncol = 1)
}




```


# OCCUPANCY: Presence/Absence
```{r fig.width=14, fig.height=5, dpi=100, echo=FALSE, message=FALSE, warning=FALSE}
PData <- prepare_block(Modeling_data, "Presence_Kōura", vars_common)

# 2) Candidate predictors for VIF (fixed effects only)
id <- "LID"
pred_fixed_occ <- setdiff(names(PData), c("Presence_Kōura", id))

# 3) VIF pruning (binomial)
if (PROTECT_FISH) {
  vif_occ <- remove_high_vif_glmBI(PData, "Presence_Kōura", pred_fixed_occ, threshold = vif_thresh, protect_vars = fish_vars)
} else {
  vif_occ <- remove_high_vif_glmBI(PData, "Presence_Kōura", pred_fixed_occ, threshold = vif_thresh, protect_vars = character(0))
}
kept_fixed_occ <- vif_occ$predictors

# 4) Build FULL (ML) formula
vars_step_occ <- c(kept_fixed_occ, if (INCLUDE_RE) id)
rhs_full_occ <- paste(vapply(vars_step_occ, function(v) build_smooth_gam(v, PData, custom_k, include_re_for_LID = INCLUDE_RE), character(1)), collapse = " + ")
form_full_occ <- as.formula(paste("Presence_Kōura ~", rhs_full_occ))

# 5) Fit FULL (ML)
full_ml_occ <- mgcv::gam(form_full_occ, data = PData, family = binomial(link="logit"), method = "ML", select = TRUE)

# 6) Stepwise (ML) by p-values
remaining_occ <- vars_step_occ
repeat {
  # fit current
  rhs_now <- paste(vapply(remaining_occ, function(v) build_smooth_gam(v, PData, custom_k, include_re_for_LID = INCLUDE_RE), character(1)), collapse = " + ")
  m_now <- mgcv::gam(as.formula(paste("Presence_Kōura ~", rhs_now)), data = PData, family = binomial(link="logit"), method = "ML", select = TRUE)
  sm <- summary(m_now)
  
  # collect p-values
  ps <- c()
  if (!is.null(sm$p.table) && nrow(sm$p.table) > 0) {
    pcol <- intersect(colnames(sm$p.table), c("Pr(>|t|)","Pr(>|z|)"))[1]
    pvec <- sm$p.table[, pcol]
    ps <- c(ps, pvec[names(pvec) != "(Intercept)"])
  }
  if (!is.null(sm$s.table) && nrow(sm$s.table) > 0) {
    pvec <- sm$s.table[, "p-value"]; names(pvec) <- rownames(sm$s.table)
    ps <- c(ps, pvec)
  }
  drop_candidates <- ps[ps > p_cutoff_ml]
  
  if (!length(drop_candidates)) { red_ml_occ <- m_now; break }
  
  ordered <- names(sort(drop_candidates, decreasing = TRUE))
  # never drop protected fish if requested
  ordered_vars <- vapply(ordered, function(x) if (grepl("^s\\(", x)) sub("^s\\(([^,]+).*\\)$", "\\1", x) else x, character(1))
  if (PROTECT_FISH) {
    ordered_vars <- setdiff(ordered_vars, fish_vars)
    if (!length(ordered_vars)) { red_ml_occ <- m_now; break }
  }
  remove_v <- ordered_vars[1]
  
  message("Dropping variable: ", remove_v)
  
  remaining_occ <- setdiff(remaining_occ, remove_v)
  if (!length(remaining_occ)) { red_ml_occ <- m_now; break }
}

# ML comparison: FULL vs REDUCED
summary(full_ml_occ)
summary(red_ml_occ)
AIC(full_ml_occ, red_ml_occ)
suppressWarnings(print(anova(full_ml_occ, red_ml_occ, test = "Chisq")))

# test the effect of LID
m_noLID <- red_ml_occ
m_withLID <- mgcv::gam(formula(Presence_Kōura ~ s(Substrate_index, bs="ts", k=10) +
                                 s(Temperature, bs="ts", k=10) + s(LID, bs="re")),data = PData, family = binomial(link="logit"), method = "ML")

anova(m_noLID, m_withLID, test = "Chisq")


# 7) Final REML
final_occ <- mgcv::gam(formula(red_ml_occ), data = PData, family = binomial(link="logit"), method = "REML", select = TRUE)
summary(final_occ)
mgcv::concurvity(final_occ)
gam.check(final_occ)

# 8) Plots
p_occ_smooths <- gratia::draw(final_occ)

# plot All significant parametric terms
p_occ_param <- parametric_panel(final_occ, PData, alpha = alpha_sig,exclude = c("LID"),family = "binomial", fish_covars = fish_vars)

# ROC curve + AUC
pred_pres <- predict(final_occ, type = "response", exclude = exclude_if_RE(final_occ))
obs_pres  <- PData$Presence_Kōura == 1
roc_obj   <- pROC::roc(obs_pres, pred_pres)
auc_val   <- as.numeric(pROC::auc(roc_obj))
roc_df    <- data.frame(tpr = roc_obj$sensitivities, fpr = 1 - roc_obj$specificities)
p_pres_roc <- ggplot(roc_df, aes(fpr, tpr)) +
  geom_path() + geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.6, y = 0.1, label = paste0("AUC = ", round(auc_val, 3))) +
  labs(x = "False positive rate", y = "True positive rate")



# Simple grouped CV-by-lake for calibration
set.seed(1)
k_occ <- min(5, nlevels(PData$LID))
lid_levels <- levels(PData$LID)
lid_folds  <- sample(rep(1:k_occ, length.out = length(lid_levels)))
lid2fold   <- setNames(lid_folds, lid_levels)
fold_vec   <- unname(lid2fold[as.character(PData$LID)])
pred_cv_occ <- rep(NA_real_, nrow(PData))

# Terms present in the final model (strip s(...))
final_terms <- {
  tl <- attr(terms(final_occ), "term.labels")
  sub("^s\\(([^,\\)]+).*$", "\\1", tl)
}

for (fold in 1:k_occ) {
  idx   <- which(fold_vec == fold)
  train <- PData[-idx, , drop = FALSE]
  test  <- PData[idx,  , drop = FALSE]
  
  # Rebuild RHS for the TRAINING DATA so k is capped to n_unique-1 in each fold
  rhs_now <- paste(
    vapply(
      final_terms,
      function(v) build_smooth_gam(v, train, custom_k, include_re_for_LID = INCLUDE_RE),
      character(1)
    ),
    collapse = " + "
  )
  fml_now <- as.formula(paste("Presence_Kōura ~", rhs_now))
  
  m <- mgcv::gam(
    fml_now,
    data   = train,
    family = binomial(link = "logit"),
    method = "REML",
    select = TRUE
  )
  
  # Align factor levels in test to those learned in train
  for (nm in names(test)) {
    if (is.factor(train[[nm]])) {
      test[[nm]] <- factor(test[[nm]], levels = levels(train[[nm]]))
    }
  }
  
  # Population-level prediction (exclude random effect)
  pred_cv_occ[idx] <- predict(
    m,
    newdata = test,
    type    = "response",
    exclude = exclude_if_RE(m)
  )
}

# Build the calibration plot
calib_data_occ <- tibble(pred = pred_cv_occ, obs = PData$Presence_Kōura) %>%
  mutate(bin = ntile(pred, 5)) %>%
  group_by(bin) %>%
  summarise(mean_pred = mean(pred, na.rm = TRUE),
            obs_rate  = mean(obs,  na.rm = TRUE),
            n = dplyr::n(), .groups = "drop")

Calibration_plot_occ <- ggplot(calib_data_occ, aes(mean_pred, obs_rate)) +
  geom_point(size = 3) + geom_line() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "Mean predicted probability", y = "Observed proportion")

final_plot_occupancy <- patchwork::wrap_plots(
  list(as_plot(p_occ_smooths), as_plot(p_occ_param), as_plot(p_pres_roc), as_plot(Calibration_plot_occ)), ncol = 4)

final_plot_occupancy

ggsave(file.path(fig_dir, "Occupancy_panels_1.png"), final_plot_occupancy, width = 14, height = 6, dpi = 300)

```

# CPUE: HURDLE
```{r fig.width=14, fig.height=5, dpi=100, echo=FALSE, message=FALSE, warning=FALSE}
CData <- prepare_block(Modeling_data, "Weighted_CPUE_Kōura", vars_common) %>%
  dplyr::mutate(CPUE_pos = as.integer(Weighted_CPUE_Kōura > 0))

# ---------- Presence component (binomial)
# VIF candidates (exclude base CPUE columns)
pred_fixed_cpue_pres <- setdiff(names(CData), c("CPUE_pos","LID","Weighted_CPUE_Kōura","Weighted_CPUE_Koura","CPUE"))
if (PROTECT_FISH) {
  vif_cpue_pres <- remove_high_vif_glmBI(CData, "CPUE_pos", pred_fixed_cpue_pres, threshold = vif_thresh, protect_vars = fish_vars)
} else {
  vif_cpue_pres <- remove_high_vif_glmBI(CData, "CPUE_pos", pred_fixed_cpue_pres, threshold = vif_thresh, protect_vars = character(0))
}
kept_fixed_cpue_pres <- vif_cpue_pres$predictors

vars_step_cpue_pres <- c(kept_fixed_cpue_pres, if (INCLUDE_RE) "LID")
rhs_cpue_pres <- paste(vapply(vars_step_cpue_pres, function(v) build_smooth_gam(v, CData, custom_k, include_re_for_LID = INCLUDE_RE), character(1)), collapse = " + ")
full_ml_cpue_pres <- mgcv::gam(as.formula(paste("CPUE_pos ~", rhs_cpue_pres)),
                               data = CData, family = binomial(link="logit"), method = "ML", select = TRUE)

# Stepwise (ML)
remaining_cpue_pres <- vars_step_cpue_pres
repeat {
  rhs_now <- paste(vapply(remaining_cpue_pres, function(v) build_smooth_gam(v, CData, custom_k, include_re_for_LID = INCLUDE_RE), character(1)), collapse = " + ")
  m_now <- mgcv::gam(as.formula(paste("CPUE_pos ~", rhs_now)), data = CData, family = binomial(link="logit"), method = "ML", select = TRUE)
  sm <- summary(m_now)
  
  ps <- c()
  if (!is.null(sm$p.table) && nrow(sm$p.table) > 0) {
    pcol <- intersect(colnames(sm$p.table), c("Pr(>|t|)","Pr(>|z|)"))[1]
    pvec <- sm$p.table[, pcol]
    ps <- c(ps, pvec[names(pvec) != "(Intercept)"])
  }
  if (!is.null(sm$s.table) && nrow(sm$s.table) > 0) {
    pvec <- sm$s.table[, "p-value"]; names(pvec) <- rownames(sm$s.table)
    ps <- c(ps, pvec)
  }
  drop_candidates <- ps[ps > p_cutoff_ml]
  if (!length(drop_candidates)) { red_ml_cpue_pres <- m_now; break }
  
  ordered <- names(sort(drop_candidates, decreasing = TRUE))
  ordered_vars <- vapply(ordered, function(x) if (grepl("^s\\(", x)) sub("^s\\(([^,]+).*\\)$", "\\1", x) else x, character(1))
  if (PROTECT_FISH) {
    ordered_vars <- setdiff(ordered_vars, fish_vars)
    if (!length(ordered_vars)) { red_ml_cpue_pres <- m_now; break }
  }
  remove_v <- ordered_vars[1]
  remaining_cpue_pres <- setdiff(remaining_cpue_pres, remove_v)
  if (!length(remaining_cpue_pres)) { red_ml_cpue_pres <- m_now; break }
}

# ML comparison: FULL vs REDUCED 
summary(full_ml_cpue_pres)
summary(red_ml_cpue_pres)
print(AIC(full_ml_cpue_pres, red_ml_cpue_pres))
suppressWarnings(print(anova(full_ml_cpue_pres, red_ml_cpue_pres, test = "Chisq")))

# test the effect of LID
mc_noLID <- red_ml_cpue_pres
mc_withLID <- mgcv::gam(formula(CPUE_pos ~ s(Substrate_index, bs = "ts", k = 10) + s(Temperature, bs = "ts", k = 10) + s(LID, bs="re")),data = CData, family = binomial(link="logit"), method = "ML")

anova(mc_noLID, mc_withLID, test = "Chisq")

# Final REML
final_cpue_pres <- mgcv::gam(formula(red_ml_cpue_pres), data = CData, family = binomial(link="logit"), method = "REML", select = TRUE)
summary(final_cpue_pres)
mgcv::concurvity(final_cpue_pres)
gam.check(final_cpue_pres)

# Plots
p_cpue_pres_smooths <- gratia::draw(final_cpue_pres)

p_cpue_pres_param <- parametric_panel(final_cpue_pres, CData, alpha = alpha_sig, exclude = c("LID"), family = "binomial", fish_covars = fish_vars)

# ROC
pred_pres <- predict(final_cpue_pres, type = "response", exclude = exclude_if_RE(final_cpue_pres))
obs_pres  <- CData$CPUE_pos
roc_obj   <- pROC::roc(obs_pres, pred_pres)
auc_val   <- as.numeric(pROC::auc(roc_obj))
roc_df    <- data.frame(tpr = roc_obj$sensitivities, fpr = 1 - roc_obj$specificities)
p_cpue_pres_roc <- ggplot(roc_df, aes(fpr, tpr)) +
  geom_path() + geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.6, y = 0.1, label = paste0("AUC = ", round(auc_val, 3))) +
  labs(x = "False positive rate", y = "True positive rate")

CPUE_panel_presence <- patchwork::wrap_plots(list(as_plot(p_cpue_pres_smooths), as_plot(p_cpue_pres_param), as_plot(p_cpue_pres_roc)), ncol = 3)
CPUE_panel_presence
ggsave(file.path(fig_dir, "CPUE_panel_presence_05.png"), CPUE_panel_presence, width = 14, height = 5, dpi = 300)

# ---------- Positives component (Gamma, log link)
Cpos <- dplyr::filter(CData, CPUE_pos == 1L)

pred_fixed_cpue_pos <- setdiff(names(Cpos), c("Weighted_CPUE_Kōura","LID"))
if (PROTECT_FISH) {
  vif_cpue_pos <- remove_high_vif_glmmTMB(Cpos, "Weighted_CPUE_Kōura", pred_fixed_cpue_pos, threshold = vif_thresh, protect_vars = fish_vars)
} else {
  vif_cpue_pos <- remove_high_vif_glmmTMB(Cpos, "Weighted_CPUE_Kōura", pred_fixed_cpue_pos, threshold = vif_thresh, protect_vars = character(0))
}
kept_fixed_cpue_pos <- vif_cpue_pos$predictors

vars_step_cpue_pos <- c(kept_fixed_cpue_pos, if (INCLUDE_RE) "LID")
rhs_cpue_pos <- paste(vapply(vars_step_cpue_pos, function(v) build_smooth_gam(v, Cpos, custom_k, include_re_for_LID = INCLUDE_RE), character(1)), collapse = " + ")
full_ml_cpue_pos <- mgcv::gam(as.formula(paste("Weighted_CPUE_Kōura ~", rhs_cpue_pos)),
                              data = Cpos, family = Gamma(link="log"), method = "ML", select = TRUE)

# Stepwise (ML)
remaining_cpue_pos <- vars_step_cpue_pos
repeat {
  rhs_now <- paste(vapply(remaining_cpue_pos, function(v) build_smooth_gam(v, Cpos, custom_k, include_re_for_LID = INCLUDE_RE), character(1)), collapse = " + ")
  m_now <- mgcv::gam(as.formula(paste("Weighted_CPUE_Kōura ~", rhs_now)), data = Cpos, family = Gamma(link="log"), method = "ML", select = TRUE)
  sm <- summary(m_now)
  
  ps <- c()
  if (!is.null(sm$p.table) && nrow(sm$p.table) > 0) {
    pcol <- intersect(colnames(sm$p.table), c("Pr(>|t|)","Pr(>|z|)"))[1]
    pvec <- sm$p.table[, pcol]
    ps <- c(ps, pvec[names(pvec) != "(Intercept)"])
  }
  if (!is.null(sm$s.table) && nrow(sm$s.table) > 0) {
    pvec <- sm$s.table[, "p-value"]; names(pvec) <- rownames(sm$s.table)
    ps <- c(ps, pvec)
  }
  drop_candidates <- ps[ps > p_cutoff_ml]
  if (!length(drop_candidates)) { red_ml_cpue_pos <- m_now; break }
  
  ordered <- names(sort(drop_candidates, decreasing = TRUE))
  ordered_vars <- vapply(ordered, function(x) if (grepl("^s\\(", x)) sub("^s\\(([^,]+).*\\)$", "\\1", x) else x, character(1))
  if (PROTECT_FISH) {
    ordered_vars <- setdiff(ordered_vars, fish_vars)
    if (!length(ordered_vars)) { red_ml_cpue_pos <- m_now; break }
  }
  remove_v <- ordered_vars[1]
  remaining_cpue_pos <- setdiff(remaining_cpue_pos, remove_v)
  if (!length(remaining_cpue_pos)) { red_ml_cpue_pos <- m_now; break }
}

# ML comparison: FULL vs REDUCED
summary(full_ml_cpue_pos)
summary(red_ml_cpue_pos)
print(AIC(full_ml_cpue_pos, red_ml_cpue_pos))
suppressWarnings(print(anova(full_ml_cpue_pos, red_ml_cpue_pos, test = "Chisq")))

# test the effect of LID
mcpos_noLID <- red_ml_cpue_pos
mcpos_withLID <- mgcv::gam(formula(Weighted_CPUE_Kōura ~ s(Slope_5m, bs = "ts", k = 10) + s(pH, bs = "ts", 
                                                                                            k = 10) + s(Emergent_Native, bs = "ts", k = 4) + s(LID, bs="re")),data = Cpos, family = Gamma(link="log"), method = "ML")

anova(mcpos_noLID, mcpos_withLID, test = "Chisq")


# Final REML
final_cpue_pos <- mgcv::gam(formula(red_ml_cpue_pos), data = Cpos, family = Gamma(link="log"), method = "REML", select = TRUE)
summary(final_cpue_pos)
mgcv::concurvity(final_cpue_pos)
gam.check(final_cpue_pos)

# Plots
p_cpue_pos_smooths <- gratia::draw(final_cpue_pos)

p_cpue_pos_param <- parametric_panel(final_cpue_pos, Cpos, alpha = alpha_sig, exclude = c("LID"),family = "gamma", fish_covars = fish_vars)

# Scatter diagnostics
pred_pos <- predict(final_cpue_pos, type = "response", exclude = exclude_if_RE(final_cpue_pos))
obs_pos  <- final_cpue_pos$y
R2   <- cor(pred_pos, obs_pos)^2
RMSE <- sqrt(mean((pred_pos - obs_pos)^2))
df_scat <- data.frame(Predicted = pred_pos, Observed = obs_pos)
p_cpue_pos_scatter <- ggplot(df_scat, aes(Predicted, Observed)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text",
           x = max(df_scat$Predicted, na.rm = TRUE) * 0.7,
           y = max(df_scat$Observed,  na.rm = TRUE) * 0.8,
           label = paste0("R² = ", round(R2, 2), "\nRMSE = ", round(RMSE, 2))) +
  labs(x = "Predicted CPUE", y = "Observed CPUE")

CPUE_panel_positives <- patchwork::wrap_plots(list(as_plot(p_cpue_pos_smooths), as_plot(p_cpue_pos_param), as_plot(p_cpue_pos_scatter)), ncol = 3)

CPUE_panel_positives

ggsave(file.path(fig_dir, "CPUE_panel_positives_05.png"), CPUE_panel_positives, width = 14, height = 5, dpi = 300)

```

# BCUE: HURDLE
```{r fig.width=10, fig.height=5, dpi=100, echo=FALSE, message=FALSE, warning=FALSE}
BData <- prepare_block(Modeling_data, "Weighted_BCUE_Kōura", vars_common) %>%
  dplyr::mutate(BCUE_pos = as.integer(Weighted_BCUE_Kōura > 0))

# ---------- Presence component (binomial)
pred_fixed_bcue_pres <- setdiff(names(BData), c("BCUE_pos","LID","Weighted_BCUE_Kōura","Weighted_BCUE_Koura","Weighted_BCUE"))
if (PROTECT_FISH) {
  vif_bcue_pres <- remove_high_vif_glmBI(BData, "BCUE_pos", pred_fixed_bcue_pres, threshold = vif_thresh, protect_vars = fish_vars)
} else {
  vif_bcue_pres <- remove_high_vif_glmBI(BData, "BCUE_pos", pred_fixed_bcue_pres, threshold = vif_thresh, protect_vars = character(0))
}
kept_fixed_bcue_pres <- vif_bcue_pres$predictors

vars_step_bcue_pres <- c(kept_fixed_bcue_pres, if (INCLUDE_RE) "LID")
rhs_bcue_pres <- paste(vapply(vars_step_bcue_pres, function(v) build_smooth_gam(v, BData, custom_k, include_re_for_LID = INCLUDE_RE), character(1)), collapse = " + ")
full_ml_bcue_pres <- mgcv::gam(as.formula(paste("BCUE_pos ~", rhs_bcue_pres)),
                               data = BData, family = binomial(link="logit"), method = "ML", select = TRUE)

remaining_bcue_pres <- vars_step_bcue_pres
repeat {
  rhs_now <- paste(vapply(remaining_bcue_pres, function(v) build_smooth_gam(v, BData, custom_k, include_re_for_LID = INCLUDE_RE), character(1)), collapse = " + ")
  m_now <- mgcv::gam(as.formula(paste("BCUE_pos ~", rhs_now)), data = BData, family = binomial(link="logit"), method = "ML", select = TRUE)
  sm <- summary(m_now)
  
  ps <- c()
  if (!is.null(sm$p.table) && nrow(sm$p.table) > 0) {
    pcol <- intersect(colnames(sm$p.table), c("Pr(>|t|)","Pr(>|z|)"))[1]
    pvec <- sm$p.table[, pcol]
    ps <- c(ps, pvec[names(pvec) != "(Intercept)"])
  }
  if (!is.null(sm$s.table) && nrow(sm$s.table) > 0) {
    pvec <- sm$s.table[, "p-value"]; names(pvec) <- rownames(sm$s.table)
    ps <- c(ps, pvec)
  }
  drop_candidates <- ps[ps > p_cutoff_ml]
  if (!length(drop_candidates)) { red_ml_bcue_pres <- m_now; break }
  
  ordered <- names(sort(drop_candidates, decreasing = TRUE))
  ordered_vars <- vapply(ordered, function(x) if (grepl("^s\\(", x)) sub("^s\\(([^,]+).*\\)$", "\\1", x) else x, character(1))
  if (PROTECT_FISH) {
    ordered_vars <- setdiff(ordered_vars, fish_vars)
    if (!length(ordered_vars)) { red_ml_bcue_pres <- m_now; break }
  }
  remove_v <- ordered_vars[1]
  remaining_bcue_pres <- setdiff(remaining_bcue_pres, remove_v)
  if (!length(remaining_bcue_pres)) { red_ml_bcue_pres <- m_now; break }
}

# ML comparison: FULL vs REDUCED (BCUE
summary(full_ml_bcue_pres)
summary(red_ml_bcue_pres)
print(AIC(full_ml_bcue_pres, red_ml_bcue_pres))
suppressWarnings(print(anova(full_ml_bcue_pres, red_ml_bcue_pres, test = "Chisq")))


# test the effect of LID
mb_noLID <- red_ml_bcue_pres
mb_withLID <- mgcv::gam(formula(BCUE_pos ~ s(Overhanging_trees, bs = "ts", k = 5) + Presence_Goldfish + 
                                  Presence_Eel + s(LID, bs="re")),data = BData, family = binomial(link="logit"), method = "ML")

anova(mb_noLID, mb_withLID, test = "Chisq")

# Final REML
final_bcue_pres <- mgcv::gam(formula(red_ml_bcue_pres), data = BData, family = binomial(link="logit"), method = "REML", select = TRUE)
summary(final_bcue_pres)
mgcv::concurvity(final_bcue_pres)
gam.check(final_bcue_pres)

# make plots
p_bcue_pres_smooths <- gratia::draw(final_bcue_pres)

p_bcue_pres_param <- parametric_panel(final_bcue_pres, BData, alpha = alpha_sig, exclude = c("LID"),family = "binomial", fish_covars = fish_vars)

pred_pres_b <- predict(final_bcue_pres, type = "response", exclude = exclude_if_RE(final_bcue_pres))
obs_pres_b  <- BData$BCUE_pos
roc_obj_b   <- pROC::roc(obs_pres_b, pred_pres_b)
auc_val_b   <- as.numeric(pROC::auc(roc_obj_b))
roc_df_b    <- data.frame(tpr = roc_obj_b$sensitivities, fpr = 1 - roc_obj_b$specificities)
p_bcue_pres_roc <- ggplot(roc_df_b, aes(fpr, tpr)) +
  geom_path() + geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.6, y = 0.1, label = paste0("AUC = ", round(auc_val_b, 3))) +
  labs(x = "False positive rate", y = "True positive rate")

BCUE_panel_presence <- patchwork::wrap_plots(list(as_plot(p_bcue_pres_smooths), as_plot(p_bcue_pres_param), as_plot(p_bcue_pres_roc)), ncol = 3)
BCUE_panel_presence
ggsave(file.path(fig_dir, "BCUE_panel_presence_05.png"), BCUE_panel_presence, width = 14, height = 5, dpi = 300)

# ---------- Positives component (Gamma, log link)
Bpos <- dplyr::filter(BData, BCUE_pos == 1L)

pred_fixed_bcue_pos <- setdiff(names(Bpos), c("Weighted_BCUE_Kōura","LID"))
if (PROTECT_FISH) {
  vif_bcue_pos <- remove_high_vif_glmmTMB(Bpos, "Weighted_BCUE_Kōura", pred_fixed_bcue_pos, threshold = vif_thresh, protect_vars = fish_vars)
} else {
  vif_bcue_pos <- remove_high_vif_glmmTMB(Bpos, "Weighted_BCUE_Kōura", pred_fixed_bcue_pos, threshold = vif_thresh, protect_vars = character(0))
}
kept_fixed_bcue_pos <- vif_bcue_pos$predictors

vars_step_bcue_pos <- c(kept_fixed_bcue_pos, if (INCLUDE_RE) "LID")
rhs_bcue_pos <- paste(vapply(vars_step_bcue_pos, function(v) build_smooth_gam(v, Bpos, custom_k, include_re_for_LID = INCLUDE_RE), character(1)), collapse = " + ")
full_ml_bcue_pos <- mgcv::gam(as.formula(paste("Weighted_BCUE_Kōura ~", rhs_bcue_pos)),
                              data = Bpos, family = Gamma(link="log"), method = "ML", select = TRUE)

remaining_bcue_pos <- vars_step_bcue_pos
repeat {
  rhs_now <- paste(vapply(remaining_bcue_pos, function(v) build_smooth_gam(v, Bpos, custom_k, include_re_for_LID = INCLUDE_RE), character(1)), collapse = " + ")
  m_now <- mgcv::gam(as.formula(paste("Weighted_BCUE_Kōura ~", rhs_now)), data = Bpos, family = Gamma(link="log"), method = "ML", select = TRUE)
  sm <- summary(m_now)
  
  ps <- c()
  if (!is.null(sm$p.table) && nrow(sm$p.table) > 0) {
    pcol <- intersect(colnames(sm$p.table), c("Pr(>|t|)","Pr(>|z|)"))[1]
    pvec <- sm$p.table[, pcol]
    ps <- c(ps, pvec[names(pvec) != "(Intercept)"])
  }
  if (!is.null(sm$s.table) && nrow(sm$s.table) > 0) {
    pvec <- sm$s.table[, "p-value"]; names(pvec) <- rownames(sm$s.table)
    ps <- c(ps, pvec)
  }
  drop_candidates <- ps[ps > p_cutoff_ml]
  if (!length(drop_candidates)) { red_ml_bcue_pos <- m_now; break }
  
  ordered <- names(sort(drop_candidates, decreasing = TRUE))
  ordered_vars <- vapply(ordered, function(x) if (grepl("^s\\(", x)) sub("^s\\(([^,]+).*\\)$", "\\1", x) else x, character(1))
  if (PROTECT_FISH) {
    ordered_vars <- setdiff(ordered_vars, fish_vars)
    if (!length(ordered_vars)) { red_ml_bcue_pos <- m_now; break }
  }
  remove_v <- ordered_vars[1]
  remaining_bcue_pos <- setdiff(remaining_bcue_pos, remove_v)
  if (!length(remaining_bcue_pos)) { red_ml_bcue_pos <- m_now; break }
}

# ML comparison: FULL vs REDUCED
summary(full_ml_bcue_pos)
summary(red_ml_bcue_pos)
print(AIC(full_ml_bcue_pos, red_ml_bcue_pos))
suppressWarnings(print(anova(full_ml_bcue_pos, red_ml_bcue_pos, test = "Chisq")))

# test the effect of LID
mbpos_noLID <- red_ml_bcue_pos
mbpos_withLID <- mgcv::gam(formula(Weighted_BCUE_Kōura ~ s(pH, bs = "ts", k = 10) + s(LID, bs="re")),data = Bpos, family = Gamma(link="log"), method = "ML")

anova(mbpos_noLID, mbpos_withLID, test = "Chisq")


# Final REML
final_bcue_pos <- mgcv::gam(formula(red_ml_bcue_pos), data = Bpos, family = Gamma(link="log"), method = "REML", select = TRUE)
summary(final_bcue_pos)
mgcv::concurvity(final_bcue_pos)
gam.check(final_bcue_pos)

# make plots
p_bcue_pos_smooths <- gratia::draw(final_bcue_pos)

p_bcue_pos_param <- parametric_panel(final_bcue_pos, Bpos, alpha = alpha_sig, exclude = c("LID"),family = "gamma", fish_covars = fish_vars)

pred_pos_b <- predict(final_bcue_pos, type = "response", exclude = exclude_if_RE(final_bcue_pos))
obs_pos_b  <- final_bcue_pos$y
R2_b   <- cor(pred_pos_b, obs_pos_b)^2
RMSE_b <- sqrt(mean((pred_pos_b - obs_pos_b)^2))
df_scat_b <- data.frame(Predicted = pred_pos_b, Observed = obs_pos_b)
p_bcue_pos_scatter <- ggplot(df_scat_b, aes(Predicted, Observed)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text",
           x = max(df_scat_b$Predicted, na.rm = TRUE) * 0.7,
           y = max(df_scat_b$Observed,  na.rm = TRUE) * 0.8,
           label = paste0("R² = ", round(R2_b, 2), "\nRMSE = ", round(RMSE_b, 2))) +
  labs(x = "Predicted BCUE", y = "Observed BCUE")

BCUE_panel_positives <- patchwork::wrap_plots(list(as_plot(p_bcue_pos_smooths),  as_plot(p_bcue_pos_param),as_plot(p_bcue_pos_scatter)), ncol = 3) 

BCUE_panel_positives

ggsave(file.path(fig_dir, "BCUE_panel_positives_05.png"), BCUE_panel_positives, width = 10, height = 5, dpi = 300)