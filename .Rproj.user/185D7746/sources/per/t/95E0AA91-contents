# Habitat + Reef data
# Explanation of this script ------------------------------------------------

# combination fo the habitat monitoring and the reef monitoring
# 
#
#


# Clean and load packages ------------------------------------------------------
cat("\014"); rm(list = ls())#; dev.off()
#sapply(.packages(), unloadNamespace)

#Set working derectory
setwd("~/PhD/Data")

# Define the list of packages
packages <- c("sf", "mgcv","tidyverse", "dplyr", "ggplot2","readxl", "writexl","readr")

# Load packages if not already installed
lapply(packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, dependencies = TRUE)
  library(pkg, character.only = TRUE)})


# Import the data sets ---------------------------------------------------------
habitat_Monitoring_CPUE_data <- read_csv("3. Natural habitat monitoring/Data_mod/Monitoring_CPUE_data.csv")
Site_info <- read_excel("3. Natural habitat monitoring/Data_raw/Natural_habitat.xlsx")
Fish_data <- read_excel("3. Natural habitat monitoring/Data_raw/Natural_habitat.xlsx", sheet = "Fish_data") %>% select(-starts_with("..."))
habitat_classification <- read_csv("3. Natural habitat monitoring/Data_mod/habitat_classification.csv")

Reef_Monitoring_CPUE_data <- read_csv("1. Reefs Rotoiti/Data_mod/Monitoring_Reef_data.csv")
Site_info_Reef <- read_excel("1. Reefs Rotoiti/Data_raw/Data_Reefs_Rotoiti.xlsx")
Fish_data_Reef <- read_excel("1. Reefs Rotoiti/Data_raw/Data_Reefs_Rotoiti.xlsx", sheet = "Fish_data") %>% select(-starts_with("..."))
habitat_classification_Reef <- read_csv("1. Reefs Rotoiti/Data_mod/habitat_classification.csv")

# Explore Monitoring data ------------------------------------------------------
All_data <- bind_rows(habitat_Monitoring_CPUE_data, Reef_Monitoring_CPUE_data)

# Filter for only Rotoiti data
Rotoiti_data <- All_data %>%
  filter(Lake == "Rotoiti" & Monitoring == 0)

names(Rotoiti_data)

# Rotoiti modeling

# Define a custom color palette
sediment_colors <- c("Bedrock" = "black","Boulders" = "gray25","Cobble" = "gray55","Gravel" = "gray80","Sand" = "gold","Mud" = "saddlebrown","Organic_matter" = "darkgreen")
weed_colors <- c("Emergent_Native" = "darkgreen","Emergent_Non_Native" = "green3","Submerged_Native" = "skyblue","Submerged_Non_Native" = "royalblue4", "Turf_Native" = "#00bfc4","Wood_cover" = "saddlebrown")
habitat_colors <- c("Rocky" = "gray25", "Sandy" = "gold", "Muddy" = "saddlebrown","Emergent Macrophyte" = "darkgreen")
habitat_order <- c("Rocky", "Sandy", "Muddy", "Emergent Macrophyte")
Rotoiti_data$Habitat_Type <- factor(Rotoiti_data$Habitat_Type, levels = habitat_order)

sediment_data <- Rotoiti_data %>%
  select(Site_ID, Site_ID_,DHT,Habitat_Type , Lake, Bedrock, Boulders, Cobble, Gravel, Sand, Mud, Organic_matter) %>%
  pivot_longer(cols = c(Bedrock, Boulders, Cobble, Gravel, Sand, Mud, Organic_matter), names_to = "Sediment_Type", values_to = "Percentage") %>%
  mutate(Sediment_Type = factor(Sediment_Type,levels = c("Bedrock", "Boulders", "Cobble", "Gravel", "Sand", "Mud", "Organic_matter")))

Predator_data <- Rotoiti_data %>%
  select(Site_ID, Habitat_Type, Weighted_CPUE_Eel, Weighted_CPUE_Trout, Weighted_CPUE_Catfish) %>%
  pivot_longer(cols = c(Weighted_CPUE_Eel, Weighted_CPUE_Trout, Weighted_CPUE_Catfish), names_to = "Fish", values_to = "CPUE")

weed_data <- Rotoiti_data %>%
  select(Site_ID, Site_ID_,DHT,Habitat_Type , Lake, Emergent_Native,Emergent_Non_Native, Submerged_Native, Submerged_Non_Native, Turf_Native,Wood_cover) %>%
  pivot_longer(cols = c(Emergent_Native,Emergent_Non_Native, Submerged_Native, Submerged_Non_Native, Turf_Native,Wood_cover), 
               names_to = "Weeds", values_to = "Percentage")

# statistics


# Plot Physical parameters
Sediment_plot <- ggplot(sediment_data, aes(factor(Site_ID), Percentage, fill = Sediment_Type)) +
  geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.2) +
  scale_fill_manual(values = sediment_colors) +
  #facet_grid( ~ Site_ID, scales = "free_x")+
  labs(x = "Site",y = "Sediment Cover (%)",fill = "Sediment Type") +
  theme_bw()+
  theme(axis.title.x = element_blank())

Substrate_index_plot <- Rotoiti_data %>%
  #mutate(Lake = factor(Lake, levels = substrate_stats$Lake, labels = substrate_stats$Label)) %>%
  ggplot(aes(factor(Site_ID), Substrate_index, fill = Habitat_Type)) +
  geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.2) +
  scale_fill_manual(values = habitat_colors) +
  #facet_grid(~ Site_ID, scales = "free_x") +
  labs(x = "Site", y = "Substrate_index", fill = "Habitat Type") +
  theme_bw()+
  theme(axis.title.x = element_blank())

Koura_plot <- Rotoiti_data %>%
  ggplot(aes(factor(Site_ID), Weighted_CPUE_Kōura)) +
  geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.2) +
  #facet_grid(~ Site_ID, scales = "free_x") +
  labs(x = "Site", y = "Kōura CPUE", fill = "Habitat Type") +
  theme_bw()+
  theme(axis.title.x = element_blank())

Predator_plot <- Predator_data %>%
  ggplot(aes(factor(Site_ID), CPUE, fill=Fish)) +
  geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.2) +
  theme_bw()+
  theme(axis.title.x = element_blank())

Weed_plot <- ggplot(weed_data, aes(factor(Site_ID),Percentage, fill = Weeds)) +
  geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.2) +
  scale_fill_manual(values = weed_colors) +
  #facet_grid( ~ Site_ID, scales = "free_x")+
  labs(x = "Site",y = "Cover (%)",fill = "Cover Type") +
  theme_bw()

Physical_plot <- Sediment_plot/Substrate_index_plot/Koura_plot/Predator_plot/Weed_plot
Physical_plot


# Remake Chemical parameters
Chemical_data <- Rotoiti_data %>%
  select(Site_ID,DHT,Habitat_Type,Lake,DO_mgl,DO_percent, Conductivity, Specific_conductivity, pH, Temperature) %>%
  pivot_longer(cols = c(DO_mgl,DO_percent, Specific_conductivity, pH, Temperature), 
               names_to = "Variable", values_to = "Values") 

# calculate statistics

# Plot Chemical parameters
Chemical_plot <- ggplot(Chemical_data, aes(as.factor(Site_ID), Values)) +
  geom_point() +
  facet_wrap(~ Variable, scales = "free", nrow = 5) +
  theme_bw() 
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

Chemical_plot


# Hierarchical Binomial Model --------------------------------------------------
library(lme4)
vars <- c("Presence_Kōura", #"Lake",
          "Slope_5m", "Riparian_vegetation", "Overhanging_trees","Wood_cover", "Substrate_index",
          "Bedrock","Boulders", "Cobble", "Gravel", "Sand", "Organic_matter",
          "Temperature", "DO_percent", "pH", "Specific_conductivity", "DO_mgl",
          "Emergent_Native","Submerged_Non_Native", "Turf_Native", "Submerged_Native","Emergent_Non_Native","Predator_Fish_Presence",
          "Presence_Kōaro", "Presence_Common_smelt","Presence_Morihana","Presence_Eel","Presence_Catfish")


MW_test <- sapply(vars, function(var) {
  wilcox.test(Rotoiti_data[[var]] ~ Rotoiti_data$Presence_Kōura, exact = FALSE)$p.value})
significant <- names(MW_test)[MW_test < 0.1]
significant

MW_test <- sapply(vars, function(var) {
  wilcox.test(Rotoiti_data[[var]] ~ Rotoiti_data$Presence_Catfish, exact = FALSE)$p.value})
significant <- names(MW_test)[MW_test < 0.1]
significant


# Random intercept for Monitoring_ID to account for repeated samples
model1 <- glmer(
Presence_Catfish   ~ Bedrock+  Emergent_Native+ Emergent_Non_Native+Predator_Fish_Presence+Presence_Kōaro+
    (1 | Presence_Kōura),
  data = Rotoiti_data,
  family = binomial(link = "logit"))

# Model summary
summary(model1)

model2 <- glmer(
    Presence_Kōura ~ Presence_Morihana  + Substrate_index + 
    (1 | Presence_Catfish),
  data = Rotoiti_data,
  family = binomial(link = "logit"))

# Model summary
summary(model2)











#
# Explore Fish data ------------------------------------------------------------
# Define the lake order
lake_order <- c("Rotorua", "Rotoiti", "Rotoehu", "Rotomā", "Ōkāreka")


### do calculations 
DF1 <- Fish_data %>% select(-Site_ID) %>%
  left_join(Site_info, by = "Monitoring_ID") %>% 
  left_join(habitat_classification, by = "Monitoring_ID") %>% 
  filter(!is.na(Monitoring_ID))

DF2 <- Fish_data_Reef %>% select(-Site_ID, -ID) %>%
  left_join(Site_info_Reef, by = "Monitoring_ID") %>% 
  left_join(habitat_classification_Reef, by = "Monitoring_ID") %>% 
  filter(!is.na(Monitoring_ID))

Fish_data_2 <- bind_rows(DF1, DF2)

names(Fish_data_2)
Fish_data_2_summary <- Fish_data_2 %>%
  filter(!is.na(Lake), !is.na(Habitat_Type)) %>%                  
  distinct(Monitoring_ID, Net_type, .keep_all = TRUE) %>% 
  group_by(Habitat_Type, Lake) %>%                                
  summarise(
    Sites_Monitored = n_distinct(Monitoring_ID),         
    Total_Nets = sum(Amount_nets, na.rm = TRUE),  
    .groups = "drop")


species_summary <- Fish_data_2 %>%
  filter(!is.na(Species)) %>%  # Remove rows where Species is NA
  group_by(Species) %>%
  summarize(
    Total_Records = n(),  
    Sites_Present = n_distinct(Site_ID),
    Avg_Length = mean(Length_mm, na.rm = TRUE),  
    Avg_Weight = mean(Weight_g, na.rm = TRUE),   
    Total_Amount = sum(Amount, na.rm = TRUE),    
    Total_Weight = sum(Weight_g, na.rm = TRUE),  
    Sites_Present = n_distinct(Site_ID),
    .groups = 'drop')
species_summary




# Koura ------------------------------------------------------------------------
Koura_data <- Fish_data_2 %>%
  filter(Species == "Kōura")

# Calculate the total count
total_count <- Koura_data %>%
  nrow()

total_count_by_lake <- Koura_data %>%
  group_by(Lake) %>%
  summarise(total_Kōura = n())

total_count_by_lake_DHT <- Koura_data %>%
  group_by(Lake, Habitat_Type) %>%
  summarise(total_Kōura = n())

Fish_data_2_summary <- Fish_data_2_summary %>%
  left_join(total_count_by_lake_DHT, by = c("Lake", "Habitat_Type"))

CPUE_by_lake <- Fish_data_2_summary %>%
  group_by(Lake) %>% 
  summarise(
    total_Kōura = sum(total_Kōura, na.rm = TRUE),
    total_sites = sum(Sites_Monitored, na.rm = TRUE),
    CPUE = total_Kōura / total_sites)


# Calculate size class counts
size_class_counts <- Koura_data %>%
  filter(!is.na(Length_mm)) %>%
  mutate(Size_Class = case_when(
    Length_mm < 23 ~ "Small",
    Length_mm >= 23 & Length_mm < 30 ~ "Medium",
    Length_mm >= 30 ~ "Large"
  )) %>%
  group_by(Size_Class) %>%
  summarise(Count = n())

# Extract counts for each size class
small_count <- size_class_counts %>% filter(Size_Class == "Small") %>% pull(Count)
medium_count <- size_class_counts %>% filter(Size_Class == "Medium") %>% pull(Count)
large_count <- size_class_counts %>% filter(Size_Class == "Large") %>% pull(Count)

# Create the histogram with the total count in the title
ggplot(Koura_data %>% 
         filter(!is.na(Length_mm)),
       aes(Length_mm)) +
  geom_histogram(binwidth = 1, color = "black", na.rm = TRUE) +
  labs(x = "Length (mm)", 
       y = "Count", 
       title = paste("Histogram of Kōura Lengths (Total Count:", total_count,")"))

ggplot(Koura_data %>% 
         filter(!is.na(Length_mm)),
       aes(Length_mm)) +
  geom_histogram(binwidth = 1, color = "black", na.rm = TRUE) +
  geom_vline(xintercept = c(23, 30), linetype = "dashed") +
  annotate("text", x = 15, y = 49, label = paste("Small (n=", small_count, ")", sep = "")) +
  annotate("text", x = 26.5, y = 49, label = paste("Medium (n=", medium_count, ")", sep = "")) +
  annotate("text", x = 35, y = 49, label = paste("Large (n=", large_count, ")", sep = "")) +
  labs(x = "Length (mm)", 
       y = "Count", 
       title = paste("Histogram of Kōura Lengths (Total Count:", sum(size_class_counts$Count), ")"))

ggplot(Koura_data %>%
         filter(!is.na(Length_mm), !is.na(Sex)) %>%
         group_by(Lake) %>%
         mutate(weight = 1 / n()),  # Scale weights for each lake
       aes(x = Length_mm, fill = Sex, weight = weight)) +
  geom_density(alpha = 0.4, position = "identity", color = "black") +
  facet_wrap(~Lake)


ggplot(Koura_data %>% 
         filter(!is.na(Length_mm), !is.na(Sex)),
       aes(Length_mm, fill = Sex)) +
  geom_histogram(binwidth = 3, color = "black", position = "dodge", na.rm = TRUE) +
  labs(x = "OCL Length (mm)", 
       y = "Count", 
       fill = "Gender", 
       title = paste("Histogram of Kōura Lengths by Gender (Total Count:",total_count,")")) +
  facet_wrap(~Net_type ) + # Create facets for Gender
  theme(legend.position = "top")






# Aggregate data by date
Koura_summary <- Koura_data %>%
  mutate(Date = as.Date(Date_Time_out)) %>%
  group_by(Date) %>%
  summarise(
    Total_Caught = sum(Amount, na.rm = TRUE),
    Soft_Shelled = sum(Soft_shelled, na.rm = TRUE),
    Berried = sum(Berried, na.rm = TRUE),
    Youngs = sum(Youngs, na.rm = TRUE))

# Reshape data for plotting
Koura_long <- Koura_summary %>%
  pivot_longer(cols = c(Total_Caught, Soft_Shelled, Berried, Youngs),
               names_to = "Category",
               values_to = "Count")

# Plot
ggplot(Koura_long, aes(as.factor(Date), Count, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge") 

# # # all species ----------------------------------------------------------------
library(ggpmisc)
filtered_data <- Fish_data_2 %>% 
  filter(!Species %in% c("Bullies", "Common_smelt", "Trout", "Eel")) %>% 
  filter(!is.na(Species)) %>% 
  filter(!is.na(Length_mm) & !is.na(Weight_g))%>%
  mutate(log_length = log(Length_mm), log_weight = log(Weight_g))

# Generate the plot with linear model equations
ggplot(filtered_data, aes(log_length, log_weight, col = Sex)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +  # Fit a linear model
  stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label.., sep = "~~~")), 
               formula = y ~ x, parse = TRUE) + # Add equation and R²
  facet_wrap(~Species, scales = "free") +
  labs(x = "Length (mm)", y = "Weight (g)", title = "Length vs. Weight Relationship by Species") +
  theme_minimal()


filtered_data <- Fish_data_2 %>% 
  filter(!Species %in% c("Bullies", "Common_smelt", "Trout", "Eel")) %>% 
  filter(!is.na(Species)) %>% 
  filter(!is.na(Length_mm) & !is.na(Weight_g)) %>%
  mutate(log_length = log(Length_mm), log_weight = log(Weight_g))

# Fit models for each species
species_params <- filtered_data %>%
  group_by(Species, Sex) %>%
  summarise(
    model = list(lm(log_weight ~ log_length)),  # Fit linear model
    .groups = "drop"
  ) %>%
  mutate(
    intercept_log_a = map_dbl(model, ~ coef(.)[1]), # Extract log(a)
    b = map_dbl(model, ~ coef(.)[2]),              # Extract b
    a = exp(intercept_log_a)                       # Back-transform to get a
  ) %>%
  select(Species, a, b)

# Display power-law formulas
species_params %>%
  mutate(formula = paste0("W = ", round(a, 3), " * L^", round(b, 3)))

