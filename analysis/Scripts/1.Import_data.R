# 1. Import_data
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
packages <- c("tidyverse", "dplyr", "ggplot2","readxl", "writexl","readr")

# Load packages if not already installed
lapply(packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, dependencies = TRUE)
  library(pkg, character.only = TRUE)})


# Import the data sets ---------------------------------------------------------
Site_info <- read_excel("Data_raw/Natural_habitat.xlsx")
Monitoring_data <- read_excel("Data_raw/Natural_habitat.xlsx", sheet = "Monitoring_data")
Weed_data <- read_excel("Data_raw/Natural_habitat.xlsx", sheet = "Weed_data")  %>% dplyr::select(-starts_with("..."))
Fish_data <- read_excel("Data_raw/Natural_habitat.xlsx", sheet = "Fish_data") %>% dplyr::select(-starts_with("..."))
Macroinvertebrates <- read_excel("Data_raw/Natural_habitat.xlsx", sheet = "Macroinvertebrates")



# Combine monitoring, cpue, and natural habitat data ---------------------------
## --- Safe helpers for missing-only vectors ---
safe_min  <- function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
safe_max  <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
safe_mean <- function(x) { m <- mean(x, na.rm = TRUE); if (is.nan(m)) NA_real_ else m }

## --- 0) Tidy base data we’ll need repeatedly ---
fish <- Fish_data %>%
  filter(!is.na(Species)) %>%                  # must have a species label
  filter(!is.na(Monitoring_ID))                # drop unusable rows early

## --- 1) Deployment (effort) table per site × net type ---
deploy <- fish %>%
  filter(!is.na(Net_type), !is.na(Amount_nets)) %>%       # avoid NA groups
  distinct(Monitoring_ID, Net_type, Amount_nets) %>%
  group_by(Monitoring_ID, Net_type) %>%
  summarise(Effort = max(Amount_nets, na.rm = TRUE), .groups = "drop") %>%
  mutate(Effort = ifelse(is.finite(Effort), Effort, 0))   # guard: no -Inf

## --- 2) Catches per site × net type × species (totals) ---
# NOTE: Amount and Weight_g may be missing depending on species/recording.
# Sum them, but keep NAs as 0 for CPUE/BCUE math.
catches <- fish %>%
  group_by(Monitoring_ID, Net_type, Species) %>%
  summarise(
    Total_Individuals = sum(Amount,   na.rm = TRUE),
    Total_Weight      = sum(Weight_g, na.rm = TRUE),
    .groups = "drop")

## --- 3) Ensure zero-catch combinations are present (complete grid) ---
# Start from deploy to make sure ALL nets deployed at a site are represented,
# then expand to all species and fill zero catches for absent combinations.
catches_full <- deploy %>%
  left_join(catches, by = c("Monitoring_ID", "Net_type")) %>%
  complete(Monitoring_ID, Net_type, Species,
           fill = list(Total_Individuals = 0, Total_Weight = 0)) %>%
  group_by(Monitoring_ID, Net_type) %>%
  mutate(Effort = first(Effort)) %>%           # re-fill Effort after complete()
  ungroup() %>%
  mutate(Effort = coalesce(Effort, 0))

## --- 4) Correct CPUE/BCUE per site × species (effort-weighted across nets) ---
CPUE_BCUE <- catches_full %>%
  group_by(Monitoring_ID, Species) %>%
  summarise(
    Total_Individuals = sum(Total_Individuals, na.rm = TRUE),
    Total_Weight      = sum(Total_Weight,      na.rm = TRUE),
    Total_Effort      = sum(Effort,            na.rm = TRUE),
    CPUE              = ifelse(Total_Effort > 0, Total_Individuals / Total_Effort, NA_real_),
    BCUE              = ifelse(Total_Effort > 0, Total_Weight      / Total_Effort, NA_real_),
    .groups = "drop")

## --- 5) True size stats (independent of net type) per site × species ---
# These are computed from raw individuals/records, not means-of-means.
size_stats <- fish %>%
  group_by(Monitoring_ID, Species) %>%
  summarise(
    Mean_Length = safe_mean(Length_mm),
    Min_Length  = safe_min(Length_mm),
    Max_Length  = safe_max(Length_mm),
    Mean_Weight = safe_mean(Weight_g),
    Min_Weight  = safe_min(Weight_g),
    Max_Weight  = safe_max(Weight_g),
    .groups = "drop")

## --- 6) Presence/absence per site × species (+ Predator_Fish_Presence) ---
species_presence <- fish %>%
  distinct(Monitoring_ID, Species) %>%
  mutate(Presence = 1L) %>%
  pivot_wider(
    names_from  = Species,
    values_from = Presence,
    values_fill = list(Presence = 0L),
    names_prefix = "Presence_")

# Ensure predator columns exist; if not, create them as 0
for (pred in c("Trout","Eel","Catfish")) {
  nm <- paste0("Presence_", pred)
  if (!nm %in% names(species_presence)) species_presence[[nm]] <- 0L
}

species_presence <- species_presence %>%
  mutate(
    Predator_Fish_Presence = pmax(
      coalesce(Presence_Trout,   0L),
      coalesce(Presence_Eel,     0L),
      coalesce(Presence_Catfish, 0L)))

## --- 7) Final LONG table (site × species) with CPUE/BCUE + size stats ---
site_species <- CPUE_BCUE %>%
  left_join(size_stats,      by = c("Monitoring_ID","Species"))

## --- 8) Optional: WIDE summary per site (one row per Monitoring_ID) ---
# Include CPUE/BCUE, totals, and size stats per species as columns;
# add richness and an abundance metric (customize exclusions as needed).
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
  # Join presence/absence (includes Predator_Fish_Presence)
  left_join(species_presence, by = "Monitoring_ID") %>%
  # Site-level richness & abundance (customize exclusions below)
  mutate(
    Richness = rowSums(across(starts_with("Total_Individuals_"), ~ . > 0), na.rm = TRUE),
    Abundance = rowSums(across(starts_with("Total_Individuals_"), ~ .),
                        na.rm = TRUE) -
      rowSums(across(matches("^Total_Individuals_(Bullies|Common_smelt)$"), ~ .), na.rm = TRUE))



# Calculate CPUE and BCUE by Site and Species and net_Type
CPUE_BCUE <- Fish_data %>%
  filter(!is.na(Species)) %>%
  group_by(Monitoring_ID, Species,Net_type) %>%
  reframe(
    Total_Individuals = sum(Amount, na.rm = T),
    Total_Weight = sum(Weight_g, na.rm = T),
    Total_Effort = first(Amount_nets),
    CPUE = Total_Individuals / Total_Effort,      
    BCUE = Total_Weight / Total_Effort,
    Mean_Length = mean(Length_mm, na.rm = T),
    Min_Length = ifelse(all(is.na(Length_mm)), NA, min(Length_mm, na.rm = T)),
    Max_Length = ifelse(all(is.na(Length_mm)), NA, max(Length_mm, na.rm = T)),
    Mean_Weight = mean(Weight_g, na.rm = T),
    Min_Weight = ifelse(all(is.na(Weight_g)), NA, min(Weight_g, na.rm = T)),
    Max_Weight = ifelse(all(is.na(Weight_g)), NA, max(Weight_g, na.rm = T)))

CPUE_BCUE_weighted <- CPUE_BCUE %>%
  group_by(Monitoring_ID, Species) %>%
  summarise(
    Total_Individuals = sum(Total_Individuals, na.rm = T),
    Total_Weight = sum(Total_Weight, na.rm = T),
    Weighted_CPUE_numerator = sum(CPUE * Total_Effort, na.rm = T),
    Weighted_BCUE_numerator = sum(BCUE * Total_Effort, na.rm = T),
    Total_Effort_sum = sum(Total_Effort, na.rm = T),
    
    Mean_Length = mean(Mean_Length, na.rm = T),
    Min_Length = ifelse(all(is.na(Min_Length)), NA, min(Min_Length, na.rm = T)),
    Max_Length = ifelse(all(is.na(Max_Length)), NA, max(Max_Length, na.rm = T)),
    
    Mean_Weight = mean(Mean_Weight, na.rm = T),
    Min_Weight = ifelse(all(is.na(Min_Weight)), NA, min(Min_Weight, na.rm = T)),
    Max_Weight = ifelse(all(is.na(Max_Weight)), NA, max(Max_Weight, na.rm = T)) ) %>%
  ungroup() %>%  
  mutate(
    Total_Effort_sum = ifelse(Monitoring_ID %in% c("96_0", "101_0", "117_1", "119_1"), 3, 4),  
    Weighted_CPUE = Weighted_CPUE_numerator / Total_Effort_sum,
    Weighted_BCUE = Weighted_BCUE_numerator / Total_Effort_sum )

# Create Presence/Absence Columns Based on Species Naming
species_presence_absence <- Fish_data %>%
  filter(!is.na(Species)) %>%
  distinct(Monitoring_ID, Species) %>% 
  mutate(Presence = 1) %>%       
  pivot_wider(
    names_from = Species, 
    values_from = Presence, 
    values_fill = list(Presence = 0), 
    names_prefix = "Presence_") %>%
  mutate(Predator_Fish_Presence = pmax(Presence_Trout, Presence_Eel, Presence_Catfish))

# Reshape CPUE_BCUE_weighted so each species becomes a column
CPUE_BCUE_weighted_summary <- CPUE_BCUE_weighted %>%
  pivot_wider(names_from = Species,
              values_from = c(Total_Individuals, Weighted_CPUE, Weighted_BCUE, Total_Weight, Mean_Length, Mean_Weight, Weighted_CPUE_numerator,Weighted_BCUE_numerator,Total_Effort_sum, Min_Length, Max_Length, Min_Weight, Max_Weight),
              names_sep = "_", values_fill = list(Total_Individuals = 0, Weighted_CPUE = 0, Weighted_BCUE = 0))%>%
  mutate(Richness = rowSums(dplyr::select(., starts_with("Total_Individuals_")) > 0), Abundance = rowSums(dplyr::select(., starts_with("Total_Individuals_") & !ends_with(c("_Bullies", "_Common_smelt"))))) 
    





# Create a metadata table with Parameter-Unit mappings
Monitoring_data <- Monitoring_data %>%  dplyr::select(-Site_ID)
unit_metadata <- Monitoring_data %>%
  dplyr::select(Parameter, Unit) %>%
  distinct()

Monitoring_summary <- Monitoring_data %>%
  dplyr::select(-Group, -Notes, -Unit) %>%
  pivot_wider(names_from = c(Parameter),values_from = Value,values_fill = list(Value = NA))

Monitoring_summary <- Monitoring_summary%>%
  mutate(across(c(Bottom_visible, Water_clarity, Depth_10m, Slope, Riparian_vegetation,Vegetation_nearby,
      Overhanging_trees, Erosion, Sructure, Bedrock, Boulders, Cobble, Gravel,
      Sand, Mud, Organic_matter, Rock_size, Temperature, DO_mgl, DO_percent,
      Conductivity, Specific_conductivity , pH, Wood_cover),~ as.numeric(.)))

# Reshape Weed_data so each Weed_Type becomes a column
Weed_summary <- Weed_data %>%
  group_by(Monitoring_ID, Weed_Type, Native_Status) %>%
  summarise(Total_Cover = sum(Percentage_Cover, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = c(Weed_Type, Native_Status), values_from = Total_Cover, values_fill = 0)

Macroinvertebrates_sum <- Macroinvertebrates %>%
  group_by(Monitoring_ID, Species) %>% 
  summarise(Total_amount = sum(Amount, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = c(Species), values_from = Total_amount, values_fill =0)


Macroinvertebrates_sum <- Macroinvertebrates %>%
  group_by(Monitoring_ID, Species) %>% 
  summarise(Total_amount = sum(Amount, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = c(Species), values_from = Total_amount, values_fill = 0) %>%
  mutate(Invertebrates_Richness = rowSums(dplyr::select(., -Monitoring_ID) > 0), 
    Invertebrates_Abundance = rowSums(dplyr::select(., -Monitoring_ID)))


# join together all data
Monitoring_CPUE_data <- Site_info %>%
  left_join(Monitoring_summary, by = "Monitoring_ID") %>%
  left_join(Weed_summary %>% dplyr::select(Monitoring_ID, Emergent_Native,Emergent_Non_Native,Submerged_Native, Submerged_Non_Native, Turf_Native), by = c("Monitoring_ID")) %>%
  left_join(species_wide, by = "Monitoring_ID")%>%
  #left_join(species_presence_absence, by = "Monitoring_ID")%>%
  left_join(Macroinvertebrates_sum, by = "Monitoring_ID")

#Calculate rock presence,  slope to 5m depth, new site-ID, dates and seasons
Monitoring_CPUE_data <- Monitoring_CPUE_data %>%
  mutate(Presence_rocks = if_else(Cobble > 1 | Boulders > 1, 1, 0),
  Slope_5m = 5/Distance_5m,
  Site_ID_ = Site_ID - 60,Monitoring_ID_ = paste0(as.numeric(sub("_.*", "", Monitoring_ID)) - 60, sub("^[^_]*", "", Monitoring_ID)), Monitoring = sub(".*?_", "", Monitoring_ID),
    Date = as.Date(Date_Time),                      
    Time = format(Date_Time, "%H:%M:%S"),
    Year = year(Date_Time),
    Month = month(Date_Time, label = TRUE),          
    Day = day(Date_Time),                             
    Season = case_when(                               
      Month %in% c("Dec", "Jan", "Feb") ~ "Summer",
      Month %in% c("Mar", "Apr", "May") ~ "Autumn",
      Month %in% c("Jun", "Jul", "Aug") ~ "Winter",
      Month %in% c("Sep", "Oct", "Nov") ~ "Spring",
      TRUE ~ NA_character_),
    Date_Time_Numeric = as.numeric(Date_Time))

# Determine the Dominant Habitat Types of each site 
habitat_classification <- Monitoring_CPUE_data %>%
  dplyr::select(Monitoring_ID, DHT, Lake, 
                Bedrock, Boulders, Cobble, Gravel, Sand, Mud, Organic_matter, 
                Emergent_Native, Emergent_Non_Native, Submerged_Native, Submerged_Non_Native, Wood_cover) %>%
  pivot_longer(cols = c(Bedrock, Boulders, Cobble, Gravel, Sand, Mud, Organic_matter, 
                        Emergent_Native, Emergent_Non_Native, Submerged_Native, 
                        Submerged_Non_Native, Wood_cover), 
               names_to = "Type", values_to = "Percentage") %>%
  group_by(Monitoring_ID) %>%
  summarise(
    Rocky_Percentage = sum(Percentage[Type %in% c("Bedrock", "Boulders", "Cobble")], na.rm = TRUE),
    Sand_Percentage  = sum(Percentage[Type == "Sand"], na.rm = TRUE),
    Mud_Percentage   = sum(Percentage[Type %in% c("Mud", "Organic_matter")], na.rm = TRUE),
    Emergent_Percentage = sum(Percentage[Type %in% c("Emergent_Native")], na.rm = TRUE),
    # Substrate index calculation
    Substrate_index = sum(
        0.08 * Percentage[Type == "Bedrock"] +
        0.07 * Percentage[Type == "Boulders"] +
        0.06 * Percentage[Type == "Cobble"] +
        0.04 * Percentage[Type == "Gravel"] +
        0.03 * Percentage[Type == "Sand"] +
        0.02 * Percentage[Type == "Organic_matter"] +
        0.01 * Percentage[Type == "Mud"],
      na.rm = TRUE),.groups = "drop") %>%
  mutate(Habitat_Type = case_when(
    Rocky_Percentage > 25 ~ "Rocky",
    Emergent_Percentage > 25 ~ "Emergent Macrophyte",
    Sand_Percentage >= Mud_Percentage ~ "Sandy",
    TRUE ~ "Muddy")) %>%
  dplyr::select(Monitoring_ID, Habitat_Type, Substrate_index)


# Prepare substrate matrix
substrates <- Monitoring_CPUE_data %>%
  mutate(Mud_OM = Mud + Organic_matter) %>%
  select(Bedrock, Boulders, Cobble, Gravel, Sand, Mud_OM)

substrate_names <- colnames(substrates)

# Find dominant and second dominant substrates
max_idx <- apply(substrates, 1, which.max)
second_max_idx <- apply(substrates, 1, function(x) {
  x_copy <- x
  x_copy[which.max(x)] <- -Inf
  which.max(x_copy)
})

# Save as CSV
write.csv(habitat_classification, "Data_mod/habitat_classification.csv", row.names = FALSE)


#Merge Habitat_Type back into Monitoring_CPUE_data
Monitoring_CPUE_data <- Monitoring_CPUE_data %>%
  left_join(habitat_classification, by = c("Monitoring_ID"))

# Save as CSV
write.csv(Monitoring_CPUE_data, "Data_mod/Monitoring_CPUE_data.csv", row.names = FALSE)

# Save as Excel file
write_xlsx(Monitoring_CPUE_data, "Data_mod/Monitoring_CPUE_data.xlsx")















