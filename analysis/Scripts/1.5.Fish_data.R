# 1.5 Fish_data
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
#Monitoring_data <- read_excel("Data_raw/Natural_habitat.xlsx", sheet = "Monitoring_data")
#Weed_data <- read_excel("Data_raw/Natural_habitat.xlsx", sheet = "Weed_data")  %>% select(-starts_with("..."))
Fish_data <- read_excel("Data_raw/Natural_habitat.xlsx", sheet = "Fish_data") %>% dplyr::select(-starts_with("..."))
habitat_classification <- read_csv("Data_mod/habitat_classification.csv")


# make the data sets -----------------------------------------------------------
# Define the lake order
lake_order <- c("Rotorua", "Rotoiti", "Rotoehu", "Rotomā", "Ōkāreka")

# combine fish data with site info and habitat classification
Fish_data_2 <- Fish_data %>% 
  dplyr::select(-Site_ID) %>%
  left_join(Site_info, by = "Monitoring_ID") %>% 
  left_join(habitat_classification, by = "Monitoring_ID") %>% 
  filter(!is.na(Monitoring_ID)) %>%
  mutate(
    Lake = factor(Lake, levels = lake_order),
    Monitoring = sub(".*?_", "", Monitoring_ID)
  ) %>%
  filter(Monitoring == 0)

Fish_data_2_summary <- Fish_data_2 %>%
  filter(!is.na(Lake), !is.na(Habitat_Type)) %>%                  
  distinct(Monitoring_ID, Net_type, .keep_all = TRUE) %>% 
  group_by(Habitat_Type, Lake) %>%                                
  summarise(Sites_Monitored = n_distinct(Monitoring_ID),         
            Total_Nets = sum(Amount_nets, na.rm = TRUE),.groups = "drop")%>%
  mutate(Lake = factor(Lake, levels = lake_order))


# Explore all Fish data-------------------------------------------------------------
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

ggplot(Fish_data_2 %>% 
         #filter(!Species %in% c("Bullies", "Common_smelt")) %>% 
         filter(!is.na(Species)) %>% 
         filter(!is.na(Length_mm)),
       aes(x = Length_mm, fill = Sex)) +
  geom_histogram(binwidth = 5, position = "dodge") +
  facet_wrap(~Species, scales = "free") +
  labs(title = "Histogram of Fish Length by Species and Sex") 

ggplot(Fish_data_2 %>% 
         filter(!Species %in% c("Bullies", "Common_smelt", "Trout", "Eel")) %>% 
         filter(!is.na(Species)) %>% 
         filter(!is.na(Length_mm) & !is.na(Weight_g)),
       aes(Length_mm, Weight_g, col = Sex)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~Species, scales = "free")



# Koura ------------------------------------------------------------------------
# Calculate the total count
Koura_data <- Fish_data_2 %>%
  filter(Species == "Kōura")

# Calculate the koura numbers in diffrent ways
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
    total_nets = sum(Total_Nets, na.rm=T),
    meanCPUE = total_Kōura / total_nets)

# Calculate mean length and weight for each habitat type and lake
mean_koura_stats <- Koura_data %>%
  group_by(Habitat_Type, Lake) %>%
  summarise(mean_length_mm = mean(Length_mm, na.rm = TRUE),sd_length_mm = sd(Length_mm, na.rm = TRUE),n_length = sum(!is.na(Length_mm)),sem_length_mm = sd_length_mm / sqrt(n_length),
    mean_weight_g = mean(Weight_g, na.rm = TRUE),sd_weight_g = sd(Weight_g, na.rm = TRUE),n_weight = sum(!is.na(Weight_g)),sem_weight_g = sd_weight_g / sqrt(n_weight),
    .groups = "drop")

# Calculate size class counts
size_class_counts <- Koura_data %>%
  #filter(Species == "Kōura", !is.na(Length_mm)) %>%
  mutate(Size_Class = case_when(
    Length_mm < 22 ~ "Small",
    Length_mm >= 22 & Length_mm < 30 ~ "Medium",
    Length_mm >= 30 ~ "Large")) %>%
  group_by(Size_Class) %>%
  summarise(Count = n())

# Extract counts for each size class
small_count <- size_class_counts %>% filter(Size_Class == "Small") %>% pull(Count)
medium_count <- size_class_counts %>% filter(Size_Class == "Medium") %>% pull(Count)
large_count <- size_class_counts %>% filter(Size_Class == "Large") %>% pull(Count)

# Create the histogram with the total count in the title
ggplot(Koura_data %>% 
         filter(Species == "Kōura") %>% 
         filter(!is.na(Length_mm)),
       aes(x = Length_mm )) +
  geom_histogram(binwidth = 1, color = "black", na.rm = TRUE) +
  labs(x = "Length (mm)", 
       y = "Count", 
       title = paste("Histogram of Kōura Lengths (Total Count:", total_count,")"))

ggplot(Koura_data %>% 
         filter(!is.na(Length_mm)),
       aes(x = Length_mm)) +
  geom_histogram(binwidth = 1, color = "black", na.rm = TRUE) +
  geom_vline(xintercept = c(22, 30), linetype = "dashed") +
  annotate("text", x = 11, y = 27, label = paste("Small (n=", small_count, ")", sep = "")) +
  annotate("text", x = 26, y = 27, label = paste("Medium (n=", medium_count, ")", sep = "")) +
  annotate("text", x = 35, y = 27, label = paste("Large (n=", large_count, ")", sep = "")) +
  labs(x = "Length (mm)", y = "Count", title = paste("Histogram of Kōura Lengths (Total Count:", sum(size_class_counts$Count), ")"))

ggplot(Koura_data %>% 
         filter(!is.na(Length_mm), !is.na(Sex)),
       aes(x = Length_mm, fill = Sex)) +
   geom_density(alpha = 0.5, color = "black", position = "identity", na.rm = TRUE) +
  #geom_histogram(binwidth = 3, color = "black", position = "dodge", na.rm = TRUE) +
   labs(x = "OCL Length (mm)", 
       y = "Count", 
       fill = "Gender", 
       title = paste("Histogram of Kōura Lengths by Gender (Total Count:",total_count,")")) +
  facet_wrap(~Net_type ) + # Create facets for Gender
  theme(legend.position = "top")

ggplot(Koura_data, aes(Length_mm,Weight_g, col=Sex))+
  geom_point()+
  geom_smooth()

plot_length <- ggplot(Koura_data %>% 
         filter(!is.na(Length_mm), !is.na(Sex)),
       aes(x = Length_mm, fill = Sex)) +
  geom_histogram(binwidth = 3, color = "black", position = "dodge", na.rm = TRUE) +
  labs(x = "OCL Length (mm)", y = "Count", fill = "Gender") +
  facet_grid(Habitat_Type ~ Lake) + 
  geom_vline(data = mean_koura_stats,aes(xintercept = mean_length_mm),linetype = "dashed") +
  geom_text(data = (Fish_data_2_summary), aes(x = Inf, y = Inf, label = paste("Sites:", Sites_Monitored,"\nTotal Kōura:", total_Kōura)),inherit.aes = FALSE,hjust = 1.1, vjust = 1.1, size = 3) + 
  theme(legend.position = "top") +
  theme_bw()

plot_weight <- ggplot(Koura_data %>% 
         filter(!is.na(Weight_g), !is.na(Sex)),
       aes(x = Weight_g, fill = Sex)) +
  geom_histogram(binwidth = 5, color = "black", position = "dodge", na.rm = TRUE) +
  labs(x = "Weight (g)", y = "Count",fill = "Gender")+
  facet_grid(Habitat_Type ~ Lake) + 
  geom_vline(data = mean_koura_stats,aes(xintercept = mean_weight_g),linetype = "dashed") +
  geom_text(data = (Fish_data_2_summary), aes(x = Inf, y = Inf, label = paste("Sites:", Sites_Monitored,"\nTotal Kōura:", total_Kōura)),inherit.aes = FALSE,hjust = 1.1, vjust = 1.1, size = 3) + 
  theme(legend.position = "top")+
  theme_bw()

plot_length / plot_weight

# other species ----------------------------------------------------------------
# Koaro
total_Catfish <- Fish_data_2 %>%
  filter(Species == "Catfish") %>%
  nrow()

ggplot(Fish_data_2 %>% 
         filter(Species == "Catfish") %>% 
         filter(!is.na(Length_mm)),
       aes(x = Length_mm, fill =Lake)) +
  geom_histogram(binwidth = 8, position = "dodge", na.rm = TRUE) +
  labs(x = "Length (mm)", 
       y = "Count", 
       title = paste("Histogram of Catfish Lengths (Total Count:", total_Catfish,")"))

ggplot(Fish_data_2 %>% 
         filter(Species == "Catfish") %>% 
         filter(!is.na(Length_mm)),
       aes(x = Length_mm)) +
  geom_histogram(binwidth = 8, color = "black", position = "dodge", na.rm = TRUE) +
  labs(x = "Length (mm)", 
       y = "Count", 
       fill = "Gender", 
       title = paste("Histogram of Catfish Lengths")) +
  facet_grid(~Habitat_Type ) + # Ensure each Lake has its own column
  geom_text(data = Fish_data_2_summary, 
            aes(x = Inf, y = Inf, label = paste("Sites:", Sites_Monitored)),
            inherit.aes = FALSE,
            hjust = 1.1, vjust = 1.1, size = 3) +
  theme(legend.position = "top")



# Koaro
total_Kōaro <- Fish_data_2 %>%
  filter(Species == "Kōaro") %>%
  nrow()

ggplot(Fish_data_2 %>% 
         filter(Species == "Kōaro") %>% 
         filter(!is.na(Length_mm)),
       aes(x = Length_mm, fill =Lake)) +
  geom_histogram(binwidth = 8, position = "dodge", na.rm = TRUE) +
  labs(x = "Length (mm)", 
       y = "Count", 
       title = paste("Histogram of Kōaro Lengths (Total Count:", total_Kōaro,")"))

ggplot(Fish_data_2 %>% 
         filter(Species == "Kōaro") %>% 
         filter(!is.na(Length_mm)),
       aes(x = Length_mm)) +
  geom_histogram(binwidth = 8, color = "black", position = "dodge", na.rm = TRUE) +
  labs(x = "Length (mm)", 
       y = "Count", 
       fill = "Gender", 
       title = paste("Histogram of Kōaro Lengths")) +
  facet_grid(~Habitat_Type ) + # Ensure each Lake has its own column
  geom_text(data = Fish_data_2_summary, 
            aes(x = Inf, y = Inf, label = paste("Sites:", Sites_Monitored)),
            inherit.aes = FALSE,
            hjust = 1.1, vjust = 1.1, size = 3) +
  theme(legend.position = "top")


# Morihana
ggplot(Fish_data_2 %>% 
         filter(Species == "Morihana") %>% 
         filter(!is.na(Length_mm)),
       aes(x = Length_mm)) +
  geom_histogram(binwidth = 3, color = "black", position = "dodge", na.rm = TRUE) +
  labs(x = "Length (mm)", 
       y = "Count", 
       fill = "Gender", 
       title = paste("Histogram of Morihana Lengths")) +
  facet_grid(DHT ~ Lake) + 
  geom_text(data = Fish_data_2_summary, 
            aes(x = Inf, y = Inf, label = paste("Sites:", Sites_Monitored)),
            inherit.aes = FALSE,
            hjust = 1.1, vjust = 1.1, size = 3) +
  theme(legend.position = "top")


# Mosquitofish
Mosquitofish_data <- Fish_data_2 %>%
  filter(Species == "Mosquitofish")



# Bullies
Bullies_data <- Fish_data_2 %>% filter(Species == "Bullies") 
summary(Bullies_data)


# Check the rows where both 'Weight_g' and 'Amount' are available
weight_and_amount_available <- Bullies_data %>%
  filter(!is.na(Weight_g) & !is.na(Amount))

# Calculate the weight per bully for each row
weight_and_amount_available <- weight_and_amount_available %>%
  mutate(weight_per_bully = Weight_g / Amount)

# Calculate the average weight per bully across all rows
Average_Weight_per_Bully_g <- mean(weight_and_amount_available$weight_per_bully, na.rm = TRUE)

# Calculate total weight, total amount, total volume, and average weight in a single summarise call
total_bullies_info <- Bullies_data %>%
  group_by(Site_ID) %>%
  summarise(
    Total_Weight_g = sum(Weight_g, na.rm = TRUE),
    Total_Amount = sum(Amount, na.rm = TRUE),
    Total_Cups = sum(Cups, na.rm =TRUE),
    Total_Volume_mL = sum(mL, na.rm = TRUE)) %>%
  left_join(average_weight_per_bully, by = "Site_ID")%>%
  mutate(Average_Weight_per_Cup = Total_Weight_g / Total_Cups)




# bullies check
Bullie_check <- M_C_data %>% select(Site_ID, DHT, Total_Individuals_Bullies, Weighted_CPUE_Bullies, 
                                    Weighted_BCUE_Bullies, Total_Weight_Bullies, Mean_Length_Bullies, 
                                    Mean_Weight_Bullies, Min_Length_Bullies, Max_Length_Bullies, 
                                    Min_Weight_Bullies, Max_Weight_Bullies)







# plot CPUE and BCUE by Site and Species and net_Type
ggplot(CPUE_BCUE_weighted, aes(x = factor(Monitoring_ID), y = Weighted_CPUE , fill = Species)) +
  geom_bar(stat = "identity", position = "dodge") 

ggplot(CPUE_BCUE_weighted, aes(x = factor(Monitoring_ID), y = Weighted_BCUE, fill = Species)) +
  geom_bar(stat = "identity", position = "dodge") 


ggplot(CPUE_BCUE_weighted %>% filter(Species %in% c("Kōura", "Kōaro")),
       aes(x = factor(Monitoring_ID), y = Weighted_CPUE, fill = Species)) +
  geom_bar(stat = "identity", position = "dodge") 


ggplot(CPUE_BCUE_weighted, aes(Weighted_BCUE, Weighted_CPUE))+
  geom_point()+
  facet_wrap(~Species, scales = "free")

