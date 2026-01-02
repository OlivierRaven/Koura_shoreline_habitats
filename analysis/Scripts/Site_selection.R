# Selection of the stratified random sites on the dominant habitat types
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
packages <- c("spsurvey","sf","readr", "tidyverse", "dplyr", "ggplot2", "ggrepel")

# Load packages if not already installed
lapply(packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, dependencies = TRUE)
  library(pkg, character.only = TRUE)})

# Set custom_colors ------------------------------------------------------------
custom_colors <- c("Muddy" = "tan", "Sandy" = "orange","Rocky" = "gray", "Raupo" = "green", "AR"= "red",  "Cliff"= "black", "Geo"= "blue") 

# Import the data sets ---------------------------------------------------------
Rotorua_DHT <- st_read("~/PhD/Pictures_Figures/Maps/Layers/Dominant_habitat_types/Rotorua/DHT_rotorua_disolved.gpkg")
Rotoiti_DHT <- st_read("~/PhD/Pictures_Figures/Maps/Layers/Dominant_habitat_types/Rotoiti/DHT_Rotoiti_Dissolve.gpkg")
Rotoehu_DHT <- st_read("~/PhD/Pictures_Figures/Maps/Layers/Dominant_habitat_types/Rotoehu/DHT_Rotoehu_Dissolved.gpkg")
Rotomā_DHT  <- st_read("~/PhD/Pictures_Figures/Maps/Layers/Dominant_habitat_types/Rotomā/DHT_Rotoma_Dissolved.gpkg")
Ōkāreka_DHT <- st_read("~/PhD/Pictures_Figures/Maps/Layers/Dominant_habitat_types/Ōkāreka/DHT_Okaroka_Dissolve.gpkg")

# Start with lake  ------------------------------------------------------
Lake_DHT <- Rotoehu_DHT %>% select(-id)

# Reproject the geometry to a projected CRS (NZGD2000 New Zealand Transverse Mercator - EPSG:2193)
Lake_DHT <- st_transform(Lake_DHT, crs = 2193)

# Calculate the length of each geometry in meters
Lake_DHT <- Lake_DHT %>%
  mutate(length_m = as.numeric(st_length(geom)))  # Convert to numeric to allow for summation

# Group by DHT and summarize the total length
dht_summary <- Lake_DHT %>%
  group_by(DHT) %>%
  summarize(Total_Length_m = sum(length_m, na.rm = TRUE))  # Summarize total lengths

# Define total number of sites
total_sites <- 12

# Calculate the proportion of sites for each DHT category
dht_summary <- dht_summary %>%
  mutate(Proportion = Total_Length_m / sum(Total_Length_m),
         Sites = round(Proportion * total_sites))
dht_summary

caty_n_Rotorua <- c(Sandy = 9, Rocky = 3, AR = 0, Geo = 0) # Rotorua_DHT
caty_n_Rotoiti <- c(Sandy = 8, Rocky = 2, Raupo = 2, AR = 0, Cliff = 0, Geo = 0)  # Rotoiti_DH
caty_n_Rotoehu <- c(Sandy = 8, Rocky = 2, Raupo = 2, Muddy = 0, AR = 0, Cliff = 0, Geo = 0)  # Rotoehu_DHT
caty_n_Rotomā  <- c(Sandy = 6, Rocky = 3, Raupo = 3, AR = 0)    # Rotomā_DHT
caty_n_Ōkāreka <- c(Sandy = 4, Rocky = 3, Raupo = 3, Muddy = 2) # Ōkāreka_DHT

set.seed(1)
# Select 12 random points, stratified by DHT
random_points <- grts(
  sframe = Lake_DHT,
  n_base = 12,          # Total number of points
  caty_var = "DHT",     # Column with habitat types
  caty_n = caty_n_Rotorua
  #projcrs_check = FALSE  # Disable the projected CRS check if needed (but it's better to reproject correctly)
)
random_points


# Extract the points from the 'grts' result
random_points_sf <- st_as_sf(random_points$sites_base, coords = c("xcoord", "ycoord"), crs = 2193)

# Plot the geometries with colors based on DHT
plot(st_geometry(Lake_DHT), 
  col = custom_colors[Lake_DHT$DHT], 
  main = "Rotoiti DHT Types with Random Points",  lwd = 4)

# Plot the selected random points
plot(st_geometry(random_points_sf), add = TRUE, pch = 21, bg="red", col = "black")
text(st_coordinates(random_points_sf), labels = 1:nrow(random_points_sf), cex = 0.8, pos = 3)

# Extract the coordinates and DHT values
coords_dht <- data.frame(
  Longitude = st_coordinates(random_points_sf)[, "X"],  # Extract longitude
  Latitude = st_coordinates(random_points_sf)[, "Y"],   # Extract latitude
  DHT = random_points_sf$DHT                        # Extract DHT values
)
coords_dht

# Create an sf object with the original CRS (NZTM)
coords_sf <- st_as_sf(coords_dht, coords = c("Longitude", "Latitude"), crs = 2193)

# Transform to WGS84
coords_wgs84 <- st_transform(coords_sf, crs = 4326)



# Rotomā_DHT ------------------------------------------------------------------
Rotomā_DHT <- st_read("~/PhD/Pictures_Figures/Maps/Layers/Dominant_habitat_types/Rotomā/DHT_RotomaII.gpkg")

# Reproject to a suitable projected CRS
DHT_proj <- st_transform(Rotomā_DHT, crs = 2193) 

unique(DHT_proj$DHT)

# Define the number of points per category
caty_n <- c(AR = 0, Sandy = 6, Rocky =4, Raupo = 2)

# Set seed for reproducibility
set.seed(1)

# Select 12 random points, stratified by DHT
random_points <- grts(
  sframe = DHT_proj,
  n_base = 12, 
  caty_var = "DHT", 
  caty_n = caty_n,
  mindis = 800
)

# Extract points and convert to sf object
random_points_sf <- st_as_sf(random_points$sites_base, coords = c("xcoord", "ycoord"), crs = 2193)

# Plot DHT geometries
plot(st_geometry(DHT_proj), col = custom_colors[DHT_proj$DHT], main = "Rotomā DHT Types with Random Points", lwd = 2)

# Plot and label random points
plot(st_geometry(random_points_sf), add = TRUE, col = "red", pch = 16)
text(st_coordinates(random_points_sf), labels = 1:nrow(random_points_sf), cex = 0.8, pos = 3)

# Extract coordinates and DHT values
coords_dht <- data.frame(
  Number = 1:nrow(random_points_sf),
  Longitude = st_coordinates(random_points_sf)[, "X"],
  Latitude = st_coordinates(random_points_sf)[, "Y"],
  DHT = random_points_sf$DHT
)
coords_dht



################################################################################
# Rororua DHT ----------------------------------------------------------------------
Rotorua_DHT <- st_read("~/PhD/Pictures_Figures/Maps/Layers/Dominant_habitat_types/Rotorua/DHT_rotorua.gpkg")

# Calculate centroids with correct column names
Rotorua_DHT_centroids <- st_centroid(Rotorua_DHT)
Rotorua_DHT_centroids <- cbind(Rotorua_DHT_centroids, st_coordinates(Rotorua_DHT_centroids))

ggplot() +
  geom_sf(data = Rotorua_DHT, aes(col = DHT), lwd = 1.5) +
  coord_sf() +
  scale_color_manual(values = custom_colors) +
  labs(x = "Longitude", y = "Latitude") + 
  #geom_text_repel(data = Rotorua_DHT_centroids, aes(label = id, x = X, y = Y), size = 3, color = "black", max.overlaps = Inf)+
  xlim(176.23, 176.28) +  
  ylim(-38.15, -38.10)  


# Rotoiti DHT ----------------------------------------------------------------------
Rotoiti_DHT <- st_read("~/PhD/Pictures_Figures/Maps/Layers/Dominant_habitat_types/Rotoiti/DHT_RotoitiIII.gpkg")

Rotoiti_DHT_centroids <- st_centroid(Rotoiti_DHT)
Rotoiti_DHT_centroids <- cbind(Rotoiti_DHT_centroids, st_coordinates(Rotoiti_DHT_centroids))

ggplot() +
  geom_sf(data = Rotoiti_DHT, aes(col = DHT), lwd = 1.5) +
  coord_sf() +
  scale_color_manual(values = custom_colors) +
  labs(x = "Longitude", y = "Latitude") +
  geom_text_repel(data = Rotoiti_DHT_centroids, aes(label = id, x = X, y = Y), size = 2, color = "black", max.overlaps = Inf)#+
#xlim(176.33, 176.40) +  
#ylim(-38.055, -38.01)

DHT_length_Rotoiti <- Rotoiti_DHT %>%
  group_by(DHT) %>%
  summarize(total_length = sum(st_length(geom)))
total_length_Rotoiti <- sum(DHT_length_Rotoiti$total_length)
DHT_length_Rotoiti
total_length_Rotoiti

#### random points
# Filter out the categories "Cliff" and "Geo"
DHT_length_RotoitiI <- DHT_length_Rotoiti %>% 
  filter(!DHT %in% c("Cliff", "Geo", "AR"))

# Transform to UTM Zone 60S for all geometries
DHT_length_Rotoiti_utm <- st_transform(DHT_length_RotoitiI, crs = 32760)  # UTM Zone 60S

# Function to generate random points for each category
set.seed(112)
generate_random_points <- function(data, n_points = 3) {
  st_sample(data, size = n_points, type = "random")
}



# Generate random points and transform back to WGS 84
random_points_list <- lapply(split(DHT_length_Rotoiti_utm, DHT_length_Rotoiti_utm$DHT), function(cat_data) {
  points_utm <- generate_random_points(cat_data)
  st_transform(points_utm, crs = st_crs(DHT_length_Rotoiti))
})

# Combine random points into one sf object with DHT attributes
all_random_points <- st_sf(
  DHT = rep(names(random_points_list), sapply(random_points_list, length)),
  geometry = st_sfc(do.call(c, random_points_list), crs = st_crs(DHT_length_RotoitiI)))

# Extract coordinates for the random points
coordinates <- st_coordinates(all_random_points)

# Add point number as labels
Rotoiti_DHT$DHT <- factor(Rotoiti_DHT$DHT, levels = names(custom_colors))
label_df <- data.frame(coordinates)
label_df$label <- seq_len(nrow(label_df))

# plot
plot_Rotoiti = ggplot() +
  geom_sf(data = Rotoiti_DHT, aes(col = DHT), lwd = 1.5) +
  geom_sf(data = all_random_points, color = "black", size = 2, shape = 21, fill = "red") +
  geom_text_repel(data = label_df, aes(x = X, y = Y, label = label), size = 3) +
  coord_sf() +
  scale_color_manual(values = custom_colors) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw()

plot_Rotoiti

#
# Rotoehu DHT ------------------------------------------------------------------
Rotoehu_DHT <- st_read("~/PhD/Pictures_Figures/Maps/Layers/Dominant_habitat_types/Rotoehu/DHT_RotoehuII.gpkg")

Rotoehu_DHT_centroids <- st_centroid(Rotoehu_DHT)
Rotoehu_DHT_centroids <- cbind(Rotoehu_DHT_centroids, st_coordinates(Rotoehu_DHT_centroids))

ggplot() +
  geom_sf(data = Rotoehu_DHT, aes(col = DHT), lwd = 1.5) +
  coord_sf() +
  scale_color_manual(values = custom_colors) +
  labs(x = "Longitude", y = "Latitude") +
  geom_text_repel(data = Rotoehu_DHT_centroids, aes(label = id, x = X, y = Y), size = 3, color = "black", max.overlaps = Inf)

DHT_length_Rotoehu <- Rotoehu_DHT %>%
  group_by(DHT) %>%
  summarize(total_length = sum(st_length(geom)))
total_length_Rotoehu <- sum(DHT_length_Rotoehu$total_length)
DHT_length_Rotoehu
total_length_Rotoehu


#### random points
# Filter out the categories "Cliff" and "Geo"
DHT_length_RotoehuI <- DHT_length_Rotoehu %>% 
  filter(!DHT %in% c("Cliff", "Geo", "AR"))

# Transform to UTM Zone 60S for all geometries
DHT_length_Rotoehu_utm <- st_transform(DHT_length_RotoehuI, crs = 32760)  # UTM Zone 60S

# Function to generate random points for each category
set.seed(654321)
generate_random_points <- function(data, n_points = 3) {
  st_sample(data, size = n_points, type = "random")
}

# Generate random points and transform back to WGS 84
random_points_list <- lapply(split(DHT_length_Rotoehu_utm, DHT_length_Rotoehu_utm$DHT), function(cat_data) {
  points_utm <- generate_random_points(cat_data)
  st_transform(points_utm, crs = st_crs(DHT_length_Rotoehu))
})

# Combine random points into one sf object with DHT attributes
all_random_points <- st_sf(
  DHT = rep(names(random_points_list), sapply(random_points_list, length)),
  geometry = st_sfc(do.call(c, random_points_list), crs = st_crs(DHT_length_RotoehuI)))

# Extract coordinates for the random points
coordinates <- st_coordinates(all_random_points)

# Add point number as labels
Rotoehu_DHT$DHT <- factor(Rotoehu_DHT$DHT, levels = names(custom_colors))
label_df <- data.frame(coordinates)
label_df$label <- seq_len(nrow(label_df))

# plot
plot_Rotoehu = ggplot() +
  geom_sf(data = Rotoehu_DHT, aes(col = DHT), lwd = 1.5) +
  geom_sf(data = all_random_points, color = "black", size = 2, shape = 21, fill = "red") +
  geom_text_repel(data = label_df, aes(x = X, y = Y, label = label), size = 3) +
  coord_sf() +
  scale_color_manual(values = custom_colors) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw()

plot_Rotoehu

# Rotomā DHT -------------------------------------------------------------------
Rotomā_DHT <- st_read("~/PhD/Pictures_Figures/Maps/Layers/Dominant_habitat_types/Rotomā/DHT_RotomaII.gpkg")

Rotomā_DHT_centroids <- st_centroid(Rotomā_DHT)
Rotomā_DHT_centroids <- cbind(Rotomā_DHT_centroids, st_coordinates(Rotomā_DHT_centroids))

ggplot() +
  geom_sf(data = Rotomā_DHT, aes(col = DHT), lwd = 1.5) +
  coord_sf() +
  scale_color_manual(values = custom_colors) +
  labs(x = "Longitude", y = "Latitude") +
  geom_text_repel(data = Rotomā_DHT_centroids, aes(label = id, x = X, y = Y), size = 3, color = "black", max.overlaps = Inf)

DHT_length_Rotomā <- Rotomā_DHT %>%
  group_by(DHT) %>%
  summarize(total_length = sum(st_length(geom)))
total_length_Rotomā <- sum(DHT_length_Rotomā$total_length)
DHT_length_Rotomā
total_length_Rotomā


#### random points
# Filter out the categories "Cliff" and "Geo"
DHT_length_RotomāI <- DHT_length_Rotomā %>% 
  filter(!DHT %in% c("Cliff", "Geo", "AR"))

# Transform to UTM Zone 60S for all geometries
DHT_length_Rotomā_utm <- st_transform(DHT_length_RotomāI, crs = 32760)  # UTM Zone 60S

# Function to generate random points for each category
set.seed(65)
generate_random_points <- function(data, n_points = 3) {
  st_sample(data, size = n_points, type = "random")
}






#
# Define the number of points to be sampled for each habitat type
points_per_habitat <- list(
  Raupo = 3,
  Rocky = 3,
  Sandy = 6
)

# Generate random points for each habitat type and transform back to WGS 84
random_points_list <- lapply(names(points_per_habitat), function(habitat) {
  cat_data <- DHT_length_Rotoiti_utm[DHT_length_Rotoiti_utm$DHT == habitat, ]
  points_utm <- generate_random_points(cat_data, n_points = points_per_habitat[[habitat]])
  st_transform(points_utm, crs = st_crs(DHT_length_Rotoiti))
})




#




# Generate random points and transform back to WGS 84
random_points_list <- lapply(split(DHT_length_Rotomā_utm, DHT_length_Rotomā_utm$DHT), function(cat_data) {
  points_utm <- generate_random_points(cat_data)
  st_transform(points_utm, crs = st_crs(DHT_length_Rotomā))
})

# Combine random points into one sf object with DHT attributes
all_random_points <- st_sf(
  DHT = rep(names(random_points_list), sapply(random_points_list, length)),
  geometry = st_sfc(do.call(c, random_points_list), crs = st_crs(DHT_length_RotomāI)))

# Extract coordinates for the random points
coordinates <- st_coordinates(all_random_points)

# Add point number as labels
Rotomā_DHT$DHT <- factor(Rotomā_DHT$DHT, levels = names(custom_colors))
label_df <- data.frame(coordinates)
label_df$label <- seq_len(nrow(label_df))

# plot
plot_Rotomā = ggplot() +
  geom_sf(data = Rotomā_DHT, aes(col = DHT), lwd = 1.5) +
  geom_sf(data = all_random_points, color = "black", size = 2, shape = 21, fill = "red") +
  geom_text_repel(data = label_df, aes(x = X, y = Y, label = label), size = 3) +
  coord_sf() +
  scale_color_manual(values = custom_colors) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw()

plot_Rotomā

# Ōkāreka DHT ------------------------------------------------------------------
Ōkāreka_DHT <- st_read("~/PhD/Pictures_Figures/Maps/Layers/Dominant_habitat_types/Ōkāreka/Okaroka_DHT.gpkg")

Ōkāreka_DHT_centroids <- st_centroid(Ōkāreka_DHT)
Ōkāreka_DHT_centroids <- cbind(Ōkāreka_DHT_centroids, st_coordinates(Ōkāreka_DHT_centroids))

# Extract the geometry of the shoreline
shoreline <- st_union(Ōkāreka_DHT$geom)

ggplot() +
  #geom_sf(data = shoreline, lwd = 2.5) +
  geom_sf(data = Ōkāreka_DHT, aes(col = DHT), lwd = 1.5) +
  coord_sf() +
  scale_color_manual(values = custom_colors) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw()+
  geom_text_repel(data = Ōkāreka_DHT_centroids, aes(label = id, x = X, y = Y), size = 3, color = "black", max.overlaps = Inf)


DHT_length_Ōkāreka <- Ōkāreka_DHT %>%
  group_by(DHT) %>%
  summarize(total_length = sum(st_length(geom)))
total_length_Ōkāreka <- sum(DHT_length_Ōkāreka$total_length)
DHT_length_Ōkāreka
total_length_Ōkāreka

#### random points
# Transform to UTM Zone 60S for all geometries
DHT_length_Ōkāreka_utm <- st_transform(DHT_length_Ōkāreka, crs = 32760)  # UTM Zone 60S

# Function to generate random points for each category
set.seed(17688)
generate_random_points <- function(data, n_points = 3) {
  st_sample(data, size = n_points, type = "random")
}

# Generate random points and transform back to WGS 84
random_points_list <- lapply(split(DHT_length_Ōkāreka_utm, DHT_length_Ōkāreka_utm$DHT), function(cat_data) {
  points_utm <- generate_random_points(cat_data)
  st_transform(points_utm, crs = st_crs(DHT_length_Ōkāreka))
})

# Combine random points into one sf object with DHT attributes
all_random_points <- st_sf(
  DHT = rep(names(random_points_list), sapply(random_points_list, length)),
  geometry = st_sfc(do.call(c, random_points_list), crs = st_crs(DHT_length_Ōkāreka)))

# Extract coordinates for the random points
coordinates <- st_coordinates(all_random_points)

# Add point number as labels
label_df <- data.frame(coordinates)
label_df$label <- seq_len(nrow(label_df))
Ōkāreka_DHT$DHT <- factor(Ōkāreka_DHT$DHT, levels = names(custom_colors))

# plot
plot_Ōkāreka = ggplot() +
  geom_sf(data = Ōkāreka_DHT, aes(col = DHT), lwd = 1.5) +
  geom_sf(data = all_random_points, color = "black", size = 2, shape = 21, fill = "red") +
  geom_text_repel(data = label_df, aes(x = X, y = Y, label = label), size = 3) +  # Add labels
  coord_sf() +
  scale_color_manual(values = custom_colors) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw()

plot_Ōkāreka

################
# Combine the DHT of the 5 lakes ####
combined_DHT <- rbind(Rotorua_DHT,Rotoiti_DHT,Rotoehu_DHT,Rotomā_DHT,Ōkāreka_DHT)

combined_DHT_centroids <- st_centroid(combined_DHT)
combined_DHT_centroids <- cbind(combined_DHT_centroids, st_coordinates(combined_DHT_centroids))


# Combine all data frames into one
combined_DHT_length <- bind_rows(
  DHT_length_Rotoiti %>% mutate(Lake = "Rotoiti"),
  DHT_length_Rotoehu %>% mutate(Lake = "Rotoehu"),
  DHT_length_Rotomā %>% mutate(Lake = "Rotomā"),
  DHT_length_Ōkāreka %>% mutate(Lake = "Ōkāreka"))

# Select and rename columns
combined_DHT_length <- combined_DHT_length %>%
  select(Lake, DHT, total_length_m = total_length)

# Print the combined data frame
print(combined_DHT_length)

#custom_colors <- c("Sandy" = "orange", "Raupo" = "green", "Rocky" = "brown", "Mud/sand" = "tan", "?" = "gray", "? clif"= "black", "AR"= "red", "Geo"= "blue")

ggplot() +
  #geom_sf(data = Lake_data) +
  geom_sf(data = combined_DHT, aes(col=DHT), lwd = 1.5) +
  scale_color_manual(values = custom_colors) +
  coord_sf() +
  labs(x="Longitude", y="Latitude")+
  theme_bw()
#geom_text_repel(data = combined_DHT_centroids, aes(label = id, x = X, y = Y), size = 2, color = "black", max.overlaps = Inf)



# Calculate the combined total length
combined_DHT_length <- combined_DHT_length %>%
  group_by(Lake) %>%
  mutate(relative_length = total_length_m / sum(total_length_m)) %>%
  ungroup()

# Plot the total length of each DHT
ggplot(combined_DHT_length, aes(x = DHT, y = relative_length, fill = DHT)) +
  scale_fill_manual(values = custom_colors) +
  geom_bar(stat = "identity", col="black") +
  labs(x = "Dominant habitat type", y = "Total length") +
  facet_wrap(~ Lake)+
  theme_bw() 

# Random select 12 point along the lakes with a fixed distance -----------------

# import the outline of the 5 lakes
# File paths for the lakes
lake_files <- list(
  "Lake Rotorua" = "~/PhD/Data/3. Natural habitat monitoring/Data_raw/Lake Rotorua.gpkg",
  "Lake Rotoiti" = "~/PhD/Data/3. Natural habitat monitoring/Data_raw/Lake Rotoiti.gpkg",
  "Lake Rotoehu" = "~/PhD/Data/3. Natural habitat monitoring/Data_raw/Lake Rotoehu.gpkg",
  "Lake Rotomā" = "~/PhD/Data/3. Natural habitat monitoring/Data_raw/Lake Rotomā.gpkg",
  "Lake Ōkāreka" = "~/PhD/Data/3. Natural habitat monitoring/Data_raw/Lake Ōkareka.gpkg")

# Import and process lake geometries
Lakes <- lapply(names(lake_files), function(name) {
  st_read(lake_files[[name]]) %>%
    st_boundary() %>%
    mutate(name = name)})

# Combine all lakes into a single data frame
Lakes_combined <- do.call(rbind, Lakes)

ggplot() + geom_sf(data = Lakes_combined)


# Rotorua ----------------------------------------------------------------------
# Extract the geometry of Lake Rotorua
Rotorua <- Lakes_combined %>% filter(name == "Lake Rotorua")

# Set seed and generate random and evenly spaced points
set.seed(123)
random_point <- st_sample(Rotorua, size = 1)
shoreline_length <- st_length(Rotorua)
num_points <- 11
spacing <- shoreline_length / (num_points + 1)
distances <- seq(spacing, shoreline_length - spacing, by = spacing)
additional_points <- st_sample(Rotorua, size = num_points, type = "regular", as_points = TRUE)

# Combine all points and convert to simple feature
all_points <- c(random_point, additional_points)
all_points_sfc <- st_sfc(all_points, crs = st_crs(Rotorua))
all_points_df <- st_sf(geometry = all_points_sfc)

# Plot Lake Rotorua with all points
ggplot() +
  geom_sf(data = Rotorua, fill = NA, color = "black") +
  geom_sf(data = all_points_df, color = "red") 
  
# Print coordinates of all points
st_coordinates(all_points_df)



# Rotoiti ----------------------------------------------------------------------
# Extract the geometry of Lake Rotoiti
Rotoiti <- Lakes_combined %>% filter(name == "Lake Rotoiti")

# Set seed and generate random and evenly spaced points
set.seed(21)
random_point <- st_sample(Rotoiti, size = 1)
shoreline_length <- st_length(Rotoiti)
num_points <- 11
spacing <- shoreline_length / (num_points + 1)
distances <- seq(spacing, shoreline_length - spacing, by = spacing)
additional_points <- st_sample(Rotoiti, size = num_points, type = "regular", as_points = TRUE)

# Combine all points and convert to simple feature
all_points <- c(random_point, additional_points)
all_points_sfc <- st_sfc(all_points, crs = st_crs(Rotoiti))
all_points_df <- st_sf(geometry = all_points_sfc)

# Plot Lake Rotoiti with all points
ggplot() +
  geom_sf(data = Rotoiti, fill = NA, color = "black") +
  geom_sf(data = all_points_df, color = "red") 

# Print coordinates of all points
st_coordinates(all_points_df)

# Rotoehu ----------------------------------------------------------------------
Rotoehu <- Lakes_combined %>% filter(name == "Lake Rotoehu")

# Set seed and generate random and evenly spaced points
set.seed(21)
random_point <- st_sample(Rotoehu, size = 1)
shoreline_length <- st_length(Rotoehu)
num_points <- 11
spacing <- shoreline_length / (num_points + 1)
distances <- seq(spacing, shoreline_length - spacing, by = spacing)
additional_points <- st_sample(Rotoehu, size = num_points, type = "regular", as_points = TRUE)

# Combine all points and convert to simple feature
all_points <- c(random_point, additional_points)
all_points_sfc <- st_sfc(all_points, crs = st_crs(Rotoehu))
all_points_df <- st_sf(geometry = all_points_sfc)

# Plot Lake Rotoehu with all points
ggplot() +
  geom_sf(data = Rotoehu, fill = NA, color = "black") +
  geom_sf(data = all_points_df, color = "red") 

# Print coordinates of all points
st_coordinates(all_points_df)

# Rotomā ----------------------------------------------------------------------
Rotomā <- Lakes_combined %>% filter(name == "Lake Rotomā")

# Set seed and generate random and evenly spaced points
set.seed(21)
random_point <- st_sample(Rotomā, size = 1)
shoreline_length <- st_length(Rotomā)
num_points <- 11
spacing <- shoreline_length / (num_points + 1)
distances <- seq(spacing, shoreline_length - spacing, by = spacing)
additional_points <- st_sample(Rotomā, size = num_points, type = "regular", as_points = TRUE)

# Combine all points and convert to simple feature
all_points <- c(random_point, additional_points)
all_points_sfc <- st_sfc(all_points, crs = st_crs(Rotomā))
all_points_df <- st_sf(geometry = all_points_sfc)

# Plot Lake Rotomā with all points
ggplot() +
  geom_sf(data = Rotomā, fill = NA, color = "black") +
  geom_sf(data = all_points_df, color = "red") 

# Print coordinates of all points
st_coordinates(all_points_df)

# Ōkāreka ----------------------------------------------------------------------
Ōkāreka <- Lakes_combined %>% filter(name == "Lake Ōkareka")

# Set seed and generate random and evenly spaced points
set.seed(12)
random_point <- st_sample(Ōkāreka, size = 1)
shoreline_length <- st_length(Ōkāreka)
num_points <- 11
spacing <- shoreline_length / (num_points + 1)
distances <- seq(spacing, shoreline_length - spacing, by = spacing)
additional_points <- st_sample(Ōkāreka, size = num_points, type = "regular", as_points = TRUE)

# Combine all points and convert to simple feature
all_points <- c(random_point, additional_points)
all_points_sfc <- st_sfc(all_points, crs = st_crs(Ōkāreka))
all_points_df <- st_sf(geometry = all_points_sfc)

# Plot Lake Ōkāreka with all points
ggplot() +
  geom_sf(data = Ōkāreka, fill = NA, color = "black") +
  geom_sf(data = all_points_df, color = "red") 

# Print coordinates of all points
st_coordinates(all_points_df)









  


# the 5 lakes random point maker -----------------------------------------------
# Function to generate and plot random and evenly spaced points
# Function to generate random and evenly spaced points for a given lake
generate_points <- function(lake_data, seed, num_points = 11) {
  set.seed(seed)
  random_point <- st_sample(lake_data, size = 1)
  additional_points <- st_sample(lake_data, size = num_points, type = "regular")
  all_points <- c(random_point, additional_points)
  st_sfc(all_points, crs = st_crs(lake_data)) %>%
    st_sf(geometry = ., lake_name = unique(lake_data$name))
}

# List of lakes and seeds for reproducibility
lakes <- list(
  "Lake Rotorua" = 2,
  "Lake Rotoiti" = 2,
  "Lake Rotoehu" = 4,
  "Lake Rotomā"  = 2,
  "Lake Ōkareka" = 2
)

# Loop over lakes to generate points, plot, and show plots
for (lake_name in names(lakes)) {
  lake_data <- Lakes_combined %>% filter(name == !!lake_name)
  points_df <- generate_points(lake_data, lakes[[lake_name]])
  
  # Plot lake with generated points
  plot<-ggplot() +
    geom_sf(data = lake_data, fill = NA, color = "black") +
    geom_sf(data = points_df, color = "red") +
    ggtitle(lake_name) +
    theme_bw()
  
  # Print coordinates with lake name
  print(plot)
  print(st_coordinates(points_df) %>% 
          as.data.frame() %>% 
          mutate(lake_name = lake_name))
}

