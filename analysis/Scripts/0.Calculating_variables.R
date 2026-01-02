# 0. Calculate Shear stress
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
packages <- c("tmap", "waver","sf","tidyverse", "dplyr", "ggplot2", "readxl", "writexl", "readr")#"bathytools","rgdal",

# Load packages if not already installed
lapply(packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, dependencies = TRUE)
  library(pkg, character.only = TRUE)})


# Import the data sets ---------------------------------------------------------
### Monitoring sites: -------------------------------------------------------------
All_sample_sites <- read_excel("~/PhD/Data/Lakes/All_lakes/Lakes_Monitoring.xlsx", 
                               sheet = "Monitoring_Sites")
# Create an sf object for All_sample_sites
All_sample_sites <- All_sample_sites %>% filter(!is.na(Lon))
All_sample_sites_sf <- st_as_sf(All_sample_sites, coords = c("Lon", "Lat"), crs = 4326)  # Assuming WGS84 for lat/lon coords
All_sample_sites_sf <- st_transform(All_sample_sites_sf, crs = 2193)  # Reproject to NZGD2000

# Filter sites for each lake
Sites_Rotorua <- All_sample_sites_sf %>% filter(Lake == "Lake Rotorua")
Sites_Rotoiti <- All_sample_sites_sf %>% filter(Lake == "Lake Rotoiti")
Sites_Rotoehu <- All_sample_sites_sf %>% filter(Lake == "Lake Rotoehu")
Sites_Rotomā  <- All_sample_sites_sf %>% filter(Lake == "Lake Rotomā")
Sites_Ōkāreka <- All_sample_sites_sf %>% filter(Lake == "Lake Ōkāreka")


### Lake perimeter: ------------------------------------------------------------
Rotorua <- st_read("Data_raw/Lakes_perimeter/Lake Rotorua.gpkg")
Rotoiti <- st_read("Data_raw/Lakes_perimeter/Lake Rotoiti.gpkg")
Rotoehu <- st_read("Data_raw/Lakes_perimeter/Lake Rotoehu.gpkg")
Rotomā  <- st_read("Data_raw/Lakes_perimeter/Lake Rotomā.gpkg")
Ōkāreka <- st_read("Data_raw/Lakes_perimeter/Lake Ōkāreka.gpkg")

# Extract boundaries
Rotorua_outline <- st_boundary(Rotorua)
Rotoiti_outline <- st_boundary(Rotoiti)
Rotoehu_outline <- st_boundary(Rotoehu)
Rotomā_outline  <- st_boundary(Rotomā)
Ōkāreka_outline <- st_boundary(Ōkāreka)

Rotorua_outline <- st_transform(Rotorua_outline, crs = 2193)
Rotoiti_outline <- st_transform(Rotoiti_outline, crs = 2193)
Rotoehu_outline <- st_transform(Rotoehu_outline, crs = 2193)
Rotomā_outline  <- st_transform(Rotomā_outline, crs = 2193)
Ōkāreka_outline <- st_transform(Ōkāreka_outline, crs = 2193)

combined_outlines <- bind_rows(Rotorua_outline, Rotoiti_outline, Rotoehu_outline, Rotomā_outline, Ōkāreka_outline)
combined_outlines <- st_transform(combined_outlines, crs = 2193)





### Lake Bathmetry: ------------------------------------------------------------
Rotlakes_bathymetry <- read_excel("~/PhD/Data/Lakes/Lakes_waterquality/Data_raw/Rotlakes_bathymetry.xls")

Bathymetry_Rotoiti <- read_excel("Data_raw/Bathymetry_Rotoiti.xlsx")

Rotoiti_bath_25_grid <- read_csv("~/PhD/0. Thesis/1. Reefs Rotoiti/Rotoiti Model/Extra/Bathmetry/Raster25.csv")

Rotoiti_bath_25 <- Rotoiti_bath_25_grid %>%
  pivot_longer(cols = -Y, names_to = "X", values_to = "Z") %>%
  mutate(Z = round(Z), Z = as.integer(Z)) %>% 
  mutate(X = as.numeric(X), Y = as.numeric(Y),Z = as.numeric(Z))%>%
  filter(!is.na(Z))

# Create an sf object for Rotoiti_bath_25
Rotoiti_bath_25_sf <- st_as_sf(Rotoiti_bath_25, coords = c("X", "Y"), crs = 27200)

ggplot()+
  geom_sf(data=Rotoiti_bath_25_sf)


### Lake Bathmetry Lernzmp: ----------------------------------------------------

#Please find attached a zip file with the lake bathymetrys as ".tif" files. 
#These can be read in to R using the "terra" package. 
#There are functions within that package to calculate slope/aspect/hillshade etc. 
#They are all in the CRS: NZGD2000 are at 2m resolution. 
#Would also recommend the "tmap" package for generating nice visualisations. 

# Define the list of packages
packages <- c("sf","tmap","terra", "readr", "readxl", "tidyverse", "dplyr", "ggplot2")

# Load packages if not already installed
lapply(packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, dependencies = TRUE)
  library(pkg, character.only = TRUE)})

# Import the data  
folder_path <- "Data_raw/bathy_oraven_Jan2025"

tif_files <- list.files(folder_path, pattern = "\\.tif$", full.names = TRUE)


# Load rasters into a list and rename them based on the filename
lake_rasters <- lapply(tif_files, function(file) {
  lake_name <- tools::file_path_sans_ext(basename(file))  # Extract filename without extension
  r <- rast(file)  # Load raster
  names(r) <- lake_name  # Rename the raster layer
  return(r)})

# Name the list elements by lake name
names(lake_rasters) <- sapply(tif_files, function(file) tools::file_path_sans_ext(basename(file)))

plot(lake_rasters[[5]])  


# Loop through the first 5 lakes and plot each raster with black contour lines
for (i in 1:5) {
  r <- lake_rasters[[i]]
  plot(r, main = paste("Lake", names(lake_rasters)[i]), col = terrain.colors(100))
  contour(r, add = TRUE, col = "black")}



# Create data frames for each of the lakes in the list
lake_data_frames <- lapply(names(lake_rasters), function(lake_name) {
  r <- lake_rasters[[lake_name]]
  r_df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)  
  r_df$lake <- lake_name  
  
  return(r_df)
})
names(lake_data_frames) <- names(lake_rasters)
















### Meteorology Data: ----------------------------------------------------------
Meteorology <- read.csv("C:/1.PhD stuff/AEM3D/Rotoiti_wall_125/run/infiles/Meteorology_sans_cloud.dat", 
                        header = TRUE, skip = 14, sep = "", stringsAsFactors = FALSE, 
                        colClasses = c("character", rep("numeric", 7)))

# Split the TIME column into date and fractional day
Meteorology$DATE <- as.Date(as.character(floor(as.numeric(Meteorology$TIME))), format = "%Y%j")
Meteorology$TIME_OF_DAY <- as.numeric(Meteorology$TIME) %% 1 * 24


# Calculate Shear stress -------------------------------------------------------

# Load geospatial data and sample sites 
lake_boundary <- Ōkāreka_outline
sites         <- Sites_Ōkāreka

ggplot()+
  geom_sf(data=Ōkāreka_outline)+
  geom_sf(data=Sites_Ōkāreka)

# Set wind direction 
# Count the frequency of each wind direction
wind_direction_count <- Meteorology %>%
  count(WIND_DIR) %>%
  arrange(desc(n))

# Get the dominant wind direction (the one with the highest count)
dominant_wind_dir <- wind_direction_count[2, ]

wind_direction <- 90  

# Convert wind direction to radians
wind_angle_rad <- (wind_direction - 90) * pi / 180  # Convert to radians

# Define the wind vector (direction of the wind)
wind_vector <- c(cos(wind_angle_rad), sin(wind_angle_rad))

# Ensure both datasets are in the same CRS (e.g., NZGD2000 / New Zealand Transverse Mercator 2000)
lake_boundary <- st_transform(lake_boundary, crs = st_crs(sites))

# Extract coordinates from lake boundary and sites
lake_boundary_coords <- st_coordinates(lake_boundary)
site_coords <- st_coordinates(sites)

# Function to calculate the fetch length for each site
calculate_fetch_for_site <- function(site_coord, lake_boundary_coords, wind_vector) {
  # For each site, calculate the projection of all lake boundary points onto the wind direction
  projections <- apply(lake_boundary_coords, 1, function(coord) {
    return(sum(coord[1:2] * wind_vector))  # Dot product for projection
  })
  
  # Calculate the projection of the site onto the wind direction
  site_projection <- sum(site_coord[1:2] * wind_vector)  # Dot product for site
  
  # Find the minimum distance from the site to the boundary along the wind direction
  fetch_distance <- max(projections) - site_projection  # Fetch length in the wind direction
  
  return(fetch_distance)
}

# Apply the fetch calculation for each site
fetch_lengths <- apply(site_coords, 1, function(site_coord) {
  calculate_fetch_for_site(site_coord, lake_boundary_coords, wind_vector)
})

# Add the fetch lengths to the `Sites_Ōkāreka` dataset
sites$fetch_length <- fetch_lengths

# View the results
print(sites[, c("Site_ID", "fetch_length")])







# Load wind data
wind_data <- read.csv("path/to/your/wind_data.csv")

# Calculate wave energy flux
mean_water_depth <- 4
wave_energy_flux <- wave_energy(wind_data, fetch_lengths, mean_water_depth)

# Log transform and aggregate wave energy flux
wave_energy_flux$log_energy <- log(wave_energy_flux$energy)
mean_log_energy <- wave_energy_flux %>%
  group_by(location) %>%
  summarize(mean_log_energy = mean(log_energy, na.rm = TRUE))

# Print the results
print(mean_log_energy)

# Slope ------------------------------------------------------------------------

#Required Data:
#  Digital Elevation Model (DEM):
#  A raster dataset that represents the elevation of the terrain.
# Shapefile of Shoreline:
#  A polygon or line shapefile representing the shoreline of the lake.

library(raster)
library(rgdal)
library(sf)
library(rgeos)


# Distance to Deep 20m ---------------------------------------------------------





# Rotoiti
ggplot()+
  #geom_sf(data=Rotoiti_bath_25_sf)+
  geom_sf(data=Rotoiti_outline)+
  geom_sf(data=Sites_Rotoiti)



# Merge Bathymetry with DEM data -----------------------------------------------
install.packages("bathytools")
library(bathytools)
library(tmap) # Used for plotting spatial data

shoreline <- readRDS(system.file("extdata/rotoma_shoreline.rds",
                                 package = "bathytools"))
catchment <- readRDS(system.file("extdata/rotoma_catchment.rds",
                                 package = "bathytools"))
point_data <- readRDS(system.file("extdata/depth_points.rds",
                                  package = "bathytools"))


