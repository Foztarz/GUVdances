# Install packages if you haven't already:
# install.packages(c("ggplot2", "sf", "dplyr", "rnaturalearth"))
library(ggplot2)
library(leaflet)

# 1. Define the coordinates
hive_coords   = c(longitude = 8.811333, latitude = 50.806306)
garden_coords = c(longitude  = 8.809010, latitude = 50.802586) #centre


# 2. Combine the data into a data frame
locations_df <- data.frame(
  Location = c("Hive", "Garden Centre"),
  Lat      = c(hive_coords["latitude"], garden_coords["latitude"]),
  Lon      = c(hive_coords["longitude"], garden_coords["longitude"])
)

# 3. Create the interactive map
leaflet(data = locations_df) %>%
  # Add the default map tiles (e.g., OpenStreetMap)
  addTiles() %>%
  # Set the map view to center on the two points
  # (The 'Garden Centre' is slightly closer to the middle)
  setView(
    lng = mean(locations_df$Lon),
    lat = mean(locations_df$Lat),
    zoom = 14 # Adjust zoom level as needed (1-18)
  ) %>%
  # Add markers for the locations
  addMarkers(
    lng   = ~Lon,
    lat   = ~Lat,
    popup = ~Location, # Show the location name when clicking
    label = ~Location # Show the location name on hover
  ) %>%
  # (Optional: Add a line connecting them)
  addPolylines(
    lng   = locations_df$Lon,
    lat   = locations_df$Lat,
    color = "red",
    weight = 2
  )