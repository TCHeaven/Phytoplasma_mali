```R
install.packages(c(
  "sf",
  "ggplot2",
  "rnaturalearth",
  "rnaturalearthdata",
  "dplyr",
  "readr",
  "ggrepel"
))

library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)
library(readr)
library(ggrepel)

# ---------------------------------------------------------
# GPS coordinates
# ---------------------------------------------------------

sites <- data.frame(
  Site = c(
    "Cagnò",
    "Barbian",
    "Schluderns",
    "Aosta",
    "Meckenheim"
  ),
  latitude = c(
    46.39502146954353,
    46.58769163277825,
    46.6641,
    45.74994662506268,
    50.62539
  ),
  longitude = c(
    11.042085235135142,
    11.51227488992364,
    10.5837,
    7.31644001133698,
    7.04392
  )
)

# Convert to spatial object
sites_sf <- st_as_sf(
  sites,
  coords = c("longitude", "latitude"),
  crs = 4326
)

# ---------------------------------------------------------
# Get countries
# ---------------------------------------------------------

world <- ne_countries(
  scale = "medium",
  returnclass = "sf"
)

countries <- world %>%
  filter(admin %in% c(
    "Italy",
    "Austria",
    "Czechia",
    "Germany",
    "Switzerland",
    "France",
    "Slovenia",
    "Luxembourg",
    "Belgium",
    "Netherlands"
  ))

# ---------------------------------------------------------
# Plot
# ---------------------------------------------------------

ggplot() +

  # Country polygons
  geom_sf(
    data = countries,
    fill = "grey90",
    color = "grey50",
    linewidth = 0.4
  ) +

  # Sampling locations
  geom_sf(
    data = sites_sf,
    size = 4,
    shape = 21,
    fill = "red",
    color = "black",
    stroke = 0.8
  ) +

  # Site labels
  geom_text_repel(
    data = sites,
    aes(
      x = longitude,
      y = latitude,
      label = Site
    ),
    size = 4,
    fontface = "bold"
  ) +

  # Map extent
  coord_sf(
    xlim = c(6, 14),
    ylim = c(44.5, 51.5),
    expand = FALSE
  ) +

  labs(
    x = "Longitude",
    y = "Latitude"
  ) +

  theme_classic() +

  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 11)
  )



```
