# ===========================================
# Geographic-data preparation functions
# ===========================================

library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

prepare_geographic_features <- function() {

  sf_use_s2(FALSE)

  # Obtain geographic features
  ocean <- ne_download(scale = 10, type = 'ocean', category = 'physical', returnclass = "sf")
  lakes <- ne_download(scale = 10, type = 'lakes', category = 'physical', returnclass = "sf")
  rivers <- ne_download(scale = 10, type = 'rivers_lake_centerlines',
                        category = 'physical', returnclass = "sf")
  coastline <- ne_coastline(scale = 10, returnclass = "sf")

  # Create study region
  study_area <- st_as_sfc(st_bbox(c(xmin = -125, xmax = -103, ymin = 29, ymax = 41),
                                  crs = st_crs(4326)))

  # Crop
  ocean_crop <- st_intersection(ocean, study_area)
  lakes_crop <- st_intersection(lakes, study_area)
  rivers_crop <- st_intersection(rivers, study_area)
  coast_crop <- st_intersection(coastline, study_area)

  sf_use_s2(TRUE)

  return(list(
    ocean = ocean_crop,
    lakes = lakes_crop,
    rivers = rivers_crop,
    coast = coast_crop,
    study_area = study_area
  ))
}
