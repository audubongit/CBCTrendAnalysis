library(SurveyCoverage)
library(tidyverse)
library (ebirdst)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

example_species <- "American Dipper"
set_ebirdst_access_key("52q326g4b3sf")
range_info <- grid_range(example_species, seasonal_range = "resident")



range_map <- range_info[["range_map"]]
grid <- range_info[["coverage_grid"]]
focus <- sf::st_bbox(range_map)
map <- ggplot()+
  geom_sf(data = grid,
          aes(fill = proportion_in_range))+
  geom_sf(data = range_map,
          colour = "darkorange",
          fill = NA)+
  coord_sf(xlim = focus[c("xmin","xmax")],
           ylim = focus[c("ymin","ymax")])+
  scale_fill_viridis_c()+
  labs(title = paste(example_species,"Breeding Range"),
       caption = paste(example_species,"breeding season range in the Americas.
                       Orange polygon includes",signif(range_info[["range_area"]],3),
                       "km^2.
                       The land area within the grid cells that intersect the
                       orange polygon includes", signif(range_info[["range_area_gridded"]],3),"km^2."))
print(map)

survey_sites <- read.csv(list.files(path = "./data", 
                                    pattern = "*modeled_records*", 
                                    full.names = T)) %>% 
  select(circle, lon, lat, count_year) %>% distinct()


basp_coverage <- overlay_range_data(range = range_info,
                                    survey_sites = survey_sites,
                                    sites = "circle",
                                    years = "count_year",
                                    x_coord = "lon",
                                    y_coord = "lat",
                                    crs_site_coordinates = 4326,
                                    add_survey_sites_to_range = TRUE)


