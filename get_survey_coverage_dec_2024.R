
# set up -----------------------------------------------------------------------
library(SurveyCoverage)
library(tidyverse)
library (ebirdst)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
# set_ebirdst_access_key("52q326g4b3sf")
# ------------------------------------------------------------------------------



# set parameters ---------------------------------------------------------------
species <- "American Dipper"
# ------------------------------------------------------------------------------



# get species range map and survey sites ---------------------------------------
range_info <- grid_range(species, seasonal_range = "resident")

# add study sites
survey_sites <- read.csv(list.files(path = "./data", 
                                    pattern = "*all_records*", 
                                    full.names = T)) %>% 
  select(circle, lon, lat, count_year, how_many) %>% 
  filter(how_many>0) %>% 
  distinct()
# ------------------------------------------------------------------------------



# extract and save total coverage ----------------------------------------------
spp_coverage <- overlay_range_data(range = range_info,
                                   survey_sites = survey_sites,
                                   sites = "circle",
                                   years = "count_year",
                                   x_coord = "lon",
                                   y_coord = "lat",
                                   crs_site_coordinates = 4326,
                                   add_survey_sites_to_range = TRUE)
overall_coverage_estimate <- spp_coverage$cumulative_coverage_estimate
write.csv(data.frame(species = species,
                     proportion_range_surveyed = overall_coverage_estimate[2]), 
          paste0("./output/", gsub(" ", "_", species), 
                              "_range_coverage_proportion.csv"), 
          na="", row.names=F)
# ------------------------------------------------------------------------------



# map it -----------------------------------------------------------------------
cumulative_coverage_map <- spp_coverage$cumulative_coverage_map
coverage_overall <- ggplot() +
  geom_sf(data = cumulative_coverage_map,
          aes(fill = coverage)) +
  scale_fill_viridis_d() +
  labs(title = paste(species,"proportion covered = ",
                     round(overall_coverage_estimate$coverage_proportion, 2)),
       fill = "Survey coverage") +
  theme_bw()
cairo_pdf(paste0("./output/", gsub(" ", "_", species), "_range_coverage_map.pdf"),  
          width = 11, height = 8.5)
print(coverage_overall)
dev.off()
# ------------------------------------------------------------------------------


