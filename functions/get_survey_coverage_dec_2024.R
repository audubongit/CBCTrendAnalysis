
get_survey_coverage <- function(){

  # set up ---------------------------------------------------------------------
  library(SurveyCoverage)
  library(tidyverse)
  library (ebirdst)
  # set_ebirdst_access_key("52q326g4b3sf")
  # ----------------------------------------------------------------------------
  
  
  
  # set parameters -------------------------------------------------------------
  species <- ebird_com_name_s
  # ----------------------------------------------------------------------------
  
  
  
  # get species range map and survey sites -------------------------------------
  range_info <- try(grid_range(species, seasonal_range = "nonbreeding"),
                    silent = TRUE)
  if(methods::is(range_info,"try-error")){
    range_info <- try(grid_range(species, seasonal_range = "resident"),
                      silent = TRUE)
  }
  
  # add study sites
  survey_sites <- read.csv(list.files(path = dir_out1, 
                                      pattern = "*filtered_modeling_data*", 
                                      full.names = T)) %>% 
    select(circle_code, longitude, latitude, count_year, number_counted) %>% 
    distinct()
  # ----------------------------------------------------------------------------
  
  
  
  # extract and save total coverage --------------------------------------------
  spp_coverage <- overlay_range_data(range = range_info,
                                     survey_sites = survey_sites,
                                     sites = "circle_code",
                                     years = "count_year",
                                     x_coord = "longitude",
                                     y_coord = "latitude",
                                     crs_site_coordinates = 4326,
                                     add_survey_sites_to_range = TRUE)
  overall_coverage_estimate <- spp_coverage$cumulative_coverage_estimate
  write.csv(data.frame(species = species,
                       proportion_range_surveyed = overall_coverage_estimate[2]), 
            file.path(dir_out1, paste0(gsub(" ", "_", species), 
                                "_range_coverage_proportion.csv")), 
            na="", row.names=F)
  # ----------------------------------------------------------------------------
  
  
  
  # map it ---------------------------------------------------------------------
  cumulative_coverage_map <- spp_coverage$cumulative_coverage_map
  coverage_overall <- ggplot() +
    geom_sf(data = cumulative_coverage_map,
            aes(fill = coverage)) +
    scale_fill_viridis_d() +
    labs(title = paste(species,"proportion covered = ",
                       round(overall_coverage_estimate$coverage_proportion, 2)),
         fill = "Survey coverage") +
    theme_bw()
  cairo_pdf(file.path(dir_out1, paste0(gsub(" ", "_", species), "_range_coverage_map.pdf")),  
            width = 11, height = 8.5)
  print(coverage_overall)
  dev.off()
  # ----------------------------------------------------------------------------
  
  
  # delete rasters
  unlink(file.path(ebirdst_data_dir(), "2022", ebird_spp_code_s), recursive = T)
  
  
} # end function


# test function
# get_survey_coverage()