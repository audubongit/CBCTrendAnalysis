
get_survey_coverage <- function(){

  # set up ---------------------------------------------------------------------
  require(SurveyCoverage)
  require(tidyverse)
  require (ebirdst)
  # set_ebirdst_access_key("52q326g4b3sf", overwrite=T)
  # ----------------------------------------------------------------------------
  
  
  
  # set parameters -------------------------------------------------------------
  species <- ebird_com_name_s
  spec <- species
  if(grepl(" or ", spec)==TRUE){
    spec <- unlist(str_split(unlist(str_split(spec, " or "))[2], " [(]"))[1]
  }
  # species_l <- gsub(" ", "_", species)
  # ----------------------------------------------------------------------------
  
  
  # if there is an ebird range map
  if(spec %in% ebirdst_runs$common_name==T){
  
    # get species range map and survey sites -------------------------------------
    season1 <- as.numeric(is.na(ebirdst_runs[
      which(ebirdst_runs$common_name==spec), "resident_start"]))
    range_season <- ifelse(season1==1, "nonbreeding", "resident")
    range_info <- grid_range(spec, seasonal_range = range_season)
    
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
    ggsave(file.path(dir_out1, paste0(gsub(" ", "_", species), "_range_coverage_map.pdf")), 
           coverage_overall, width = 11, height = 8.5, units="in")
    # ----------------------------------------------------------------------------
    
    
    # delete rasters
    unlink(file.path(ebirdst_data_dir(), "2022", ebird_spp_code_s), recursive = T)
  }
  
  
  # if no ebird range map ----------------------------------------------------------
  if(spec %in% ebirdst_runs$common_name==F){
    write.csv(data.frame(species = species,
                       proportion_range_surveyed = NA), 
            file.path(dir_out1, paste0(gsub(" ", "_", species), 
                                       "_range_coverage_proportion.csv")), 
            na="", row.names=F)
  }
  # -------------------------------------------------------------------------------
  
} # end function


# test function
# get_survey_coverage()