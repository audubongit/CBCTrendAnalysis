

# add estimate quality function
add_estimate_quality <- function(){
  
  # set up
  require(tidyverse)
  species <- ebird_com_name_s
  
  # get data
  idx1 <- read.csv(file.path(dir_out1, 
                             paste0(gsub(" ", "_", species), 
                                    "_relative_abundance_indices.csv")))
  trds1 <- read.csv(file.path(dir_out1, 
                              paste0(gsub(" ", "_", species), 
                                     "_relative_abundance_trends.csv")))
  cov1 <- read.csv(file.path(dir_out1, 
                             paste0(gsub(" ", "_", species), 
                                    "_range_coverage_proportion.csv")))
  
  # add quality info to indices
  idx2 <- idx1 %>% 
    mutate(estimate_precision=index_q0.975-index_q0.025,
           survey_suitability=survey_suitability_s,
           range_coverage=cov1[1,2]) %>% 
    mutate(survey_suitability=ifelse(survey_suitability=="", NA, 
                                     survey_suitability),
           range_coverage=ifelse(region_type!="continent", NA, range_coverage)) %>% 
    rename(common_name=species) %>% 
    mutate(scientific_name=ebird_sci_name_s,
           species_code=ebird_spp_code_s, .after=common_name)
  
  # add quality info to trends
  trds2 <- trds1 %>% 
    mutate(estimate_precision=trend_q0.975-trend_q0.025,
           survey_suitability=survey_suitability_s,
           range_coverage=cov1[1,2]) %>% 
    mutate(survey_suitability=ifelse(survey_suitability=="", NA, 
                                     survey_suitability),
           range_coverage=ifelse(region_type!="continent", NA, range_coverage)) %>% 
    rename(common_name=species) %>% 
    mutate(scientific_name=ebird_sci_name_s,
           species_code=ebird_spp_code_s, .after=common_name)
  
  # rewrite tables
  write.csv(idx2, file.path(dir_out1, paste0(gsub(" ", "_", species), 
                                                    "_relative_abundance_indices.csv")), 
            na="", row.names=F)
  write.csv(trds2, file.path(dir_out1, paste0(gsub(" ", "_", species), 
                                                   "_relative_abundance_trends.csv")),
            na="", row.names=F)
  
} # end function


# test function
# add_estimate_quality()

