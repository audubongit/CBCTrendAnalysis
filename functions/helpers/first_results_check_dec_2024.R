


# Script for checking the quality of results from the first model run.
# Tim Meehan, Quantitative Ecologist, National Audubon Society
# 21 April 2025 



# set up -----------------------------------------------------------------------
library(tidyverse)

# set paths for files and results
git_data_path <- "C:/Users/tmeehan/Documents/GitHub/CBCTrendAnalysis"
results_root_directory <- "Z:/7_CommunityScience/CBCAnalysisResults/cbc_results_v2023.0"
species_list_file <- paste(git_data_path, "data/taxon_key_dec_2024.csv", 
                           sep="/")
options(stringsAsFactors=F)
setwd(results_root_directory)
# ------------------------------------------------------------------------------



# check results per species ----------------------------------------------------
# get species file
species_list <- read.csv(species_list_file) %>% 
  mutate(ebird_com_name=gsub("Â", "", ebird_com_name)) %>% 
  mutate(ebird_com_name=gsub("/", " or ", ebird_com_name))

# start data frame for results
res1 <- data.frame(ebird_com_name=species_list[,2], 
                   dir_exists=0, 
                   got_data=0, 
                   got_coverage=0, 
                   mapped_neighbors=0,
                   got_draws = 0,
                   got_par_sums = 0,
                   made_trend_maps = 0, 
                   got_indices = 0,
                   got_trends = 0, 
                   prop_bad_rh = NA,
                   prop_bad_ess = NA,
                   trend_median = NA,
                   trend_q0.025 = NA,
                   trend_q0.975 = NA)

# loop through species
i <- 1
for(i in 1:nrow(species_list)){
  
  # set species code
  species_code <- gsub(" ", "_", species_list[i, 2])
  
  # set wd
  setwd(results_root_directory)
  
  # status
  print(paste("trying", species_code))
  
  # if dir exists
  if(dir.exists(species_code)) {
    
    # verify directory
    res1[i, 2] <- 1
    
    # go to that directory
    setwd(file.path(results_root_directory, species_code))
    
    # check existence of data
    if(file.exists(paste0(species_code, 
                          "_zero_filled_and_filtered_modeling_data.csv")) == 
       TRUE) {
      res1[i, 3] <- 1
    }
    
    # check existence of coverage proportions
    if(file.exists(paste0(species_code, "_range_coverage_proportion.csv")) == 
       TRUE) {
      res1[i, 4] <- 1
    }
    
    # check existence of neighborhood map
    if(file.exists(paste0(species_code, "_strata_map.pdf")) == 
       TRUE) {
      res1[i, 5] <- 1
    }
    
    # check existence of posterior draws
    if(file.exists(paste0("posterior_draws_", species_code, 
                          "_CBC_spatial_first_diff.rds")) == 
       TRUE) {
      res1[i, 6] <- 1
    }
    
    # check existence of parameter summaries
    if(file.exists(paste0(species_code, "_parameter_estimate_summaries.csv")) == 
       TRUE) {
      res1[i, 7] <- 1
    }
    
    # check existence of stratum trend maps
    if(file.exists(paste0(species_code, "_stratum_trend_map.pdf")) == 
       TRUE) {
      res1[i, 8] <- 1
    }
    
    # check existence of abundance indices
    if(file.exists(paste0(species_code, "_relative_abundance_indices.csv")) == 
       TRUE) {
      res1[i, 9] <- 1
    }
    
    # check existence of abundance trends
    if(file.exists(paste0(species_code, "_relative_abundance_trends.csv")) == 
       TRUE) {
      res1[i, 10] <- 1
    }
    
    # look in parameter summary file
    if(file.exists(paste0(species_code, "_parameter_estimate_summaries.csv")) == 
       TRUE) {
      
      # load error file
      fname <- paste(results_root_directory, species_code, 
                     paste0(species_code,
                     "_parameter_estimate_summaries.csv"), 
                     sep="/")
      log_file <- read.csv(fname)
      
      # get rhat
      rh1 <- log_file$rhat
      prop_bad_rh <- sum(as.numeric(rh1 > 1.05)) / length(rh1)
      res1[i, 11] <- prop_bad_rh
      
      # get ess
      ess1 <- log_file$ess_bulk
      prop_bad_ess <- sum(as.numeric(ess1 < 200)) / length(ess1)
      res1[i, 12] <- prop_bad_ess
    }
    
    # grab abundance trend
    if(file.exists(paste0(species_code, "_relative_abundance_trends.csv")) == 
       TRUE) {
      
      trd_file <- read.csv(paste0(species_code, 
                                  "_relative_abundance_trends.csv"))
      trds <- trd_file %>% 
        filter(region_type=="continent", trend_type=="long-term") %>% 
        select(trend_median, trend_q0.025, trend_q0.975) %>% 
        as.numeric()
      res1[i, 13:15] <- trds

    }
    
  } else {next}
  
}
res1$sum_success <- rowSums(res1[,3:8])
# ------------------------------------------------------------------------------



# save output tables -----------------------------------------------------------
write.csv(res1, 
          paste0(git_data_path, "/data/", "first_results_check.csv"), 
          row.names=F)
# ------------------------------------------------------------------------------



