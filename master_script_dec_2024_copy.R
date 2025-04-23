# # first change the temp directory to something that can handle cmdstan output files
# usethis::edit_r_environ()
# add something like this: TMPDIR = "C:\Users\tmeehan\Desktop\test_data"
# .rs.restartR()
tempdir()


# load some libraries
library(SurveyCoverage)
library(ebirdst)
library(cmdstanr)
library(bbsBayes2)
library(tidyverse)

# define some directories
code_dir <- "C:/Users/tmeehan/Documents/GitHub/CBCTrendAnalysis"
results_dir <- "Z:/7_CommunityScience/CBCAnalysisResults/cbc_results_v2023.0"
# results_dir <- "C:/Users/tmeehan/Desktop/test_data"

# set some stratum selection settings
number_years_per_circle_threshold <- 5 # minimum
nonzero_circles_threshold <- 3 # minimum
stratification1 <- "bbs_cws"

# # set some stan model settings for testing
# refresh1 <- 0
# sig_figs1 <- 4
# chains1 <- 4
# iter_sampling1 <- 20
# iter_warmup1 <- 20
# parallel_chains1 <- 4
# adapt_delta1 <- 0.8
# max_treedepth1 <- 12
# init1 <- 1

# set some stan model settings for actual analysis
refresh1 <- 250
sig_figs1 <- 4
chains1 <- 4
iter_sampling1 <- 1000
iter_warmup1 <- 1000
parallel_chains1 <- 4
adapt_delta1 <- 0.8
max_treedepth1 <- 12
init1 <- 1

# input three major tables
setwd(code_dir)
species_table <- read.csv(file.path(code_dir, "/data/taxon_key_dec_2024.csv"),
                          encoding="latin1") %>% 
  # mutate(ebird_com_name=gsub("�", "", ebird_com_name)) %>% 
  # mutate(historic_cbc_com_name=gsub("�", "", historic_cbc_com_name)) %>% 
  mutate(ebird_com_name=gsub("/", " or ", ebird_com_name))
count_table <- read.csv("./output/count_table_dec_2024.csv") # %>% 
  # mutate(common_name=gsub("�", "", common_name))
site_table <- read.csv("./data/site_table_dec_2024.csv")
# Encoding(species_table$ebird_com_name) <- "UTF-8"
# Encoding(species_table$historic_cbc_com_name) <- "UTF-8"
# Encoding(count_table$common_name) <- "UTF-8"

# identify which worker to use
# worker_number <- 1
# species_table <- species_table %>% arrange(desc(total_counted)) %>%
#   mutate(worker_id=rep(1:10, length.out=nrow(species_table))) %>%
#   filter(worker_id==worker_number) %>%
#   sample_n(size=nrow(.), replace = FALSE)

# define vectors for looping
ebird_spp_codes <- species_table$ebird_spp_code
ebird_com_names <- species_table$ebird_com_name
ebird_sci_names <- species_table$ebird_sci_name
historic_cbc_com_names <- species_table$historic_cbc_com_name
lon_filters <- species_table$lon_filter
lat_filters <- species_table$lat_filter
prov_state_filters <- species_table$prov_state_filter
bcr_filters <- species_table$bcr_filter
add_nocturnals <- species_table$add_nocturnal
add_feeders <- species_table$add_feeder
survey_suitabilities <- species_table$survey_suitability

# loop through species
s <- 121
for(s in 1:nrow(species_table)){ # start for loop
  
  # define species variables
  ebird_spp_code_s <- ebird_spp_codes[s]
  ebird_com_name_s <- ebird_com_names[s]
  ebird_sci_name_s <- ebird_sci_names[s]
  historic_cbc_com_name_s <- historic_cbc_com_names[s]
  lon_filter_s <- lon_filters[s]
  lat_filter_s <- lat_filters[s]
  prov_state_filter_s <- prov_state_filters[s]
  bcr_filter_s <- bcr_filters[s]
  add_nocturnal_s <- add_nocturnals[s]
  add_feeder_s <- add_feeders[s]
  survey_suitability_s <- survey_suitabilities[s]
  
  # define output location
  dir_out1 <- file.path(results_dir, gsub(" ", "_", ebird_com_name_s))
  
  # skip if the species has been finished
  if(file.exists(paste0(dir_out1, "/", gsub(" ", "_", ebird_com_name_s), 
                        "_stratum_trend_map.pdf"))){
    next
  }

  # otherwise continue
  if(!file.exists(paste0(dir_out1, "/", gsub(" ", "_", ebird_com_name_s), 
                        "_stratum_trend_map.pdf"))){
    
    # if not create directory
    dir.create(dir_out1)
    
    # set up progress log file ---------------------------------------------------
    capture.output(print(paste0("Create progress log file. ", Sys.time())),
                   file=file.path(dir_out1, "analysis_progress_log.txt"), 
                   append=F)
    
    
    # get data -------------------------------------------------------------------
    source("functions/get_and_filter_count_data_dec_2024.R")
    tryCatch({
      get_and_filter_count_data()
    }, error=function(e){})
    
    print(" "); print(" "); print(" ") 
    print(paste("got data for", ebird_com_name_s))
    print(" "); print(" "); print(" ") 
    
    
    # get coverage ---------------------------------------------------------------
    source("functions/get_survey_coverage_dec_2024.R")
    tryCatch({
      get_survey_coverage()
      capture.output(print(paste0("Finished getting coverage. ", Sys.time())),
                     file=file.path(dir_out1, "analysis_progress_log.txt"), append=T)
    }, error=function(e){})
    
    try(unlink(file.path(ebirdst_data_dir(), "2022", ebird_spp_code_s), 
               recursive = T), silent=TRUE)
    
    print(" "); print(" "); print(" ") 
    print(paste("tried coverage for", ebird_com_name_s))
    print(" "); print(" "); print(" ") 
    
    
    # set up and run model -------------------------------------------------------
    source("functions/neighbors_define_dec_2024.R")
    source("functions/prep_and_fit_model_dec_2024.R")
    tryCatch({
      prep_and_fit_model()
      capture.output(print(paste0("Finished running model. ", Sys.time())),
                     file=file.path(dir_out1, "analysis_progress_log.txt"), append=T)
    }, error=function(e){})
    
    print(" "); print(" "); print(" ") 
    print(paste("tried model for", ebird_com_name_s))
    print(" "); print(" "); print(" ") 
    
    
    # process model results ------------------------------------------------------
    source("functions/post_processing_functions_dec_2024.R")
    source("functions/process_model_results_dec_2024.R")
    tryCatch({
      draws1 <- readRDS(file.path(dir_out1, 
                                  paste0("posterior_draws_",
                                         gsub(" ", "_", ebird_com_name_s),
                                         "_CBC_spatial_first_diff.rds")))
      process_model_results()
      rm(draws1)
      capture.output(print(paste0("Finished processing results. ", Sys.time())),
                     file=file.path(dir_out1, "analysis_progress_log.txt"), append=T)
    }, error=function(e){})
    
    print(" "); print(" "); print(" ") 
    print(paste("tried to process results for", ebird_com_name_s))
    print(" "); print(" "); print(" ") 
    
    
    # add quality info -----------------------------------------------------------
    source("functions/add_estimate_quality_dec_2024.R")
    tryCatch({
      add_estimate_quality()
      capture.output(print(paste0("Finished adding quality. ", Sys.time())),
                     file=file.path(dir_out1, "analysis_progress_log.txt"), append=T)
    }, error=function(e){})
    
    print(" "); print(" "); print(" ") 
    print(paste("tried to finish for", ebird_com_name_s))
    print(" "); print(" "); print(" ") 
  }

} # end for loop


# done. that was easy.


# # now change back the temp directory
# usethis::edit_r_environ()
# remove something like this: TMPDIR = "C:\Users\tmeehan\Desktop\test_data"
# .rs.restartR()
tempdir()






