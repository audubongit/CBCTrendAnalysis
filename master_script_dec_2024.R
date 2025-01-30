
# loadsome libraries
library(cmdstanr)
library(bbsBayes2)
library(parallel)
library(doParallel)
library(foreach)
library(tidyverse)

# define some directories
code_dir <- "C:/Users/tmeehan/Documents/GitHub/CBCTrendAnalysis"
results_dir <- "C:/Users/tmeehan/Desktop/test_results"

# set some analysis settings
number_years_per_circle_threshold <- 5 # minimum
nonzero_circles_threshold <- 3 # minimum
stratification1 <- "bbs_cws"

# input three major tables
setwd(code_dir)
species_table <- read.csv(file.path(code_dir, "/data/taxon_key_dec_2024.csv"))
count_table <- read_csv("./output/count_table_dec_2024.csv")
site_table <- read_csv("./data/site_table_dec_2024.csv")

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

# set up cluster
#n_cores <- floor((parallel::detectCores()-1)/4) # requires 4 cores per species



# for loop to be converted to parallel loop
s <- 1
#for(s in 1:nrow(species_table)){

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
  survey_reliability_s <- survey_suitabilities[s]
  
  # define output location
  dir_out1 <- file.path(results_dir, gsub(" ", "_", ebird_com_name_s))
  dir.create(dir_out1)
  
  # get data
  source("functions/get_and_filter_count_data_dec_2024.R")
  get_and_filter_count_data()
  
  # set up and run model
  source("functions/prep_and_fit_model_dec_2024.R")
  
  
#}












