  

# prep and fit model function
prep_and_fit_model <- function(){
  
  # libraries
  require(patchwork)
  require(cmdstanr)
  require(bbsBayes2)
  require(tidyverse)
  
  
  # select species and load data
  species <- species_f <- ebird_com_name_s
  data_1 <- read.csv(file.path(results_dir, gsub(" ", "_", ebird_com_name_s),
                               paste0(gsub(" ", "_", ebird_com_name_s), 
                          "_zero_filled_and_filtered_modeling_data.csv"))) %>% 
    mutate(field_hours=tot_hours)
  
  # get, tweak and view strata df 
  strat_df <- data_1 %>% 
    select(strata_name,
           circles_per_stratum, nonzero_circles) %>% 
    distinct() %>% 
    mutate(non_zero = nonzero_circles / circles_per_stratum)
  strata_map <- bbsBayes2::load_map(stratify_by = stratification1) %>% 
    inner_join(., strat_df,
               by = "strata_name") %>% 
    mutate(strata_vec = as.integer(factor(strata_name))) %>% 
    arrange(strata_vec)
  nstrata <- max(strata_map$strata_vec)
  nonzeroweight <- as.numeric(strata_map$non_zero)
  
  # build neighbor relationships
  source("functions/neighbors_define_dec_2024.R")
  neighbours <- neighbours_define(strata_map,
                                  species = species,
                                  strat_indicator = "strata_vec",
                                  plot_dir = dir_out1,
                                  plot_file = "_strata_map")
  N_edges <- neighbours$N_edges
  node1 <- neighbours$node1
  node2 <- neighbours$node2
  strata_information <- strata_map %>% sf::st_set_geometry(., NULL)
  
  # join strata back to data table and drop unnecessary columns
  data_prep <- data_1 %>% 
    select(circle_code, count_year,
           longitude, latitude, number_counted, field_hours,
           scaled_effort, strata_name) %>% 
    left_join(., strata_information,
              by = c("strata_name")) %>% 
    arrange(strata_name, circle_code, count_year) %>% 
    mutate(circle_vec = as.integer(factor(circle_code)),
           year_vec = count_year - (min(count_year) - 1))
  
  # save for later
  saveRDS(data_prep, paste0("data/data_prep_", gsub(" ", "_", species), ".rds"))
  # ------------------------------------------------------------------------------
  
  
  
  # circle per strata work -------------------------------------------------------
  nsites <- max(data_prep$circle_vec)
  
  # list of site and strat combos
  sites_df <- data_prep %>% 
    select(strata_vec, circle_vec) %>% 
    distinct() %>% 
    arrange(strata_vec, circle_vec) 
  
  # number of sites in each stratum
  nsites_strata <- data_prep %>% 
    select(strata_vec, circle_vec) %>% 
    distinct() %>% 
    arrange(strata_vec, circle_vec) %>% 
    group_by(strata_vec) %>% 
    summarise(nsites = n())
  
  nsites_strata <- as.integer(nsites_strata$nsites)
  maxnsites_strata <- max(nsites_strata)
  
  # matrix of which sites are in which strata
  ste_mat <- matrix(data = 0, nrow = nstrata, ncol = maxnsites_strata)
  for(i in 1:nstrata){
    ste_mat[i, 1:nsites_strata[i]] <- sites_df[which(sites_df$strata_vec == i),
                                               "circle_vec"]
  }
  # ------------------------------------------------------------------------------
  
  
  
  # data list for stan -----------------------------------------------------------
  ncounts <- nrow(data_prep)
  nyears <- max(data_prep$year_vec)
  count <- data_prep$number_counted
  strat <- data_prep$strata_vec
  year <- data_prep$year_vec
  site <- data_prep$circle_vec
  hours <- data_prep$scaled_effort
  field_hours <- data_prep$field_hours
  
  stan_data <- list(# scalar indicators
    nsites = nsites,
    nstrata = nstrata,
    ncounts = ncounts,
    nyears = nyears,
    
    # basic data
    count = count,
    strat = strat,
    year = year,
    site = site,
    
    # spatial structure
    N_edges = N_edges,
    node1 = node1,
    node2 = node2,
    
    # effort information
    hours = hours,
    
    # ragged array information to link sites to strata
    nsites_strata = nsites_strata,
    maxnsites_strata = maxnsites_strata,
    ste_mat = ste_mat,
    
    # weights
    nonzeroweight = nonzeroweight
  )
  # ------------------------------------------------------------------------------
  
  
  
  # set up spatial first difference model ----------------------------------------
  stan_data[["N_edges"]] <- N_edges
  stan_data[["node1"]] <- node1
  stan_data[["node2"]] <- node2
  stan_data[["fixed_year"]] <- floor(stan_data$nyears / 2)
  stan_data[["zero_betas"]] <- rep(0, stan_data$nstrata)
  stan_data[["Iy1"]] <- c((stan_data$fixed_year - 1):1)
  stan_data[["nIy1"]] <- length(stan_data[["Iy1"]])
  stan_data[["Iy2"]] <- c((stan_data$fixed_year + 1):stan_data$nyears)
  stan_data[["nIy2"]] <- length(stan_data[["Iy2"]])
  
  # add additional effort vizualisation variables
  stan_data[["n_effort_preds"]] <- 100
  stan_data[["effort_preds"]] <- c(seq(from = min(stan_data$hours), 
                                       to = max(stan_data$hours),
                                       length.out = 100))
  
  # get model file
  #mod.file <- "models/first_diff_spatial_CBC_nonspatial_effort.stan"
  mod_file <- "models/first_diff_spatial_CBC.stan"
  
  # compile model
  #set_cmdstan_path(path = "D:/Users/tmeehan/Documents/.cmdstan/cmdstan-2.35.0")
  model <- cmdstan_model(mod_file, stanc_options = list("Oexperimental"))
  # ------------------------------------------------------------------------------
  
  
  
  # run stan model ---------------------------------------------------------------
  ## initial Values (not used?) 
  # init_def <- function(){ list(strata_raw = rnorm(nstrata, 0, 0.1),
  #                              STRATA = 0,
  #                              sdstrata = runif(1, 0.01, 0.1),
  #                              ste_raw = rnorm(nsites, 0, 0.1),
  #                              sdnoise = runif(1, 0.01, 0.2),
  #                              sdb = runif(1, 0.01, 0.1),
  #                              sdp = runif(1, 0.01, 0.1),
  #                              b_raw = rnorm(nstrata, 0, 0.01),
  #                              p_raw = rnorm(nstrata, 0, 0.01),
  #                              B = 0,
  #                              P = 0,
  #                              sdste = runif(1, 0.01, 0.2),
  #                              sdbeta = runif(1, 0.01, 0.1),
  #                              sdBETA = runif(1, 0.01, 0.1),
  #                              BETA_raw = rnorm(nyears - 1, 0, 0.1),
  #                              beta_raw = matrix(rnorm((nyears - 1) * 
  #                                nstrata, 0, 0.01), nrow = nstrata, 
  #                                ncol = nyears - 1))}
  
  # start sampling
  stanfit <- model$sample(
    data = stan_data,
    refresh = 200,
    chains = 4, 
    iter_sampling = 1000,
    iter_warmup = 1000,
    parallel_chains = 4,
    #pars = parms,
    adapt_delta = 0.8,
    max_treedepth = 11,
    #seed = 123,
    init = 1,
    show_exceptions = FALSE)
  
  
  # save cmdstan objects
  summ <- stanfit$summary()
  stanfit$save_object(paste0("output/fit_", gsub(" ", "_", species), "_CBC_spatial_first_diff.rds"))
  saveRDS(stan_data, paste0("output/datalist_", gsub(" ", "_", species), "_CBC_spatial_first_diff.rds"))
  saveRDS(summ, paste0("output/parameter_summary_", gsub(" ", "_", species), "_CBC_spatial_first_diff.rds"))
  
  # stanfit$save_object(paste0("output/fit_",species,"_CBC_spatial_first_diff_nonspat_eff.rds"))
  # saveRDS(stan_data, paste0("output/datalist_",species,"_CBC_spatial_first_diff_nonspat_eff.rds"))
  # saveRDS(summ, paste0("output/parameter_summary_",species,"_CBC_spatial_first_diff_nonspat_eff.rds"))
  # ------------------------------------------------------------------------------
  
  
  
  
  
  
  
}

  
