


# cbc annual index function ----------------------------------------------------
index_function <- function(cmdstan_draws = draws1,
                           quant = 0.95, # credible interval
                           year_1 = 1966, # first year
                           parameter = "n", # parameter name from samples
                           year = "Year", # not sure
                           first_dim = "s", # "s" stratum, "y" year, for samples
                           # if weights_df supplied, composite summary estimated
                           weights_df = NULL, # or "cntry_wts"
                           strat = NULL, # or "stratum", match weights_df
                           area = NULL, # or "Area", match weights_df
                           # optional summary_regions column to get region summaries 
                           # summarizes across all strata if weights_df is supplied
                           summary_regions = NULL, # or "country", match weights_df
                           # to_summarise is TRUE, summarizes indices across regions
                           to_summarise = TRUE # or FALSE, just strata estimates
){
  
  # set credible limits
  lu <- ((1-(quant))/2)
  uu <- 1-((1-(quant))/2)

  ############################# checks and tweaks ##############################
  
  if(is.null(strat)){ # when strat IS null 
    dims <- year 
  } else { # when strat NOT null
    
    # and when weights ARE given
    if(!is.null(weights_df)){ 
      if(is.null(strat)){
        stop("Weights were supplied but no stratum dimension")
        return(NULL)
      }
      
      if(!any(grepl(pattern = strat, x = names(weights_df)))){
        stop("weights_df must include a column with a name that matches strat")
        return(NULL)
      }
      
      if(!to_summarise){
        to_summarise <- TRUE
        warning("weights were supplied but to_summarise is FALSE, changing to_summaries to TRUE")
      }
    }
    
    # and when to_summarise IS true
    if(to_summarise){
      if(is.null(weights_df) | is.null(strat)){
        stop("to_summarise is TRUE but no weights were supplied or strat is NULL")
        return(NULL)
      }
    }
    
    # and when first dim is s
    if(tolower(substr(first_dim, 1, 1)) == "s"){
      dims <- c(strat, year) # set dims to strat then year
    }
    
    # and when first dim is y
    if(tolower(substr(first_dim, 1, 1)) == "y"){
      dims <- c(year, strat) # set dims to year then strat
    }
    
  } #### end checks and tweaks
  
  # get initial posterior samples following checks and tweaks
  smpls <- posterior_samples(cmdstan_draws = draws1, parm = parameter, dims = dims)
  
  #################### processing when weights ARE given #######################
  
  if(!is.null(weights_df)){ # weights ARE given
    
    # check for lengths
    nstrata_fit <- length(unique(smpls[[strat]]))
    nstrata_w <- nrow(weights_df)
    if(nstrata_fit != nstrata_w){ 
      stop("Lengths of strata and weights are different")
      return(NULL)
    }
    
    if(!is.null(summary_regions)){ # and summary regions ARE given
      
      # sum weights across regions
      weights_df <- weights_df %>% 
        rename_with(., ~gsub(pattern = area, replacement = "stratum_area", 
                             x = .x, fixed = TRUE)) %>% 
        rename_with(., ~gsub(pattern = summary_regions, 
                             replacement = "summary_region", .x, fixed = TRUE))
      tmp <- weights_df %>% 
        group_by(summary_region) %>% 
        summarise(regional_area_sum = sum(stratum_area), .groups = "keep")
      weights_df <- weights_df %>% 
        left_join(., tmp, by = "summary_region") %>% 
        mutate(stratum_weight = stratum_area / regional_area_sum)
      
    } else { # and summary regions NOT given
      
      # teak weights
      weights_df <- weights_df %>% 
        rename_with(., ~gsub(pattern = area, replacement = "stratum_area", 
                             x = .x, fixed = TRUE)) %>% 
        mutate(stratum_weight = stratum_area / sum(stratum_area),
               summary_region = "Survey_wide") 
      
      # set summary regions
      summary_regions <- "Survey_wide"
      
    }
    
    # join samples and weights when weights ARE given
    smpls <- smpls %>% 
      left_join(., weights_df, by = strat) %>% 
      ungroup()
    
    # calc indices when weights ARE given
    inds <- smpls %>% 
      mutate(.value = .value * stratum_weight) %>% 
      rename_with(., ~gsub(pattern = year, replacement = "yyy", .x, 
                           fixed = TRUE)) %>% 
      group_by(yyy, summary_region, .draw) %>% 
      summarise(.vsum = sum(.value), .groups = "keep") %>% 
      group_by(yyy, summary_region) %>% 
      summarise(mean = mean(.vsum),
                median = median(.vsum),
                lci = quantile(.vsum, lu),
                uci = quantile(.vsum, uu),
                q0.05 = quantile(.vsum, probs = 0.05),
                q0.95 = quantile(.vsum, probs = 0.95),
                q0.165 = quantile(.vsum, probs = 0.165),
                q0.835 = quantile(.vsum, probs = 0.835),
                .groups = "keep") %>% 
      mutate(true_year = yyy + year_1 - 1) %>% 
      ungroup() %>% 
      rename_with(., ~gsub(replacement = summary_regions, 
                           pattern = "summary_region", .x, fixed = TRUE)) %>% 
      rename_with(., ~gsub(replacement = year,pattern = "yyy", .x, 
                           fixed = TRUE))
    
    # calc samples when weights ARE given
    smpls <- smpls %>% 
      mutate(.value = .value * stratum_weight) %>% 
      rename_with(., ~gsub(pattern = year, replacement = "yyy", .x,
                           fixed = TRUE)) %>% 
      group_by(yyy, summary_region, .draw) %>% 
      summarise(.value = sum(.value), .groups = "keep") %>%
      rename_with(., ~gsub(pattern = year, replacement = "yyy", .x,
                           fixed = TRUE)) %>% 
      rename_with(., ~gsub(replacement = summary_regions,
                           pattern = "summary_region", .x, fixed = TRUE)) %>% 
      mutate(true_year = yyy + year_1 - 1) %>% 
      rename_with(., ~gsub(pattern = "yyy", replacement = year, .x,
                           fixed = TRUE)) %>% 
      ungroup()
    
    # calc weights for output when weights ARE given
    weights_df <- weights_df %>% 
      rename_with(., ~gsub(replacement = area, pattern = "stratum_area", .x,
                           fixed = TRUE)) %>% 
      rename_with(., ~gsub(replacement = summary_regions, 
                           pattern = "summary_region", .x, fixed = TRUE))
    
  } else { #### end processing when weights ARE given
    
  ##################### processing when weights NOT given ######################
    
    if(!is.null(strat)){ # and strata ARE given
      
      # calc indices for strata when weights NOT given
      inds <- smpls %>% 
        rename_with(., ~gsub(pattern = strat, replacement = "sss", .x,
                             fixed = TRUE)) %>% 
        rename_with(., ~gsub(pattern = year, replacement = "yyy", .x,
                             fixed = TRUE)) %>% 
        group_by(sss, yyy) %>% 
        summarise(mean = mean(.value),
                  median = median(.value),
                  lci = quantile(.value, lu),
                  uci = quantile(.value, uu),
                  q0.05 = quantile(.value, probs = 0.05),
                  q0.95 = quantile(.value, probs = 0.95),
                  q0.165 = quantile(.value, probs = 0.165),
                  q0.835 = quantile(.value, probs = 0.835),
                  .groups = "keep") %>%
        ungroup() %>% 
        mutate(true_year = yyy + year_1 - 1) %>% 
        rename_with(., ~gsub(replacement = strat, pattern = "sss", .x,
                             fixed = TRUE)) %>% 
        rename_with(., ~gsub(replacement = year, pattern = "yyy", .x,
                             fixed = TRUE)) %>% 
        mutate(region_type = "stratum") %>% 
        select(Year, stratum, mean, median, lci, uci, 
               q0.05, q0.95, q0.165, q0.835, 
               true_year, region_type)
      
      # assign summary region when weights NOT given
      summary_regions <- strat
      
      # calc samples for strata when weights NOT given
      smpls <- smpls %>% 
        rename_with(., ~gsub(pattern = year, replacement = "yyy", .x,
                             fixed = TRUE)) %>% 
        rename_with(., ~gsub(replacement = strat, pattern = "sss", .x,
                             fixed = TRUE)) %>% 
        mutate(true_year = yyy + year_1 - 1) %>% 
        rename_with(., ~gsub(pattern = "yyy", replacement = year, .x,
                             fixed = TRUE)) %>% 
        select(Year, stratum, .draw, .value, true_year)
      
    } 
  } #### end processing when weights NOT given
  
  # set final summary region
  inds <- inds %>% 
    mutate(region_type = summary_regions)
  
  # return index list
  return(list(indices = inds,
              samples = smpls,
              parameter = parameter,
              strat = strat,
              year = year,
              dims = dims,
              quant = quant,
              weights_df = weights_df,
              area = area,
              summary_regions = summary_regions,
              to_summarise = to_summarise,
              year_1 = year_1
  ))
}
# ------------------------------------------------------------------------------





# cbc smoothed index function --------------------------------------------------
get_smooth_ns <- function(idx_list){
  
  # establish quantiles
  lu <- ((1-(as.numeric(idx_list$quant)))/2)
  uu <- 1-((1-(as.numeric(idx_list$quant)))/2)
  
  # get n samples and smooth n
  dat_i <- idx_list$samples %>% # ungroup() %>% 
    group_by(across(c(2, 3))) %>% 
    mutate(.value = smooth_sample(x=Year, y=.value))

  # make summaries of smooth n
  sum_smooth <- dat_i %>% # ungroup() %>% 
    group_by(across(c(1, 2, 5))) %>% 
    summarise(mean = mean(.value),
              median = median(.value),
              lci = quantile(.value, lu),
              uci = quantile(.value, uu),
              q0.05 = quantile(.value, probs = 0.05),
              q0.95 = quantile(.value, probs = 0.95),
              q0.165 = quantile(.value, probs = 0.165),
              q0.835 = quantile(.value, probs = 0.835)) %>% 
    select(Year, 2, mean, median, lci, uci, q0.05, q0.95, q0.165, q0.835, 
           true_year) %>% 
    mutate(region_type=idx_list$summary_regions)
  
  # return list with two dataframes
  return(list(samples_smooth=dat_i, indices_smooth=sum_smooth))
}
# ------------------------------------------------------------------------------





# cbc trends mapping function --------------------------------------------------
map_function <- function(trds = bcr_lt_trds, spunit = "bcr"){
  require(ggplot2)
  require(ggpattern)
  require(bbsBayes2)
  require(sf)
  require(dplyr)
  species <- ebird_com_name_s
  breaks1 <- c(-7, -4, -2, -1, -0.5, 0.5, 1, 2, 4, 7)
  labls1 <- c(paste0("< ", breaks1[1]),
              paste0(breaks1[-c(length(breaks1))],":", breaks1[-c(1)]),
              paste0("> ",breaks1[length(breaks1)]))
  labls1 <- paste0(labls1, "%")
  cols1 <- c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", 
             "#ffffbf", "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", 
             "#313695")
  map1 <- bbsBayes2::load_map(stratify_by = spunit) %>% 
    mutate(stratum = strata_name)
  if(spunit == "bcr"){
    map1 <- bbsBayes2::load_map(stratify_by = spunit, type = "strata") %>% 
      mutate(stratum = as.numeric(str_extract(strata_name, pattern = "(\\d)+")))
  }
  if(spunit == "prov_state"){
    map1 <- bbsBayes2::load_map(stratify_by = spunit, type = "strata") %>% 
      mutate(stratum = strata_name)
  }
  map_ext1 <- sf::st_bbox(right_join(map1, trds)) * 1
  trds_map1 <- left_join(map1, trds, by = "strata_name") %>% 
    mutate(t_plot = cut(trend, breaks = c(-Inf, breaks1, Inf),
                        labels = labls1, ordered_result = TRUE)) %>% 
    mutate(t_plot = factor(t_plot, levels = labls1)) %>% 
    mutate(sig = as.character(sig)) %>% 
    mutate(sig = ifelse(is.na(sig), "0", sig))
  pal1 <- setNames(cols1, levels(trds_map1$t_plot))
  title1 <- paste0(species, ": ", 
                   unique(na.omit(trds_map1$start_year)), " through ", 
                   unique(na.omit(trds_map1$end_year)))
  tp <- ggplot() + 
    geom_sf(data = map1, fill = "grey80") + 
    geom_sf(data = trds_map1 %>% filter(!is.na(t_plot)), 
            aes(fill = t_plot), show.legend = TRUE) + 
    # geom_sf_pattern(data = trds_map1,
    #                 aes(fill = t_plot, pattern_type = sig),
    #                 pattern = 'magick',
    #                 pattern_scale = 0.5,
    #                 pattern_fill = "black",
    #                 show.legend = F) +
    #scale_pattern_type_discrete(choices = c("gray100", "left45")) +
    scale_fill_manual(values = pal1,
                      drop = F,
                      na.translate = F,
                      na.value = "grey80",
                      guide = guide_legend(title = "Trend (%/yr)", 
                                           reverse = TRUE)) +
    coord_sf(xlim = map_ext1[c("xmin", "xmax")],
             ylim = map_ext1[c("ymin", "ymax")]) +
    labs(title = title1)+
    theme_bw() +
    guides(pattern_type = "none") + 
    theme(plot.margin = unit(rep(1, 4),"mm"),
          axis.text = element_text(size = 8),
          legend.text = element_text(size = 8),
          title = element_text(size = 11),
          plot.title = element_text(hjust = 0.5))
  return(tp)
}
# ------------------------------------------------------------------------------





# cbc trends function ----------------------------------------------------------
trends_function <- function(ind_list,
                            start_year = 1966,
                            end_year = 2019,
                            quant = 0.95){
  samples <- ind_list$samples
  samples_smooth <- ind_list$samples_smooth
  year <- ind_list$year
  summary_regions <- ind_list$summary_regions

  # fill in null years
  if(is.null(end_year)){
    end_year <- max(samples$true_year)
  }
  if(is.null(start_year)){
    start_year <- min(samples$true_year)
  }
  # get period and quantiles
  nyrs <- end_year-start_year
  lu <- ((1 - (quant)) / 2)
  uu <- 1 - ((1 - (quant)) / 2)
  
  # get samples
  indt <- samples %>% 
    ungroup() %>% 
    filter(true_year %in% c(start_year, end_year)) %>% 
    select(-matches(match = year, ignore.case = FALSE)) %>% 
    pivot_wider(names_from = true_year,
                values_from = .value,
                names_prefix = "Y") %>% 
    rename_with(., ~gsub(replacement = "start",
                         pattern = paste0("Y",start_year), .x,
                         fixed = TRUE))%>% 
    rename_with(., ~gsub(replacement = "end",
                         pattern = paste0("Y",end_year), .x,
                         fixed = TRUE))%>% 
    rename_with(., ~gsub(replacement = "stratum_trend",
                         pattern = summary_regions, .x,
                         fixed = TRUE))
  
  # compute endpoint trends
  tt <- indt %>% 
    group_by(.draw, stratum_trend) %>% 
    summarise(end = sum(end),
              start = sum(start),
              t = texp(end / start,ny = nyrs),
              ch = chng(end / start),
              .groups = "keep") %>% 
    group_by(stratum_trend) %>% 
    summarise(trend = mean(t),
              lci = quantile(t, lu, names = FALSE),
              uci = quantile(t, uu, names = FALSE),
              q0.5 = quantile(t, 0.5),
              q0.05 = quantile(t, 0.05),
              q0.95 = quantile(t, 0.95),
              q0.165 = quantile(t, 0.165),
              q0.835 = quantile(t, 0.835),
              percent_change = median(ch),
              p_ch_lci = quantile(ch, lu, names = FALSE),
              p_ch_uci = quantile(ch, uu, names = FALSE),
              prob_decline = prob_dec(ch, 0),
              prob_decline_gt30 = prob_dec(ch, -30),
              prob_decline_gt50 = prob_dec(ch, -50),
              prob_decline_gt70 = prob_dec(ch, -70))%>% 
    rename_with(., ~gsub(replacement = summary_regions,
                         pattern = "stratum_trend", .x,
                         fixed = TRUE)) %>% 
    mutate(region_type = summary_regions,
           start_year = start_year,
           end_year = end_year)
  
  # save to list
  out_lst <- list(endpoint_trends = tt)
  
  # get smooth samples
  indt <- samples_smooth %>%
    ungroup() %>% 
    filter(true_year %in% c(start_year, end_year)) %>% 
    select(-matches(match = year, ignore.case = FALSE)) %>% 
    pivot_wider(names_from = true_year,
                values_from = .value,
                names_prefix = "Y") %>% 
    rename_with(., ~gsub(replacement = "start",
                         pattern = paste0("Y",start_year), .x,
                         fixed = TRUE))%>% 
    rename_with(., ~gsub(replacement = "end",
                         pattern = paste0("Y",end_year), .x,
                         fixed = TRUE))%>% 
    rename_with(., ~gsub(replacement = "stratum_trend",
                         pattern = summary_regions, .x,
                         fixed = TRUE))
  
  # compute smooth endpoint trends
  tt <- indt %>% 
    group_by(.draw, stratum_trend) %>% 
    summarise(end = sum(end),
              start = sum(start),
              t = texp(end / start, ny = nyrs),
              ch = chng(end / start),
              .groups = "keep") %>% 
    group_by(stratum_trend) %>% 
    summarise(trend = mean(t),
              lci = quantile(t, lu, names = FALSE),
              uci = quantile(t, uu, names = FALSE),
              q0.5 = quantile(t, 0.5),
              q0.05 = quantile(t, 0.05),
              q0.95 = quantile(t, 0.95),
              q0.165 = quantile(t, 0.165),
              q0.835 = quantile(t, 0.835),
              percent_change = median(ch),
              p_ch_lci = quantile(ch, lu, names = FALSE),
              p_ch_uci = quantile(ch, uu, names = FALSE),
              prob_decline = prob_dec(ch, 0),
              prob_decline_gt30 = prob_dec(ch, -30),
              prob_decline_gt50 = prob_dec(ch, -50),
              prob_decline_gt70 = prob_dec(ch, -70))%>% 
    rename_with(., ~gsub(replacement = summary_regions,
                         pattern = "stratum_trend", .x,
                         fixed = TRUE)) %>% 
    mutate(region_type = summary_regions,
           start_year = start_year,
           end_year = end_year)
  
  # save to list
  out_lst$smooth_endpoint_trends <- tt
  
  # return list
  return(out_lst)
}
# ------------------------------------------------------------------------------






# helper functions --------------------------------------------------------
# add smooth indices to index list
smooth_sample <- function(y, x, min_interval_yrs = 5){
  require(mgcv)
  ny <- length(x)
  kts <- floor(ny / min_interval_yrs)
  yhat <- exp(fitted(gam(log(y) ~ 1 + s(x, k = kts, bs="tp"),
                         family = "gaussian")))
  return(yhat)
}

# calculate endpoint trends
texp <- function(x,ny = 2019-1966){
  (x^(1 / ny) - 1) * 100
}

# calculate percent change
chng <- function(x){
  (x - 1) * 100
}

# summarize trends
prob_dec <- function(ch,thresh){
  length(which(ch < thresh)) / length(ch)
}

### function to extract the dimension values from an bayesian model fit
### works within the gather_samples function
# this function is used in posterior_samples()
dim_ext <- function(dim = 1,
                    var = "",
                    cl = "Parameter",
                    dat = NULL){
  # function to extract the indicator values from cmdstanr output
  require(stringr)
  
  pat <- paste0("(?<=", var, "\\[")
  
  if(dim > 1){
    for(j in 1:(dim - 1)){
      pat2 <- paste0(pat, ")[:digit:]+")
      cl2 <- str_extract(unlist(dat[, cl]), pattern = pat2)
      d <- max(nchar(cl2))
      pat <- paste0(pat, "[:digit:]{1,", d, "}[:punct:]")
    }
  }
  
  pat <- paste0(pat, ")[:digit:]+")
  dds <- as.integer(str_extract(unlist(dat[,cl]), pattern = pat))
  return(dds)
}

### function to generate the same tidy output as gather-draws in tidybayes package
## dims should be a character vector defining the dimensions of the parameter
## e.g., parm = "nsmooth", dims = c("stratum","year"),
## function works with cmdstanr output by default and rstan fits
## with is_rstan == TRUE
# this function is used in index_function()
posterior_samples <- function(cmdstan_draws = draws1,
                              parm = "n",
                              dims = NULL){
  # require(posterior) # not loaded to avoid masking
  require(tidyverse)
  
  if(length(dims) > 0){
    parm_ex <- paste0(parm,"[")
  } else {
    parm_ex <- parm
  }
  
  # get subset of samples
  samples <- posterior::as_draws_df(
    posterior::subset_draws(draws1, variable = parm))

  plong <- suppressWarnings(samples %>% pivot_longer(
    cols = starts_with(parm_ex, ignore.case = FALSE),
    names_to = c(".variable"),
    values_to = ".value",
    values_drop_na = TRUE
  )) 
  
  for(dn in 1:length(dims)){
    dd <- dims[dn]
    plong[,dd] <- dim_ext(dim = dn, var = parm, cl = ".variable", dat = plong)
  }
  
  plong <- plong %>% mutate(.variable = parm)
  return(plong)
}


# posterior summaries function
posterior_sums <- function(samples = n_samples,
                           quantiles = c(0.025, 0.05, 0.165, 0.5, 0.835, 
                                         0.95, 0.975),
                           ci = 0.95,
                           dims = NULL){
  
  cis = c((1 - ci) / 2, 1 - ((1 - ci) / 2))
  
  if(!is.null(dims)){
    for(i in 1:length(dims)){
      samples <- samples %>% 
        rename_with(~gsub(dims[i], paste0("d",i), .x, fixed = TRUE))
    }
    
    if(length(dims) == 1){
      sums <- samples %>% 
        group_by(d1) %>% 
        summarise(mean = mean(.value),
                  median = median(.value),
                  sd = sd(.value),
                  lci = quantile(.value, cis[1]),
                  uci = quantile(.value, cis[2]),
                  .groups = "keep") %>% 
        rename_with(~gsub("d1", dims[1], .x, fixed = TRUE))
      
      if(!is.null(quantiles)){
        qs <- paste0("Q_", gsub(quantiles, pattern = "0.", replacement = "", 
                               fix = TRUE))
        for(i in 1:length(quantiles)){
          qq <- quantiles[i]
          qn <- qs[i]
          sumt <- samples %>% 
            group_by(d1) %>% 
            summarise(tt = as.numeric(quantile(.value, qq)),
                      .groups = "keep") %>% 
            rename_with(~gsub(pattern = "tt", replacement = qn, .x, 
                              fixed = TRUE)) %>% 
            rename_with(~gsub("d1", dims[1], .x, fixed = TRUE)) 
          sums <- left_join(sums, sumt, by = dims)
        }
      }
    }
    
    if(length(dims) == 2){
      sums <- samples %>% 
        group_by(d1, d2) %>% 
        summarise(mean = mean(.value),
                  median = median(.value),
                  sd = sd(.value),
                  lci = quantile(.value, cis[1]),
                  uci = quantile(.value, cis[2]),
                  .groups = "keep") %>% 
        rename_with(~gsub("d1", dims[1], .x, fixed = TRUE)) %>% 
        rename_with(~gsub("d2", dims[2], .x, fixed = TRUE))
      
      if(!is.null(quantiles)){
        qs <- paste0("Q_", gsub(quantiles, pattern = "0.", replacement = "",
                               fixed = TRUE))
        for(i in 1:length(quantiles)){
          qq <- quantiles[i]
          qn <- qs[i]
          sumt <- samples %>% 
            group_by(d1, d2) %>% 
            summarise(tt = as.numeric(quantile(.value, qq)),
                      .groups = "keep") %>% 
            rename_with(~gsub(pattern = "tt", replacement = qn, .x, 
                              fixed = TRUE)) %>% 
            rename_with(~gsub("d1", dims[1], .x, fixed = TRUE))  %>% 
            rename_with(~gsub("d2", dims[2], .x, fixed = TRUE))
          sums <- left_join(sums, sumt, by = dims)
        }
      }
    }
    
    if(length(dims) == 3){
      sums <- samples %>% 
        group_by(d1, d2, d3) %>% 
        summarise(mean = mean(.value),
                  median = median(.value),
                  sd = sd(.value),
                  lci = quantile(.value, cis[1]),
                  uci = quantile(.value, cis[2]),
                  .groups = "keep") %>% 
        rename_with(~gsub("d1", dims[1], .x, fixed = TRUE)) %>% 
        rename_with(~gsub("d2", dims[2], .x, fixed = TRUE)) %>% 
        rename_with(~gsub("d3", dims[3], .x, fixed = TRUE))
      
      if(!is.null(quantiles)){
        qs <- paste0("Q_", gsub(quantiles, pattern = "0.", replacement = "",
                               fix = TRUE))
        for(i in 1:length(quantiles)){
          qq <- quantiles[i]
          qn <- qs[i]
          sumt <- samples %>% 
            group_by(d1, d2, d3) %>% 
            summarise(tt = as.numeric(quantile(.value, qq)),
                      .groups = "keep") %>% 
            rename_with(~gsub(pattern = "tt",replacement = qn, .x, 
                              fixed = TRUE)) %>% 
            rename_with(~gsub("d1", dims[1], .x, fixed = TRUE))  %>% 
            rename_with(~gsub("d2", dims[2], .x, fixed = TRUE))  %>% 
            rename_with(~gsub("d3", dims[3], .x, fixed = TRUE))
          sums <- left_join(sums, sumt, by = dims)
        }
      }
    }
    
    if(length(dims) == 4){
      sums <- samples %>% 
        group_by(d1, d2, d3, d4) %>% 
        summarise(mean = mean(.value),
                  median = median(.value),
                  sd = sd(.value),
                  lci = quantile(.value, cis[1]),
                  uci = quantile(.value, cis[2]),
                  .groups = "keep") %>% 
        rename_with(~gsub("d1", dims[1], .x, fixed = TRUE)) %>% 
        rename_with(~gsub("d2", dims[2], .x, fixed = TRUE)) %>% 
        rename_with(~gsub("d3", dims[3], .x, fixed = TRUE)) %>% 
        rename_with(~gsub("d4", dims[4], .x, fixed = TRUE))
  
      if(!is.null(quantiles)){
        qs <- paste0("Q_", gsub(quantiles, pattern = "0.", replacement = "",
                               fixed = TRUE))
        for(i in 1:length(quantiles)){
          qq <- quantiles[i]
          qn <- qs[i]
          sumt <- samples %>% 
            group_by(d1, d2, d3, d4) %>% 
            summarise(tt = as.numeric(quantile(.value, qq)),
                      .groups = "keep") %>% 
            rename_with(~gsub(pattern = "tt", replacement = qn, .x,
                              fixed = TRUE)) %>% 
            rename_with(~gsub("d1", dims[1], .x, fixed = TRUE))  %>% 
            rename_with(~gsub("d2", dims[2], .x, fixed = TRUE))  %>% 
            rename_with(~gsub("d3", dims[3], .x, fixed = TRUE))  %>% 
            rename_with(~gsub("d4", dims[4], .x, fixed = TRUE))
          sums <- left_join(sums, sumt, by = dims)
        }
      }
    }
    
    if(length(dims) > 4){stop("ERROR this function cannot handle more than 4 dimensions")}
    
  }else{
    
    sums <- samples %>% 
      summarise(mean = mean(.value),
                median = median(.value),
                sd = sd(.value),
                lci = quantile(.value, cis[1]),
                uci = quantile(.value, cis[2]),
                .groups = "keep")
    
    if(!is.null(quantiles)){
      qs <- paste0("Q_", gsub(quantiles, pattern = "0.", replacement = "",
                             fix = TRUE))
      for(i in 1:length(quantiles)){
        qq <- quantiles[i]
        qn <- qs[i]
        sumt <- samples %>% 
          summarise(tt = as.numeric(quantile(.value, qq)),
                    .groups = "keep") %>% 
          rename_with(~gsub(pattern = "tt", replacement = qn, .x, fixed = TRUE))
        sums <- bind_cols(sums, sumt)
      }
    }
  }
  
  return(sums)
  
}
# ------------------------------------------------------------------------------

