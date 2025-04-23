

# get and prepare count data function
get_and_filter_count_data <- function(){

  # libraries
  require(tidyverse)
  
  # filter count table based on taxon key
  potential_common_names <- as.character(unlist(str_split(historic_cbc_com_name_s, ",")))
  # Encoding(potential_common_names) <- "UTF-8"
  # Encoding(count_table$common_name) <- "UTF-8"
  
  
 # test <- c("Woodhouse's Scrub-Jay","California Scrub-Jay","Florida Scrub-Jay",
 #           "Island Scrub-Jay","Western Scrub-Jay","Western Scrub-Jay (Coastal)",
 #           "Western Scrub-Jay (Woodhouse's)")
 # potential_common_names <- stri_enc_toutf8(potential_common_names)
 # 
 # potential_common_names <- stri_encode(potential_common_names, "ASCII", "UTF-8")
 # 
 # table(Encoding(count_table$common_name))
 # table(Encoding(species_table$ebird_com_name))
 # library(stringi)
 # stri_enc_mark(count_table$common_name)
 # all(stri_enc_isutf8(count_table$common_name))
 # stri_enc_mark(species_table$ebird_com_name)
 # all(stri_enc_isutf8(species_table$ebird_com_name))
 # count_table$common_name <- stri_enc_toutf8(count_table$common_name)
 # species_table$ebird_com_name <- stri_enc_toutf8(species_table$ebird_com_name)
  
  count1 <- count_table %>% filter(common_name %in% potential_common_names)
  # count1 <- count_table %>% filter(common_name %in% test)
  
  # join count and site tables for zero filling
  dat1 <- site_table %>% left_join(count1, by="join_code") %>% 
    mutate(common_name=ebird_com_name_s,
           ebird_species_code=ebird_spp_code_s,
           ebird_common_name=ebird_com_name_s,
           ebird_sci_name=ebird_sci_name_s,
           number_counted=ifelse(is.na(number_counted), 0, number_counted))
  
  # filter for current species geography
  lon_filter <- as.numeric(unlist(str_split(lon_filter_s, pattern=",")))
  lat_filter <- as.numeric(c(str_sub(lat_filter_s, 1, 2), 
                  str_sub(lat_filter_s, 4, nchar(lat_filter_s))))
  prov_state_filter <- unlist(str_split(prov_state_filter_s, ","))
  bcr_filter <- str_replace(unlist(str_split(bcr_filter_s, ",")), "^0+", "")
  dat2 <- dat1 %>% filter(abs(longitude)*-1 <= lon_filter[1] & abs(longitude)*-1 >= lon_filter[2])
  dat2 <- dat2 %>% filter(latitude <= lat_filter[1] & latitude >= lat_filter[2])
  if(prov_state_filter[1]=="!ZZ") {
    dat2 <- dat2
  }
  if(prov_state_filter[1]!="!ZZ" & !grepl("!", prov_state_filter[1])) {
    dat2 <- dat2 %>% filter(prov_state %in% prov_state_filter)
  }
  if(prov_state_filter[1]!="!ZZ" & grepl("!", prov_state_filter[1])) {
    dat2 <- dat2 %>% filter(!prov_state %in% gsub("!", "", prov_state_filter))
  }
  if(bcr_filter[1]=="!ZZ") {
    dat2 <- dat2
  }
  if(bcr_filter[1]!="!ZZ" & !grepl("!", bcr_filter[1])) {
    dat2 <- dat2 %>% filter(bcr %in% bcr_filter)
  }
  if(bcr_filter[1]!="!ZZ" & grepl("!", bcr_filter[1])) {
    dat2 <- dat2 %>% filter(!bcr %in% gsub("!", "", bcr_filter))
  }
  
  # create hours column for modeling
  add_nocturnal <- as.numeric(add_nocturnal_s)
  add_feeder <- as.numeric(add_feeder_s)
  dat2 <- dat2 %>% mutate(tot_hours=field_hours) 
  if(add_nocturnal==1) dat2$tot_hours <- dat2$field_hours + dat2$nocturnal_hours
  if(add_feeder==1) dat2$tot_hours <- dat2$field_hours + dat2$feeder_hours
  if(add_nocturnal==1 & add_feeder==1) dat2$tot_hours <- dat2$field_hours + 
    dat2$nocturnal_hours + dat2$feeder_hours
  
  # select strata for inclusion and calculate z proportion
  strat_select_df <- dat2 %>% 
    select(strata_name, circle_code, count_year, number_counted) %>% 
    mutate(detected=ifelse(number_counted>0, 1, 0)) %>% 
    select(-number_counted) %>% 
    group_by(strata_name, circle_code) %>% 
    summarise(number_years_detected=sum(detected)) %>% 
    ungroup() %>% 
    group_by(strata_name) %>% 
    summarise(circles_per_stratum=length(unique(circle_code)),
           nonzero_circles=sum(number_years_detected >=
                                 number_years_per_circle_threshold)) %>% 
    ungroup() %>% 
    filter(nonzero_circles>=nonzero_circles_threshold) %>% 
    mutate(nonzero_ratio_z=nonzero_circles/circles_per_stratum)
  selected_strata <- strat_select_df %>% pull(strata_name)
  
  # filter count data for selected strata
  dat3 <- dat2 %>% filter(strata_name %in% selected_strata) %>% 
    left_join(strat_select_df, by="strata_name") %>% 
    mutate(ebird_species_code=ebird_spp_code_s,
           ebird_common_name=ebird_com_name_s,
           ebird_sci_name=ebird_sci_name_s,
           scaled_effort=tot_hours/mean(tot_hours),
           year_vec=as.numeric(factor(count_year)),
           strata_vec=as.numeric(factor(strata_name)),
           circle_vec=as.numeric(factor(circle_code))) %>% 
    select(strata_name, area_sq_km, strata_vec, 
           country, country_code, prov_state, bcr,
      circles_per_stratum, nonzero_circles, nonzero_ratio_z, 
      circle_code, circle_name, longitude, latitude, circle_vec,
      count_number, count_year, year_vec,
      field_hours, feeder_hours, nocturnal_hours, tot_hours,scaled_effort,
      ebird_common_name, ebird_sci_name, ebird_species_code, number_counted)
  
  # save data for modeling
  fname1 <- file.path(results_dir, gsub(" ", "_", ebird_com_name_s),
                      paste0(gsub(" ", "_", ebird_com_name_s), 
                             "_zero_filled_and_filtered_modeling_data.csv"))
  write.csv(dat3, fname1, na="", row.names=F)
  
} # end function


# test function
# get_and_filter_count_data()
