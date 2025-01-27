# setup ------------------------------------------------------------------------
library(tidyverse)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
getwd()
# ------------------------------------------------------------------------------


# select species current common name
species <- "Sharp-shinned Hawk"

# grab three data files
taxon0 <- read_csv("./data/taxon_key_dec_2024.csv")
count0 <- read_csv("./output/count_table_dec_2024.csv")
site0 <- read_csv("./data/site_table_dec_2024.csv")

# filter count table based on taxon key
target_row <- taxon0 %>% filter(ebird_com_name==species)
potential_common_names <- unlist(str_split(target_row$historic_cbc_com_name, 
                                           pattern=","))
count1 <- count0 %>% filter(common_name %in% potential_common_names)

# join count and site tables
dat1 <- site0 %>% left_join(count1, by="join_code") %>% 
  mutate(common_name=species,
         ebird_species_code=target_row$ebird_spp_code,
         ebird_common_name=target_row$ebird_com_name,
         ebird_sci_name=target_row$ebird_sci_name,
         number_counted=ifelse(is.na(number_counted), 0, number_counted))

# filter for current species geography
head(as.data.frame(target_row))
lon_filter <- as.numeric(unlist(str_split(target_row$lon_filter, pattern=",")))
lat_filter <- as.numeric(c(str_sub(target_row$lat_filter, 1, 2), 
                str_sub(target_row$lat_filter, 3, nchar(target_row$lat_filter))))
prov_state_filter <- target_row$prov_state_filter
bcr_filter <- target_row$bcr_filter






# stop here. go back and figure out how to assign bcr to islands!
dat1 <- dat1 %>% filter(longitude <= lon_filter[1] & longitude >= lon_filter[2])
dat1 <- dat1 %>% filter(latitude <= lat_filter[1] & latitude >= lat_filter[2])
if(prov_state_filter!="!ZZ") dat1 <- dat1 %>% filter(prov_state %in% prov_state_filter)
if(bcr_filter!="!ZZ") dat1 <- dat1 %>% filter(bcr_state %in% bcr_filter)







# then work on this part
# create hours column for modeling
add_nocturnal <- as.numeric(target_row$add_nocturnal)
add_feeder <- as.numeric(target_row$add_feeder)





# select strata for inclusion and calculate z proportion
number_years_per_circle_threshold <- 5 # minimum
nonzero_circles_threshold <- 3 # minimum
strat_select_df <- dat1 %>% 
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
dat2 <- dat1 %>% filter(strata_name %in% selected_strata) %>% 
  left_join(strat_select_df, by="strata_name") %>% 
  mutate(ebird_species_code=target_row$ebird_spp_code,
         ebird_common_name=target_row$ebird_com_name,
         ebird_sci_name=target_row$ebird_sci_name) %>% 
  select(strata_name, area_sq_km, country, country_code, prov_state, bcr,
    circles_per_stratum, nonzero_circles, nonzero_ratio_z, 
    circle_code, circle_name, longitude, latitude, 
    count_number, count_year,
    field_hours, feeder_hours, nocturnal_hours, 
    ebird_common_name, ebird_sci_name, ebird_species_code, number_counted)






# save data for modeling
fname1 <- paste0("./output/", gsub(" ", "_", species), 
                 "_zero_filled_and_filtered_modeling_data.csv")
write.csv(dat2, fname1, na="", row.names=F)




