# this script has chunks to make three distinct types of files. first is the 
# sample or circle file that has all the circles along with information on
# location, strata id, and three types of effort.
#
# the second is a count file that has species, count and a join field that joins
# with the sample or circle file. it has all count records of all species, and 
# taxonomy is historic not present taxonomy.
#
# the third file is a taxonomic key.it maps the present taxonmic id to historic
# ones in the cbc database.
#
# eventually, the taxon key is used to query the count file for all counts for
# a species. then the species count file is joined with the sample file, and the
# sample file is filtered based on current species geography. then stratum 
# selection and Z proportion is calculated per stratum and some strata are 
# removed before modeling.
#
# tim meehan, january 2025



# setup ------------------------------------------------------------------------
# install the required packages
print("Set up query_cbc_data function.")
library(RODBC)
library(sf)
library(bbsBayes2)
library(tidyverse)
options(scipen = 9999999)
options(stringsAsFactors=F)

#set working directory
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
getwd()
# ------------------------------------------------------------------------------




# make sample file from cbc database -------------------------------------------
# generic query stuff
first_count_num = 67
last_count_num = 124
query_text = "SELECT loc_circle.abbrev,
  loc_circle.subnational_code, loc_circle.name, loc_circle.Latitude,
  loc_circle.Longitude, cnt_submission.count_yr, ref_transportation.description, 
  cnt_effort_time_distance.distance, ref_dist_unit.description, 
  cnt_effort_time_distance.hours,
  cnt_effort_feeder_night.Feeder_hrs,                 
  cnt_effort_feeder_night.Nocturnal_hrs FROM 
  cnt_submission FULL JOIN
  loc_circle ON loc_circle.circle_id = cnt_submission.circle_id FULL JOIN
  cnt_effort_time_distance ON cnt_submission.submission_id =
  cnt_effort_time_distance.submission_id FULL JOIN 
  cnt_effort_feeder_night ON 
  cnt_submission.submission_id = cnt_effort_feeder_night.submission_id FULL JOIN 
  ref_transportation ON cnt_effort_time_distance.trans_id =
  ref_transportation.trans_id FULL JOIN 
  ref_dist_unit ON 
  cnt_effort_time_distance.distance_unit_id = ref_dist_unit.unit_id WHERE"

# establish a connection
channel <- odbcConnect(dsn = "CBC database", uid = "azurecbc",
                       pwd = "bcVirtualM@chine")

sqlTables(channel, tableType = "TABLE")
res <- sqlFetch(channel, "CNT_EFFORT_FEEDER_NIGHT", max=5)

# sql query database
dat1 <- sqlQuery(channel, paste(query_text, 
                                "CNT_SUBMISSION.COUNT_YR BETWEEN",
                                first_count_num, "AND", last_count_num))

# close the connection
odbcClose(channel)

# clean data
dat2 <- dat1 %>% 
  # keep us and ca
  filter(grepl("CA-|US-", subnational_code)) %>% 
  # remove hi
  filter(subnational_code!="US-HI") %>%
  # recode dc
  rename(prov_state=subnational_code) %>% 
  mutate(prov_state=str_trim(str_split(prov_state, "-", simplify=T)[,2])) %>% 
  mutate(prov_state=ifelse(prov_state=="DC", "MD", prov_state)) %>% 
  # deal with dirty values
  mutate(hours=ifelse(hours<=0, NA, hours),
         Feeder_hrs=ifelse(Feeder_hrs<=0, NA, Feeder_hrs),
         Nocturnal_hrs=ifelse(Nocturnal_hrs<=0, NA, Nocturnal_hrs)) %>% 
  # collapse by mode
  group_by(prov_state, circle_code=abbrev, circle_name=name, 
           count_number=count_yr) %>% 
  summarise(longitude=mean(Longitude, na.rm=T), latitude=mean(Latitude, na.rm=T),
            field_hours=sum(hours, na.rm=T),
            feeder_hours=mean(Feeder_hrs, na.rm=T),
            nocturnal_hours=mean(Nocturnal_hrs)) %>%
  ungroup() %>% 
  # zero out nocturnal and feeder
  mutate(feeder_hours=ifelse(is.na(feeder_hours), 0, feeder_hours),
         nocturnal_hours=ifelse(is.na(nocturnal_hours), 0, nocturnal_hours)) %>%
  # clean up 
  mutate(count_year=count_number+1899) %>% 
  select(prov_state, circle_code, circle_name, count_number, count_year, 
  longitude, latitude, field_hours, feeder_hours, nocturnal_hours) %>% 
  mutate(field_hours=ifelse(field_hours<=0, NA, field_hours)) %>% 
  # remove no effort data
  filter(!is.na(field_hours)) %>% 
  mutate(join_code=paste(circle_code, count_year, sep="-")) %>% 
  arrange(circle_code, count_number)

# save data
getwd()
write.csv(dat2, "cbc_sample_data_Dec_1966_to_Jan_2024_cont_USA_CAN.csv", 
          na="", row.names=F)

# assign analytical strata to sample file 
map1 <- bbsBayes2::load_map("bbs_cws")
plot(map1$geom, col="red")

# samples
samps0 <- read_csv("cbc_sample_data_Dec_1966_to_Jan_2024_cont_USA_CAN.csv") 
samps1 <- samps0 %>% st_as_sf(coords = c("longitude", "latitude"), 
                              remove = F, crs=4326) %>% 
  st_transform(crs=st_crs(map1))

# spatial join
samps2 <- samps1 %>% select(-prov_state) %>% st_join(map1) %>% 
  select(-bcr_by_country) %>% st_drop_geometry()

# save as circle file
write.csv(samps2, "temporary_circle_table_v2.csv", na="", row.names=F)
# ------------------------------------------------------------------------------





# make count data file from cbc database --------------------------------------
# generic query stuff
query_text = "SELECT loc_circle.abbrev,
  loc_circle.subnational_code, loc_circle.name, loc_circle.Latitude,
  loc_circle.Longitude, cnt_submission.count_yr,
  ref_species.com_name, cnt_observation.how_many, ref_transportation.description,
  cnt_effort_time_distance.distance, ref_dist_unit.description,
  cnt_effort_time_distance.hours FROM cnt_submission FULL JOIN
  loc_circle ON loc_circle.circle_id = cnt_submission.circle_id FULL JOIN
  cnt_observation ON cnt_submission.submission_id =
  cnt_observation.submission_id FULL JOIN ref_species ON
  cnt_observation.species_id = ref_species.species_id FULL JOIN
  cnt_effort_time_distance ON cnt_submission.submission_id =
  cnt_effort_time_distance.submission_id FULL JOIN cnt_effort_feeder_night ON
  cnt_submission.submission_id = cnt_effort_feeder_night.submission_id
  FULL JOIN ref_transportation ON cnt_effort_time_distance.trans_id =
  ref_transportation.trans_id FULL JOIN ref_dist_unit  ON
  cnt_effort_time_distance.distance_unit_id = ref_dist_unit.unit_id WHERE"

# establish a connection
channel <- odbcConnect(dsn = "CBC database", uid = "tmeehan",
                       pwd = "Aud2MPitWEB@R5")

# sql query database
dat1 <- sqlQuery(channel, paste(query_text, 
                                "CNT_SUBMISSION.COUNT_YR BETWEEN",
                                first_count_num, "AND", last_count_num))

# close the connection
odbcClose(channel)

# clean data
dat2 <- dat1 %>% 
  # keep us and ca
  filter(grepl("CA-|US-", subnational_code)) %>% 
  # remove hi
  filter(subnational_code!="US-HI") %>%
  # recode dc
  rename(state=subnational_code) %>% 
  mutate(state=str_trim(str_split(state, "-", simplify=T)[,2])) %>% 
  mutate(state=ifelse(state=="DC", "MD", state)) %>% 
  # sum hours
  filter(!is.na(how_many)) %>% 
  filter(how_many > 0) %>% 
  # collapse by mode
  group_by(state, circle_code=abbrev, circle_name=name, 
           count_year=count_yr, com_name) %>% 
  summarise(lon=mean(Longitude), lat=mean(Latitude),
            count=max(how_many, na.rm=T), hours=sum(hours, na.rm=T)) %>% 
  ungroup() 

# format for joining
names(dat2)
dat3 <- dat2 %>% 
  filter(count>0, hours>0) %>% 
  select(state_province=state, circle_code, circle_name, 
         count_number=count_year, common_name=com_name, number_counted=count) %>% 
  mutate(count_year=count_number+1899, 
         join_code=paste(circle_code, count_year, sep="-")) %>% 
  select(5,6,8) %>% 
  arrange(join_code, common_name)

# save
write.csv(dat3, "cbc_count_data_Dec_1966_to_Jan_2024_cont_USA_CAN.csv", 
          na="", row.names=F)

# final count data manipulation
dat4 <- read.csv("cbc_count_data_Dec_1966_to_Jan_2024_cont_USA_CAN.csv")
write.csv(dat4, "temporary_count_table_v2.csv", na="", row.names=F)
# ------------------------------------------------------------------------------





# make taxon map ---------------------------------------------------------------
# counts
cnts1 <- read_csv("cbc_count_data_Dec_1966_to_Jan_2024_cont_USA_CAN.csv")
spp_tab1 <- cnts1 %>% arrange(common_name, number_counted) %>% 
  filter(!grepl("sp.", common_name, fixed=T)) %>% 
  filter(!grepl("hybrid", common_name, fixed=T)) %>% 
  filter(!grepl("NA", common_name, fixed=T)) %>% 
  filter(!is.na(common_name)) %>% 
  group_by(common_name) %>% summarise(total_counted=sum(number_counted))

# current ebird taxonomy
ebird_tax <- 
  read_csv("eBird-Clements-v2024-integrated-checklist-October-2024-rev.csv") %>% 
  rename_with(., ~tolower(gsub(" ", "_", .x, fixed = TRUE))) %>% 
  select(species_code, common_name=english_name, scientific_name) %>% 
  arrange(common_name)

# initial taxon join
spp_tab1 <- spp_tab1 %>% left_join(ebird_tax)

# write out for manual tweaking
write_csv(spp_tab1, "temporary_species_table.csv")
# ------------------------------------------------------------------------------





# bring files back in and merge ------------------------------------------------
dir()
tst1 <- read.csv("temporary_species_table.csv")
ost1 <- read.csv("cbc_trends_spp_list_66_to_121.csv")

# clean query columns
# in the next round, no need to do this as we won't query database for each 
# species
head(tst1)
tst2 <- tst1 %>% mutate_if(grepl("q", names(tst1)), 
                   ~paste0("ref_species.com_name = ", "'", ., "'")) %>% 
  mutate(q2=ifelse(q2=="ref_species.com_name = ''", NA, q2),
         q3=ifelse(q3=="ref_species.com_name = ''", NA, q3),
         q4=ifelse(q4=="ref_species.com_name = ''", NA, q4),
         q5=ifelse(q5=="ref_species.com_name = ''", NA, q5),
         q6=ifelse(q6=="ref_species.com_name = ''", NA, q6),
         q7=ifelse(q7=="ref_species.com_name = ''", NA, q7),
         q8=ifelse(q8=="ref_species.com_name = ''", NA, q8),
         q9=ifelse(q9=="ref_species.com_name = ''", NA, q9),
         q10=ifelse(q10=="ref_species.com_name = ''", NA, q10),
         q11=ifelse(q11=="ref_species.com_name = ''", NA, q11),
         q12=ifelse(q12=="ref_species.com_name = ''", NA, q12),
         q13=ifelse(q13=="ref_species.com_name = ''", NA, q13)) %>% 
  tidyr::unite(col="q0", names(tst1)[grep("q", names(tst1))], sep=" OR ", 
               remove=T, na.rm=T) %>% 
  mutate(q0=str_trim(q0))

# join and export again for more manual tweaking
jst1 <- tst2 %>% full_join(ost1, by = join_by(ebird_spp_code))
write_csv(jst1, "temporary_species_table_v2.csv")

# 2024 only, remove ref species stuff
dat1 <- read_csv("./data/temporary_species_table_v2.csv")
dat2 <- dat1 %>% 
  rename(historic_cbc_com_name=query_text,
         prov_state_filter=state_filter) %>% 
  mutate(historic_cbc_com_name=gsub("ref_species.com_name = ", "", 
                                    historic_cbc_com_name)) %>% 
  mutate(historic_cbc_com_name=gsub("\'", "", historic_cbc_com_name),
    lon_filter=gsub("\"", "", lon_filter), 
    lat_filter=gsub("\"", "", lat_filter), 
    prov_state_filter=gsub("\"", "", prov_state_filter), 
    bcr_filter=gsub("\"", "", bcr_filter)) %>% 
  mutate(historic_cbc_com_name=gsub(" OR ", ",", historic_cbc_com_name))
write.csv(dat2, "temporary_species_table_v2.csv")
# ------------------------------------------------------------------------------







