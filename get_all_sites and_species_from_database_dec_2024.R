

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




# make sample file from cbc database ###########################################
# generic query stuff
first_count_num = 67
last_count_num = 124
query_text = "SELECT loc_circle.abbrev,
  loc_circle.subnational_code, loc_circle.name, loc_circle.Latitude,
  loc_circle.Longitude, cnt_submission.count_yr, ref_transportation.description, 
  cnt_effort_time_distance.distance, ref_dist_unit.description, 
  cnt_effort_time_distance.hours FROM 
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
  # collapse by mode
  group_by(state_province=state, circle_code=abbrev, circle_name=name, 
           count_number=count_yr) %>% 
  summarise(longitude=mean(Longitude, na.rm=T), latitude=mean(Latitude, na.rm=T),
            field_party_hours=sum(hours, na.rm=T)) %>%
  ungroup() %>% 
  # clean up 
  mutate(count_year=count_number+1899) %>% 
  select(1,2,3,4,8,5,6,7) %>% 
  mutate(field_party_hours=ifelse(field_party_hours<=0, NA, 
                                  field_party_hours)) %>% 
  # remove no effort data
  filter(!is.na(field_party_hours)) %>% 
  mutate(join_code=paste(circle_code, count_year, sep="-")) %>% 
  arrange(circle_code, count_number)

# save data
getwd()
write.csv(dat2, "cbc_sample_data_Dec_1966_to_Jan_2024_cont_USA_CAN.csv", 
          na="", row.names=F)
#-------------------------------------------------------------------------------




# make count data file from  cbc database ######################################
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

write.csv(dat3, "cbc_count_data_Dec_1966_to_Jan_2024_cont_USA_CAN.csv", 
          na="", row.names=F)
#-------------------------------------------------------------------------------




# make taxon map ---------------------------------------------------------------
# counts
cnts1 <- read_csv("cbc_count_data_Dec_1966_to_Jan_2024_cont_USA_CAN.csv")
spp_tab1 <- cnts1 %>% arrange(common_name, number_counted) %>% 
  filter(!grepl("sp.", common_name, fixed=T)) %>% 
  filter(!grepl("hybrid", common_name, fixed=T)) %>% 
  filter(!grepl("NA", common_name, fixed=T)) %>% 
  filter(!is.na(common_name)) %>% 
  group_by(common_name) %>% summarise(total_counted=sum(number_counted))
ebird_tax <- 
  read_csv("eBird-Clements-v2024-integrated-checklist-October-2024-rev.csv") %>% 
  rename_with(., ~tolower(gsub(" ", "_", .x, fixed = TRUE))) %>% 
  select(species_code, common_name=english_name, scientific_name) %>% 
  arrange(common_name)

spp_tab1 <- spp_tab1 %>% left_join(ebird_tax)

write_csv(spp_tab1, "temporary_species_table.csv")
# ------------------------------------------------------------------------------




# modify manually in excel -----------------------------------------------------
# ------------------------------------------------------------------------------




# bring files back in and merge ------------------------------------------------
dir()
tst1 <- read.csv("temporary_species_table.csv")
ost1 <- read.csv("cbc_trends_spp_list_66_to_121.csv")

# clean query columns
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

# join and export again for manual tweaking
jst1 <- tst2 %>% full_join(ost1, by = join_by(ebird_spp_code))
write_csv(jst1, "temporary_species_table_v2.csv")
# ------------------------------------------------------------------------------




# import bcr and sample files --------------------------------------------------
# bcrs
library(sf)

map1 <- bbsBayes2::load_map("bbs_cws")
plot(map1$geom, col="red")

# samples
samps0 <- read_csv("cbc_sample_data_Dec_1966_to_Jan_2024_cont_USA_CAN.csv") 

samps1 <- samps0 %>% 
  group_by(state_province, circle_code, circle_name, count_number,
           count_year,longitude,latitude,field_party_hours,join_code)










# ------------------------------------------------------------------------------










