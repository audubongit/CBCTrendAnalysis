

# setup ------------------------------------------------------------------------

# install the required packages
print("Set up query_cbc_data function.")
suppressPackageStartupMessages(library(RODBC))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
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
channel <- odbcConnect(dsn = "CBC database", uid = "tmeehan",
                       pwd = "Aud2MPitWEB@R5")

# sql query database
dat1 <- sqlQuery(channel, paste(query_text, 
                                "CNT_SUBMISSION.COUNT_YR BETWEEN",
                                first_count_num, "AND", last_count_num))

# close the connection
odbcClose(channel)

# clean data
library(stringr)
library(sf)
library(dplyr)
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
write.csv(dat2, "./output/cbc_sample_data_Dec_1966_to_Jan_2024_cont_USA_CAN.csv", 
          na="", row.names=F)
#------------------------------------------------------------------------------




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

write.csv(dat3, "./output/cbc_count_data_Dec_1966_to_Jan_2024_cont_USA_CAN.csv", 
          na="", row.names=F)
#-------------------------------------------------------------------------------





# extract subset of data #######################################################
# wd
setwd("D:/Users/tmeehan/Box/-.Science Team shared/2_Projects/6_CommunityScience/CBC/data_requests_filled_by_tim")
dir()

# clean data
library(sf)
library(dplyr)

# open and join data
samp_dat <- read.csv("cbc_sample_data_Dec_1930_to_Jan_2023_cont_USA_CAN.csv")
count_dat <- read.csv("cbc_count_data_Dec_1930_to_Jan_2023_cont_USA_CAN.csv")

# select species
taxa <- unique(count_dat$common_name)
focal_taxa <- grep("Grebe", taxa, fixed=T, value = T)
focal_taxa <- c(grep("Western", focal_taxa, fixed=T, value = T),
                grep("Clark's", focal_taxa, fixed=T, value = T))

# subset
focal_dat <- count_dat %>% filter(common_name %in% focal_taxa)

# join
all_count_dat <- focal_dat %>% left_join(samp_dat)
all_samp_dat <- samp_dat %>% left_join(focal_dat) %>% 
  mutate(common_name="Western/Clark's Grebe", 
         number_counted=ifelse(is.na(number_counted), 0, number_counted)) %>% 
  group_by(pick(-number_counted)) %>% 
  summarise(number_counted=sum(number_counted))

# save data
write.csv(all_samp_dat, "cbc_western-clark-grebe_data_Dec_1930_to_Jan_2023_cont_USA_CAN.csv", na="", row.names=F)
#-------------------------------------------------------------------------------

