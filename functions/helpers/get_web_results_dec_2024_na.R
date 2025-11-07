
# setup ------------------------------------------------------------------------
library(tidyverse)
output_dir <- "D:/Users/tim.meehan/Box/-.Tim.Meehan individual/CBCForWeb"
trends_dir <- "Z:/7_CommunityScience/CBCAnalysisResults/cbc_results_v2023.0_na"
code_dir <- "D:/Users/tim.meehan/Documents/GitHub/CBCTrendAnalysis/functions/helpers"
# ------------------------------------------------------------------------------



# for trend viewer -------------------------------------------------------------
setwd(code_dir)
spp_tab <- read.csv("../../data/taxon_key_dec_2024.csv") %>%
  dplyr::select(2,3,1) 

# get dir names
setwd(trends_dir)
dir_lst <- list.dirs(recursive = F)

# loop through dirs to get web app format data
out_all <- c()
i <- 1
for(i in 1:length(dir_lst)){
  spp_i <- gsub("./", "", as.character(dir_lst[i]))
  spp_i <- gsub("_", " ", as.character(spp_i))
  print(paste("trying", spp_i))
  tryCatch({
  setwd(trends_dir)
  setwd(dir_lst[i])
  fi <- list.files(pattern="relative_abundance_trends")
  f2 <- read.csv(fi) 
  # get trend estimates
  f3 <- f2 %>% filter(trend_type %in% c("long-term", "ten-year")) %>%
    select(1:6, 11:14)  
  # get abundance indices
  fi <- list.files(pattern="relative_abundance_indices")
  f4 <- read.csv(fi) 
  f5 <- f4 %>% filter(region_type=="continent", 
                      index_type=="first_diff") %>%
    select(year_start, index_median) %>%
    pivot_wider(names_from=year_start, values_from = index_median,
                names_prefix="F") 
  # get global trends
  f6 <- f3 %>% filter(trend_type=="long-term", region_type=="continent") %>%
    select(7:9) 
  names(f6) <- c("region_estimate_median", "region_estimate_lcl", 
                 "region_estimate_ucl")
  # put it together and stack
  f6 <- cbind(f6, f5) %>% slice(rep(1:n(), each = nrow(f3)))
  f7 <- cbind(f3, f6) 
  out_all <- rbind(out_all, f7)
  
  # back out and rerun
  setwd("../")
  print(paste("finished", spp_i))
}, error=function(e){})
}

# save
write.csv(out_all, file.path(output_dir, "na_viewer/cbc_trends_version_2023.0_na.csv"),
          na="", row.names=F)
# ------------------------------------------------------------------------------



# for download -----------------------------------------------------------------
out_all_trends <- c()
out_all_indices <- c()
i <- 1
for(i in 1:length(dir_lst)){
  spp_i <- gsub("./", "", as.character(dir_lst[i]))
  spp_i <- gsub("_", " ", as.character(spp_i))
  print(paste("trying", spp_i))
  tryCatch({
  setwd(trends_dir)
  setwd(dir_lst[i])
  # get trend estimates
  fi <- list.files(pattern="relative_abundance_trends")
  f2 <- read.csv(fi) %>% 
    filter(trend_type %in% c("long-term", "medium-term", "ten-year")) %>% 
    select(-trend_q0.165, -trend_q0.835, -length_3_gens,
           -trend_q0.05, -trend_q0.95, -trend_mean) %>% 
    rename(percent_change_median=percent_change_mean)
  # take out indices
  fi <- list.files(pattern="relative_abundance_indices")
  f3 <- read.csv(fi) %>% filter(index_type=="first_diff") %>% 
    select(-index_q0.165, -index_q0.835,
           -index_q0.05, -index_q0.95, -index_mean)
    # put it together and stack
  out_all_trends <- rbind(out_all_trends, f2)
  out_all_indices <- rbind(out_all_indices, f3)
  # back out and rerun
  setwd("../")
  print(paste("finished", dir_i))
}, error=function(e){})
}

# select and rename columns
names(out_all_indices)
out_all_indices <- out_all_indices %>% 
  select(common_name, scientific_name, 
         region, region_type, 
         count_number, year_start, year_end, 
         annual_index=index_median, 
         index_lcl=index_q0.025,  
         index_ucl=index_q0.975,   
         survey_suitability,
         range_proportion=range_coverage)
write.csv(out_all_indices, 
          file = file.path(output_dir, 
                           "na_viewer/cbc_indices_version_2023.0_na_web_download.csv"),
          na = "", row.names=F)

names(out_all_trends)
out_all_trends <- out_all_trends %>% 
  select(scientific_name, 
         region, region_type, 
         trend_type,
         year_start, year_end, 
         annual_percent_change=trend_median, 
         annual_change_lcl=trend_q0.025,  
         annual_change_ucl=trend_q0.975, 
         total_percent_change=percent_change_median,
         total_change_lcl=percent_change_q0.025,
         total_change_ucl=percent_change_q0.975,
         prob_decline_gt0p=prob_decline,
         prob_decline_gt30p=prob_decline_gt30,
         prob_decline_gt50p=prob_decline_gt50,
         prob_decline_gt70p=prob_decline_gt70,
         survey_suitability,
         range_proportion=range_coverage)
write.csv(out_all_trends, 
          file = file.path(output_dir, 
                           "na_viewer/cbc_trends_version_2023.0_na_web_download.csv"),
          na = "", row.names=F)
# ------------------------------------------------------------------------------




# quick comparison -------------------------------------------------------------
# dir()
# setwd(output_dir)
# d1 <- read.csv("./archive/cbc_trend_2022_v2.csv", stringsAsFactors = F) %>% 
#   filter(parameter=="RegressionTrend1970On") %>%
#   select(scientific_name=ebird_sci_name, region=stratum, 
#          med_2022=estimate_median) %>%
#   arrange(scientific_name, region)
# 
# d2 <- read.csv(file.path(output_dir, 
#                          "./na_viewer/cbc_trends_version_2023.0_na_web_download.csv"), 
#                stringsAsFactors = F) %>% 
#    filter(trend_type=="long-term",
#           region_type!="stratum") %>%
#   select(scientific_name, region, med_2024=trend_median)
# 
# d3 <- left_join(d2, d1)
# plot(d3$med_2022, d3$med_2024)
# summary(mod1 <- lm(d3$med_2024~d3$med_2022))
# plot(mod1)
# ------------------------------------------------------------------------------



# map stuff --------------------------------------------------------------------
# library(sf)
# map1 <- read_sf("cbc_strata_version_2023.0_na.gpkg") %>% 
#   mutate(bcr=paste0("BCR", str_pad(bcr, width=2, side = "left", pad=0)))
# plot(map1$geom)
# write_sf(map1, "cbc_strata_version_2023.0_na.shp")
