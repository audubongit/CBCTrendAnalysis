
# setup ------------------------------------------------------------------------
library(tidyverse)
output_dir <- "D:/Users/tim.meehan/Box/-.Tim.Meehan individual/CBCForBirdsCanada"
trends_dir <- "Z:/7_CommunityScience/CBCAnalysisResults/cbc_results_v2023.0_na"
code_dir <- "D:/Users/tim.meehan/Documents/GitHub/CBCTrendAnalysis/functions/helpers"
# ------------------------------------------------------------------------------



# for birds canada -------------------------------------------------------------
setwd(code_dir)
spp_tab <- read.csv("../../data/taxon_key_dec_2024.csv") %>%
  dplyr::select(2,3,1) 

# get dir names
setwd(trends_dir)
dir_lst <- list.dirs(recursive = F)

# loop through species 
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
  
  trds1 <- f2 %>% mutate(results_code="CBC", version=2023,
                         season="non-breeding", index_type="individuals",
                         model_type="spatial-first-difference",
                         prob_increase_0=1-prob_decline) %>% 
    select(results_code, version, area_code=region, season,
           period=trend_type, species_code, 
           year_start, year_end, trnd=trend_median, 
           lower_ci=trend_q0.025, upper_ci=trend_q0.975, index_type, 
           model_type, percent_change=percent_change_median,
           percent_change_low=percent_change_q0.025, 
           percent_change_high=percent_change_q0.975,
           prob_decrease_0=prob_decline, 
           prob_decrease_30=prob_decline_gt30,
           prob_decrease_50=prob_decline_gt50,
           prob_increase_0,
           suitability=survey_suitability,
           precision_number=estimate_precision,
           coverage_num=range_coverage
    )
  
  # get indices
  fi <- list.files(pattern="relative_abundance_indices")
  f3 <- read.csv(fi)
  idxs1 <- f3 %>% filter(index_type=="first_diff") %>% 
    mutate(results_code="CBC", version=2023, 
           season="non-breeding", period="long-term") %>% 
    select(results_code, version, area_code=region, 
           season, period, species_code, year=year_start,
           index=index_median, lower_ci=index_q0.025, upper_ci=index_q0.975)
  idxs2 <- f3 %>% filter(index_type=="smooth_first_diff") %>% 
    mutate(results_code="CBC", version=2023, 
           season="non-breeding", period="long-term") %>% 
    select(results_code, version, area_code=region, 
           season, period, species_code, year=year_start,
           LOESS_index=index_median, trend_index=index_median, 
           trend_index_lower_ci=index_q0.025, 
           trend_index_upper_ci=index_q0.975)
  idxs1_wide <- idxs1 %>% left_join(idxs2)  
    
  # put it together and stack
  out_all_trends <- rbind(out_all_trends, trds1)
  out_all_indices <- rbind(out_all_indices, idxs1_wide)
  
  # back out and rerun
  setwd("../")
  print(paste("finished", dir_i))
}, error=function(e){})
}

# tests
length(unique(out_all_indices$species_code))
length(unique(out_all_trends$species_code))
summary(out_all_indices)
summary(out_all_trends)

# write csvs
write.csv(out_all_indices, 
          file = file.path(output_dir, 
                           "na_cbc_indices_version_2023.0_birds_canada.csv"),
          na = "", row.names=F)


write.csv(out_all_trends, 
          file = file.path(output_dir, 
                           "na_cbc_trends_version_2023.0_birds_canada.csv"),
          na = "", row.names=F)
# ------------------------------------------------------------------------------


