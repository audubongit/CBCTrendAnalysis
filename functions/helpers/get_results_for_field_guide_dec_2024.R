library(dplyr)
library(tidyr)
library(stringr)
output_dir <- "D:/Users/tim.meehan/Box/-.Tim.Meehan individual/cbc_results_for_field_guide"
trends_dir <- "Z:/7_CommunityScience/CBCAnalysisResults/cbc_results_v2023.0_na"
code_dir <- "D:/Users/tim.meehan/Documents/GitHub/CBCTrendAnalysis/data"

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
  
  # get abundance indices
  fi <- list.files(pattern="relative_abundance_indices")
  f1 <- read.csv(fi) 
  f2 <- f1 %>% 
    filter(region_type %in% c("continent", "prov_state"),
           index_type=="smooth_first_diff") %>%
    select(common_name, scientific_name, region, 
           year_start, year_end, index_median, index_q0.025, index_q0.975) %>% 
    arrange(common_name, scientific_name, region, year_start)
  
  # get abundance trends
  fi <- list.files(pattern="relative_abundance_trends")
  f1 <- read.csv(fi) 
  f3 <- f1 %>% 
    filter(region_type %in% c("continent", "prov_state"),
           trend_type %in% c("long-term", "ten-year")) %>%
    select(common_name, scientific_name, region, trend_type,
           trend_diff_zero_95) %>%
    pivot_wider(names_from=trend_type, values_from=trend_diff_zero_95,
                names_prefix = "sig_", values_fn=mean) %>% 
    rename(sig_long_term="sig_long-term",
           sig_ten_year="sig_ten-year") %>% 
    ungroup()
  f4 <- f2 %>% left_join(f3)

  # stack it
  out_all <- rbind(out_all, f4)
  
  # back out and rerun
  setwd("../")
  print(paste("finished", spp_i))
  
  }, error=function(e){})
}


# save
write.csv(out_all, file.path(output_dir, "cbc_time_series_for_field_guide.csv"), 
          row.names=F, na="")







