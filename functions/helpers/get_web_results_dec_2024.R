library(dplyr)
library(tidyr)
library(stringr)
#library(ebirdst)
#library(ggplot2)
output_dir <- "D:/Users/tmeehan/Box/-.Tim.Meehan individual/CBCForWeb"
trends_dir <- "Z:/7_CommunityScience/CBCAnalysisResults/CBC_Draft_4.0_1966-2021_2021/analysis_output"
code_dir <- "D:/Users/tmeehan/Documents/GitHub/QuantitativeMetrics/cbc_analysis_library/code_2021/helpers/get_web_results"

# get species names
setwd(code_dir)
spp_tab <- read.csv("cbc_trends_spp_list_66_to_121.csv") %>%
  dplyr::select(1:3, 6, 8:12) 
spp_tab <- spp_tab %>% dplyr::select(1,3) %>%
  mutate(ebird_sci_name=gsub("_", " ", species_code)) %>% select(3,2)

# get old web file
setwd(output_dir)
web_tab <- read.csv("cbc_trend_2020.csv")
strat_key <- web_tab %>% select(3,2,4,5) %>% distinct()
spp_key <- web_tab %>% select(6,7,70,71) %>% distinct() %>%
  mutate(ebird_sci_name=gsub("/"," ",ebird_sci_name))
spp_key <- spp_tab %>% left_join(spp_key)

# get dir names
setwd(trends_dir)
dir_lst <- list.dirs(recursive = F)

# loop through dirs to get web app format data
out_all <- c()
i <- 538
for(i in 1:length(dir_lst)){
  spp_i <- gsub("./", "", as.character(dir_lst[i]))
  spp_i <- gsub("_", " ", as.character(spp_i))
  print(paste("trying", spp_i))
  tryCatch({
  setwd(trends_dir)
  setwd(dir_lst[i])
  fi <- list.files(pattern="aggregate_estimates_summary")
  f2 <- read.csv(fi) 
  # get trend estimates
  f3 <- f2 %>% filter(grepl("Regression", parameter)) %>%
    mutate(ebird_sci_name=spp_i) %>% 
    select(ebird_sci_name, 1,2,4,5,6) %>%
    mutate(estimate_median=(estimate_median - 1) * 100,
           estimate_lcl=(estimate_lcl - 1) * 100,
           estimate_ucl=(estimate_ucl - 1) * 100) # %>%
    # mutate(Geography=NA,
    #        Geography=ifelse(stratum=="USACAN", "Region", Geography),
    #        Geography=ifelse(grepl("BCR", stratum), "BCR", Geography),
    #        Geography=ifelse(str_count(stratum)==2, "State/Province", Geography),
    #        Geography=ifelse(stratum=="USA", "Country", Geography),
    #        Geography=ifelse(stratum=="CAN", "Country", Geography))
  
  # get abundance indices
  f4 <- f2 %>% filter(parameter=="AbundanceIndex", stratum=="USACAN") %>%
    select(count_year, estimate_median) %>%
    pivot_wider(names_from=count_year, values_from = estimate_median,
                names_prefix="F") %>%
    mutate(ValType="EstMed")
  
  # get global trends
  f5 <- f2 %>% filter(parameter=="RegressionTrend1970On", stratum=="USACAN") %>%
    select(4,5,6) %>%
    mutate_if(is.numeric, function(x) (x - 1) * 100)
  names(f5) <- c("region_estimate_median", "region_estimate_lcl", 
                 "region_estimate_ucl")
  
  # put it together and stack
  f6 <- cbind(f5, f4) %>% slice(rep(1:n(), each = nrow(f3)))
  f7 <- cbind(f3, f6) 
  out_all <- rbind(out_all, f7)
  
  # back out and rerun
  setwd("../")
  print(paste("finished", dir_i))
}, error=function(e){})
}

# add strat key
f8 <- left_join(out_all, strat_key)

# add spp key
f9 <- left_join(f8, spp_key)

# save
write.csv(f9, file.path(output_dir, "cbc_trend_2022_v2.csv"))


# loop through dirs to get downloadable data format
out_all_dl <- c()
i <- 538
for(i in 1:length(dir_lst)){
  spp_i <- gsub("./", "", as.character(dir_lst[i]))
  spp_i <- gsub("_", " ", as.character(spp_i))
  print(paste("trying", spp_i))
  tryCatch({
  setwd(trends_dir)
  setwd(dir_lst[i])
  fi <- list.files(pattern="aggregate_estimates_summary")
  f2 <- read.csv(fi) 
  # get trend estimates
  f3 <- f2 %>% 
    mutate(ebird_sci_name=spp_i) 
  
  # take out indices
  f4 <- f3 %>% filter(parameter=="AbundanceIndex" | 
                        parameter=="AbundanceScalingFactor")
  
  
  f5 <- f3 %>% filter(parameter!="AbundanceIndex", 
                        parameter!="AbundanceScalingFactor") %>%
    mutate(estimate_mean=(estimate_mean - 1) * 100,
           estimate_median=(estimate_median - 1) * 100,
           estimate_lcl=(estimate_lcl - 1) * 100,
           estimate_ucl=(estimate_ucl - 1) * 100,
           quantile050=(quantile050 - 1) * 100,
           quantile165=(quantile165 - 1) * 100,     
           quantile835=(quantile835 - 1) * 100,    
           quantile950=(quantile950 - 1) * 100) 
  f6 <- rbind(f5, f4) %>%
    select(ebird_sci_names=ebird_sci_name, everything(), 
           -estimate_sd)
  
  # put it together and stack
  out_all_dl <- rbind(out_all_dl, f6)
  
  # back out and rerun
  setwd("../")
  print(paste("finished", dir_i))
}, error=function(e){})
}

# add spp key
f7 <- left_join(out_all_dl, spp_tab %>% rename(ebird_sci_names=ebird_sci_name)) %>%
  select(ebird_com_name, ebird_sci_name=ebird_sci_names, everything()) %>%
  rename(quantile_050=quantile050, quantile_165=quantile165, 
         quantile_835=quantile835, quantile_950=quantile950)

# save
write.csv(f7, file.path(output_dir, "cbc_trends_and_abundance_indices_v4.1_web_download_02Dec2022.csv"))

# quick check
dir()
setwd(output_dir)
d1 <- read.csv("cbc_trend_2020.csv", stringsAsFactors = F) %>% 
  filter(parameter=="RegressionTrend1970On") %>%
  select(ebird_sci_name, stratum, med_2020=estimate_median) %>%
  arrange(ebird_sci_name, stratum) %>%
  slice(-1) %>%
  mutate(mf=paste(ebird_sci_name, stratum)) %>%
  select(mf, med_2020)
d2 <- read.csv("cbc_trend_2022_v2.csv", stringsAsFactors = F) %>% 
   filter(parameter=="RegressionTrend1970On") %>%
  select(ebird_sci_name, stratum, med_2022=estimate_median,
         lcl_2022=estimate_lcl,
         ucl_2022=estimate_ucl) %>%
  arrange(ebird_sci_name, stratum) %>%
  mutate(mf=paste(ebird_sci_name, stratum)) %>%
  select(mf, med_2022, lcl_2022, ucl_2022)
d3 <- left_join(d2, d1)
d3$species <- paste(str_split(d3$mf, pattern=" ", simplify=T)[,1],
                    str_split(d3$mf, pattern=" ", simplify=T)[,2])
library(ggplot2)
library(lme4)
ggplot(d3, aes(x=med_2020, y=med_2022)) + geom_point() + geom_abline() +
  geom_smooth(method="lm")
ggplot(d3 %>% filter(grepl("USACAN", mf)), aes(x=med_2020, y=med_2022)) + 
  geom_point() + geom_abline() +
  geom_smooth(method="lm")
summary(lmer(med_2022~med_2020+(1|species) + (0+med_2020|species), d3))
summary(lm(med_2022~med_2020, d3 %>% filter(grepl("USACAN", mf))))
