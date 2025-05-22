


# Script for checking the quality of results from the first model run.
# Tim Meehan, Quantitative Ecologist, National Audubon Society
# 21 April 2025 



# set up -----------------------------------------------------------------------
library(tidyverse)

# set paths for files and results
git_data_path <- "D:/Users/tim.meehan/Documents/GitHub/CBCTrendAnalysis"
results_root_directory <- "Z:/7_CommunityScience/CBCAnalysisResults/cbc_results_v2023.0"
species_list_file <- paste(git_data_path, "data/taxon_key_dec_2024.csv", 
                           sep="/")
options(stringsAsFactors=F)
setwd(results_root_directory)

# get species file
species_list <- read.csv(species_list_file) %>% 
  mutate(ebird_com_name=gsub("Â", "", ebird_com_name)) %>% 
  mutate(ebird_com_name=gsub("/", " or ", ebird_com_name))
# ------------------------------------------------------------------------------



# check results per species ----------------------------------------------------
# start data frame for results
res1 <- data.frame(ebird_com_name=species_list[,2], 
                   dir_exists=0, 
                   got_data=0, 
                   got_coverage=0, 
                   mapped_neighbors=0,
                   got_draws = 0,
                   got_par_sums = 0,
                   made_trend_maps = 0, 
                   got_indices = 0,
                   got_trends = 0, 
                   prop_bad_rh = NA,
                   prop_bad_ess = NA,
                   trend_median = NA,
                   trend_q0.025 = NA,
                   trend_q0.975 = NA,
                   trend_median_10yr = NA,
                   trend_q0.025_10yr = NA,
                   trend_q0.975_10yr = NA,
                   ebird_sci_name=gsub(" ", "_", species_list[,3]))

# loop through species
i <- 1
for(i in 1:nrow(species_list)){
  
  # set species code
  species_code <- gsub(" ", "_", species_list[i, 2])
  
  # set wd
  setwd(results_root_directory)
  
  # status
  print(paste("trying", species_code))
  
  # if dir exists
  if(dir.exists(species_code)) {
    
    # verify directory
    res1[i, 2] <- 1
    
    # go to that directory
    setwd(file.path(results_root_directory, species_code))
    
    # check existence of data
    if(file.exists(paste0(species_code, 
                          "_zero_filled_and_filtered_modeling_data.csv")) == 
       TRUE) {
      res1[i, 3] <- 1
    }
    
    # check existence of coverage proportions
    if(file.exists(paste0(species_code, "_range_coverage_proportion.csv")) == 
       TRUE) {
      res1[i, 4] <- 1
    }
    
    # check existence of neighborhood map
    if(file.exists(paste0(species_code, "_strata_map.pdf")) == 
       TRUE) {
      res1[i, 5] <- 1
    }
    
    # check existence of posterior draws
    if(file.exists(paste0("posterior_draws_", species_code, 
                          "_CBC_spatial_first_diff.rds")) == 
       TRUE) {
      res1[i, 6] <- 1
    }
    
    # check existence of parameter summaries
    if(file.exists(paste0(species_code, "_parameter_estimate_summaries.csv")) == 
       TRUE) {
      res1[i, 7] <- 1
    }
    
    # check existence of stratum trend maps
    if(file.exists(paste0(species_code, "_stratum_trend_map.pdf")) == 
       TRUE) {
      res1[i, 8] <- 1
    }
    
    # check existence of abundance indices
    if(file.exists(paste0(species_code, "_relative_abundance_indices.csv")) == 
       TRUE) {
      res1[i, 9] <- 1
    }
    
    # check existence of abundance trends
    if(file.exists(paste0(species_code, "_relative_abundance_trends.csv")) == 
       TRUE) {
      res1[i, 10] <- 1
    }
    
    # look in parameter summary file
    if(file.exists(paste0(species_code, "_parameter_estimate_summaries.csv")) == 
       TRUE) {
      
      # load error file
      fname <- paste(results_root_directory, species_code, 
                     paste0(species_code,
                     "_parameter_estimate_summaries.csv"), 
                     sep="/")
      log_file <- read.csv(fname)
      
      # get rhat
      rh1 <- log_file$rhat
      prop_bad_rh <- sum(as.numeric(rh1 > 1.05)) / length(rh1)
      res1[i, 11] <- prop_bad_rh
      
      # get ess
      ess1 <- log_file$ess_bulk
      prop_bad_ess <- sum(as.numeric(ess1 < 200)) / length(ess1)
      res1[i, 12] <- prop_bad_ess
    }
    
    # grab abundance trend
    if(file.exists(paste0(species_code, "_relative_abundance_trends.csv")) == 
       TRUE) {
      
      trd_file <- read.csv(paste0(species_code, 
                                  "_relative_abundance_trends.csv"))
      trds <- trd_file %>% 
        filter(region_type=="continent", trend_type=="long-term") %>% 
        select(trend_median, trend_q0.025, trend_q0.975) %>% 
        as.numeric()
      res1[i, 13:15] <- trds
      trds2 <- trd_file %>% 
        filter(region_type=="continent", trend_type=="ten-year") %>% 
        select(trend_median, trend_q0.025, trend_q0.975) %>% 
        colMeans() %>% 
        as.numeric()
      res1[i, 16:18] <- trds2

    }
    
  } else {next}
  
}

# summary metric
res1$sum_success <- rowSums(res1[,3:8])

# save
write.csv(res1, 
          paste0(git_data_path, "/data/", "first_results_check.csv"), 
          row.names=F)
# ------------------------------------------------------------------------------



# compare trends with old ones -------------------------------------------------
# set up directories
old_results <- "E:/7_CommunityScience/CBCAnalysisResults/archived_results/cbc_results_v2021.1/analysis_output"

# new results
res1 <- read.csv(paste0(git_data_path, "/data/", "first_results_check.csv"))
res2 <- res1 %>% mutate(old_trend_median=NA, old_trend_q0.025=NA, 
                        old_trend_q0.975=NA,
                        old_trend_median_10yr=NA, old_trend_q0.025_10yr=NA, 
                        old_trend_q0.975_10yr=NA)

# loop through species
for(i in 1:nrow(res2)){
  print(paste("trying", res2$ebird_com_name[i]))
  tryCatch({
    dir_i <- file.path(old_results, res1$ebird_sci_name[i])
    file_i <- list.files(dir_i, pattern="aggregate_estimates_summary.csv", 
                         full.names = T)
    dat_i <- read.csv(file_i) 
    dat1 <- dat_i %>% filter(parameter=="RegressionTrend1970On", stratum=="USACAN") 
    dat2 <- dat_i %>% filter(parameter=="RegressionTrend10Year", stratum=="USACAN") 
    res2[i, 21] <- dat1$estimate_median
    res2[i, 22] <- dat1$estimate_lcl
    res2[i, 23] <- dat1$estimate_ucl
    res2[i, 24] <- dat2$estimate_median
    res2[i, 25] <- dat2$estimate_lcl
    res2[i, 26] <- dat2$estimate_ucl
    print(paste("finished", res2$ebird_com_name[i]))
  }, error=function(e){})
}

# convert trends
res3 <- res2 %>% 
  mutate(old_trend_median=(old_trend_median-1)*100,
  old_trend_q0.025=(old_trend_q0.025-1)*100,
  old_trend_q0.975=(old_trend_q0.975-1)*100,
  old_trend_median_10yr=(old_trend_median_10yr-1)*100,
  old_trend_q0.025_10yr=(old_trend_q0.025_10yr-1)*100,
  old_trend_q0.975_10yr=(old_trend_q0.975_10yr-1)*100) %>% 
  mutate(new_prec=trend_q0.975-trend_q0.025,
         old_prec=old_trend_q0.975-old_trend_q0.025,
         prec_diff=new_prec-old_prec,
         trend_diff=trend_median-old_trend_median,
         new_prec_10yr=trend_q0.975_10yr-trend_q0.025_10yr,) %>% 
  mutate(trend_sig = ifelse(trend_q0.025 > 0 | trend_q0.975 < 0, 1, 0),
         trend10_sig = ifelse(trend_q0.025_10yr > 0 | trend_q0.975_10yr < 0, 1, 0))

# save
write.csv(res3, 
          paste0(git_data_path, "/data/", "first_results_check.csv"), 
          row.names=F)
# ------------------------------------------------------------------------------




# quick summaries --------------------------------------------------------------
# plot stuff
res3 <- read.csv(paste0(git_data_path, "/data/", "first_results_check.csv"))
names(res3)
res3 %>% ggplot() +
  geom_text(aes(x=trend_median, y=old_trend_median, label=ebird_com_name)) +
  geom_abline(intercept=0, slope=1)
res3 %>% ggplot() +
  geom_point(aes(x=trend_median, y=old_trend_median)) +
  geom_abline(intercept=0, slope=1)

# correlation
cor(cbind(res3$trend_median, res3$old_trend_median), use="pairwise.complete",
    method="pearson")

# new trend
summary(res3$trend_median)
summary(res3$trend_median_10yr)
res3 %>% filter(trend_sig==1) %>% select(trend_median) %>% summary()
res3 %>% filter(trend_sig==1) %>% ggplot() +
  geom_histogram(aes(x=trend_median))
res3 %>% filter(trend_sig==1) %>% select(trend_median_10yr) %>% summary()
res3 %>% filter(trend_sig==1) %>% ggplot() +
  geom_histogram(aes(x=trend_median_10yr))

# trend diffs
summary(res3$trend_diff)
res3 %>% ggplot() +
  geom_histogram(aes(x=trend_diff))

# prec diffs
summary(res3$prec_diff)
res3 %>% ggplot() +
  geom_histogram(aes(x=prec_diff))

# lt increases and decreases
res3 %>% filter(!is.na(trend_median)) %>% nrow()
res3 %>% filter(trend_sig==0) %>% nrow()
res3 %>% filter(trend_sig==1, trend_median>0) %>% nrow()
res3 %>% filter(trend_sig==1, trend_median<0) %>% nrow()

# st increases and decreases
res3 %>% filter(!is.na(trend_median_10yr)) %>% nrow()
res3 %>% filter(trend10_sig==0) %>% nrow()
res3 %>% filter(trend10_sig==1, trend_median_10yr>0) %>% nrow()
res3 %>% filter(trend10_sig==1, trend_median_10yr<0) %>% nrow()
# ------------------------------------------------------------------------------



# more divergence analysis -----------------------------------------------------
names(res3)
res3 %>% mutate(perc_diff_trend=(trend_diff/old_trend_median)*100) %>% 
  select(ebird_com_name, ebird_sci_name, trend_median, 
         old_trend_median, trend_diff, perc_diff_trend) %>% 
  arrange(desc(trend_diff)) %>% View()
# ------------------------------------------------------------------------------


# bring in acad ----------------------------------------------------------------
acad0 <- read.csv("D:/Users/tim.meehan/Box/-.Science Team shared/2_Projects/0_Climate/Climate 2.0/ACAD PIF/Analysis/ACAD_SBD_join.csv")
acad1 <- acad0 %>% 
  select(ebird_com_name=Common.Name,family, order, 
          cbd_group=Group.x, 
          breed_biome=Breeding.Biome,
          winter_biome=Nonbreeding.Biome,
          breed_habitat=Primary.Breeding.Habitat,
          winter_habitat=Primary.Nonbreeding.Habitat,
          mig_status=Mig.Status,
          glob_pop_score=PS.g,
          breed_range_score=BD.g,
          winter_range_score=ND.g,
          breed_threats_score=TB.c,
          winter_threats_score=TN.c,
          breed_vuln_score=CV.b,
          winter_vuln_score=CV.w,
          breed_trend_score=PT.c,
          breed_tot_score=CCS.b,
          winter_tot_score=CCS.n,
          urban=Urban,
          ag=Agriculture) %>% 
  mutate(breed_habitat=str_split(breed_habitat, ": ", simplify = T)[,2],
         winter_habitat=str_split(winter_habitat, ": ", simplify = T)[,2]) %>% 
  mutate(urban=ifelse(urban=="yes", 1, 0),
         ag=ifelse(ag %in% c("w", "b,w", "b"), 1, 0))

# join
d1 <- res3 %>% left_join(acad1)
summary(d1)
# ------------------------------------------------------------------------------




# models stats -----------------------------------------------------------------
# library(MuMIn)
library(cowplot)
library(car)
library(performance)
library(glmmTMB)
library(emmeans)
setwd(git_data_path)

############### bird group 
# lt models
names(d1)
m1 <- glmmTMB(trend_median ~ cbd_group + (1|order),  
              data=d1)
summary(m1)
Anova(m1, type="III")
performance(m1)
EMM <- emmeans(m1, "cbd_group"); EMM
summary(EMM, type = "response", adjust="tukey")
p1 <- plot(EMM, type = "response", comparisons = TRUE) + 
  labs(x="Predicted long-term trend",
       y="Bird group") +
  theme_bw() +
  geom_vline(xintercept=0, col="darkred"); p1

# st models
m1 <- glmmTMB(trend_median_10yr ~ cbd_group + (1|order),
              data=d1)
summary(m1)
Anova(m1, type="III")
performance(m1) 
EMM <- emmeans(m1, "cbd_group"); EMM
summary(EMM, type = "response", adjust="tukey")
p2 <- plot(EMM, type = "response", comparisons = TRUE) + 
  labs(x="Predicted short-term trend",
       y="Bird group") +
  theme_bw() +
  geom_vline(xintercept=0, col="darkred")

# plots
plot_grid(p1, p2)
ggsave("./output/bird_group_and_trends.pdf", width=220, height=140, units = "mm", dpi=300)
getwd()


############### migration 
# lt models
names(d1)
summary(d1)
m1 <- glmmTMB(trend_median ~ mig_status + (1|order),  
              data=d1)
summary(m1)
Anova(m1, type="III")
performance(m1)
EMM <- emmeans(m1, "mig_status"); EMM
summary(EMM, type = "response", adjust="tukey")
p1 <- plot(EMM, type = "response", comparisons = TRUE) + 
  labs(x="Predicted long-term trend",
       y="Migration status") +
  theme_bw() +
  geom_vline(xintercept=0, col="darkred"); p1

# st models
m1 <- glmmTMB(trend_median_10yr ~ mig_status + (1|order),
              data=d1)
summary(m1)
Anova(m1, type="III")
performance(m1) 
EMM <- emmeans(m1, "mig_status"); EMM
summary(EMM, type = "response", adjust="tukey")
p2 <- plot(EMM, type = "response", comparisons = TRUE) + 
  labs(x="Predicted short-term trend",
       y="Migration status") +
  theme_bw() +
  geom_vline(xintercept=0, col="darkred"); p2

# plots
plot_grid(p1, p2)
ggsave("./output/migration_status_and_trends.pdf", width=220, height=140, units = "mm", dpi=300)


############### breed trend 
# lt models
names(d1)
summary(d1)
m1 <- glmmTMB(trend_median ~ factor(breed_trend_score) + (1|order),  
              data=d1)
summary(m1)
Anova(m1, type="III")
performance(m1)
EMM <- emmeans(m1, "breed_trend_score"); EMM
summary(EMM, type = "response", adjust="tukey")
p1 <- plot(EMM, type = "response", comparisons = TRUE) + 
  labs(x="Predicted long-term trend",
       y="Breeding trend score") +
  theme_bw() +
  geom_vline(xintercept=0, col="darkred"); p1

# st models
m1 <- glmmTMB(trend_median_10yr ~ factor(breed_trend_score) + (1|order),
              data=d1)
summary(m1)
Anova(m1, type="III")
performance(m1) 
EMM <- emmeans(m1, "breed_trend_score"); EMM
summary(EMM, type = "response", adjust="tukey")
p2 <- plot(EMM, type = "response", comparisons = TRUE) + 
  labs(x="Predicted short-term trend",
       y="Breeding trend score") +
  theme_bw() +
  geom_vline(xintercept=0, col="darkred"); p2

# plots
plot_grid(p1, p2)
ggsave("./output/breeding_trend_score_and_trends.pdf", width=220, height=140, units = "mm", dpi=300)


############### breed climate vuln 
# lt models
names(d1)
summary(d1)
m1 <- glmmTMB(trend_median ~ factor(breed_vuln_score) + (1|order),  
              data=d1)
summary(m1)
Anova(m1, type="III")
performance(m1)
EMM <- emmeans(m1, "breed_vuln_score"); EMM
summary(EMM, type = "response", adjust="tukey")
p1 <- plot(EMM, type = "response", comparisons = TRUE) + 
  labs(x="Predicted long-term trend",
       y="Breeding climate vulnerability score") +
  theme_bw() +
  geom_vline(xintercept=0, col="darkred"); p1

# st models
m1 <- glmmTMB(trend_median_10yr ~ factor(breed_vuln_score) + (1|order),
              data=d1)
summary(m1)
Anova(m1, type="III")
performance(m1) 
EMM <- emmeans(m1, "breed_vuln_score"); EMM
summary(EMM, type = "response", adjust="tukey")
p2 <- plot(EMM, type = "response", comparisons = TRUE) + 
  labs(x="Predicted short-term trend",
       y="Breeding climate vulnerability score") +
  theme_bw() +
  geom_vline(xintercept=0, col="darkred"); p2

# plots
plot_grid(p1, p2)
ggsave("./output/breeding_climate_score_and_trends.pdf", width=220, height=140, units = "mm", dpi=300)


############### breed threats
# lt models
names(d1)
summary(d1)
m1 <- glmmTMB(trend_median ~ factor(breed_threats_score) + (1|order),  
              data=d1)
summary(m1)
Anova(m1, type="III")
performance(m1)
EMM <- emmeans(m1, "breed_threats_score"); EMM
summary(EMM, type = "response", adjust="tukey")
p1 <- plot(EMM, type = "response", comparisons = TRUE) + 
  labs(x="Predicted long-term trend",
       y="Breeding threats score") +
  theme_bw() +
  geom_vline(xintercept=0, col="darkred"); p1

# st models
m1 <- glmmTMB(trend_median_10yr ~ factor(breed_threats_score) + (1|order),
              data=d1)
summary(m1)
Anova(m1, type="III")
performance(m1) 
EMM <- emmeans(m1, "breed_threats_score"); EMM
summary(EMM, type = "response", adjust="tukey")
p2 <- plot(EMM, type = "response", comparisons = TRUE) + 
  labs(x="Predicted short-term trend",
       y="Breeding threats score") +
  theme_bw() +
  geom_vline(xintercept=0, col="darkred"); p2

# plots
plot_grid(p1, p2)
ggsave("./output/breed_threats_score_and_trends.pdf", width=220, height=140, units = "mm", dpi=300)
