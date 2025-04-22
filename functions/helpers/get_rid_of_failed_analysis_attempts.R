# define some directories
code_dir <- "C:/Users/tmeehan/Documents/GitHub/CBCTrendAnalysis"
results_dir <- "Z:/7_CommunityScience/CBCAnalysisResults/cbc_results_v2023.0"

# input species table
setwd(code_dir)
species_table <- read.csv(file.path(code_dir, "/data/taxon_key_dec_2024.csv")) %>% 
  mutate(ebird_com_name=gsub("Â", "", ebird_com_name)) %>% 
  mutate(ebird_com_name=gsub("/", " or ", ebird_com_name))

# define vector for looping
ebird_com_names <- species_table$ebird_com_name

# loop through species
# s <- 1
for(s in 1:nrow(species_table)){ # start for loop
  
  # define species variables
  ebird_com_name_s <- ebird_com_names[s]
  
  # define output location
  dir_out1 <- file.path(results_dir, gsub(" ", "_", ebird_com_name_s))
  
  # skip if the species has been finished
  if(file.exists(paste0(dir_out1, "/", gsub(" ", "_", ebird_com_name_s), 
                        "_stratum_trend_map.pdf"))){
    next
  }
    
  # if not delete dir
  if(!file.exists(paste0(dir_out1, "/", gsub(" ", "_", ebird_com_name_s), 
                         "_stratum_trend_map.pdf"))){
    unlink(dir_out1, recursive = T, force = T)
  }
}
  
