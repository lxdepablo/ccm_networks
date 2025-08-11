# load libraries
library(tidyverse)
library(odin)
#library(GPEDM) # for empirical dynamic modeling (devtools::install_github("tanyalrogers/GPEDM"))
library(rEDM)
library(janitor)

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("/projects/lude8513/ccm_networks/code/")
print(getwd())

# load helper functions
source("edm_utils.R")

# load GOM data
gom_raw <- read_csv("../data/experimental_data_stacked.csv")

# iterate over every site
all_sites_xmaps <- bind_rows(lapply(unique(gom_raw$site), function(s){
  # wrangle data
  curr_site_wide <- gom_raw %>%
    # use only control plots
    filter(site == s, plot == "C") %>%
    # exclude double-counted species
    filter(!(metric_type == "count" & species == "MYED")) %>%
    # combine same-species observations within survey groups
    group_by(survey_group, species) %>%
    summarize(value_scaled = sum(value_scaled)) %>%
    # drop unneeded columns
    select(c(survey_group, species, value_scaled)) %>%
    # make wide
    pivot_wider(id_cols = survey_group, names_from = species, values_from = value_scaled) %>%
    # add time column
    mutate(date = my(survey_group))
  
  # now do CCM
  ccm_data <- curr_site_wide %>%
    ungroup() %>%
    select(-survey_group) %>%
    clean_names() %>%
    #rename(FUDI_canopy = "FUDI (canopy)", FUDI_attachment = "FUDI (attachment)") %>%
    arrange(date) %>%
    relocate(date)
  
  # do CCM on every pair of species
  all_xmaps_long <- par_calc_all_xmaps(ccm_data) %>%
    mutate(site = s)
}))


# write CSV's
write_csv(all_sites_xmaps, "../data/xmaps.csv")







