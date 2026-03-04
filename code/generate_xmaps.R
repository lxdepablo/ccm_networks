# load libraries
library(tidyverse)
library(odin)
#library(GPEDM) # for empirical dynamic modeling (devtools::install_github("tanyalrogers/GPEDM"))
library(rEDM)
library(janitor)

# set working directory
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("/projects/lude8513/ccm_networks/code/")
print(getwd())

# load helper functions
source("edm_utils.R")

# load GOM data
gom_raw <- read_csv("../data/experimental_data_stacked.csv")
# load SST anomaly data
sst <- read_csv("../data/sst/sst.csv")

# visualize temperature data
ggplot(data = sst, aes(x = as.Date(date), y = sst)) +
  geom_line()

# iterate over every site
all_sites_xmaps <- bind_rows(lapply(unique(gom_raw$site), function(s){
  print(paste0("site ", s))
  
  # prep timeseries data
  curr_site_wide <- gom_raw %>%
    # use only control plots
    filter(site == s, plot == "C") %>%
    # exclude double-counted species
    filter(!(metric_type == "count" & species == "MYED")) %>%
    # combine same-species observations within survey groups
    group_by(survey_group, species) %>%
    summarize(value_scaled = sum(value_scaled, na.rm = TRUE), .groups = "drop") %>%
    # add time column
    mutate(date = my(survey_group)) %>%
    arrange(date) %>%
    select(-survey_group) %>%
    # drop species with < n observations
    group_by(species) %>%
    filter(sum(value_scaled != 0, na.rm = TRUE) >= 3) %>%
    ungroup() %>%
    # backfill zeroes before scaling: make all species x dates explicit
    complete(date, species, fill = list(value_scaled = 0)) %>%
    # z score each species
    group_by(species) %>%
    mutate(value_z = as.numeric(scale(value_scaled))) %>%
    ungroup() %>%
    # guard against species with zero variance
    mutate(value_z = replace_na(value_z, 0)) %>%
    select(date, species, value_z) %>%
    # make wide (use value_z)
    pivot_wider(
      id_cols     = date,
      names_from  = species,
      values_from = value_z
    ) %>%
    # backfill missing species observations with 0's
    replace(is.na(.), 0)
  
  # bring in temperature data
  sst_z <- sst %>%                       # <-- rename to your SST dataframe object
    mutate(date = as.Date(date)) %>%        # ensure Date class
    select(date, sst) %>%
    arrange(date) %>%
    mutate(
      sst = as.numeric(sst),
      sst_z = as.numeric(scale(sst)),
      sst_z = ifelse(is.na(sst_z), 0, sst_z)
    ) %>%
    select(date, sst = sst_z)
  
  curr_site_wide <- curr_site_wide %>%
    left_join(sst_z, by = "date") %>%
    mutate(sst = replace_na(sst, 0))
  
  # # plot
  # curr_site_long <- curr_site_wide %>%
  #   pivot_longer(cols = -date)
  # ggplot(data = curr_site_long, aes(x = date, y = value, col = name)) +
  #   geom_line()
  
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







