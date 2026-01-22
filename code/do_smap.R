# load libraries
library(tidyverse)
library(igraph)
library(rEDM)
library(janitor)
library(furrr)

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# source helper functions ----
source("edm_utils.R")

# read in data ----
gom_raw <- read_csv("../data/experimental_data_stacked.csv")
edge_lists <- read_csv("../data/edge_lists.csv")

# iterate over all sites
all_smap_coefs <- bind_rows(lapply(unique(edge_lists$site), function(s){
  print(paste0("site ", s))
  
  # fit multivariate S-map to one target species
  this_site <- filter(edge_lists, site == s)
  
  # prep timeseries data
  this_site_wide <- gom_raw %>%
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
  
  # now do CCM
  smap_data <- this_site_wide %>%
    ungroup() %>%
    clean_names() %>%
    arrange(date) %>%
    relocate(date)
  
  # iterate over species
  plan(multisession)
  this_site_coefs <- bind_rows(future_map(seq_along(unique(this_site$sp2)), function(sp){
    # pick target species and get its predictors
    this_target <- unique(this_site$sp2)[sp]
    # check if target is in data
    valid_spp <- names(smap_data)
    # if target was dropped, skip
    if (!(this_target %in% valid_spp)) {
      return(NULL)
    }
    these_drivers <- filter(this_site, sp2 == this_target)$sp1
    valid_drivers <- these_drivers[these_drivers %in% valid_spp]
    # if no valid drivers, skip
    if (length(these_drivers) == 0) {
      return(NULL)
    }
    # fit smap
    do_smap(smap_data, this_target, valid_drivers) %>%
      mutate(site = s)
  }))
  plan(sequential)
  
  this_site_coefs
}))

# plot interactions strengths
ggplot(data=filter(all_smap_coefs, site == 1), aes(x=date, y = value, col = name)) +
  geom_line() +
  labs(x = "Date", y = "Interaction Strength") +
  theme_classic() +
  theme(legend.position = "none")
  #facet_wrap(~site)





