# load libraries
library(igraph)
library(tidyverse)

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(dirname("/projects/lude8513/ccm_networks/code"))

# source helper functions
source("edm_utils.R")

# read in data
edge_lists <- read_csv("../data/edge_lists.csv")
gom_raw <- read_csv("../data/experimental_data_stacked.csv")

network_stats <- bind_rows(lapply(unique(gom_raw$site), function(s){
  # wrangle dynamics data
  gom_clean <- gom_raw %>%
    # just 1 site and only control plots
    filter(site == s, plot == "C") %>%
    # exclude mobile species
    filter(metric_type == "percent cover") %>%
    # combine same-species observations within survey groups (eg epibiont and encrusting)
    group_by(survey_group, species) %>%
    summarize(value_scaled = sum(value_scaled)) %>%
    # drop unneeded columns
    select(c(survey_group, species, value_scaled)) %>%
    # sum cover of all species
    group_by(survey_group) %>%
    summarize(total_cover = sum(value_scaled))
  
  # calculate temporal stability for current site
  temporal_stability <- gom_clean %>%
    summarize(cv = sd(total_cover)/mean(total_cover))
  
  # build network from edge list
  spp_network <- graph_from_edgelist(as.matrix(filter(edge_lists, site == s)[, 1:2]),
                                     directed = TRUE)
  
  # calculate network stats
  network_stats <- calc_network_stats(spp_network) %>%
    cbind(temporal_stability) %>%
    mutate(site = s)
}))

# visualize networks
spp_network <- graph_from_edgelist(as.matrix(filter(edge_lists, site == 4)[, 1:2]),
                                   directed = TRUE)
plot(spp_network)

# fit models ----
model <- lm(cv ~ mean_in_deg + connectance + relative_ascendancy, data = network_stats, na.action = na.omit)
summary(model)






  
  
  
