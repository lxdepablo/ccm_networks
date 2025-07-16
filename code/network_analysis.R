# load libraries
library(igraph)
library(tidyverse)

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# source helper functions
source("edm_utils.R")

# read in data
edge_list <- read_csv("edge_list.csv")
gom_raw <- read_csv("experimental_data_stacked.csv")

# wrangle dynamics data
gom_clean <- gom_raw %>%
  # just 1 site and only control plots
  filter(site == 1, plot == "C") %>%
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

temporal_stability <- gom_clean %>%
  summarize(cv = sd(total_cover)/mean(total_cover))

# build network from edge list
spp_network <- graph_from_edgelist(as.matrix(edge_list[, 1:2]),
                                   directed = TRUE)
plot(spp_network)

# calculate network stats
network_stats <- calc_network_stats(spp_network) %>%
  cbind(temporal_stability)





  
  
  