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
metaweb_raw <- read_csv("../data/metaweb_edges.csv")

# wrangle metaweb
metaweb_clean <- metaweb_raw %>%
  dplyr::select(c(species_1.accepted.name, species_2.accepted.name)) %>%
  unique()

metaweb <- reverse_edges(graph_from_edgelist(as.matrix(filter(metaweb_clean)[, 1:2]),
                               directed = TRUE))

# calculate stats for causal networks by site
network_stats <- bind_rows(lapply(unique(gom_raw$site), function(s){
  # wrangle dynamics data
  gom_clean <- gom_raw %>%
    # just 1 site and only control plots
    filter(site == s, plot == "C") %>%
    # exclude mobile species
    #filter(metric_type == "percent cover") %>%
    # combine same-species observations within survey groups (eg epibiont and encrusting)
    group_by(survey_group, species) %>%
    summarize(value_scaled = sum(value_scaled)) %>%
    # drop unneeded columns
    dplyr::select(c(survey_group, species, value_scaled)) %>%
    # sum cover of all species
    group_by(survey_group) %>%
    summarize(total_cover = sum(value_scaled, na.rm = T))
  
  # calculate temporal stability for current site
  temporal_stability <- gom_clean %>%
    summarize(cv = sd(total_cover, na.rm=T)/mean(total_cover, na.rm = T))
  
  # build network from edge list
  spp_network <- graph_from_edgelist(as.matrix(filter(edge_lists, site == s)[, 1:2]),
                                     directed = TRUE)
  
  # calculate network stats
  network_stats <- calc_network_stats(spp_network) %>%
    cbind(temporal_stability) %>%
    mutate(site = s)
}))

# build causal metaweb
causal_metaweb <- unique(dplyr::select(edge_lists, c(sp1, sp2)))

# visualize networks
spp_network <- graph_from_edgelist(as.matrix(causal_metaweb[, 1:2]),
                                   directed = TRUE)

# calculate stats for metaweb and causal web
causal_metaweb_stats <- calc_network_stats(spp_network)
trophic_metaweb_stats <- calc_network_stats(metaweb)

# compare keystone species
# by degree
causal_degrees <- sort(degree(spp_network))
trophic_degrees <- sort(degree(metaweb))
# by page rank
causal_pageranks <- sort(page_rank(spp_network)$vector)
trophic_pageranks <- sort(page_rank(metaweb)$vector)

# Kamada–Kawai layout with extra spacing
layout <- layout_with_kk(spp_network) *1.5

# plot
plot(
  spp_network,
  layout = layout,
  vertex.size = 20,
  vertex.label = NA,
  vertex.color = "skyblue",
  edge.arrow.size = 0.1,
  edge.color = adjustcolor("grey40", alpha.f = 0.5),
  edge.curved = 0.15,
  margin = 0.05,
  asp = 0,
  rescale = F,
  xlim = c(-8,8),
  ylim = c(-7.2,7.2)
)

# plot subgraph
# pick the focal node
focal <- "fudi_canopy"

# get all neighbors (both in and out)
nbrs <- neighbors(spp_network, focal, mode = "all")

# include the focal node itself
nodes_to_keep <- c(focal, V(spp_network)[nbrs]$name)

# make subgraph
sub_net <- induced_subgraph(spp_network, vids = nodes_to_keep)

# layout: keep it tidy with KK
kk_layout <- layout_with_kk(sub_net)
#kk_layout <- norm_coords(kk_layout, xmin = -2, xmax = 2, ymin = -2, ymax = 2)

# plot
plot(
  sub_net,
  layout = kk_layout,
  vertex.size = 20,
  vertex.label = NA,
  vertex.color = ifelse(V(sub_net)$name == focal, "tomato", "skyblue"), # highlight focal
  edge.arrow.size = 0.1,
  edge.color = adjustcolor("grey40", alpha.f = 0.5),
  edge.curved = 0.15,
  margin = 0,
  asp = 0,
  rescale = F,
  xlim = c(-4,4.1),
  ylim = c(-4.5,4.5)
)

# look at metaweb
# define focal node
focal <- "Fucus distichus"

# get all neighbors (both in and out)
nbrs <- neighbors(metaweb, focal, mode = "all")

# include the focal node itself
nodes_to_keep <- c(focal, V(metaweb)[nbrs]$name)

# make subgraph
sub_net <- induced_subgraph(metaweb, vids = nodes_to_keep)

# layout: keep it tidy with KK
kk_layout <- layout_with_kk(sub_net)
#kk_layout <- norm_coords(kk_layout, xmin = -2, xmax = 2, ymin = -2, ymax = 2)

# plot
plot(
  sub_net,
  layout = kk_layout,
  vertex.size = 20,
  vertex.label = NA,
  vertex.color = ifelse(V(sub_net)$name == focal, "tomato", "skyblue"), # highlight focal
  edge.arrow.size = 0.1,
  edge.color = adjustcolor("grey40", alpha.f = 0.5),
  edge.curved = 0.15,
  margin = 0,
  asp = 0,
  rescale = F,
  xlim = c(-4,4.1),
  ylim = c(-4.5,4.5)
)

# fit models ----
model <- lm(cv ~ mean_in_deg + connectance + relative_ascendancy, data = network_stats, na.action = na.omit)
summary(model)






  
  
  
