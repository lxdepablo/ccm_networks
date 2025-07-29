# load libraries
library(tidyverse)
library(igraph)

# set working directory
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("/projects/lude8513/ccm_networks/code/")

# source helper functions ----
source("edm_utils.R")

# read in data ----
xmaps_raw <- read_csv("../data/xmaps.csv")

# visualize xmaps
ggplot(data = filter(xmaps_raw, site == 1), aes(x = LibSize, y = skill, col = xmap)) + 
  geom_line() +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = "L", y = "Correlation")

# filter valid xmaps ----
valid_xmaps <- bind_rows(lapply(unique(xmaps_raw$site), function(s){
  curr_site <- filter(xmaps_raw, site == s)
  
  curr_site_valid <- filter_xmaps(curr_site) %>%
    mutate(site = s)
}))

# plot correlations again
ggplot(data = filter(valid_xmaps, site == 1), aes(x = LibSize, y = skill, col = xmap)) + 
  geom_line() +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = "L", y = "Correlation")

# build network --------
# get correlation from largest library size
edge_lists <- bind_rows(lapply(unique(valid_xmaps$site), function(s){
  curr_site <- filter(valid_xmaps, site == s)
  
  final_cors <- bind_rows(lapply(unique(curr_site$xmap), function(x){
    curr_xmap <- curr_site %>%
      filter(xmap == x)
    final_cor <- curr_xmap %>%
      filter(LibSize == max(unique(LibSize)))
  }))
  
  # build edge list
  edge_list <- bind_rows(lapply(1:nrow(final_cors), function(i){
    curr_row <- final_cors[i, ]
    # separate xmap column into two nodes
    split <- strsplit(curr_row$xmap, split = ":")[[1]]
    node_1 <- split[1]
    node_2 <- split[2]
    
    data.frame(sp1 = node_1, sp2 = node_2, weight = curr_row$skill)
  })) %>%
    # remove self edges, which show species temporal autocorrelation (not meaningful)
    filter(sp1 != sp2) %>%
    mutate(site = s)
}))

# write edgelist to csv
write_csv(edge_lists, "../data/edge_lists.csv")








