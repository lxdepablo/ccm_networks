# load libraries
library(tidyverse)
library(igraph)

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
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

# use S-map to approximate interaction strengths ----


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
    # separate xmap column into two nodes. CCM()'s "node_1:node_2" column is
    # the skill of using node_1's manifold to predict node_2 - CCM's classic
    # (and here empirically verified, see log.md) result is that a high score
    # for that direction means node_2's dynamics are recoverable from
    # node_1's manifold because node_2 causally drives node_1, i.e. the score
    # supports node_2 -> node_1, not node_1 -> node_2. So the edge is drawn
    # cause (node_2) -> effect (node_1), reversed from the raw column name.
    split <- strsplit(curr_row$xmap, split = ":")[[1]]
    node_1 <- split[1]
    node_2 <- split[2]

    data.frame(sp1 = node_2, sp2 = node_1, weight = curr_row$skill)
  })) %>%
    # remove self edges, which show species temporal autocorrelation (not meaningful)
    filter(sp1 != sp2) %>%
    mutate(site = s)
}))

# write edgelist to csv
write_csv(edge_lists, "../data/edge_lists.csv")

# prune false positives with a block-permutation significance test ----
# (see edm_utils.R's bootstrap_ccm_significance()/prune_nonsignificant_edges();
# see log.md for why this is a permutation test on to_col's blocks, not a
# resample of the pair, and for the empirical check that motivated that).
# This is the "Approach A" deliverable: pairwise CCM + false-positive pruning.
#
# NOTE ON RUNTIME: this reruns Simplex n_boot times per surviving edge, per
# site - do not run this over the full dataset locally (see claude.md). Use
# a reduced n_boot for local sanity checks; use the defaults (or larger) on
# the HPC.
edge_lists_bootstrap <- bind_rows(lapply(unique(edge_lists$site), function(s){
  print(paste0("bootstrap pruning site ", s))
  ccm_data <- read_csv(paste0("../data/ccm_data/site_", s, ".csv"), show_col_types = FALSE)
  curr_site_edges <- filter(edge_lists, site == s)

  prune_nonsignificant_edges(ccm_data, curr_site_edges, n_boot = 200, block_size = 10, alpha = 0.05) %>%
    filter(!pruned)
}))

write_csv(edge_lists_bootstrap, "../data/edge_lists_bootstrap.csv")

# multivariate cross mapping (MXMap-style) ----
# "Approach B": phase 1 (the bootstrap-pruned pairwise graph above) + phase 2,
# pruning edges that are fully explained by an indirect (2-hop) path using
# multiPCM (see edm_utils.R's multi_pcm()/prune_indirect_edges()).
#
# gamma_threshold is the paper's own default (0.45) but should be
# recalibrated for this system - see log.md. Per-site series here are far
# shorter than the paper's own validation data (~26 timepoints vs. their
# L=3500), which limits how reliable this pruning step can be; see log.md
# for the validation this is based on and that limitation.
edge_lists_multivariate <- bind_rows(lapply(unique(edge_lists_bootstrap$site), function(s){
  print(paste0("multivariate (multiPCM) pruning site ", s))
  ccm_data <- read_csv(paste0("../data/ccm_data/site_", s, ".csv"), show_col_types = FALSE)
  curr_site_edges <- filter(edge_lists_bootstrap, site == s)

  prune_indirect_edges(ccm_data, curr_site_edges, gamma_threshold = 0.45, max_conds = 3, knn = 10) %>%
    filter(!pruned)
}))

write_csv(edge_lists_multivariate, "../data/edge_lists_multivariate.csv")








