# compare the pairwise-only, bootstrap-pruned, and multivariate (MXMap-style)
# causal networks produced by xmap_analysis.R
#
# bootstrap pruning and multiPCM pruning catch *different* kinds of false
# positives (see log.md): bootstrap removes edges that aren't reliably
# distinguishable from an uncoupled null; multiPCM removes edges that are
# real but fully explained by an indirect path through another node. This
# script reports both the network-level stats and which specific edges each
# step removed, rather than just edge counts.

# load libraries
library(igraph)
library(tidyverse)

# set working directory
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("/projects/lude8513/ccm_networks/code/")

# source helper functions
source("edm_utils.R")

# read in the three edge lists produced by xmap_analysis.R
edges_raw          <- read_csv("../data/edge_lists.csv", show_col_types = FALSE)
edges_bootstrap    <- read_csv("../data/edge_lists_bootstrap.csv", show_col_types = FALSE)
edges_multivariate <- read_csv("../data/edge_lists_multivariate.csv", show_col_types = FALSE)

edge_key <- function(df) paste(df$site, df$sp1, df$sp2, sep = "||")

# edges removed at each step, by site ----
removed_by_bootstrap <- edges_raw %>%
  mutate(key = edge_key(edges_raw)) %>%
  filter(!(key %in% edge_key(edges_bootstrap))) %>%
  select(site, sp1, sp2, weight)

removed_by_multipcm <- edges_bootstrap %>%
  mutate(key = edge_key(edges_bootstrap)) %>%
  filter(!(key %in% edge_key(edges_multivariate))) %>%
  select(site, sp1, sp2, weight)

cat("Edges removed by bootstrap significance pruning (Approach A):", nrow(removed_by_bootstrap), "\n")
print(removed_by_bootstrap)

cat("\nEdges additionally removed by multiPCM indirect-path pruning (Approach B):", nrow(removed_by_multipcm), "\n")
print(removed_by_multipcm)

# network-level stats per site, for each of the three networks ----
network_stats_by_approach <- bind_rows(lapply(unique(edges_raw$site), function(s){
  build_stats <- function(edge_df, approach) {
    site_edges <- filter(edge_df, site == s)
    if (nrow(site_edges) == 0) {
      return(data.frame(site = s, approach = approach, S = 0, L = 0))
    }
    net <- graph_from_edgelist(as.matrix(site_edges[, c("sp1", "sp2")]), directed = TRUE)
    calc_network_stats(net) %>%
      mutate(site = s, approach = approach)
  }

  bind_rows(
    build_stats(edges_raw, "pairwise (raw)"),
    build_stats(edges_bootstrap, "pairwise + bootstrap"),
    build_stats(edges_multivariate, "multivariate (MXMap-style)")
  )
}))

cat("\nNetwork stats by approach and site:\n")
print(network_stats_by_approach)

write_csv(network_stats_by_approach, "../data/network_stats_comparison.csv")
