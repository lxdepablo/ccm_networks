# generate and save figures summarizing the pipeline's results: CCM
# convergence, the three causal networks (raw pairwise / bootstrap-pruned /
# multivariate-pruned) and how their edges/stats compare, S-map interaction
# strengths, and the causal metaweb vs. trophic metaweb. Reads the CSVs
# written by generate_xmaps.R, xmap_analysis.R, do_smap.R,
# network_analysis.R, and compare_networks.R - run this after those, not
# instead of them.

# load libraries
library(tidyverse)
library(igraph)

# set working directory
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("/projects/lude8513/ccm_networks/code/")

dir.create("../figures", showWarnings = FALSE)

# 1. CCM convergence curves, by site (diagnostic: are skill curves actually
#    converging with library size, or just noisy?) ----
xmaps_raw <- read_csv("../data/xmaps.csv", show_col_types = FALSE)

p_convergence <- ggplot(xmaps_raw, aes(x = LibSize, y = skill, group = xmap)) +
  geom_line(alpha = 0.25, linewidth = 0.3) +
  facet_wrap(~site) +
  theme_minimal() +
  labs(title = "CCM cross-map skill vs. library size, by site",
       subtitle = "one line per species pair/direction; convergence (flattening at high L) supports a real cross-map link",
       x = "Library size (L)", y = "Cross-map skill (rho)")
ggsave("../figures/xmap_convergence_by_site.pdf", p_convergence, width = 12, height = 8)

# 2. edge count by pruning approach and site - the headline "what did each
#    method remove" summary ----
edges_raw          <- read_csv("../data/edge_lists.csv", show_col_types = FALSE)
edges_bootstrap    <- read_csv("../data/edge_lists_bootstrap.csv", show_col_types = FALSE)
edges_multivariate <- read_csv("../data/edge_lists_multivariate.csv", show_col_types = FALSE)

approach_levels <- c("pairwise (raw)", "+ bootstrap", "+ multivariate (MXMap)")
edge_counts <- bind_rows(
  count(edges_raw, site)          %>% mutate(approach = approach_levels[1]),
  count(edges_bootstrap, site)    %>% mutate(approach = approach_levels[2]),
  count(edges_multivariate, site) %>% mutate(approach = approach_levels[3])
) %>%
  mutate(approach = factor(approach, levels = approach_levels))

p_edge_counts <- ggplot(edge_counts, aes(x = factor(site), y = n, fill = approach)) +
  geom_col(position = "dodge") +
  theme_minimal() +
  labs(title = "Edge count by pruning approach and site",
       x = "Site", y = "Number of edges", fill = "Approach")
ggsave("../figures/edge_counts_by_approach.pdf", p_edge_counts, width = 9, height = 5)

# 3. the three networks, side by side, for every site ----
plot_network <- function(edge_df, s, title) {
  site_edges <- filter(edge_df, site == s)
  if (nrow(site_edges) == 0) {
    plot.new()
    title(title)
    return(invisible(NULL))
  }
  net <- graph_from_edgelist(as.matrix(site_edges[, c("sp1", "sp2")]), directed = TRUE)
  plot(
    net,
    layout = layout_with_kk(net),
    vertex.size = 8,
    vertex.label = NA,
    vertex.color = "skyblue",
    edge.arrow.size = 0.2,
    edge.color = adjustcolor("grey40", alpha.f = 0.5),
    edge.curved = 0.15,
    main = title
  )
}

for (s in sort(unique(edges_raw$site))) {
  pdf(paste0("../figures/networks_site_", s, ".pdf"), width = 15, height = 5)
  par(mfrow = c(1, 3))
  plot_network(edges_raw, s, paste0("Site ", s, ": pairwise (raw)"))
  plot_network(edges_bootstrap, s, paste0("Site ", s, ": + bootstrap"))
  plot_network(edges_multivariate, s, paste0("Site ", s, ": + multivariate (MXMap)"))
  dev.off()
}

# 4. network stats (S, L, degree, connectance, reciprocity, ascendancy, ...)
#    by approach and site, from compare_networks.R's output ----
network_stats_comparison <- read_csv("../data/network_stats_comparison.csv", show_col_types = FALSE)

stats_long <- network_stats_comparison %>%
  mutate(approach = factor(approach, levels = approach_levels)) %>%
  pivot_longer(
    cols = c(S, L, mean_in_deg, mean_out_deg, connectance, reciprocity,
             feedback_edge_frac, relative_ascendancy),
    names_to = "metric", values_to = "value"
  )

p_stats <- ggplot(stats_long, aes(x = factor(site), y = value, fill = approach)) +
  geom_col(position = "dodge") +
  facet_wrap(~metric, scales = "free_y") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(title = "Network stats by pruning approach and site", x = "Site", y = NULL, fill = "Approach")
ggsave("../figures/network_stats_comparison.pdf", p_stats, width = 14, height = 10)

# 5. S-map interaction strengths, all sites ----
smap_coefs <- read_csv("../data/smap_coefs.csv", show_col_types = FALSE)

p_smap <- ggplot(smap_coefs, aes(x = date, y = value, col = name)) +
  geom_line() +
  facet_wrap(~site, scales = "free_x") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(title = "S-map interaction strengths by site (final pruned edges)",
       x = "Date", y = "Interaction strength")
ggsave("../figures/smap_interaction_strengths.pdf", p_smap, width = 14, height = 10)

# 6. causal metaweb vs. trophic metaweb stats ----
metaweb_stats <- read_csv("../data/metaweb_vs_causal_stats.csv", show_col_types = FALSE)

metaweb_long <- metaweb_stats %>%
  pivot_longer(cols = -network, names_to = "metric", values_to = "value")

p_metaweb <- ggplot(metaweb_long, aes(x = metric, y = value, fill = network)) +
  geom_col(position = "dodge") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Causal metaweb vs. trophic metaweb", x = NULL, y = "Value", fill = "Network")
ggsave("../figures/metaweb_vs_causal_stats.pdf", p_metaweb, width = 8, height = 6)

cat("Figures written to ../figures/\n")
