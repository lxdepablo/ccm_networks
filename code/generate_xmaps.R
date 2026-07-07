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

# make sure the per-site ccm_data directory exists (xmap_analysis.R reads
# these back for bootstrap/multivariate pruning, without needing to redo
# this data prep or re-run CCM)
dir.create("../data/ccm_data", showWarnings = FALSE)

# iterate over every site
all_sites_xmaps <- bind_rows(lapply(unique(gom_raw$site), function(s){
  print(paste0("site ", s))

  # prep timeseries data
  ccm_data <- prep_site_ccm_data(gom_raw, sst, s)

  # save this site's wide time series for reuse (bootstrap/multivariate pruning)
  write_csv(ccm_data, paste0("../data/ccm_data/site_", s, ".csv"))

  # do CCM on every pair of species
  all_xmaps_long <- par_calc_all_xmaps(ccm_data) %>%
    mutate(site = s)
}))


# write CSV's
write_csv(all_sites_xmaps, "../data/xmaps.csv")







