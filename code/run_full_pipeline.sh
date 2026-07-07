#!/bin/bash
#SBATCH -p amilan # Partition or queue - run.sh's atesting/testing is a short debug-only queue; switch to your cluster's standard compute partition/qos for a full run (adjust if this isn't Alpine or your allocation uses a different partition name)
#SBATCH --qos=normal
#SBATCH --job-name=ccm_full_pipeline # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=luis.depablo@colorado.edu
#SBATCH --nodes=1 # Only use a single node
#SBATCH --ntasks=8 # Run on n CPUs
#SBATCH --mem=32gb # Bumped up from run.sh's 8gb - bootstrap/multiPCM pruning hold extra reconstructions per edge; adjust based on observed usage
#SBATCH --time=24:00:00 # Full pipeline (pairwise CCM across all sites/species, then bootstrap + multiPCM pruning on top) is much heavier than the 20-min test run in run.sh and hasn't been benchmarked at this scale - treat as a starting point and adjust after your first run
#SBATCH --output=log_full_pipeline.out # Standard output and error log
#SBATCH --error=log_full_pipeline.err # %j inserts job number

pwd; hostname; date
echo "You've requested $SLURM_CPUS_ON_NODE core(s)."
date

module purge
module load miniforge
mamba activate ccm_env

# 1. pairwise CCM for every species pair at every site.
#    Writes data/xmaps.csv (raw cross-map skill curves) and
#    data/ccm_data/site_<n>.csv (each site's prepped wide time series -
#    needed by steps 2 and 3 below so they don't redo this data prep or
#    re-run CCM).
Rscript generate_xmaps.R
date

# 2. screen xmaps for convergence, then prune false positives:
#      - block-permutation significance test (bootstrap pruning)
#      - multiPCM indirect-edge pruning (MXMap-style multivariate CCM)
#    Writes data/edge_lists.csv (raw), data/edge_lists_bootstrap.csv
#    (+ bootstrap pruning), data/edge_lists_multivariate.csv (+ multiPCM
#    pruning on top of that) - see log.md for what each approach removes
#    and why they're not directly interchangeable.
Rscript xmap_analysis.R
date

# 3. multivariate S-map interaction strengths along the final (multiPCM-
#    pruned) edge list. Writes data/smap_coefs.csv.
Rscript do_smap.R
date

# 4. network stats/plots: causal network(s) vs. the trophic metaweb. Writes
#    data/network_stats.csv (per-site causal network stats) and
#    data/metaweb_vs_causal_stats.csv (causal metaweb vs. trophic metaweb).
Rscript network_analysis.R
date

# 5. compare the pairwise/bootstrap/multivariate networks against each
#    other - which edges each pruning step removed, and network-level stats
#    for all three. Writes data/network_stats_comparison.csv.
Rscript compare_networks.R
date

# 6. save figures summarizing everything above (CCM convergence, the three
#    networks per site, edge counts and network stats by approach, S-map
#    interaction strengths, causal vs. trophic metaweb) to figures/.
Rscript make_figures.R
date
