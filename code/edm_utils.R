library(furrr)
library(tictoc)

# prep one site's data into the wide, z-scored ccm_data format used by CCM,
# do_smap, and (new) the bootstrap/multivariate pruning steps. Pulled out of
# generate_xmaps.R so xmap_analysis.R can rebuild the same per-site time
# series without duplicating this wrangling (and without re-running CCM).
# Fixes a latent bug in the original inline version: survey_group is dropped
# by the earlier `select(-survey_group)` (after `date` is derived from it),
# so the later, duplicate `select(-survey_group)` would error - removed here.
prep_site_ccm_data <- function(gom_raw, sst, site) {
  curr_site_wide <- gom_raw %>%
    # use only control plots
    filter(site == !!site, plot == "C") %>%
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

  # bring in temperature data
  sst_z <- sst %>%
    mutate(date = as.Date(date)) %>%
    select(date, sst) %>%
    arrange(date) %>%
    mutate(
      sst = as.numeric(sst),
      sst_z = as.numeric(scale(sst)),
      sst_z = ifelse(is.na(sst_z), 0, sst_z)
    ) %>%
    select(date, sst = sst_z)

  curr_site_wide <- curr_site_wide %>%
    left_join(sst_z, by = "date") %>%
    mutate(sst = replace_na(sst, 0))

  curr_site_wide %>%
    ungroup() %>%
    clean_names() %>%
    arrange(date) %>%
    relocate(date)
}

# generalized lotka volterra
gLV_model <- odin::odin({
  # initial conditions
  initial(B[]) <- B0[i]
  
  # differential equations
  # dx/dt = x*(r+AX)
  AB[, ] <- A[i, j] * B[j]
  deriv(B[]) <- B[i] * (r[i] + sum(AB[i, ]))
  
  # user parameters
  n_spp <- user() # number of species
  r[] <- user() # intrinsic growth rates vector
  A[, ] <- user() # adjacency matrix
  B0[] <- user() # initial state
  
  # define lengths
  dim(r) <- n_spp
  dim(B) <- n_spp
  dim(B0) <- n_spp
  dim(A) <- c(n_spp, n_spp)
  dim(AB) <- c(n_spp, n_spp)
})

# function to create bootstrapped data, preserving temporal autocorrelation
block_bootstrap <- function(df, block_size = 10) {
  N <- nrow(df)
  # Number of blocks needed
  n_blocks <- ceiling(N / block_size)
  # Randomly choose starting indices for each block
  start_indices <- sample(seq(1, N - block_size + 1), 
                          size = n_blocks, replace = TRUE)
  
  # Build the bootstrapped data by stacking blocks
  new_data <- data.frame()
  for (s in start_indices) {
    end_i <- s + block_size - 1
    if (end_i > N) end_i <- N
    new_data <- rbind(new_data, df[s:end_i, ])
  }
  # If we overshoot, trim back to original length
  new_data <- new_data[1:N, ]
  
  return(new_data)
}

# Core network stats for a single directed graph
calc_network_stats <- function(network){
  S <- vcount(network)
  L <- ecount(network)
  
  mean_in_degree  <- if (S > 0) mean(degree(network, mode = "in"))  else NA_real_
  mean_out_degree <- if (S > 0) mean(degree(network, mode = "out")) else NA_real_
  
  # Connectance (directed, excludes self-loops)
  denom <- S * (S - 1)
  connectance <- if (denom > 0) L / denom else NA_real_
  
  # Reciprocity (fraction of edges that are mutual)
  reciprocity_val <- if (L > 0) reciprocity(network, mode = "default") else NA_real_
  
  # Feedback/cycle edge fraction: edges whose endpoints are both in SCCs of size > 1
  feedback_edge_frac <- if (L == 0) NA_real_ else {
    scc <- components(network, mode = "strong")
    cycV <- which(scc$csize[scc$membership] > 1)
    ed <- as_edgelist(network, names = FALSE)
    mean(ed[,1] %in% cycV & ed[,2] %in% cycV)
  }
  
  # Relative ascendancy (simple Shannon-style version on normalized adjacency)
  P <- as.matrix(as_adjacency_matrix(network, sparse = FALSE))
  tot <- sum(P)
  relative_ascendancy <- if (tot > 0 && S > 1) {
    Pn <- P / tot
    Pn[Pn == 0] <- NA
    ascendancy <- -sum(Pn * log(Pn), na.rm = TRUE)
    capacity   <- log(S^2)
    if (capacity > 0) ascendancy / capacity else NA_real_
  } else NA_real_
  
  data.frame(
    S = S,
    L = L,
    mean_in_deg  = mean_in_degree,
    mean_out_deg = mean_out_degree,
    connectance  = connectance,
    reciprocity  = reciprocity_val,
    feedback_edge_frac = feedback_edge_frac,
    relative_ascendancy = relative_ascendancy,
    row.names = NULL
  )
}

# centrality alignment between two graphs on the shared node set (spearman)
centrality_alignment <- function(g1, g2, centrality_fun = function(g) degree(g, mode = "out")){
  nms <- intersect(V(g1)$name, V(g2)$name)
  if (length(nms) < 2) return(NA_real_)
  g1s <- induced_subgraph(g1, vids = nms)
  g2s <- induced_subgraph(g2, vids = nms)
  x <- centrality_fun(g1s); x <- x[nms]
  y <- centrality_fun(g2s); y <- y[nms]
  suppressWarnings(cor(x, y, method = "spearman", use = "complete.obs"))
}

# parallelized function to calculate xmaps for every pair of species in a dataframe
par_calc_all_xmaps <- function(ccm_data, ncols = (ncol(ccm_data)-1)){
  test_cols <- colnames(ccm_data)[2:ncols]
  
  # set maximum possible libsize for any CCM run
  max_possible_libsize <- nrow(ccm_data) - 2
  
  # set up parallel session
  plan(multisession)
  
  # try CCM on every pair of species
  all_xmaps <- bind_cols(future_map(1:(length(test_cols)-1), function(i){
    bind_cols(lapply((i+1):length(test_cols), function(j){
      # find optimal embedding dimension
      e_df <- EmbedDimension(
        dataFrame = ccm_data,
        lib = c(1, nrow(ccm_data) - 2),
        pred = c(1, nrow(ccm_data) - 2),
        columns = test_cols[i],
        target = test_cols[j],
        maxE = 10
      )
      
      # select E with highest rho
      e_opt <- e_df$E[which.max(e_df$rho)]
      max_libsize <- nrow(ccm_data) - e_opt - 1
      lib_sizes <- seq(1, max_libsize, by = 1)
      
      # run CCM
      curr_xmap <- CCM(
        dataFrame = ccm_data,
        E = e_opt,
        columns = test_cols[i],
        target = test_cols[j],
        libSizes = lib_sizes,
        sample = 100,
        random = TRUE
      )
      
      # drop libSize column
      skill_vals <- curr_xmap[, -1]
      
      # pad with NA if needed
      if (nrow(skill_vals) < max_possible_libsize) {
        pad_n <- max_possible_libsize - nrow(skill_vals)
        padding <- matrix(NA, nrow = pad_n, ncol = 2)
        colnames(padding) <- colnames(skill_vals)
        skill_vals <- rbind(skill_vals, padding)
      }
      
      skill_vals
    }))
  }))
  
  plan(sequential)
  
  # reshape to long format
  all_xmaps_long <- all_xmaps[ , !grepl("\\.1$", names(all_xmaps))] %>%
    mutate(LibSize = 1:nrow(.)) %>%
    relocate(LibSize) %>%
    pivot_longer(cols = -LibSize, names_to = "xmap", values_to = "skill")
}

calc_all_xmaps <- function(ccm_data, ncols = (ncol(ccm_data)-1)){
  all_cols <- colnames(ccm_data)[2:47]
  test_cols <- all_cols[1:ncols]
  
  # set maximum possible libsize for any CCM run
  max_possible_libsize <- nrow(ccm_data) - 2
  
  # try CCM on every pair of species
  all_xmaps <- bind_cols(lapply(1:(length(test_cols)-1), function(i){
    bind_cols(lapply((i+1):length(test_cols), function(j){
      # find optimal embedding dimension
      e_df <- EmbedDimension(
        dataFrame = ccm_data,
        lib = c(1, nrow(ccm_data) - 2),
        pred = c(1, nrow(ccm_data) - 2),
        columns = test_cols[i],
        target = test_cols[j],
        maxE = 10
      )
      
      # select E with highest rho
      e_opt <- e_df$E[which.max(e_df$rho)]
      max_libsize <- nrow(ccm_data) - e_opt - 1
      lib_sizes <- seq(1, max_libsize, by = 1)
      
      # run CCM
      curr_xmap <- CCM(
        dataFrame = ccm_data,
        E = e_opt,
        columns = test_cols[i],
        target = test_cols[j],
        libSizes = lib_sizes,
        sample = 100,
        random = TRUE
      )
      
      # drop libSize column
      skill_vals <- curr_xmap[, -1]
      
      # pad with NA if needed
      if (nrow(skill_vals) < max_possible_libsize) {
        pad_n <- max_possible_libsize - nrow(skill_vals)
        padding <- matrix(NA, nrow = pad_n, ncol = 2)
        colnames(padding) <- colnames(skill_vals)
        skill_vals <- rbind(skill_vals, padding)
      }
      
      skill_vals
    }))
  }))

  # reshape to long format
  all_xmaps_long <- all_xmaps[ , !grepl("\\.1$", names(all_xmaps))] %>%
    mutate(LibSize = 1:nrow(.)) %>%
    relocate(LibSize) %>%
    pivot_longer(cols = -LibSize, names_to = "xmap", values_to = "skill")
}

# function to screen out valid xmaps using second derivative
filter_xmaps <- function(xmaps){
  # filter out xmaps with negative prediction skill
  pos_cor_list <- xmaps %>%
    group_by(xmap, site) %>%
    summarize(skill = mean(skill, na.rm = T)) %>%
    filter(skill > 0)
  xmaps_filtered <- xmaps %>%
    filter(xmap %in% pos_cor_list$xmap)

  # take second derivatives of each xmap
  second_derivs <- bind_rows(lapply(unique(xmaps_filtered$xmap), function(x){
    curr_xmap <- filter(xmaps_filtered, xmap == x)
    
    # downsample data to smooth out curves
    downsample_indices <- seq(1, nrow(curr_xmap), by = 3)
    curr_xmap <- curr_xmap[downsample_indices,]
    
    first_derivs <- lapply(2:nrow(curr_xmap), function(l){
      curr_deriv <- curr_xmap[[l, 3]] - curr_xmap[[l-1, 3]]
    })
    
    second_derivs <- bind_rows(lapply(2:length(first_derivs), function(i){
      curr_deriv <- first_derivs[[i]] - first_derivs[[i-1]]
      data.frame(xmap = x, second_deriv = curr_deriv)
    }))
  })) 
  
  avg_derivs <- second_derivs %>%
    group_by(xmap) %>%
    summarize(med_deriv = median(second_deriv, na.rm = T),
              mean_deriv = mean(second_deriv, na.rm = T))
  
  # include only relationships where second derivative is negative
  valid_xmap_list <- avg_derivs %>%
    filter(mean_deriv < -0.005 | med_deriv < -0.005) %>%
    select(xmap)
  
  valid_xmaps <- xmaps_filtered %>%
    # if correlation > 0 and ddy/dx < 0, xmap is valid
    filter(xmap %in% valid_xmap_list$xmap)
}

# ---------------------------------------------------------------------------
# Bootstrapping to prune false-positive xmaps (block_bootstrap() above was
# defined but never wired in - this is that wiring).
#
# NOTE ON EDGE ORIENTATION: par_calc_all_xmaps()/CCM() name their output
# columns "A:B" for CCM(columns = A, target = B), i.e. "use A's manifold to
# predict B". A high "A:B" score means B's dynamics leave a footprint in A -
# i.e. it is the *effect* (A) reconstructing the *cause* (B) - so the causal
# claim actually supported by a high "A:B" score is B -> A (verified
# empirically, see log.md). xmap_analysis.R now draws the edge accordingly:
# sp1 (= B = cause) -> sp2 (= A = effect), i.e. edge_lists.csv's sp1/sp2
# already reflect the direction of causality, not the raw "columns:target"
# order.
#
# Everything below therefore needs *two* different roles for an edge
# sp1 -> sp2: sp1/sp2 as literal causal-graph endpoints (cause/effect) for
# finding intermediate/mediator nodes by graph topology, versus from_col/
# to_col for the actual Simplex calls, which must reconstruct the cause
# (to_col) from the effect's manifold (from_col) - i.e. from_col = sp2,
# to_col = sp1. Getting this backwards (from_col = sp1, to_col = sp2) would
# silently retest the wrong direction now that sp1/sp2 mean cause/effect
# instead of columns/target.
# ---------------------------------------------------------------------------

# get optimal embedding dimension for one variable pair (mirrors the E search
# already used in par_calc_all_xmaps, exposed standalone so it can be re-run
# for just the small set of edges that survive filtering)
get_optimal_E <- function(data, from_col, to_col, maxE = 10, lib = NULL, pred = NULL) {
  if (is.null(lib)) lib <- c(1, nrow(data) - 2)
  if (is.null(pred)) pred <- lib

  e_df <- EmbedDimension(
    dataFrame = data,
    lib = lib,
    pred = pred,
    columns = from_col,
    target = to_col,
    maxE = maxE
  )

  e_df$E[which.max(e_df$rho)]
}

# univariate cross-map reconstruction: reconstruct `to_col` using the delay
# embedding of `from_col` (embedded = FALSE lets rEDM build the E-dim lag
# embedding itself). Tp = 0 (contemporaneous) matches classic CCM.
# knn defaults to rEDM's own E+1 when left at 0; MXMap's own experiments use
# knn = 10 (see log.md), which is noticeably more robust for the chained
# reconstructions multi_pcm() does below, at the cost of needing a large
# enough library to have 10 neighbors to draw on.
xmap_reconstruct <- function(data, from_col, to_col, E, tau = 1, lib = NULL, pred = NULL, knn = 0) {
  lib_str <- if (is.null(lib)) paste(1, nrow(data)) else paste(lib, collapse = " ")
  pred_str <- if (is.null(pred)) lib_str else paste(pred, collapse = " ")

  Simplex(
    dataFrame = data,
    lib = lib_str,
    pred = pred_str,
    E = E,
    tau = -abs(tau),
    Tp = 0,
    columns = from_col,
    target = to_col,
    embedded = FALSE,
    knn = knn
  )
}

# multivariate cross-map reconstruction: `from_cols` must already be columns
# present in `data` (e.g. reconstructed intermediate series) forming the
# state vector as-is, one dimension per column (embedded = TRUE means rEDM
# does not build any further lags - see multi_pcm()'s "short series"
# adaptation note in log.md).
multi_xmap_reconstruct <- function(data, from_cols, to_col, lib = NULL, pred = NULL, knn = 0) {
  lib_str <- if (is.null(lib)) paste(1, nrow(data)) else paste(lib, collapse = " ")
  pred_str <- if (is.null(pred)) lib_str else paste(pred, collapse = " ")

  Simplex(
    dataFrame = data,
    lib = lib_str,
    pred = pred_str,
    E = length(from_cols),
    Tp = 0,
    columns = paste(from_cols, collapse = " "),
    target = to_col,
    embedded = TRUE,
    knn = knn
  )
}

# partial correlation of x and y, controlling for z
partial_cor <- function(x, y, z) {
  ok <- complete.cases(x, y, z)
  x <- x[ok]; y <- y[ok]; z <- z[ok]
  if (length(x) < 4) return(NA_real_)

  rxy <- suppressWarnings(cor(x, y))
  rxz <- suppressWarnings(cor(x, z))
  ryz <- suppressWarnings(cor(y, z))
  denom <- sqrt((1 - rxz^2) * (1 - ryz^2))
  if (is.na(denom) || denom <= 0) return(NA_real_)

  (rxy - rxz * ryz) / denom
}

# Multivariate Partial Cross Mapping (multiPCM), after Zhang et al. 2025
# (MXMap), Section 3.2 / Eq. 7.
#
# Tests whether an edge is direct or fully explained by an indirect path
# through `conds` (the other node(s) already sitting on a 2-hop path in the
# current graph). Note the parameter names here are in *Simplex* terms, not
# graph-edge terms: to_col is the variable being reconstructed (the cause,
# in causal-graph terms) and from_col is the manifold/columns variable whose
# reconstruction is being tested (the effect) - see prune_indirect_edges()
# for how a graph edge cause -> effect maps onto this call.
#
# This is a genuine two-hop composition (NOT "reconstruct to_col directly
# from the true conds values", which would over-prune real direct links
# whenever conds also happens to correlate with to_col):
#   1. apparent:    to_col reconstructed from from_col's manifold
#   2. conds "seen" from_col's manifold: each intermediate reconstructed from
#      from_col (this is the X2 -> Conds hop in Eq. 7)
#   3. conditioned:  to_col reconstructed from that *reconstructed* conds
#      block (the Conds -> X1 hop)
#   4. rho_direct = partial correlation of (to_col, apparent-reconstruction)
#      controlling for (conditioned-reconstruction)
#
# ADAPTATION FOR SHORT ECOLOGICAL SERIES (see log.md): the paper stacks a
# full E-dim lag embedding per conditioning variable (multiSSR, Eq. 6). With
# ~26 timepoints per site that blows up the conditioned manifold's
# dimensionality far past what the library can support, so here each
# conditioning variable contributes a single (unlagged) reconstructed value
# to the conditioned block instead of a full E-dim lag stack. This keeps the
# conditioned embedding dimension equal to length(conds) regardless of E.
#
# CALIBRATION CAVEAT (see log.md): validated on simulated chain systems, this
# two-hop composition reliably ranks a genuinely indirect edge's ratio below
# a genuinely direct edge's ratio, but the *absolute* ratio values (and thus
# how well the paper's own gamma_threshold = 0.45 transfers) depend on
# coupling strength, series length, and knn - the paper itself notes the
# threshold needs re-calibrating per system (Appendix C). knn defaults to 10
# (the paper's own choice) since chained reconstruction is noise-sensitive.
multi_pcm <- function(data, from_col, to_col, conds, E = NULL, tau = 1,
                       lib = NULL, pred = NULL, knn = 10) {
  conds <- setdiff(unique(conds), c(from_col, to_col))
  if (length(conds) == 0) {
    return(list(rho_all = NA_real_, rho_direct = NA_real_, ratio = NA_real_,
                decision = "keep", reason = "no intermediate nodes"))
  }

  if (is.null(E)) E <- max(3, length(conds))

  # 1. apparent cross map: from_col -> to_col
  # (Simplex()'s output always names its first column after whatever the
  # index column in `data` was called, e.g. "date" - not literally "time" -
  # so we address it positionally via [[1]] rather than by name.)
  apparent <- xmap_reconstruct(data, from_col, to_col, E = E, tau = tau, lib = lib, pred = pred, knn = knn)
  rho_all <- abs(suppressWarnings(cor(apparent$Observations, apparent$Predictions, use = "complete.obs")))

  # 2. reconstruct each intermediate from from_col's manifold (X2 -> Conds)
  conds_recon <- lapply(conds, function(cvar) {
    xmap_reconstruct(data, from_col, cvar, E = E, tau = tau, lib = lib, pred = pred, knn = knn)$Predictions
  })
  names(conds_recon) <- conds

  # Simplex treats a dataFrame's *first* column as the time/index column, so
  # that has to come first here too - otherwise it silently swallows one of
  # the conds columns as if it were the index.
  recon_df <- data.frame(time_idx = seq_along(apparent[[1]]))
  recon_df <- cbind(recon_df, as.data.frame(conds_recon))
  recon_df$to_col_true <- apparent$Observations

  # 3. reconstruct to_col from the *reconstructed* conds block (Conds -> X1)
  conditioned <- multi_xmap_reconstruct(recon_df, from_cols = conds, to_col = "to_col_true", knn = knn)

  # 4. partial correlation, aligned on time (both apparent/conditioned were
  # built from the same recon_df row order, so a positional index lines them
  # up correctly even though Simplex may drop leading/trailing NA rows)
  merged <- merge(
    data.frame(time_idx = seq_along(apparent[[1]]), obs = apparent$Observations, pred_direct = apparent$Predictions),
    data.frame(time_idx = conditioned[[1]], pred_conditioned = conditioned$Predictions),
    by = "time_idx"
  )
  rho_direct <- abs(partial_cor(merged$obs, merged$pred_direct, merged$pred_conditioned))

  ratio <- if (is.na(rho_direct) || is.na(rho_all) || rho_all == 0) NA_real_ else rho_direct / rho_all

  list(rho_all = rho_all, rho_direct = rho_direct, ratio = ratio,
       decision = NA_character_, reason = NA_character_)
}

# Phase 2 of MXMap: prune edges in a (single-site) edge list that are fully
# explained by an indirect path, using multi_pcm().
#
# ADAPTATION: MXMap's Algorithm 1 conditions on *all* intermediate nodes
# along any path between a parent/child pair. For short series we cap the
# number of conditioning variables (max_conds) since the conditioned
# embedding dimension grows with the number of intermediates - see log.md.
#
# gamma_threshold defaults to the paper's own value (0.45) but should be
# recalibrated for this system/dataset (see log.md and Appendix C of the
# paper, which makes the same point) - it is exposed as a parameter, not a
# fixed constant, for exactly that reason.
prune_indirect_edges <- function(ccm_data, edge_df, E_default = NULL, tau = 1,
                                  gamma_threshold = 0.45, max_conds = 3, knn = 10) {
  edge_df$ratio <- NA_real_
  edge_df$rho_all <- NA_real_
  edge_df$rho_direct <- NA_real_
  edge_df$pruned <- FALSE
  edge_df$n_conds <- 0L

  for (i in seq_len(nrow(edge_df))) {
    # sp1/sp2 are true cause/effect (edge_lists.csv is now drawn in the
    # direction of causality - see the NOTE ON EDGE ORIENTATION above), so
    # graph topology (finding mediators) uses cause/effect directly...
    cause <- edge_df$sp1[i]
    effect <- edge_df$sp2[i]

    # intermediates: nodes k with a cause -> k edge AND a k -> effect edge
    # already present in this site's graph (2-hop path)
    cause_children <- edge_df$sp2[edge_df$sp1 == cause]
    effect_parents <- edge_df$sp1[edge_df$sp2 == effect]
    conds <- intersect(cause_children, effect_parents)
    conds <- setdiff(conds, c(cause, effect))
    if (length(conds) > max_conds) conds <- conds[seq_len(max_conds)]

    if (length(conds) == 0) next

    edge_df$n_conds[i] <- length(conds)

    # ...but multi_pcm()'s Simplex calls must reconstruct the cause from the
    # effect's manifold (from_col = effect, to_col = cause), matching how
    # this edge's apparent CCM score was actually computed.
    res <- tryCatch(
      multi_pcm(ccm_data, from_col = effect, to_col = cause, conds, E = E_default, tau = tau, knn = knn),
      error = function(e) list(rho_all = NA_real_, rho_direct = NA_real_, ratio = NA_real_)
    )

    edge_df$rho_all[i] <- res$rho_all
    edge_df$rho_direct[i] <- res$rho_direct
    edge_df$ratio[i] <- res$ratio
    edge_df$pruned[i] <- !is.na(res$ratio) && res$ratio < gamma_threshold
  }

  edge_df
}

# ---------------------------------------------------------------------------
# Block-permutation significance test for a single pairwise cross-map edge.
#
# H0: from_col and to_col are not dynamically coupled.
#
# An earlier version of this test block-*resampled* (with replacement) the
# joint (from_col, to_col) pair via block_bootstrap(). That's wrong for a
# null test in two ways, confirmed empirically against a known-independent
# pair before this version was written (see log.md): (1) resampling the pair
# jointly keeps from_col[i] paired with to_col[i] in every row, which
# preserves - rather than breaks - any real coupling, so it structurally
# cannot test "not coupled"; and (2) sampling blocks *with replacement*
# creates duplicate blocks, and a duplicated point is its own nearest
# neighbor in the simplex projection (distance 0, matching target), which
# trivially inflates cross-map skill - worse with larger blocks, which is
# exactly the monotonic inflation that was observed.
#
# This version instead permutes to_col's blocks *without* replacement
# (block_permute() - every original block is used exactly once, just
# reordered) while leaving from_col in its original order. That breaks the
# temporal alignment between the two series (destroying real coupling, if
# any) while preserving each series' own within-block autocorrelation, and
# introduces no duplicate points. The edge is significant if the observed
# rho exceeds most of this null distribution.
# ---------------------------------------------------------------------------

# permute a single vector's blocks (without replacement - a reordering, not
# a resample) so autocorrelation within each block survives but the block
# sequence itself is scrambled
block_permute <- function(x, block_size = 10) {
  N <- length(x)
  block_id <- rep(seq_len(ceiling(N / block_size)), each = block_size, length.out = N)
  blocks <- split(seq_len(N), block_id)
  new_order <- unlist(blocks[sample(length(blocks))], use.names = FALSE)
  x[new_order]
}

bootstrap_ccm_significance <- function(data, from_col, to_col, E = NULL, tau = 1,
                                        n_boot = 200, block_size = 10,
                                        alpha = 0.05, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (is.null(E)) E <- get_optimal_E(data, from_col, to_col)

  obs <- xmap_reconstruct(data, from_col, to_col, E = E, tau = tau)
  rho_obs <- suppressWarnings(cor(obs$Observations, obs$Predictions, use = "complete.obs"))

  null_rhos <- vapply(seq_len(n_boot), function(b) {
    permuted <- data.frame(
      time_idx = seq_len(nrow(data)),
      x = data[[from_col]],
      y = block_permute(data[[to_col]], block_size = block_size)
    )
    names(permuted)[2:3] <- c(from_col, to_col)
    out <- tryCatch(
      xmap_reconstruct(permuted, from_col, to_col, E = E, tau = tau),
      error = function(e) NULL
    )
    if (is.null(out)) return(NA_real_)
    suppressWarnings(cor(out$Observations, out$Predictions, use = "complete.obs"))
  }, numeric(1))
  null_rhos <- null_rhos[!is.na(null_rhos)]

  p_value <- if (length(null_rhos) >= 10) {
    (1 + sum(null_rhos >= rho_obs)) / (1 + length(null_rhos))
  } else {
    NA_real_
  }

  list(
    rho_obs = rho_obs,
    p_value = p_value,
    n_null_valid = length(null_rhos),
    significant = !is.na(p_value) && p_value < alpha
  )
}

# Apply the bootstrap significance test to every edge in a (single-site)
# edge list and flag/prune the ones that don't hold up. Meant to run
# alongside/after filter_xmaps() (convergence-based screening) as an
# additional, independent false-positive check - not a replacement.
prune_nonsignificant_edges <- function(ccm_data, edge_df, E_default = NULL, tau = 1,
                                        n_boot = 200, block_size = 10,
                                        alpha = 0.05, seed = NULL) {
  edge_df$rho_obs <- NA_real_
  edge_df$p_value <- NA_real_
  edge_df$significant <- NA
  edge_df$pruned <- FALSE

  for (i in seq_len(nrow(edge_df))) {
    # sp1/sp2 are cause/effect (see NOTE ON EDGE ORIENTATION above); the
    # Simplex call needs to reconstruct the cause from the effect's manifold,
    # i.e. from_col = effect (sp2), to_col = cause (sp1).
    cause <- edge_df$sp1[i]
    effect <- edge_df$sp2[i]

    res <- tryCatch(
      bootstrap_ccm_significance(ccm_data, from_col = effect, to_col = cause, E = E_default, tau = tau,
                                  n_boot = n_boot, block_size = block_size,
                                  alpha = alpha, seed = seed),
      error = function(e) list(rho_obs = NA_real_, p_value = NA_real_, significant = NA)
    )

    edge_df$rho_obs[i] <- res$rho_obs
    edge_df$p_value[i] <- res$p_value
    edge_df$significant[i] <- res$significant
    edge_df$pruned[i] <- isFALSE(res$significant)
  }

  edge_df
}

# find best theta value for S-map
find_theta <- function(smap_data, target, predictors, lib, pred){
  theta_grid <- seq(0, 8, by = 0.5)
  
  theta_fit <- map_dfr(theta_grid, function(th) {
    pn <- PredictNonlinear(
      dataFrame = smap_data,
      lib       = lib,
      pred      = pred,
      embedded  = TRUE,
      columns   = paste(predictors, collapse = " "),
      target    = target,
      E         = length(predictors),
      Tp        = 1,
      theta     = th
    )
    
    tibble(theta = th, rho = pn$rho)
  })
  
  best_theta <- theta_fit %>%
    filter(!is.na(rho)) %>%
    slice_max(rho, n = 1, with_ties = FALSE) %>%
    pull(theta)
}

# do S-map for one target species and its set of predictors
do_smap <- function(smap_data, target, predictors){
  # define library/pred ranges
  N <- nrow(smap_data)
  lib  <- paste0("1 ", N)
  pred <- paste0("1 ", N)
  
  # find optimal theta value
  theta <- find_theta(smap_data, target, predictors, lib, pred)
  
  # fit smap
  smap_fit <- SMap(
    dataFrame = smap_data,
    lib       = lib,
    pred      = pred,
    embedded  = TRUE,
    columns   = paste(predictors, collapse = " "),
    target    = target,
    E         = length(predictors),
    Tp        = 1,
    theta     = theta
  )
  
  # extract interaction strengths
  coefs <- smap_fit$coefficients
  
  coefs_long <- coefs %>%
    select(-C0) %>%
    pivot_longer(cols = -date)
}












