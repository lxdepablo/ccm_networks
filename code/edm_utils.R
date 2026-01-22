library(furrr)
library(tictoc)

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












