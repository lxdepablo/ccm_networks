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

# function to calculate network stats
calc_network_stats <- function(network){
  # number of nodes and edges
  S <- vcount(network)
  L <- ecount(network)
  
  # mean in-degree and out-degree
  mean_in_degree <- mean(degree(network, mode = "in"))
  mean_out_degree <- mean(degree(network, mode = "out"))
  
  # connectance: proportion of possible links realized
  connectance <- L / (S * (S - 1))  # excludes self-loops
  
  # relative ascendancy (simplified version)
  # from Ulanowicz: A/C = ascendancy / capacity
  P <- as_adjacency_matrix(network, attr = NULL, sparse = FALSE)
  P <- P / sum(P)  # normalize to make it a flow matrix
  P[P == 0] <- NA
  ascendancy <- -sum(P * log(P), na.rm = TRUE)
  capacity <- log(S^2)
  rel_asc <- ascendancy / capacity
  
  data.frame(
    mean_in_deg = mean_in_degree,
    mean_out_deg = mean_out_degree,
    connectance = connectance,
    relative_ascendancy = rel_asc
  )
}

# function to calculate xmaps for every pair of species in a dataframe
par_calc_all_xmaps <- function(ccm_data, ncols = (ncol(ccm_data)-1)){
  test_cols <- colnames(ccm_data)[1:ncols]
  
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
















