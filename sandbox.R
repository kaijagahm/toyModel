# Parameters

# n = number of individuals
# m = number of individuals who die (mortality)
# r = proportion of dyads removed per day
# c = proportion of dyads created per day
# Removal and connection should depend just on the previous day (maybe later we can allow memory beyond that)

library(igraph)

# Loss and gain of edges --------------------------------------------------
generateAMs <- function(n, density, seed, r, c, days){
  # Set the seed
  set.seed(seed)
  initialConnections <- sample(0:1, size = n*n, replace = T, prob = c(1-density, density))
  
  # Create the initial adjacency matrix
  am_init <- matrix(data = initialConnections, nrow = n, ncol = n)
  
  # For each day, remove r connections and create c connections
  ams <- vector(mode = "list", length = days) # create a list to store the adjacency matrices
  ams[[1]] <- am_init
  for(i in 2:days){
    mat_prev <- ams[[i-1]]
    mat_new <- mat_prev
    for(col in 1:ncol(mat_prev)){
      for(row in 1:nrow(mat_prev)){
        if(mat_prev[row, col] == 0){
          mat_new[row, col] <- sample(0:1, size = 1, replace = TRUE, prob = c(1-c, c))
        }else{
          mat_new[row, col] <- sample(0:1, size = 1, replace = TRUE, prob = c(1-r, r))
        }
      }
    }
    ams[[i]] <- mat_new # assign new matrix to the corresponding place in the list.
  }
  
  return(ams)
}

ams <- generateAMs(n = 20, density = 0.1, seed = 3, r = 0.1, c = 0.05, days = 10)
graphs <- lapply(ams, graph_from_adjacency_matrix, weighted = NULL, mode = "undirected", diag = FALSE)
lapply(graphs, plot)
